// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#include "surface.hh"

#include "piecewise.hh"
#include "mesh_utils.hh"

#include "plane.hh"
#include "liming.hh"
#include "iloft.hh"

#include <dc.hh>
#include <mc.h>

#include <algorithm>
#include <cmath>
#include <execution>
#include <numeric>
#include <cstring>
#include <fstream>

using namespace Geometry;

namespace I_patch
{

using BoundingBox = std::array<Point3D, 2>;
using SharedCurve = std::shared_ptr<CurveType>;
using SharedSurface = std::shared_ptr<const ImplicitSurface>;

constexpr size_t boundary_resolution = 10;
constexpr double step_size = 0.023;       // multiplied by the bounding box axis
constexpr double eval_step_size = 0.0023; // multiplied by the bounding box axis
constexpr size_t trace_sampling = 50;
constexpr size_t trace_max_iterations = 2000;
constexpr size_t trace_bisection_iteration = 10;

static BoundingBox computeBoundingBox(const PointVector &p) {
	BoundingBox bbox;
	bbox[0] = bbox[1] = p[0];
	for (size_t i = 1; i < p.size(); ++i)
		for (int j = 0; j < 3; ++j) {
			if (p[i][j] < bbox[0][j])
				bbox[0][j] = p[i][j];
			if (p[i][j] > bbox[1][j])
				bbox[1][j] = p[i][j];
		}
	return bbox;
}

Surface::Surface(C0Coons coons_patch, Ipatch i_patch, std::vector<CornerPoint> corner_points):
		corners(corner_points), coons(std::make_unique<C0Coons>(coons_patch)), patch(std::make_unique<Ipatch>(i_patch))
{
	Geometry::PointVector p(corner_points.size());
	std::transform(corner_points.begin(), corner_points.end(), p.begin(), [](CornerPoint& cp){ return cp.first; });
	setBBox(computeBoundingBox(p));
}

Surface::Surface(TriMesh inputMesh, const PointVector &p, const VectorVector &n,
				 const VectorVector &B1, const VectorVector &B2, const PointVector &mps,
				 OptimizeType optimize, ApproxError err, RibbonType ribbonType,
				 BoundingType boundingType, double offset)
{
	boundingBox = computeBoundingBox(p);
	size_t N = p.size();
	for(size_t i = 0; i < N; ++i)
		corners.emplace_back(p[i], n[i]);

	// Create ribbons, boundings & boundary curves

	std::vector<SharedSurface> ribbons(N);
	std::vector<SharedSurface> boundings(N);
	std::vector<SharedCurve> boundCurvesMesh;
	std::vector<SharedCurve> boundCurvesI;

	for (size_t i = 0; i < N; ++i) {
		std::cerr << "Ribbon #" << i << std::endl;
		auto p1 = p[i], p2 = p[(i + 1) % N];
		auto n1 = n[i], n2 = n[(i + 1) % N];
		Point3D b1(0,0,0), b2(0,0,0);
		if(B2.size() > 0) b1 = B2[i];
		if(B1.size() > 0) b2 = B1[(i + 1) % N];
		Point3D mp;
		if(boundingType != BoundingType::ByNormal) mp = mps[i];
		auto nn = ((n1+n2)/2).normalize();

		Plane bounding1 = Plane(p1, (p2-p1)^nn);
		Plane bounding2 = Plane(p1, (p1 - mp) ^ (p2 - mp));
		if (bounding1(p[(i + 2) % N]) < 0)
			bounding1.n *= -1;
		if (bounding2(p[(i + 2) % N]) < 0)
			bounding2.n *= -1;

		bool curvedBounding = b1.norm() > 1e-10 && b2.norm() > 1e-10;
		
		SharedCurve boundary1, boundary2;
		if(boundingType != BoundingType::ByMidpoint)
			try {
				boundary1 = boundaryFromPlane(curvedBounding, inputMesh,
											  p1, p2, b1 ^ n1, b2 ^ n2, bounding1);
			} catch (std::exception& e) {
				std::cerr << "Error while tracing a boundary: " << e.what() << std::endl;
			}
		if(boundingType != BoundingType::ByNormal)
			try {
				boundary2 = boundaryFromPlane(curvedBounding, inputMesh,
											  p1, p2, b1 ^ n1, b2 ^ n2, bounding2);
			} catch (std::exception& e) {
				std::cerr << "Error while tracing a boundary: " << e.what() << std::endl;
			}

		SharedCurve boundary;
		if (boundingType == BoundingType::ByNormal)
			boundary = boundary1;
		else if (boundingType == BoundingType::ByMidpoint)
			boundary = boundary2;
		else if (boundingType == BoundingType::Better)
			boundary = selectBetter(boundary1, boundary2, inputMesh, bounding1.n, bounding2.n);
		boundCurvesMesh.push_back(boundary);

		if(boundary.get() == nullptr)
			throw std::runtime_error("boundary could not be traced");

		auto bounding = (boundary == boundary1) ? bounding1 : bounding2;

		PointVector approxPoints;
		for (size_t j = 1; j < boundary_resolution; ++j)
			approxPoints.push_back(boundary->eval((double)j / boundary_resolution));

		auto corner1 = Plane(p1, n1), corner2 = Plane(p2, n2);
		const auto &midPoint = approxPoints[(boundary_resolution-1)/2];
		double distance = (p1 - p2).norm();

		bool concave = ((corner1.p + corner1.n * distance / 10) -
						(corner2.p + corner2.n * distance / 10)).norm() < distance;

		auto mn = roughMeshNormal(inputMesh, midPoint);
		auto limingNormal_primary = ((p2 - p1) ^ mn) ^ (p2 - p1);
		if(limingNormal_primary * nn < 0) limingNormal_primary *= -1;
		std::cerr << limingNormal_primary << std::endl;
		auto limingPlane_primary = Plane(p1, limingNormal_primary);

		Liming lr({{corner1, corner2, limingPlane_primary}}, midPoint);
		double limingError = lr.fit(approxPoints, true, err);
		if (!concave)
			lr.invert();

		bool limingBad = (lr.grad(p1) * n1) * (lr.grad(p2) * n2) < 0 ||
			lr(midPoint + nn * (p1 - p2).norm() / 5) < 0 ||
            lr(midPoint - nn * (p1 - p2).norm() / 5) > 0;

		Plane pb1 = B1.size() > 0 && B1[i].norm() > 1e-10 && Plane(p1, B1[i])(p[(i - 1 + N) % N]) > 0 ?
			Plane(p1, B1[i]) :
			Plane(p1, (p[i] - p[(i - 1 + N) % N]) ^ (n[(i - 1 + N) % N] + n[i]));
		Plane pb2 = B2.size() > 0 && B2[(i+1)%N].norm() > 1e-10 && Plane(p2, B2[(i + 1) % N])(p[(i + 2) % N]) ?
			Plane(p2, B2[(i+1)%N]) :
			Plane(p2, (p[(i + 1) % N] - p[(i + 2) % N]) ^ (n[(i + 1) % N] + n[(i + 2) % N]));

		Vector3D sweep = boundary == boundary2 ? bounding.n : (mn ^ (p2 - p1));
		Vector3D b1n = nn ^ sweep;
		b1n *= (b1n * (p1 - p2) > 0 ? 1 : -1);
		Vector3D b2n = nn ^ sweep;
		b2n *= (b2n * (p2 - p1) > 0 ? 1 : -1);
		Iloft ir({ corner1, corner2 },
				 { curvedBounding ? pb1 : Plane(p1, b1n), curvedBounding ? pb2 : Plane(p2, b2n) }, midPoint);
		double iribbonError = ir.fit(approxPoints, true, err);

		std::cerr << "Liming error: " << limingError << std::endl;
		std::cerr << "I-ribbon error: " << iribbonError << std::endl;

		if (limingBad || lr.getLambda() < 0 || ribbonType == RibbonType::Iloft || (ribbonType != RibbonType::Liming && iribbonError < limingError)) {
			ribbons[i] = std::make_shared<Iloft>(ir);
			std::cerr << "I-ribbon" << std::endl;
			ir.printLog();
		} else {
			ribbons[i] = std::make_shared<Liming>(lr);
			std::cerr << "Liming-ribbon" << std::endl;
			lr.printLog();
		}
		std::cerr << "\tcenter: " << midPoint << std::endl;

		std::cerr << "Bounding #" << i << std::endl;
		if (curvedBounding) {
			Plane pl1(p1, b1), pl2(p2, b2), pl3(p1, ((p2 - p1) ^ nn).normalize());
			Liming lb({ pl1, pl2, pl3 }, midPoint);
			double limErr = lb.fit(approxPoints, true, err);
			if (((p1 + b1 * distance / 10) - (p2 + b2 * distance / 10)).norm() < distance)
				lb.invert();

			Iloft ib({ pl1,pl2 }, { pb1,pb2 }, midPoint);
			double ibErr = ib.fit(approxPoints, true, err);

			std::cerr << "Liming error: " << limErr << std::endl;
			std::cerr << "I-loft error: " << ibErr << std::endl;

			bool l = true;
			for (auto q : p)
				if (lb(q) < 0) l = false;
			if (!(ribbonType == RibbonType::Iloft) && l && ibErr > limErr) {
				boundings[i] = std::make_shared<Liming>(lb);
				std::cerr << "Curved bounding (Liming)" << std::endl;
				lb.printLog();
			}
			else {
				boundings[i] = std::make_shared<Iloft>(ib);
				std::cerr << "Curved bounding (I-loft)" << std::endl;
				ib.printLog();
			}
		}
		else {
			boundings[i] = std::make_shared<Plane>(bounding);
			std::cerr << "Planar bounding" << std::endl;
		}

		PointVector points = traceBoundaryCurve(ribbons[i], p1 + n1*offset, p2 + n2*offset, offset, boundary == boundary2 ? (p2-p1)^bounding.n : nn);

		if (curvedBounding && offset == 0.0)
			boundCurvesI.push_back(boundary);
		else
			boundCurvesI.push_back(std::make_shared<PiecewiseCurve>(points));
		std::cerr << std::endl;
	}

	coons = std::make_unique<C0Coons>(boundCurvesI);
	patch = std::make_unique<Ipatch>(ribbons, boundings);

	auto pointCloud = filterPoints(inputMesh.points());
	std::cerr << "Number of points: " << pointCloud.size() << std::endl;
	Point3D centerPoint = findRefPoint(pointCloud, boundCurvesMesh);
	std::cerr << "Center point: " << centerPoint << std::endl;
	patch->setInitialCoeffs(centerPoint);

	std::cerr << "Values in center point" << std::endl;
	for (size_t i = 0; i < patch->N(); ++i)
	{
		std::cerr << "R_" << i << "(c) = " << patch->evalR(i, centerPoint)
			<< " (normalized";
		if (!Iloft::getNormalize())
			std::cerr << ": " << patch->evalR(i, centerPoint) / patch->getR(i)->grad(centerPoint).norm();
		std::cerr << ")" << std::endl;

		std::cerr << "B_" << i << "(c) = " << patch->evalB(i, centerPoint)
			<< " (normalized";
		if(!Iloft::getNormalize())
			std::cerr << ": " << patch->evalB(i, centerPoint) / patch->getB(i)->grad(centerPoint).norm();
		std::cerr << ")" << std::endl;
	}

	std::cerr << "Initial coeffs: " << std::endl;
	for (size_t i = 0; i <= patch->N(); ++i)
	{
		std::cerr << patch->weight(i) << " ";
	}
	std::cerr << std::endl;

	if(optimize != OptimizeType::None)
		patch->setOptimizedCoeffs(pointCloud, boundingBoxSize() * step_size, optimize == OptimizeType::Algebraic);
}

Vector3D HSVtoRGB(double H, double S, double V) {
	if (H > 360 || H < 0 || S>1 || S < 0 || V>1 || V < 0) {
		throw std::runtime_error("The given HSV values are not in valid range");
	}
	double C = S * V;
	double X = C * (1 - std::abs(std::fmod(H / 60.0, 2) - 1));
	double m = V - C;
	double r, g, b;
	if (H >= 0 && H < 60) {
		r = C, g = X, b = 0;
	}
	else if (H >= 60 && H < 120) {
		r = X, g = C, b = 0;
	}
	else if (H >= 120 && H < 180) {
		r = 0, g = C, b = X;
	}
	else if (H >= 180 && H < 240) {
		r = 0, g = X, b = C;
	}
	else if (H >= 240 && H < 300) {
		r = X, g = 0, b = C;
	}
	else {
		r = C, g = 0, b = X;
	}
	double R = (r + m);
	double G = (g + m);
	double B = (b + m);
	return Vector3D(R, G, B);
}

void Surface::exportColoredMesh(std::string filename)
{
	TriMesh mesh = eval(50, 0.0);
	std::ofstream file(filename);
	file << "ply" << std::endl <<
		"format ascii 1.0" << std::endl <<
		"element vertex " << mesh.points().size() << std::endl <<
		"property float x" << std::endl <<
		"property float y" << std::endl <<
		"property float z" << std::endl <<
		"property uchar red" << std::endl <<
		"property uchar green" << std::endl <<
		"property uchar blue" << std::endl <<
		"element face " << mesh.triangles().size() << std::endl << 
		"property list uchar int vertex_index" << std::endl <<
		"end_header" << std::endl;
	std::vector<Vector3D> colors(patch->N());
	for (size_t i = 0; i < patch->N(); ++i)
	{
		int step = (patch->N() % 2 != 0) ? (patch->N() / 2) : ((patch->N() % 3 != 0) ? (patch->N() / 3) : 1);
		colors[i] = HSVtoRGB((step*i % patch->N()) * 360.0 / patch->N(), 1., 1.);
	}
	for (auto it : mesh.points())
	{
		std::vector<double> blends;
		for (size_t i = 0; i < patch->N(); ++i)
		{
			blends.push_back(patch->weight(i) / std::pow(patch->evalB(i, it), 2));
		}
		double sum = std::accumulate(blends.begin(), blends.end(), 0.0);
		std::for_each(blends.begin(), blends.end(), [sum](double& x) { x /= sum; });
		Vector3D color = Vector3D(0, 0, 0);
		for (size_t i = 0; i < patch->N(); ++i)
		{
			color += colors[i] * blends[i];
		}
		std::array<int, 3> col;
		std::transform(color.data(), color.data() + 3, col.begin(), [](double d) -> int { return static_cast<int>(255 * d); });
		file << it[0] << " " << it[1] << " " << it[2] << " " << col[0] << " " << col[1] << " " << col[2] << std::endl;
	}
	for (auto it : mesh.triangles())
	{
		file << 3 << " " << it[0] << " " << it[1] << " " << it[2] << std::endl;
	}
	file.close();
}

PointVector Surface::filterPoints(const PointVector &points) const {
	PointVector result;
	std::copy_if(points.begin(), points.end(), std::back_inserter(result),
				 [this](const Point3D &p) {
					 if ((p - corners[0].first).norm() > boundingBoxSize())
						 return false;
					 for (size_t i = 0; i < patch->N(); ++i)
						 if (patch->evalB(i, p) <= 0)
							 return false;
					 return true;
				 });
	return result;
}

static std::vector<double> wachspress(const Point2DVector &domain, const Point2D &uv) {
  size_t n = domain.size();
  Vector2DVector vectors; vectors.reserve(n);
  std::transform(domain.begin(), domain.end(), std::back_inserter(vectors),
                 [uv](const Point2D &p) { return uv - p; });

  DoubleVector areas; areas.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    const Vector2D &si = vectors[i];
    const Vector2D &si1 = vectors[(i+1)%n];
    areas.push_back((si[0] * si1[1] - si[1] * si1[0]) / 2.0);
  }

  DoubleVector l; l.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    size_t i_1 = (i + n - 1) % n, i1 = (i + 1) % n;
    double Ai = 1.0, Ai_1 = 1.0, Ai_1i = 1.0;
    for (size_t j = 0; j < n; ++j) {
      if (j == i)
        Ai_1 *= areas[j];
      else if (j == i_1)
        Ai *= areas[j];
      else {
        Ai_1 *= areas[j];
        Ai *= areas[j];
        Ai_1i *= areas[j];
      }
    }
    const Vector2D &si_1 = vectors[i_1];
    const Vector2D &si1 = vectors[i1];
    double Bi = (si_1[0] * si1[1] - si_1[1] * si1[0]) / 2.0;
    l.push_back(Ai_1 + Ai - Bi * Ai_1i);
  }

  double sum = std::accumulate(l.begin(), l.end(), 0.0);
  std::transform(l.begin(), l.end(), l.begin(), [sum](double x) { return x / sum; });
  return l;
}

TriMesh Surface::eval(size_t resolution, double isolevel) const
{
	TriMesh mesh = coons->meshTopology(resolution);
	Point2DVector uvs = coons->parameters(resolution);
	PointVector points(uvs.size());
	std::vector<bool> on_edge; on_edge.reserve(uvs.size());

	for (size_t i = 0; i < uvs.size(); ++i)
		on_edge.push_back(coons->onEdge(resolution, i));

	static constexpr auto policy = std::execution::seq;
	std::transform(policy, uvs.begin(), uvs.end(), on_edge.begin(), points.begin(),
		[this, isolevel](const Point2D& uv, bool on_edge) {
			Point3D p = coons->eval(uv);
			if (on_edge)
				if(std::isnan((*patch)(p)) || (*patch)(p) == 0)
					return p;
			Vector3D dir = { 0.0, 0.0, 0.0 };
			auto w = wachspress(coons->domain(), uv);
			for (size_t i = 0; i < patch->N(); ++i)
				dir += corners[i].second;
			dir *= -(patch->evaluateFunction(p) - isolevel);
			return patch->projectToSurf(p, dir.normalize(), boundingBoxSize() * eval_step_size, isolevel);
		}
	);
	mesh.setPoints(points);
	return mesh;
}


TriMesh Surface::evalMC(int resolution, double isolevel) const
{
	std::function<double(Vector3D)> imp = [this](Point3D p){ return (*patch)(p); };
	Vector3D diag = boundingBox[1] - boundingBox[0];
	return IMC::marching_cubes(imp, isolevel, { boundingBox[0] - diag, boundingBox[1] + diag }, { {resolution, resolution, resolution} }, false);
}

static TriMesh quadToTriMesh(const DualContouring::QuadMesh &mesh) {
	TriMesh result;
	size_t n = mesh.points.size();
	result.resizePoints(n);
	for (size_t i = 0; i < n; ++i) {
		const auto &p = mesh.points[i];
		result[i] = Point3D(p[0], p[1], p[2]);
	}
	for (const auto &q : mesh.quads) {
		result.addTriangle(q[0] - 1, q[1] - 1, q[2] - 1);
		result.addTriangle(q[2] - 1, q[3] - 1, q[0] - 1);
	}
	return result;
}

static TriMesh wrapDC(const std::function<double(Vector3D)> &f, double isolevel,
					  const Point3D &min, const Point3D &max, size_t resolution) {
	auto wrap = [&f](const DualContouring::Point3D &p) { return f({ p[0], p[1], p[2] }); };
	std::array<DualContouring::Point3D, 2> bbox =
		{ { { min[0], min[1], min[2] }, { max[0], max[1], max[2] } } };
	auto mesh =
		DualContouring::isosurface(wrap, isolevel, bbox, { {resolution, resolution, resolution} });
	return quadToTriMesh(mesh);
}

TriMesh Surface::evalDC(size_t resolution, double isolevel) const
{
	auto imp = [this](Point3D p){ return (*patch)(p); };
	Vector3D diag = boundingBox[1] - boundingBox[0];
	return wrapDC(imp, isolevel, boundingBox[0] - diag, boundingBox[1] + diag, resolution);
}

TriMesh Surface::evalRibbons(size_t res, double width)
{
	TriMesh ribbons;
	size_t N = patch->N();
	for (size_t i = 0; i < N; ++i) {
		auto dp = coons->domain()[i];
		auto nextdp = coons->domain()[(i + 1) % N];
		TriMesh ribbon;
		PointVector points;
		auto p1 = coons->eval(dp);
		auto p2 = coons->eval(nextdp);

		auto n1 = corners[(i + 1) % N].second;
		auto n2 = corners[(i + 2) % N].second;

		auto d1 = n1 ^ (p2 - p1);
		auto d2 = n2 ^ (p2 - p1);

		for (size_t j = 0; j <= res; ++j) {
			double r = (double)j / res;
			points.push_back(coons->eval(dp * (1 - r) + nextdp * r));
			auto d = d1 * (1 - r) + d2 * r;

			auto p = points.back() + d.normalize() * width;

			double step = boundingBoxSize() * 0.0023;
			auto point = patch->getR((i+1)%N)->projectToSurf(p, step);
			points.push_back(point);
		}
		ribbon.setPoints(points);
		for (size_t j = 0; j < res; ++j) {
			ribbon.addTriangle(2 * j, 2 * j + 1, 2 * (j + 1));
			ribbon.addTriangle(2 * (j + 1), 2 * j + 1, 2 * (j + 1) + 1);
		}
		ribbons.append(ribbon);
	}
	return ribbons;
}

TriMesh Surface::evalRibbonMC(size_t ind, int res, double isolevel)
{
	if(ind < 0 || ind >= N())
		return TriMesh();
	auto imp = [this, ind](Point3D p){ return patch->evalR(ind, p); };
	Vector3D diag = boundingBox[1] - boundingBox[0];
	return IMC::marching_cubes(imp, isolevel, { boundingBox[0] - diag, boundingBox[1] + diag }, { {res, res, res} }, false);
}

TriMesh Surface::evalRibbonDC(size_t ind, size_t res, double isolevel)
{
	if(ind < 0 || ind >= N())
		return TriMesh();
	auto imp = [this, ind](Point3D p){ return patch->evalR(ind, p); };
	Vector3D diag = boundingBox[1] - boundingBox[0];
	return wrapDC(imp, isolevel, boundingBox[0] - diag, boundingBox[1] + diag, res);
}

TriMesh Surface::evalBoundings(size_t res, double width)
{
	TriMesh boundings;
	size_t N = patch->N();
	for (size_t i = 0; i < N; ++i) {
		auto dp = coons->domain()[i];
		auto nextdp = coons->domain()[i + 1 == N ? 0 : i + 1];
		TriMesh bounding;
		PointVector points;

		auto n1 = corners[(i + 1) % N].second;
		auto n2 = corners[(i + 2) % N].second;
		auto nn = (n1 + n2) / 2;

		for (size_t j = 0; j <= res; ++j) {
			double r = (double)j / res;
			points.push_back(coons->eval(dp * (1 - r) + nextdp * r));
			points.push_back(coons->eval(dp * (1 - r) + nextdp * r) + nn * width);
		}
		bounding.setPoints(points);
		for (size_t j = 0; j < res; ++j) {
			bounding.addTriangle(2 * j, 2 * j + 1, 2 * (j + 1));
			bounding.addTriangle(2 * (j + 1), 2 * j + 1, 2 * (j + 1) + 1);
		}
		boundings.append(bounding);
	}
	return boundings;
}

TriMesh Surface::evalBoundingMC(size_t ind, int res, double isolevel)
{
	if(ind < 0 || ind >= N())
		return TriMesh();
	auto imp = [this, ind](Point3D p){ return patch->evalB(ind, p); };
	Vector3D diag = boundingBox[1] - boundingBox[0];
	auto mesh1 = IMC::marching_cubes(imp, isolevel, {boundingBox[0] - diag, boundingBox[1] + diag}, { {res, res, res} }, false);
	return mesh1;
}

TriMesh Surface::evalBoundingDC(size_t ind, size_t res, double isolevel)
{
	if(ind < 0 || ind >= N())
		return TriMesh();
	auto imp = [this, ind](Point3D p){ return patch->evalB(ind, p); };
	Vector3D diag = boundingBox[1] - boundingBox[0];
	return wrapDC(imp, isolevel, boundingBox[0] - diag, boundingBox[1] + diag, res);
}

PointVector traceBoundaryCurve(const SharedSurface &ribbon,
								const Point3D &p1, const Point3D &p2, double isoval, const Vector3D &nn)
{
	PointVector points;
	auto func = [ribbon, isoval](const Point3D& p){ return ribbon->evaluateFunction(p) - isoval; };
	points.push_back(p1);
	for (size_t j = 1; j < trace_sampling; ++j) {
		double r = (double)j / trace_sampling;
		auto p = p1 * (1 - r) + p2 * r, d = nn;
		double side = func(p);
		double step = step_size * (p2 - p1).norm() / trace_sampling;
		Point3D q;
		for (size_t u = 1; u < trace_max_iterations; ++u) {
			if (side * func(p + d * u * step) < 0) {
				q = p + d * (u - 1) * step;
				break;
			}
			if (side * func(p - d * u * step) < 0) {
				q = p - d * u * step;
				break;
			}
		}
		auto a = q, b = q + d * step, c = q;
		for (size_t t = 0; t < trace_bisection_iteration; ++t) {
			c = (a + b) / 2;
			if (func(a) * func(c) > 0)
				a = c;
			else
				b = c;
		}
		points.push_back(c);
	}
	points.push_back(p2);
	return points;
}

}