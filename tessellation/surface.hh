// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#pragma once

#include "c0coons.hh"
#include "ipatch.hh"

namespace I_patch
{

enum class OptimizeType { None, Projection, Algebraic };
enum class BoundingType { ByNormal = 1, ByMidpoint = 2, Better = 3 };
enum class RibbonType   { Liming = 1, Better = 2, Iloft = 3 };

// Represents a multi-sided portion of an I-patch isosurface
// `corners` are the points and normal vectors in the corners
// `coons` holds the boundary curves
class Surface
{
public:
	using CornerPoint = std::pair<Geometry::Point3D, Geometry::Vector3D>;

	// creates object from already created implicit surface and boundaries
	// only computes bounding box
	Surface(C0Coons coons_patch, Ipatch i_patch, std::vector<CornerPoint> corner_points);

	// creates patch by approximating onto a supplied mesh (see Sipos et al. (2022))
	// boundary_normals* only has to be supplied if point is we want a curved bounding surface (see Sipos et al. (2020), 4.2--4.3)
	// midpoints only has to be supplied if boundingType is not ByNormal (i.e. we want a bounding surface based on them)
	Surface(Geometry::TriMesh inputMesh, const Geometry::PointVector &corners,
			const Geometry::VectorVector &corner_normals,
			const Geometry::VectorVector &boundary_normals1,
			const Geometry::VectorVector &boundary_normals2,
			const Geometry::PointVector &midpoints, OptimizeType optimize, ApproxError err,
			RibbonType ribbonType, BoundingType boundingType, double offset = 0.0);

	size_t N() { return patch->N(); }
	double boundingBoxSize() const { return (boundingBox[0] - boundingBox[1]).norm(); }

	Geometry::TriMesh eval(size_t res, double isolevel) const;
	Geometry::TriMesh evalMC(int res, double isolevel) const;
	Geometry::TriMesh evalDC(size_t res, double isolevel) const;
	Geometry::TriMesh evalCoons(size_t res) const { return coons->eval(res); }

	Geometry::TriMesh evalRibbons(size_t res, double width);
	Geometry::TriMesh evalRibbonMC(size_t ind, int res, double isolevel);
	Geometry::TriMesh evalRibbonDC(size_t ind, size_t res, double isolevel);
	Geometry::TriMesh evalBoundings(size_t res, double width);
	Geometry::TriMesh evalBoundingMC(size_t ind, int res, double isolevel);
	Geometry::TriMesh evalBoundingDC(size_t ind, size_t res, double isolevel);

	double implicitEval(Geometry::Point3D p){ return patch->evaluateFunction(p); }

	void setBBox(const std::array<Geometry::Vector3D, 2>& bbox){ boundingBox = bbox; }

	void exportColoredMesh(std::string filename);

private:
	Geometry::PointVector filterPoints(const Geometry::PointVector &points) const;

	std::vector<CornerPoint> corners;

	std::unique_ptr<C0Coons> coons;
	std::unique_ptr<Ipatch> patch;

	std::array<Geometry::Point3D, 2> boundingBox; // AABB
};

Geometry::PointVector traceBoundaryCurve(const std::shared_ptr<const ImplicitSurface>& ribbon,
	const Geometry::Point3D& p1, const Geometry::Point3D& p2, double isoval, const Geometry::Vector3D& nn);

}