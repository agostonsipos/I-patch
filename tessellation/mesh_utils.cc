// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#include "mesh_utils.hh"
#include "iloft.hh"

#include <limits>
#include <set>
#include <cstring>

using namespace Geometry;

namespace I_patch
{

using SharedCurve = std::shared_ptr<CurveType>;

bool inside(const Ipatch &ipatch, const Point3D &p) {
	for (size_t i = 0; i < ipatch.N(); ++i)
	{
		if (ipatch.evalB(i, p) < 0)
			return false;
	}
	return true;
}

Point3D findRefPoint(const PointVector &points, const std::vector<SharedCurve> &boundaries)
{
	Point3D cp = Point3D(0, 0, 0);
	int N = 50;
	for (auto it : boundaries)
	{
		for(int i = 0; i < N; ++i)
			cp += it->eval(i / (double)N);
	}
	cp /= (boundaries.size() * N);
	double min = std::numeric_limits<double>::infinity();
	Point3D rp;
	for (auto p : points)
	{
		if ((p - cp).norm() < min)
		{
			rp = p;
			min = (p - cp).norm();
		}
	}
	return rp;
}

std::map<size_t, std::vector<TriMesh::Triangle>> triangleMapFromMesh(const TriMesh& mesh)
{
	std::map<size_t, std::vector<TriMesh::Triangle>> triangleMap;
	for(auto it : mesh.triangles())
	{
		for(int i = 0; i < 3; ++i)
		{
			triangleMap[it[i]].push_back(it);
		}
	}
	return triangleMap;
}

bool segmentPlaneIntersection(Plane p, Point3D p1, Point3D p2, Point3D& outP)
{
	double d1 = p(p1), d2 = p(p2);

	if (d1 * d2 > 0)
		return false;

	double t = d1 / (d1 - d2);
	outP = p1 + (p2 - p1) * t;

	return true;
}

std::vector<Point3D> trianglePlaneIntersection(Plane p,
	const TriMesh& mesh, TriMesh::Triangle tr,
	std::vector<std::pair<size_t, size_t>>& edges)
{
	Point3D triA = mesh.points()[tr[0]], triB = mesh.points()[tr[1]], triC = mesh.points()[tr[2]];
	edges.resize(0);

	std::vector<Point3D> intersections;
	Point3D isp;

	if (segmentPlaneIntersection(p, triA, triB, isp)) {
		edges.push_back({ tr[0],tr[1] });
		intersections.push_back(isp);
	}

	if (segmentPlaneIntersection(p, triB, triC, isp)) {
		edges.push_back({ tr[1],tr[2] });
		intersections.push_back(isp);
	}

	if (segmentPlaneIntersection(p, triC, triA, isp)) {
		edges.push_back({ tr[2],tr[0] });
		intersections.push_back(isp);
	}
	if (intersections.size() == 0) throw std::runtime_error("no intersection");
	return intersections;
}

bool nextTriangle(const std::map<size_t, std::vector<TriMesh::Triangle>>& triangleMap,
	TriMesh::Triangle tr, TriMesh::Triangle& nexttr,
	size_t a, size_t b)
{
	for(auto it : triangleMap.at(a))
		if(it != tr)
			for(int i = 0; i < 3; ++i)
				if(it[i] == b)
				{
					nexttr = it;
					return true;
				}
	return false;
}

bool isBoundaryPoint(const std::map<size_t, std::vector<TriMesh::Triangle>>& triangleMap,
	size_t ind, std::vector<size_t> &others)
{
	others.resize(0);
	std::multiset<size_t> neighbours;
	for (auto tr : triangleMap.at(ind))
	{
		for(int j = 0; j < 3; ++j)
			if (tr[j] == ind)
			{
				neighbours.insert(tr[(j+1)%3]);
				neighbours.insert(tr[(j+2)%3]);
			}
	}
	bool l = false;
	for (auto it : neighbours)
	{
		if (neighbours.count(it) == 1 )
		{
			l = true;
			others.push_back(it);
		}
	}
	return l;
}

std::pair<int, int> startEdge(const TriMesh& mesh,
	const std::map<size_t, std::vector<TriMesh::Triangle>>& triangleMap,
	TriMesh::Triangle tr, Vector3D dir)
{
	std::pair<int, int> startedge;
	double bestdir = -1;
	dir.normalize();
	for (int i = 0; i < 3; ++i)
	{
		std::vector<size_t> nexts;
		if(!isBoundaryPoint(triangleMap, tr[i], nexts))
			continue;
		double dir1 = ((mesh.points()[nexts[0]] - mesh.points()[tr[i]]).normalize() * dir);
		double dir2 = ((mesh.points()[nexts[1]] - mesh.points()[tr[i]]).normalize() * dir);
		if(dir1 > bestdir){
			bestdir = dir1;
			startedge = std::make_pair(tr[i], nexts[0]);
		}
		if(dir2 > bestdir){
			bestdir = dir2;
			startedge = std::make_pair(tr[i], nexts[1]);
		}
	}
	if(bestdir == -1)
		throw std::runtime_error("Triangle not on edge with curved bounding");
	return startedge;
}


std::pair<int, int> nextEdge(const std::map<size_t, std::vector<TriMesh::Triangle>>& triangleMap, std::pair<int, int> edge)
{
	size_t a = edge.first, b = edge.second;
	std::multiset<size_t> nexts;
	for (auto tr : triangleMap.at(b))
	{
		if (tr[0] == b)
		{
			nexts.insert(tr[1]);
			nexts.insert(tr[2]);
		}
		if (tr[1] == b)
		{
			nexts.insert(tr[0]);
			nexts.insert(tr[2]);
		}
		if (tr[2] == b)
		{
			nexts.insert(tr[1]);
			nexts.insert(tr[0]);
		}
	}
	for (auto it : nexts)
	{
		if (nexts.count(it) == 1 && it != a)
			return std::make_pair(b, it);
	}
	throw std::runtime_error("Could not trace mesh boundary (this should not happen)");
}

PiecewiseCurve traceMeshPlanarBoundary(const TriMesh &mesh,
												const Point3D &p1, const Point3D &p2,
												const Plane &bounding)
{
	std::map<size_t, std::vector<TriMesh::Triangle>> triangleMap = triangleMapFromMesh(mesh);

	auto firstTr = mesh.closestTriangle(p1);
	auto lastTr = mesh.closestTriangle(p2);
	auto start = mesh.projectToTriangle(p1, firstTr), end = mesh.projectToTriangle(p2, lastTr);
	start -= bounding.n * bounding(start);
	end -= bounding.n * bounding(end);

	PointVector points{ start };
	auto tr = firstTr;
	std::vector<std::pair<size_t, size_t>> edges;
	PointVector ips = trianglePlaneIntersection(bounding, mesh, tr, edges);
	bool which = (ips[0] - start) * (end - start) > 0;
	Point3D ip = which ? ips[0] : ips[1];
	std::pair<size_t, size_t> edge = which ? edges[0] : edges[1];

	TriMesh::Triangle nexttr;
	try{
		bool onmesh = nextTriangle(triangleMap, tr, nexttr, edge.first, edge.second);
		if(!onmesh)
			throw std::runtime_error("Tracing going off mesh");
	} catch(std::runtime_error& e)
	{
		if (strcmp(e.what(), "Tracing going off mesh") == 0) {
			ip = !which ? ips[0] : ips[1];
			std::pair<size_t, size_t> edge = !which ? edges[0] : edges[1];
			bool onmesh = nextTriangle(triangleMap, tr, nexttr, edge.first, edge.second);
			if(!onmesh)
				throw std::runtime_error("Tracing going off mesh");
		}
		else throw;
	}
	points.push_back(ip);
	tr = nexttr;

	while (tr != lastTr) {
		std::vector<std::pair<size_t, size_t>> edges;
		PointVector ips = trianglePlaneIntersection(bounding, mesh, tr, edges);
		bool which = (ips[0] - ip).norm() > (ips[1] - ip).norm();
		Point3D next = which ? ips[0] : ips[1];
		std::pair<size_t, size_t> edge = which ? edges[0] : edges[1];

		bool onmesh = nextTriangle(triangleMap, tr, nexttr, edge.first, edge.second);
		if(!onmesh)
			throw std::runtime_error("Tracing going off mesh");
		ip = next;
		points.push_back(ip);
		tr = nexttr;
	}
	points.push_back(end);
	return PiecewiseCurve(points);
}

PiecewiseCurve traceMeshCurvedBoundary(const TriMesh &mesh,
												const Point3D &p1, const Point3D &p2,
												const Vector3D &t1, const Vector3D &t2)
{
	std::map<size_t, std::vector<TriMesh::Triangle>> triangleMap = triangleMapFromMesh(mesh);
	
	auto firstTr = mesh.closestTriangle(p1);
	auto lastTr = mesh.closestTriangle(p2);
	
	auto startE = startEdge(mesh, triangleMap, firstTr, t1);
	auto endE = startEdge(mesh, triangleMap, lastTr, t2);
	endE = std::make_pair(endE.second, endE.first);

	PointVector points{ p1 };
	points.push_back(mesh.points()[startE.second]);
	auto edge = nextEdge(triangleMap, startE);
	while (edge != endE && edge != std::make_pair(endE.second, endE.first))
	{
		points.push_back(mesh.points()[edge.second]);
		edge = nextEdge(triangleMap, edge);
	}
	points.push_back(p2);
	return PiecewiseCurve(points);
}

SharedCurve boundaryFromPlane(bool& curvedBounding, const TriMesh& mesh, const Point3D& p1, const Point3D& p2, const Vector3D& t1, const Vector3D& t2, const Plane& bounding)
{
	SharedCurve boundary;
	try {
		boundary = std::make_shared<PiecewiseCurve>(curvedBounding
			? traceMeshCurvedBoundary(mesh, p1, p2, t1, t2)
			: traceMeshPlanarBoundary(mesh, p1, p2, bounding));
	}
	catch (std::runtime_error & e) {
		if (strcmp(e.what(), "Triangle not on edge with curved bounding") == 0) {
			boundary = std::make_shared<PiecewiseCurve>(traceMeshPlanarBoundary(mesh, p1, p2, bounding));
			curvedBounding = false;
		}
		else if (strcmp(e.what(), "Tracing going off mesh") == 0 || strcmp(e.what(), "no intersection") == 0) {
			boundary = std::make_shared<PiecewiseCurve>(traceMeshCurvedBoundary(mesh, p1, p2, t1, t2));
			curvedBounding = true;
		}
		else
			throw;
	}
	return boundary;
}

SharedCurve selectBetter(const SharedCurve &c1, const SharedCurve &c2,
						 const TriMesh &mesh, const Vector3D &b1n, const Vector3D &b2n)
{
	if (!c1) return c2;
	if (!c2) return c1;
	double x1 = std::abs(b1n * roughMeshNormal(mesh, c1->eval(0.5)));
	double x2 = std::abs(b2n * roughMeshNormal(mesh, c2->eval(0.5)));
	return x1 < x2 ? c1 : c2;
}


bool isEdge(const BSCurve& c, const PointVector& meshEdge)
{
	PointVector curvePoints {c.eval(0), c.eval(0.5), c.eval(1)};
	std::vector<bool> onEdge{false, false, false};
	double chord = (*curvePoints.rbegin() - *curvePoints.begin()).norm();
	for(size_t i = 0; i < meshEdge.size(); ++i)
	{
		const Point3D& a = meshEdge[i];
		const Point3D& b = meshEdge[(i+1)%meshEdge.size()];
		for(size_t j = 0; j < curvePoints.size(); ++j)
		{
			if(onEdge[j])
				continue;
			Point3D& p = curvePoints[j];
			Vector3D d = (b - a).normalize();
			Point3D proj = a + d*((p-a)*d);
			if((p - proj).norm() < chord / 10)
			{
				onEdge[j] = true;
			}
		}
	}
	for(auto it : onEdge)
	{
		if(!it)
			return false;
	}
	std::cerr << std::endl;
	return true;
}

PointVector collectMeshEdge(const TriMesh& m)
{
	PointVector meshEdge;
	auto triangleMap = triangleMapFromMesh(m);
	std::set<size_t> boundPointIndices;
	std::map<size_t, std::set<size_t>> boundNeighbours;
	for(size_t i = 0; i < m.points().size(); ++i)
	{
		std::multiset<size_t> neighbours;
		for (auto tr : triangleMap.at(i))
		{
			for(int j = 0; j < 3; ++j)
				if (tr[j] == i)
				{
					neighbours.insert(tr[(j+1)%3]);
					neighbours.insert(tr[(j+2)%3]);
				}
		}
		for(auto n : neighbours)
		{
			if(neighbours.count(n) == 1)
			{
				boundPointIndices.insert(i);
				boundNeighbours[i].insert(n);
			}
		}
	}
	size_t start = *boundPointIndices.begin();
	meshEdge.push_back(m.points()[start]);
	size_t ind = *boundNeighbours[start].begin(), prev = start;
	while(ind != start)
	{
		meshEdge.push_back(m.points()[ind]);
		for(auto it : boundNeighbours[ind])
		{
			if(it != prev)
			{
				prev = ind;
				ind = it;
				break;
			}
		}
	}
	return meshEdge;
}

Vector3D roughMeshNormal(const TriMesh& m, Point3D p)
{
	auto tr = m.closestTriangle(p);
	return (m.points()[tr[1]] - m.points()[tr[0]]) ^ (m.points()[tr[2]] - m.points()[tr[0]]);
}

}