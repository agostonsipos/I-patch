// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#pragma once

#include "geometry.hh"

#include "plane.hh"
#include "ipatch.hh"
#include "piecewise.hh"

#include <map>

namespace I_patch
{

Geometry::Point3D
findRefPoint(const Geometry::PointVector &points,
			 const std::vector<std::shared_ptr<Geometry::CurveType>> &ipatch);

Geometry::PiecewiseCurve
traceMeshPlanarBoundary(const Geometry::TriMesh &mesh,
						const Geometry::Point3D &p1, const Geometry::Point3D &p2,
						const Plane &bounding);

Geometry::PiecewiseCurve
traceMeshCurvedBoundary(const Geometry::TriMesh &mesh,
						const Geometry::Point3D &p1, const Geometry::Point3D &p2,
						const Geometry::Vector3D &t1, const Geometry::Vector3D &t2);

std::shared_ptr<Geometry::CurveType>
boundaryFromPlane(bool &curvedBounding, const Geometry::TriMesh& mesh,
				  const Geometry::Point3D& p1, const Geometry::Point3D& p2,
				  const Geometry::Vector3D& t1, const Geometry::Vector3D& t2,
				  const Plane& bounding);

std::shared_ptr<Geometry::CurveType>
selectBetter(const std::shared_ptr<Geometry::CurveType> &c1,
			 const std::shared_ptr<Geometry::CurveType> &c2,
			 const Geometry::TriMesh &mesh,
			 const Geometry::Vector3D &b1n, const Geometry::Vector3D &b2n);

std::map<size_t, std::vector<Geometry::TriMesh::Triangle>>
triangleMapFromMesh(const Geometry::TriMesh& mesh);

bool isEdge(const Geometry::BSCurve& c, const Geometry::PointVector& meshEdge);

Geometry::PointVector collectMeshEdge(const Geometry::TriMesh& m);

Geometry::Vector3D roughMeshNormal(const Geometry::TriMesh& m, Geometry::Point3D p);

}