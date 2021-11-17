// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#pragma once

#include <iostream>
#include "autodiff/forward/dual.hpp"

#include "geometry.hh"

namespace I_patch
{

enum class ApproxError { Sum = 1, Max = 2 };

class Ipatch;

class ImplicitSurface
{
public:
	using Point3D = Geometry::Point3D;
	using Vector3D = Geometry::Vector3D;
	using PointVector = Geometry::PointVector;

	double operator()(const Point3D &p) const { return evaluateFunction(p); }
	Vector3D grad(const Point3D &p) const { return evaluateGradient(p); }
	Point3D projectToSurf(const Point3D &p, double step = 1) const;
	Point3D projectToSurf(const Point3D &p, const Vector3D &dir, double step = 1, double isolevel = 0.0) const;

	double evaluateFunction(const Point3D &p) const;
	Vector3D evaluateGradient(const Point3D &p) const;

protected:
	friend Ipatch;

	using Dual = autodiff::dual;
	virtual Dual eval(const Dual &x, const Dual &y, const Dual &z) const = 0;
	virtual double eval(double x, double y, double z) const = 0;
};

}