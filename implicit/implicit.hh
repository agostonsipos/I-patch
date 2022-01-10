// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#pragma once

#include <iostream>
#include "autodiff/forward/dual.hpp"

#include "geometry.hh"

namespace I_patch
{

enum class ApproxError { Sum = 1, Max = 2 };

class ImplicitSurface
{
public:
	using Point3D = Geometry::Point3D;
	using Vector3D = Geometry::Vector3D;
	using Matrix3x3 = Geometry::Matrix3x3;
	using PointVector = Geometry::PointVector;

	double operator()(const Point3D &p) const { return evaluateFunction(p); }
	Vector3D grad(const Point3D &p) const { return evaluateGradient(p); }
	Point3D projectToSurf(const Point3D &p, double step = 1) const;
	Point3D projectToSurf(const Point3D &p, const Vector3D &dir, double step = 1, double isolevel = 0.0) const;

	double evaluateFunction(const Point3D &p) const;
	Vector3D evaluateGradient(const Point3D &p) const;
	Matrix3x3 evaluateHessian(const Point3D &p) const;

	using Dual = autodiff::dual;
	using Dual2 = autodiff::dual2nd;
	virtual Dual eval(const Dual& x, const Dual& y, const Dual& z) const = 0;
	virtual Dual2 eval(const Dual2& x, const Dual2& y, const Dual2& z) const = 0;
	virtual double eval(double x, double y, double z) const = 0;
};

}