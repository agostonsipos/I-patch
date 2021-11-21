// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#include "implicit.hh"

#include <exception>

using namespace Geometry;

namespace I_patch
{

constexpr double projection_tolerance = 1e-10;
constexpr size_t projection_max_iteration = 20000;
constexpr size_t projection_bisection_iteraton = 20;

Point3D ImplicitSurface::projectToSurf(const Point3D &p, double step) const
{
	return projectToSurf(p, (-grad(p)*evaluateFunction(p)).normalize(), step);
}


Point3D ImplicitSurface::projectToSurf(const Point3D &p, const Vector3D &dir, double step, double isolevel) const
{
	auto func = [this, isolevel](const Point3D& p){ return evaluateFunction(p) - isolevel; };
	// Check if it is already within tolerance
	double val = func(p);
	auto grad = evaluateGradient(p);
	if (std::abs(val) <= grad.norm() * projection_tolerance)
		return p;

	// Search for a point on the other side
	double t = 0;
	bool found = false;
	for (size_t i = 1; i < projection_max_iteration; ++i, t += step)
		if (func(p + dir * t) * val <= 0) {
			found = true;
			break;
		}
	if (!found)
		throw std::runtime_error("point cannot be projected");

	// Search for a zero point by bisection
	double lower = t - step, upper = t;
	for (size_t i = 0; i < projection_bisection_iteraton; ++i)
	{
		double middle = (lower + upper) / 2;
		if (func(p + dir * middle) * val >= 0)
			lower = middle;
		else
			upper = middle;
	}
	return p + dir * (lower + upper) / 2;
}

double ImplicitSurface::evaluateFunction(const Point3D &p) const {
	return eval(p[0], p[1], p[2]);
}

Vector3D ImplicitSurface::evaluateGradient(const Point3D &p) const {
	using namespace autodiff;
	auto f = [&](const dual &x, const dual &y, const dual &z) { return eval(x, y, z); };
	dual x = p[0], y = p[1], z = p[2];
	return {
		derivative(f, wrt(x), at(x, y, z)),
		derivative(f, wrt(y), at(x, y, z)),
		derivative(f, wrt(z), at(x, y, z))
	};
}

Matrix3x3 ImplicitSurface::evaluateHessian(const Point3D& p) const {
	using namespace autodiff;
	auto f = [&](const dual2nd& x, const dual2nd& y, const dual2nd& z) { return eval(x, y, z); };
	dual2nd x = p[0], y = p[1], z = p[2];
	std::vector<double> data{
		derivative<2>(f, wrt(x, x), at(x, y, z)),
		derivative<2>(f, wrt(x, y), at(x, y, z)),
		derivative<2>(f, wrt(x, z), at(x, y, z)),
		derivative<2>(f, wrt(y, x), at(x, y, z)),
		derivative<2>(f, wrt(y, y), at(x, y, z)),
		derivative<2>(f, wrt(y, z), at(x, y, z)),
		derivative<2>(f, wrt(z, x), at(x, y, z)),
		derivative<2>(f, wrt(z, y), at(x, y, z)),
		derivative<2>(f, wrt(z, z), at(x, y, z)),
	};
	return Matrix3x3(data.data());
}

}