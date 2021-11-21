// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#pragma once

#include "plane.hh"

namespace I_patch
{

class Liming : public ImplicitSurface
{
public:
	Liming(std::array<Plane, 3> P, double l) : planes(std::move(P)), lambda(l) { }
	Liming(std::array<Plane, 3> P, const Point3D &point);

	double fit(const PointVector &approx, bool algebraicError, ApproxError errType);

	const Plane &plane(size_t i) const { return planes[i]; }

	void invert(){ coeff *= -1; }

	void printLog();

	static void setNormalize(bool norm) { normalize = norm; }
	double getLambda() { return lambda; }

private:
	double eval(double x, double y, double z) const override;
	Dual eval(const Dual &x, const Dual &y, const Dual &z) const override;
	Dual2 eval(const Dual2& x, const Dual2& y, const Dual2& z) const override;

	template<typename T>
	T evalImpl(const T &x, const T &y, const T &z) const;

	std::array<Plane, 3> planes;
	double coeff = 1; // multiplier, basically used as a sign
	double lambda;
	static bool normalize;
};

#include "liming.hpp"

}