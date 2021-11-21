// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#pragma once

#include "plane.hh"

namespace I_patch
{

class Iloft : public ImplicitSurface
{
public:
	Iloft(std::array<Plane, 2> S, std::array<Plane, 2> B, std::array<double, 3> C)
		: primaries(std::move(S)), boundings(std::move(B)), coeffs(std::move(C))
	{
	}

	Iloft(std::array<Plane, 2> S, std::array<Plane, 2> B, const Point3D &midpoint);

	double fit(const PointVector &approx, bool algebraicError, ApproxError errType);
	
	void printLog();
	static void setNormalize(bool norm) { normalize = norm; }
	static bool getNormalize() { return normalize; }

	static double maxRatio;
private:
	double eval(double x, double y, double z) const override;
	Dual eval(const Dual &x, const Dual &y, const Dual &z) const override;
	Dual2 eval(const Dual2& x, const Dual2& y, const Dual2& z) const override;

	template<typename T>
	T evalImpl(const T &x, const T &y, const T &z) const;

	std::array<Plane, 2> primaries;
	std::array<Plane, 2> boundings;
	std::array<double, 3> coeffs;
	static bool normalize;
};

#include "iloft.hpp"

}