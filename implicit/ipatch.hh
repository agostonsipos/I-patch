// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#pragma once
#include "implicit.hh"

namespace I_patch
{

class Ipatch : public ImplicitSurface
{
public:
	Ipatch(std::vector<std::shared_ptr<const ImplicitSurface>> R,
		   std::vector<std::shared_ptr<const ImplicitSurface>> B,
		   std::vector<double> w = {})
		: R(std::move(R)), B(std::move(B)), w(std::move(w))
	{
	}

	void setInitialCoeffs(const Point3D &p);
	void setInitialCoeffs(const Point3D &p, const std::vector<double> &d);
	void setOptimizedCoeffs(const PointVector &pv, double step = 1, bool algebraic = false);

	size_t N() const { return R.size(); }
	std::shared_ptr<const ImplicitSurface> Ribbon(size_t i) const { return R[i]; }
	std::shared_ptr<const ImplicitSurface> Bounding(size_t i) const { return B[i]; }
	double weight(size_t i) const { return w[i]; }

	double evalR(size_t i, const Point3D &p) const { return R[i]->evaluateFunction(p); }
	double evalB(size_t i, const Point3D &p) const { return B[i]->evaluateFunction(p); }

	std::shared_ptr<const ImplicitSurface> getR(size_t i) const { return R[i]; }
	std::shared_ptr<const ImplicitSurface> getB(size_t i) const { return B[i]; }

	Point3D center;

	static double maxRatio;

	double eval(double x, double y, double z) const override;
	Dual eval(const Dual &x, const Dual &y, const Dual &z) const override;
	Dual2 eval(const Dual2& x, const Dual2& y, const Dual2& z) const override;

protected:
	template<typename T>
	T evalImpl(const T &x, const T &y, const T &z) const;

	const static int exponent;

	std::vector<std::shared_ptr<const ImplicitSurface>> R;
	std::vector<std::shared_ptr<const ImplicitSurface>> B;
	std::vector<double> w;
	std::vector<double> _d;
};

#include "ipatch.hpp"

}