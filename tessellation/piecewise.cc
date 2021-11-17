// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#include "piecewise.hh"

namespace Geometry{

PiecewiseCurve::PiecewiseCurve(const PointVector &pt) : pts(pt)
{
	DoubleVector kn = { 0 };
	for (size_t i = 1; i < pt.size(); ++i)
	{
		kn.push_back(kn.back() + (pt[i] - pt[i - 1]).norm());
	}
	double l = kn.back();
	for (auto& it : kn)
	{
		it /= l;
	}
	knots = kn;
}

Point3D PiecewiseCurve::eval(double u) const
{
	if (u == 0) return pts.front();
	if (u == 1) return pts.back();
	for (size_t i = 0; i < knots.size(); ++i)
	{
		if (u < knots[i])
		{
			double t = (u - knots[i - 1]) / (knots[i] - knots[i - 1]);
			return pts[i - 1] * (1 - t) + pts[i] * t;
		}
	}
	return pts.back();
}

Vector3D PiecewiseCurve::evalDerivative(double u) const
{
	if (u == 1) return (pts.back() - *(pts.end() - 1)) * (pts.size() - 1);
	u *= (pts.size() - 1);
	int i = u;
	return (pts[i + 1] - pts[i]) * (pts.size() - 1);
}

}
