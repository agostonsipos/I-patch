// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#pragma once

#include "curves.hh"

namespace Geometry{

class PiecewiseCurve : public CurveType {
public:
	PiecewiseCurve(const PointVector& pts);
	virtual Point3D eval(double u) const override;
	virtual Vector3D evalDerivative(double u) const override;
private:
	PointVector pts;
	DoubleVector knots;
};

}
