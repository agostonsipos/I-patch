// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-

#include "plane.hh"

namespace I_patch
{

double Plane::eval(double x, double y, double z) const {
	return evalImpl<double>(x, y, z);
}

Plane::Dual Plane::eval(const Dual &x, const Dual &y, const Dual &z) const {
	return evalImpl<Dual>(x, y, z);
}

Plane::Dual2 Plane::eval(const Dual2& x, const Dual2& y, const Dual2& z) const {
	return evalImpl<Dual2>(x, y, z);
}

}