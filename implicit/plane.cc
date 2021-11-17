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

}