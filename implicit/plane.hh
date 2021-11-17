// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#pragma once

#include "implicit.hh"

namespace I_patch
{

class Iloft;
class Liming;

struct Plane : public ImplicitSurface {
	Plane() : p(0, 0, 0), n(1, 0, 0) { }
	Plane(const Point3D &p, const Vector3D &n) : p(p), n(n) { this->n.normalize(); }

	Point3D p;
	Vector3D n;

private:
	friend Iloft;
	friend Liming;

	double eval(double x, double y, double z) const override;
	Dual eval(const Dual &x, const Dual &y, const Dual &z) const override;

	template<typename T>
	T evalImpl(const T &x, const T &y, const T &z) const;
};

#include "plane.hpp"

}