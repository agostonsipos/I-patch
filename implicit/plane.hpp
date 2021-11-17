// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-

template<typename T>
T Plane::evalImpl(const T &x, const T &y, const T &z) const {
	return n[0] * (x - p[0]) + n[1] * (y - p[1]) + n[2] * (z - p[2]);
}
