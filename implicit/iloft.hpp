// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-

template<typename T>
T Iloft::evalImpl(const T &x, const T &y, const T &z) const {
	using std::pow;
	T v =
		coeffs[0] * primaries[0].eval(x, y, z) * pow(boundings[1].eval(x, y, z), 2) +
		coeffs[1] * primaries[1].eval(x, y, z) * pow(boundings[0].eval(x, y, z), 2) +
		coeffs[2] * pow(boundings[0].eval(x, y, z) * boundings[1].eval(x, y, z), 2);
	v /= (pow(boundings[0].eval(x, y, z), 2) * coeffs[1] + pow(boundings[1].eval(x, y, z), 2) * coeffs[0]);
	
	if (normalize) {
		if constexpr(std::is_same<T, double>::value)
			v /= evaluateGradient({ x, y, z }).norm();
	}
	return v;
}
