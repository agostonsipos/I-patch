// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-

template<typename T>
T Liming::evalImpl(const T &x, const T &y, const T &z) const {
	auto cut = planes[2].eval(x, y, z);
	T v = coeff * (planes[0].eval(x, y, z) * planes[1].eval(x, y, z) - lambda * cut * cut);
	//v /= -(planes[0].eval(x,y,z) + planes[1].eval(x,y,z) - 2*sqrt(lambda)*planes[2].eval(x,y,z));
	if (normalize) {
		if constexpr(std::is_same<T, double>::value)
			v /= evaluateGradient({ x, y, z }).norm();
	}
	return v;
}
