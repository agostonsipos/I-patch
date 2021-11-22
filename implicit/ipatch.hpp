// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-

template<typename T>
T Ipatch::evalImpl(const T &x, const T &y, const T &z) const {
	using std::abs, std::pow;
	T val = 0.0, norm = 0.0;
	for (size_t i = 0; i < N(); ++i) {
		T part = w[i];
		for (size_t j = 0; j < N(); ++j)
			if (i != j)
				part *= pow(B[j]->eval(x, y, z), exponent);
		val += R[i]->eval(x, y, z) * part;
		norm += part;
	}
	T part = w[N()];
	for (size_t j = 0; j < N(); ++j)
		part = part * pow(B[j]->eval(x, y, z), exponent);
	val += part;
	val /= norm;
	return val;
}
