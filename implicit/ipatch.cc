// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#include "ipatch.hh"

#include "nelder-mead.hh"

#include <cmath>
#include <functional>
#include <numeric>
#include <algorithm>

using namespace Geometry;

namespace I_patch
{

const int Ipatch::exponent = 2;

double Ipatch::maxRatio;

void Ipatch::setInitialCoeffs(const Point3D &center, bool use_p_i)
{
	setInitialCoeffs(center, std::vector<double>(N() + 1, 1.0), use_p_i);
}

void Ipatch::setInitialCoeffs(const Point3D &center, const std::vector<double> &d, bool use_p_i)
{
	_d = d;
	w.resize(N() + 1);
	double lastCoeff = 0;

	for (size_t i = 0; i < N(); ++i)
	{
		double pb2 = (*(R[i]))(center) / std::abs(pow((*(B[i]))(center), exponent));
		double b2 = 1 / std::abs(pow((*(B[i]))(center), exponent));
		if(use_p_i){
			w[i] = std::abs(d[i] / pb2);
			lastCoeff -= d[i] * ((pb2 >= 0) ? 1 : -1);
		}
		else{
			w[i] = std::abs(d[i] / b2);
			lastCoeff -= d[i] * ((b2 >= 0) ? 1 : -1) * (*(R[i]))(center);
		}
		//lastCoeff -= d[i] * ((pb2 >= 0) ? 1 : -1);
	}
	if (d.size() > N())
		lastCoeff *= d[N()];
	w[N()] = lastCoeff;
	this->center = center;
}

void Ipatch::setOptimizedCoeffs(const PointVector &approx, double step, bool algebraic)
{
	auto costFunc = [this, &approx, step, algebraic](const std::vector<double> &args) {
		setInitialCoeffs(center, args);

		double x = (*std::max_element(args.begin(), args.end())) / (*std::min_element(args.begin(), args.end()));
		double pen = std::max(x, 1 / x);
		if ((pen > maxRatio && maxRatio != 0.0) || pen < 0) return std::numeric_limits<double>::infinity();

		/*std::cerr << "Values in iteration: ";
		for (auto it : args)
			std::cerr << it << " ";
		std::cerr << std::endl;*/

		double sumDistance = 0;
		double max = 0;

		for (auto p : approx)
		{
			double x;
			if(algebraic)
			{
				x = pow(evaluateFunction(p), 2);
			}
			else{
				try{
					x = (p - projectToSurf(p, step)).normSqr();
				} catch(std::runtime_error &) {
					return std::numeric_limits<double>::infinity();
				}
			}
			sumDistance += x;
			if (x > max)
				max = x;
		}
		//std::cerr << "sum of distance: " << sumDistance << std::endl;
		//std::cerr << "max distance: " << max << std::endl;
		//std::cerr << std::endl;
		return sumDistance;
	};

	std::vector<double> args(_d);

	std::cerr << "Initial sumDistance: " << costFunc(args) << std::endl;

	/*bool success = */NelderMead::optimize(costFunc, args, 500, 1e-4, 0.5);

	//std::cerr << (success ? std::string("Iteration successful") : std::string("Reached max iterations")) << std::endl;

	setInitialCoeffs(center, args);

	std::cerr << "Optimized coeffs: " << std::endl;
	for (size_t i = 0; i <= N(); ++i)
	{
		std::cerr << weight(i) << " ";
	}
	std::cerr << std::endl;

	std::cerr << "Final sumDistance: " << costFunc(args) << std::endl;
}

double Ipatch::eval(double x, double y, double z) const {
	return evalImpl<double>(x, y, z);
}

Ipatch::Dual Ipatch::eval(const Dual& x, const Dual& y, const Dual& z) const {
	return evalImpl<Dual>(x, y, z);
}

Ipatch::Dual2 Ipatch::eval(const Dual2& x, const Dual2& y, const Dual2& z) const {
	return evalImpl<Dual2>(x, y, z);
}

}