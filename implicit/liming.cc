// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#include "liming.hh"

#include "nelder-mead.hh"

#include <cmath>
#include <string>

using namespace Geometry;

namespace I_patch
{

bool Liming::normalize = false;

Liming::Liming(std::array<Plane, 3> P, const Point3D &point) : planes(std::move(P))
{
	lambda = planes[0](point) * planes[1](point) / pow(planes[2](point), 2);
}

double Liming::fit(const PointVector &approx, bool algebraicError, ApproxError errType)
{
	const auto &P = planes;
	auto midPoint = approx[4];
	std::vector<double> lambdas = { P[0](midPoint) * P[1](midPoint) / pow(P[2](midPoint), 2) };
	auto costFunc = [this, approx, algebraicError, errType](std::vector<double> args) {
		lambda = args[0];
		if(lambda < 0) return std::numeric_limits<double>::infinity();
		double x = 0.0, m = 0.0;
		for (auto p : approx)
		{
			double z;
			if(algebraicError)
			{
				z = pow(evaluateFunction(p), 2);
			}
			else
			{
				Point3D p0;
				try{
					p0 = projectToSurf(p, (planes[0].p - planes[1].p).norm() / 1000);
				} catch(std::runtime_error&) {
					return std::numeric_limits<double>::infinity();
				}
				z = (p - p0).normSqr();
			}
			x += z;
			if (z > m) m = z;
		}
		switch (errType) {
		case ApproxError::Sum:
			return x;
		case ApproxError::Max:
			return m;
		default:
			throw std::runtime_error("unknown approximation error type");
		}
	};
	/*bool success = */NelderMead::optimize(costFunc, lambdas, 50, 1e-10, 1/10.0);
	/*std::cerr << (success
				  ? std::string("Iteration successful")
				  : std::string("Reached max iterations"))
			  << std::endl;*/
	lambda = lambdas[0];

	double d = 0.1 * (approx.front() - approx.back()).norm();
	double val0 = d / evaluateFunction(planes[0].p + planes[0].n * d);
	double val1 = d / evaluateFunction(planes[1].p + planes[1].n * d);

	coeff *= -(val0 + val1) / 2;

	return costFunc(lambdas);
}

double Liming::eval(double x, double y, double z) const {
	return evalImpl<double>(x, y, z);
}

Liming::Dual Liming::eval(const Dual& x, const Dual& y, const Dual& z) const {
	return evalImpl<Dual>(x, y, z);
}

Liming::Dual2 Liming::eval(const Dual2& x, const Dual2& y, const Dual2& z) const {
	return evalImpl<Dual2>(x, y, z);
}

void Liming::printLog()
{
	std::cerr << "\tp1: " << planes[0].p << "\tn1: " << planes[0].n << "\n"
		<< "\tp2: " << planes[1].p << "\tn2: " << planes[1].n << "\n"
		<< "\tlambda: " << lambda << std::endl;
}

}