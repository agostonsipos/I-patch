// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#include "iloft.hh"

#include "nelder-mead.hh"

#include <cmath>
#include <string>

using namespace Geometry;

namespace I_patch
{

bool Iloft::normalize = false;

double Iloft::maxRatio;

Iloft::Iloft(std::array<Plane, 2> S, std::array<Plane, 2> B, const Point3D &midpoint)
	: primaries(std::move(S)), boundings(std::move(B))
{
	double P0 = primaries[0](midpoint);
	double P1 = primaries[1](midpoint);
	double B0_2 = std::pow(boundings[0](midpoint), 2);
	double B1_2 = std::pow(boundings[1](midpoint), 2);
	double c0 = 1 / B1_2;
	double c1 = 1 / B0_2;
	double coeff = -(c0 * P0 * B1_2 + c1 * P1 * B0_2) / (B0_2 * B1_2);
	coeffs[0] = c0 / coeff;
	coeffs[1] = c1 / coeff;
	coeffs[2] = 1;
}

double Iloft::fit(const PointVector &approx, bool algebraicError, ApproxError errType)
{
	auto sgn = [](double x) {return (x > 0) - (x < 0); };
	std::array<double, 2> startCoeffs = { coeffs[0], coeffs[1] };
	auto costFunc = [this, approx, startCoeffs, algebraicError, errType, sgn](const std::vector<double>& args) {
		if (args[0] * args[1] <= 0) return std::numeric_limits<double>::infinity();
		if (std::max(args[0] / args[1], args[1] / args[0]) > maxRatio && maxRatio != 0) return std::numeric_limits<double>::infinity();
		double w1 = args[0] * startCoeffs[0], w2 = args[1] * startCoeffs[1];
		coeffs = std::array<double, 3>{w1, w2, 1};
		double x = 0.0, m = 0.0;
		for (auto p : approx)
		{
			Point3D p0;
			double z;
			if(algebraicError)
			{
				z = pow(evaluateFunction(p), 2);
			}
			else
			{
				try{
					p0 = projectToSurf(p, (primaries[0].p - primaries[1].p).norm() / 1000);
				} catch(std::runtime_error&) {
					return std::numeric_limits<double>::infinity();
				}
				z = (p - p0).normSqr();
				if ((p - approx[4]).norm() < 1e-10 && args[0] == args[1])
					z = 0;
			}
			x += z;
			if (z > m) m = z;
		}
		if(std::isnan(x)) return std::numeric_limits<double>::infinity();
		switch (errType) {
		case ApproxError::Sum:
			return x;
		case ApproxError::Max:
			return m;
		default:
			throw std::runtime_error("unknown approximation error type");
		}
	};
	std::vector<double> args{ 1, 1 };
	/*bool success = */NelderMead::optimize(costFunc, args, 50, 1e-10, (args[0] + args[1])/10);
	/*std::cerr << (success 
				  ? std::string("Iteration successful")
				  : std::string("Reached max iterations")) 
			  << std::endl;*/
	coeffs = {args[0] * coeffs[0], args[1] * coeffs[1], 1};
	return costFunc(args);
}

double Iloft::eval(double x, double y, double z) const {
	return evalImpl<double>(x, y, z);
}

Iloft::Dual Iloft::eval(const Dual &x, const Dual &y, const Dual &z) const {
	return evalImpl<Dual>(x, y, z);
}

void Iloft::printLog()
{
	std::cerr << "\tp1: " << primaries[0].p << "\tn1: " << primaries[0].n << "\n"
		<< "\tp2: " << primaries[1].p << "\tn2: " << primaries[1].n << "\n"
		<< "\tcoeffs: " << coeffs[0] << " " << coeffs[1] << " " << coeffs[2] << std::endl;
}

}