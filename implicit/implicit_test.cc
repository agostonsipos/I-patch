// -*- c-basic-offset: 4; indent-tabs-mode: t; tab-width: 4 -*-
#include "plane.hh"
#include "liming.hh"
#include "iloft.hh"
#include "ipatch.hh"

#include <cmath>
#define CATCH_CONFIG_MAIN
#include <../misc/catch.hpp>

using namespace I_patch;
using namespace Geometry;
using namespace std;

TEST_CASE( "Iloft test", "[iloft]" )
{
	Plane p1 (Point3D(0,0,0), Vector3D(-1,1,0)), p2 (Point3D(1,0,0), Vector3D(1,1,0));
	Plane b1 (Point3D(0,0,0), Vector3D(-1,0,0)), b2 (Point3D(1,0,0), Vector3D(1,0,0));

	Point3D mid (0.5, 0.25, 0);
	Iloft il ({p1, p2}, {b1, b2}, mid);

	REQUIRE(il(mid) == 0.0);
	REQUIRE(il(Point3D(0.5,0.3,0)) > 0.0);
	REQUIRE(il(Point3D(0.5,0.2,0)) < 0.0);

	REQUIRE(il(Point3D(0,2,10)) == p1(Point3D(0,2,10)));
	REQUIRE(il(Point3D(1,-2,7)) == p2(Point3D(1,-2,7)));
	REQUIRE(il(Point3D(0,0,10)) == 0.0);
}

TEST_CASE( "Ipatch test", "[ipatch]" )
{
	double sqh = 1/sqrt(2);
	Plane p11 {Point3D(1,0,0), Vector3D(1,0,0)};
	Plane p12 {Point3D(0,1,0), Vector3D(0,1,0)};
	Plane p13 {Point3D(1,0,0), Vector3D(sqh, sqh, 0)};
	Liming l1({{p11, p12, p13}}, Point3D(sqh, sqh, 0));
	l1.invert();
	using namespace Catch::literals;
	REQUIRE( l1(Point3D(sqh, sqh, 0)) == Approx(0).margin(1e-10) );  // as it is on the surface
	REQUIRE( l1(Point3D(sqh, -sqh, 0)) == Approx(0).margin(1e-10) ); // as it is a cylinder
	REQUIRE( l1(Point3D(0,0,0)) < 0 );                               // inside should be negative

	
	Plane p21 {Point3D(0,1,0), Vector3D(0,1,0)};
	Plane p22 {Point3D(0,0,1), Vector3D(0,0,1)};
	Plane p23 {Point3D(0,1,0), Vector3D(0, sqh, sqh)};
	Liming l2({{p21, p22, p23}}, Point3D(0, sqh, sqh));
	l2.invert();
	REQUIRE( l2(Point3D(0,0,0)) < 0 );                               // inside should be negative
	
	
	Plane p31 {Point3D(0,0,1), Vector3D(0,0,1)};
	Plane p32 {Point3D(1,0,0), Vector3D(1,0,0)};
	Plane p33 {Point3D(0,0,1), Vector3D(sqh, 0, sqh)};
	Liming l3({{p31, p32, p33}}, Point3D(sqh, 0, sqh));
	l3.invert();
	REQUIRE( l3(Point3D(0,0,0)) < 0 );                               // inside should be negative
	
	Plane b1 {Point3D(0,0,0), Vector3D(0,0,1)};
	Plane b2 {Point3D(0,0,0), Vector3D(1,0,0)};
	Plane b3 {Point3D(0,0,0), Vector3D(0,1,0)};
	
	Ipatch I({make_shared<Liming>(l1), make_shared<Liming>(l2), make_shared<Liming>(l3)},
			{make_shared<Plane>(b1), make_shared<Plane>(b2), make_shared<Plane>(b3)});
	
	double sqt = 1/sqrt(3);
	I.setInitialCoeffs(Point3D(sqt, sqt, sqt));
	
	for(int i = 0; i < 3; ++i)
		REQUIRE(I.weight(i) == 2_a);
	
	REQUIRE( I(Point3D(sqt, sqt, sqt)) == Approx(0).margin(1e-10) ); // as it is a sphere
	auto d = I.grad(Point3D(sqt, sqt, sqt)).normalize();
	REQUIRE( d[0] == Approx(sqt) );
	REQUIRE( d[0] == Approx(d[1]) );
	REQUIRE( d[1] == Approx(d[2]) );
	REQUIRE( d[2] == Approx(d[0]) );
	
	REQUIRE( I(Point3D(-sqt, sqt, -sqt)) == Approx(0).margin(1e-10) ); // as it is a sphere
}
