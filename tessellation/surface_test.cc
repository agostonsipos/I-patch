#include "surface.hh"

#include "plane.hh"
#include "liming.hh"

#define CATCH_CONFIG_MAIN
#include <../misc/catch.hpp>

using namespace I_patch;
using namespace Geometry;
using namespace std;

TEST_CASE("Surface test 1", "[surface]")
{
	Point3D p1(1., 0., 0.), p2(0., 1., 0.), p3(0., 0., 1.);

	RationalBezierCurve c1({ p1, p1 + p2, p2 }, { 1,1 / sqrt(2),1 }),
		c2({ p2, p2 + p3, p3 }, { 1,1 / sqrt(2),1 }),
		c3({ p3, p3 + p1, p1 }, { 1,1 / sqrt(2),1 });

	C0Coons coons({
		make_shared<RationalBezierCurve>(c1),
		make_shared<RationalBezierCurve>(c2),
		make_shared<RationalBezierCurve>(c3)
	});

	Vector3D n1 = p1, n2 = p2, n3 = p3;
	Plane pl1(p1, n1), pl2(p2, n2), pl3(p3, n3);

	Liming l1({ pl1, pl2, Plane(p1, (p1 + p2).normalize()) }, c1.eval(0.5)); l1.invert();
	Liming l2({ pl2, pl3, Plane(p2, (p2 + p3).normalize()) }, c2.eval(0.5)); l2.invert();
	Liming l3({ pl3, pl1, Plane(p3, (p3 + p1).normalize()) }, c3.eval(0.5)); l3.invert();

	Plane b1(p1, n3), b2(p2, n1), b3(p3, n2);

	Ipatch patch(
		{ make_shared<Liming>(l1), make_shared<Liming>(l2), make_shared<Liming>(l3) }, 
		{ make_shared<Plane>(b1), make_shared<Plane>(b2), make_shared<Plane>(b3) }
	);

	patch.setInitialCoeffs(Vector3D(1, 1, 1).normalize());

	Geometry::PointVector sphere_points{ {1,1,1},{2,1,-3},{-1,2,5},{3,0,3} };
	for (auto& it : sphere_points) it = it.normalize();
	
	for (const auto& it : sphere_points)
		REQUIRE(patch(it) < 1e-3);

	Surface surf(coons, patch, { {p1,n1},{p2,n2},{p3,n3} });

	TriMesh result = surf.eval(30, 0.0);

	for (auto it : result.points())
	{
		REQUIRE(((it[0] >= -1e-5) && (it[1] >= -1e-5) && (it[2] >= -1e-5)));
		double val = surf.implicitEval(it);
		if (!std::isnan(val))
			REQUIRE(std::abs(val) < 1e-3);
	}
}

TEST_CASE("Surface test 2", "[surface]")
{
	TriMesh sphere = TriMesh::readOBJ("sphere.obj");

	Point3D p1(1., 0.1, 0.1), p2(0.1, 1., 0.1), p3(0.1, 0.1, 1.);
	p1.normalize(); p2.normalize(); p3.normalize();
	PointVector p {Point3D(p1), Point3D(p2), Point3D(p3)};
	VectorVector n {Vector3D(p1), Vector3D(p2), Vector3D(p3)};

	Surface surf(sphere, p, n, {}, {}, {}, OptimizeType::Algebraic, ApproxError::Sum, RibbonType::Better, BoundingType::ByNormal);

	TriMesh result = surf.eval(30, 0.0);

	for(auto it : result.points())
	{
		REQUIRE(((it[0] >= -1e-5) && (it[1] >= -1e-5) && (it[2] >= -1e-5)));
		double val = surf.implicitEval(it);
		if(!std::isnan(val))
			REQUIRE(std::abs(val) < 1e-3);
	}

	result.writeOBJ("patch.obj");

	surf.setBBox({{Point3D(-1,-1,-1), Point3D(2,2,2)}});

	TriMesh resultMC = surf.evalMC(30, 0.0);

	resultMC.writeOBJ("patchMC.obj");

	TriMesh resultDC = surf.evalDC(30, 0.0);

	resultDC.writeOBJ("patchDC.obj");

	TriMesh ribbons = surf.evalRibbons(30, 0.5);

	ribbons.writeOBJ("ribbons.obj");

	TriMesh ribbon1 = surf.evalRibbonMC(0, 30, 0.0);

	ribbon1.writeOBJ("ribbon1.obj");

	TriMesh ribbon2 = surf.evalRibbonMC(1, 30, 0.0);

	ribbon2.writeOBJ("ribbon2.obj");

	TriMesh ribbon3 = surf.evalRibbonMC(2, 30, 0.0);

	ribbon3.writeOBJ("ribbon3.obj");

	surf.exportColoredMesh("colored1.ply");
}

TEST_CASE( "Surface test 3", "[surface]" )
{
	TriMesh cagd86 = TriMesh::readOBJ("cagd86.obj");

	Point3D p0(-114.506,51.3331,13.558);
	Point3D p1(-65.0471,51.0248,105.198);
	Point3D p2(-14.0368,28.5036,108.462);
	Point3D p3(-10.1002,4.08701,44.0614);
	Point3D p4(-100.937,28.9021,2.73092);

	Vector3D n0(-0.790815,-0.581245,0.191745);
	Vector3D n1(-0.55589,-0.659113,0.506514);
	Vector3D n2(-0.15918,-0.892538,0.421945);
	Vector3D n3(-0.0775253,-0.892662,0.444008);
	Vector3D n4(-0.789165,-0.572749,0.221759);

	PointVector p {p0,p1,p2,p3,p4};
	VectorVector n {n0,n1,n2,n3,n4};

	Surface surf(cagd86, p, n, {}, {}, {}, OptimizeType::None, ApproxError::Sum, RibbonType::Better, BoundingType::ByNormal, 0.0);

	TriMesh result = surf.eval(30, 0.0);

	for(auto it : result.points())
	{
		double val = surf.implicitEval(it);
		if(!std::isnan(val))
			REQUIRE(std::abs(val) < 1e-3);
	}

	result.writeOBJ("cagd86_patch.obj");

	surf.exportColoredMesh("colored2.ply");
}
