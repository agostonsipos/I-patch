cmake_minimum_required(VERSION 3.1.0)

project(tessellation)

if(CMAKE_VERSION VERSION_LESS "3.7.0")
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif(CMAKE_COMPILER_IS_GNUCXX)

set(GEOM_DIR "../libgeom")
include_directories("${GEOM_DIR}")

set(IMC_DIR "../implicit_marching_cubes")
include_directories("${IMC_DIR}")

set(DC_DIR "../dual-contouring")
include_directories("${DC_DIR}")

set(C0C_DIR "../c0coons")
include_directories("${C0C_DIR}")

include_directories("../implicit")

add_library(tessellation STATIC
    surface.hh
    surface.cc
    piecewise.hh
    piecewise.cc
    mesh_utils.hh
    mesh_utils.cc
)

add_executable(tessellation_test
    surface_test.cc
)

set_property(TARGET tessellation PROPERTY CXX_STANDARD 17)
set_property(TARGET tessellation PROPERTY CXX_STANDARD_REQUIRED ON)

set_property(TARGET tessellation_test PROPERTY CXX_STANDARD 17)
set_property(TARGET tessellation_test PROPERTY CXX_STANDARD_REQUIRED ON)

target_link_libraries(tessellation geom)
target_link_libraries(tessellation dualcontour)
target_link_libraries(tessellation marching_cubes)
target_link_libraries(tessellation c0coons)
target_link_libraries(tessellation implicit)
if(UNIX)
	target_link_libraries(tessellation tbb)
endif()

target_link_libraries(tessellation_test tessellation)
target_link_libraries(tessellation_test dualcontour)
target_link_libraries(tessellation_test marching_cubes)
