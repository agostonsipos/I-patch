cmake_minimum_required(VERSION 3.1.0)

project(marching_cubes)

if(CMAKE_VERSION VERSION_LESS "3.7.0")
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif(CMAKE_COMPILER_IS_GNUCXX)

set(GEOM_DIR "../libgeom") # https://github.com/salvipeter/libgeom
include_directories("${GEOM_DIR}")

set(DFO_DIR "../dfo") # https://github.com/salvipeter/libgeom
include_directories("${DFO_DIR}")

link_directories("${GEOM_DIR}/build")


add_library(implicit STATIC
    implicit.hh
    implicit.cc
    plane.hh
    plane.hpp
    plane.cc
    liming.hh
    liming.hpp
    liming.cc
    iloft.hh
    iloft.hpp
    iloft.cc
    ipatch.hh
    ipatch.hpp
    ipatch.cc
)

add_executable(implicit_test
	implicit_test.cc
)

set_property(TARGET implicit PROPERTY CXX_STANDARD 17)
set_property(TARGET implicit PROPERTY CXX_STANDARD_REQUIRED ON)

set_property(TARGET implicit_test PROPERTY CXX_STANDARD 17)
set_property(TARGET implicit_test PROPERTY CXX_STANDARD_REQUIRED ON)

target_link_libraries(implicit geom)
target_link_libraries(implicit dfo)

target_link_libraries(implicit_test implicit)
