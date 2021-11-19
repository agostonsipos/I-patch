cmake_minimum_required(VERSION 3.9)

project(I_patch)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif(CMAKE_COMPILER_IS_GNUCXX)


set(LIBGEOM_DIR "libgeom") # https://github.com/salvipeter/libgeom
include_directories("${LIBGEOM_DIR}")

set(AD_DIR "autodiff") # https://github.com/autodiff/autodiff
include_directories("${AD_DIR}")

add_subdirectory(libgeom)
add_subdirectory(c0coons)
add_subdirectory(dfo)
add_subdirectory(dual-contouring)
add_subdirectory(implicit_marching_cubes)
add_subdirectory(implicit)
add_subdirectory(tessellation)

file(COPY misc/assets/sphere.obj DESTINATION ${CMAKE_BINARY_DIR})
file(COPY misc/assets/cagd86.obj DESTINATION ${CMAKE_BINARY_DIR})
file(COPY misc/assets/tests DESTINATION ${CMAKE_BINARY_DIR})
