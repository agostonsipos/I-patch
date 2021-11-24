# I-patch
Library for creating, evaluating and tessellating I-patches

Based on papers:
* Várady, T., Benkő, P., Kós, G., & Rockwood, A. (2001). Implicit surfaces revisited—I-patches. In Geometric Modelling (pp. 323-335). Springer, Vienna. [doi.org/10.1007/978-3-7091-6270-5_19](https://doi.org/10.1007/978-3-7091-6270-5_19)
* Sipos, Á., Várady, T., Salvi, P., & Vaitkus, M. (2020). Multi-sided implicit surfacing with I-patches. Computers & Graphics, Vol. 90, 29-42. [doi.org/10.1016/j.cag.2020.05.009](https://doi.org/10.1016/j.cag.2020.05.009)
* Sipos, Á., Várady, T., & Salvi, P. (2022). Approximating Triangular Meshes by Implicit, Multi-Sided Surfaces. Computer-Aided Design and Applications, Vol. 19 (to be published)

## Introduction

I-patches are a family of multi-sided, implicit patches. They are defined by a number of side interpolant (also called ribbon) surfaces and the same number of bounding surfaces. The resulting surface is defined inside the part of space delimited by the bounding surfaces, and interpolates the ribbon surfaces with a given order of geometric continuity.

![ribbons & boundings](https://user-images.githubusercontent.com/25045084/143272235-4164b982-9029-490c-b4dc-7603b03b0bc0.png) ![patch](https://user-images.githubusercontent.com/25045084/143272242-a8bf97b9-2fa0-426f-afcc-ab1b4553b56c.png)

![example2](https://user-images.githubusercontent.com/25045084/143278556-9b6e08b5-25d4-4876-a7e0-7267d8808254.png)

Implicit surfaces are efficient for a handful of operations, including inside/outside classification, distance calculation or line intersection (raycasting). I-patches seek to overcome their weaknesses in producing shape problems (like singularities, self-intersections, and disconnected parts), by limiting the surface to a portion of the whole space, while retaining general topology.

![polyhedral](http://3dgeo.iit.bme.hu/images/implicit/polyhedral2-shaded.png)

## Compiling

All dependencies are included as git submodules except [Eigen](https://gitlab.com/libeigen/eigen). Building the code should be straightforward with CMake, e.g.:
```
git submodule update --init --recursive
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make
```

On Windows, you may use CMake to generate project files for the IDE of your choice.

## Usage

### As implicit functions

`class Ipatch`, defined in `implicit/ipatch.hh` represents an implicit function (R<sup>3</sup>→R), whose 0-isosurface is of our interest. Usage examples (creating, evaluating) can be seen in `implicit_test.cc`

### As multi-sided patches
`class Surface`, defined in `tessellation/surface.hh` represents a finite portion of the implicit surface, enclosed by a boundary curve loop. (The curve loop should consist of the intersection curves of the ribbon and bounding surfaces.) A `Surface` object can be created two ways.
* From an existing implicit function and curveloop
* By supplying a mesh and a number of corner points on it; the surface is automatically approximated

A `Surface` object is capable of creating a triangle mesh of the area inside the curveloop (see Sipos, Á., & Salvi, P. (2021). [Creating good quality meshes from smooth implicit surfaces](http://3dgeo.iit.bme.hu/papers/implicit/meshing.pdf)). However, it can also create the isosurface meshes of the whole implicit surface, inside a bounding box, with either Marching Cubes, or Dual Contouring.

Examples are in `surface_test.cc`.