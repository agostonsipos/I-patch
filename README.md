# I-patch
Library for creating, evaluating and tessellating I-patches

Based on papers:
* Várady, T., Benkő, P., Kós, G., & Rockwood, A. (2001). Implicit surfaces revisited—I-patches. In Geometric Modelling (pp. 323-335). Springer, Vienna.
* Sipos, Á., Várady, T., Salvi, P., & Vaitkus, M. (2020). Multi-sided implicit surfacing with I-patches. Computers & Graphics, Vol. 90, 29-42.
* Sipos, Á., Várady, T., & Salvi, P. (2022). Approximating Triangular Meshes by Implicit, Multi-Sided Surfaces. Computer-Aided Design and Applications, Vol. 19

## Compiling

All dependencies are included as git submodules except [Eigen](https://gitlab.com/libeigen/eigen). Building the code should be straightforward with CMake, e.g.:
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make
```

On Windows, you may use CMake to generate project files for the IDE of your choice.

## Usage

### As implicit functions

`class Ipatch`, defined in `implicit/ipatch.hh` represents an implicit function ($$\mathbb{R}^3 \rightarrow \mathbb{R}$$), whose 0-isosurface is of our interest. Usage examples (creating, evaluating) can be seen in `implicit_test.cc`

### As multi-sided patches
`class Surface`, defined in `tessellation/surface.hh` represents a finite portion of the implicit surface, enclosed by a boundary curve loop. A `Surface` object can be created two ways.
* From an existing implicit function and curveloop
* By supplying a mesh and a number of corner points on it; the surface is automatically approximated

Examples are in `surface_test.cc`.