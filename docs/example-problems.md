Erstellt mit KI
# Example problems

The program solves `-Δu = f` with zero Dirichlet values on the mesh boundary.
For each example below, replace the `f` lambda in `src/main.cc`, rebuild, and
run the matching mesh. The listed exact solution is useful for checking the
shape and sign of the computed result in ParaView.

Generate the disk and L-shaped meshes with Gmsh (the default output is MSH 4.x,
which the reader supports):

```bash
gmsh -2 meshes/disk.geo -o meshes/disk.msh
gmsh -2 meshes/l_shape.geo -o meshes/l_shape.msh
```

## Unit square

Mesh: `meshes/square.msh` (generated from `meshes/square.geo`).

Exact solution: `u(x,y) = sin(πx) sin(πy)`.

```cpp
const double pi = std::acos(-1.0);
auto f = [pi](double x, double y) {
  return 2.0 * pi * pi * std::sin(pi * x) * std::sin(pi * y);
};
```

## Unit disk

Mesh: `meshes/disk.msh`.

Exact solution: `u(x,y) = 1 - x² - y²`.

```cpp
auto f = [](double, double) { return 4.0; };
```

## L-shaped domain

Mesh: `meshes/l_shape.msh`. The domain is `[0,1]²` with
`[0,0.5] × [0.5,1]` removed.

Define `p(t) = t(1-t)(t-0.5)` and use the exact solution
`u(x,y) = p(x)p(y)`. This vanishes on every boundary segment.

```cpp
const auto p = [](double t) { return t * (1.0 - t) * (t - 0.5); };
const auto p_second = [](double t) { return 3.0 - 6.0 * t; };
auto f = [p, p_second](double x, double y) {
  return -(p_second(x) * p(y) + p(x) * p_second(y));
};
```

Run any case and view its solution with ParaView:

```bash
cmake --build build --target fem2d
./build/fem2d meshes/disk.msh results/disk.vtk
paraview results/disk.vtk
```
