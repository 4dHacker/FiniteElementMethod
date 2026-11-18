# fem2d

2D finite element method (P1 triangles) using [HDNUM](https://parcomp-git.iwr.uni-heidelberg.de/Teaching/hdnum), Gmsh and ParaView.

Pipeline: `Gmsh .geo → .msh → Mesh{points, triangles} → assemble sparse matrix (HDNUM) → solve → .vtk → ParaView`

The driver assembles and solves the Poisson problem with P1 triangles and
homogeneous Dirichlet boundary conditions. The right-hand side is selected in
`src/main.cc`; see the example problems for matching functions and meshes.

## Setup

```bash
git clone --recursive <repo-url>      # or: git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Run

```bash
gmsh -2 meshes/square.geo -o meshes/square.msh
./build/fem2d meshes/square.msh results/solution.vtk
paraview results/solution.vtk
```

The mesh reader accepts ASCII Gmsh MSH 2.x and 4.x files.
See [example problems](docs/example-problems.md) for meshes and matching
right-hand-side lambdas.

## Layout

| Path | Purpose |
|---|---|
| `include/fem/mesh.hh` | `Mesh`: `points` (N×2), `triangles` (M×3, 0-based), `on_boundary` |
| `include/fem/mesh_reader.hh` | Gmsh MSH 2.x/4.x → `Mesh` |
| `include/fem/vtk_writer.hh` | solution → legacy VTK |
| `src/main.cc` | driver, FEM assembly, boundary conditions and solve |
| `tests/` | small checks with known answers |
| `meshes/` | `.geo` sources (commit) and small `.msh` files |
| `docs/` | math derivation |
| `external/hdnum` | HDNUM git submodule |

## Team conventions

- Work on feature branches, merge via pull requests; keep `main` building.
- Commit `.geo` files, not large generated meshes or results.
- Agree on interface changes (`Mesh`, `assemble_*` signatures) before changing them.
