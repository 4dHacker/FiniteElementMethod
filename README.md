# fem2d

2D finite element method (P1 triangles) using [HDNUM](https://parcomp-git.iwr.uni-heidelberg.de/Teaching/hdnum), Gmsh and ParaView.

Pipeline: `Gmsh .geo → .msh → Mesh{points, triangles} → assemble sparse matrix (HDNUM) → solve → .vtk → ParaView`

## Setup

```bash
git clone --recursive <repo-url>      # or: git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Run

```bash
gmsh -2 -format msh2 meshes/square.geo -o meshes/square.msh   # ASCII format 2 is required
./build/fem2d meshes/square.msh results/solution.vtk
paraview results/solution.vtk
```

## Layout

| Path | Purpose |
|---|---|
| `include/fem/mesh.hh` | `Mesh`: `points` (N×2), `triangles` (M×3, 0-based), `on_boundary` |
| `include/fem/mesh_reader.hh` | Gmsh msh2 → `Mesh` |
| `include/fem/problem.hh` | PDE data (`f`, `g`); edit for your equation |
| `include/fem/assembly.hh` | local → global assembly into HDNUM `SparseMatrix`, Dirichlet BCs |
| `include/fem/vtk_writer.hh` | solution → legacy VTK |
| `src/main.cc` | driver |
| `tests/` | small checks with known answers |
| `meshes/` | `.geo` sources (commit) and small `.msh` files |
| `docs/` | math derivation |
| `external/hdnum` | HDNUM git submodule |

## Team conventions

- Work on feature branches, merge via pull requests; keep `main` building.
- Commit `.geo` files, not large generated meshes or results.
- Agree on interface changes (`Mesh`, `assemble_*` signatures) before changing them.
