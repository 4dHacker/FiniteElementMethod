# fem2d

2D finite element method (P1 triangles) using [HDNUM](https://parcomp-git.iwr.uni-heidelberg.de/Teaching/hdnum), Gmsh and ParaView.

Pipeline: `Gmsh .geo → .msh → Mesh{points, triangles} → assemble sparse matrix (HDNUM) → solve → .vtk → ParaView`

## AI branch: reader milestone

This branch stops at the 2D mesh reader; assembly is still a placeholder.
See [AI reader scope, provenance and validation](docs/ai-gmsh-reader.md).
On the configured Windows machine, run `./scripts/build-reader.ps1` in PowerShell.
With a configured CMake build, use `cmake --build build --target read_mesh test_mesh_reader`,
then `ctest --test-dir build -R '^mesh_reader$' --output-on-failure` and
`./build/read_mesh meshes/square.msh` (Windows: `read_mesh.exe`).

## AI local-element milestone

Local P1 Poisson element calculations are available independently of assembly.
See [the hand-verification worksheet](docs/element-hand-check.md) and the separate
[reference answers](docs/element-reference-results.md).
On the configured Windows machine, run `./scripts/build-element.ps1` in PowerShell.
After doing the hand calculation, run `./build/element_demo.exe reference` or
`./build/element_demo.exe skew` with the MinGW runtime on PATH.
The full assembly/solver remains unfinished.

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
## Comparing the implementations

The AI track on this branch and the team implementation on main started at
ba6f294. Compare each branch against that common commit for attribution.
On 2026-11-22 the team snapshot was 5007728. Its commit messages identify AI
assistance on some tasks. Record assistance per task when reporting the comparison.

Verified on Windows: 39 reader checks and 221 element checks passed. The shared
square mesh yielded 142 points, 242 triangles and 40 boundary points. The AI
milestone covers the ASCII MSH 2.2 reader and local P1 elements. Global assembly
is still a placeholder on this branch. Human hand verification is pending.
