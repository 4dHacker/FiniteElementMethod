# AI group: Gmsh reader

## Provenance and scope (2026-11-10)

This branch is the AI-assisted contribution, separate from the human group.
Starting commit: `ba6f2946c629e1b981387dca7ea68f7a120fdd26` on `origin/main`.
Local branch: `ai/gmsh-reader`. No commits or changes have been pushed.

The shared repository already supplied the Mesh interface and a basic reader.
The AI adapted its previous standalone reader (written during the earlier project
work) to this interface, added validation, reader-only execution, and tests.
This is not a from-scratch implementation starting on this date; earlier reader
development and environment setup are part of the AI work history.

The user requested migration to this repository and work only up to the 2D
triangle reader, identifying themselves as the AI group and the two teammates
as the human group. AI generated the code, tests, setup scripts, and this record.
The conversation supplies the prompt/tool history; this document is a concise
record, not a complete transcript or a measurement of human development time.

## Shared interface

Existing `fem::Mesh` and `read_gmsh_msh2(const std::string&)` are preserved.
An input-stream overload supports tests without temporary files.
The reader returns points, zero-based triangle connectivity, and an on_boundary
flag per point. Like the original implementation, type-1 line endpoints set the
flag; internal tagged curves would also be flagged. Boundary topology is not
inferred from triangles. Physical names and tags are not returned because the
team's Mesh has no corresponding fields. Supporting distinct boundary-condition
groups later requires an agreed interface extension.

Supported input: MSH 2.2 ASCII, first-order triangles (type 2), line segments
(type 1), and optional point elements (type 15, validated then discarded).
Coordinates must be finite and z=0. Sparse/unsorted IDs and either triangle
orientation work. Unknown data sections are skipped. Binary, MSH 4.x, quads,
higher-order elements and 3D elements are rejected explicitly.

Compared with the initial reader, invalid counts, truncated sections, duplicate
IDs, unknown node references, zero-area triangles and zero-length lines now
produce diagnostics. This intentionally rejects unsupported elements instead of
silently returning a partial mesh. It does not check complete mesh topology,
overlapping cells, or mesh quality.

## Run on this Windows machine

Open the new repository folder, then in PowerShell:

```powershell
.\scripts\build-reader.ps1
```

This configures the project's C++20 CMake build with the existing MSYS2 GCC
15.1.0, builds only reader targets, runs reader tests, then reads the committed
square mesh. It uses the project's pinned HDNUM submodule, not the old checkout.
No new software or global PATH settings were installed.

For an already configured build (ensure MinGW runtime DLLs are on PATH):

```powershell
$env:PATH = "C:\msys64\mingw64\bin;$env:PATH"
.\build\read_mesh.exe .\meshes\square.msh
ctest --test-dir build -R '^mesh_reader$' --output-on-failure
```

On Linux/macOS, use the existing README's CMake setup and the `read_mesh` and
`test_mesh_reader` targets. The Windows-only script is optional. These platforms
have not been tested in this session.

VS Code: open this repository and use **Terminal > Run Build Task** for the
reader build/test/demo. Compiling a single file with Code Runner does not supply
the project's include directories or required program arguments.

## Validation performed

- All project targets compile with GCC 15.1.0 in C++20 Release mode.
- Reader validation: 39 checks passed. Includes sparse IDs, orientation, boundary
  flags with an interior point, total area, optional tags/names, point elements,
  CRLF, additional sections, and malformed/unsupported inputs.
- Committed `meshes/square.msh`: 142 points, 242 triangles, 40 boundary points.
  This validates against the team's existing mesh; Gmsh was not run in this session.
- Full CTest run: reader passes; the existing assembly test segfaults. Its
  `assemble_poisson` currently returns an empty matrix. Neither assembly nor its
  test, PDE data, solver or VTK writer were changed in this contribution.
  Do not describe the full solver or full test suite as working yet.

The test cases were also written by AI and are not an independent assessment.
For the final comparison, the professor/team should supply common held-out meshes
and acceptance cases under the same agreed requirements.
