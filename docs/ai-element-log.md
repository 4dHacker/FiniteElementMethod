# AI local-element implementation record

Date: 2026-11-11. Working branch: ai/gmsh-reader (uncommitted local milestones).

User request: implement the local element calculation, then provide a detailed
English hand-verification workflow with given, computed and unknown quantities.
Global assembly and solving are not part of this change.

Implementation:

- Added include/fem/element.hh for K=1 P1 triangles, using the supplied handout's
  reference map, inverse-transpose gradient transformation and three interior
  quadrature points with reference weights 1/6.
- Preserved input vertex order and used abs(det(B)) for integration, so clockwise
  and counter-clockwise cells give matching results after local permutation.
- Used std::array for tiny fixed-size local data; HDNUM remains the global
  sparse linear-algebra dependency. No new numerical library is required.
- Added finite-value and relative-degeneracy checks with explicit arithmetic
  range limitations.
- Added element_demo to print all intermediate quantities for two fixed cases,
  plus an English worksheet and separate exact reference answers.
- Added focused CTest checks and a Windows build script; retained reader tests.

Verification covers exact local matrix/load values, reference quadrature
polynomial moments, a sheared element with affine source, all six vertex
permutations, rotation/scaling/translation, affine-field energy, zero source,
and invalid/non-finite/unrepresentable inputs.

All automated checks and reference answers in this milestone were produced by
AI. Human hand verification is pending and must not be described as completed.
Retain the user's independent calculations and discrepancies as separate evidence.
The existing assembly placeholder and global solver behavior were not changed.
