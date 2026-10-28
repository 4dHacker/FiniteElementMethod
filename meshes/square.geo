// Unit square. Generate the mesh with:
//   gmsh -2 -format msh2 meshes/square.geo -o meshes/square.msh
h = 0.1;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Curve("boundary") = {1, 2, 3, 4};
Physical Surface("domain") = {1};
