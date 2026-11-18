// Unit disk centered at the origin.
h = 0.2;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {0, 1, 0, h};
Point(4) = {-1, 0, 0, h};
Point(5) = {0, -1, 0, h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Curve("boundary") = {1, 2, 3, 4};
Physical Surface("domain") = {1};
