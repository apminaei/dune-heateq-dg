// Gmsh project created on Fri Sep 11 15:40:58 2020
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1000, 0, 0, 1.0};
//+
Point(3) = {1000, -800, 0, 1.0};
//+
Point(4) = {0, -800, 0, 1.0};
//+
Point(5) = {0, -500, 0, 1.0};
//+
Point(6) = {1000, -500, 0, 1.0};
//+
Point(7) = {1000, -400, 0, 1.0};
//+
Point(8) = {0, -400, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 7};
//+
Line(3) = {7, 6};
//+
Line(4) = {6, 3};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 5};
//+
Line(7) = {5, 8};
//+
Line(8) = {8, 1};
//+
Line(9) = {1, 8};
//+
Line(10) = {7, 8};
//+
Line(11) = {5, 6};
//+
Physical Curve("isBaseGHSZ") = {11};
//+
Physical Curve("isTopGHSZ") = {10};
//+
Physical Curve("isTopBoundary") = {1};
//+
Physical Curve("isBottomBoundary") = {5};
//+
Physical Curve("isLeftBoundary") = {8, 7, 6};
//+
Physical Curve("isRightBoundary") = {2, 3, 4};
//+
Transfinite Curve {6, 4} = 10 Using Progression 1;
//+
Transfinite Curve {3, 7} = 15 Using Progression 1;
//+
Transfinite Curve {8, 2} = 20 Using Progression 1;
//+
Transfinite Curve {1, 10, 11, 5} = 15 Using Progression 1;

//+
Curve Loop(1) = {1, 2, 10, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, -7, 11, -3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 4, 5, 6};
//+
Plane Surface(3) = {3};
//+
Transfinite Surface {1} = {1, 2, 7, 8};
//+
Transfinite Surface {2} = {8, 7, 6, 5};
//+
Transfinite Surface {3} = {5, 6, 3, 4};
//+
Recombine Surface {1, 2, 3};
//+
Recombine Surface {1, 2, 3};
