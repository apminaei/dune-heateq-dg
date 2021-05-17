// Gmsh project created on Mon March 2 2020

// characteristic length
char_X = 1000.0;

//######################
// DOMAIN AND STRUCTURES
//######################
// origin
O_x = 0.0/char_X;
O_y = 0.0/char_X;
// Domain X-length, Y-length
L_x = 1000.0/char_X;
L_y = 800.0/char_X;
// Depth of hydrate BSR below sea floor
hyd_bsr = 360.0/char_X;
// free-gas-pocket: location below BSR, depth, width
fgp_offset = 0.0/char_X;
fgp_depth = 100.0/char_X;
fgp_width = 500.0/char_X;
// Fracture diameter and depth below sea floor
frac_diameter = 50.0/char_X;
frac_depth = 300.0/char_X;

//######################
// MESH PARAMETERS
//######################
mesh_param = 0.25/char_X;
A1_y = 1.25;
N1_y = 15;
A3_y = 1;
N3_y = 10;
A2_x = 1;
N2_x = 10;
A4_x = 1.125;
N4_x = 18;
A5_x = 1;
N5_x = 10;

//######################
// MESH
// !! DO NOT TOUCH THIS PART !!
//######################
// full domain
Point(1) = {O_x	, O_y	, 0, mesh_param};
Point(2) = {L_x	, O_y	, 0, mesh_param};
Point(3) = {L_x	, L_y	, 0, mesh_param};
Point(4) = {O_x	, L_y	, 0, mesh_param};

// fracture
frac_O_x = L_x*0.5-frac_diameter/2;
frac_O_y = L_y-frac_depth;
frac_L_x = L_x*0.5+frac_diameter/2;
frac_L_y = L_y;
Point(5) = {frac_O_x	, frac_O_y	, 0, mesh_param};
Point(6) = {frac_L_x	, frac_O_y	, 0, mesh_param};
Point(7) = {frac_L_x	, frac_L_y	, 0, mesh_param};
Point(8) = {frac_O_x	, frac_L_y	, 0, mesh_param};

// free gas pocket
fgp_O_x = L_x*0.5-fgp_width/2;
fgp_O_y = L_y-hyd_bsr-fgp_offset-fgp_depth;
fgp_L_x = L_x*0.5+fgp_width/2;
fgp_L_y = L_y-hyd_bsr-fgp_offset;
Point(9)  = {fgp_O_x	, fgp_O_y	, 0, mesh_param};
Point(10) = {fgp_L_x	, fgp_O_y	, 0, mesh_param};
Point(11) = {fgp_L_x	, fgp_L_y	, 0, mesh_param};
Point(12) = {fgp_O_x	, fgp_L_y	, 0, mesh_param};


Point(13) = {L_x	, fgp_O_y 	, 0, mesh_param};
Point(14) = {O_x	, fgp_O_y 	, 0, mesh_param};
Point(15) = {L_x	, frac_O_y	, 0, mesh_param};
Point(16) = {O_x	, frac_O_y	, 0, mesh_param};
Point(17) = {frac_O_x	, fgp_L_y	, 0, mesh_param};
Point(18) = {frac_L_x	, fgp_L_y	, 0, mesh_param};
Point(19) = {frac_O_x	, fgp_O_y	, 0, mesh_param};
Point(20) = {frac_L_x	, fgp_O_y	, 0, mesh_param};
Point(21) = {frac_O_x	, O_y		, 0, mesh_param};
Point(22) = {frac_L_x	, O_y		, 0, mesh_param};
Point(23) = {fgp_O_x	, O_y		, 0, mesh_param};
Point(24) = {fgp_L_x	, O_y		, 0, mesh_param};
Point(25) = {fgp_O_x	, frac_O_y	, 0, mesh_param};
Point(26) = {fgp_L_x	, frac_O_y	, 0, mesh_param};
Point(27) = {fgp_O_x	, L_y		, 0, mesh_param};
Point(28) = {fgp_L_x	, L_y		, 0, mesh_param};
Point(29) = {O_x	, fgp_L_y	, 0, mesh_param};
Point(30) = {L_x	, fgp_L_y	, 0, mesh_param};


//+
Line(1) = {1, 23};
//+
Line(2) = {23, 21};
//+
Line(3) = {21, 22};
//+
Line(4) = {22, 24};
//+
Line(5) = {24, 2};
//+
Line(6) = {2, 13};
//+
Line(7) = {13, 10};
//+
Line(8) = {10, 20};
//+
Line(9) = {20, 19};
//+
Line(10) = {19, 9};
//+
Line(11) = {9, 14};
//+
Line(12) = {14, 1};
//+
Line(13) = {13, 30};
//+
Line(14) = {30, 11};
//+
Line(15) = {11, 18};
//+
Line(16) = {18, 17};
//+
Line(17) = {17, 12};
//+
Line(18) = {12, 29};
//+
Line(19) = {29, 14};
//+
Line(20) = {30, 15};
//+
Line(21) = {15, 26};
//+
Line(22) = {26, 6};
//+
Line(23) = {6, 5};
//+
Line(24) = {5, 25};
//+
Line(25) = {25, 16};
//+
Line(26) = {16, 29};
//+
Line(27) = {15, 3};
//+
Line(28) = {3, 28};
//+
Line(29) = {28, 7};
//+
Line(30) = {7, 8};
//+
Line(31) = {8, 27};
//+
Line(32) = {27, 4};
//+
Line(33) = {4, 16};
//+
Line(34) = {28, 26};
//+
Line(35) = {26, 11};
//+
Line(36) = {11, 10};
//+
Line(37) = {10, 24};
//+
Line(38) = {22, 20};
//+
Line(39) = {20, 18};
//+
Line(40) = {18, 6};
//+
Line(41) = {6, 7};
//+
Line(42) = {8, 5};
//+
Line(43) = {5, 17};
//+
Line(44) = {17, 19};
//+
Line(45) = {19, 21};
//+
Line(46) = {23, 9};
//+
Line(47) = {9, 12};
//+
Line(48) = {12, 25};
//+
Line(49) = {25, 27};


//+
Transfinite Line {12, -46, 45, -38, 37, -6} = N1_y Using Progression A1_y;
//+
Transfinite Line {19, 47, 44, 39, 36, 13, 26, 48, 43, 40, 35, 20} = N3_y Using Progression A3_y;
Transfinite Line {33, 49, 42, 41, 34, 27} = 25 Using Progression 1;
//+
Transfinite Line {1, 11, 18, 25, 32, 5, 7, 14, 21, 28} = N2_x Using Progression A2_x;
//+
Transfinite Line {-2, 10, 17, 24, 31} = N4_x Using Progression A4_x;
//+
Transfinite Line {4, -8, -15, -22, -29} = N4_x Using Progression A4_x;
//+
Transfinite Line {3, 9, 16, 23, 30} = N5_x Using Progression A5_x;


//+
Line Loop(1) = {12, 1, 46, 11};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {46, -10, 45, -2};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {45, 3, 38, 9};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {38, -8, 37, -4};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {37, 5, 6, 7};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {19, -11, 47, 18};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {17, -47, -10, -44};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {16, 44, -9, 39};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {15, -39, -8, -36};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {14, 36, -7, 13};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {21, 35, -14, 20};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {22, -40, -15, -35};
//+
Plane Surface(12) = {12};
//+
Line Loop(13) = {23, 43, -16, 40};
//+
Plane Surface(13) = {13};
//+
Line Loop(14) = {24, -48, -17, -43};
//+
Plane Surface(14) = {14};
//+
Line Loop(15) = {25, 26, -18, 48};
//+
Plane Surface(15) = {15};
//+
Line Loop(16) = {32, 33, -25, 49};
//+
Plane Surface(16) = {16};
//+
Line Loop(17) = {31, -49, -24, -42};
//+
Plane Surface(17) = {17};
//+
Line Loop(18) = {30, 42, -23, 41};
//+
Plane Surface(18) = {18};
//+
Line Loop(19) = {29, -41, -22, -34};
//+
Plane Surface(19) = {19};
//+
Line Loop(20) = {28, 34, -21, 27};
//+
Plane Surface(20) = {20};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {10};
//+
Transfinite Surface {9};
//+
Transfinite Surface {8};
//+
Transfinite Surface {7};
//+
Transfinite Surface {6};
//+
Transfinite Surface {15};
//+
Transfinite Surface {14};
//+
Transfinite Surface {13};
//+
Transfinite Surface {12};
//+
Transfinite Surface {11};
//+
Transfinite Surface {20};
//+
Transfinite Surface {19};
//+
Transfinite Surface {18};
//+
Transfinite Surface {17};
//+
Transfinite Surface {16};
//+
Recombine Surface {1};
//+
Recombine Surface {2};
//+
Recombine Surface {3};
//+
Recombine Surface {4};
//+
Recombine Surface {5};
//+
Recombine Surface {10};
//+
Recombine Surface {9};
//+
Recombine Surface {8};
//+
Recombine Surface {7};
//+
Recombine Surface {6};
//+
Recombine Surface {15};
//+
Recombine Surface {14};
//+
Recombine Surface {13};
//+
Recombine Surface {12};
//+
Recombine Surface {11};
//+
Recombine Surface {20};
//+
Recombine Surface {19};
//+
Recombine Surface {18};
//+
Recombine Surface {17};
//+
Recombine Surface {16};
//+

//######################
// PHYSICAL GROUPS
//######################
// FRACTURE
//Physical Surface(1) = {18};
// FREE GAS POCKET
//Physical Surface(2) = {7, 8, 9};
// TOP (SEA FLOOR)
//Physical Line(3) = {32, 31, 30, 29, 28};
// BOTTOM (FAR END)
//Physical Line(4) = {1, 2, 3, 4, 5};
