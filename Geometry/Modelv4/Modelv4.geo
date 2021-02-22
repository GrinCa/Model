// Gmsh project created on Sun Nov 29 10:54:50 2020
SetFactory("OpenCASCADE");

plate_thickness=0.004;
L=0.4;
sizemesh = 0.01;

Point(1) = {-L/2, -L/2, 0, sizemesh};
//+
Point(2) = {-L/2, L/2, 0, sizemesh};
//+
Point(3) = {L/2, L/2, 0, sizemesh};
//+
Point(4) = {L/2, -L/2, 0, sizemesh};
//+
Point(5) = {-L/2, -L/2, plate_thickness, sizemesh};
//+
Point(6) = {-L/2, L/2, plate_thickness, sizemesh};
//+
Point(7) = {L/2, L/2, plate_thickness, sizemesh};
//+
Point(8) = {L/2, -L/2, plate_thickness, sizemesh};



//+
Line(1) = {6, 7};
//+
Line(2) = {7, 8};
//+
Line(3) = {8, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 4};
//+
Line(8) = {4, 1};
//+
Line(9) = {1, 2};
//+
Line(10) = {7, 3};
//+
Line(11) = {8, 4};
//+
Line(12) = {5, 1};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {6, 7, 8, 9};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {6, -10, -1, 5};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {7, -11, -2, 10};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {11, 8, -12, -3};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {12, 9, -5, -4};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 3, 2, 4, 5, 6};
//+
Volume(1) = {1};
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
Transfinite Surface {6};
//+

Physical Surface("inplate", 2) = {1};
//+
Physical Surface("extplate", 3) = {2};
//+
Physical Surface("embedding", 4) = {3, 4, 5, 6};
//+
Physical Volume("plate", 1) = {1};


