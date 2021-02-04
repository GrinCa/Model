// Gmsh project created on Sun Nov 29 10:54:50 2020
SetFactory("OpenCASCADE");

plate_thickness=0.004;
L=0.4;
sizemesh = 0.02;

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
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Line(5) = {2, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {3, 7};
//+
Line(8) = {7, 8};
//+
Line(9) = {8, 4};
//+
Line(10) = {8, 5};
//+
Line(11) = {5, 1};
//+
Line(12) = {5, 6};
//+
Line Loop(1) = {6, 8, 10, 12};
//+


Plane Surface(1) = {1};
//+
Line Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {10, 11, -3, -9};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {2, -9, -8, -7};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {1, 7, -6, -5};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {12, -5, -4, -11};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 5, 2, 3, 6, 4};
//+
Volume(1) = {1};


//+
Physical Surface("inplate",2) = {2};
//+
Physical Surface("extplate",3) = {1};
//+
Physical Surface("embedding",4) = {4, 5, 6, 3};
//+
Physical Volume("plate",1) = {1};

