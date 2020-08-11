// Gmsh project created on Mon Apr 27 08:53:20 2020
//+

//+

SetFactory("OpenCASCADE");

sizemesh = 0.15;

l = 0.2;
L = 0.3;
plate_thickness = 0.005;

lpml = 0.5;

coeff = 5;


Point(1) = {0, l, l, sizemesh/coeff};
//+
Point(2) = {0, -l, l, sizemesh/coeff};
//+
Point(3) = {0, -l, -l, sizemesh/coeff};
//+
Point(4) = {0, +l, -l, sizemesh/coeff};
//+
Point(13) = {0, L, L, sizemesh};
//+
Point(14) = {0, -L, L, sizemesh};
//+
Point(15) = {0, -L, -L, sizemesh};
//+
Point(16) = {0, L, -L, sizemesh};

Point(17) = {2*L, L, L, sizemesh};
//+
Point(18) = {2*L, -L, L, sizemesh};
//+
Point(19) = {2*L, -L, -L, sizemesh};
//+
Point(20) = {2*L, L, -L, sizemesh};
//+
Point(21) = {-plate_thickness, l, l, sizemesh/coeff};
//+
Point(22) = {-plate_thickness, -l, l, sizemesh/coeff};
//+
Point(23) = {-plate_thickness, -l, -l, sizemesh/coeff};
//+
Point(24) = {-plate_thickness, l, -l, sizemesh/coeff};
//+


//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {13, 14};
//+
Line(6) = {14, 15};
//+
Line(7) = {15, 16};
//+
Line(8) = {16, 13};
//+

Line(9) = {13, 17};
//+
Line(10) = {17, 20};
//+
Line(11) = {20, 16};
//+
Line(12) = {14, 18};
//+
Line(13) = {18, 19};
//+
Line(14) = {19, 15};
//+
Line(15) = {18, 17};
//+
Line(16) = {19, 20};
//+
Line(17) = {1, 21};
//+
Line(18) = {2, 22};
//+
Line(19) = {3, 23};
//+
Line(20) = {4, 24};
//+
Line(21) = {24, 23};
//+
Line(22) = {23, 22};
//+
Line(23) = {22, 21};
//+
Line(24) = {21, 24};


//+
Line Loop(1) = {9, 10, 11, 8};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {5, 12, 15, -9};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {13, 16, -10, -15};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {11, -7, -14, 16};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {5, 6, 7, 8};
//+
Line Loop(6) = {1, 2, 3, 4};
//+
Plane Surface(5) = {5, 6};
//+
Line Loop(7) = {12, 13, 14, -6};
//+
Plane Surface(6) = {7};
//+
Plane Surface(7) = {6};
//+
Line Loop(8) = {22, 23, 24, 21};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {17, -23, -18, -1};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {18, -22, -19, -2};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {19, -21, -20, -3};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {4, 17, 24, -20};
//+
Plane Surface(12) = {12};
//+
Surface Loop(1) = {1, 2, 5, 6, 3, 4, 7};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {8, 10, 9, 12, 11, 7};
//+
Volume(2) = {2};
//+



//+
Point(25) = {0, -(L+lpml), -(L+lpml), sizemesh};
//+
Point(26) = {0, -(L+lpml), (L+lpml), sizemesh};
//+
Point(27) = {0, (L+lpml), -(L+lpml), sizemesh};
//+
Point(28) = {0, (L+lpml), (L+lpml), sizemesh};
//+
Point(29) = {0, -(L+lpml), -L, sizemesh};
//+
Point(30) = {0, -(L+lpml), L, sizemesh};
//+
Point(31) = {0, (L+lpml), -L, sizemesh};
//+
Point(32) = {0, (L+lpml), L, sizemesh};
//+
Point(33) = {0, -L, -(L+lpml), sizemesh};
//+
Point(34) = {0, -L, (L+lpml), sizemesh};
//+
Point(35) = {0, L, -(L+lpml), sizemesh};
//+
Point(36) = {0, L, (L+lpml), sizemesh};
//+

Point(37) = {2*L, -(L+lpml), -(L+lpml), sizemesh};
//+
Point(38) = {2*L, -(L+lpml), (L+lpml), sizemesh};
//+
Point(39) = {2*L, (L+lpml), -(L+lpml), sizemesh};
//+
Point(40) = {2*L, (L+lpml), (L+lpml), sizemesh};
//+
Point(41) = {2*L, -(L+lpml), -L, sizemesh};
//+
Point(42) = {2*L, -(L+lpml), L, sizemesh};
//+
Point(43) = {2*L, (L+lpml), -L, sizemesh};
//+
Point(44) = {2*L, (L+lpml), L, sizemesh};
//+
Point(45) = {2*L, -L, -(L+lpml), sizemesh};
//+
Point(46) = {2*L, -L, (L+lpml), sizemesh};
//+
Point(47) = {2*L, L, -(L+lpml), sizemesh};
//+
Point(48) = {2*L, L, (L+lpml), sizemesh};
//+

Point(49) = {2*L+lpml, -(L+lpml), -(L+lpml), sizemesh};
//+
Point(50) = {2*L+lpml, -(L+lpml), (L+lpml), sizemesh};
//+
Point(51) = {2*L+lpml, (L+lpml), -(L+lpml), sizemesh};
//+
Point(52) = {2*L+lpml, (L+lpml), (L+lpml), sizemesh};
//+
Point(53) = {2*L+lpml, -(L+lpml), -L, sizemesh};
//+
Point(54) = {2*L+lpml, -(L+lpml), L, sizemesh};
//+
Point(55) = {2*L+lpml, (L+lpml), -L, sizemesh};
//+
Point(56) = {2*L+lpml, (L+lpml), L, sizemesh};
//+
Point(57) = {2*L+lpml, -L, -(L+lpml), sizemesh};
//+
Point(58) = {2*L+lpml, -L, (L+lpml), sizemesh};
//+
Point(59) = {2*L+lpml, L, -(L+lpml), sizemesh};
//+
Point(60) = {2*L+lpml, L, (L+lpml), sizemesh};
//+

Point(61) = {2*L+lpml, -L, -L, sizemesh};
//+
Point(62) = {2*L+lpml, -L, L, sizemesh};
//+
Point(63) = {2*L+lpml, L, -L, sizemesh};
//+
Point(64) = {2*L+lpml, L, L, sizemesh};



//+
Line(25) = {13, 32};
//+
Line(26) = {32, 28};
//+
Line(27) = {28, 36};
//+
Line(28) = {36, 13};
//+
Line(29) = {36, 34};
//+
Line(30) = {34, 14};
//+
Line(31) = {14, 30};
//+
Line(32) = {30, 26};
//+
Line(33) = {26, 34};
//+
Line(34) = {26, 38};
//+
Line(35) = {38, 50};
//+
Line(36) = {50, 54};
//+
Line(37) = {54, 42};
//+
Line(38) = {42, 38};
//+
Line(39) = {42, 30};
//+
Line(40) = {34, 46};
//+
Line(41) = {46, 18};
//+
Line(42) = {18, 62};
//+
Line(43) = {62, 58};
//+
Line(44) = {58, 46};
//+
Line(45) = {46, 38};
//+
Line(46) = {54, 62};
//+
Line(47) = {18, 42};
//+
Line(48) = {50, 58};
//+
Line(49) = {58, 60};
//+
Line(50) = {60, 52};
//+
Line(51) = {40, 52};
//+
Line(52) = {40, 48};
//+
Line(53) = {48, 17};
//+
Line(54) = {17, 44};
//+
Line(55) = {44, 40};
//+
Line(56) = {44, 56};
//+
Line(57) = {56, 52};
//+
Line(58) = {17, 64};
//+
Line(59) = {64, 60};
//+
Line(60) = {60, 48};
//+
Line(61) = {64, 56};
//+
Line(62) = {48, 36};
//+
Line(63) = {28, 40};
//+
Line(64) = {44, 32};
//+
Line(65) = {48, 46};
//+
Line(66) = {30, 29};
//+
Line(67) = {29, 25};
//+
Line(68) = {25, 33};
//+
Line(69) = {33, 15};
//+
Line(70) = {15, 29};
//+
Line(71) = {33, 35};
//+
Line(72) = {35, 27};
//+
Line(73) = {27, 31};
//+
Line(74) = {31, 16};
//+
Line(75) = {16, 35};
//+
Line(76) = {35, 47};
//+
Line(77) = {47, 20};
//+
Line(78) = {20, 43};
//+
Line(79) = {43, 39};
//+
Line(80) = {39, 47};
//+
Line(81) = {39, 27};
//+
Line(82) = {31, 43};
//+
Line(83) = {43, 55};
//+
Line(84) = {55, 51};
//+
Line(85) = {51, 39};
//+
Line(86) = {20, 63};
//+
Line(87) = {63, 59};
//+
Line(88) = {59, 47};
//+
Line(89) = {63, 55};
//+
Line(90) = {63, 64};
//+
Line(91) = {56, 55};
//+
Line(92) = {43, 44};
//+
Line(93) = {31, 32};
//+
Line(94) = {63, 61};
//+
Line(95) = {61, 19};
//+
Line(96) = {41, 53};
//+
Line(97) = {53, 49};
//+
Line(98) = {49, 37};
//+
Line(99) = {37, 41};
//+
Line(100) = {41, 19};
//+
Line(101) = {61, 53};
//+
Line(102) = {61, 57};
//+
Line(103) = {57, 49};
//+
Line(104) = {37, 45};
//+
Line(105) = {45, 57};
//+
Line(106) = {57, 59};
//+
Line(107) = {59, 51};
//+
Line(108) = {47, 45};
//+
Line(109) = {45, 33};
//+
Line(110) = {25, 37};
//+
Line(111) = {42, 41};
//+
Line(112) = {53, 54};
//+
Line(113) = {62, 61};
//+
Line(114) = {19, 45};
//+
Line(115) = {62, 64};
//+
Line(116) = {41, 29};
//+
Line Loop(13) = {110, 99, 116, 67};
//+
Plane Surface(13) = {13};
//+
Line Loop(14) = {99, 96, 97, 98};
//+
Plane Surface(14) = {14};
//+
Line Loop(15) = {96, 112, 37, 111};
//+
Plane Surface(15) = {15};
//+
Line Loop(16) = {37, 38, 35, 36};
//+
Plane Surface(16) = {16};
//+
Line Loop(17) = {34, -38, 39, 32};
//+
Plane Surface(17) = {17};
//+
Line Loop(18) = {66, -116, -111, 39};
//+
Plane Surface(18) = {18};
//+
Line Loop(19) = {109, 69, -14, 114};
//+
Plane Surface(19) = {19};
//+
Line Loop(20) = {105, -102, 95, 114};
//+
Plane Surface(20) = {20};
//+
Line Loop(21) = {95, -13, 42, 113};
//+
Plane Surface(21) = {21};
//+
Line Loop(22) = {41, 42, 43, 44};
//+
Plane Surface(22) = {22};
//+
Line Loop(23) = {12, -41, -40, 30};
//+
Plane Surface(23) = {23};
//+
Line Loop(24) = {76, 77, 11, 75};
//+
Plane Surface(24) = {24};
//+
Line Loop(25) = {77, 86, 87, 88};
//+
Plane Surface(25) = {25};
//+
Line Loop(26) = {10, 86, 90, -58};
//+
Plane Surface(26) = {26};
//+
Line Loop(27) = {53, 58, 59, 60};
//+
Plane Surface(27) = {27};
//+
Line Loop(28) = {9, -53, 62, 28};
//+
Plane Surface(28) = {28};
//+
Line Loop(29) = {73, 82, 79, 81};
//+
Plane Surface(29) = {29};
//+
Line Loop(30) = {79, -85, -84, -83};
//+
Plane Surface(30) = {30};
//+
Line Loop(31) = {92, 56, 91, -83};
//+
Plane Surface(31) = {31};
//+
Line Loop(32) = {55, 51, -57, -56};
//+
Plane Surface(32) = {32};
//+
Line Loop(33) = {64, 26, 63, -55};
//+
Plane Surface(33) = {33};
//+
Line Loop(34) = {93, -64, -92, -82};
//+
Plane Surface(34) = {34};
//+
Line Loop(35) = {81, -72, 76, -80};
//+
Plane Surface(35) = {35};
//+
Line Loop(36) = {80, -88, 107, 85};
//+
Plane Surface(36) = {36};
//+
Line Loop(37) = {88, 108, 105, 106};
//+
Plane Surface(37) = {37};
//+
Line Loop(38) = {105, 103, 98, 104};
//+
Plane Surface(38) = {38};
//+
Line Loop(39) = {109, -68, 110, 104};
//+
Plane Surface(39) = {39};
//+
Line Loop(40) = {76, 108, 109, 71};
//+
Plane Surface(40) = {40};
//+
Line Loop(41) = {82, -78, 11, -74};
//+
Plane Surface(41) = {41};
//+
Line Loop(42) = {83, -89, -86, 78};
//+
Plane Surface(42) = {42};
//+
Line Loop(43) = {86, 94, 95, 16};
//+
Plane Surface(43) = {43};
//+
Line Loop(44) = {95, -100, 96, -101};
//+
Plane Surface(44) = {44};
//+
Line Loop(45) = {14, 70, -116, 100};
//+
Plane Surface(45) = {45};
//+
Line Loop(46) = {12, 47, 39, -31};
//+
Plane Surface(46) = {46};
//+
Line Loop(47) = {42, -46, 37, -47};
//+
Plane Surface(47) = {47};
//+
Line Loop(48) = {15, 58, -115, -42};
//+
Plane Surface(48) = {48};
//+
Line Loop(49) = {54, 56, -61, -58};
//+
Plane Surface(49) = {49};
//+
Line Loop(50) = {64, -25, 9, 54};
//+
Plane Surface(50) = {50};
//+
Line Loop(51) = {63, 52, 62, -27};
//+
Plane Surface(51) = {51};
//+
Line Loop(52) = {51, -50, 60, -52};
//+
Plane Surface(52) = {52};
//+
Line Loop(53) = {60, 65, -44, 49};
//+
Plane Surface(53) = {53};
//+
Line Loop(54) = {44, 45, 35, 48};
//+
Plane Surface(54) = {54};
//+
Line Loop(55) = {40, 45, -34, 33};
//+
Plane Surface(55) = {55};
//+
Line Loop(56) = {40, -65, 62, 29};
//+
Plane Surface(56) = {56};
//+
Line Loop(57) = {6, 70, -66, -31};
//+
Plane Surface(57) = {57};
//+
Line Loop(58) = {30, 31, 32, 33};
//+
Plane Surface(58) = {58};
//+
Line Loop(59) = {28, 5, -30, -29};
//+
Plane Surface(59) = {59};
//+
Line Loop(60) = {26, 27, 28, 25};
//+
Plane Surface(60) = {60};
//+
Line Loop(61) = {93, -25, -8, -74};
//+
Plane Surface(61) = {61};
//+
Line Loop(62) = {73, 74, 75, 72};
//+
Plane Surface(62) = {62};
//+
Line Loop(63) = {75, -71, 69, 7};
//+
Plane Surface(63) = {63};
//+
Line Loop(64) = {69, 70, 67, 68};
//+
Plane Surface(64) = {64};
//+
Line Loop(65) = {41, 47, 38, -45};
//+
Plane Surface(65) = {65};
//+
Line Loop(66) = {41, 15, -53, 65};
//+
Plane Surface(66) = {66};
//+
Line Loop(67) = {55, 52, 53, 54};
//+
Plane Surface(67) = {67};
//+
Line Loop(68) = {92, -54, 10, 78};
//+
Plane Surface(68) = {68};
//+
Line Loop(69) = {79, 80, 77, 78};
//+
Plane Surface(69) = {69};
//+
Line Loop(70) = {77, -16, 114, -108};
//+
Plane Surface(70) = {70};
//+
Line Loop(71) = {114, -104, 99, 100};
//+
Plane Surface(71) = {71};
//+
Line Loop(72) = {13, -100, -111, -47};
//+
Plane Surface(72) = {72};
//+
Line Loop(73) = {43, -48, 36, 46};
//+
Plane Surface(73) = {73};
//+
Line Loop(74) = {49, -59, -115, 43};
//+
Plane Surface(74) = {74};
//+
Line Loop(75) = {57, -50, -59, 61};
//+
Plane Surface(75) = {75};
//+
Line Loop(76) = {91, -89, 90, 61};
//+
Plane Surface(76) = {76};
//+
Line Loop(77) = {84, -107, -87, 89};
//+
Plane Surface(77) = {77};
//+
Line Loop(78) = {87, -106, -102, -94};
//+
Plane Surface(78) = {78};
//+
Line Loop(79) = {102, 103, -97, -101};
//+
Plane Surface(79) = {79};
//+
Line Loop(80) = {113, 101, 112, 46};
//+
Plane Surface(80) = {80};
//+
Line Loop(81) = {115, -90, 94, -113};
//+
Plane Surface(81) = {81};
//+
Surface Loop(3) = {61, 34, 1, 50, 68, 41};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {33, 60, 51, 50, 67, 28};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {35, 29, 62, 24, 41, 69};
//+
Volume(5) = {5};
//+
Surface Loop(6) = {30, 36, 77, 25, 42, 69};
//+
Volume(6) = {6};
//+
Surface Loop(7) = {31, 76, 68, 26, 49, 42};
//+
Volume(7) = {7};
//+
Surface Loop(8) = {27, 32, 52, 75, 67, 49};
//+
Volume(8) = {8};
//+
Surface Loop(9) = {56, 59, 2, 23, 28, 66};
//+
Volume(9) = {9};
//+
Surface Loop(10) = {22, 53, 74, 66, 48, 27};
//+
Volume(10) = {10};
//+
Surface Loop(11) = {21, 48, 81, 3, 26, 43};
//+
Volume(11) = {11};
//+
Surface Loop(12) = {21, 80, 15, 72, 44, 47};
//+
Volume(12) = {12};
//+
Surface Loop(13) = {16, 54, 73, 65, 22, 47};
//+
Volume(13) = {13};
//+
Surface Loop(14) = {17, 55, 58, 46, 23, 65};
//+
Volume(14) = {14};
//+
Surface Loop(15) = {79, 38, 14, 20, 44, 71};
//+
Volume(15) = {15};
//+
Surface Loop(16) = {13, 39, 64, 45, 19, 71};
//+
Volume(16) = {16};
//+
Surface Loop(17) = {6, 46, 18, 57, 72, 45};
//+
Volume(17) = {17};
//+
Surface Loop(18) = {20, 37, 78, 70, 43, 25};
//+
Volume(18) = {18};
//+
Surface Loop(19) = {40, 63, 70, 4, 19, 24};
//+
Volume(19) = {19};
//+

Physical Volume("cavity",1) = {1};
//+
Physical Volume("plate",2) = {2};
//+
Physical Surface("embedding",3) = {9, 12, 11, 10};
//+
Physical Surface("coupling",4) = {7};
//+
Physical Surface("ext_pressure", 5) = {8};
//+
Physical Volume("PML",6) = {5, 4, 3, 8, 7, 6, 18, 11, 10, 9, 14, 13, 12, 17, 16, 15, 19};




//+
