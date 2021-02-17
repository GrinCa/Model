

plate_thickness = 0.004;
l_ROOMfield = 0.8;
L_PML = 0.5;
l_plate = 0.4;
l1 = 0.8;
coeff = 8;
sizemesh = 0.02;

Point(1) = {l1/2, l1/2, 0, sizemesh*coeff};
//+
Point(2) = {-l1/2, l1/2, 0, sizemesh*coeff};
//+
Point(3) = {-l1/2, -l1/2, 0, sizemesh*coeff};
//+
Point(4) = {l1/2, -l1/2, 0, sizemesh*coeff};
//+
Point(9) = {l1/2+L_PML, l1/2, 0, sizemesh*coeff};
//+
Point(10) = {-l1/2-L_PML, l1/2, 0, sizemesh*coeff};
//+
Point(11) = {-l1/2-L_PML, -l1/2, 0, sizemesh*coeff};
//+
Point(12) = {l1/2+L_PML, -l1/2, 0, sizemesh*coeff};
//+
Point(13) = {l1/2, l1/2+L_PML, 0, sizemesh*coeff};
//+
Point(14) = {-l1/2, l1/2+L_PML, 0, sizemesh*coeff};
//+
Point(15) = {-l1/2, -l1/2-L_PML, 0, sizemesh*coeff};
//+
Point(16) = {l1/2, -l1/2-L_PML, 0, sizemesh*coeff};
//+
Point(17) = {l1/2+L_PML, l1/2+L_PML, 0, sizemesh*coeff};
//+
Point(18) = {-l1/2-L_PML, l1/2+L_PML, 0, sizemesh*coeff};
//+
Point(19) = {-l1/2-L_PML, -l1/2-L_PML, 0, sizemesh*coeff};
//+
Point(20) = {l1/2+L_PML, -l1/2-L_PML, 0, sizemesh*coeff};
//*
Point(21) = {l1/2+L_PML, l1/2, l_ROOMfield, sizemesh*coeff};
//+
Point(22) = {-l1/2-L_PML, l1/2, l_ROOMfield, sizemesh*coeff};
//+
Point(23) = {-l1/2-L_PML, -l1/2, l_ROOMfield, sizemesh*coeff};
//+
Point(24) = {l1/2+L_PML, -l1/2, l_ROOMfield, sizemesh*coeff};
//+
Point(25) = {l1/2, l1/2+L_PML, l_ROOMfield, sizemesh*coeff};
//+
Point(26) = {-l1/2, l1/2+L_PML, l_ROOMfield, sizemesh*coeff};
//+
Point(27) = {-l1/2, -l1/2-L_PML, l_ROOMfield, sizemesh*coeff};
//+
Point(28) = {l1/2, -l1/2-L_PML, l_ROOMfield, sizemesh*coeff};
//+
Point(29) = {l1/2+L_PML, l1/2+L_PML, l_ROOMfield, sizemesh*coeff};
//+
Point(30) = {-l1/2-L_PML, l1/2+L_PML, l_ROOMfield, sizemesh*coeff};
//+
Point(31) = {-l1/2-L_PML, -l1/2-L_PML, l_ROOMfield, sizemesh*coeff};
//+
Point(32) = {l1/2+L_PML, -l1/2-L_PML, l_ROOMfield, sizemesh*coeff};
//*
Point(33) = {l1/2+L_PML, l1/2, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(34) = {-l1/2-L_PML, l1/2, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(35) = {-l1/2-L_PML, -l1/2, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(36) = {l1/2+L_PML, -l1/2, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(37) = {l1/2, l1/2+L_PML, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(38) = {-l1/2, l1/2+L_PML, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(39) = {-l1/2, -l1/2-L_PML, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(40) = {l1/2, -l1/2-L_PML, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(41) = {l1/2+L_PML, l1/2+L_PML, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(42) = {-l1/2-L_PML, l1/2+L_PML, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(43) = {-l1/2-L_PML, -l1/2-L_PML, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(44) = {l1/2+L_PML, -l1/2-L_PML, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(45) = {l1/2, -l1/2, l_ROOMfield, sizemesh*coeff};
//+
Point(46) = {-l1/2, -l1/2, l_ROOMfield, sizemesh*coeff};
//+
Point(47) = {-l1/2, l1/2, l_ROOMfield, sizemesh*coeff};
//+
Point(48) = {l1/2, l1/2, l_ROOMfield, sizemesh*coeff};
//+
Point(49) = {l1/2, -l1/2, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(50) = {-l1/2, -l1/2, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(51) = {-l1/2, l1/2, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(52) = {l1/2, l1/2, l_ROOMfield+L_PML, sizemesh*coeff};
//+
Point(53) = {l_plate/2, l_plate/2, -plate_thickness, sizemesh};
//+
Point(54) = {l_plate/2, -l_plate/2, -plate_thickness, sizemesh};
//+
Point(55) = {-l_plate/2, -l_plate/2, -plate_thickness, sizemesh};
//+
Point(56) = {-l_plate/2, l_plate/2, -plate_thickness, sizemesh};
//+
Point(57) = {l_plate/2, l_plate/2, 0, sizemesh};
//+
Point(58) = {l_plate/2, -l_plate/2, 0, sizemesh};
//+
Point(59) = {-l_plate/2, -l_plate/2, 0, sizemesh};
//+
Point(60) = {-l_plate/2, l_plate/2, 0, sizemesh};

//+
Line(1) = {41, 37};
//+
Line(2) = {37, 25};
//+
Line(3) = {25, 29};
//+
Line(4) = {29, 21};
//+
Line(5) = {21, 33};
//+
Line(6) = {33, 41};
//+
Line(7) = {41, 29};
//+
Line(8) = {29, 17};
//+
Line(9) = {17, 13};
//+
Line(10) = {13, 25};
//+
Line(11) = {9, 1};
//+
Line(15) = {2, 10};
//+
Line(16) = {2, 14};
//+
Line(17) = {14, 18};
//+
Line(18) = {18, 10};
//+
Line(19) = {2, 1};
//+
Line(20) = {1, 13};
//+
Line(21) = {13, 14};
//+
Line(22) = {14, 26};
//+
Line(23) = {26, 38};
//+
Line(24) = {38, 42};
//+
Line(25) = {42, 34};
//+
Line(26) = {34, 22};
//+
Line(27) = {22, 47};
//+
Line(28) = {47, 51};
//+
Line(29) = {51, 38};
//+
Line(30) = {51, 34};
//+
Line(31) = {26, 47};
//+
Line(32) = {2, 47};
//+
Line(33) = {22, 10};
//+
Line(34) = {15, 19};
//+
Line(35) = {19, 11};
//+
Line(36) = {3, 11};
//+
Line(37) = {3, 15};
//+
Line(38) = {3, 46};
//+
Line(39) = {46, 50};
//+
Line(40) = {50, 35};
//+
Line(41) = {35, 23};
//+
Line(42) = {23, 46};
//+
Line(43) = {11, 23};
//+
Line(44) = {23, 31};
//+
Line(45) = {19, 31};
//+
Line(46) = {31, 43};
//+
Line(47) = {43, 39};
//+
Line(48) = {39, 27};
//+
Line(49) = {27, 15};
//+
Line(50) = {31, 27};
//+
Line(51) = {39, 50};
//+
Line(52) = {43, 35};
//+
Line(53) = {46, 27};
//+
Line(54) = {23, 22};
//+
Line(55) = {35, 34};
//+
Line(56) = {2, 3};
//+
Line(57) = {11, 10};
//+
Line(58) = {47, 46};
//+
Line(59) = {50, 51};
//+
Line(60) = {22, 30};
//+
Line(61) = {30, 18};
//+
Line(62) = {30, 42};
//+
Line(63) = {30, 26};
//+
Line(64) = {26, 25};
//+
Line(65) = {37, 38};
//+
Line(66) = {51, 33};
//+
Line(67) = {21, 47};
//+
Line(68) = {48, 52};
//+
Line(69) = {48, 1};
//+
Line(70) = {9, 17};
//+
Line(71) = {21, 9};
//+
Line(72) = {9, 12};
//+
Line(73) = {12, 20};
//+
Line(74) = {20, 16};
//+
Line(75) = {16, 28};
//+
Line(76) = {28, 40};
//+
Line(77) = {40, 49};
//+
Line(78) = {49, 45};
//+
Line(79) = {45, 24};
//+
Line(85) = {4, 1};
//+
Line(86) = {4, 3};
//+
Line(87) = {16, 15};
//+
Line(88) = {28, 27};
//+
Line(89) = {39, 40};
//+
Line(90) = {44, 40};
//+
Line(91) = {44, 36};
//+
Line(92) = {36, 33};
//+
Line(93) = {52, 49};
//+
Line(94) = {36, 49};
//+
Line(95) = {45, 28};
//+
Line(96) = {24, 36};
//+
Line(97) = {24, 32};
//+
Line(98) = {32, 28};
//+
Line(99) = {44, 32};
//+
Line(100) = {45, 48};
//+
Line(101) = {24, 21};
//+
Line(102) = {24, 12};
//+
Line(103) = {12, 4};
//+
Line(104) = {4, 45};
//+
Line(105) = {4, 16};
//+
Line(106) = {49, 50};
//+
Line(107) = {46, 45};
//+
Line(108) = {37, 52};
//+
Line(109) = {32, 20};
//+
Line(115) = {57, 53};
//+
Line(116) = {57, 60};
//+
Line(117) = {60, 56};
//+
Line(118) = {56, 53};
//+
Line(119) = {53, 54};
//+
Line(120) = {54, 58};
//+
Line(121) = {58, 59};
//+
Line(122) = {59, 55};
//+
Line(123) = {55, 54};
//+
Line(124) = {58, 57};
//+
Line(125) = {60, 59};
//+
Line(126) = {55, 56};


Line(110) = {48, 25};
//+
Recursive Delete {
  Line{67}; 
}
//+
Line(111) = {21, 48};
//+
Line(112) = {48, 47};
//+
Line Loop(3) = {112, 58, 107, 100};
//+
Plane Surface(3) = {3};
//+
Recursive Delete {
  Line{66}; 
}
//+
Line(113) = {33, 52};
//+
Line(114) = {52, 51};
//+
Line Loop(4) = {114, -59, -106, -93};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {87, -37, -86, 105};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {88, -53, 107, 95};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {89, 77, 106, -51};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {19, 20, 21, -16};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {112, -31, 64, -110};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {108, 114, 29, -65};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {36, 57, -15, 56};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {27, 58, -42, 54};
//+
Plane Surface(12) = {12};
//+
Line Loop(13) = {40, 55, -30, -59};
//+
Plane Surface(13) = {13};
//+
Line Loop(14) = {11, -85, -103, -72};
//+
Plane Surface(14) = {14};
//+
Line Loop(15) = {101, 111, -100, 79};
//+
Plane Surface(15) = {15};
//+
Line Loop(16) = {92, 113, 93, -94};
//+
Plane Surface(16) = {16};
//+
Line Loop(17) = {32, 58, -38, -56};
//+
Plane Surface(17) = {17};
//+
Line Loop(18) = {43, 54, 33, -57};
//+
Plane Surface(18) = {18};
//+
Line Loop(19) = {69, -85, 104, 100};
//+
Plane Surface(19) = {19};
//+
Line Loop(20) = {71, 72, -102, 101};
//+
Plane Surface(20) = {20};
//+
Line Loop(21) = {10, -64, -22, -21};
//+
Plane Surface(21) = {21};
//+
Line Loop(22) = {32, -112, 69, -19};
//+
Plane Surface(22) = {22};
//+
Line Loop(23) = {104, -107, -38, -86};
//+
Plane Surface(23) = {23};
//+
Line Loop(24) = {49, -87, 75, 88};
//+
Plane Surface(24) = {24};
//+
Line Loop(25) = {74, -105, -103, 73};
//+
Plane Surface(25) = {25};
//+
Line Loop(26) = {98, -95, 79, 97};
//+
Plane Surface(26) = {26};
//+
Line Loop(27) = {90, 77, -94, -91};
//+
Plane Surface(27) = {27};
//+
Line Loop(28) = {34, 35, -36, 37};
//+
Plane Surface(28) = {28};
//+
Line Loop(29) = {42, 53, -50, -44};
//+
Plane Surface(29) = {29};
//+
Line Loop(30) = {47, 51, 40, -52};
//+
Plane Surface(30) = {30};
//+
Line Loop(31) = {18, -15, 16, 17};
//+
Plane Surface(31) = {31};
//+
Line Loop(32) = {63, 31, -27, 60};
//+
Plane Surface(32) = {32};
//+
Line Loop(33) = {30, -25, -24, -29};
//+
Plane Surface(33) = {33};
//+
Line Loop(34) = {9, -20, -11, 70};
//+
Plane Surface(34) = {34};
//+
Line Loop(35) = {3, 4, 111, 110};
//+
Plane Surface(35) = {35};
//+
Line Loop(36) = {1, 108, -113, 6};
//+
Plane Surface(36) = {36};
//+
Line Loop(37) = {9, 10, 3, 8};
//+
Plane Surface(37) = {37};
//+
Line Loop(38) = {2, 3, -7, 1};
//+
Plane Surface(38) = {38};
//+
Line Loop(39) = {64, -2, 65, -23};
//+
Plane Surface(39) = {39};
//+
Line Loop(40) = {62, -24, -23, -63};
//+
Plane Surface(40) = {40};
//+
Line Loop(41) = {61, -17, 22, -63};
//+
Plane Surface(41) = {41};
//+
Line Loop(42) = {33, -15, 32, -27};
//+
Plane Surface(42) = {42};
//+
Line Loop(43) = {27, 28, 30, 26};
//+
Plane Surface(43) = {43};
//+
Line Loop(44) = {28, -114, -68, 112};
//+
Plane Surface(44) = {44};
//+
Line Loop(45) = {68, -113, -5, 111};
//+
Plane Surface(45) = {45};
//+
Line Loop(46) = {69, -11, -71, 111};
//+
Plane Surface(46) = {46};
//+
Line Loop(47) = {5, 6, 7, 4};
//+
Plane Surface(47) = {47};
//+
Line Loop(48) = {96, -91, 99, -97};
//+
Plane Surface(48) = {48};
//+
Line Loop(49) = {90, -76, -98, -99};
//+
Plane Surface(49) = {49};
//+
Line Loop(50) = {89, -76, 88, -48};
//+
Plane Surface(50) = {50};
//+
Line Loop(51) = {47, 48, -50, 46};
//+
Plane Surface(51) = {51};
//+
Line Loop(52) = {45, 50, 49, 34};
//+
Plane Surface(52) = {52};
//+
Plane Surface(53) = {28};
//+
Line Loop(53) = {38, 53, 49, -37};
//+
Plane Surface(54) = {53};
//+
Line Loop(54) = {104, 79, 102, 103};
//+
Plane Surface(55) = {54};
//+
Line Loop(55) = {98, -75, -74, -109};
//+
Plane Surface(56) = {55};
//+
Line Loop(56) = {94, 78, 79, 96};
//+
Plane Surface(57) = {56};
//+
Line Loop(57) = {78, -107, 39, -106};
//+
Plane Surface(58) = {57};
//+
Line Loop(58) = {38, -42, -43, -36};
//+
Plane Surface(59) = {58};
//+
Line Loop(59) = {39, 40, 41, 42};
//+
Plane Surface(60) = {59};
//+
Line Loop(60) = {48, -53, 39, -51};
//+
Plane Surface(61) = {60};
//+
Line Loop(61) = {45, -44, -43, -35};
//+
Plane Surface(62) = {61};
//+
Line Loop(62) = {60, 61, 18, -33};
//+
Plane Surface(63) = {62};
//+
Line Loop(63) = {60, 62, 25, 26};
//+
Plane Surface(64) = {63};
//+
Line Loop(64) = {29, -23, 31, 28};
//+
Plane Surface(65) = {64};
//+
Line Loop(65) = {68, -108, 2, -110};
//+
Plane Surface(66) = {65};
//+
Line Loop(66) = {20, 10, -110, 69};
//+
Plane Surface(67) = {66};
//+
Line Loop(67) = {71, 70, -8, 4};
//+
Plane Surface(68) = {67};
//+
Line Loop(68) = {41, 44, 46, 52};
//+
Plane Surface(69) = {68};
//+
Surface Loop(1) = {40, 64, 33, 43, 32, 65};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {47, 45, 36, 38, 35, 66};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {39, 10, 9, 65, 44, 66};
//+
Volume(3) = {3};
//+
Line Loop(69) = {31, -32, 16, 22};
//+
Plane Surface(70) = {69};
//+
Surface Loop(4) = {41, 63, 31, 32, 42, 70};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {21, 8, 22, 9, 70, 67};
//+
Volume(5) = {5};
//+
Surface Loop(6) = {34, 37, 68, 35, 67, 46};
//+
Volume(6) = {6};
//+
Line Loop(70) = {68, 93, 78, 100};
//+
Plane Surface(71) = {70};
//+
Line Loop(71) = {28, -59, -39, -58};
//+
Plane Surface(72) = {71};
//+
Line Loop(72) = {101, 5, -92, -96};
//+
Plane Surface(73) = {72};
//+
Surface Loop(7) = {45, 71, 57, 15, 73, 16};
//+
Volume(7) = {7};
//+
Line Loop(73) = {26, -54, -41, 55};
//+
Plane Surface(74) = {73};
//+
Surface Loop(8) = {13, 74, 60, 72, 12, 43};
//+
Volume(8) = {8};
//+
Line Loop(74) = {77, 78, 95, 76};
//+
Plane Surface(75) = {74};
//+
Surface Loop(9) = {27, 49, 48, 26, 57, 75};
//+
Volume(9) = {9};
//+
Surface Loop(10) = {58, 7, 50, 6, 61, 75};
//+
Volume(10) = {10};
//+
Surface Loop(11) = {4, 3, 71, 58, 44, 72};
//+
Volume(11) = {11};
//+
Surface Loop(13) = {20, 14, 55, 15, 19, 46};
//+
Volume(13) = {13};
//+
Line Loop(75) = {105, 75, -95, -104};
//+
Plane Surface(76) = {75};
//+
Surface Loop(14) = {54, 5, 24, 76, 23, 6};
//+
Volume(14) = {14};
//+
Surface Loop(15) = {52, 62, 28, 54, 59, 29};
//+
Volume(15) = {15};
//+
Surface Loop(16) = {60, 61, 51, 30, 69, 29};
//+
Volume(16) = {16};
//+
Surface Loop(17) = {42, 17, 11, 18, 59, 12};
//+
Volume(17) = {17};
//+
Line Loop(76) = {73, -109, -97, 102};
//+
Plane Surface(77) = {76};
//+
Surface Loop(18) = {25, 56, 77, 55, 26, 76};
//+
Volume(18) = {18};
//+
Line Loop(77) = {116, 117, 118, -115};
//+
Plane Surface(78) = {77};
//+
Line Loop(78) = {117, -126, -122, -125};
//+
Plane Surface(79) = {78};
//+
Line Loop(79) = {123, 120, 121, 122};
//+
Plane Surface(80) = {79};
//+
Line Loop(80) = {119, 120, 124, 115};
//+
Plane Surface(81) = {80};
//+
Line Loop(81) = {118, 119, -123, 126};
//+
Plane Surface(82) = {81};
//+
Line Loop(82) = {116, 125, -121, 124};
//+
Plane Surface(85) = {82};
//+
Line Loop(83) = {19, -85, 86, -56};
//+
Plane Surface(86) = {82, 83};
//+
Surface Loop(19) = {82, 78, 79, 80, 81, 85};
//+
Volume(19) = {19};
//+
Surface Loop(20) = {19, 23, 17, 22, 3, 86, 82, 78, 79, 80, 81};
//+
Surface Loop(21) = {23, 19, 22, 3, 17, 85, 86};
//+
Volume(20) = {21};


//+
Physical Surface("embedding",1) = {78, 81, 80, 79};
//+
Physical Surface("in_plate",2) = {82};
//+
Physical Surface("ext_plate",3) = {2};
//+
Physical Surface("ROOM_PML",4) = {23, 17, 22, 19, 3};
//+
Physical Volume("ROOM",6) = {12};
//+
Physical Volume("PML",7) = {14, 15, 16, 10, 9, 18, 7, 13, 11, 8, 17, 4, 1, 5, 3, 2, 6};






