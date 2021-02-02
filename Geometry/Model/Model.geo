// Gmsh project created on Mon Apr 27 08:53:20 2020
//+

sizemesh = 1;

Lcavity = 1;

Lbg = 0.75;
lbg = 0.5;
lpml = 0.5;

Lplate = 0.5;
plate_thickness = 0.05;

coeff = 7;

Point(1) = {plate_thickness, Lplate, Lplate, sizemesh/coeff};
//+
Point(2) = {plate_thickness, -Lplate, Lplate, sizemesh/coeff};
//+
Point(3) = {plate_thickness, -Lplate, -Lplate, sizemesh/coeff};
//+
Point(4) = {plate_thickness, +Lplate, -Lplate, sizemesh/coeff};
//+
Point(13) = {plate_thickness, Lcavity, Lcavity, sizemesh};
//+
Point(14) = {plate_thickness, -Lcavity, Lcavity, sizemesh};
//+
Point(15) = {plate_thickness, -Lcavity, -Lcavity, sizemesh};
//+
Point(16) = {plate_thickness, Lcavity, -Lcavity, sizemesh};

Point(17) = {2*Lcavity, Lcavity, Lcavity, sizemesh};
//+
Point(18) = {2*Lcavity, -Lcavity, Lcavity, sizemesh};
//+
Point(19) = {2*Lcavity, -Lcavity, -Lcavity, sizemesh};
//+
Point(20) = {2*Lcavity, Lcavity, -Lcavity, sizemesh};
//+
Point(21) = {0, Lplate, Lplate, sizemesh/coeff};
//+
Point(22) = {0, -Lplate, Lplate, sizemesh/coeff};
//+
Point(23) = {0, -Lplate, -Lplate, sizemesh/coeff};
//+
Point(24) = {0, Lplate, -Lplate, sizemesh/coeff};
//+
Point(25) = {0, Lbg, -Lbg, sizemesh};
//+
Point(26) = {0, Lbg, Lbg, sizemesh};
//+
Point(27) = {0, -Lbg, Lbg, sizemesh};
//+
Point(28) = {0, -Lbg, -Lbg, sizemesh};
//+
Point(29) = {-2*Lbg, -Lbg, -Lbg, sizemesh};
//+
Point(30) = {-2*Lbg, Lbg, -Lbg, sizemesh};
//+
Point(31) = {-2*Lbg, Lbg, Lbg, sizemesh};
//+
Point(32) = {-2*Lbg, -Lbg, Lbg, sizemesh};
//+
Point(33) = {0, Lbg+lpml, Lbg+lpml, sizemesh};
//+
Point(34) = {0, Lbg+lpml, -(Lbg+lpml), sizemesh};
//+
Point(35) = {0, -(Lbg+lpml), -(Lbg+lpml), sizemesh};
//+
Point(36) = {0, -(Lbg+lpml), Lbg+lpml, sizemesh};
//+
Point(37) = {-(2*Lbg+lpml), Lbg+lpml, Lbg+lpml, sizemesh};
//+
Point(38) = {-(2*Lbg+lpml), Lbg+lpml, -(Lbg+lpml), sizemesh};
//+
Point(39) = {-(2*Lbg+lpml), -(Lbg+lpml), -(Lbg+lpml), sizemesh};
//+
Point(40) = {-(2*Lbg+lpml), -(Lbg+lpml), Lbg+lpml, sizemesh};
//+
Point(41) = {0, Lbg+lpml, Lbg, sizemesh};
//+
Point(42) = {0, Lbg+lpml, -Lbg, sizemesh};
//+
Point(43) = {0, -(Lbg+lpml), -Lbg, sizemesh};
//+
Point(44) = {0, -(Lbg+lpml), Lbg, sizemesh};
//+
Point(45) = {-(2*Lbg+lpml), Lbg+lpml, Lbg, sizemesh};
//+
Point(46) = {-(2*Lbg+lpml), Lbg+lpml, -Lbg, sizemesh};
//+
Point(47) = {-(2*Lbg+lpml), -(Lbg+lpml), -Lbg, sizemesh};
//+
Point(48) = {-(2*Lbg+lpml), -(Lbg+lpml), Lbg, sizemesh};
//+
Point(49) = {-(2*Lbg), Lbg+lpml, Lbg, sizemesh};
//+
Point(50) = {-(2*Lbg), Lbg+lpml, -Lbg, sizemesh};
//+
Point(51) = {-(2*Lbg), -(Lbg+lpml), -Lbg, sizemesh};
//+
Point(52) = {-(2*Lbg), -(Lbg+lpml), Lbg, sizemesh};
//+
Point(53) = {-(2*Lbg), Lbg+lpml, Lbg+lpml, sizemesh};
//+
Point(54) = {-(2*Lbg), Lbg+lpml, -(Lbg+lpml), sizemesh};
//+
Point(55) = {-(2*Lbg), -(Lbg+lpml), -(Lbg+lpml), sizemesh};
//+
Point(56) = {-(2*Lbg), -(Lbg+lpml), Lbg+lpml, sizemesh};
//+
Point(57) = {-2*Lbg, -Lbg, -(Lbg+lpml), sizemesh};
//+
Point(58) = {-2*Lbg, Lbg, -(Lbg+lpml), sizemesh};
//+
Point(59) = {-2*Lbg, Lbg, Lbg+lpml, sizemesh};
//+
Point(60) = {-2*Lbg, -Lbg, Lbg+lpml, sizemesh};
//+
Point(61) = {-(2*Lbg+lpml), -Lbg, -Lbg, sizemesh};
//+
Point(62) = {-(2*Lbg+lpml), Lbg, -Lbg, sizemesh};
//+
Point(63) = {-(2*Lbg+lpml), Lbg, Lbg, sizemesh};
//+
Point(64) = {-(2*Lbg+lpml), -Lbg, Lbg, sizemesh};
//+
Point(65) = {0, -Lbg, -(Lbg+lpml), sizemesh};
//+
Point(66) = {0, Lbg, -(Lbg+lpml), sizemesh};
//+
Point(67) = {0, Lbg, Lbg+lpml, sizemesh};
//+
Point(68) = {0, -Lbg, Lbg+lpml, sizemesh};
//+
Point(69) = {-(2*Lbg+lpml), -Lbg, -(Lbg+lpml), sizemesh};
//+
Point(70) = {-(2*Lbg+lpml), Lbg, -(Lbg+lpml), sizemesh};
//+
Point(71) = {-(2*Lbg+lpml), Lbg, Lbg+lpml, sizemesh};
//+
Point(72) = {-(2*Lbg+lpml), -Lbg, Lbg+lpml, sizemesh};

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
Line(25) = {25, 30};
//+
Line(26) = {30, 31};
//+
Line(27) = {31, 32};
//+
Line(28) = {32, 27};
//+
Line(29) = {27, 26};
//+
Line(30) = {26, 31};
//+
Line(31) = {26, 25};
//+
Line(32) = {30, 29};
//+
Line(33) = {29, 32};
//+
Line(34) = {29, 28};
//+
Line(35) = {28, 25};
//+
Line(36) = {28, 27};
//+
Line(48) = {39, 47};
//+
Line(49) = {47, 48};
//+
Line(50) = {48, 40};
//+
Line(51) = {40, 56};
//+
Line(52) = {56, 36};
//+
Line(53) = {36, 44};
//+
Line(54) = {44, 43};
//+
Line(55) = {43, 35};
//+
Line(56) = {35, 55};
//+
Line(57) = {55, 39};
//+
Line(58) = {47, 51};
//+
Line(59) = {51, 55};
//+
Line(60) = {51, 52};
//+
Line(61) = {48, 52};
//+
Line(62) = {52, 56};
//+
Line(63) = {52, 44};
//+
Line(64) = {43, 51};
//+
Line(65) = {50, 30};
//+
Line(66) = {29, 51};
//+
Line(67) = {32, 52};
//+
Line(68) = {54, 38};
//+
Line(69) = {38, 46};
//+
Line(70) = {46, 50};
//+
Line(71) = {54, 50};
//+
Line(72) = {54, 34};
//+
Line(73) = {34, 42};
//+
Line(74) = {42, 50};
//+
Line(75) = {42, 41};
//+
Line(76) = {41, 33};
//+
Line(77) = {33, 53};
//+
Line(78) = {49, 53};
//+
Line(79) = {53, 37};
//+
Line(80) = {37, 45};
//+
Line(81) = {45, 49};
//+
Line(82) = {45, 46};
//+
Line(83) = {50, 49};
//+
Line(84) = {49, 41};
//+
Line(85) = {70, 58};
//+
Line(86) = {58, 30};
//+
Line(87) = {30, 62};
//+
Line(88) = {62, 70};
//+
Line(89) = {70, 38};
//+
Line(90) = {46, 62};
//+
Line(91) = {69, 39};
//+
Line(92) = {47, 61};
//+
Line(93) = {61, 69};
//+
Line(94) = {69, 70};
//+
Line(95) = {62, 61};
//+
Line(96) = {69, 57};
//+
Line(97) = {57, 29};
//+
Line(98) = {29, 61};
//+
Line(99) = {61, 64};
//+
Line(100) = {64, 72};
//+
Line(101) = {72, 60};
//+
Line(102) = {60, 32};
//+
Line(103) = {32, 64};
//+
Line(104) = {72, 40};
//+
Line(105) = {48, 64};
//+
Line(106) = {56, 60};
//+
Line(107) = {72, 71};
//+
Line(108) = {63, 71};
//+
Line(109) = {71, 59};
//+
Line(110) = {31, 59};
//+
Line(111) = {63, 31};
//+
Line(112) = {62, 63};
//+
Line(113) = {63, 45};
//+
Line(114) = {37, 71};
//+
Line(115) = {59, 53};
//+
Line(116) = {64, 63};
//+
Line(117) = {59, 60};
//+
Line(118) = {59, 67};
//+
Line(119) = {67, 33};
//+
Line(120) = {60, 68};
//+
Line(121) = {68, 36};
//+
Line(122) = {44, 27};
//+
Line(123) = {27, 68};
//+
Line(124) = {68, 67};
//+
Line(125) = {26, 67};
//+
Line(126) = {26, 41};
//+
Line(127) = {66, 25};
//+
Line(128) = {66, 34};
//+
Line(129) = {58, 66};
//+
Line(130) = {65, 35};
//+
Line(131) = {65, 57};
//+
Line(132) = {65, 28};
//+
Line(133) = {43, 28};
//+
Line(134) = {65, 66};
//+
Line(135) = {55, 57};
//+
Line(136) = {58, 57};
//+
Line(137) = {25, 42};
//+
Line(138) = {31, 49};
//+
Line(139) = {58, 54};




//+
Line Loop(1) = {11, 8, 9, 10};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {9, -15, -12, -5};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {13, 14, -6, 12};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {14, 7, -11, -16};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {15, 10, -16, -13};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {8, 5, 6, 7};
//+
Line Loop(7) = {4, 1, 2, 3};
//+
Plane Surface(6) = {6, 7};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {24, -20, 4, 17};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {3, 20, 21, -19};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {2, 19, 22, -18};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {23, -17, 1, 18};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {24, 21, 22, 23};
//+
Plane Surface(12) = {12};
//+
Line Loop(13) = {31, -35, 36, 29};
//+
Plane Surface(13) = {12, 13};
//+
Line Loop(14) = {133, 36, -122, 54};
//+
Plane Surface(14) = {14};
//+
Line Loop(15) = {29, 125, -124, -123};
//+
Plane Surface(15) = {15};
//+
Line Loop(16) = {31, 137, 75, -126};
//+
Plane Surface(16) = {16};
//+
Line Loop(17) = {35, -127, -134, 132};
//+
Plane Surface(17) = {17};
//+
Line Loop(18) = {132, -133, 55, -130};
//+
Plane Surface(18) = {18};
//+
Line Loop(19) = {122, 123, 121, 53};
//+
Plane Surface(19) = {19};
//+
Line Loop(20) = {125, 119, -76, -126};
//+
Plane Surface(20) = {20};
//+
Line Loop(21) = {73, -137, -127, 128};
//+
Plane Surface(21) = {21};
//+
Line Loop(22) = {26, 27, -33, -32};
//+
Plane Surface(22) = {22};
//+
Line Loop(23) = {32, -97, -136, 86};
//+
Plane Surface(23) = {23};
//+
Line Loop(24) = {33, 67, -60, -66};
//+
Plane Surface(24) = {24};
//+
Line Loop(25) = {27, -102, -117, -110};
//+
Plane Surface(25) = {25};
//+
Line Loop(26) = {26, 138, -83, 65};
//+
Plane Surface(26) = {26};
//+
Line Loop(27) = {86, -65, -71, -139};
//+
Plane Surface(27) = {27};
//+
Line Loop(28) = {97, 66, 59, 135};
//+
Plane Surface(28) = {28};
//+
Line Loop(29) = {67, 62, 106, 102};
//+
Plane Surface(29) = {29};
//+
Line Loop(30) = {138, 78, -115, -110};
//+
Plane Surface(30) = {30};
//+
Line Loop(31) = {112, -116, -99, -95};
//+
Plane Surface(31) = {31};
//+
Line Loop(32) = {99, -105, -49, 92};
//+
Plane Surface(32) = {32};
//+
Line Loop(33) = {116, 108, -107, -100};
//+
Plane Surface(33) = {33};
//+
Line Loop(34) = {112, 113, 82, 90};
//+
Plane Surface(34) = {34};
//+
Line Loop(35) = {88, -94, -93, -95};
//+
Plane Surface(35) = {35};
//+
Line Loop(36) = {93, 91, 48, 92};
//+
Plane Surface(36) = {36};
//+
Line Loop(37) = {50, -104, -100, -105};
//+
Plane Surface(37) = {37};
//+
Line Loop(38) = {113, -80, 114, -108};
//+
Plane Surface(38) = {38};
//+
Line Loop(39) = {88, 89, 69, 90};
//+
Plane Surface(39) = {39};
//+
Line Loop(40) = {25, 32, 34, 35};
//+
Plane Surface(40) = {40};
//+
Line Loop(41) = {74, 65, -25, 137};
//+
Plane Surface(41) = {41};
//+
Line Loop(42) = {34, -133, 64, -66};
//+
Plane Surface(42) = {42};
//+
Line Loop(43) = {98, -95, -87, 32};
//+
Plane Surface(43) = {43};
//+
Line Loop(44) = {98, -92, 58, -66};
//+
Plane Surface(44) = {44};
//+
Line Loop(45) = {70, 65, 87, -90};
//+
Plane Surface(45) = {45};
//+
Line Loop(46) = {129, -134, 131, -136};
//+
Plane Surface(46) = {46};
//+
Line Loop(47) = {72, -128, -129, 139};
//+
Plane Surface(47) = {47};
//+
Line Loop(48) = {131, -135, -56, -130};
//+
Plane Surface(48) = {48};
//+
Line Loop(49) = {85, 136, -96, 94};
//+
Plane Surface(49) = {49};
//+
Line Loop(50) = {96, -135, 57, -91};
//+
Plane Surface(50) = {50};
//+
Line Loop(51) = {68, -89, 85, 139};
//+
Plane Surface(51) = {51};
//+
Line Loop(52) = {30, 27, 28, 29};
//+
Plane Surface(52) = {52};
//+
Line Loop(53) = {84, -126, 30, 138};
//+
Plane Surface(53) = {53};
//+
Line Loop(54) = {28, -122, -63, -67};
//+
Plane Surface(54) = {54};
//+
Line Loop(55) = {111, 27, 103, 116};
//+
Plane Surface(55) = {55};
//+
Line Loop(56) = {103, -105, 61, -67};
//+
Plane Surface(56) = {56};
//+
Line Loop(57) = {81, -138, -111, 113};
//+
Plane Surface(57) = {57};
//+
Line Loop(58) = {118, -124, -120, -117};
//+
Plane Surface(58) = {58};
//+
Line Loop(59) = {77, -115, 118, 119};
//+
Plane Surface(59) = {59};
//+
Line Loop(60) = {120, 121, -52, 106};
//+
Plane Surface(60) = {60};
//+
Line Loop(61) = {109, 117, -101, 107};
//+
Plane Surface(61) = {61};
//+
Line Loop(62) = {79, 114, 109, 115};
//+
Plane Surface(62) = {62};
//+
Line Loop(63) = {101, -106, -51, -104};
//+
Plane Surface(63) = {63};
//+
Line Loop(64) = {64, 60, 63, 54};
//+
Plane Surface(64) = {64};
//+
Line Loop(65) = {56, -59, -64, 55};
//+
Plane Surface(65) = {65};
//+
Line Loop(66) = {63, -53, -52, -62};
//+
Plane Surface(66) = {66};
//+
Line Loop(67) = {60, -61, -49, 58};
//+
Plane Surface(67) = {67};
//+
Line Loop(68) = {61, 62, -51, -50};
//+
Plane Surface(68) = {68};
//+
Line Loop(69) = {57, 48, 58, 59};
//+
Plane Surface(69) = {69};
//+
Line Loop(70) = {34, 36, -28, -33};
//+
Plane Surface(70) = {70};
//+
Line Loop(71) = {131, 97, 34, -132};
//+
Plane Surface(71) = {71};
//+
Line Loop(72) = {28, 123, -120, 102};
//+
Plane Surface(72) = {72};
//+
Line Loop(73) = {98, 99, -103, -33};
//+
Plane Surface(73) = {73};
//+
Line Loop(74) = {96, 97, 98, 93};
//+
Plane Surface(74) = {74};
//+
Line Loop(75) = {103, 100, 101, 102};
//+
Plane Surface(75) = {75};
//+
Line Loop(76) = {25, 26, -30, 31};
//+
Plane Surface(76) = {76};
//+
Line Loop(77) = {129, 127, 25, -86};
//+
Plane Surface(77) = {77};
//+
Line Loop(78) = {30, 110, 118, -125};
//+
Plane Surface(78) = {78};
//+
Line Loop(79) = {87, 112, 111, -26};
//+
Plane Surface(79) = {79};
//+
Line Loop(80) = {85, 86, 87, 88};
//+
Plane Surface(80) = {80};
//+
Line Loop(81) = {111, 110, -109, -108};
//+
Plane Surface(81) = {81};
//+
Line Loop(82) = {74, 83, 84, -75};
//+
Plane Surface(82) = {82};
//+
Line Loop(83) = {72, 73, 74, -71};
//+
Plane Surface(83) = {83};
//+
Line Loop(84) = {84, 76, 77, -78};
//+
Plane Surface(84) = {84};
//+
Line Loop(85) = {70, 83, -81, 82};
//+
Plane Surface(85) = {85};
//+
Line Loop(86) = {68, 69, 70, -71};
//+
Plane Surface(86) = {86};
//+
Line Loop(87) = {81, 78, 79, 80};
//+
Plane Surface(87) = {87};
//+
Surface Loop(1) = {4, 3, 5, 2, 1, 6, 7};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {9, 8, 11, 10, 7, 12};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {40, 76, 22, 52, 70, 13, 12};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {46, 17, 23, 40, 77, 71};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {82, 16, 41, 76, 26, 53};
//+
Volume(5) = {5};
//+
Surface Loop(6) = {58, 15, 78, 52, 25, 72};
//+
Volume(6) = {6};
//+
Surface Loop(7) = {54, 24, 70, 14, 64, 42};
//+
Volume(7) = {7};
//+
Surface Loop(8) = {31, 73, 22, 79, 43, 55};
//+
Volume(8) = {8};
//+
Surface Loop(9) = {84, 20, 59, 53, 78, 30};
//+
Volume(9) = {9};
//+
Surface Loop(10) = {33, 61, 25, 55, 75, 81};
//+
Volume(10) = {10};
//+
Surface Loop(11) = {19, 60, 66, 29, 72, 54};
//+
Volume(11) = {11};
//+
Surface Loop(12) = {77, 41, 83, 47, 21, 27};
//+
Volume(12) = {12};
//+
Surface Loop(13) = {35, 49, 43, 23, 74, 80};
//+
Volume(13) = {13};
//+
Surface Loop(14) = {18, 65, 48, 42, 71, 28};
//+
Volume(14) = {14};
//+
Surface Loop(15) = {79, 57, 26, 85, 34, 45};
//+
Volume(15) = {15};
//+
Surface Loop(16) = {44, 67, 32, 56, 24, 73};
//+
Volume(16) = {16};
//+
Surface Loop(17) = {62, 87, 38, 81, 30, 57};
//+
Volume(17) = {17};
//+
Surface Loop(18) = {63, 68, 37, 56, 29, 75};
//+
Volume(18) = {18};
//+
Surface Loop(19) = {36, 50, 69, 28, 44, 74};
//+
Volume(19) = {19};
//+
Surface Loop(20) = {45, 86, 51, 39, 80, 27};
//+
Volume(20) = {20};




Physical Volume("cavity", 1) = {1};
//+
Physical Volume("plate", 2) = {2};
//+
Physical Volume("background", 3) = {3};
//+
Physical Volume("PML", 4) = {14, 19, 16, 18, 10, 17, 9, 6, 11, 8, 15, 5, 12, 4, 13, 7, 20};
//+
Physical Surface("PlateCavity", 5) = {7};
//+
Physical Surface("PlateBG", 6) = {12};
//+
Physical Surface("embedding", 7) = {9, 8, 11, 10};
//+
Physical Surface("BGPML",8) = {70, 22, 76, 52, 40};
