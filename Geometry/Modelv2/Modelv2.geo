// Gmsh project created on Mon Apr 27 08:53:20 2020
//+

sizemesh = 1;



Lbg = 0.75;
lbg = 0.5;
lpml = 0.5;

Lplate = 0.5;
plate_thickness = 0.05;

Lcavity = Lplate+0.1;

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
Point(73) = {plate_thickness, -(Lcavity+lpml), -(Lcavity+lpml), sizemesh};
//+
Point(74) = {plate_thickness, -(Lcavity+lpml), Lcavity+lpml, sizemesh};
//+
Point(75) = {plate_thickness, Lcavity+lpml, -(Lcavity+lpml), sizemesh};
//+
Point(76) = {plate_thickness, Lcavity+lpml, Lcavity+lpml, sizemesh};
//+
Point(77) = {2*Lcavity, -(Lcavity+lpml), -(Lcavity+lpml), sizemesh};
//+
Point(78) = {2*Lcavity, -(Lcavity+lpml), Lcavity+lpml, sizemesh};
//+
Point(79) = {2*Lcavity, Lcavity+lpml, -(Lcavity+lpml), sizemesh};
//+
Point(80) = {2*Lcavity, Lcavity+lpml, Lcavity+lpml, sizemesh};
//+
Point(81) = {2*Lcavity+lpml, -(Lcavity+lpml), -(Lcavity+lpml), sizemesh};
//+
Point(82) = {2*Lcavity+lpml, -(Lcavity+lpml), Lcavity+lpml, sizemesh};
//+
Point(83) = {2*Lcavity+lpml, Lcavity+lpml, -(Lcavity+lpml), sizemesh};
//+
Point(84) = {2*Lcavity+lpml, Lcavity+lpml, Lcavity+lpml, sizemesh};
//+
Point(85) = {2*Lcavity+lpml, -(Lcavity+lpml), -Lcavity, sizemesh};
//+
Point(86) = {2*Lcavity+lpml, -(Lcavity+lpml), Lcavity, sizemesh};
//+
Point(87) = {2*Lcavity+lpml, Lcavity+lpml, -Lcavity, sizemesh};
//+
Point(88) = {2*Lcavity+lpml, Lcavity+lpml, Lcavity, sizemesh};
//+
Point(89) = {2*Lcavity+lpml, -Lcavity, -(Lcavity+lpml), sizemesh};
//+
Point(90) = {2*Lcavity+lpml, -Lcavity, Lcavity+lpml, sizemesh};
//+
Point(91) = {2*Lcavity+lpml, Lcavity, -(Lcavity+lpml), sizemesh};
//+
Point(92) = {2*Lcavity+lpml, Lcavity, Lcavity+lpml, sizemesh};
//+
Point(93) = {2*Lcavity, -Lcavity, -(Lcavity+lpml), sizemesh};
//+
Point(94) = {2*Lcavity, -Lcavity, Lcavity+lpml, sizemesh};
//+
Point(95) = {2*Lcavity, Lcavity, -(Lcavity+lpml), sizemesh};
//+
Point(96) = {2*Lcavity, Lcavity, Lcavity+lpml, sizemesh};
//+
Point(97) = {2*Lcavity, -(Lcavity+lpml), -Lcavity, sizemesh};
//+
Point(98) = {2*Lcavity, -(Lcavity+lpml), Lcavity, sizemesh};
//+
Point(99) = {2*Lcavity, Lcavity+lpml, -Lcavity, sizemesh};
//+
Point(100) = {2*Lcavity, Lcavity+lpml, Lcavity, sizemesh};
//+
Point(101) = {2*Lcavity+lpml, -Lcavity, -Lcavity, sizemesh};
//+
Point(102) = {2*Lcavity+lpml, -Lcavity, Lcavity, sizemesh};
//+
Point(103) = {2*Lcavity+lpml, Lcavity, -Lcavity, sizemesh};
//+
Point(104) = {2*Lcavity+lpml, Lcavity, Lcavity, sizemesh};
//+
Point(105) = {plate_thickness, -(Lcavity+lpml), -Lcavity, sizemesh};
//+
Point(106) = {plate_thickness, -(Lcavity+lpml), Lcavity, sizemesh};
//+
Point(107) = {plate_thickness, Lcavity+lpml, -Lcavity, sizemesh};
//+
Point(108) = {plate_thickness, Lcavity+lpml, Lcavity, sizemesh};
//+
Point(109) = {plate_thickness, -Lcavity, -(Lcavity+lpml), sizemesh};
//+
Point(110) = {plate_thickness, -Lcavity, Lcavity+lpml, sizemesh};
//+
Point(111) = {plate_thickness, Lcavity, -(Lcavity+lpml), sizemesh};
//+
Point(112) = {plate_thickness, Lcavity, Lcavity+lpml, sizemesh};





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
Line(140) = {76, 108};
//+
Line(141) = {108, 107};
//+
Line(142) = {107, 75};
//+
Line(143) = {75, 79};
//+
Line(144) = {79, 83};
//+
Line(145) = {83, 87};
//+
Line(146) = {87, 88};
//+
Line(147) = {88, 84};
//+
Line(149) = {108, 100};
//+
Line(150) = {100, 99};
//+
Line(151) = {99, 107};
//+
Line(152) = {99, 87};
//+
Line(153) = {99, 79};
//+
Line(154) = {100, 88};
//+
Line(155) = {100, 80};
//+
Line(156) = {76, 112};
//+
Line(157) = {112, 110};
//+
Line(158) = {110, 74};
//+
Line(159) = {74, 106};
//+
Line(160) = {106, 105};
//+
Line(161) = {105, 73};
//+
Line(162) = {73, 109};
//+
Line(163) = {109, 111};
//+
Line(164) = {111, 75};
//+
Line(165) = {106, 14};
//+
Line(166) = {14, 110};
//+
Line(167) = {112, 13};
//+
Line(168) = {13, 108};
//+
Line(169) = {112, 96};
//+
Line(170) = {96, 92};
//+
Line(171) = {92, 104};
//+
Line(172) = {104, 103};
//+
Line(173) = {103, 91};
//+
Line(174) = {91, 95};
//+
Line(175) = {95, 111};
//+
Line(176) = {83, 91};
//+
Line(177) = {95, 79};
//+
Line(178) = {87, 103};
//+
Line(179) = {20, 99};
//+
Line(180) = {20, 103};
//+
Line(181) = {20, 95};
//+
Line(182) = {17, 104};
//+
Line(183) = {17, 96};
//+
Line(184) = {92, 84};
//+
Line(185) = {80, 96};
//+
Line(186) = {104, 88};
//+
Line(187) = {90, 92};
//+
Line(188) = {90, 82};
//+
Line(189) = {82, 86};
//+
Line(190) = {86, 85};
//+
Line(191) = {85, 81};
//+
Line(192) = {81, 77};
//+
Line(193) = {77, 73};
//+
Line(194) = {105, 97};
//+
Line(195) = {97, 98};
//+
Line(198) = {78, 82};

//+
Line(199) = {74, 78};
//+
Line(200) = {98, 106};
//+
Line(201) = {98, 86};
//+
Line(202) = {98, 78};
//+
Line(203) = {98, 18};
//+
Line(204) = {18, 102};
//+
Line(205) = {102, 86};
//+
Line(206) = {94, 78};
//+
Line(207) = {94, 18};
//+
Line(208) = {94, 90};
//+
Line(209) = {90, 102};
//+
Line(210) = {110, 94};
//+
Line(211) = {102, 101};
//+
Line(212) = {101, 89};
//+
Line(213) = {89, 93};
//+
Line(214) = {93, 109};
//+
Line(215) = {19, 93};
//+
Line(216) = {101, 19};
//+
Line(217) = {101, 85};
//+
Line(218) = {81, 89};
//+
Line(219) = {93, 77};

//+
Line(220) = {89, 91};
//+
Line(221) = {95, 93};
//+
Line(222) = {101, 103};
//+
Line(223) = {104, 102};
//+
Line(224) = {94, 96};
//+
Line(225) = {17, 100};
//+
Line(226) = {16, 107};
//+
Line(227) = {16, 111};
//+
Line(228) = {97, 77};
//+
Line(229) = {97, 85};
//+
Line(230) = {15, 105};
//+
Line(231) = {15, 109};
//+
Line(232) = {19, 97};
//+
Line(233) = {76, 80};
//+
Line(234) = {80, 84};



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




Line Loop(88) = {141, -226, 8, 168};
//+
Plane Surface(88) = {88};
//+
Line Loop(89) = {142, -164, -227, 226};
//+
Plane Surface(89) = {89};
//+
Line Loop(90) = {227, -163, -231, 7};
//+
Plane Surface(90) = {90};
//+
Line Loop(91) = {231, -162, -161, -230};
//+
Plane Surface(91) = {91};
//+
Line Loop(92) = {6, 230, -160, 165};
//+
Plane Surface(92) = {92};
//+
Line Loop(93) = {165, 166, 158, 159};
//+
Plane Surface(93) = {93};
//+
Line Loop(94) = {166, -157, 167, 5};
//+
Plane Surface(94) = {94};
//+
Line Loop(95) = {140, -168, -167, -156};
//+
Plane Surface(95) = {95};

//+
Line Loop(96) = {150, -179, -10, 225};
//+
Plane Surface(96) = {96};
//+
Line Loop(97) = {153, -177, -181, 179};
//+
Plane Surface(97) = {97};
//+
Line Loop(98) = {181, 221, -215, 16};
//+
Plane Surface(98) = {98};

//+
Line Loop(99) = {215, 219, -228, -232};
//+
Plane Surface(99) = {99};
//+
Line Loop(100) = {13, 232, 195, 203};
//+
Plane Surface(100) = {100};
//+
Line Loop(101) = {202, -206, 207, -203};
//+
Plane Surface(101) = {101};
//+
Line Loop(102) = {224, -183, -15, -207};
//+
Plane Surface(102) = {102};
//+
Line Loop(103) = {155, 185, -183, 225};
//+
Plane Surface(103) = {103};
//+
Line Loop(104) = {172, -222, -211, -223};
//+
Plane Surface(104) = {104};
//+
Line Loop(105) = {146, -186, 172, -178};
//+
Plane Surface(105) = {105};
//+
Line Loop(106) = {145, 178, 173, -176};
//+
Plane Surface(106) = {106};
//+
Line Loop(107) = {173, -220, -212, 222};
//+
Plane Surface(107) = {107};
//+
Line Loop(108) = {212, -218, -191, -217};
//+
Plane Surface(108) = {108};
//+
Line Loop(109) = {211, 217, -190, -205};
//+
Plane Surface(109) = {109};
//+
Line Loop(110) = {209, 205, -189, -188};
//+
Plane Surface(110) = {110};
//+
Line Loop(111) = {209, -223, -171, -187};
//+
Plane Surface(111) = {111};
//+
Line Loop(112) = {184, -147, -186, -171};
//+
Plane Surface(112) = {112};

//+
Line Loop(113) = {233, 185, -169, -156};
//+
Plane Surface(113) = {113};
//+
Line Loop(114) = {234, -184, -170, -185};
//+
Plane Surface(114) = {114};
//+
Line Loop(115) = {170, -187, -208, 224};
//+
Plane Surface(115) = {115};
//+
Line Loop(116) = {208, 188, -198, -206};
//+
Plane Surface(116) = {116};
//+
Line Loop(117) = {210, 206, -199, -158};
//+
Plane Surface(117) = {117};
//+
Line Loop(118) = {169, -224, -210, -157};
//+
Plane Surface(118) = {118};
//+
Line Loop(119) = {168, 149, -225, -9};
//+
Plane Surface(119) = {119};
//+
Line Loop(120) = {154, -186, -182, 225};
//+
Plane Surface(120) = {120};
//+
Line Loop(121) = {223, -204, 15, 182};
//+
Plane Surface(121) = {121};
//+
Line Loop(122) = {204, 205, -201, 203};
//+
Plane Surface(122) = {122};
//+
Line Loop(123) = {200, 165, 12, -203};
//+
Plane Surface(123) = {123};
//+
Line Loop(124) = {151, -226, -11, 179};
//+
Plane Surface(124) = {124};
//+
Line Loop(125) = {152, 178, -180, 179};
//+
Plane Surface(125) = {125};
//+
Line Loop(126) = {180, -222, 216, 16};
//+
Plane Surface(126) = {126};
//+
Line Loop(127) = {216, 232, 229, -217};
//+
Plane Surface(127) = {127};
//+
Line Loop(128) = {14, 230, 194, -232};
//+
Plane Surface(128) = {128};
//+
Line Loop(129) = {193, 162, -214, 219};
//+
Plane Surface(129) = {129};
//+
Line Loop(130) = {213, 219, -192, 218};
//+
Plane Surface(130) = {130};
//+
Line Loop(131) = {220, 174, 221, -213};
//+
Plane Surface(131) = {131};
//+
Line Loop(132) = {144, 176, 174, 177};
//+
Plane Surface(132) = {132};
//+
Line Loop(133) = {143, -177, 175, 164};
//+
Plane Surface(133) = {133};
//+
Line Loop(134) = {175, -163, -214, -221};
//+
Plane Surface(134) = {134};
//+
Line Loop(135) = {194, 228, 193, -161};
//+
Plane Surface(135) = {135};
//+
Line Loop(136) = {229, 191, 192, -228};
//+
Plane Surface(136) = {136};
//+
Line Loop(137) = {201, 190, -229, 195};
//+
Plane Surface(137) = {137};
//+
Line Loop(138) = {198, 189, -201, 202};
//+
Plane Surface(138) = {138};
//+
Line Loop(139) = {199, -202, 200, -159};
//+
Plane Surface(139) = {139};
//+
Line Loop(140) = {200, 160, 194, 195};
//+
Plane Surface(140) = {140};
//+
Line Loop(141) = {231, -214, -215, 14};
//+
Plane Surface(141) = {141};
//+
Line Loop(142) = {216, 215, -213, -212};
//+
Plane Surface(142) = {142};
//+
Line Loop(143) = {204, 211, 216, -13};
//+
Plane Surface(143) = {143};
//+
Line Loop(144) = {208, 209, -204, -207};
//+
Plane Surface(144) = {144};
//+
Line Loop(145) = {210, 207, -12, 166};
//+
Plane Surface(145) = {145};
//+
Line Loop(146) = {227, -175, -181, 11};
//+
Plane Surface(146) = {146};
//+
Line Loop(147) = {180, 173, 174, -181};
//+
Plane Surface(147) = {147};
//+
Line Loop(148) = {182, 172, -180, -10};
//+
Plane Surface(148) = {148};
//+
Line Loop(149) = {170, 171, -182, 183};
//+
Plane Surface(149) = {149};
//+
Line Loop(150) = {169, -183, -9, -167};
//+
Plane Surface(150) = {150};
//+
Line Loop(151) = {142, 143, -153, 151};
//+
Plane Surface(151) = {151};
//+
Line Loop(152) = {152, -145, -144, -153};
//+
Plane Surface(152) = {152};
//+
Line Loop(153) = {154, -146, -152, -150};
//+
Plane Surface(153) = {153};
//+
Line Loop(154) = {234, -147, -154, 155};
//+
Plane Surface(154) = {154};
//+
Line Loop(155) = {233, -155, -149, -140};
//+
Plane Surface(155) = {155};
//+
Line Loop(156) = {149, 150, 151, -141};
//+
Plane Surface(156) = {156};
//+
Surface Loop(21) = {151, 89, 133, 146, 124, 97};
//+
Volume(21) = {21};
//+
Surface Loop(22) = {132, 152, 106, 125, 147, 97};
//+
Volume(22) = {22};
//+
Surface Loop(23) = {131, 107, 98, 126, 147, 142};
//+
Volume(23) = {23};
//+
Surface Loop(24) = {108, 130, 136, 142, 99, 127};
//+
Volume(24) = {24};
//+
Surface Loop(25) = {129, 135, 91, 99, 128, 141};
//+
Volume(25) = {25};
//+
Surface Loop(26) = {141, 134, 90, 4, 98, 146};
//+
Volume(26) = {26};
//+
Surface Loop(27) = {109, 137, 122, 127, 143, 100};
//+
Volume(27) = {27};
//+
Surface Loop(28) = {110, 138, 116, 122, 144, 101};
//+
Volume(28) = {28};
//+
Surface Loop(29) = {115, 111, 149, 121, 144, 102};
//+
Volume(29) = {29};
//+
Surface Loop(30) = {149, 114, 154, 112, 120, 103};
//+
Volume(30) = {30};
//+
Surface Loop(31) = {139, 117, 93, 123, 101, 145};
//+
Volume(31) = {31};
//+
Surface Loop(32) = {118, 94, 102, 150, 2, 145};
//+
Volume(32) = {32};
//+
Surface Loop(33) = {155, 113, 95, 103, 119, 150};
//+
Volume(33) = {33};
//+
Surface Loop(34) = {104, 121, 5, 148, 143, 126};
//+
Volume(34) = {34};
//+
Surface Loop(35) = {140, 92, 123, 3, 128, 100};
//+
Volume(35) = {35};
//+
Surface Loop(36) = {120, 105, 153, 125, 148, 96};
//+
Volume(36) = {36};
//+
Surface Loop(37) = {119, 156, 88, 96, 124, 1};
//+
Volume(37) = {37};


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
//+
Physical Volume("PML2",9) = {33, 37, 21, 22, 36, 30, 34, 26, 23, 24, 25, 27, 28, 35, 31, 32, 29};
