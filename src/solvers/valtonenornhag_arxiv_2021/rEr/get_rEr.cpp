// Copyright (c) 2021 Marcus Valtonen Örnhag
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <Eigen/Dense>
#include <cmath>  // max
#include <vector>
#include <complex>
#include "get_valtonenornhag_arxiv_2021.hpp"
#include "quartic.hpp"
#include "relpose.hpp"
#include "normalize2dpts.hpp"

namespace DronePoseLib {
namespace ValtonenOrnhagArxiv2021Extra {
    inline Eigen::Matrix<double, 5, 1> calculate_coefficients(const Eigen::Matrix<double, 22, 1> &x);
    inline Eigen::Matrix3d create_M(const double r, const Eigen::Matrix<double, 22, 1> &x);

    std::vector<RelPose> get_rEr(
        const Eigen::MatrixXd &p1,
        const Eigen::MatrixXd &p2,
        const Eigen::Matrix3d &R1,
        const Eigen::Matrix3d &R2,
        double focal_length
    ) {
        // This is a 3-point method
        int nbr_pts = 3;

        // We expect inhomogenous input data, i.e. p1 and p2 are 2x3 matrices
        assert(p1.rows() == 2);
        assert(p2.rows() == 2);
        assert(p1.cols() == nbr_pts);
        assert(p2.cols() == nbr_pts);

        // Compute normalization matrix
        double scale1 = normalize2dpts(p1);
        double scale2 = normalize2dpts(p2);
        double scale = std::max(scale1, scale2);

        // Make a copy
        Eigen::MatrixXd pp1 = p1;
        Eigen::MatrixXd pp2 = p2;
        pp1 *= scale;
        pp2 *= scale;
        focal_length *= scale;

        // Compute relative rotation
        Eigen::Matrix3d R = R2 * R1.transpose();

        Eigen::Matrix<double, 22, 1> x;
        x << pp1.col(0), pp1.col(1), pp1.col(2),
             pp2.col(0), pp2.col(1), pp2.col(2),
             R.col(0), R.col(1), R.col(2), focal_length;

        Eigen::Matrix<double, 5, 1> coeffs = calculate_coefficients(x);

        Eigen::Vector4cd r  = solve_quartic(
            coeffs(1) / coeffs(0),
            coeffs(2) / coeffs(0),
            coeffs(3) / coeffs(0),
            coeffs(4) / coeffs(0));

        // Construct putative output
        double thresh = 1e-12;
        Eigen::Matrix3d M;
        std::vector<RelPose> output;
        RelPose relpose;
        Eigen::Matrix3d Q;
        Eigen::Vector3d kinv;
        Eigen::DiagonalMatrix<double, 3> Kinv;
        Eigen::Matrix3d skew_t;

        // Householder is just 500 ns slower than null3, but is significantly more stable
        for (int i=0; i < r.size(); i++) {
            if (std::abs(std::imag(r[i])) < thresh) {
                M = create_M(std::real(r[i]), x);
                auto qr = M.transpose().colPivHouseholderQr();
                Q = qr.householderQ();
                relpose.t = Q.col(2);
                relpose.f = focal_length / scale;
                relpose.r = std::real(r[i]) * std::pow(scale, 2);

                // Compute fundamental matrix
                kinv << 1.0 / relpose.f, 1.0 / relpose.f, 1.0;
                Kinv = kinv.asDiagonal();
                skew_t << 0, -relpose.t(2), relpose.t(1),
                          relpose.t(2), 0, -relpose.t(0),
                         -relpose.t(1), relpose.t(0), 0;
                relpose.F = Kinv * skew_t * R * Kinv;

                output.push_back(relpose);
            }
        }

        return output;
    }

    inline Eigen::Matrix<double, 5, 1> calculate_coefficients(const Eigen::Matrix<double, 22, 1> &x) {
        Eigen::Matrix<double, 5, 1> coeffs;
        double t607 = x[11];
        double t585 = std::pow(t607, 2);
        double t608 = x[10];
        double t586 = std::pow(t608, 2);
        double t1141 = t585+t586;
        double t615 = x[3];
        double t593 = std::pow(t615, 2);
        double t616 = x[2];
        double t594 = std::pow(t616, 2);
        double t564 = t593+t594;
        double t617 = x[1];
        double t595 = std::pow(t617, 2);
        double t618 = x[0];
        double t596 = std::pow(t618, 2);
        double t565 = t595+t596;
        double t1143 = t564+t565;
        double t499 = t564*t565;
        double t250 = t1141*t1143+t499;
        double t597 = x[21];
        double t1209 = t250*t597;
        double t601 = x[17];
        double t575 = t615*t601;
        double t604 = x[14];
        double t576 = t616*t604;
        double t1142 = t575+t576;
        double t961 = t604*t618;
        double t981 = t601*t617;
        double t533 = t961+t981;
        double t351 = t533-t1142;
        double t1130 = t1141*t351;
        double t1062 = t1142*t596;
        double t330 = t1142*t617-t564*t601;
        double t229 = t330*t617-t564*t961+t1062;
        double t108 = t229-t1130;
        double t614 = x[4];
        double t592 = std::pow(t614, 2);
        double t1193 = t592+t595;
        double t541 = t592+t564;
        double t912 = t533*t541;
        double t158 = t1142*t1193+t1062-t912;
        double t36 = t1141*t158+t229*t592;
        double t613 = x[5];
        double t591 = std::pow(t613, 2);
        double t1208 = t108*t591+t36;
        double t602 = x[16];
        double t895 = t615-t617;
        double t1160 = t602*t895;
        double t605 = x[13];
        double t894 = -t616+t618;
        double t352 = t605*t894-t1160;
        double t1150 = t1141*t352;
        double t957 = t605*t616;
        double t972 = t602*t615;
        double t531 = t957+t972;
        double t1056 = t531*t596;
        double t331 = t531*t617-t564*t602;
        double t1176 = t331*t617+t1056;
        double t955 = t605*t618;
        double t230 = -t564*t955+t1176;
        double t109 = t230-t1150;
        double t610 = x[8];
        double t572 = t604*t614;
        double t220 = t1143*t572-t912;
        double t609 = x[9];
        double t587 = std::pow(t609, 2);
        double t588 = std::pow(t610, 2);
        double t899 = -t588-t587;
        double t741 = -t220-(t533-t572)*t899;
        double t85 = t158-t1130;
        double t1207 = t608*t741+t610*t85;
        double t603 = x[15];
        double t1159 = t603*t895;
        double t606 = x[12];
        double t353 = t606*t894-t1159;
        double t240 = t352*t608+t353*t607;
        double t1194 = t899*t240;
        double t1155 = t1143*t614;
        double t971 = t602*t617;
        double t534 = t955+t971;
        double t911 = t534*t541;
        double t221 = t1155*t605-t911;
        double t573 = t605*t614;
        double t421 = t534-t573;
        double t754 = t1141*t421;
        double t124 = t221-t754;
        double t950 = t606*t618;
        double t965 = t603*t617;
        double t535 = t950+t965;
        double t910 = t535*t541;
        double t222 = t1155*t606-t910;
        double t574 = t606*t614;
        double t422 = t535-t574;
        double t737 = t1141*t422;
        double t126 = t222-t737;
        double t1183 = t592*t531;
        double t160 = t531*t595+t1056+t1183-t911;
        double t952 = t606*t616;
        double t966 = t603*t615;
        double t532 = t952+t966;
        double t498 = t532*t596;
        double t161 = t1193*t532+t498-t910;
        double t68 = t160*t608+t161*t607;
        double t1206 = -t124*t610-t126*t609+t1194+t68;
        double t1028 = t564*t607;
        double t358 = t533*t1028;
        double t1205 = (t108*t609+t358)*t591+t36*t609;
        double t763 = t564+t1141;
        double t1203 = t533*t763;
        double t1202 = t597*t763;
        double t1029 = t564*t603;
        double t1045 = t532*t617;
        double t332 = -t1029+t1045;
        double t315 = t332*t617;
        double t231 = -t564*t950+t315+t498;
        double t90 = t230*t608+t231*t607;
        double t1201 = -t90-t1194;
        double t1152 = t899*t352;
        double t1199 = -t124*t609+(t160+t1152)*t607;
        double t960 = t605*t603;
        double t979 = t602*t606;
        double t518 = -t960+t979;
        double t1074 = t518*t594;
        double t956 = t605*t617;
        double t435 = t532*t956;
        double t1014 = t596*t605;
        double t436 = t532*t1014;
        double t437 = t532*t955;
        double t106 = t436+(-t518*t593-t1074+t435)*t617+(-t1045*t602-t437)*t614;
        double t444 = t531*t574;
        double t469 = t518*t615;
        double t477 = t616*t518;
        double t791 = t605*t574;
        double t940 = t608*t614;
        double t111 = t352*t605*t940-t607*((t791-t469)*t618+(t574*t602+t477)*t617-t444);
        double t438 = t532*t573;
        double t473 = t614*t518;
        double t486 = t602*t532;
        double t1136 = -t438+t437+(t473+t486)*t617;
        double t598 = x[20];
        double t566 = t598*t607;
        double t147 = t160*t566;
        double t929 = t614*t618;
        double t646 = -t929+t565;
        double t930 = t614*t617;
        double t320 = -t602*t930+t605*t646;
        double t951 = t606*t617;
        double t441 = t531*t951;
        double t1013 = t596*t606;
        double t442 = t531*t1013;
        double t953 = t606*t615;
        double t549 = t605*t953;
        double t746 = t564*t606;
        double t1016 = t594*t606;
        double t785 = t605*t1016;
        double t839 = t614*t477;
        double t959 = t605*t608;
        double t58 = t230*t959-(t442+(-t785-(-t473+t549)*t615)*t618+(-t602*t746+t441-t839)*t617)*t607;
        double t1006 = t598*t608;
        double t416 = t532-t574;
        double t1133 = t565*t416;
        double t1026 = t564*t614;
        double t633 = -t1026*t606+t532*t592;
        double t178 = t633+t1133;
        double t923 = t1006*t221+t178*t566;
        double t1198 = (-t106*t614+t1136*t1141+t923)*t609-t899*t111-(t147+(-t320*t614+t754)*t531)*t610-t614*t58;
        double t782 = t606*t929;
        double t693 = t531*t782;
        double t931 = t614*t616;
        double t789 = t605*t931;
        double t700 = t603*t789;
        double t967 = t603*t614;
        double t103 = t442-t693+t617*(t441+t1074-t700+(-t602*t967+t469)*t615);
        double t947 = t607*t614;
        double t112 = -t606*t353*t947+((t573*t603-t477)*t617-(-t791-t469)*t618-t438)*t608;
        double t443 = t531*t950;
        double t487 = t603*t531;
        double t1137 = -t444+t443+(-t473+t487)*t617;
        double t148 = t161*t1006;
        double t321 = -t603*t930+t606*t646;
        double t747 = t564*t605;
        double t954 = t606*t607;
        double t59 = (t436+(-t785-(t473+t549)*t615)*t618+(-t603*t747+t435+t839)*t617)*t608-t231*t954;
        double t415 = t531-t573;
        double t1131 = t565*t415;
        double t177 = -t1026*t605+t1131+t1183;
        double t924 = t1006*t177+t222*t566;
        double t1196 = t899*t112+(-t103*t614+t1137*t1141+t924)*t610-(t148+(-t321*t614+t737)*t532)*t609+t614*t59;
        double t611 = x[7];
        double t589 = std::pow(t611, 2);
        double t612 = x[6];
        double t590 = std::pow(t612, 2);
        double t898 = t589+t590;
        double t756 = t1141*t415;
        double t482 = t534*t610;
        double t483 = t535*t609;
        double t1144 = t483+t482;
        double t1191 = -t240+t1144;
        double t1032 = t1141*t564;
        double t1118 = t588*t90;
        double t823 = t610*t1032;
        double t1190 = -t534*t823-(t1032*t535+t609*t90)*t609-t1118;
        double t673 = t1026*t321+t1141*t222;
        double t958 = t605*t615;
        double t973 = t602*t614;
        double t1066 = (t958-t973)*t615;
        double t726 = t594-t931;
        double t317 = t605*t726+t1066;
        double t293 = t614*t317;
        double t1189 = -(-t756+t293)*t534-t147;
        double t1073 = t518*t618;
        double t1134 = t1073-t487;
        double t977 = t602*t608;
        double t542 = t598*t977;
        double t970 = t603*t607;
        double t543 = t598*t970;
        double t901 = t543-t542;
        double t919 = t901*t1143;
        double t1188 = t1134*t1141+t919;
        double t1025 = t564*t618;
        double t356 = t518*t1025;
        double t425 = t596*t486;
        double t227 = t486*t595+t356+t425;
        double t519 = -t970+t977;
        double t663 = t230+t1152;
        double t906 = t1073+t486;
        double t1185 = -(-t1141*t906-t227-t919)*t609+t663*t519;
        double t414 = -t572+t1142;
        double t1184 = t414*t565;
        double t1181 = t1141*t353;
        double t1179 = (t353*t899+t231)*t519;
        double t599 = x[19];
        double t999 = t599*t601;
        double t510 = t598*t602-t999;
        double t1178 = t1142*t510;
        double t1007 = t598*t603;
        double t600 = x[18];
        double t993 = t600*t601;
        double t511 = -t993+t1007;
        double t1177 = t1142*t511;
        double t1064 = t1142*t591;
        double t176 = -t1026*t604+t1142*t592+t1184;
        double t1172 = t414*t898+t1064+t176;
        double t1021 = t565*t614;
        double t982 = t601*t614;
        double t1067 = (t604*t615-t982)*t615;
        double t316 = t604*t726+t1067;
        double t778 = t316*t1021;
        double t983 = t601*t613;
        double t1171 = -t1064*(t565+t898)+t778-(-t1143*t898-t499)*t983;
        double t1164 = t534*t603;
        double t926 = t616*t617;
        double t412 = t518*t926;
        double t413 = t596*t477;
        double t796 = t603*t955;
        double t488 = t564*t796;
        double t980 = t602*t603;
        double t815 = t564*t980;
        double t189 = (t412-t815)*t617+t413-t488;
        double t1170 = t1141*(-t477+t1164)-t189-t919;
        double t1162 = t535*t602;
        double t799 = t602*t950;
        double t489 = t564*t799;
        double t190 = t413+t489+(t412+t815)*t617;
        double t1169 = t1141*(t477+t1162)+t190+t919;
        double t638 = t519*t598;
        double t896 = t614-t618;
        double t758 = t531*t896;
        double t1163 = t535*t1141;
        double t998 = t599*t608;
        double t997 = t599*t610;
        double t757 = t608*t899;
        double t975 = t602*t612;
        double t1157 = t613*(t603*t611+t975);
        double t1048 = t532*t609;
        double t1054 = t531*t610;
        double t650 = -t1048+t1054;
        double t1156 = t617*t650;
        double t1154 = t898*t609;
        double t1151 = t899*t416;
        double t319 = -t601*t930+t604*t646;
        double t878 = t319*t1026;
        double t622 = (-t1143*t899+t499)*t983+t878-t899*t220;
        double t383 = t1143+t1141;
        double t662 = -t383*t898-t250;
        double t344 = t607*t352;
        double t480 = t534*t609;
        double t1149 = -t344+t480;
        double t840 = t518*t576;
        double t963 = t604*t603;
        double t986 = t601*t606;
        double t516 = -t963+t986;
        double t850 = t516*t972;
        double t275 = t840+t850;
        double t964 = t604*t602;
        double t987 = t601*t605;
        double t514 = -t964+t987;
        double t409 = t514*t926;
        double t382 = t603*t409;
        double t1148 = t275*t618-t382;
        double t513 = t598*t605-t599*t604;
        double t909 = t510*t617+t513*t618;
        double t515 = t598*t606-t600*t604;
        double t908 = t511*t617+t515*t618;
        double t474 = t615*t510;
        double t296 = t513*t616+t474;
        double t475 = t615*t511;
        double t297 = t515*t616+t475;
        double t512 = t599*t603-t600*t602;
        double t517 = t599*t606-t600*t605;
        double t1147 = t512*t615+t517*t616;
        double t540 = t598*t565;
        double t1146 = -t480+t540;
        double t1145 = t482-t483;
        double t399 = t534-t566;
        double t1140 = -t231+t1181;
        double t1044 = t532*t618;
        double t1070 = t519*t609;
        double t927 = t615*t618;
        double t536 = -t926+t927;
        double t809 = t534*t1006;
        double t708 = t564*t809;
        double t774 = t531*t1073;
        double t1030 = t564*t598;
        double t829 = t534*t1030;
        double t1075 = t518*t536;
        double t843 = t588*t1075;
        double t891 = -(t586*t774-t602*t708+t607*(t603*t829+t607*t774))*t610+t519*t843+t609*(-t564*t535*t638+(t1044*t1141+t1070*t536)*t518);  // NOLINT [whitespace/line_length]
        double t476 = t616*t516;
        double t394 = t605*t476;
        double t841 = t518*t575;
        double t277 = t394+t841;
        double t396 = t514*t952;
        double t279 = t396-t841;
        double t410 = t516*t926;
        double t381 = t602*t410;
        double t849 = t607*t476;
        double t1083 = t514*t616;
        double t859 = t608*t1083;
        double t1135 = t609*((t279*t618+t381)*t608-t535*t849)-(t534*t859-t607*(t277*t618+t382))*t610;
        double t1023 = t565*t601;
        double t1051 = t532*t601;
        double t423 = t596*t1051;
        double t798 = t603*t961;
        double t704 = t564*t798;
        double t863 = t514*t1025;
        double t739 = (t1023*t531+t863)*t608-(t332*t981+t423-t704)*t607;
        double t1102 = t353*t608;
        double t484 = t535*t610;
        double t267 = t484-t1102;
        double t1132 = t1141*t416;
        double t1129 = -t592*t90+t68*t899;
        double t1128 = (t633+t1132)*t535-t148;
        double t1065 = (t953-t967)*t615;
        double t318 = t606*t726+t1065;
        double t284 = t318*t955;
        double t928 = t615*t617;
        double t400 = t518*t928;
        double t1012 = t596*t615;
        double t401 = t518*t1012;
        double t115 = t401+t284+(t318*t602+t400)*t617;
        double t932 = t614*t615;
        double t408 = t518*t932;
        double t1127 = -t1141*(t416*t534-t408)+t614*t115+t923;
        double t787 = t317*t950;
        double t116 = t401-t787+(-t317*t603+t400)*t617;
        double t1126 = -t1141*(t415*t535+t408)-t614*t116+t924;
        double t629 = t1026*t320+t1141*t221;
        double t1060 = t1142*t603;
        double t427 = t596*t1060;
        double t943 = t608*t609;
        double t688 = t943*t1051;
        double t551 = t601*t950;
        double t697 = t564*t551;
        double t801 = t601*t955;
        double t701 = t564*t801;
        double t777 = t596*t1083;
        double t989 = t601*t602;
        double t817 = t564*t989;
        double t1124 = -t565*t688+((t777+t701+(t409+t817)*t617)*t608+(t330*t965+t427-t697)*t607)*t610;
        double t1043 = t533*t588;
        double t1123 = t591*(t533*t587+t1043)-t622;
        double t1031 = t1141*t603;
        double t1071 = t519*t588;
        double t428 = t596*t487;
        double t916 = t428-t356;
        double t228 = t487*t595+t916;
        double t711 = t565*t542;
        double t357 = t564*t711;
        double t813 = t565*t1007;
        double t387 = t564*t813;
        double t1122 = ((t228*t607-t387)*t607+t228*t586+t357)*t610-t609*(t1031*t532*t565+t1070*t231)-t231*t1071;
        double t1112 = t229*t609;
        double t1009 = t597*t614;
        double t814 = t565*t1009;
        double t687 = t317*t814;
        double t755 = t1141*t177;
        double t1121 = (t597*t755+t607*t878-t687)*t899+t1141*(t1112*t592+t687);
        double t1024 = t565*t597;
        double t1120 = (-t358+(t565+t1141)*t597*t531)*t899-t1141*(t1024*t531-t1112);
        double t1119 = t176*t898-t1171;
        double t1111 = t229*t610;
        double t1110 = t230*t592;
        double t1109 = t231*t598;
        double t1088 = t513*t614;
        double t1108 = (t909-t1088)*t608;
        double t1081 = t515*t614;
        double t1107 = (t908-t1081)*t607;
        double t1036 = t535*t607;
        double t479 = t534*t608;
        double t1106 = t650*(t479-t1036);
        double t1101 = t383*t598;
        double t539 = t565+t585;
        double t507 = t539+t586;
        double t1092 = t507*t532;
        double t1091 = t510*t608;
        double t1090 = t511*t607;
        double t1089 = t512*t531;
        double t462 = t512*t617;
        double t1087 = t514*t536;
        double t1086 = t514*t594;
        double t1085 = t514*t614;
        double t1084 = t514*t615;
        double t1082 = t514*t618;
        double t1080 = t516*t536;
        double t1079 = t516*t594;
        double t1078 = t516*t615;
        double t1077 = t516*t618;
        double t1076 = t517*t614;
        double t1072 = t519*t536;
        double t1069 = (-t954+t959)*t536;
        double t969 = t603*t609;
        double t976 = t602*t610;
        double t521 = t969+t976;
        double t1068 = t521*t613;
        double t1061 = t1142*t602;
        double t1059 = t1142*t608;
        double t1058 = t1142*t609;
        double t1057 = t531*t586;
        double t1055 = t531*t601;
        double t478 = t531*t607;
        double t1053 = t531*t611;
        double t348 = t532*t535;
        double t1052 = t532*t598;
        double t1050 = t532*t604;
        double t1049 = t532*t607;
        double t1047 = t532*t611;
        double t1046 = t532*t612;
        double t1042 = t533*t591;
        double t1041 = t533*t607;
        double t1040 = t533*t608;
        double t1039 = t534*t586;
        double t481 = t534*t588;
        double t1038 = t534*t607;
        double t1037 = t535*t586;
        double t485 = t535*t588;
        double t1035 = t535*t611;
        double t1034 = t536*t586;
        double t1033 = t536*t607;
        double t1027 = t564*t608;
        double t1022 = t565*t607;
        double t1020 = t586*t602;
        double t1019 = t588*t607;
        double t1018 = t591*t598;
        double t1017 = t594*t605;
        double t1015 = t596*t604;
        double t1011 = t597*t607;
        double t1010 = t597*t608;
        double t1008 = t598*t601;
        double t1005 = t598*t609;
        double t1004 = t598*t610;
        double t1003 = t598*t614;
        double t1002 = t598*t615;
        double t1001 = t598*t618;
        double t1000 = t599*t600;
        double t996 = t599*t614;
        double t995 = t599*t615;
        double t994 = t599*t616;
        double t992 = t600*t607;
        double t991 = t600*t611;
        double t990 = t600*t615;
        double t988 = t601*t603;
        double t985 = t601*t607;
        double t984 = t601*t608;
        double t978 = t602*t607;
        double t974 = t602*t613;
        double t968 = t603*t613;
        double t962 = t604*t617;
        double t949 = t607*t609;
        double t948 = t607*t610;
        double t946 = t607*t615;
        double t945 = t607*t618;
        double t944 = t608*t531;
        double t941 = t608*t611;
        double t939 = t608*t618;
        double t937 = t609*t176;
        double t936 = t609*t611;
        double t935 = t609*t612;
        double t934 = t609*t614;
        double t933 = t610*t611;
        double t472 = t614*t516;
        double t920 = t1037*t532+t348*t585;
        double t393 = t515*t957;
        double t544 = t600*t964;
        double t806 = t598*t960;
        double t446 = -t544+t806;
        double t918 = t446*t615+t393;
        double t395 = t513*t952;
        double t545 = t599*t963;
        double t807 = t598*t979;
        double t448 = -t545+t807;
        double t917 = t448*t615+t395;
        double t907 = t517*t618+t462;
        double t905 = t488-t425;
        double t904 = -t763-t587;
        double t538 = t587+t564;
        double t902 = -t538-t588;
        double t571 = t609*t590;
        double t900 = t589*t609+t571;
        double t897 = -t611+t609;
        double t214 = t598*t230;
        double t780 = t518*t1033;
        double t833 = t531*t1038;
        double t834 = t531*t1039;
        double t26 = t588*t780+(-t834-t607*(t214+t833))*t610+t609*(t609*t780-t708+(t1039+t607*(t540+t1038))*t532);
        double t889 = t229*t948;
        double t888 = t229*t943;
        double t887 = t230*t1004;
        double t886 = t230*t1003;
        double t885 = t230*t949;
        double t884 = t1141*t1109;
        double t883 = t250*t1005;
        double t882 = t1145*t607*t1142;
        double t881 = t316*t950;
        double t880 = t565*t293;
        double t879 = t318*t1021;
        double t876 = t321*t1030;
        double t875 = t598*t344;
        double t874 = t353*t1006;
        double t873 = t763*t480;
        double t872 = t763*t484;
        double t871 = t507*t1054;
        double t870 = t507*t1048;
        double t869 = t598*t462;
        double t868 = t512*t981;
        double t867 = t513*t574;
        double t866 = t602*t1087;
        double t865 = t605*t1087;
        double t864 = t514*t1033;
        double t862 = t598*t1083;
        double t861 = t514*t994;
        double t860 = t608*t1084;
        double t858 = t514*t608*t617;
        double t857 = t514*t939;
        double t856 = t514*t932;
        double t855 = t514*t931;
        double t854 = t515*t573;
        double t853 = t516*t1025;
        double t852 = t598*t476;
        double t851 = t600*t476;
        double t848 = t516*t607*t617;
        double t847 = t516*t945;
        double t846 = t516*t932;
        double t845 = t614*t476;
        double t844 = t517*t572;
        double t842 = t586*t469;
        double t838 = t518*t930;
        double t837 = t1142*t948;
        double t836 = t1142*t943;
        double t835 = t531*t1041;
        double t832 = t565*t944;
        double t831 = t598*t944;
        double t830 = t533*t949;
        double t828 = t534*t985;
        double t827 = t535*t946;
        double t826 = t1141*t899*t565;
        double t825 = t1141*t499;
        double t824 = t1141*t1030;
        double t822 = t1141*t540;
        double t821 = t1141*t1004;
        double t820 = t1141*t976;
        double t819 = t564*t1023;
        double t818 = t564*t1006;
        double t816 = t564*t988;
        double t812 = t565*t1006;
        double t811 = t565*t1004;
        double t810 = t250*t1004;
        double t808 = t598*t984;
        double t805 = t598*t940;
        double t804 = t598*t932;
        double t803 = t599*t986;
        double t802 = t600*t987;
        double t800 = t602*t961;
        double t795 = t614*t576;
        double t794 = t604*t929;
        double t793 = t316*t955;
        double t792 = t604*t1017;
        double t790 = t1142*t573;
        double t788 = t605*t929;
        double t786 = t604*t1016;
        double t784 = t1142*t574;
        double t783 = t606*t931;
        double t781 = t230*t1019;
        double t779 = t608*t1075;
        double t776 = t534*t477;
        double t775 = t535*t477;
        double t550 = t604*t955;
        double t552 = t604*t950;
        double t553 = t605*t950;
        double t773 = t190*t586-t357+t607*(t190*t607+t387);
        double t206 = t231*t1006;
        double t772 = t206-t920;
        double t771 = t293-t1131;
        double t770 = t318*t614-t1133;
        double t201 = t331*t965+t428-t489;
        double t765 = t1143+t898;
        double t764 = -t477-t901;
        double t215 = t592*t231;
        double t760 = t161*t899-t215;
        double t759 = t1142*t896;
        double t753 = t1141*t532;
        double t748 = t564*t604;
        double t744 = t896*t532;
        double t343 = t353*t587;
        double t742 = t353*t588-t231+t343;
        double t740 = t189*t586-t357;
        double t391 = t514*t966;
        double t276 = -t391+t840;
        double t218 = -t276*t618-t381;
        double t388 = t512*t575;
        double t445 = -t545+t802;
        double t736 = t445*t616-t388;
        double t447 = -t544+t803;
        double t735 = t447*t616+t388;
        double t389 = t602*t475;
        double t734 = t446*t616+t389;
        double t449 = -t802+t807;
        double t733 = t449*t616+t389;
        double t390 = t603*t474;
        double t732 = t448*t616+t390;
        double t450 = t803-t806;
        double t731 = -t450*t616+t390;
        double t392 = t517*t576;
        double t730 = -t445*t615+t392;
        double t729 = t447*t615+t392;
        double t728 = t449*t615+t393;
        double t727 = -t450*t615+t395;
        double t725 = t231*t805;
        double t724 = t1143*t598*t985;
        double t723 = t1143*t808;
        double t722 = t514*t606*t927;
        double t721 = t516*t605*t927;
        double t720 = (-t933+t935)*t1142*t1042;
        double t719 = t612*t830;
        double t718 = t933*t1040;
        double t717 = t612*t826;
        double t716 = t597*t825;
        double t715 = t613*t825;
        double t714 = t564*t821;
        double t713 = t1141*t811;
        double t712 = t565*t820;
        double t710 = t230*t821;
        double t709 = t320*t818;
        double t346 = t531*t812;
        double t707 = t601*t412;
        double t706 = t564*t800;
        double t705 = t602*t795;
        double t703 = t603*t795;
        double t702 = t604*t518*t927;
        double t699 = t564*t550;
        double t698 = t1142*t788;
        double t696 = t564*t552;
        double t695 = t564*t553;
        double t694 = t1142*t782;
        double t692 = t610*t864;
        double t691 = t948*t1055;
        double t689 = t943*t1080;
        double t686 = t534*t408;
        double t685 = t535*t408;
        double t682 = t1141*t1080;
        double t681 = t231*t757;
        double t678 = t396+t850;
        double t670 = t602*t716;
        double t669 = t603*t716;
        double t666 = t770-t1132;
        double t658 = (t478-t1030)*t607+t1057;
        double t664 = -(t535*t658+t346)*t610+(t609*t779+t206+t920)*t609;
        double t83 = (t1176*t601-t706)*t608-(t1051*t595+t423+t853)*t607;
        double t657 = t351*t610-t1040;
        double t653 = -t594*t604-t1067;
        double t652 = t1017+t1066;
        double t651 = t1016+t1065;
        double t649 = -t531*t612-t1047;
        double t136 = -t532*t847+(t618*t678-t707)*t608;
        double t33 = -t1102*t588+t608*(-t343+t161)-t126*t610;
        double t384 = t1143-t899;
        double t433 = t1142*t956;
        double t434 = t1142*t1014;
        double t546 = t604*t958;
        double t101 = (-t601*t747+t433-t855)*t617-(t792+(t546-t1085)*t615)*t618+t434;
        double t636 = -t160*t899+t1110;
        double t635 = t763*t899-t1032;
        double t19 = -t1071*t230-(t230*t1070+t227*t586-t357+(t227*t607+t387)*t607)*t609+t531*t712;
        double t278 = t394+t391;
        double t133 = -(t278*t618+t707)*t607+t531*t857;
        double t632 = (t531*t813+t607*t776)*t607+t586*t776;
        double t630 = t1141*t178-t879;
        double t628 = -t1090-t477+t1091;
        double t627 = t535-t1006;
        double t411 = t596*t476;
        double t424 = t596*t1061;
        double t62 = -t565*t691+((t330*t971+t424-t701)*t608+t607*(t411+t697+(t410+t816)*t617))*t609;
        double t95 = -(t532*t813+t607*t775)*t607-t586*t775+t532*t711;
        double t347 = t531*t534;
        double t81 = (-t347+t875)*t610-(t809-t532*(t534+t566))*t609+t780;
        double t398 = -t531+t566;
        double t625 = -(-t348+t874)*t609-(-t398*t535+t831)*t610+t779;
        double t624 = t177+t756;
        double t623 = t178+t1132;
        double t584 = std::pow(t603, 2);
        double t583 = std::pow(t602, 2);
        double t582 = std::pow(t600, 2);
        double t581 = std::pow(t599, 2);
        double t580 = std::pow(t598, 2);
        double t579 = std::pow(t597, 2);
        double t578 = t597*t579;
        double t577 = std::pow(t579, 2);
        double t548 = t604*t953;
        double t547 = t601*t952;
        double t529 = t574-t966;
        double t527 = t573-t972;
        double t525 = t572-t575;
        double t522 = t943+t948;
        double t505 = -t584*t599+t600*t980;
        double t504 = -t584*t598+t600*t988;
        double t503 = -t583*t598+t599*t989;
        double t440 = t1142*t1013;
        double t439 = t1142*t951;
        double t432 = t532*t1015;
        double t431 = t532*t962;
        double t430 = t531*t1015;
        double t429 = t531*t962;
        double t419 = t952+t1159;
        double t418 = t957+t1160;
        double t417 = t1142-t981;
        double t407 = t516*t1012;
        double t406 = t516*t928;
        double t405 = t514*t1012;
        double t404 = t514*t928;
        double t403 = t586*t1060;
        double t402 = t1142*t1020;
        double t397 = t531+t566;
        double t378 = t507-t899;
        double t355 = t414*t950;
        double t354 = t414*t955;
        double t350 = t422*t566;
        double t349 = t421*t1006;
        double t329 = t540-t1038;
        double t327 = t478+t1030;
        double t326 = t384+t1141;
        double t309 = t944-t1049;
        double t285 = t317*t961;
        double t283 = t318*t961;
        double t273 = t318*t540;
        double t270 = t329*t1049;
        double t268 = t317*t812;
        double t264 = t1147-t1076;
        double t262 = t297-t1081;
        double t260 = t296-t1088;
        double t217 = -t277*t617+t722;
        double t216 = -t279*t617+t721;
        double t210 = t907-t1147;
        double t209 = t908-t297;
        double t208 = t909-t296;
        double t199 = t332*t971-t905;
        double t194 = t1059*t565-t1202*t535;
        double t192 = t442-t695+(-t564*t960+t441)*t617;
        double t191 = t436-t695+(-t564*t979+t435)*t617;
        double t170 = t793+(t601*t652-t705)*t617;
        double t169 = t285-(t601*t789+t602*t653)*t617;
        double t168 = t283-(t601*t783+t603*t653)*t617;
        double t167 = t881+(t601*t651-t703)*t617;
        double t163 = t440-t532*t794+(-t532*t982+t439)*t617;
        double t162 = t434-t531*t794+(-t531*t982+t433)*t617;
        double t155 = t1209*t603-t608*t819;
        double t152 = -t1022*t1142+t1202*t534;
        double t146 = t1209*t602-t607*t819;
        double t143 = t214-t833;
        double t142 = t189*t607+t387;
        double t139 = -t270+(t532*t586+t818)*t534;
        double t138 = t346+(-t327*t607-t1057)*t535;
        double t134 = (-t278*t617+t702)*t608+t532*t848;
        double t132 = t531*t858-t607*(t617*t678+t702);
        double t128 = -t383*t899+t250;
        double t102 = t440+(-t786-(-t472+t548)*t615)*t618+(-t601*t746+t439-t845)*t617;
        double t94 = t346*t602-t632;
        double t91 = -t143*t607+t834;
        double t75 = t1163*t318+t725;
        double t74 = -t415*t899+t624;
        double t72 = t317*t1039+t607*(t1038*t317+t886);
        double t71 = t535*t823-t681;
        double t69 = t781+(t1032*t534+t885)*t609;
        double t60 = -t321*t823-t614*t681;
        double t57 = t614*t781+(-t1032*t320+t614*t885)*t609;
        double t51 = (t430-t699+(-t564*t987+t429)*t617)*t608-t607*(t432+(-t516*t593-t1079+t431)*t617+(-t1045*t601-t532*t961)*t614);  // NOLINT [whitespace/line_length]
        double t50 = (t430+(-t514*t593-t1086+t429)*t617+(-t1055*t617-t531*t961)*t614)*t608-(t432-t696+(-t564*t986+t431)*t617)*t607;  // NOLINT [whitespace/line_length]
        double t48 = t142*t607+t740;
        double t47 = t535*t842+(t436-t693+(-t531*t967+t435)*t617)*t1006-t607*(t1052*t321-t518*t827);
        double t46 = t534*t842+t320*t831-t607*(-t518*t534*t946+t598*(t442-t532*t788+(-t532*t973+t441)*t617));
        double t42 = t155*t899-t669;
        double t41 = (-t534*t860+t607*(t275*t617+t722))*t610+t609*((-t276*t617+t721)*t608-t516*t827);
        double t39 = t597*t673-t608*t778;
        double t38 = t597*t629-t607*t778;
        double t35 = (t1024*t532-t1111)*t1141-t899*(-t1027*t533+t1092*t597);
        double t29 = ((t405-t793+(-t316*t602+t404)*t617)*t608-t607*(t440-t696+(-t564*t963+t439)*t617))*t610+(t432-t694+(-t1142*t967+t431)*t617)*t943;  // NOLINT [whitespace/line_length]
        double t28 = -(t430-t698+(-t1142*t973+t429)*t617)*t948+t609*((t434-t699+(-t564*t964+t433)*t617)*t608-t607*(t407-t881+(-t316*t603+t406)*t617));  // NOLINT [whitespace/line_length]
        double t27 = t588*t779+t664;
        double t22 = t116*t586+t268+t607*(t116*t607-t876);
        double t21 = t115*t586+t709+(t115*t607-t273)*t607;
        double t18 = t614*t1118-t320*t823+t609*(-t1032*t321+t90*t934);
        double t12 = t1141*t770-t666*t899+t879;
        double t11 = t880-t899*(-t756+t771)+t1141*t771;
        double t8 = ((t317*t479-t607*(t787+(t602*t651-t700)*t617))*t610+t609*(-(t284-t617*(t602*t783-t603*t652))*t608+t318*t1036))*t598+(-t1069*t899+t1141*t1156)*t518;  // NOLINT [whitespace/line_length]
        double t6 = (-t1111*t592-t318*t814)*t1141-t899*(t597*t630+t608*t878);
        double t5 = t58*t588+(-t320*t1057-(t320*t478-t886)*t607)*t610+t609*(t58*t609+t106*t586-t709+t607*(t106*t607+t273));  // NOLINT [whitespace/line_length]
        double t4 = t59*t588+(-t103*t586-t268-t607*(t103*t607-t876))*t610+t609*(t321*t753+t59*t609-t725);

        coeffs <<-t577*(((-t69*t590+(-t531*t826-t611*t69)*t611)*t591-t602*(-t571+t611*(-t899-t936))*t715-t614*(t57*t590+t611*(-t317*t826+t57*t611)))*t582+((-t1190*t590+t531*t717-(t1190*t611-t532*t826)*t611)*t591+(-t521*t590-t899*t975+t611*(-t521*t611-t603*t899))*t715+t614*(t18*t590-t317*t717+(t18*t611-t318*t826)*t611))*t1000-((t532*t717+t71*t898)*t591+(-t610*t898-t612*t899)*t603*t715+t614*(-t318*t717+t60*t898))*t581),  // NOLINT [whitespace/line_length]
        (((t1120*t611+t597*t69+t898*(t1011*t663+t152*t609))*t591+(-t897*t670+(-t611*t899-t1154)*t146)*t613+t1121*t611+t57*t1009+t898*(t1011*t636-t38*t609))*t582+(((-t1120*t612+t35*t611+t597*t1190+t898*(t1201*t597-t152*t610+t194*t609))*t591+((t146*t899-t670)*t612+t42*t611+t521*t716+t898*(t146*t610+t155*t609))*t613-t1121*t612+t6*t611-t18*t1009+t898*(t1129*t597+t38*t610+t39*t609))*t599+(t26*t590+(t91*t588-(t534*t824-t609*t91)*t609)*t612+t611*(t26*t611+t138*t588+t531*t713+t609*(t138*t609-t884)))*t591+(t19*t590+(t1146*t609-t481)*t975*t1032+(-t1030*t712+t19*t611-t773*t899)*t611)*t613-t614*(t5*t590+(t72*t588+t609*(-t320*t824+t609*t72))*t612+(t5*t611+t22*t588+t317*t713+t609*(t22*t609+t614*t884))*t611))*t600-(((t35*t612-t597*t71+t898*(t1010*t742+t194*t610))*t591+(t42*t612+(t155*t898+t669)*t610)*t613+t6*t612-t60*t1009+t898*(t1010*t760+t39*t610))*t599+(t27*t590+(t139*t588+t710-t609*(-t139*t609+t532*t822))*t612+t611*(t27*t611+t535*t714-t772*t899))*t591+(t1122*t590+(t48*t588+(t1141*t387+t48*t609)*t609)*t612-t611*(-t1122*t611+(-t535*t587-t485+t811)*t1141*t1029))*t613-t614*(t4*t590+(t21*t588-t614*t710+(t21*t609-t318*t822)*t609)*t612+(t321*t714+t4*t611+t75*t899)*t611))*t599)*t578,  // NOLINT [whitespace/line_length]
        ((-((t1144*t763+t1191*t898+t378*t649-t1201)*t591+t11*t612+t12*t611-t629*t610-t673*t609+t898*t1206+t128*t1157+t662*t1068-t1129)*t1000+((-t1053*t378+t1149*t898+t873)*t591+(t128*t611+t609*t662)*t974+t11*t611-t629*t609+(t591*t663+t636)*t607+t898*t1199)*t582+((-t1046*t378+t267*t898+t872)*t591+(t128*t612+t610*t662)*t968+t12*t612-t673*t610+(-t591*t742-t760)*t608+t898*t33)*t581)*t579+(((t81*t898+t26)*t591+(t19+t898*(t602*t871-t1185))*t613-t614*t5+t898*t1198)*t600+t599*((-t625*t898-t664)*t591+(-t1122+t898*((t1188-t228)*t610+t603*t870+t1179))*t613+t614*t4+t898*t1196)+((-t1205*t599+(-t598*t873+t834-t899*(t347+t875))*t591+(t534*t635+t883)*t974+t629*t1005-t614*t72-t899*t1189)*t600+t599*((-t1039*t532+t598*t870+t270-t899*(-t399*t532-t809))*t591+(-t1170*t899-t603*t883-t740)*t613+t630*t1005+t614*t21-t899*t1127)+(t1208*t599+(-t109*t591-t1141*t160-t1110)*t598)*t997)*t612+(t1205*t582+((t1140*t1005-t1037*t531+t598*t871+t346-t899*(-t397*t535+t831))*t591+(-t1169*t899-t602*t810+t773)*t613+(t755-t880)*t1004-(t1141*t161+t215)*t1005-t614*t22-t1208*t997-t899*t1126)*t600+t599*((-t598*t872-t899*(t348+t874)-t772)*t591+(t535*t635+t810)*t968+t673*t1004-t75*t614-t899*t1128))*t611+(-t591*t843+(-t176*t590+(t316*t614-t1184)*t589+t1171)*t997+(t1042*t902+t622)*t991+(((t533*t538+t1043)*t591-t622)*t599-t591*t829)*t612)*t998+((-t327*t535*t591*t600+t1123*t582)*t611+((-t1123*t599-t143*t591)*t600-t599*t142*t613)*t612-t1119*t582*t609)*t607+t1119*t522*t1000)*t597+((t590*t882+(-t588*t835+(-t531*t830+t1142*(t329*t607-t1039))*t609)*t612-(-t611*t882+t609*(t1006*t229-t1034*t516-t1080*t585)+(t1049*t899-t610*t658)*t533)*t611)*t591+(t62*t590-(-t588*t828+(-t533*t1020+(t1146*t601-t533*t978)*t607)*t609)*t564*t612+t611*(t62*t611+t83*t588-(t402+t607*(-t1008*t564+t1142*t978))*t565*t610+(t1031*t229+t609*t83)*t609))*t613+t614*(t28*t590+(t169*t1019-t609*(-t170*t586+(-t169*t609-t170*t607+t316*t540)*t607))*t612+(t28*t611+t51*t588+(-t162*t586-t607*(-t1030*t319+t162*t607))*t610+t609*(t102*t1141-t229*t805+t51*t609))*t611))*t600+(((t514*t1034-t607*(t229*t598-t864))*t610*t612+(t753*t935+((t531*t588+(t531*t609-t1030)*t609)*t612+t899*t1047)*t608)*t533+(-t1145*t608*t590+(-t1145*t941+(t812-t1163)*t610)*t611)*t1142)*t591+(t1124*t590+(-t739*t588+t229*t820-t609*(t739*t609+(t1060*t585-t564*t808+t403)*t565))*t612-t611*(-t1124*t611+((-t1031*t533+t565*t808)*t610+t601*t535*t757)*t564))*t613-t614*(t29*t590+(t50*t588+(-t101*t586-t607*(-t1003*t229+t101*t607))*t610+t609*(t1141*t163-t319*t818+t50*t609))*t612+(t29*t611+(-t1141*t167+t316*t812)*t610+t168*t757)*t611))*t599+(t590*t1106+(-t309*t481+t607*t887-(t309*t480+t598*(t1022*t532-t1027*t534))*t609)*t612+(t611*t1106+t309*t485-(-t1028*t535+t832)*t1004+t609*(t309*t483+t206))*t611)*t1018+(-t891*t590+(t94*t588-t602*t710+t609*(t94*t609+t598*(t199*t586-t357+t607*(t199*t607+t387))))*t612-(t891*t611+t95*t588-(t201*t586+t357+t607*(t201*t607-t387))*t1004+(t603*t884+t609*t95)*t609)*t611)*t613-t614*(t8*t590+(-t46*t588+(t586*t605+t607*(t605*t607-t1003))*t887-(t46*t609+t598*(t191*t586-t709+t607*(t191*t607+t273)))*t609)*t612+t611*(t8*t611+t47*t588-(t192*t586+t268+(t192*t607-t876)*t607)*t1004+(t47*t609+(t1141*t606-t805)*t1109)*t609)))*t579,  // NOLINT [whitespace/line_length]
        -((((t1053-t1149)*t591+(-t326*t611+t383*t609+t900)*t974+t74*t611+t898*(-t421*t609+t344)-t1199)*t582-((-t649-t1191)*t591+t74*t612+(t623-t1151)*t611-t326*t1157+(t383+t898)*t1068+t898*(-t421*t610-t422*t609+t240)-t1206)*t1000-((t267-t1046)*t591+(t326*t612+(-t765-t1141)*t610)*t968+(t666+t1151)*t612+t33+t898*(t422*t610-t1102))*t581)*t578+(((t351*t591-t85)*t936+((-t533*t611+t1058)*t591-(t1143*t609-t384*t611+t900)*t983-t741*t611+t937+t414*t1154)*t607)*t582+((((-t351*t609+t1041)*t612-t657*t611)*t591+(t607*t741+t609*t85)*t612+t1207*t611+(-t607*t612-t941)*t384*t983+(t765*t983-t1172)*t522)*t599+(t397*t1035-t347*t612+(t1149*t612+(-t353*t609-t1054-t944)*t611)*t598-t81)*t591+((-t899*(t764-t1162)-t1169)*t611+(-(t1005*t383+t534*t904-t481)*t612+(t1101*t611-t507*t531)*t610)*t602+t898*((-t583*t615-t602*t957)*t610+(t901+t906)*t609-t519*t352)+t1185)*t613+(-t124*t1005-t899*(-t415*t534-t875)-t1189)*t612+((-t624*t610+(t161-t1181)*t609)*t598-t899*(t415*t627+t350+t408)-t1126)*t611+t898*((t421*t531-t875)*t610+(-t416*t566-t1136+t349)*t609-t111)-t1198)*t600+(((t591*t657-t1207)*t612+((t384*t612-t610*t765)*t983+t1172*t610)*t608)*t599+(t399*t1046-t348*t611+((-t352*t610-t1048+t479)*t612+t267*t611)*t598+t625)*t591+((-t899*(-t764-t1164)-t1170)*t612+(-t1188+t916)*t610+(-(t1004*t383+t535*t904-t485)*t611+t595*t1054+(t1101*t612-t1092)*t609)*t603+t898*((-t901-t1134)*t610+(-t584*t615-t603*t952)*t609+t519*t353)-t1179)*t613+(((t160-t1150)*t610-t623*t609)*t598-t899*(t399*t416+t349-t408)-t1127)*t612+(-t126*t1004-t899*(-t416*t535-t874)-t1128)*t611+t898*((-t1006*t415-t1137+t350)*t610+(t422*t532-t874)*t609+t112)-t1196)*t599)*t579+((((t1058*t399+t835)*t612+((-t1006*t351-t1080)*t609+(t398*t610-t1049)*t533)*t611-t882)*t591+(((-t1203*t602+t724)*t609+(-t564+t899)*t828)*t612+((t1061*t539+t402-t724)*t610-t108*t969-t899*((-t418*t601+t800)*t608+t607*(t1051+t1077))-t83)*t611-t62+t898*(t691+((-t417*t602+t801)*t608-t607*(t601*t965+t476+t551))*t609))*t613+((-t170*t614+t1141*(t354-(t527*t601-t576*t602)*t617))*t609+(-t598*t937-t614*t169-t899*(t415*t961+(-t525*t602+t601*t957)*t617))*t607)*t612+((-t220*t566+t162*t614+t1141*(-t531*t533+t790))*t610+(t158*t1006-t102*t614+t1141*((t574*t604-t1078)*t618+(t574*t601+t476)*t617-t784))*t609-t614*t51-t899*((-t531*t604+t601*t956+t550)*t940-t607*((t575*t603+t472+t547)*t617-t896*t1050)))*t611-t614*t28+t898*(-(t1142*t534-t531*t572)*t948+t609*((-t1142*t605+t602*t962+t550)*t940+t607*(t414*t965+t355+t846))))*t600+((((-t351*t566-t1087)*t610+((-t532+t1006)*t609-t944)*t533)*t612+(t1142*t610*t627+t1040*t532)*t611+t1145*t1059)*t591+((-t108*t976+(t1060*t539+t403-t723)*t609-t899*((t1055+t1082)*t608+(-t419*t601+t798)*t607)+t739)*t612+((-t1203*t603+t723)*t610+t902*t535*t984)*t611-t1124+t898*(((-t534*t601-t1083)*t608+(-t417*t603+t551)*t607)*t610+t688))*t613+((t158*t566-t101*t614+t1141*((t573*t604-t1084)*t618+(t573*t601+t1083)*t617-t790))*t610+(-t220*t1006+t163*t614+t1141*(-t532*t533+t784))*t609+t614*t50-t899*(((-t1055-t1085)*t617+t604*t758)*t608+(t601*t951-t1050+t552)*t947))*t612+((-t167*t614+t1141*(t355-(t529*t601-t576*t603)*t617))*t610+(-t176*t1004-t168*t614-t899*(t416*t961+(-t525*t603+t547)*t617))*t608)*t611+t614*t29+t898*(((t414*t971+t354+t856)*t608+(-t1142*t606+t603*t962+t552)*t947)*t610-(t1142*t535-t532*t572)*t943))*t599+(-t1106+(t534*t612-t1035)*t309+((t352*t948-(t479-t1049)*t609)*t612+(-(-t944+t1036)*t610+t353*t943)*t611)*t598)*t1018+((-t899*(-t531*t638+t776)+((t1141*(-t419*t602+t796)+t905-t919)*t609+(t109*t610-t315*t609-t832)*t602)*t598+t632)*t612+((-(t1141*(t418*t603-t799)+t201-t919)*t610-t1140*t969)*t598-t899*(t1052*t519-t775)+t95)*t611+t891+t898*((-t618*t650+t1072)*t518+t1145*t638))*t613+(-t899*(-t686+(t421*t944-(t532*t971+t437-t444)*t607)*t598)+(-t147*t610+t609*t923)*t598+((-t191+t1141*(-t532*t605+t602*t951+t553))*t1005-t46+t109*t605*t1004)*t614)*t612+(((-t192+t1141*(-t531*t606+t603*t956+t553))*t614+t924)*t1004-(t1140*t574+t148)*t1005+t614*t47-t899*(t685+(-(t531*t965-t438+t443)*t608+t422*t1049)*t598))*t611+t614*t8+t898*((t1069+t1156)*t473+((-t415*t479+(t415*t950-(t529*t602-t603*t957)*t617)*t607)*t610+((t416*t955+(-t527*t603+t602*t952)*t617)*t608-t416*t1036)*t609)*t598))*t597+(-t607*t720+(t564*t719+(-t565*t837+t888)*t611)*t983+t614*(t316*t719+(-t319*t837+t604*t888)*t611))*t600+(t608*t720-((t565*t836-t889)*t612-t564*t718)*t983-t614*((t319*t836-t604*t889)*t612-t316*t718))*t599+((t692-t609*(-t1041*t532+t1142*t479))*t612+((-t1036*t1142+t533*t944)*t610+t689)*t611)*t1018+(t1135*t590+(t133*t588+(-t586*t866-(t1008*t230+t602*t864)*t607)*t610+t609*(t133*t609+t218*t586+(t1061*t595+t424-t863)*t1006+(t218*t607+t598*(t411-t704+t617*(t410-t816)))*t607))*t612-(-t1135*t611+t136*t588+(-t1148*t586-(t777-t706+(t409-t817)*t617)*t1006-t607*(t1148*t607+t598*(t1060*t595+t427-t853)))*t610+(t136*t609+t231*t808+t603*t682)*t609)*t611)*t613-t614*(t41*t590+(t132*t588+(t586*t865+t607*(t605*t864+t598*(t430+(-t792-(t546+t1085)*t615)*t618+(-t602*t748+t429+t855)*t617)))*t610-(-t132*t609+t216*t586+(t434-t698+(t433+t1086-t705+(-t601*t973+t1084)*t615)*t617)*t1006-t607*(-t216*t607+t598*(t407+t283+(t318*t601+t406)*t617)))*t609)*t612+t611*(t41*t611+t134*t588+(-t217*t586+(t405+t285+t617*(t317*t601+t404))*t1006-(t217*t607+t598*(t440-t694+(t439+t1079-t703+(-t601*t967+t1078)*t615)*t617))*t607)*t610+(t134*t609+(t432+(-t786-(t472+t548)*t615)*t618+t617*(-t603*t748+t431+t845))*t1006+t606*t682)*t609)))*t597,  // NOLINT [whitespace/line_length]
        (t512*(t897*t600+(t612-t610)*t599)*t613-t210*(-t992+t998)+(-t600*t609+t997)*(t907-t1076)+(-t599*t612+t991)*t264)*t577+(((t628*t599+(t534-t1005)*t512)*t612+(t505*t617+(-t950+t1004)*t512-t628*t600)*t611+(t581*t984-t1089+(t1073-t542+t1090)*t599)*t610+(t582*t985+t512*t532+(-t1073-t543+t1091)*t600)*t609-t210*t519)*t613+((t1147*t607-t517*t934)*t598-t534*t264+(t208*t610+t297*t609-t1108+(-t515*t607+t469)*t614)*t599+(t607-t609)*t600*t909)*t612+((t1076*t610-t1147*t608)*t598+t535*t264+(t296*t610+t209*t609-t1107+(-t513*t608-t469)*t614)*t600+(t608-t610)*t599*t908)*t611+(t260*t998+(-t296*t600+t515*t996-t869)*t607+(-t473*t599+t1089)*t617+(-t598*t945-t758)*t517)*t610+((t1088*t600-t475*t599-t515*t994+t869)*t608+t262*t992+(t473*t600+t505*t615-t512*t952)*t617+(t598*t939+t744)*t517)*t609+((-t517*t573-t518*t995)*t618+(t477*t599-t512*t573)*t617+t1147*t573)*t608+((t517*t574+t518*t990)*t618+(-t477*t600+t512*t574)*t617-t1147*t574)*t607)*t578+(((t208*t976+((-t580*t602+t598*t999)*t608-t446*t618+(-t971+t566)*t511+t732)*t609+(-t1082*t599+t510*t531)*t608+(t445*t618-t731-t868)*t607-t776)*t612+(((-t580*t603+t598*t993)*t607-t448*t618+(-t965+t1006)*t510+t734)*t610+t209*t969+(t447*t618-t733+t868)*t608+(-t1077*t600+t511*t532)*t607+t775)*t611+((t503*t617-t510*t955+t861)*t608+(-t450*t618+t510*t965-t736)*t607+t774)*t610+((t449*t618+t511*t971-t735)*t608+(t504*t617-t511*t950+t851)*t607-t518*t1044)*t609-t518*t1072)*t613+((-t208*t566+(t513*t573+t514*t995)*t618+(t510*t573-t861)*t617-t296*t573)*t610+((t730-t854)*t618+(-t449*t614+t735)*t617+t614*t727+(-t262*t607+t1108)*t598)*t609+((t503*t615-t510*t957+t514*t996)*t617+t513*t758)*t608+((-t844+t918)*t618+(-t447*t614+t733)*t617-t614*t917)*t607+t686)*t612+(((-t729-t867)*t618+(t450*t614+t736)*t617+t614*t728+(-t260*t608+t1107)*t598)*t610+(-t209*t1006+(t515*t574+t516*t990)*t618+(t511*t574-t851)*t617-t297*t574)*t609+((t844+t917)*t618+(-t445*t614+t731)*t617-t614*t918)*t608+((t472*t600+t504*t615-t511*t952)*t617+t515*t744)*t607-t685)*t611+((t260*t534-t599*t856)*t608+((-t728+t867)*t618+(t448*t614-t734)*t617+t614*t729)*t607-t531*t838)*t610+(((-t727+t854)*t618+(t446*t614-t732)*t617-t614*t730)*t608+(t262*t535-t600*t846)*t607+t532*t838)*t609-t473*t1069)*t579+((((-t208*t985+t866)*t610+((t1001*t514-t1178)*t608+(t511*t533-t852)*t607-t218)*t609-t133)*t612+(((t510*t533-t862)*t608+(t1001*t516-t1177)*t607-t1148)*t610+(t1080*t603-t209*t984)*t609+t136)*t611-t1135)*t613+((((-t1002*t514-t513*t572)*t618+(-t510*t572+t862)*t617+t296*t572)*t607+t614*t865)*t610+(((-t1003*t514+t1178)*t617-t513*t759)*t608+(-t262*t533+t516*t804)*t607-t614*t216)*t609+t614*t132)*t612+(((-t260*t533+t514*t804)*t608+((-t472*t598+t1177)*t617-t515*t759)*t607-t614*t217)*t610+(((-t1002*t516-t515*t572)*t618+(-t511*t572+t852)*t617+t297*t572)*t608+t606*t536*t472)*t609+t614*t134)*t611+t41*t614)*t597+((-t601*t692+t609*(t1142*t857-t533*t849))*t612-((-t1142*t847+t533*t859)*t610+t601*t689)*t611)*t613-t614*((t604*t692-t609*(t516*t533*t946-t1142*t858))*t612+t611*((t1142*t848-t533*t860)*t610+t604*t689));  // NOLINT [whitespace/line_length]
        return coeffs;
    }

    inline Eigen::Matrix3d create_M(const double r, const Eigen::Matrix<double, 22, 1> &x) {
        Eigen::Matrix3d M;
        double t1 = std::pow(x[6], 2);
        double t2 = std::pow(x[7], 2);
        double t5 = 1+r*(t1+t2);
        double t7 = std::pow(x[0], 2);
        double t8 = std::pow(x[1], 2);
        double t9 = t7+t8;
        double t11 = t9*r+1;
        double t12 = std::pow(x[21], 2);
        double t13 = t12*t11;
        double t15 = x[0]*x[13];
        double t16 = x[1]*x[16];
        double t17 = t15+t16;
        double t19 = t9*x[20];
        double t29 = x[0]*x[14]+x[1]*x[17];
        double t32 = std::pow(x[8], 2);
        double t33 = std::pow(x[9], 2);
        double t36 = 1+r*(t32+t33);
        double t38 = std::pow(x[2], 2);
        double t39 = std::pow(x[3], 2);
        double t40 = t38+t39;
        double t42 = t40*r+1;
        double t43 = t12*t42;
        double t45 = x[2]*x[13];
        double t46 = x[3]*x[16];
        double t47 = t45+t46;
        double t49 = t40*x[20];
        double t59 = x[2]*x[14]+x[3]*x[17];
        double t62 = std::pow(x[10], 2);
        double t63 = std::pow(x[11], 2);
        double t66 = 1+r*(t62+t63);
        double t68 = std::pow(x[4], 2);
        double t69 = std::pow(x[5], 2);
        double t70 = t68+t69;
        double t73 = t12*(t70*r+1);
        double t75 = x[4]*x[13];
        double t76 = x[5]*x[16];
        double t77 = t75+t76;
        double t79 = t70*x[20];
        double t89 = x[4]*x[14]+x[5]*x[17];
        double t94 = x[0]*x[12];
        double t95 = x[1]*x[15];
        double t96 = -t94-t95;
        double t109 = x[2]*x[12];
        double t110 = x[3]*x[15];
        double t111 = -t109-t110;
        double t124 = x[4]*x[12];
        double t125 = x[5]*x[15];
        double t126 = -t124-t125;

        M << t13*t5*x[19]+x[21]*(r*(t1*t17+t17*t2-t19*x[7])-x[20]*x[7]+t15+t16)-t29*x[7],
             t43*t36*x[19]+x[21]*(r*(t32*t47+t33*t47-t49*x[9])-x[20]*x[9]+t45+t46)-t59*x[9],
             t73*t66*x[19]+x[21]*(r*(t62*t77+t63*t77-t79*x[11])-x[20]*x[11]+t75+t76)-t89*x[11],
             -t13*t5*x[18]+x[21]*(r*(t1*t96+t19*x[6]+t2*t96)+x[20]*x[6]-t94-t95)+t29*x[6],
             -t43*t36*x[18]+x[21]*(r*(t111*t32+t111*t33+t49*x[8])+x[20]*x[8]-t109-t110)+t59*x[8],
             -t73*t66*x[18]+x[21]*(r*(t126*t62+t126*t63+t79*x[10])+x[20]*x[10]-t124-t125)+t89*x[10],
             x[21]*t11*(-x[6]*x[19]+x[7]*x[18])-x[6]*t17-t96*x[7],
             x[21]*t42*(-x[8]*x[19]+x[9]*x[18])-x[8]*t47-t111*x[9],
             x[21]*(t68*r+t69*r+1)*(-x[10]*x[19]+x[11]*x[18])-x[10]*t77-t126*x[11];

        M.transposeInPlace();
        return M;
    }
}  // namespace ValtonenOrnhagArxiv2021Extra
}  // namespace DronePoseLib
