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
#include <catch2/catch.hpp>
#include "get_valtonenornhag_arxiv_2021.hpp"
#include "relpose.hpp"

TEST_CASE("Valtonen Ornhag Arxiv 2021 - fEf") {
    Eigen::MatrixXd p1(2, 3);
    Eigen::MatrixXd p2(2, 3);
    Eigen::Matrix3d R1, R2;

    p1 << -34.998621048798569, -0.064090328787399, -8.738571263872560,
          -26.916377657772212, 10.374453893285320, -8.028592688400703;

    p2 << -164.2306957039563, -4.85637790217750, -44.47452317786740,
           18.59197557487190,  33.7875137840784,  17.86327225402930;

    R1 << -0.180619577025105,  0.225943086758464, -0.957249335304720,
          -0.982474033611525,  0.004129208291573,  0.186353757456038,
           0.046058025081098,  0.974131752477529,  0.221237399959154;

    R2 << -0.460139184327122,  0.450642405815460, -0.764979315490049,
          -0.848730576230669, -0.476191421985198,  0.229995953440211,
          -0.260630658246353,  0.755091485654930,  0.601588321257572;

    std::vector<DronePoseLib::RelPose> poses;
    poses = DronePoseLib::ValtonenOrnhagArxiv2021::get_fEf(p1, p2, R1, R2);
    double tol = 1e-12;

    // Test size
    REQUIRE(poses.size() == 4);

    // Test focal length
    REQUIRE(poses[0].f == Approx(5.998701610439166).margin(tol));
    REQUIRE(poses[1].f == Approx(-10.455953354056966).margin(tol));
    REQUIRE(poses[2].f == Approx(51.610462716016606).margin(tol));
    REQUIRE(poses[3].f == Approx(-37.468067324800963).margin(tol));

    // Test fundamental matrix
    tol = 1e-7;
    Eigen::Matrix3d expected;

    // Normalize putative fundamental matrices
    for (int i=0; i < poses.size(); i++) {
        poses[i].F = poses[i].F / poses[i].F(2, 2);
    }

    expected <<  0.203735821138506,   -0.29174022910316,   -2.31756557066407,
                 0.185333085040082, -0.0455980893449054,   -0.95829243942079,
                 2.588053239528700,    2.11210721694911,                   1;
    REQUIRE(poses[0].F.isApprox(expected, tol));

    expected << -0.00817766746226784,  0.00704331784387064,   -0.223273361781265,
                -0.01004655401635660,    0.013021333724335,   -0.303980972622073,
                 0.16778861460682200,    0.382529463130022,                    1;
    REQUIRE(poses[1].F.isApprox(expected, tol));

    expected << -0.0156135226001851,   0.0307665764025506,    0.984399006404908,
                -0.0261403533412491, 0.000573829535995767,    0.577929719099143,
                -1.4063570500754600,    -1.40143545655041,                    1;
    REQUIRE(poses[2].F.isApprox(expected, tol));

    expected << -0.000117461773838000,  7.96421781084678e-05,   -0.0125025810022051,
                -0.000579112765750187,  0.000676899187049556,   -0.0574476205394655,
                 0.000701435507820614,    0.0622551897616618,                     1;
    REQUIRE(poses[3].F.isApprox(expected, tol));
}

TEST_CASE("Valtonen Ornhag Arxiv 2021 - frEfr") {
    Eigen::MatrixXd p1(2, 4);
    Eigen::MatrixXd p2(2, 4);
    Eigen::Matrix3d R1, R2;

    p1 <<  99.0825859985542, 1136.84241396289, -1031.66650596755, -117.325418998838,
          -301.923289351466, 1760.62028612233, -533.989983528509,  566.900954605729;

    p2 <<  1829.78818884974, 15378.7866880232,  612.159309750213,  2756.44403323769,
           474.180433404958, 4677.92468041337,  1092.76420021176,  1874.89973953780;

    R1 <<  0.983417707845482,  0.013453875580959,  0.180855204867843,
          -0.060089831339915,  0.965084992860598,  0.254951306575981,
          -0.171110560940808, -0.261591188282616,  0.949890112669571,

    R2 <<  0.556837962037774,  0.329755275145440,  0.762360113428931,
          -0.787076103923440,  0.502747650438714,  0.357429722618378,
          -0.265410419287405, -0.799065866178820,  0.539491474299248;

    std::vector<DronePoseLib::RelPose> poses;
    bool use_fast_solver = false;
    poses = DronePoseLib::ValtonenOrnhagArxiv2021::get_frEfr(p1, p2, R1, R2, use_fast_solver);
    double tol = 1e-12;

    // Test size
    REQUIRE(poses.size() == 3);

    // Test focal length
    REQUIRE(poses[0].f == Approx(1081.280695559530).margin(tol));
    REQUIRE(poses[1].f == Approx(7.142069186813480).margin(tol));
    REQUIRE(poses[2].f == Approx(96.42064345373980).margin(tol));

    // Test radial distortion coefficient
    REQUIRE(poses[0].r == Approx(-0.00000005366375).margin(tol));
    REQUIRE(poses[1].r == Approx(-0.00003572050787).margin(tol));
    REQUIRE(poses[2].r == Approx(-0.00000634720100).margin(tol));

    // Test fundamental matrix
    tol = 1e-7;
    Eigen::Matrix3d expected;

    // Normalize putative fundamental matrices
    for (int i=0; i < poses.size(); i++) {
        poses[i].F = poses[i].F / poses[i].F(2, 2);
    }
    expected << -1.36694183091702e-06,  2.28782882252314e-07,   0.00166084307538912,
                -3.37831263794585e-07,  2.00768808300091e-06,  -0.00469955999290139,
                 -0.00438068604065931,   0.00238365873365962,                     1;
    REQUIRE(poses[0].F.isApprox(expected, tol));

    expected << -0.00191334166323602,   0.0127558338972733,  -0.0630718336078865,
                 0.00168816158196316,   0.0183267091593544,    -0.21347550674579,
                  -0.102757406138353,    0.202821018650091,                    1;
    REQUIRE(poses[1].F.isApprox(expected, tol));

    expected << 0.000164734835612738, 3.07215616172707e-05,  -0.0228114435936927,
                0.000114053882939565, 5.49251999934715e-05,   0.0129104230967059,
                  0.0185530116486357,  0.00764472822544315,                    1;
    REQUIRE(poses[2].F.isApprox(expected, tol));
}

TEST_CASE("Valtonen Ornhag Arxiv 2021 - frEfr - fast option") {
    Eigen::MatrixXd p1(2, 4);
    Eigen::MatrixXd p2(2, 4);
    Eigen::Matrix3d R1, R2;

    p1 <<  99.0825859985542, 1136.84241396289, -1031.66650596755, -117.325418998838,
          -301.923289351466, 1760.62028612233, -533.989983528509,  566.900954605729;

    p2 <<  1829.78818884974, 15378.7866880232,  612.159309750213,  2756.44403323769,
           474.180433404958, 4677.92468041337,  1092.76420021176,  1874.89973953780;

    R1 <<  0.983417707845482,  0.013453875580959,  0.180855204867843,
          -0.060089831339915,  0.965084992860598,  0.254951306575981,
          -0.171110560940808, -0.261591188282616,  0.949890112669571,

    R2 <<  0.556837962037774,  0.329755275145440,  0.762360113428931,
          -0.787076103923440,  0.502747650438714,  0.357429722618378,
          -0.265410419287405, -0.799065866178820,  0.539491474299248;

    std::vector<DronePoseLib::RelPose> poses;
    bool use_fast_solver = true;
    poses = DronePoseLib::ValtonenOrnhagArxiv2021::get_frEfr(p1, p2, R1, R2, use_fast_solver);
    double tol = 1e-12;

    // Test size
    REQUIRE(poses.size() == 3);

    // Test focal length
    REQUIRE(poses[0].f == Approx(7.1420691867328010).margin(tol));
    REQUIRE(poses[1].f == Approx(96.420643453768420).margin(tol));
    REQUIRE(poses[2].f == Approx(1081.2806955548920).margin(tol));

    // Test radial distortion coefficient
    REQUIRE(poses[0].r == Approx(-0.00003572050787).margin(tol));
    REQUIRE(poses[1].r == Approx(-0.00000634720100).margin(tol));
    REQUIRE(poses[2].r == Approx(-0.00000005366375).margin(tol));

    // Test fundamental matrix
    tol = 1e-7;
    Eigen::Matrix3d expected;

    // Normalize putative fundamental matrices
    for (int i=0; i < poses.size(); i++) {
        poses[i].F = poses[i].F / poses[i].F(2, 2);
    }

    expected << -0.00191334166327502,    0.012755833897505,   -0.063071833608412,
                 0.00168816158199461,   0.0183267091596996,   -0.213475506747847,
                  -0.102757406139415,    0.202821018651993,                    1;
    REQUIRE(poses[0].F.isApprox(expected, tol));
    expected <<  0.000164734835612805, 3.07215616172335e-05,  -0.0228114435936991,
                 0.000114053882939612, 5.49251999934633e-05,   0.0129104230967174,
                   0.0185530116486443,  0.00764472822544001,                    1;
    REQUIRE(poses[1].F.isApprox(expected, tol));
    expected << -1.36694183127657e-06,  2.28782882538085e-07,   0.00166084307560333,
                -3.37831264124548e-07,  2.00768808281651e-06,  -0.00469955999324961,
                 -0.00438068604074091,   0.00238365873368763,                     1;
    REQUIRE(poses[2].F.isApprox(expected, tol));
}

TEST_CASE("Valtonen Ornhag Arxiv 2021 Extra - rEr") {
    Eigen::MatrixXd p1(2, 3);
    Eigen::MatrixXd p2(2, 3);
    Eigen::Matrix3d R1, R2;

    p1 <<  99.0825859985542, 1136.84241396289, -1031.66650596755,
          -301.923289351466, 1760.62028612233, -533.989983528509;

    p2 <<  1829.78818884974, 15378.7866880232,  612.159309750213,
           474.180433404958, 4677.92468041337,  1092.76420021176;

    R1 <<  0.983417707845482,  0.013453875580959,  0.180855204867843,
          -0.060089831339915,  0.965084992860598,  0.254951306575981,
          -0.171110560940808, -0.261591188282616,  0.949890112669571,

    R2 <<  0.556837962037774,  0.329755275145440,  0.762360113428931,
          -0.787076103923440,  0.502747650438714,  0.357429722618378,
          -0.265410419287405, -0.799065866178820,  0.539491474299248;

    std::vector<DronePoseLib::RelPose> poses;
    double focal_length = 1.0;
    poses = DronePoseLib::ValtonenOrnhagArxiv2021Extra::get_rEr(p1, p2, R1, R2, focal_length);
    double tol = 1e-12;

    // Test size
    REQUIRE(poses.size() == 4);

    // Test focal length
    REQUIRE(poses[0].f == Approx(1.0).margin(tol));
    REQUIRE(poses[1].f == Approx(1.0).margin(tol));
    REQUIRE(poses[2].f == Approx(1.0).margin(tol));
    REQUIRE(poses[3].f == Approx(1.0).margin(tol));

    // Test radial distortion coefficient
    REQUIRE(poses[0].r == Approx(-4.2874876051071e-05).margin(tol));
    REQUIRE(poses[1].r == Approx(-0.000246913930531072).margin(tol));
    REQUIRE(poses[2].r == Approx(0.00184604041266755).margin(tol));
    REQUIRE(poses[3].r == Approx(-0.00109924194761612).margin(tol));

    // Test fundamental matrix
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected <<  0.264934547706507,   0.145240095976332,  -0.465106205345824,
                0.0309556703247373,   -0.39180514557381,   0.776482268850591,
                 0.924467577138297,  -0.282908100262259, -0.0178797053794708;
    REQUIRE(poses[0].F.isApprox(expected, tol));

    expected <<  0.0598818526297369,  -0.332627369542942,   0.217582108961897,
                -0.0388012497427252,  -0.481754588366569,   0.797147303075138,
                  0.393684983857834,  -0.744417584343487,  -0.510146031080832;
    REQUIRE(poses[1].F.isApprox(expected, tol));

    expected <<  0.642219648136942, -0.678026024845748, -0.217842639644248,
                 0.703919384125803,  0.629690865155572,  0.297361918557707,
                -0.295697506917121,  0.109786187781648, 0.0212602749320976;
    REQUIRE(poses[2].F.isApprox(expected, tol));

    expected << 0.591572411289522,  -0.130509314604598,  -0.636919856945032,
                0.410441925689088,  0.0482225493952967,   0.753793979143006,
                0.687481917761856, -0.0525530771809502,   0.113824127649439;
    REQUIRE(poses[3].F.isApprox(expected, tol));
}
