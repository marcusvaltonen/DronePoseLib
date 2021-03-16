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

TEST_CASE("Valtonen Ornhag Arxiv 2021") {
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

    expected << -0.009941290102839,  0.014235465496328,  0.113085620101434,
                -0.009043328530750,  0.002224958928576,  0.046759882922909,
                -0.126284066846820, -0.103060279015522, -0.048795003486796;
    REQUIRE(poses[0].F.isApprox(expected, tol));

    expected <<  0.001822823155453, -0.001569973701701,  0.049768207832915,
                 0.002239402785454, -0.002902488849921,  0.067758142314926,
                -0.037400514674608, -0.085266803309565, -0.222902577521412;
    REQUIRE(poses[0].F.isApprox(expected, tol));

    expected <<  0.000133718765037, -0.000263493940881, -0.008430680430705,
                 0.000223873616219, -0.000004914443643, -0.004949558808401,
                 0.012044452283588,  0.012002302320058, -0.008564291893685;
    REQUIRE(poses[0].F.isApprox(expected, tol));

    expected <<  0.000046254516421, -0.000031361781069,  0.004923310957907,
                 0.000228045091248, -0.000266551777142,  0.022621928996723,
                -0.000276213777084, -0.024515070759081, -0.393783568131982;
    REQUIRE(poses[0].F.isApprox(expected, tol));
}
