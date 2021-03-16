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

    p1 << -0.170344517767220,  0.108864416575805,  0.661903907097742,
           1.542599674500269,  0.956408496791784,  2.441464813306650;
    p2 <<  0.729243664675995,  0.914374886276414,  3.761356659359217,
          -0.571732334703832, -1.329303225626401, -1.507683466738961;
    R1 <<  0.775811502326779, -0.525223981869245, -0.349651657692168,
           0.591756227378766,  0.413366286756765,  0.692064216912979,
          -0.218954516317697, -0.743819925682513,  0.631498881979806;
    R2 <<  0.904829024395140,  0.344671896790600,  0.249971438718325,
          -0.216541885342809,  0.878022117084280, -0.426833426295341,
          -0.366597938488913,  0.332081986072128,  0.869095797954443;

    std::vector<DronePoseLib::RelPose> poses;
    poses = DronePoseLib::ValtonenOrnhagArxiv2021::get_fEf(p1, p2, R1, R2);
    double tol = 1e-12;

    // Test distortion parameter and focal length
    REQUIRE(poses[0].f == Approx(1.194848331672758).margin(tol));

    // Test homography
    Eigen::Matrix3d expected;
    expected << 0.73948604786043, -1.78006293093289,  2.34478746256462,
               -1.62833579813011,  1.72202744447341, -2.52052334684554,
                2.09573004733515, -1.58339929944571,  2.07673762020455;
    REQUIRE(poses[0].F.isApprox(expected, tol));
}
