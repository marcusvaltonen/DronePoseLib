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
#include "distortion.hpp"
#include "scene_and_pose_generation.hpp"
#include "relpose.hpp"
#include "triangulate.hpp"

TEST_CASE("triangulate") {
    DronePoseLib::Camera pose;
    pose.R <<
       0.114117755567683,   0.720696770083189,  -0.683793319253599,
      -0.416731265253450,   0.659533743021984,   0.625579966411706,
       0.901838228522418,   0.213568273409398,   0.375601387335664;
    pose.t <<
          394.963979386373,
          310.693571838401,
         -1461.80694350949;
    pose.focal = 0.673139392770936;
    pose.dist_params.push_back(-0.00425062473859951);

    Eigen::Vector2d p1;
    p1 << 2.76452368277409,
          13.8755682729254;
    Eigen::Vector2d p2;
    p2 << 0.984409768200559,
          7.54701555819437;

    // Triangulate
    Eigen::Vector3d X;
    bool succ = triangulate(pose, p1, p2, &X);

    // Test output
    double tol = 1e-12;
    Eigen::Vector3d expected;
    expected << 324.525725801189,
                5837.4883872621,
                2298.86986335411;

    REQUIRE(succ);
    REQUIRE(X.isApprox(expected, tol));
}

TEST_CASE("distortion") {
    int N = 5;
    double lambda = -0.01;
    Eigen::Matrix<double, 2, Eigen::Dynamic> x(2, N);
    x.row(0).setLinSpaced(N, -1.0, 1.0);
    x.row(1).setLinSpaced(N, 1.0, 5.0);

    double tol = 1e-12;
    Eigen::Matrix<double, 2, Eigen::Dynamic> y(2, N);
    Eigen::Matrix<double, 2, Eigen::Dynamic> expected(2, N);

    // Distort
    expected << -0.980762113533166, -0.480384603759981,                  0,  0.437728089025404,  0.823626318670327,
                 0.980762113533166,   1.92153841503992,   2.76983964948433,   3.50182471220324,   4.11813159335163;
    DronePoseLib::inverse_1param_division_model(lambda, x, &y);
    REQUIRE(y.isApprox(expected, tol));

    // Undistort
    expected << -1.02040816326531, -0.522193211488251,                  0,  0.597014925373134,   1.35135135135135,
                 1.02040816326531,     2.088772845953,    3.2967032967033,   4.77611940298507,   6.75675675675676;
    DronePoseLib::forward_1param_division_model(lambda, x, &y);
    REQUIRE(y.isApprox(expected, tol));
}

TEST_CASE("add_parameters") {
    DronePoseLib::Camera pose;
    double focal = 2.0;
    double dist_param = -0.01;

    Eigen::Matrix<double, 2, Eigen::Dynamic> x(2, 3);
    x << 1, 2, 3,
         4, 5, 6;

    add_focal(focal, &pose, &x);

    REQUIRE(focal == pose.focal);

    Eigen::Matrix<double, 2, Eigen::Dynamic> expected(2, 3);
    expected << 2, 4, 6,
                8, 10, 12;
    double tol = 1e-12;
    REQUIRE(x.isApprox(expected, tol));

    add_distortion(dist_param, &pose, &x);
    REQUIRE(dist_param == pose.dist_params[0]);

    expected << 1.36577963558616, 2.37046278863376, 3.10594035442545,
                5.46311854234465, 5.92615697158441,  6.2118807088509;
    REQUIRE(x.isApprox(expected, tol));
}
