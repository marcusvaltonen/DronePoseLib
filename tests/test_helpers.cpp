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
