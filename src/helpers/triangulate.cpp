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

#include "triangulate.hpp"
#include <Eigen/Dense>
#include "distortion.hpp"
#include "relpose.hpp"

//DEBUG
#include <iostream>

namespace DronePoseLib {
bool triangulate(const Camera& pose, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, Eigen::Vector3d *t) {

    // Normalize points
    Eigen::Matrix<double, 2, Eigen::Dynamic> x1, x2;
    DronePoseLib::forward_1param_division_model(pose.dist_params[0], p1, &x1);
    DronePoseLib::forward_1param_division_model(pose.dist_params[0], p2, &x2);
    x1 /= pose.focal;
    x2 /= pose.focal;

    std::cout << "x1 = \n" << x1 << std::endl;
    std::cout << "x2 = \n" << x2 << std::endl;

    // First pose is assumed to be the identity
    // TODO: Fixed 6x6
    Eigen::MatrixXd M(6, 6);
    M.setZero();
    M.topLeftCorner(3, 3).setIdentity();
    M.bottomLeftCorner(3, 3) = pose.R;
    M.block(3, 3, 3, 1) = pose.t;
    M.block(0, 4, 2, 1) = -x1;
    M(2, 4) = -1;
    M.block(3, 5, 2, 1) = -x2;
    M(5, 5) = -1;

    std::cout << "M = \n" << M << std::endl;

    // Extract solution via SVD
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullV);
    Eigen::Matrix<double, 6, 6> Q = svd.matrixV();
    std::cout << "Q = \n" << Q << std::endl;
	(*t) = Q.block(0, 5, 3, 1) / Q(3, 5);

    // Check if point is in front of the camera
    bool succ = Q(4,5) / Q(3,5) > 0 && Q(5,5) / Q(3,5) > 0;

    return succ;
}
}
