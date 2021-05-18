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

#ifndef INCLUDES_DRONEPOSELIB_RELPOSE_HPP_
#define INCLUDES_DRONEPOSELIB_RELPOSE_HPP_

#include <vector>
#include <Eigen/Dense>

namespace DronePoseLib {
struct RelPose {
    Eigen::Matrix3d F;
    Eigen::Vector3d t;
    double f;
    double r;
};
struct Camera {
    Camera() : focal(1.0) {}
    Camera(Eigen::Matrix3d rot, Eigen::Vector3d trans, double f) : R(rot), t(trans), focal(f) {};
    Camera(Eigen::Matrix3d rot, Eigen::Vector3d trans) : R(rot), t(trans), focal(1.0) {};
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    double focal;
    std::vector<double> dist_params;
};
}  // namespace DronePoseLib

#endif  // INCLUDES_DRONEPOSELIB_RELPOSE_HPP_
