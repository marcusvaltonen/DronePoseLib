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

#ifndef SRC_RANSAC_REFINEMENT_HPP_
#define SRC_RANSAC_REFINEMENT_HPP_

#include <Eigen/Dense>
#include "relpose.hpp"

namespace DronePoseLib {
void refinement_undist_with_structure(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x1,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x2,
    DronePoseLib::Camera &p,
    Eigen::Matrix<double, 3, Eigen::Dynamic> &X,
    const DronePoseLib::RefinementSettings &settings
);
};  // DronePoseLib
#endif  // SRC_RANSAC_REFINEMENT_HPP_
