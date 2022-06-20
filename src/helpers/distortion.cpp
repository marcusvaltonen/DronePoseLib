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

#include "distortion.hpp"
#include <Eigen/Dense>

void DronePoseLib::inverse_1param_division_model(
    double lambda,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x0,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *x1) {

    if (lambda == 0.0) {
        *x1 = x0;
        return;
    }

    Eigen::Array<double, 1, Eigen::Dynamic> ru2 = x0.colwise().squaredNorm();
    Eigen::Array<double, 1, Eigen::Dynamic> ru = ru2.sqrt();

    double dist_sign = 1.0;
    if (lambda < 0) {
        dist_sign = -1.0;
    }

    Eigen::Array<double, 1, Eigen::Dynamic> rd;
    rd.resizeLike(ru);
    rd = (0.5 / lambda) / ru - dist_sign * ((0.25 / (lambda * lambda)) / ru2 - 1 / lambda).sqrt();

    rd /= ru;
    x1->resizeLike(x0);
    x1->row(0) = x0.row(0).array() * rd;
    x1->row(1) = x0.row(1).array() * rd;
}

void DronePoseLib::forward_1param_division_model(
    double lambda,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x0,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *x1) {

    Eigen::Array<double, 1, Eigen::Dynamic> r2 = x0.colwise().squaredNorm();
    r2 *= lambda;
    r2 += 1.0;

    x1->resizeLike(x0);
    x1->row(0) = x0.row(0).array() / r2;
    x1->row(1) = x0.row(1).array() / r2;
}
