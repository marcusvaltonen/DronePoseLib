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
#include <chrono>  // NOLINT [build/c++11]
#include <iostream>
#include "get_valtonenornhag_arxiv_2021.hpp"
#include "relpose.hpp"

int main() {
    /* Timing experiments */
    int nbr_iter = 1e4;

    // Test fEf
    int N = 3;
    Eigen::MatrixXd x1(2, N);
    x1 = Eigen::MatrixXd::Random(2, N);
    Eigen::MatrixXd x2(2, N);
    x2 = Eigen::MatrixXd::Random(2, N);
    Eigen::Matrix3d R1;
    R1 = Eigen::Matrix3d::Random(3, 3);
    Eigen::Matrix3d R2;
    R2 = Eigen::Matrix3d::Random(3, 3);

    auto start = std::chrono::steady_clock::now();
    std::cout << "=== DronePoseLib timings ===" << std::endl;

    std::vector<DronePoseLib::RelPose> poses;

    for (int i = 0; i < nbr_iter; i++) {
        poses = DronePoseLib::ValtonenOrnhagArxiv2021::get_fEf(x1, x2, R1, R2);
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time (fEf): "
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / nbr_iter
        << " ns" << std::endl;

    // Test rEr
    double focal_length = 1.0;
    start = std::chrono::steady_clock::now();
    for (int i = 0; i < nbr_iter; i++) {
        poses = DronePoseLib::ValtonenOrnhagArxiv2021Extra::get_rEr(x1, x2, R1, R2, focal_length);
    }
    end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time (rEr): "
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / nbr_iter
        << " ns" << std::endl;

    // Test frEfr
    N = 4;
    x1 = Eigen::MatrixXd::Random(2, N);
    x2 = Eigen::MatrixXd::Random(2, N);
    R1 = Eigen::Matrix3d::Random(3, 3);
    R2 = Eigen::Matrix3d::Random(3, 3);

    bool use_fast_solver = false;

    start = std::chrono::steady_clock::now();
    for (int i = 0; i < nbr_iter; i++) {
        poses = DronePoseLib::ValtonenOrnhagArxiv2021::get_frEfr(x1, x2, R1, R2, use_fast_solver);
    }
    end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time (frEfr): "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / nbr_iter
        << " μs" << std::endl;

    use_fast_solver = true;

    start = std::chrono::steady_clock::now();
    for (int i = 0; i < nbr_iter; i++) {
        poses = DronePoseLib::ValtonenOrnhagArxiv2021::get_frEfr(x1, x2, R1, R2, use_fast_solver);
    }
    end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time (frEfr) fast: "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / nbr_iter
        << " μs" << std::endl;

    return 0;
}
