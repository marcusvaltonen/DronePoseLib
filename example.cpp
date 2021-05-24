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


#define _USE_MATH_DEFINES
#include <Eigen/Geometry>
#include <limits>
#include <RansacLib/ransac.h>

#include <cmath>
#include <numeric>
#include "ransac_estimator.hpp"
#include "ransac_valtonenornhag_arxiv_2021.hpp"
#include "scene_and_pose_generation.hpp"




int main() {
    /* Timing experiments */
    int nbr_iter = 1e2;

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


    // Test RANSAC
	std::cout << "Running RANSAC tests...\n\n";

	Eigen::Matrix<double, 2, Eigen::Dynamic> xx1;
	Eigen::Matrix<double, 2, Eigen::Dynamic> xx2;
	DronePoseLib::Camera pose_gt;

	double dist_param = -0.0000001;

    DronePoseLib::ValtonenOrnhagArxiv2021::Solver estimator;

	generate_scene_and_image(100, 2, 10, 70, false, &pose_gt, &xx1, &xx2, 1.0);
	add_focal(2000.0, &pose_gt, &xx1);
	add_focal(2000.0, &pose_gt, &xx2);
	add_distortion_1pdiv(dist_param, &pose_gt, &xx1);
	add_distortion_1pdiv(dist_param, &pose_gt, &xx2);
	add_noise(0.5, &xx1);
	add_noise(0.5, &xx2);

    std::cout << "===  GT pose ===" << std::endl;
    debug_print_pose(pose_gt);
    std::cout << "xx1 = \n" << xx1 << std::endl;
    std::cout << "xx2 = \n" << xx2 << std::endl;

    // Outliers
    for (int i = 0; i < 20; ++i) {
		Vector2d n;
        n.setRandom();
        n *= 0.2 * pose_gt.focal;
		xx1.col(i) += n;
        n.setRandom();
        n *= 0.2 * pose_gt.focal;
		xx2.col(i) += n;
	}

    // We consider the relative pose problem
    R1.setIdentity();
    R2 = pose_gt.R;

	DronePoseLib::RansacEstimator<DronePoseLib::ValtonenOrnhagArxiv2021::Solver> solver(xx1, xx2, R1, R2, estimator);

	ransac_lib::LORansacOptions options;
	options.squared_inlier_threshold_ = 1;
    std::random_device rand_dev;
    options.random_seed_ = rand_dev();

	ransac_lib::LocallyOptimizedMSAC<
        DronePoseLib::Camera,
		std::vector<DronePoseLib::Camera>,
		DronePoseLib::RansacEstimator<DronePoseLib::ValtonenOrnhagArxiv2021::Solver>
    > lomsac;
	ransac_lib::RansacStatistics ransac_stats;

	DronePoseLib::Camera best_model;
	int num_ransac_inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);

	std::cout << "   ... LOMSAC found " << num_ransac_inliers
		<< " inliers in " << ransac_stats.num_iterations
		<< " iterations with an inlier ratio of "
		<< ransac_stats.inlier_ratio << std::endl;

    std::cout << "GT\n";
    debug_print_pose(pose_gt, true);
    std::cout << "Est\n";
    debug_print_pose(best_model, true);

    return 0;
}
