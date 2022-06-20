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

#include <stdlib.h>  // srand
#include <RansacLib/ransac.h>
#include <Eigen/Dense>
#include <catch2/catch.hpp>
#include "get_valtonenornhag_arxiv_2021.hpp"
#include "relpose.hpp"
#include "ransac_estimator.hpp"
#include "ransac_valtonenornhag_arxiv_2021.hpp"
#include "scene_and_pose_generation.hpp"


TEST_CASE("RANSAC frEfr - intergration test") {
    Eigen::Matrix3d R1;
    Eigen::Matrix3d R2;
    // Test RANSAC
    Eigen::Matrix<double, 2, Eigen::Dynamic> xx1;
    Eigen::Matrix<double, 2, Eigen::Dynamic> xx2;
    DronePoseLib::Camera pose_gt;

    double dist_param = -0.0000001;

    DronePoseLib::ValtonenOrnhagArxiv2021::Solver estimator;

    // Eigen uses the std library seed
    srand((unsigned int) 0);

    generate_scene_and_image(100, 2, 10, 70, &pose_gt, &xx1, &xx2, 1.0);
    add_focal(2000.0, &pose_gt, &xx1);
    add_focal(2000.0, &pose_gt, &xx2);
    add_distortion(dist_param, &pose_gt, &xx1);
    add_distortion(dist_param, &pose_gt, &xx2);
    add_noise(0.5, &xx1);
    add_noise(0.5, &xx2);

    debug_print_pose(pose_gt);

    // Outliers
    for (int i = 0; i < 20; ++i) {
        Eigen::Vector2d n;
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

    // Refinement settings
    DronePoseLib::RansacEstimator<DronePoseLib::ValtonenOrnhagArxiv2021::Solver> solver(xx1, xx2, R1, R2, estimator);

    // RANSAC parmeters (seed default 0)
    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 1;
    options.success_probability_ = 0.9999;
    options.min_num_iterations_ = 100u;
    options.max_num_iterations_ = 200u;

    ransac_lib::LocallyOptimizedMSAC<
        DronePoseLib::Camera,
        std::vector<DronePoseLib::Camera>,
        DronePoseLib::RansacEstimator<DronePoseLib::ValtonenOrnhagArxiv2021::Solver>
    > lomsac;
    ransac_lib::RansacStatistics ransac_stats;

    DronePoseLib::Camera best_model;
    int num_ransac_inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);

    double tol = 1e-10;
    REQUIRE(num_ransac_inliers == 80);
    REQUIRE(ransac_stats.num_iterations == 100);
    REQUIRE(ransac_stats.best_num_inliers == 80);
    REQUIRE(ransac_stats.inlier_ratio == 0.8);
    for (int i=0; i <ransac_stats.inlier_indices.size(); i++)
        REQUIRE(ransac_stats.inlier_indices[i] == 20 + i);
    REQUIRE(ransac_stats.number_lo_iterations == 3);

    Eigen::Matrix3d R_expected;
    R_expected <<  0.954859221838922,  -0.208314564297171,   0.211775609480524,
                   0.189925545753350,   0.976276785415069,   0.103980408394627,
                  -0.228412244723131, -0.0590650536383497,   0.971771148933348;

    Eigen::Vector3d t_expected;
    t_expected <<  0.01012099831878770,
                 -4.69487630234731e-05,
                   0.00204335953247533;

    debug_print_pose(best_model);

    REQUIRE(best_model.R.isApprox(R_expected, tol));
    REQUIRE(best_model.t.isApprox(t_expected, tol));
    REQUIRE(best_model.focal == Approx(1993.685303806515).margin(tol));
    REQUIRE(best_model.dist_params[0] == Approx(-1.00907826044204e-07).margin(tol));
    REQUIRE(ransac_stats.best_model_score == Approx(26.6962845521).margin(tol));
}
