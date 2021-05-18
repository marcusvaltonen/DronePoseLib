#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <limits>
#include <RansacLib/ransac.h>

#include <time.h>
#include <cmath>
#include <iostream>
#include <numeric>
#include "radialpose.h"
#include "misc/ransac_estimator.h"
#include "misc/unit_test_misc.h"



int main() {
	std::cout << "Running tests...\n\n";

	Matrix<double, 2, Dynamic> x1;
	Matrix<double, 2, Dynamic> x2;
	Camera pose_gt;

	std::vector<double> params2 = { -0.12, 0.034 };

    DronePoseLib::ValtonenOrnhagArxiv2021::Solver<true, true> estimator;

	generate_scene_and_image(100, 2, 20, 70, false, &pose_gt, &x1, &x1, 1.0);
	add_rational_distortion(params2, 2, 0, &pose_gt, &x1);
	add_rational_distortion(params2, 2, 0, &pose_gt, &x2);
	add_focal(2000.0, &pose_gt, &x1);
	add_focal(2000.0, &pose_gt, &x2);
	add_noise(0.5, &x1);
	add_noise(0.5, &x2);

	DronePoseLib::RansacEstimator<
        DronePoseLib::ValtonenOrnhagArxiv2021::Solver<true, true>
    > solver(x1, x2, estimator);

	ransac_lib::LORansacOptions options;
	options.squared_inlier_threshold_ = 4;

	ransac_lib::LocallyOptimizedMSAC<Camera,
		std::vector<Camera>,
		DronePoselib::RansacEstimator<DronePoseLib::ValtonenOrnhagArxiv2021::Solver<true, true>>> lomsac;
	ransac_lib::RansacStatistics ransac_stats;

	Camera best_model;
	int num_ransac_inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);

	std::cout << "   ... LOMSAC found " << num_ransac_inliers
		<< " inliers in " << ransac_stats.num_iterations
		<< " iterations with an inlier ratio of "
		<< ransac_stats.inlier_ratio << std::endl;

    return 0;
}
