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

#ifndef SRC_RANSAC_RANSAC_ESTIMATOR_HPP_
#define SRC_RANSAC_RANSAC_ESTIMATOR_HPP_

#include <Eigen/Dense>
#include "relpose.hpp"
#include "pose_estimator.hpp"
#include "triangulate.hpp"

// TODO: No const& declaration as input to class? Check guidelines
// NOTE: Points2D is optimal for up to 8 (non-minimal).


//DEBUG
#include <iostream>

namespace DronePoseLib {
template<class Solver>
class RansacEstimator {
public:
	RansacEstimator(
        Eigen::Matrix<double, 2, Eigen::Dynamic> p1,
        Eigen::Matrix<double, 2, Eigen::Dynamic> p2,
        Eigen::Matrix3d R1,
        Eigen::Matrix3d R2,
        Solver est) {
        // TODO: Make sure p1 and p2 are of the same length (and larger than minimal sample size)
		image_points1 = p1;
		image_points2 = p2;
        rotation_matrix1 = R1;
        rotation_matrix2 = R2;
		solver = est;
	}

	inline int min_sample_size() const {
		return solver.minimal_sample_size();
	}
	inline int non_minimal_sample_size() const {
		return solver.minimal_sample_size() * 2;
	}
	inline int num_data() const {
		return image_points1.cols();
	}

	int MinimalSolver(const std::vector<int>& sample,
		std::vector<Camera>* poses) const {
        std::cout << "called MinimalSolver()" << std::endl;
		Points2D p1(2, sample.size());
		Points2D p2(2, sample.size());

		for (int i = 0; i < sample.size(); i++) {
			p1.col(i) = image_points1.col(sample[i]);
			p2.col(i) = image_points2.col(sample[i]);
		}
		solver.estimate(p1, p2, rotation_matrix1, rotation_matrix2, poses);

		return poses->size();
	}

	// Returns 0 if no model could be estimated and 1 otherwise.
	int NonMinimalSolver(const std::vector<int>& sample,
		Camera* pose) const {
        std::cout << "called NonMinimalSolver()" << std::endl;
		if (!use_non_minimal)
			return 0;

		Eigen::Matrix<double, 2, Eigen::Dynamic> p1(2, sample.size());
		Eigen::Matrix<double, 2, Eigen::Dynamic> p2(2, sample.size());

		for (int i = 0; i < sample.size(); i++) {
			p1.col(i) = image_points1.col(sample[i]);
			p2.col(i) = image_points2.col(sample[i]);
		}

		// Call minimal solver
		std::vector<Camera> poses;
		Points2D p1s = p1.block(0, 0, 2, min_sample_size());
		Points2D p2s = p2.block(0, 0, 2, min_sample_size());
		solver.estimate(p1s, p2s, rotation_matrix1, rotation_matrix2, &poses);

		// For all pose candidates compute score
		double best_score = std::numeric_limits<double>::max();
		int best_idx = -1;

		for (int i = 0; i < poses.size(); ++i) {
			double score = 0;
			for (int j = 0; j < sample.size(); ++j) {
				score += EvaluateModelOnPoint(poses[i], sample[j]);
            }
			if (score < best_score) {
				best_score = score;
				best_idx = i;
			}
		}

		if (best_idx != -1) {
			*pose = poses[best_idx];
			return 1;
		} else {
			return 0;
		}
	}

	// Evaluates the line on the i:th data point
	double EvaluateModelOnPoint(const Camera& pose, int i) const {
        std::cout << "called EvaluateModelOnPoint(): " << i << std::endl;
        Eigen::Vector3d X;
        bool succ = DronePoseLib::triangulate(pose, image_points1.col(i), image_points2.col(i), &X);

        std::cout << "x1 =\n" << image_points1.col(i) << std::endl;
        std::cout << "x2 =\n" << image_points2.col(i) << std::endl;
        std::cout << "triangulate: succ" << succ << std::endl;
        std::cout << "triangulate: X = \n" << X << std::endl;

		// Compute reprojected point in first camera
		Eigen::Matrix<double, 2, Eigen::Dynamic> z1(2, 1);
		z1 << X(0) / X(2), X(1) / X(2);
		solver.distort(pose.dist_params, z1, &z1);
		z1 = pose.focal * z1;

		// Compute reprojected point in second camera
		Eigen::Vector3d Z2 = pose.R * X + pose.t;
		Eigen::Matrix<double, 2, Eigen::Dynamic> z2(2, 1);
		z1 << Z2(0) / Z2(2), Z2(1) / Z2(2);
		solver.distort(pose.dist_params, z2, &z2);
		z2 = pose.focal * z2;

        // TODO: Should we compute it both ways?
        // Maybe Take mean/squared sum between z1-im1 and z2-im2
		return (z2 - image_points2.col(i)).squaredNorm();
	}

	// Linear least squares solver. Calls NonMinimalSolver.
	inline void LeastSquares(const std::vector<int>& sample,
		DronePoseLib::Camera* p) const {
        std::cout << "called LeastSquares()" << std::endl;
		if (!use_local_opt)
			return;
		Eigen::Matrix<double, 2, Eigen::Dynamic> p1(2, sample.size());
		Eigen::Matrix<double, 2, Eigen::Dynamic> p2(2, sample.size());

		for (int i = 0; i < sample.size(); i++) {
			p1.col(i) = image_points1.col(sample[i]);
			p2.col(i) = image_points2.col(sample[i]);
		}

        // TODO: Refine both structure and motion?
        // Right now this does not work..
        // FIXME: This is just a dummy - refinement disabled
		//solver.refine(*p, p1, p2);
		Eigen::Matrix<double, 3, Eigen::Dynamic> dummy(3, sample.size());
		solver.refine(*p, p1, dummy);
	}
	bool use_non_minimal = true;
	bool use_local_opt = false;  // FIXME: This is temporarily disabled
private:
	Solver solver;
	Eigen::Matrix<double, 2, Eigen::Dynamic> image_points1;
	Eigen::Matrix<double, 2, Eigen::Dynamic> image_points2;
    Eigen::Matrix3d rotation_matrix1;
    Eigen::Matrix3d rotation_matrix2;
};
}
#endif  // SRC_RANSAC_RANSAC_ESTIMATOR_HPP_
