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
#include "distortion.hpp"
#include "relpose.hpp"
#include "pose_estimator.hpp"
#include "triangulate.hpp"

// TODO: No const& declaration as input to class? Check guidelines
// NOTE: Points2D is optimal for up to 8 (non-minimal).


//DEBUG
#include <iostream>
#include "scene_and_pose_generation.hpp"

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
        image_points1 = p1;
        image_points2 = p2;
        rotation_matrix1 = R1;
        rotation_matrix2 = R2;
        solver = est;

        // Use default settings
        DronePoseLib::RefinementSettings s;
        settings = s;
    }
    RansacEstimator(
        Eigen::Matrix<double, 2, Eigen::Dynamic> p1,
        Eigen::Matrix<double, 2, Eigen::Dynamic> p2,
        Eigen::Matrix3d R1,
        Eigen::Matrix3d R2,
        Solver est,
        DronePoseLib::RefinementSettings s) : RansacEstimator(p1, p2, R1, R2) {
        this->settings = s;
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
        Eigen::Vector3d X;
        bool succ = DronePoseLib::triangulate(pose, image_points1.col(i), image_points2.col(i), &X);

        if (!succ) {
            return std::numeric_limits<double>::max();
        }

        // Measure distance in distorted space
        // TODO: Make distortion accept Vector2d as well
        Eigen::Matrix<double, 2, Eigen::Dynamic> x1_p = pose.focal * X.hnormalized();
        DronePoseLib::inverse_1param_division_model(pose.dist_params[0], x1_p, &x1_p);

        X = pose.R * X + pose.t;
        Eigen::Matrix<double, 2, Eigen::Dynamic> x2_p = pose.focal * X.hnormalized();
        DronePoseLib::inverse_1param_division_model(pose.dist_params[0], x2_p, &x2_p);

        double val = (x1_p - image_points1.col(i)).squaredNorm() + (x2_p - image_points2.col(i)).squaredNorm();
        return val;
    }

    // Linear least squares solver. Calls NonMinimalSolver.
    inline void LeastSquares(const std::vector<int>& sample,
        DronePoseLib::Camera* p) const {
        if (!use_local_opt)
            return;
        Eigen::Matrix<double, 2, Eigen::Dynamic> p1(2, sample.size());
        Eigen::Matrix<double, 2, Eigen::Dynamic> p2(2, sample.size());
        Eigen::Matrix<double, 3, Eigen::Dynamic> X(3, sample.size());
        Eigen::Vector3d Xi;

        for (int i = 0; i < sample.size(); i++) {
            p1.col(i) = image_points1.col(sample[i]);
            p2.col(i) = image_points2.col(sample[i]);

            bool succ = DronePoseLib::triangulate(*p, p1.col(i), p2.col(i), &Xi);
            X.col(i) = Xi;

            if (!succ) {
                // TODO: Do nothing?
            }
        }
        solver.refine(*p, p1, p2, X, settings);
    }

    bool use_non_minimal = true;
    bool use_local_opt = true;
private:
    Solver solver;
    DronePoseLib::RefinementSettings settings;
    Eigen::Matrix<double, 2, Eigen::Dynamic> image_points1;
    Eigen::Matrix<double, 2, Eigen::Dynamic> image_points2;
    Eigen::Matrix3d rotation_matrix1;
    Eigen::Matrix3d rotation_matrix2;
};
}
#endif  // SRC_RANSAC_RANSAC_ESTIMATOR_HPP_
