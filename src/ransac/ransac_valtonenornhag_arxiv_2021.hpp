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

#ifndef SRC_RANSAC_RANSAC_VALTONENORNHAG_ARXIV_2021_HPP_
#define SRC_RANSAC_RANSAC_VALTONENORNHAG_ARXIV_2021_HPP_

#include <Eigen/Dense>
#include <vector>
#include "pose_estimator.hpp"
#include "distortion.hpp"
#include "refinement.hpp"

// TODO(marcusvaltonen): Change name to reflect frEfr
namespace DronePoseLib {
namespace ValtonenOrnhagArxiv2021 {
class Solver : public PoseEstimator<Solver> {
 public:
    Solver() = default;
    int solve(
        const DronePoseLib::Points2D &x1,
        const DronePoseLib::Points2D &x2,
        const Eigen::Matrix3d &R1,
        const Eigen::Matrix3d &R2,
        std::vector<DronePoseLib::Camera> *poses) const;

    int minimal_sample_size() const {
        return 4;
    }

    inline void distort(
        const std::vector<double> &dist_params,
        const Eigen::Matrix<double, 2, Eigen::Dynamic> xu,
        Eigen::Matrix<double, 2, Eigen::Dynamic> *xd) const {
        inverse_1param_division_model(dist_params[0], xu, xd);
    }

    inline void refine(
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &x1,
        const Eigen::Matrix<double, 2, Eigen::Dynamic> &x2,
        DronePoseLib::Camera &pose,
        Eigen::Matrix<double, 3, Eigen::Dynamic> &X,
        const DronePoseLib::RefinementSettings &settings) const {
        // TODO(marcusvaltonen): Can (or should we) avoid this situation?
        // pose->dist_params[0] = pose->dist_params[0] * pose->focal * pose->focal;
        pose.dist_params[0] = pose.dist_params[0] * pose.focal * pose.focal;
        DronePoseLib::refinement_undist_with_structure(x1, x2, pose, X, settings);
        // pose->dist_params[0] = pose->dist_params[0] / pose->focal / pose->focal;
        pose.dist_params[0] = pose.dist_params[0] / pose.focal / pose.focal;
    }
};
};  // namespace ValtonenOrnhagArxiv2021
};  // namespace DronePoseLib

#endif  // SRC_RANSAC_RANSAC_VALTONENORNHAG_ARXIV_2021_HPP_
