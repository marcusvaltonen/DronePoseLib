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

#ifndef SRC_RANSAC_VALTONENOERNHAGARXIV2021_HPP_
#define SRC_RANSAC_VALTONENOERNHAGARXIV2021_HPP_

#include <Eigen/Dense>
#include "pose_estimator.hpp"
#include "distortion.hpp"
#include "refinement.hpp"


// TODO: Change name to reflect frEfr
namespace DronePoseLib {
namespace ValtonenOrnhagArxiv2021 {
class Solver : public PoseEstimator<Solver> {
    public:
        Solver() = default;
        int solve(
            DronePoseLib::Points2D &x1,
            DronePoseLib::Points2D &x2,
            const Eigen::Matrix3d &R1,
            const Eigen::Matrix3d &R2,
            std::vector<DronePoseLib::Camera> *poses
        ) const;
        int minimal_sample_size() const {
            return 4;
        }
        inline void distort(
            const std::vector<double>& dist_params, const Eigen::Matrix<double, 2, Eigen::Dynamic> xu, Eigen::Matrix<double, 2, Eigen::Dynamic>* xd) const {
            inverse_1param_division_model(dist_params[0], xu, xd);
        }

        inline void refine(Camera &pose, const Eigen::Matrix<double, 2, Eigen::Dynamic> &x, const Eigen::Matrix<double, 3, Eigen::Dynamic> &X) const {
            pose.dist_params[0] = pose.dist_params[0] * pose.focal * pose.focal;
            DronePoseLib::refinement_undist_with_structure(x, X, pose, 0, 1);
            pose.dist_params[0] = pose.dist_params[0] / pose.focal / pose.focal;
        }
    };
};  // ValtonenOrnhagArxiv2021
};  // DronePoseLib

#endif  // SRC_RANSAC_VALTONENOERNHAGARXIV2021_HPP_
