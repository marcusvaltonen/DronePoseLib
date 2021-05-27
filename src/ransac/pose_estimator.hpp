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

#ifndef SRC_RANSAC_POSE_ESTIMATOR_HPP_
#define SRC_RANSAC_POSE_ESTIMATOR_HPP_

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "relpose.hpp"

namespace DronePoseLib {
// TODO(marcusvaltonen): Consider moving Points2D
static const int MAX_SAMPLE_SIZE = 8;
typedef Eigen::Matrix<double, 2, Eigen::Dynamic, 0, 2, MAX_SAMPLE_SIZE> Points2D;

// We use CRTP here for the solvers.
template<class Solver>
class PoseEstimator {
 public:
  int estimate(
    const DronePoseLib::Points2D &image_points1,
    const DronePoseLib::Points2D &image_points2,
    const Eigen::Matrix3d &rotation_matrix1,
    const Eigen::Matrix3d &rotation_matrix2,
    std::vector<DronePoseLib::Camera> *poses) const;

  inline int minimal_sample_size() const {
      return static_cast<const Solver*>(this)->minimal_sample_size();
  }

 protected:
  PoseEstimator() = default;
};
};  // namespace DronePoseLib

template<class Solver>
int DronePoseLib::PoseEstimator<Solver>::estimate(
    const DronePoseLib::Points2D &image_points1,
    const DronePoseLib::Points2D &image_points2,
    const Eigen::Matrix3d &rotation_matrix1,
    const Eigen::Matrix3d &rotation_matrix2,
    std::vector<DronePoseLib::Camera> *poses) const {
    DronePoseLib::Points2D x1 = image_points1;
    DronePoseLib::Points2D x2 = image_points2;
    Eigen::Matrix3d R1 = rotation_matrix1;
    Eigen::Matrix3d R2 = rotation_matrix2;

    // Call solver implementation
    poses->clear();
    int n_sols = static_cast<const Solver*>(this)->solve(x1, x2, R1, R2, poses);

    return n_sols;
}

#endif  // SRC_RANSAC_POSE_ESTIMATOR_HPP_
