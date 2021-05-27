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

#ifndef SRC_HELPERS_SCENE_AND_POSE_GENERATION_HPP_
#define SRC_HELPERS_SCENE_AND_POSE_GENERATION_HPP_

#include <Eigen/Dense>
#include <vector>
#include "relpose.hpp"

void set_random_pose(DronePoseLib::Camera *pose, double translation_scaling = 1.0);

void generate_scene_and_image(
    int N,
    double min_depth,
    double max_depth,
    double h_fov,
    DronePoseLib::Camera *pose,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points1,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points2,
    double translation_scaling = 1.0
);

/* Adds focal length to the camera and image points.
  Note that the order of add_focal and add_distortion* matters! */
void add_focal(
    double focal,
    DronePoseLib::Camera *pose,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points
);

/* Adds 1 param. division model to the camera and image points.
  Note that the order of add_focal and add_distortion* matters! */
void add_distortion(
    double lambda,
    DronePoseLib::Camera *pose,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points
);

/* Adds random noise to the image points. */
void add_noise(
    double sigma,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points
);

void debug_print_pose(DronePoseLib::Camera pose, bool relative = false);

#endif  // SRC_HELPERS_SCENE_AND_POSE_GENERATION_HPP_
