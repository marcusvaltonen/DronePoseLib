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
//
// Based on an implementation by Viktor Larsson, see radialpose.

#include <stdlib.h>   // rand
#define _USE_MATH_DEFINES
#include <math.h>
#include <Eigen/Dense>
#include <random>
#include <iostream>
#include <iomanip>
#include "scene_and_pose_generation.hpp"
#include "relpose.hpp"
#include "distortion.hpp"

void set_random_pose(DronePoseLib::Camera *pose, double translation_scaling) {
    // Randomize angles, not to large, so that the camera is in front of the camera.
    Eigen::Matrix3d R;
    double xv = (20.0 * rand() / RAND_MAX - 10.0) * M_PI / 180.0;
    double yv = (20.0 * rand() / RAND_MAX - 10.0) * M_PI / 180.0;
    double zv = (20.0 * rand() / RAND_MAX - 10.0) * M_PI / 180.0;

    std::cout << "x = " << xv << " y = " << yv << " z = " << zv << std::endl;

    R = Eigen::AngleAxisd(xv, Eigen::Vector3d::UnitX())
      * Eigen::AngleAxisd(yv, Eigen::Vector3d::UnitY())
      * Eigen::AngleAxisd(zv, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d t;
    t.setRandom();
    t *= translation_scaling;
    pose->R = R;
    if (pose->R.determinant() < 0)
        pose->R *= -1.0;
    pose->t = t;
}

void generate_scene_and_image(
    int N,
    double min_depth,
    double max_depth,
    double h_fov,
    DronePoseLib::Camera *pose,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points1,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points2,
    double translation_scaling) {

    // Make sure image points are of the correct size
    image_points1->resize(2, N);
    image_points2->resize(2, N);

    // Use horizontal field of view to determine maximal coordinate
    double max_coord = tan(h_fov / 2 * M_PI / 180);

    // Create world points
    Eigen::Matrix<double, 3, Eigen::Dynamic> world_points;
    world_points.resize(3, N);
    world_points.block(0, 0, 2, N).setRandom();
    world_points.block(0, 0, 2, N) *= max_coord;
    world_points.block(2, 0, 1, N).setOnes();

    // Add depth (third coordinate)
    Eigen::Array<double, 1, Eigen::Dynamic> depths(1, N);
    depths.setRandom();
    depths = (max_depth - min_depth) * (depths + 1.0) / 2.0 + min_depth;

    for (int i = 0; i < N; ++i) {
        world_points.col(i) *= depths(i);
    }

    // Get a random (relative) pose
    set_random_pose(pose, translation_scaling);

    // First camera is assumed to be [I 0]
    image_points1->row(0) = world_points.row(0).array() / world_points.row(2).array();
    image_points1->row(1) = world_points.row(1).array() / world_points.row(2).array();

    // Second camera is assumed to be [R t]
    world_points = pose->R * world_points;
    world_points.colwise() += pose->t;
    image_points2->row(0) = world_points.row(0).array() / world_points.row(2).array();
    image_points2->row(1) = world_points.row(1).array() / world_points.row(2).array();
}

/* Adds focal length to the camera and image points.
  Note that the order of add_focal and add_distortion* matters! */
void add_focal(
    double focal,
    DronePoseLib::Camera *pose,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points) {
    pose->focal = focal;
    (*image_points) *= focal;
}

/* Adds 1 param. division model to the camera and image points.
  Note that the order of add_focal and add_distortion* matters! */
void add_distortion(
    double lambda,
    DronePoseLib::Camera *pose,
    Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points) {
    pose->dist_params.clear();
    pose->dist_params.push_back(lambda);
    DronePoseLib::inverse_1param_division_model(lambda, *image_points, image_points);
}

/* Adds random noise to the image points. */
void add_noise(double sigma, Eigen::Matrix<double, 2, Eigen::Dynamic> *image_points) {
    Eigen::Matrix<double, 2, Eigen::Dynamic> noise;
    noise.resizeLike(*image_points);
    // TODO(marcusvaltonen): Normal distribution, but without using distribution, since
    // it is compiler dependent (makes integration tests harder to implement).
    /*
    std::default_random_engine gen;
    std::normal_distribution<double> d(0.0, sigma);
    for (int i = 0; i < noise.cols(); i++) {
        noise(0, i) = d(gen);
        noise(1, i) = d(gen);
    }
    */
    noise.setRandom();
    noise *= sigma;
    *image_points += noise;
}

void debug_print_pose(DronePoseLib::Camera pose, bool relative) {
    Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
    std::cout << "R: " << pose.R.format(HeavyFmt) << std::endl;
    Eigen::Vector3d t = pose.t;

    if (relative) {
        t /= t.norm();
    }

    std::cout << "t:" << t.format(HeavyFmt) << std::endl;
    std::cout << "f:" << std::setprecision(16) << pose.focal << std::endl;
    for (int j = 0; j < pose.dist_params.size(); j++) {
        std::cout << "d[" << j << "]:" << std::setprecision(16) << pose.dist_params[j] << std::endl;
    }
}
