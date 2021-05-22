#pragma once
#define _USE_MATH_DEFINES
#include "relpose.hpp"
#include <Eigen/Dense>
#include "distortion.hpp"

#define FLIPPYBOOL(boolname) bool boolname = false; for(int o_##boolname = 0; o_##boolname < 2; ++o_##boolname, boolname = true)

//static const double TOL_POSE = 1.0e-6;

// Passing tests return true.
#define TEST(FUNC) if(!FUNC()) { std::cout << #FUNC"\033[1m\033[31m FAILED!\033[0m\n"; } else { std::cout << #FUNC"\033[1m\033[32m PASSED!\033[0m\n"; passed++;} num_tests++;

//TODO: Cleanup using namespace
using namespace Eigen;
using namespace DronePoseLib;


double pose_distance(Camera ref, Camera pose, bool compare_focal = true, bool compare_distortion = true);

double minimum_pose_distance(Camera pose_gt, std::vector<Camera> poses, bool compare_focal = true, bool compare_distortion = true);

void set_random_pose(Camera* pose, double translation_scaling = 1.0);

void project_3d_points_to_plane(Matrix<double, 3, Dynamic> *X);

void generate_scene_and_image(int N, double min_depth, double max_depth, double h_fov, bool planar, Camera* pose, Matrix<double, 2, Dynamic>* image_points1, Matrix<double, 2, Dynamic>* image_points2, double translation_scaling = 1.0);

/* Adds focal length to the camera and image points.
  Note that the order of add_focal and add_distortion* matters! */
void add_focal(double focal, Camera* pose, Matrix<double, 2, Dynamic>* image_points);

/* Adds 1 param. division model to the camera and image points.
  Note that the order of add_focal and add_distortion* matters! */
void add_distortion_1pdiv(double lambda, Camera* pose, Matrix<double, 2, Dynamic>* image_points);

/* Adds random noise to the image points. */
void add_noise(double sigma, Matrix<double, 2, Dynamic>* image_points);

void debug_print_poses(Camera pose_gt, std::vector<Camera> poses);
void debug_print_pose(Camera pose);
