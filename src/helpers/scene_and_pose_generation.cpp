#define _USE_MATH_DEFINES
#include "scene_and_pose_generation.hpp"
#include "relpose.hpp"
#include <Eigen/Dense>
#include "distortion.hpp"

//#define FLIPPYBOOL(boolname) bool boolname = false; for(int o_##boolname = 0; o_##boolname < 2; ++o_##boolname, boolname = true)

static const double TOL_POSE = 1.0e-6;

// Passing tests return true.
//#define TEST(FUNC) if(!FUNC()) { std::cout << #FUNC"\033[1m\033[31m FAILED!\033[0m\n"; } else { std::cout << #FUNC"\033[1m\033[32m PASSED!\033[0m\n"; passed++;} num_tests++;

#include <iostream>

using namespace Eigen;
using namespace DronePoseLib;
using namespace std;


double pose_distance(Camera ref, Camera pose, bool compare_focal, bool compare_distortion) {

	double dist_R = (ref.R - pose.R).norm();
	double dist_t = (ref.t - pose.t).norm();
	double dist_f = std::abs(ref.focal - pose.focal) / std::abs(ref.focal);
	double dist_d = 0.0;

	if (ref.dist_params.size() != pose.dist_params.size()) {
		dist_d = std::numeric_limits<double>::infinity();
	} else {
		for (int i = 0; i < ref.dist_params.size(); ++i) {
			dist_d += std::abs(ref.dist_params[i] - pose.dist_params[i]);
		}
	}

	double dist = dist_R + dist_t;
	if (compare_focal)
		dist += dist_f;
	if (compare_distortion)
		dist += dist_d;

	return dist;
}

double minimum_pose_distance(Camera pose_gt, std::vector<Camera> poses, bool compare_focal, bool compare_distortion) {
	double min_dist = std::numeric_limits<double>::max();

	for (int i = 0; i < poses.size(); ++i) {
		double dist = pose_distance(pose_gt, poses[i], compare_focal, compare_distortion);
		min_dist = std::min(min_dist, dist);
	}
	return min_dist;
}


void set_random_pose(Camera* pose, double translation_scaling) {
	//Quaterniond qq = Quaterniond::UnitRandom();
    Matrix3d R;
    // TODO: Randomize
    R = Eigen::AngleAxisd(0.025*M_PI, Eigen::Vector3d::UnitX())
      * Eigen::AngleAxisd(0.05*M_PI,  Eigen::Vector3d::UnitY())
      * Eigen::AngleAxisd(0.033*M_PI, Eigen::Vector3d::UnitZ());
	Vector3d t;
	t.setRandom();
	t *= translation_scaling;
	//pose->R = qq.toRotationMatrix();
	pose->R = R;
	if (pose->R.determinant() < 0)
		pose->R *= -1.0;
	pose->t = t;
}

void project_3d_points_to_plane(Matrix<double, 3, Dynamic> *X) {

	Vector3d t = X->rowwise().mean();

	// subtract mean
	X->colwise() -= t;

	// project to rank 2
	JacobiSVD<Matrix<double, 3, Dynamic>> svd(*X, ComputeThinU | ComputeThinV);
	Matrix<double, 3, 3> D = svd.singularValues().asDiagonal();
	D(2, 2) = 0.0;
	*X = svd.matrixU() * D * svd.matrixV().transpose();

	// add mean again
	X->colwise() += t;
}

void generate_scene_and_image(int N, double min_depth, double max_depth, double h_fov, bool planar, Camera* pose, Matrix<double, 2, Dynamic>* image_points1, Matrix<double, 2, Dynamic>* image_points2, double translation_scaling) {

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
	Array<double, 1, Dynamic> depths(1, N);
	depths.setRandom();
	depths = (max_depth - min_depth) * (depths + 1.0) / 2.0 + min_depth;

	for (int i = 0; i < N; ++i) {
		world_points.col(i) *= depths(i);
	}

    /*
	if (planar) {
		project_3d_points_to_plane(world_points);

		// reproject the now planar 3D points
		image_points1->row(0) = world_points.row(0).array() / world_points.row(2).array();
		image_points1->row(1) = world_points.row(1).array() / world_points.row(2).array();
	}
    */

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
void add_focal(double focal, Camera* pose, Matrix<double, 2, Dynamic>* image_points) {
	pose->focal = focal;
	(*image_points) *= focal;
}

/* Adds 1 param. division model to the camera and image points.
  Note that the order of add_focal and add_distortion* matters! */
void add_distortion_1pdiv(double lambda, Camera* pose, Matrix<double, 2, Dynamic>* image_points) {
	pose->dist_params.clear();
	pose->dist_params.push_back(lambda);
	inverse_1param_division_model(lambda, *image_points, image_points);
}

/* Adds random noise to the image points. */
void add_noise(double sigma, Matrix<double, 2, Dynamic>* image_points) {

	Matrix<double, 2, Dynamic> noise;
	noise.resizeLike(*image_points);
	noise.setRandom(); // TODO Normal distribution instead...
	noise *= sigma;

	*image_points += noise;
}

void debug_print_poses(Camera pose_gt, std::vector<Camera> poses) {
    std::cout << "GROUND TRUTH:\n";
    debug_print_pose(pose_gt);
	for (int i = 0; i < poses.size(); ++i) {
		std::cout << "----- POSE no. " << i << std::endl;
        debug_print_pose(poses[i]);
	}
}

void debug_print_pose(Camera pose, bool relative) {
	std::cout << "R: " << pose.R << "\n";
    Vector3d t = pose.t;

    if (relative) {
        t /= t.norm();
    }

	std::cout << "t:" << t << "\n";
	std::cout << "f:" << pose.focal << "\n";
	for (int j = 0; j < pose.dist_params.size(); j++) {
		std::cout << "d[" << j << "]:" << pose.dist_params[j] << "\n";
	}
}