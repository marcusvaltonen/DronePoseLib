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

	static const int MAX_SAMPLE_SIZE = 8;
	typedef Eigen::Matrix<double, 2, Eigen::Dynamic, 0, 2, MAX_SAMPLE_SIZE> Points2D;

	// We use CRTP here for the solvers.
	template<class Solver>
	class PoseEstimator {
	public:
		int estimate(const Points2D &image_points1, const Points2D &image_points2, std::vector<Camera> *poses) const;

		inline int minimal_sample_size() const {
			return static_cast<const Solver*>(this)->minimal_sample_size();
		}

		// Options
		bool normalize_image_coord = false;  // TODO: Think this over
		bool center_world_coord = false;  // TODO: Think this over
		bool normalize_world_coord = false;  // TODO: Think this over
		bool check_chirality = false; // TODO: Think this over
		bool check_reprojection_error = false;
		double reprojection_threshold = 10.0; // in pixels

	protected:
		PoseEstimator() = default;
	private:
		inline void filter_chirality(const Points3D &world_points, std::vector<Camera> &poses) const;
		inline void filter_reprojection_error(const Points2D &image_points, const Points3D &world_points, double tol, std::vector<Camera> &poses) const;
	};
};


template<class Solver>
int DronePoseLib::PoseEstimator<Solver>::estimate(const Points2D &image_points1, const Points2D &image_points2, std::vector<Camera> *poses) const
{
	Points2D x1 = image_points1;
	Points2D x2 = image_points2;

	// Rescale image plane
	double f0 = 1.0;
	if (normalize_image_coord) {
		f0 = x.colwise().norm().mean() / std::sqrt(2.0);
		x /= f0;
	}

	// Rescale and translate world coordinate system
	Eigen::Vector3d t0(0.0, 0.0, 0.0);
	double s0 = 1.0;

	if (center_world_coord) {
		t0 = X.rowwise().mean();
		X.colwise() -= t0;
	}
	if(normalize_world_coord) {
		s0 = X.colwise().norm().mean();
		X /= s0;
	}

	// Call solver implementation
	poses->clear();
	int n_sols = static_cast<const Solver*>(this)->solve(x, X, poses);

	if (check_chirality)
		filter_chirality(X, *poses);

	if (check_reprojection_error)
		filter_reprojection_error(x, X, reprojection_threshold / f0, *poses);

	// Revert image coordinate scaling
	if (normalize_image_coord) {
		for (int i = 0; i < poses->size(); ++i) {
			(*poses)[i].focal *= f0;
		}
	}

	// Revert world coordinate changes
	if (normalize_world_coord) {
		for (int i = 0; i < poses->size(); ++i)
			(*poses)[i].t = (*poses)[i].t * s0;
	}
	if (center_world_coord) {
		for (int i = 0; i < poses->size(); ++i)
			(*poses)[i].t -= (*poses)[i].R * t0;
	}


	return n_sols;
}


template<class Solver>
void DronePoseLib::PoseEstimator<Solver>::filter_chirality(const Points3D& world_points, std::vector<Camera>& poses) const {
	auto any_point_behind = [&world_points](const Camera& p) { return ((p.R.block<1, 3>(2, 0) * world_points).array() + p.t(2)).minCoeff() < 0;  };
	poses.erase(std::remove_if(poses.begin(), poses.end(), any_point_behind), poses.end());
}


template<class Solver>
void DronePoseLib::PoseEstimator<Solver>::filter_reprojection_error(const Points2D& image_points, const Points3D& world_points, double tol, std::vector<Camera>& poses) const {
	Eigen::Matrix<double, 2, Eigen::Dynamic> x;
	Eigen::Matrix<double, 3, Eigen::Dynamic> X;
	x.resizeLike(image_points);
	X.resizeLike(world_points);


	for (auto it = poses.begin(); it != poses.end();) {

		// Computer reprojection error
		X = it->R * world_points;
		X.colwise() += it->t;
		x.row(0) = X.row(0).array() / X.row(2).array();
		x.row(1) = X.row(1).array() / X.row(2).array();
		static_cast<const Solver*>(this)->distort(it->dist_params, x, &x);

		x *= it->focal;

		double max_error = std::sqrt((x - image_points).colwise().squaredNorm().maxCoeff());

		if (max_error < tol) {
			++it;
		} else {
			it = poses.erase(it);
		}
	}
}  // SRC_RANSAC_POSE_ESTIMATOR_HPP_
