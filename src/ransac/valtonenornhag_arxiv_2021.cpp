#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <complex>

static const double SMALL_NUMBER = 1e-8;
static const double DAMP_FACTOR = 1e-8;

template<bool unknown_focal, bool unknown_radial_distortion_coeff>
int DronePoseLib::ValtonenOrnhagArxiv2021::Solver<unknown_focal, unknown_radial_distortion_coeff>::solve(const Points2D& image_points1, const Points2D& image_points2, std::vector<Camera>* poses) const
{

	std::vector<Camera> initial_poses;
	std::vector<double> t3;

	if (use_radial_solver) {
		kukelova_iccv13::Radial1DSolver::p5p_radial_impl(image_points, world_points, &initial_poses);
	} else {
		initial_poses.push_back(Camera(Matrix3d::Identity(), Vector3d::Zero()));
	}
	Matrix<double, 3, Dynamic> X;

	for (int k = 0; k < initial_poses.size(); k++) {
		t3.clear();

		X = initial_poses[k].R * world_points;
		X.colwise() += initial_poses[k].t;

		//std::cout << "Initial pose " << k + 1 << "/" << initial_poses.size() << "\n";
		//std::cout << "R=" << initial_poses[k].R << "\nt=" << initial_poses[k].t << "\n";
		//std::cout << "X=" << X << "\n";

		double t0 = 0;
		if (use_precond) {
			t0 = simple_preconditioner(image_points, X);
			X.row(2).array() += t0;
		}

		solver_impl(image_points, X, &t3);

		for (int i = 0; i < t3.size(); ++i) {
			Camera pose;
			pose.R = initial_poses[k].R;
			pose.t = initial_poses[k].t;
			pose.t(2) = t3[i];

			if (DistortionModel) {
				linsolve_known_pose_dist(image_points, X, t3[i], Np, Nd, damp_factor, &pose);
				if(root_refinement)
					radial_refinement_dist<Np, Nd>(pose, image_points, X, damp_factor);
			} else {
				linsolve_known_pose_undist(image_points, X, t3[i], Np, Nd, damp_factor, &pose);
				if(root_refinement)
					radial_refinement_undist<Np, Nd>(pose, image_points, X, damp_factor);
			}

			if (pose.focal < 0) {
				// flipped solution
				pose.focal = -pose.focal;
				pose.R.row(0) = -pose.R.row(0);
				pose.R.row(1) = -pose.R.row(1);
				pose.t(0) = -pose.t(0);
				pose.t(1) = -pose.t(1);
			}


			// Revert precond
			pose.t(2) += t0;

			//std::cout << "solution[" << i << "], t3=" << pose.t(2) << "\n";
			poses->push_back(pose);
		}
	}
	return poses->size();
}

// Template instantiations
template class DronePoseLib::ValtonenOrnhagArxiv2021::Solver<true, true>;
template class DronePoseLib::ValtonenOrnhagArxiv2021::Solver<true, false>;
template class DronePoseLib::ValtonenOrnhagArxiv2021::Solver<false, true>;

template class DronePoseLib::PoseEstimator<DronePoseLib::ValtonenOrnhagArxiv2021::Solver<true, true>>;
template class DronePoseLib::PoseEstimator<DronePoseLib::ValtonenOrnhagArxiv2021::Solver<true, false>>;
template class DronePoseLib::PoseEstimator<DronePoseLib::ValtonenOrnhagArxiv2021::Solver<false, true>>;
