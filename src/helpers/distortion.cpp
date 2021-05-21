#include "distortion.hpp"

using namespace Eigen;

//#include <iostream>

static const double TOL_ITERATIVE = 1e-10;

void DronePoseLib::inverse_1param_division_model(double lambda, const Eigen::Matrix<double, 2, Eigen::Dynamic> &x0, Eigen::Matrix<double, 2, Eigen::Dynamic>* x1) {

	if (lambda == 0.0) {
		*x1 = x0;
		return;
	}


	Array<double, 1, Dynamic> ru2 = x0.colwise().squaredNorm();
	Array<double, 1, Dynamic> ru = ru2.sqrt();

	double dist_sign = 1.0;
	if (lambda < 0)
		dist_sign = -1.0;

	Array<double, 1, Dynamic> rd;
	rd.resizeLike(ru);
	rd = (0.5 / lambda) / ru - dist_sign * ((0.25 / (lambda * lambda)) / ru2 - 1 / lambda).sqrt();
//  rd = 1 / 2 / kappa. / ru - sign(kappa) * sqrt(1 / 4 / kappa ^ 2. / ru2 - 1 / kappa);

	rd /= ru;
	x1->resizeLike(x0);
	x1->row(0) = x0.row(0).array() * rd;
	x1->row(1) = x0.row(1).array() * rd;

}

void DronePoseLib::forward_1param_division_model(double lambda, const Eigen::Matrix<double, 2, Eigen::Dynamic> &x0, Eigen::Matrix<double, 2, Eigen::Dynamic>* x1) {

	Array<double, 1, Dynamic> r2 = x0.colwise().squaredNorm();
	r2 *= lambda;
	r2 += 1.0;

	x1->resizeLike(x0);
	x1->row(0) = x0.row(0).array() / r2;
	x1->row(1) = x0.row(1).array() / r2;
}
