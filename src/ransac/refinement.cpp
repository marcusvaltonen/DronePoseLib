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
// Modified from Viktor Larsson's radialpose implementation
// https://github.com/vlarsson/radialpose

#include "refinement.hpp"
#include <Eigen/Dense>


inline void drot(const Eigen::Matrix3d &R, Eigen::Matrix3d *dr1, Eigen::Matrix3d *dr2, Eigen::Matrix3d *dr3) {
    // skew = [0 -v(2) v(1); v(2) 0 -v(0); -v(1) v(0) 0]

    // dr1 = [0 0 0; 0 0 -1; 0 1 0]*R
    dr1->row(0).setZero();
    dr1->row(1) = -R.row(2);
    dr1->row(2) = R.row(1);

    // dr2 = [0 0 1; 0 0 0; -1 0 0]*R
    dr2->row(0) = R.row(2);
    dr2->row(1).setZero();
    dr2->row(2) = -R.row(0);

    // dr3 = [0 -1 0; 1 0 0; 0 0 0]*R
    dr3->row(0) = -R.row(1);
    dr3->row(1) = R.row(0);
    dr3->row(2).setZero();
}

inline void update_rot(Eigen::Vector3d &v, Eigen::Matrix3d &rot, const DronePoseLib::RefinementSettings &settings) {
    double stheta = v.norm();
    if (stheta < settings.SMALL_NUMBER)
        return;
    v /= stheta;
    double theta = asin(stheta);
    Eigen::Matrix3d K;
    K << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
    Eigen::Matrix3d deltaR = Eigen::Matrix3d::Identity() + stheta * K + (1 - cos(theta))*K*K;
    rot = deltaR * rot;
}


// TODO: Consider returning a flag on success, break etc.
void DronePoseLib::refinement_undist_with_structure(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x1,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x2,
    Eigen::Matrix<double, 3, Eigen::Dynamic> &X,
    DronePoseLib::Camera &p,
    const DronePoseLib::RefinementSettings &settings) {
    int n_pts = x2.cols();

    // One for each view
    int n_res = 2 * 2 * n_pts;

    // First camera is fixed P1 = [I 0]
    int n_params = 6 + 1 + 1 + 3 * n_pts;

    Eigen::Matrix3d dr1, dr2, dr3;
    double lm_damp = settings.INITIAL_LM_DAMP;

    // Order for Jacobian is: rotation, translation, focal, dist_params[0], 3D-points
    Eigen::MatrixXd J(n_res, n_params);
    J.setZero();
    Eigen::VectorXd res(n_res, 1);
    Eigen::VectorXd dx(n_params, 1);
    res.setZero();

    Eigen::Matrix<double, 3, Eigen::Dynamic> Z = X;

    // Change of variables to simplify equations
    p.dist_params[0] /= p.focal;

    Eigen::MatrixXd H;
    Eigen::VectorXd g;
    int iter;
    for (iter = 0; iter < settings.MAX_ITER; ++iter) {
        // Z = R*X + t
        Z = p.R * X;
        Z.colwise() += p.t;

        drot(p.R, &dr1, &dr2, &dr3);

        for (int i = 0; i < n_pts; ++i) {
            double d = Z(2, i);
            double r12 = x1.col(i).squaredNorm();
            double r22 = x2.col(i).squaredNorm();

            double denom1 = p.focal + p.dist_params[0] * r12;
            double denom2 = p.focal + p.dist_params[0] * r22;

            double factor1 = 1.0 / denom1;
            double factor2 = 1.0 / denom2;

            double dfactor_df1 = -factor1 / denom1;
            double dfactor_df2 = -factor2 / denom2;

            // Residuals (undistorted space) - (x, y) coord for first cam, (x, y) coord for second cam
            res(4 * i + 0) = factor1 * x1(0, i) - X(0, i) / X(2, i);
            res(4 * i + 1) = factor1 * x1(1, i) - X(1, i) / X(2, i);
            res(4 * i + 2) = factor2 * x2(0, i) - Z(0, i) / d;
            res(4 * i + 3) = factor2 * x2(1, i) - Z(1, i) / d;

            // Rotation
            // d(-Z(1)/Z(3)) = dZ(1)/Z(3) - Z(1)*dZ(3)/Z(3)^2
            double Z2i2 = Z(2, i) * Z(2, i);

            Eigen::Vector3d dZ_dr1 = dr1 * X.col(i);
            Eigen::Vector3d dZ_dr2 = dr2 * X.col(i);
            Eigen::Vector3d dZ_dr3 = dr3 * X.col(i);

            J(4 * i + 2, 0) = - dZ_dr1(0) / Z(2, i) + Z(0, i) * dZ_dr1(2) / Z2i2;
            J(4 * i + 2, 1) = - dZ_dr2(0) / Z(2, i) + Z(0, i) * dZ_dr2(2) / Z2i2;
            J(4 * i + 2, 2) = - dZ_dr3(0) / Z(2, i) + Z(0, i) * dZ_dr3(2) / Z2i2;

            J(4 * i + 3, 0) = -dZ_dr1(1) / Z(2, i) + Z(1, i) * dZ_dr1(2) / Z2i2;
            J(4 * i + 3, 1) = -dZ_dr2(1) / Z(2, i) + Z(1, i) * dZ_dr2(2) / Z2i2;
            J(4 * i + 3, 2) = -dZ_dr3(1) / Z(2, i) + Z(1, i) * dZ_dr3(2) / Z2i2;

            // tx
            J(4 * i + 2, 3) = - 1.0 / d;

            // ty
            J(4 * i + 3, 4) = - 1.0 / d;

            // tz
            J(4 * i + 2, 5) = Z(0, i) / (d*d);
            J(4 * i + 3, 5) = Z(1, i) / (d*d);

            // focal
            J(4 * i + 0, 6) = dfactor_df1 * x1(0, i);
            J(4 * i + 1, 6) = dfactor_df1 * x1(1, i);
            J(4 * i + 2, 6) = dfactor_df2 * x2(0, i);
            J(4 * i + 3, 6) = dfactor_df2 * x2(1, i);

            // dist_params[0]
            double dfactor_dlambda1 = -r12 / (denom1 * denom1);
            double dfactor_dlambda2 = -r22 / (denom2 * denom2);
            J(4 * i + 0, 7) = dfactor_dlambda1 * x1(0, i);
            J(4 * i + 1, 7) = dfactor_dlambda1 * x1(1, i);
            J(4 * i + 2, 7) = dfactor_dlambda2 * x2(0, i);
            J(4 * i + 3, 7) = dfactor_dlambda2 * x2(1, i);

            // 3D-points
            Eigen::Matrix3d R = p.R;
            Eigen::Vector3d Xi = X.col(i);

            // X(0)
            J(4 * i + 0, 8 + 3 * i + 0) = -1 / Xi(2);
            J(4 * i + 1, 8 + 3 * i + 0) = 0;
            J(4 * i + 2, 8 + 3 * i + 0) = -R(0, 0) / Z(2, i) + R(2, 0) * Z(0, i) / Z2i2;
            J(4 * i + 3, 8 + 3 * i + 0) = -R(1, 0) / Z(2, i) + R(2, 0) * Z(1, i) / Z2i2;

            // X(1)
            J(4 * i + 0, 8 + 3 * i + 1) = 0;
            J(4 * i + 1, 8 + 3 * i + 1) = -1 / Xi(2);
            J(4 * i + 2, 8 + 3 * i + 1) = -R(0, 1) / Z(2, i) + R(2, 1) * Z(0, i) / Z2i2;
            J(4 * i + 3, 8 + 3 * i + 1) = -R(1, 1) / Z(2, i) + R(2, 1) * Z(1, i) / Z2i2;

            // X(2)
            J(4 * i + 0, 8 + 3 * i + 2) = Xi(0) / (Xi(2) * Xi(2));
            J(4 * i + 1, 8 + 3 * i + 2) = Xi(1) / (Xi(2) * Xi(2));
            J(4 * i + 2, 8 + 3 * i + 2) = -R(0, 2) / Z(2, i) + R(2, 2) * Z(0, i) / Z2i2;
            J(4 * i + 3, 8 + 3 * i + 2) = -R(1, 2) / Z(2, i) + R(2, 2) * Z(1, i) / Z2i2;
        }

        if (res.norm() < settings.TOL_CONVERGENCE)
            break;

        H = J.transpose() * J;
        H.diagonal().array() += lm_damp;  // LM dampening
        g = -J.transpose() * res;


        if (g.cwiseAbs().maxCoeff() < settings.TOL_CONVERGENCE)
            break;

        dx = H.ldlt().solve(g);

        // TODO: add check that cost decreases in LM
        Eigen::Vector3d dx_r = dx.head(3);
        update_rot(dx_r, p.R, settings);
        p.t(0) += dx(3);
        p.t(1) += dx(4);
        p.t(2) += dx(5);
        p.focal += dx(6);
        p.dist_params[0] += dx(7);

        for (int k = 0; k < n_pts; k++) {
            X.col(k) += dx.segment(8 + 3 * k, 3);
        }

        if (dx.array().abs().maxCoeff() < settings.SMALL_NUMBER)
            break;

        lm_damp = std::max(1e-8, lm_damp / settings.DECREASE_FACTOR);
    }

    // Revert change of variables
    p.dist_params[0] *= p.focal;
}
