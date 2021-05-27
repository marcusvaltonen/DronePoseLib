//TODO: License here

#include "refinement.hpp"

using namespace Eigen;
using namespace DronePoseLib;

// TODO: Probably do not want these static declared..
static const double SMALL_NUMBER = 1e-8;
static const double TOL_CONVERGENCE = 1e-10;
static const double INITIAL_LM_DAMP = 1e-6;
static const int MAX_ITER = 10;


/* NOTE: These are Viktor's notes!
 TODO: add check that cost decreases in LM
 TODO: move radial refinement from larsson_iccv19.cc to here
 TODO: r_pow[] thing for refinement_dist as well...
*/

#include <iostream>  // DEBUG

inline void drot(const Matrix3d &R, Matrix3d *dr1, Matrix3d *dr2, Matrix3d *dr3) {
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

inline void update_rot(Vector3d &v, Matrix3d &rot) {
    double stheta = v.norm();
    if (stheta < SMALL_NUMBER)
        return;
    v /= stheta;
    double theta = asin(stheta);
    Matrix3d K;
    K << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
    Matrix3d deltaR = Matrix3d::Identity() + stheta * K + (1 - cos(theta))*K*K;
    rot = deltaR * rot;
}

void DronePoseLib::refinement_dist_with_structure(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x1,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x2,
    Eigen::Matrix<double, 3, Eigen::Dynamic> &X,
    DronePoseLib::Camera &p) {

    int n_pts = x1.cols();
    int n_res = 2 * 2 * n_pts;
    int n_params = 6 + 1 + 1 + 3 * n_pts;
    double lm_damp = INITIAL_LM_DAMP;

    Matrix3d dr1, dr2, dr3;

    // Order for jacobian is: rotation, translation, focal, dist_params[0], 3D-points
    MatrixXd J(n_res, n_params);
    J.setZero();
    VectorXd res(n_res, 1);
    VectorXd dx(n_params, 1);
    res.setZero();

    Matrix<double, 3, Dynamic> Z = X;

    MatrixXd H;
    VectorXd g;

    int iter;
    for (iter = 0; iter < MAX_ITER; ++iter) {
        // Z = R * X + t
        Z = p.R * X;
        Z.colwise() += p.t;

        drot(p.R, &dr1, &dr2, &dr3);

        for (int i = 0; i < n_pts; ++i) {
            double d1 = X(2, i);
            double d2 = Z(2, i);
            double d12 = d1 * d1;
            double d22 = d2 * d2;
            double r12 = X.block<2, 1>(0, i).squaredNorm();
            double r22 = Z.block<2, 1>(0, i).squaredNorm();

            double num = p.focal;
            double dnum_dt3 = 0.0;

            double r12d12 = r12 / d12;
            double r22d22 = r22 / d22;
            double denom1 = d1 + d1 * p.dist_params[0] * r12d12;
            double denom2 = d2 + d2 * p.dist_params[0] * r22d22;
            double ddenom_dt3 = 1.0 - p.dist_params[0] * r22d22;

            double factor1 = num / denom1;
            double factor2 = num / denom2;
            double dfactor_dt3 = (dnum_dt3 * denom2 - num * ddenom_dt3) / (denom2 * denom2);
            double dfactor_df1 = 1.0 / denom1;
            double dfactor_df2 = 1.0 / denom2;

            // Residual in distorted space
            std::cout << "factor1 = " << factor1 << std::endl;
            std::cout << "X = \n" << X.col(i) << std::endl;
            std::cout << "x1 = \n" << x1.col(i) << std::endl;
            std::cout << "factor2 = " << factor2 << std::endl;
            std::cout << "Z = \n" << Z.col(i) << std::endl;
            std::cout << "x2 = \n" << x2.col(i) << std::endl;

            res(4 * i + 0) = factor1 * X(0, i) - x1(0, i);
            res(4 * i + 1) = factor1 * X(1, i) - x1(1, i);
            res(4 * i + 2) = factor2 * Z(0, i) - x2(0, i);
            res(4 * i + 3) = factor2 * Z(1, i) - x2(1, i);

            // rotation
            // dfactor_dz = [0 0 dfactor_dt3]
            // dfactor_dr = dfactor_dz * dz_dr
            // d(factor*Z(0,i))_dr = dfactor_dr*Z(0,i) + factor * dZ(0,i)_dr

            // dZ_dr1
            Vector3d dZ_dr1 = dr1 * X.col(i);
            Vector3d dZ_dr2 = dr2 * X.col(i);
            Vector3d dZ_dr3 = dr3 * X.col(i);

            double dfactor_dr1 = dfactor_dt3 * dZ_dr1(2);
            double dfactor_dr2 = dfactor_dt3 * dZ_dr2(2);
            double dfactor_dr3 = dfactor_dt3 * dZ_dr3(2);

            J(4 * i + 2, 0) = dfactor_dr1 * Z(0, i) + factor2 * dZ_dr1(0);
            J(4 * i + 2, 1) = dfactor_dr2 * Z(0, i) + factor2 * dZ_dr2(0);
            J(4 * i + 2, 2) = dfactor_dr3 * Z(0, i) + factor2 * dZ_dr3(0);

            J(4 * i + 3, 0) = dfactor_dr1 * Z(1, i) + factor2 * dZ_dr1(1);
            J(4 * i + 3, 1) = dfactor_dr2 * Z(1, i) + factor2 * dZ_dr2(1);
            J(4 * i + 3, 2) = dfactor_dr3 * Z(1, i) + factor2 * dZ_dr3(1);

            // t_x
            J(4 * i + 2, 3) = factor2;

            // t_y
            J(4 * i + 3, 4) = factor2;

            // t_z
            J(4 * i + 2, 5) = dfactor_dt3 * Z(0, i);
            J(4 * i + 3, 5) = dfactor_dt3 * Z(1, i);

            // focal
            J(4 * i + 0, 6) = dfactor_df1 * X(0, i);
            J(4 * i + 1, 6) = dfactor_df1 * X(1, i);
            J(4 * i + 2, 6) = dfactor_df2 * Z(0, i);
            J(4 * i + 3, 6) = dfactor_df2 * Z(1, i);

            // dist_params[0]
            double dfactor_dlambda1 = -d1 * r12d12 * num / (denom1*denom1);
            double dfactor_dlambda2 = -d2 * r22d22 * num / (denom2*denom2);
            J(4 * i + 1, 7) = dfactor_dlambda1 * X(0, i);
            J(4 * i + 1, 7) = dfactor_dlambda1 * X(1, i);
            J(4 * i + 2, 7) = dfactor_dlambda2 * Z(0, i);
            J(4 * i + 3, 7) = dfactor_dlambda2 * Z(1, i);

            // 3D-points
            Matrix3d R = p.R;
            Vector3d Xi = X.col(i);
            Vector3d t = p.t;
            double f = p.focal;
            double lam = p.dist_params[0];
            double X2i2 = X(2, i) * X(2, i);
            double Z2i2 = Z(2, i) * Z(2, i);

        // TODO: Redo these computations with pen and paper to make it more legible.... chain-rule style
            double s1 = 1.0 / (1.0 + lam*(Z(0, i)*Z(0, i)/Z2i2 + Z(1, i)*Z(1, i)/Z2i2));
            double s2 = 1 / Z(2, i);  // "1/(t(2) + R(2, 0)*Xi(0) + R(2, 1)*Xi(1) + R(2, 2)*Xi(2))"
            double s3 = Z(1, i);  // "t(1) + R(1, 0)*Xi(0) + R(1, 1)*Xi(1) + R(1, 2)*Xi(2)";
            double s4 = Z(0, i); // "t(0) + R(0, 0)*Xi(0) + R(0, 1)*Xi(1) + R(0, 2)*Xi(2)";
            double s5 = s2 * s2;
            double s6 = s1 * s1;
            double s7 = s2*s5*s3*s3;
            double s8 = s2*s5*s4*s4;

            double t1 = 1.0 / (1.0 + lam * (X(0, i)*X(0, i)/X2i2 + X(1, i)*X(1, i)/X2i2));
            double t2 = 1 / X2i2;
            double t3 = 2*t2*(X(0, i)*X(0, i) + X(1, i)*X(1, i));
            double t4 = 1 / X(2, i);
            double t5 = -2*X(0, i)*X(1, i)*f*lam*X(0, i)*X(0, i)*t2;
            double t6 = t1 * t1;

            // X(0)
            J(4 * i + 0, 8 + 3 * i + 0) = f*t(0)*t4 - 2*Xi(0)*Xi(0)*f*lam*t(1)*t6;
            J(4 * i + 1, 8 + 3 * i + 0) = t5;
            J(4 * i + 2, 8 + 3 * i + 0) = R(0, 0)*f*s1*s2 - R(2, 0)*f*s1*s4*s5 + f*lam*s6*s2*s4*(2*R(2, 0)*s7 + 2*R(2, 0)*s8 - 2*R(0, 0)*s4*s5 - 2*R(1, 0)*s3*s5);
            J(4 * i + 3, 8 + 3 * i + 0) = R(1, 0)*f*s1*s2 - R(2, 0)*f*s1*s3*s5 + f*lam*s6*s2*s3*(2*R(2, 0)*s7 + 2*R(2, 0)*s8 - 2*R(0, 0)*s4*s5 - 2*R(1, 0)*s3*s5);

            // X(1)
            J(4 * i + 0, 8 + 3 * i + 1) = t5;
            J(4 * i + 1, 8 + 3 * i + 1) = f*t(0)*t4 - 2*Xi(1)*Xi(1)*f*lam*t(1)*t6;
            J(4 * i + 2, 8 + 3 * i + 1) = R(0, 1)*f*s1*s2 - R(2, 1)*f*s1*s4*s5 + f*lam*s6*s2*s4*(2*R(2, 1)*s7 + 2*R(2, 1)*s8 - 2*R(0, 1)*s4*s5 - 2*R(1, 1)*s3*s5);
            J(4 * i + 3, 8 + 3 * i + 1) = R(1, 1)*f*s1*s2 - R(2, 1)*f*s1*s3*s5 + f*lam*s6*s2*s3*(2*R(2, 1)*s7 + 2*R(2, 1)*s8 - 2*R(0, 1)*s4*s5 - 2*R(1, 1)*s3*s5);

            // X(2)
            J(4 * i + 0, 8 + 3 * i + 2) = Xi(0)*f*lam*t(2)*t4*t6 - Xi(0)*f*t(0)*t4*t4;
            J(4 * i + 1, 8 + 3 * i + 2) = Xi(1)*f*lam*t(2)*t4*t6 - Xi(1)*f*t(0)*t4*t4;
            J(4 * i + 2, 8 + 3 * i + 2) = R(0, 2)*f*s1*s2 - R(2, 2)*f*s1*s4*s5 + f*lam*s6*s2*s4*(2*R(2, 2)*s7 + 2*R(2, 2)*s8 - 2*R(0, 2)*s4*s5 - 2*R(1, 2)*s3*s5);
            J(4 * i + 3, 8 + 3 * i + 2) = R(1, 2)*f*s1*s2 - R(2, 2)*f*s1*s3*s5 + f*lam*s6*s2*s3*(2*R(2, 2)*s7 + 2*R(2, 2)*s8 - 2*R(0, 2)*s4*s5 - 2*R(1, 2)*s3*s5);
        }

        if (res.norm() < TOL_CONVERGENCE)
            break;

        H = J.transpose()*J;
        H.diagonal().array() += lm_damp; // LM dampening
        g = -J.transpose()*res;


        if (g.cwiseAbs().maxCoeff() < TOL_CONVERGENCE)
            break;

        std::cout << "iter = " << iter << " res = " << res.squaredNorm() << ", g = "<< g.squaredNorm() << "\n";
        std::cout << "res =\n" << res << "\n";

        dx = H.ldlt().solve(g);

        // std::cout << "dx =\n" << dx << "\n";

        Vector3d dx_r = dx.head(3);
        update_rot(dx_r, p.R);
        p.t(0) += dx(3);
        p.t(1) += dx(4);
        p.t(2) += dx(5);
        p.focal += dx(6);
        p.dist_params[0] += dx(7);

        for (int k = 0; k < n_pts; k++) {
            X.col(k) += dx.segment(8 + 3 * k, 3);
        }

        if (dx.array().abs().maxCoeff() < SMALL_NUMBER)
            break;
        lm_damp = std::max(1e-8, lm_damp / 10.0);
    }
     std::cout << "Local opt finished. iter=" << iter << ", res=" << res.norm() << ", g=" << g.norm() << "\n";
}

void DronePoseLib::refinement_undist_with_structure(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x1,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x2,
    Eigen::Matrix<double, 3, Eigen::Dynamic> &X,
    DronePoseLib::Camera &p) {
    int n_pts = x2.cols();

    // One for each view
    int n_res = 2 * 2 * n_pts;

    // First camera is fixed P1 = [I 0]
    int n_params = 6 + 1 + 1 + 3 * n_pts;

    Matrix3d dr1, dr2, dr3;
    double lm_damp = INITIAL_LM_DAMP;

    // Order for jacobian is: rotation, translation, focal, dist_params[0], 3D-points
    Matrix<double, Dynamic, Dynamic> J(n_res, n_params);
    J.setZero();
    VectorXd res(n_res, 1);
    VectorXd dx(n_params, 1);
    res.setZero();

    Matrix<double, 3, Dynamic> Z = X;

    // Change of variables to simplify equations
    p.dist_params[0] /= p.focal;

    MatrixXd H;
    VectorXd g;
    int iter;
    for (iter = 0; iter < MAX_ITER; ++iter) {
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

            Vector3d dZ_dr1 = dr1 * X.col(i);
            Vector3d dZ_dr2 = dr2 * X.col(i);
            Vector3d dZ_dr3 = dr3 * X.col(i);

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
            Matrix3d R = p.R;
            Vector3d Xi = X.col(i);

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

        if (res.norm() < TOL_CONVERGENCE)
            break;

        H = J.transpose() * J;
        H.diagonal().array() += lm_damp; // LM dampening
        g = -J.transpose() * res;


        if (g.cwiseAbs().maxCoeff() < TOL_CONVERGENCE)
            break;

        std::cout << "iter=" << iter << " res.squaredNorm()=" << res.squaredNorm() << ", g="<< g.squaredNorm() << "\n";
        // std::cout << "res =\n" << res << std::endl;
        // std::cout << "J =\n" << J << "\n";

        dx = H.ldlt().solve(g);
        // std::cout << "dx =\n"<< dx << "\n";

        Vector3d dx_r = dx.head(3);
        update_rot(dx_r, p.R);
        p.t(0) += dx(3);
        p.t(1) += dx(4);
        p.t(2) += dx(5);
        p.focal += dx(6);
        p.dist_params[0] += dx(7);

        for (int k = 0; k < n_pts; k++) {
            X.col(k) += dx.segment(8 + 3 * k, 3);
        }

        if (dx.array().abs().maxCoeff() < SMALL_NUMBER)
            break;

        lm_damp = std::max(1e-8, lm_damp / 10.0);
    }

    // std::cout << "Local opt finished. iter=" << iter << ", res=" << res.norm() << ", g=" << g.norm() << "\n";

    // Revert change of variables
    p.dist_params[0] *= p.focal;
}
