// TODO: License

// TODO: Use header guards instead
#pragma once
#include <Eigen/Dense>
#include "relpose.hpp"

namespace DronePoseLib {
void refinement_dist(const Eigen::Matrix<double, 2, Eigen::Dynamic> &x, const Eigen::Matrix<double, 3, Eigen::Dynamic> &X, DronePoseLib::Camera &p, int Np, int Nd);
void refinement_undist(const Eigen::Matrix<double, 2, Eigen::Dynamic> &x, const Eigen::Matrix<double, 3, Eigen::Dynamic> &X, DronePoseLib::Camera &p, int Np, int Nd);
void refinement_undist_with_structure(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x1,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x2,
    const Eigen::Matrix<double, 3, Eigen::Dynamic> &X,
    DronePoseLib::Camera &p,
    int Np,
    int Nd
);
}
