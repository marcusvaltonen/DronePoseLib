#include <Eigen/Dense>
#include <chrono>  // NOLINT [build/c++11]
#include <iostream>
#include "get_valtonenornhag_arxiv_2021.hpp"
#include "relpose.hpp"


int main() {
    /* Timing experiments */
    int nbr_iter = 1e4;

    // Test fEf
    int N = 3;
    Eigen::MatrixXd x1(2, N);
    x1 = Eigen::MatrixXd::Random(2, N);
    Eigen::MatrixXd x2(2, N);
    x2 = Eigen::MatrixXd::Random(2, N);
    Eigen::Matrix3d R1;
    R1 = Eigen::Matrix3d::Random(3, 3);
    Eigen::Matrix3d R2;
    R2 = Eigen::Matrix3d::Random(3, 3);

    auto start = std::chrono::steady_clock::now();
    std::cout << "test" << std::endl;

    std::vector<DronePoseLib::RelPose> poses;

    for (int i = 0; i < nbr_iter; i++) {
        poses = DronePoseLib::ValtonenOrnhagArxiv2021::get_fEf(x1, x2, R1, R2);
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time (fEf): "
        << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / nbr_iter
        << " ns" << std::endl;

    // Test frEfr
    N = 4;
    x1 = Eigen::MatrixXd::Random(2, N);
    x2 = Eigen::MatrixXd::Random(2, N);
    R1 = Eigen::Matrix3d::Random(3, 3);
    R2 = Eigen::Matrix3d::Random(3, 3);

    start = std::chrono::steady_clock::now();
    for (int i = 0; i < nbr_iter; i++) {
        poses = DronePoseLib::ValtonenOrnhagArxiv2021::get_frEfr(x1, x2, R1, R2);
    }
    end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time (frEfr): "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / nbr_iter
        << " μs" << std::endl;

    return 0;
}
