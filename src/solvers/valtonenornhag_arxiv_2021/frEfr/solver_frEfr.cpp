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

#include "solver_frEfr.hpp"
#include <Eigen/Dense>
#include "sturm.h"
#include "charpoly.h"
#include "coeffs_frEfr.hpp"

namespace DronePoseLib {
namespace ValtonenOrnhagArxiv2021 {
inline void fast_eigenvector_solver(double *eigv,
    int neig,
    const Eigen::Matrix<double, 11, 11> &AM,
    Eigen::MatrixXcd *sols);
Eigen::MatrixXcd solver_frEfr(const Eigen::VectorXd &data, const bool use_fast_solver) {
    // Compute coefficients
    Eigen::VectorXd coeffs = DronePoseLib::ValtonenOrnhagArxiv2021::coeffs_frEfr(data);

    // Setup elimination template
    static const int coeffs0_ind[] = {0, 15, 1, 0, 15, 16, 30, 45, 2, 1, 16, 17, 31, 46, 3, 15, 0, 18, 30, 45, 5, 3,
        18, 16, 1, 20, 31, 33, 48, 46, 8, 18, 3, 23, 33, 48, 6, 21, 36, 51, 4, 2, 17, 19, 32, 47, 6, 4, 19, 21, 34,
        49, 9, 24, 21, 6, 36, 39, 54, 51};
    static const int coeffs1_ind[] = {7, 5, 20, 17, 2, 22, 32, 35, 50, 47, 9, 7, 22, 19, 4, 24, 34, 37, 52, 49, 10, 8,
        23, 20, 5, 25, 35, 38, 53, 50, 11, 10, 25, 22, 7, 26, 37, 40, 55, 52, 12, 23, 8, 27, 38, 53, 11, 26, 24, 9,
        39, 41, 56, 54, 13, 12, 27, 25, 10, 28, 40, 42, 57, 55, 13, 28, 26, 11, 41, 43, 58, 56, 14, 27, 12, 29, 42,
        57, 14, 29, 28, 13, 43, 44, 59, 58, 29, 14, 44, 59};
    static const int C0_ind[] = {0, 5, 10, 11, 12, 15, 17, 18, 20, 21, 22, 25, 27, 28, 30, 33, 34, 35, 36, 39, 40, 41,
        42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 54, 55, 56, 59, 61, 62, 67, 68, 70, 71, 72, 75, 77, 78, 80, 81, 82, 85,
        87, 88, 91, 92, 93, 94, 96, 97, 98, 99};
    static const int C1_ind[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
        24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 43, 44, 45, 46, 49, 51, 52, 53, 54, 56,
        57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 71, 72, 73, 74, 76, 77, 78, 79, 80, 83, 84, 85, 86, 89, 91,
        92, 93, 94, 96, 97, 98, 99, 103, 104, 106, 109};

    Eigen::Matrix<double, 10, 10> C0; C0.setZero();
    Eigen::Matrix<double, 10, 11> C1; C1.setZero();
    for (int i = 0; i < 60; i++) {
        C0(C0_ind[i]) = coeffs(coeffs0_ind[i]);
    }
    for (int i = 0; i < 90; i++) {
        C1(C1_ind[i]) = coeffs(coeffs1_ind[i]);
    }

    Eigen::Matrix<double, 10, 11> C12 = C0.partialPivLu().solve(C1);

    // Setup action matrix
    Eigen::Matrix<double, 14, 11> RR;
    RR << -C12.bottomRows(3), Eigen::Matrix<double, 11, 11>::Identity(11, 11);

    static const int AM_ind[] = {0, 1, 3, 4, 5, 2, 6, 8, 9, 10, 12};
    Eigen::Matrix<double, 11, 11> AM;
    for (int i = 0; i < 11; i++) {
        AM.row(i) = RR.row(AM_ind[i]);
    }

    // Solve eigenvalue problem
    Eigen::MatrixXcd sols(2, 11);
    sols.setZero();
    if (use_fast_solver) {
        double p[1+11];
        Eigen::Matrix<double, 11, 11> AMp = AM;
        charpoly_danilevsky_piv(AMp, p);
        double roots[11];
        int nroots;
        find_real_roots_sturm(p, 11, roots, &nroots, 8, 0);
        fast_eigenvector_solver(roots, nroots, AM, &sols);

        // Remove unneccessary columns
        sols.conservativeResize(2, nroots);
    } else {
        Eigen::EigenSolver<Eigen::Matrix<double, 11, 11> > es(AM);
        Eigen::ArrayXcd D = es.eigenvalues();
        Eigen::ArrayXXcd V = es.eigenvectors();
        V = (V / V.row(10).array().replicate(11, 1)).eval();

        sols.row(0) = D.transpose().array();
        sols.row(1) = V.row(8).array() / (sols.row(0).array());
    }
    return sols;
}

inline void fast_eigenvector_solver(
    double *eigv,
    int neig,
    const Eigen::Matrix<double, 11, 11> &AM,
    Eigen::MatrixXcd *sols) {
    static const int ind[] = {0, 1, 5};
    // Truncated action matrix containing non-trivial rows
    Eigen::Matrix<double, 3, 11> AMs;
    double zi[5];

    for (int i = 0; i < 3; i++) {
        AMs.row(i) = AM.row(ind[i]);
    }
    for (int i = 0; i < neig; i++) {
        zi[0] = eigv[i];
        for (int j = 1; j < 5; j++) {
            zi[j] = zi[j - 1] * eigv[i];
        }
        Eigen::Matrix3d AA;
        AA.col(0) = zi[3] * AMs.col(0) + zi[2] * AMs.col(2) + zi[1] * AMs.col(4);
        AA.col(1) = zi[3] * AMs.col(1) + zi[2] * AMs.col(3) + zi[1] * AMs.col(6) + zi[0] * AMs.col(8);
        AA.col(2) = zi[2] * AMs.col(5) + zi[1] * AMs.col(7) + zi[0] * AMs.col(9) + AMs.col(10);
        AA(0, 0) = AA(0, 0) - zi[4];
        AA(1, 1) = AA(1, 1) - zi[4];
        AA(2, 2) = AA(2, 2) - zi[3];

        Eigen::Vector2d s = AA.leftCols(2).colPivHouseholderQr().solve(-AA.col(2));
        (*sols)(0, i) = zi[0];
        (*sols)(1, i) = s(1);
    }
}
}  // namespace ValtonenOrnhagArxiv2021
}  // namespace DronePoseLib
