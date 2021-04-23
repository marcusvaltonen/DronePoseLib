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

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "get_valtonenornhag_arxiv_2021.hpp"
#include "normalize2dpts.hpp"
#include "radial.hpp"
#include "relpose.hpp"
#include "solver_frEfr.hpp"

namespace DronePoseLib {
namespace ValtonenOrnhagArxiv2021 {
    inline Eigen::Vector3d extract_translation(
        const double f,
        const double r,
        const Eigen::Matrix3d &R1,
        const Eigen::Matrix3d &R2,
        const Eigen::Matrix2d &x1,
        const Eigen::Matrix2d &x2);

    std::vector<RelPose> get_frEfr(
        const Eigen::MatrixXd &p1,
        const Eigen::MatrixXd &p2,
        const Eigen::Matrix3d &R1,
        const Eigen::Matrix3d &R2,
        const bool use_fast_solver
    ) {
        // This is a 4-point method
        const int nbr_pts = 4;

        // We expect inhomogenous input data, i.e. p1 and p2 are 2x3 matrices
        assert(p1.rows() == 2);
        assert(p2.rows() == 2);
        assert(p1.cols() == nbr_pts);
        assert(p2.cols() == nbr_pts);

        // Compute normalization matrix
        double scale1 = normalize2dpts(p1);
        double scale2 = normalize2dpts(p2);
        double scale = std::max(scale1, scale2);
        Eigen::Vector3d s;
        s << scale, scale, 1.0;
        Eigen::DiagonalMatrix<double, 3> S = s.asDiagonal();

        // Normalize data
        Eigen::Matrix<double, 3, nbr_pts> x1;
        Eigen::Matrix<double, 3, nbr_pts> x2;
        x1 = p1.colwise().homogeneous();
        x2 = p2.colwise().homogeneous();
        x1 = S * x1;
        x2 = S * x2;

        Eigen::Matrix<double, 2, nbr_pts> x1t;
        Eigen::Matrix<double, 2, nbr_pts> x2t;
        x1t << x1.colwise().hnormalized();
        x2t << x2.colwise().hnormalized();

        // Compute relative rotation
        Eigen::Matrix3d R = R2 * R1.transpose();

        // Wrap input data to expected format
        Eigen::VectorXd input(25);
        input << Eigen::Map<Eigen::VectorXd>(x1t.data(), 8),
                 Eigen::Map<Eigen::VectorXd>(x2t.data(), 8),
                 Eigen::Map<Eigen::VectorXd>(R.data(), 9);

        // Extract solution
        Eigen::MatrixXcd sols = DronePoseLib::ValtonenOrnhagArxiv2021::solver_frEfr(input, use_fast_solver);

        // Pre-processing: Remove complex-valued solutions
        double thresh = 1e-12;
        Eigen::ArrayXd real_sols(11);
        real_sols = sols.imag().cwiseAbs().colwise().sum();
        int nbr_real_sols = (real_sols <= thresh).count();

        // Construct putative output
        // Eigen::Vector3d t;
        std::vector<RelPose> output;
        RelPose relpose;
        double f, r;
        Eigen::Vector3d kinv;
        Eigen::DiagonalMatrix<double, 3> Kinv;
        Eigen::Matrix3d skew_t;

        // Loop over real solutions
        for (int i=0; i < real_sols.size(); i++) {
            if (real_sols(i) < thresh) {
                f = sols(0, i).real();
                r = sols(1, i).real();

                // Extract translation
                relpose.t = extract_translation(f, r, R1, R2, x1t.leftCols<2>(), x2t.leftCols<2>());
                relpose.f = f / scale;
                relpose.r = r * std::pow(scale, 2);

                // Compute fundamental matrix
                kinv << 1.0 / relpose.f, 1.0 / relpose.f, 1.0;
                Kinv = kinv.asDiagonal();
                skew_t << 0, -relpose.t(2), relpose.t(1),
                          relpose.t(2), 0, -relpose.t(0),
                         -relpose.t(1), relpose.t(0), 0;
                relpose.F = Kinv * skew_t * R * Kinv;

                // Add
                output.push_back(relpose);
            }
        }

        return output;
    }

    inline Eigen::Vector3d extract_translation(
        const double f,
        const double r,
        const Eigen::Matrix3d &R1,
        const Eigen::Matrix3d &R2,
        const Eigen::Matrix2d &x1,
        const Eigen::Matrix2d &x2
    ) {
        Eigen::Vector3d fmat;
        fmat << 1.0 / f, 1.0 / f, 1.0;
        Eigen::DiagonalMatrix<double, 3> Kinv = fmat.asDiagonal();

        // Transform points
        Eigen::Matrix<double, 3, 2> y1, y2;
        y1 = R1.transpose() * Kinv * DronePoseLib::radialundistort(x1, r).colwise().homogeneous();
        y2 = R2.transpose() * Kinv * DronePoseLib::radialundistort(x2, r).colwise().homogeneous();

        // Extract translation
        Eigen::Vector3d t;
        t = R2 * y1.col(0).cross(y2.col(0)).cross(y1.col(1).cross(y2.col(1)));
        return t;
    }
}  // namespace ValtonenOrnhagArxiv2021
}  // namespace DronePoseLib
