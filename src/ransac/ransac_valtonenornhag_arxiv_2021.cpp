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

#include "ransac_valtonenornhag_arxiv_2021.hpp"
#include <Eigen/Dense>
#include <vector>
#include "pose_estimator.hpp"
#include "get_valtonenornhag_arxiv_2021.hpp"
#include "relpose.hpp"
// #include "distortion.hpp"
// #include "refinement.hpp"



// DEBUG
#include <iostream>
#include "scene_and_pose_generation.hpp"


// TODO: Change name to reflect frEfr
int DronePoseLib::ValtonenOrnhagArxiv2021::Solver::solve(
    DronePoseLib::Points2D &x1,
    DronePoseLib::Points2D &x2,
    const Eigen::Matrix3d &R1,
    const Eigen::Matrix3d &R2,
    std::vector<DronePoseLib::Camera>* poses) const {

    std::cout << "calling solve()" << std::endl;
    std::cout << "R1 = " << R1 << std::endl;
    std::cout << "R2 = " << R2 << std::endl;
    std::cout << "x1 = " << x1 << std::endl;
    std::cout << "x2 = " << x2 << std::endl;

    // TODO: Consider replaceing RelPose with Camera, and just call a computerFundamentalMatrix function.
    std::vector<DronePoseLib::RelPose> relpose = DronePoseLib::ValtonenOrnhagArxiv2021::get_frEfr(x1, x2, R1, R2, false);

    std::cout << "nbr poses = " << relpose.size() << std::endl;

    // TODO: Consider saving the relative pose early on and send it in instead
    for (int i=0; i < relpose.size(); i++) {
        DronePoseLib::Camera p;
        p.R = R2 * R1.transpose();
        p.t = relpose[i].t;
        p.focal = relpose[i].f;
        p.dist_params.push_back(relpose[i].r);

        std::cout << "nbr  = " << i << std::endl;
        debug_print_pose(p);
        poses->push_back(p);
    }

    return poses->size();
}

template class DronePoseLib::PoseEstimator<DronePoseLib::ValtonenOrnhagArxiv2021::Solver>;
