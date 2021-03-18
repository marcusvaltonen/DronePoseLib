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
#include "get_valtonenornhag_arxiv_2021.hpp"
#include "relpose.hpp"

#ifdef MATLAB_MEX_FILE
#include "mex.h"  // NOLINT [build/include_subdir]
#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names))
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 4 || nrhs > 5) {
        mexErrMsgIdAndTxt("get_valtonenornhag_arxiv_2021_frEfr:nrhs", "Four or five inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("get_valtonenornhag_arxiv_2021_frEfr:nlhs", "One output required.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
        !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
        !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("get_valtonenornhag_arxiv_2021_frEfr:notDouble", "Input data must be type double.");
    }
    if (mxGetNumberOfElements(prhs[0]) != 8) {
        mexErrMsgIdAndTxt("get_valtonenornhag_arxiv_2021_frEfr:incorrectSize1",
                          "Incorrect input size of first argument");
    }
    if (mxGetNumberOfElements(prhs[1]) != 8) {
        mexErrMsgIdAndTxt("get_valtonenornhag_arxiv_2021_frEfr:incorrectSize2",
                          "Incorrect input size of second argument");
    }
    if (mxGetNumberOfElements(prhs[2]) != 9) {
        mexErrMsgIdAndTxt("get_valtonenornhag_arxiv_2021_frEfr:incorrectSize3",
                          "Incorrect input size of third argument");
    }
    if (mxGetNumberOfElements(prhs[3]) != 9) {
        mexErrMsgIdAndTxt("get_valtonenornhag_arxiv_2021_frEfr:incorrectSize4",
                          "Incorrect input size of fourth argument");
    }
    if (nrhs == 5 && (!mxIsLogical(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1)) {
        mexErrMsgIdAndTxt("get_valtonenornhag_arxiv_2021_frEfr:notLogical", "Last input must be type logical.");
    }
    bool use_fast_solver = false;
    if (nrhs == 5) {
        bool *var = mxGetLogicals(prhs[4]);
        use_fast_solver = var[0];
    }

    // Convert to expected input
    // TODO(marcusvaltonen): Cast directly to MatrixXd?
    Eigen::VectorXd x1_tmp = Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[0]), 8);
    Eigen::VectorXd x2_tmp = Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[1]), 8);
    Eigen::MatrixXd x1 = Eigen::Map<Eigen::MatrixXd>(x1_tmp.data(), 2, 4);
    Eigen::MatrixXd x2 = Eigen::Map<Eigen::MatrixXd>(x2_tmp.data(), 2, 4);
    Eigen::VectorXd R1_tmp = Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[2]), 9);
    Eigen::VectorXd R2_tmp = Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[3]), 9);
    Eigen::MatrixXd R1 = Eigen::Map<Eigen::MatrixXd>(R1_tmp.data(), 3, 3);
    Eigen::MatrixXd R2 = Eigen::Map<Eigen::MatrixXd>(R2_tmp.data(), 3, 3);

    // Compute output
    std::vector<DronePoseLib::RelPose> posedata =
        DronePoseLib::ValtonenOrnhagArxiv2021::get_frEfr(x1, x2, R1, R2, use_fast_solver);

    // Wrap it up to Matlab compatible output
    std::size_t NUMBER_OF_STRUCTS = posedata.size();
    const char *field_names[] = {"t", "f", "r", "F"};
    mwSize dims[2] = {1, NUMBER_OF_STRUCTS };
    int t_field, f_field, r_field, F_field;
    mwIndex i;

    plhs[0] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names);

    t_field = mxGetFieldNumber(plhs[0], "t");
    f_field = mxGetFieldNumber(plhs[0], "f");
    r_field = mxGetFieldNumber(plhs[0], "r");
    F_field = mxGetFieldNumber(plhs[0], "F");

    double* zr;
    for (i = 0; i < NUMBER_OF_STRUCTS; i++) {
        mxArray *field_value;

        // Create t
        field_value = mxCreateDoubleMatrix(3, 1, mxREAL);
        zr = mxGetPr(field_value);
        for (Eigen::Index j = 0; j < 3; j++) {
            zr[j] = posedata[i].t(j);
        }
        mxSetFieldByNumber(plhs[0], i, t_field, field_value);

        // Create f
        field_value = mxCreateDoubleMatrix(1, 1, mxREAL);
        zr = mxGetPr(field_value);
        zr[0] = posedata[i].f;
        mxSetFieldByNumber(plhs[0], i, f_field, field_value);

        // Create r
        field_value = mxCreateDoubleMatrix(1, 1, mxREAL);
        zr = mxGetPr(field_value);
        zr[0] = posedata[i].r;
        mxSetFieldByNumber(plhs[0], i, r_field, field_value);

        // Create F
        field_value = mxCreateDoubleMatrix(3, 3, mxREAL);
        zr = mxGetPr(field_value);
        for (Eigen::Index j = 0; j < 9; j++) {
            zr[j] = posedata[i].F(j);
        }
        mxSetFieldByNumber(plhs[0], i, F_field, field_value);
    }
}
#endif  // MATLAB_MEX_FILE
