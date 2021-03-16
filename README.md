# DronePoseLib

[![Build Status](https://travis-ci.com/marcusvaltonen/DronePoseLib.svg?branch=main)](https://travis-ci.com/marcusvaltonen/DronePoseLib)
![GitHub](https://img.shields.io/github/license/marcusvaltonen/DronePoseLib)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/marcusvaltonen/DronePoseLib.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/marcusvaltonen/DronePoseLib/context:cpp)
[![codecov](https://codecov.io/gh/marcusvaltonen/DronePoseLib/branch/main/graph/badge.svg)](https://codecov.io/gh/marcusvaltonen/DronePoseLib)

Library for Visual-Inertial Odometry (VIO) using minimal solvers utilizing a common reference
direction (obtained from IMU data). The code is related to the ArXiV paper (coming soon).

## Dependencies
The implementation uses Eigen 3 (older versions are not compatible), which is
a C++ template library for linear algebra: matrices, vectors,
numerical solvers, and related algorithms.

Installation for Ubuntu/Debian:

```bash
    $ apt-get install libeigen3-dev
```
Tested on version 3.3.4.

## Using the solver in MATLAB
It is possible to MEX-compile the solver and use it in MATLAB.
In order to do so you may use the file `matlab/compile_mex.m`, which
should result in an executable.

Note that your local Eigen path may be different, e.g. `/usr/local/include/eigen3`.
Tested on version R2020a Linux (64-bit).

## Using the solver in Python
This is unfortunately not supported yet, but will be in the future.

## Additional information
These solvers were generated using the automatic generator proposed by
Larsson et al. "Efficient Solvers for Minimal Problems by Syzygy-based
Reduction" (CVPR 2017)
