# DronePoseLib

[![Build Status](https://travis-ci.com/marcusvaltonen/DronePoseLib.svg?branch=main)](https://travis-ci.com/marcusvaltonen/DronePoseLib)
![GitHub](https://img.shields.io/github/license/marcusvaltonen/DronePoseLib)

Library for Visual-Inertial Odometry (VIO) using minimal solvers utilizing a common reference
direction (obtained from IMU data). The code is related to the ArXiV paper [[link](https://arxiv.org/abs/2103.08286)]:

```
@misc{valtonenornhag-etal-2021-arxiv,
      title={Trust Your IMU: Consequences of Ignoring the IMU Drift},
      author={Marcus {Valtonen~Örnhag} and Patrik Persson and Mårten Wadenbäck and Kalle Åström and Anders Heyden},
      year={2021},
      eprint={2103.08286},
      archivePrefix={arXiv},
      primaryClass={cs.CV}
}
```

## Dependencies
The implementation uses Eigen 3 (older versions are not compatible), which is
a C++ template library for linear algebra: matrices, vectors,
numerical solvers, and related algorithms.

Installation for Ubuntu/Debian:

```bash
    $ apt-get install libeigen3-dev
```
The source code has been compiled and tested on Ubuntu 18.04 (Bionic Beaver) with clang++-7 to clang++-9.
Furthermore, it is tested on OSX with Xcode 12.

### Compiling with g++
It is possible to compile it using g++ (tested on version 7.5.0) but requires more
memory (tested succesfully with 16 GB internal memory).
We are working on reducing the memory footprint.

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
