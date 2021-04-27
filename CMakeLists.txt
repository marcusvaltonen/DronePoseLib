cmake_minimum_required(VERSION 3.2)
project(DronePoseLib)

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

include(CTest)

SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
set(EIGEN_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

include_directories(
    includes/DronePoseLib
    src/helpers
    src/solvers/valtonenornhag_arxiv_2021/fEf
    src/solvers/valtonenornhag_arxiv_2021/frEfr
    src/solvers/valtonenornhag_arxiv_2021/rEr
)

set(SRCS
    src/helpers/normalize2dpts.cpp
    src/helpers/quartic.cpp
    src/helpers/radial.cpp
    src/solvers/valtonenornhag_arxiv_2021/fEf/get_fEf.cpp
    src/solvers/valtonenornhag_arxiv_2021/frEfr/coeffs_frEfr.cpp
    src/solvers/valtonenornhag_arxiv_2021/frEfr/get_frEfr.cpp
    src/solvers/valtonenornhag_arxiv_2021/frEfr/solver_frEfr.cpp
    src/solvers/valtonenornhag_arxiv_2021/rEr/get_rEr.cpp
)

set(EXAMPLE_SRCS example.cpp)

add_executable(example ${SRCS} ${EXAMPLE_SRCS})
target_link_libraries(example PRIVATE Eigen3::Eigen)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(TEST_SRCS
        tests/test_main.cpp
        tests/test_valtonenornhag_arxiv_2021.cpp
    )
    add_executable(tests ${SRCS} ${TEST_SRCS})
    target_include_directories(tests PRIVATE tests)
    target_link_libraries(tests PRIVATE Eigen3::Eigen)
    if(CMAKE_COMPILER_IS_GNUCXX)
        LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")
        include(CodeCoverage)
        APPEND_COVERAGE_COMPILER_FLAGS()
    endif()
    add_test(NAME tests COMMAND tests)
    enable_testing()
endif()
