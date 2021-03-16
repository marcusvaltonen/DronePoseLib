#ifndef RELPOSEH
#define RELPOSEH

#include <Eigen/Dense>

namespace DronePoseLib {
struct RelPose {
    Eigen::Matrix3d F;
    Eigen::Vector3d t;
    double f;
    double r;
};
}  // namespace DronePoseLib

#endif /* POSEDATAH */
