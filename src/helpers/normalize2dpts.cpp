#include <math.h>
#include "normalize2dpts.hpp"

using namespace Eigen;

double normalize2dpts(MatrixXd &pts)
{
    double scale = sqrt(2) / pts.colwise().norm().mean();
    return scale;
}
