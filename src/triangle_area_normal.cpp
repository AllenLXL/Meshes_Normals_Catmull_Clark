#include "triangle_area_normal.h"
#include <Eigen/Geometry>

using namespace Eigen;
Eigen::RowVector3d triangle_area_normal(
  const Eigen::RowVector3d & a, 
  const Eigen::RowVector3d & b, 
  const Eigen::RowVector3d & c)
{
  // magnitude of cross product is area of parallelogram
  RowVector3d res = (b-a).cross(c-a).normalized();
  return ((b-a).cross(c-a).norm()/2.0) * res;
}
