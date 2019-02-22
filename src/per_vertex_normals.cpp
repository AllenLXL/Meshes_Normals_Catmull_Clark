#include "per_vertex_normals.h"
#include "triangle_area_normal.h"
#include <unordered_map>
#include <vector>
//#include <Eigen/Dense>
//#include <Eigen/Geometry>

using namespace Eigen;
using namespace std;

void per_vertex_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
  N.resize(V.rows(),3);
//  N = MatrixXd::Zero(V.rows(),3);
  // pt to its adj face
  unordered_map<int, vector<int>> ptToAdjFc;

  for (int i=0;i<F.rows();i++){
    ptToAdjFc[F(i,0)].push_back(i);

    ptToAdjFc[F(i,1)].push_back(i);

    ptToAdjFc[F(i,2)].push_back(i);
  }


  for (int i=0;i<V.rows();i++){
    RowVector3d res(0.0,0.0,0.0);
//    RowVector3d perFcNml;
    for (int j=0;j<ptToAdjFc[i].size();j++){
      int triPt0 = F.row(ptToAdjFc[i][j])(0);
      int triPt1 = F.row(ptToAdjFc[i][j])(1);
      int triPt2 = F.row(ptToAdjFc[i][j])(2);

      res += triangle_area_normal(V.row(triPt0),V.row(triPt1),V.row(triPt2));
    }
    N.row(i) << res.normalized();
  }

}
