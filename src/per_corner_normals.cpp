#include "per_corner_normals.h"
#include "triangle_area_normal.h"
// Hint:
#include "vertex_triangle_adjacency.h"
#include <iostream>

#include <unordered_map>

using namespace std;
using namespace Eigen;

void per_corner_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double corner_threshold,
  Eigen::MatrixXd & N)
{
  N = Eigen::MatrixXd::Zero(F.rows()*3,3);
  vector<vector<int> > VF;
  vertex_triangle_adjacency(F, static_cast<const int>(V.rows()), VF);

  int cur=0;
  for (int i=0; i<F.rows(); i++) {
    RowVector3d curFcNml = triangle_area_normal(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2))).normalized();

    for (int j=0;j<F.cols();j++){
      int vtxIdx = F(i,j);
      vector<int> adjVtx = VF[vtxIdx];

      RowVector3d res(0.0,0.0,0.0);
      for (int k=0; k<adjVtx.size();k++){
        RowVector3d triPt0 = V.row(F(adjVtx[k],0));
        RowVector3d triPt1 = V.row(F(adjVtx[k],1));
        RowVector3d triPt2 = V.row(F(adjVtx[k],2));

        RowVector3d adjNml = triangle_area_normal(triPt0, triPt1,triPt2).normalized();

        if (curFcNml.dot(adjNml) > cos(corner_threshold)){
          res += triangle_area_normal(triPt0, triPt1,triPt2);
        }
      }
      N.row(cur) << res.normalized();
      cur++;
    }

  }

}
