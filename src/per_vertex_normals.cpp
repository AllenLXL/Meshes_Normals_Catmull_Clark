#include "per_vertex_normals.h"
#include "triangle_area_normal.h"

using namespace Eigen;
using namespace std;

void per_vertex_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
  N = Eigen::MatrixXd::Zero(V.rows(),3);

  MatrixXi edgeMat = MatrixXi::Constant(V.rows(),V.rows(),-1);
  for (int i=0;i<F.rows();i++){
    int pt0=F(i,0);
    int pt1=F(i,1);
    int pt2=F(i,2);
    int pt3=F(i,3);

    edgeMat(pt0,pt1) = i;
    edgeMat(pt1,pt2) = i;
    edgeMat(pt2,pt3) = i;
    edgeMat(pt3,pt0) = i;
  }

  for (int i=0; i<F.rows(); i++) {
    double sumArea = 0.0;
//    RowVector3d triNml = triangle_area_normal(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));

    for (int j=0; j<V.rows();j++) {

      if (edgeMat(i,j)!=-1){
        int faceIdx = edgeMat(i,j);
        sumArea += triangle_area_normal(V.row(F(faceIdx,0)),V.row(F(faceIdx,1)),V.row(F(faceIdx,2))).norm();
        N.row(i) += triangle_area_normal(V.row(F(faceIdx, 0)), V.row(F(faceIdx, 1)), V.row(F(faceIdx, 2)));
      }
    }
    N.row(i) = (N.row(i)/sumArea).normalized();
  }

}
