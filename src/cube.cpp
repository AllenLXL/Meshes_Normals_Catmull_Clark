#include "cube.h"

void cube(
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & UV,
  Eigen::MatrixXi & UF,
  Eigen::MatrixXd & NV,
  Eigen::MatrixXi & NF)
{

  V.resize(8,3);
  F.resize(6,4);
  UV.resize(14,2);
  UF.resize(6,4);
  NV.resize(6,3);
  NF.resize(6,4);

  V << 0,0,0,
      0,1,0,
      1,1,0,
      1,0,0,
      0,0,1,
      0,1,1,
      1,1,1,
      1,0,1;

  // NOTE face points has to be consecutive
  F << 0,1,2,3,
      0,3,7,4,
      0,1,5,4,
      4,5,6,7,
      2,3,7,6,
      1,2,6,5;

  UV << 0,2,
       0,1,
       1,3,
       1,2,
       1,1,
       1,0,
       2,3,
       2,2,
       2,1,
       2,0,
       3,2,
       3,1,
       4,2,
       4,1;

  UV/=4;

  UF << 4,3,7,8,
       4,8,9,5,
       4,3,0,1,
       13,12,10,11,
       7,8,11,10,
       3,7,6,2;

  NV << 0,0,-1,
       0,0,1,
       1,0,0,
       -1,0,0,
       0,1,0,
       0,-1,0;

  NF << 0,0,0,0,
       5,5,5,5,
       3,3,3,3,
       1,1,1,1,
       2,2,2,2,
       4,4,4,4;

}
