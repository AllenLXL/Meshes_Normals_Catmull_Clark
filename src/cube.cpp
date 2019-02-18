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

   F << 0,1,2,3,
        0,3,4,7,
        0,1,4,5,
        4,5,6,7,
        2,3,6,7,
        1,2,5,6;

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

   UF << 4,3,7,8,
         4,8,5,9,
         4,3,1,0,
         13,12,10,11,
         7,8,10,11,
         3,7,2,6;

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
