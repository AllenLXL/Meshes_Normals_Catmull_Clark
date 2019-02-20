#include "sphere.h"
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace Eigen;
using namespace std;

//https://stackoverflow.com/questions/7840429/calculate-the-xyz-point-of-a-sphere-given-a-uv-coordinate-of-its-texture

void sphere(
  const int num_faces_u, // vertical line num longitudinal  xx
  const int num_faces_v, // horizontal line num latitudinal yy
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & UV,
  Eigen::MatrixXi & UF,
  Eigen::MatrixXd & NV,
  Eigen::MatrixXi & NF)
{
  V.resize((num_faces_u+1)*(num_faces_v+1), 3);  // #V by 3
  F.resize(num_faces_u*num_faces_v, 4);          // #F by 4
  UV.resize((num_faces_u+1)*(num_faces_v+1), 2); // #V by 2 (xyz)->(uv)
  UF.resize(num_faces_u*num_faces_v, 4);         // #F by 4

  NV.resize((num_faces_u+1)*(num_faces_v+1), 3);         // #NV by 3 #NV==#V
  NF.resize(num_faces_u*num_faces_v, 4);         // #F by 4

  double xLen=1.0/num_faces_u, yLen=1.0/num_faces_v;
  double radius=1;

  int cur=0;
  for (int i=0;i<num_faces_v+1;i++){
    for (int j=0;j<num_faces_u+1;j++){
      double u = j*xLen, v=i*yLen;

      double the = (2.0*M_PI)*u;
      double phi = (M_PI)*v;

      double x = cos(the) * sin(phi) * radius;
      double y = sin(the) * sin(phi) * radius;
      double z = -cos(phi) * radius;

      V.row(cur) << x,y,z;
      UV.row(cur) << u,v;
      NV.row(cur) << x,y,z;

      cur++;
    }
  }

  cur=0;
  for (int i=0;i<num_faces_v;i++){
    for (int j=0;j<num_faces_u;j++){

      int topLeft = i*(num_faces_u+1)+j;
      int topRight = topLeft+1;
      int botRight = (i+1)*(num_faces_u+1)+(j+1);
      int botLeft = botRight-1;

      F.row(cur) << topLeft, topRight, botRight, botLeft;
      UF.row(cur) << topLeft, topRight, botRight, botLeft;
      NF.row(cur) << topLeft, topRight, botRight, botLeft;

      cur++;
    }
  }



}
