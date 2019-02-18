#include "write_obj.h"
#include <fstream>
#include <cassert>
#include <iostream>

#include <Eigen/Core>

using namespace std;
bool write_obj(
  const std::string & filename,
  const Eigen::MatrixXd & V,     // Double
  const Eigen::MatrixXi & F,     // Int
  const Eigen::MatrixXd & UV,    // Double
  const Eigen::MatrixXi & UF,    // Int
  const Eigen::MatrixXd & NV,    // Double
  const Eigen::MatrixXi & NF)    // Int
{
  assert((F.size() == 0 || F.cols() == 3 || F.cols() == 4) && "F must have 3 or 4 columns");
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
  ofstream objFile;
  objFile.open(filename);
  if (objFile.is_open()){ // successfully opened
    // List of geometric vertices, with (x, y, z)
    for (int i=0;i<V.rows();i++){
        objFile << "v " << V(i,0) <<' ' << V(i,1) <<' ' << V(i,2) << endl;
    }
    // List of texture coordinates,
    for (int i=0;i<UV.rows();i++){
        objFile << "vt " << UV(i,0) <<' ' << UV(i,1) << endl;
    }
    // List of vertex normals
    for (int i=0;i<NV.rows();i++){
      objFile << "vn " << NV(i,0) <<' ' << NV(i,1) <<' ' << NV(i,2) << endl;
    }
    // no vp
    // Vertex indices/Vertex texture coordinate indices/Vertex normal indices
    // F/UF/NF (NOTE same num of rows)
    for (int i=0;i<F.rows();i++){
      objFile<<"f ";
      if (F.cols() == 3){
        objFile<<F(i,0)+1<<'/'<<UF(i,0)+1<<'/'<<NF(i,0)+1<<' ';
        objFile<<F(i,1)+1<<'/'<<UF(i,1)+1<<'/'<<NF(i,1)+1<<' ';
        objFile<<F(i,2)+1<<'/'<<UF(i,2)+1<<'/'<<NF(i,2)+1<<endl;
      } else{
        objFile<<F(i,0)+1<<'/'<<UF(i,0)+1<<'/'<<NF(i,0)+1<<' ';
        objFile<<F(i,1)+1<<'/'<<UF(i,1)+1<<'/'<<NF(i,1)+1<<' ';
        objFile<<F(i,2)+1<<'/'<<UF(i,2)+1<<'/'<<NF(i,2)+1<<' ';
        objFile<<F(i,3)+1<<'/'<<UF(i,3)+1<<'/'<<NF(i,3)+1<<endl;
      }
    }
    objFile.close();
    return true;
  } else{  // error opening file
    return false;
  }
  ////////////////////////////////////////////////////////////////////////////
}
