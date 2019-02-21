#include "catmull_clark.h"
#include <unordered_map>
#include <utility>
#include <functional>

#include <vector>

#include <string>
#include <iostream>

using namespace Eigen;
using namespace std;

void catmull_clark(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_iters,
  Eigen::MatrixXd & SV,
  Eigen::MatrixXi & SF) {
  //////////////////////////////////////////////////////////////////////////

  if (num_iters==0){
    return;
  }

  int numF = (int) F.rows();
  int numV = (int) V.rows();

  SV=MatrixXd::Zero(numV*3+numF,3);
  SF=MatrixXi::Zero(numF*4,4);

  // facePoint #F * 3
  MatrixXd facePoint = MatrixXd::Zero(numF,3);

  for(int i=0;i<numF;i++){
    // op for original point
    int op0 = F(i,0), op1 = F(i,1), op2 = F(i,2), op3 = F(i,3);

    facePoint(i,0) = (V(op0,0), V(op1,0), V(op2,0), V(op3,0))/4.0;
    facePoint(i,1) = (V(op0,1), V(op1,1), V(op2,1), V(op3,1))/4.0;
    facePoint(i,2) = (V(op0,2), V(op1,2), V(op2,2), V(op3,2))/4.0;
  }

  /*
  First we collect all the information that we need to run the algorithm.
  Each point must know its adjacent edges and faces.
  Each face must know its edges and points.
  Each edge must know its adjacent faces and points.
  We collect all this information in this loop.
   */
  // ij entry stores face idx into F
  MatrixXi edgeMat = MatrixXi::Constant(numV,numV,-1);
  for (int i=0;i<numF;i++){
    int pt0=F(i,0);
    int pt1=F(i,1);
    int pt2=F(i,2);
    int pt3=F(i,3);
    edgeMat(pt0,pt1) = i;
//    edgeMat(pt1,pt0) = i;

    edgeMat(pt1,pt2) = i;
//    edgeMat(pt2,pt1) = i;

    edgeMat(pt2,pt3) = i;
//    edgeMat(pt3,pt2) = i;

    edgeMat(pt3,pt0) = i;
//    edgeMat(pt0,pt3) = i;
  }

  MatrixXd edgePtMat = MatrixXd::Zero(numF*2,3);

  // ij entry store idx x => edge ij has edge point at edgePtMat.row(x)
//  cout << "74=======================DEBUG=======================" << endl;
  MatrixXi edgeToEdgePt = edgeMat;
  int cur=0;

  for (int i=0;i<numV;i++){
    RowVector3d avgFacePt;
    for (int j=0;j<numV;j++){
      if (edgeToEdgePt(i,j)!=-1 && i<j){
        int faceId1, faceId2;

        faceId1=edgeToEdgePt(i,j);
        faceId2=edgeToEdgePt(j,i);

        avgFacePt = facePoint.row(faceId1) + facePoint.row(faceId2) + V.row(i) + V.row(j);
        avgFacePt /= 4.0;

        edgeToEdgePt(i, j) = cur;
        edgeToEdgePt(j, i) = cur;

        edgePtMat.row(cur) << avgFacePt;
        cur++;
      }
//      cout <<"i j is " << i << " " << j << endl;
    }
  }


  // fill SF matrix
  cur=0;
  for (int i=0;i<numF;i++){
    int facePtId = i;
    int oriPtId0 = F(i,0), oriPtId1 = F(i,1), oriPtId2 = F(i,2), oriPtId3 = F(i,3);
    // TODO counter clock wise attempt first
    // 1st sub face
    SF(cur,0) = oriPtId0;
    SF(cur,1) = edgeToEdgePt(oriPtId0,oriPtId1)+numV;
    SF(cur,2) = i+(numV*3);
    SF(cur,3) = edgeToEdgePt(oriPtId3,oriPtId0)+numV;
    cur++;

    // 2nd sub face
    SF(cur,0) = oriPtId1;
    SF(cur,1) = edgeToEdgePt(oriPtId1,oriPtId2)+numV;
    SF(cur,2) = i+(numV*3);
    SF(cur,3) = edgeToEdgePt(oriPtId0,oriPtId1)+numV;
    cur++;

    // 3rd sub face
    SF(cur,0) = oriPtId2;
    SF(cur,1) = edgeToEdgePt(oriPtId2,oriPtId3)+numV;
    SF(cur,2) = i+(numV*3);
    SF(cur,3) = edgeToEdgePt(oriPtId1,oriPtId2)+numV;
    cur++;

    // 4th sub face
    SF(cur,0) = oriPtId3;
    SF(cur,1) = edgeToEdgePt(oriPtId3,oriPtId0)+numV;
    SF(cur,2) = i+(numV*3);
    SF(cur,3) = edgeToEdgePt(oriPtId2,oriPtId3)+numV;
    cur++;
  }


  MatrixXd newV = MatrixXd::Zero(numV,3);
  for (int i=0;i<numV;i++){
    vector<int> nbrPt;
    for (int j=0;j<numV;j++){
      if (edgeMat(i,j)!=-1){
        nbrPt.push_back(j);
      }
//      cout <<"i j is " << i << " " << j << endl;
    }
    int N = static_cast<int>(nbrPt.size());
    RowVector3d FSum(0.0,0.0,0.0), RSum(0.0,0.0,0.0);
    for (int k=0;k<N;k++){
      RowVector3d temp = facePoint.row(edgeMat(i,nbrPt[k]));
      FSum += temp;
      temp = edgePtMat.row(edgeToEdgePt(i,nbrPt[k]));
      RSum += temp;
    }
//    cout << "152=======================DEBUG=======================" << endl;
    RowVector3d FF = FSum/N;
    RowVector3d RR = RSum/N;
    RowVector3d PP = V.row(i);

    RowVector3d newP = (FF+2.0*RR+(N-1)*PP)/N;
    newV.row(i) << newP;

  }
//  cout << "=======================DEBUG=======================" << endl;


  // new SV
  SV<<newV,edgePtMat,facePoint;
  catmull_clark(MatrixXd(SV), MatrixXi(SF), num_iters-1, SV, SF);
}

