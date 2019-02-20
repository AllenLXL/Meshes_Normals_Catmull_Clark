#include "catmull_clark.h"
#include <unordered_map>
#include <utility>
#include <functional>
#include <vector>

//#include <string>
//#include <iostream>

using namespace Eigen;
using namespace std;

void catmull_clark(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_iters,
  Eigen::MatrixXd & SV,
  Eigen::MatrixXi & SF)
{
  ////////////////////////////////////////////////////////////////////////////
//   Replace with your code here:
  if (num_iters==0){
    return;
  }

  int numF = (int) F.rows();
  int numV = (int) V.rows();

  SV.resize(numV*3+numF,3);
  SF.resize(numF*4,4);

  // facePoint #F * 3
  MatrixXd facePoint = MatrixXd::Zero(numF,3);

  for(int i=0;i<numF;i++){
    // op for original point
    int op0 = F(i,0), op1 = F(i,1), op2 = F(i,2), op3 = F(i,3);

    facePoint(i,0) = (V(op0,0), V(op1,0), V(op2,0), V(op3,0))/4;
    facePoint(i,1) = (V(op0,1), V(op1,1), V(op2,1), V(op3,1))/4;
    facePoint(i,2) = (V(op0,2), V(op1,2), V(op2,2), V(op3,2))/4;
  }
  /*
  First we collect all the information that we need to run the algorithm.
  Each point must know its adjacent edges and faces.
  Each face must know its edges and points.
  Each edge must know its adjacent faces and points.
  We collect all this information in this loop.
   */
  // ij entry stores face idx into F
  MatrixXi edgeMat = MatrixXi::Zero(numV,numV);
  for (int i=0;i<numF;i++){
    int pt0=F(i,0);
    int pt1=F(i,1);
    int pt2=F(i,2);
    int pt3=F(i,3);
    edgeMat(pt0,pt1) = i;
    edgeMat(pt1,pt2) = i;
    edgeMat(pt2,pt3) = i;
    edgeMat(pt3,pt0) = i;
  }

  MatrixXi adjVtxMat = MatrixXi::Zero(numV,4);
  int cur;
  for (int i=0;i<numV;i++){
    cur=0;
    for (int j=0;j<numV;j++){
      if (edgeMat(i,j)!=0){
        adjVtxMat(i,cur)=j;
        cur++;
      }
    }
  }

  MatrixXi vtxToFc = MatrixXi::Zero(numV,4);
  cur=0;
  for (int i=0;i<numV;i++){
    cur=0;
    for (int j=0;j<numV;j++){
      if (edgeMat(i,j)!=0){
        adjVtxMat(i,cur)=edgeMat(i,j);
        cur++;
      }
    }
  }


  MatrixXd edgePtMat = MatrixXd::Zero(numV*2,3);

  // below matrix only half filled
  MatrixXi edgeToEdgePt = MatrixXd::Zero(numV,numV);
  cur=0;

  for (int i=0;i<numV;i++){
    RowVector3d avgFacePt;
    for (int j=0;j<numV;j++){

      if (edgeMat(i,j)!=0 && i<j){
        int faceId1, faceId2;

        edgeToEdgePt(i, j) = cur;
        edgeToEdgePt(j, i) = cur;

        faceId1=adjVtxMat(i,j); faceId2=adjVtxMat(j,i);
        avgFacePt = facePoint.row(faceId1) + facePoint.row(faceId2) + V.row(i) + V.row(j);
        avgFacePt /= 4;

        edgePtMat.row(cur) << avgFacePt;
        cur++;
      }
    }
  }
  // fill SF matrix
  for (int i=0;i<numF;i++){
    int facePtId = i;
    int oriPtId0 = F(i,0), oriPtId1 = F(i,1), oriPtId2 = F(i,2), oriPtId3 = F(i,3);
    for (int j=0;j<4;j++){
      // TODO counter clock wise attempt first

    }
  }

  MatrixXd newV = MatrixXd::Zero(numV,3);
  for (int i=0;i<numV;i++){
    int nbrFcId0, nbrFcId1, nbrFcId2, nbrFcId3;
    nbrFcId0 = vtxToFc(i,0);
    nbrFcId1 = vtxToFc(i,1);
    nbrFcId2 = vtxToFc(i,2);
    nbrFcId3 = vtxToFc(i,3);

    RowVector3d facePt0, facePt1, facePt2, facePt3;
    facePt0 = facePoint.row(nbrFcId0);
    facePt1 = facePoint.row(nbrFcId1);
    facePt2 = facePoint.row(nbrFcId2);
    facePt3 = facePoint.row(nbrFcId3);
    RowVector3d F = (facePt0+facePt1+facePt2+facePt3)/4.0;

    int nbrPt0, nbrPt1, nbrPt2, nbrPt3;
    nbrPt0 = adjVtxMat(i,0);
    nbrPt1 = adjVtxMat(i,1);
    nbrPt2 = adjVtxMat(i,2);
    nbrPt3 = adjVtxMat(i,3);

    RowVector3d P = V.row(i);
    RowVector3d edgeMid0, edgeMid1, edgeMid2, edgeMid3;
    edgeMid0 = (P+V.row(nbrPt0))/2.0;
    edgeMid1 = (P+V.row(nbrPt1))/2.0;
    edgeMid2 = (P+V.row(nbrPt2))/2.0;
    edgeMid3 = (P+V.row(nbrPt3))/2.0;

    RowVector3d R = (edgeMid0+edgeMid1+edgeMid2+edgeMid3)/4.0;

    // n=4 since quad mesh
    RowVector3d newP = (F+2*R+P)/4;

    newV.row(i) << newP;
  }

  // trying to get all edgePt

  // new SV
  SV<<V,edgePtMat,facePoint;




  ////////////////////////////////////////////////////////////////////////////
}

