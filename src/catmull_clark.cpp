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

    facePoint.row(i) << (V.row(op0)+V.row(op1)+V.row(op2)+V.row(op3))/4.0;
  }
  // ij entry stores face idx into F

  // pt to its adj pts
  unordered_map<int, vector<int>> ptToAdjPt;
  // pt to its adj face
  unordered_map<int, vector<int>> ptToAdjFc;

  for (int i=0;i<numF;i++){
    ptToAdjPt[F(i,0)].push_back(F(i,1));
    ptToAdjFc[F(i,0)].push_back(i);

    ptToAdjPt[F(i,1)].push_back(F(i,2));
    ptToAdjFc[F(i,1)].push_back(i);

    ptToAdjPt[F(i,2)].push_back(F(i,3));
    ptToAdjFc[F(i,2)].push_back(i);

    ptToAdjPt[F(i,3)].push_back(F(i,0));
    ptToAdjFc[F(i,3)].push_back(i);
  }


  MatrixXd edgePtMat = MatrixXd::Zero(numF*2,3);

  // ij entry store idx x => edge ij has edge point at edgePtMat.row(x)
  MatrixXi edgeToEdgePt = MatrixXi::Constant(numV,numV,-1);

  int cur=0;
  for (int i=0;i<numV;i++){
    int ctrPtId = i;
    vector<int> surroundingPt = ptToAdjPt[ctrPtId];
    vector<int> surroundingFc = ptToAdjFc[ctrPtId];
    for (int j=0;j<surroundingPt.size();j++){
      int endPtId = surroundingPt[j];
      if (ctrPtId < endPtId){
        int faceId1;
        faceId1=surroundingFc[j];
        // want face id2
        int faceId2;
        for (int k=0;k<ptToAdjPt[endPtId].size();k++){
          if (ptToAdjPt[endPtId][k]==ctrPtId){
            faceId2 = ptToAdjFc[endPtId][k];
            break;
          }
        }
        edgeToEdgePt(ctrPtId, endPtId) = cur;
        edgeToEdgePt(endPtId, ctrPtId) = cur;

        RowVector3d avgFacePt = (facePoint.row(faceId1) + facePoint.row(faceId2) + V.row(ctrPtId) + V.row(endPtId))/4.0;
        edgePtMat.row(cur) << avgFacePt;
        cur++;
      }
    }
  }
  SV<<V,edgePtMat,facePoint;

  // fill SF matrix
  cur=0;

  for (int i=0;i<numF;i++){
    int facePtId = i;
    int oriPtId0 = F(i,0), oriPtId1 = F(i,1), oriPtId2 = F(i,2), oriPtId3 = F(i,3);
    // counter clock wise attempt first
    // 1st sub face
    RowVector4i temp0(oriPtId0, edgeToEdgePt(oriPtId0,oriPtId1)+numV, i+(numV*3), edgeToEdgePt(oriPtId3,oriPtId0)+numV);
    SF.row(cur) << temp0;
    cur++;

    // 2nd sub face
    RowVector4i temp1(oriPtId1, edgeToEdgePt(oriPtId1,oriPtId2)+numV, i+(numV*3), edgeToEdgePt(oriPtId0,oriPtId1)+numV);
    SF.row(cur) << temp1;
    cur++;

    // 3rd sub face
    RowVector4i temp2(oriPtId2, edgeToEdgePt(oriPtId2,oriPtId3)+numV, i+(numV*3), edgeToEdgePt(oriPtId1,oriPtId2)+numV);
    SF.row(cur) << temp2;
    cur++;

    // 4th sub face
    RowVector4i temp3(oriPtId3, edgeToEdgePt(oriPtId3,oriPtId0)+numV, i+(numV*3), edgeToEdgePt(oriPtId3,oriPtId2)+numV);
    SF.row(cur) << temp3;
    cur++;
  }


  unordered_map<int, vector<int>> ptToNewAdjFc;

  for (int i=0;i<SF.rows();i++){
    ptToNewAdjFc[SF(i,0)].push_back(i);

    ptToNewAdjFc[SF(i,1)].push_back(i);

    ptToNewAdjFc[SF(i,2)].push_back(i);

    ptToNewAdjFc[SF(i,3)].push_back(i);

  }

  MatrixXd newFacePoint = MatrixXd::Zero(SF.rows(),3);
  for(int i=0;i<SF.rows();i++){
    // op for original point
    int op0 = SF(i,0), op1 = SF(i,1), op2 = SF(i,2), op3 = SF(i,3);
    newFacePoint.row(i) << (SV.row(op0)+SV.row(op1)+SV.row(op2)+SV.row(op3))/4.0;
  }


  MatrixXd newV = MatrixXd::Zero(numV,3);
  for (int i=0;i<numV;i++){
    vector<int> nbrPt = ptToAdjPt[i];
    int N = static_cast<int>(nbrPt.size());
    RowVector3d PP = V.row(i);

    RowVector3d FSum(0.0,0.0,0.0), RSum(0.0,0.0,0.0);
    RowVector3d temp;

    for (int j=0;j<N;j++){
      temp = (V.row(nbrPt[j]) + PP)/2.0;
      RSum += temp;
    }

    for (int k=0;k<ptToNewAdjFc[i].size();k++){

      temp = newFacePoint.row(ptToNewAdjFc[i][k]);
      FSum += temp;
    }

    RowVector3d FF = FSum/N;
    RowVector3d RR = RSum/N;


    RowVector3d newP = (FF+2.0*RR+(N-3)*PP)/N;
    newV.row(i) << newP;
  }
  SV << newV,edgePtMat,facePoint;

  catmull_clark(MatrixXd(SV), MatrixXi(SF), num_iters-1, SV, SF);
}
//
