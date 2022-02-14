/*
 * AliFemtoDreamEventHist.h
 *
 *  Created on: Nov 22, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMEVENTHIST_H_
#define ALIFEMTODREAMEVENTHIST_H_
#include "Rtypes.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
class AliFemtoDreamEventHist {
 public:
  AliFemtoDreamEventHist();
  AliFemtoDreamEventHist(bool centVsMultPlot);
  AliFemtoDreamEventHist(const AliFemtoDreamEventHist& hists);
  AliFemtoDreamEventHist& operator=(const AliFemtoDreamEventHist& hists);
  virtual ~AliFemtoDreamEventHist();
  void FillEvtCounter(int iBin) {
    fEvtCounter->Fill(iBin);
  }
  ;
  void FillCuts(int iBin, float val) {
    fCutConfig->Fill(iBin, val);
  }
  ;
  void FillEvtNCont(int i, float val) {
    fEvtNCont[i]->Fill(val);
  }
  ;
  void FillEvtVtxX(int i, float val) {
    fEvtVtxX[i]->Fill(val);
  }
  ;
  void FillEvtVtxY(int i, float val) {
    fEvtVtxY[i]->Fill(val);
  }
  ;
  void FillEvtVtxZ(int i, float val) {
    fEvtVtxZ[i]->Fill(val);
  }
  ;
  void FillEvtSpher(int i, float val) {
    fEvtSpher[i]->Fill(val);
  }
  ;
  void FillEvtSphero(int i, float val) {
    fEvtSphero[i]->Fill(val);
  }
  ;
  void FillMultSPD(int i, float val) {
    fMultDistSPD[i]->Fill(val);
  }
  ;
  void FillMultV0A(int i, float val) {
    fMultDistV0A[i]->Fill(val);
  }
  ;
  void FillMultV0C(int i, float val) {
    fMultDistV0C[i]->Fill(val);
  }
  ;
  void FillMultV0M(int i, float val) {
    fMultDistV0M[i]->Fill(val);
  }
  ;
  void FillMultPercentV0(int i, float val) {
    fMultPercentV0[i]->Fill(val);
  }
  void FillMultRef08(int i, float val) {
    fMultDistRef08[i]->Fill(val);
  }
  ;
  void FillSPDTrackletsLyVsCluster(int i, int ly, int spdTrkl, int spdCls) {
    if (ly == 0) {
      fSPDTrklClsLy0[i]->Fill(spdTrkl, spdCls);
    } else if (ly == 1) {
      fSPDTrklClsLy1[i]->Fill(spdTrkl, spdCls);
    }
  }
  void FillSPDTrackletsVsCluster(int i, int spdTrkl, int spdCls) {
    fSPDTrklClsLySum[i]->Fill(spdTrkl, spdCls);
  }
  void FillEvtVtxZTrackvsSPD(int i, float zVtxSPD, float zVtxTracks) {
    fSPDTrackZVtx[i]->Fill(zVtxSPD, zVtxTracks);
    fSPDTrkZVtxDispl[i]->Fill(TMath::Abs(zVtxSPD - zVtxTracks));
  }
  void FillMagneticField(int i, float bField) {
    fBField[i]->Fill(bField);
  }
  void FillVZEROtiming(int i, float tV0A, float tV0C) {
    fPileUpVZEROTime[i]->Fill(tV0A - tV0C, tV0A + tV0C);
  }
  void FillCentVsMultV0A(float cent, float mult) {
    if (fCentVsMultPlots)
      fCentVsV0A->Fill(cent, mult);
  }
  void FillCentVsMultV0M(float cent, float mult) {
    if (fCentVsMultPlots)
      fCentVsV0M->Fill(cent, mult);
  }
  void FillCentVsMultV0C(float cent, float mult) {
    if (fCentVsMultPlots)
      fCentVsV0C->Fill(cent, mult);
  }
  void FillCentVsMultRef(float cent, float mult) {
    if (fCentVsMultPlots)
      fCentVsRefMult->Fill(cent, mult);
  }

  TList *GetHistList() {
    return fEventCutList;
  }
  ;
 private:
  TList *fEventCutList;     //!
  TList *fEvtCutQA[2];      //!
  TH1F *fEvtCounter;        //!
  TProfile *fCutConfig;     //!
  TH1F *fEvtNCont[2];       //!
  TH1F *fEvtVtxX[2];        //!
  TH1F *fEvtVtxY[2];        //!
  TH1F *fEvtVtxZ[2];        //!
  TH1F *fMultDistSPD[2];    //!
  TH1F *fMultDistV0A[2];    //!
  TH1F *fMultDistV0C[2];    //!
  TH1F *fMultDistV0M[2];    //!
  TH1F *fMultPercentV0[2];  //!
  TH1F *fMultDistRef08[2];  //!
  TH2F *fSPDTrklClsLy0[2];     //!
  TH2F *fSPDTrklClsLy1[2];     //!
  TH2F *fSPDTrklClsLySum[2];     //!
  TH2F *fSPDTrackZVtx[2];   //!
  TH1F *fSPDTrkZVtxDispl[2];   //!
  TH1F *fBField[2];         //!
  TH1F *fEvtSpher[2];         //!
  TH2F* fPileUpVZEROTime[2]; //!
  bool fCentVsMultPlots;    //!
  TH1F *fEvtSphero[2];      //!
  TH2F *fCentVsV0A;         //!
  TH2F *fCentVsV0M;         //!
  TH2F *fCentVsV0C;         //!
  TH2F *fCentVsRefMult;     //!
ClassDef(AliFemtoDreamEventHist,7)
};

#endif /* ALIFEMTODREAMEVENTHIST_H_ */
