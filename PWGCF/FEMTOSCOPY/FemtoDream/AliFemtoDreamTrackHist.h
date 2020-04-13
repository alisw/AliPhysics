/*
 * AliFemtoDreamTrackHist.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMTRACKHIST_H_
#define ALIFEMTODREAMTRACKHIST_H_
#include "AliPIDResponse.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TProfile.h"

class AliFemtoDreamTrackHist {
 public:
  AliFemtoDreamTrackHist();
  AliFemtoDreamTrackHist(bool DCADist, bool CombSig, bool TOFM, float pTmin = 0.5, float pTmax = 4.05);
  AliFemtoDreamTrackHist(TString MinimalBooking);
  virtual ~AliFemtoDreamTrackHist();
  void FillConfig(int iBin, float val) {
    if (!fMinimalBooking)
      fConfig->Fill(iBin, val);
  }
  ;
  void FillTrackCounter(int iBin) {
    if (!fMinimalBooking)
      fCutCounter->Fill(iBin);
  }
  ;
  void FillpTCut(int i, float pT) {
    if ((!fMinimalBooking) || (i == 1))
      fpTDist[i]->Fill(pT);
  }
  ;
  void FillpCut(int i, float p) {
    if (!fMinimalBooking)
      fpDist[i]->Fill(p);
  }
  ;
  void FillpTPCCut(int i, float pTPC) {
    if (!fMinimalBooking)
      fpTPCDist[i]->Fill(pTPC);
  }
  ;
  void FilletaCut(int i, float eta) {
    if (!fMinimalBooking)
      fetaDist[i]->Fill(eta);
  }
  ;
  void FillphiCut(int i, float phi) {
    if (!fMinimalBooking)
      fphiDist[i]->Fill(phi);
  }
  ;
  void FillTPCclsCut(int i, float nCls) {
    if (!fMinimalBooking)
      fTPCCls[i]->Fill(nCls);
  }
  ;
  void FillDCAxyCut(int i, float pT, float dcaxy) {
    if (!fMinimalBooking)
      fDCAxy[i]->Fill(pT, dcaxy);
  }
  ;
  void FillDCAzCut(int i, float pT, float dcaz) {
    if (!fMinimalBooking)
      fDCAz[i]->Fill(pT, dcaz);
  }
  ;
  void FillDCAxyPropCut(int i, float pT, float dcaxy) {
    if (!fMinimalBooking)
      fDCAxyProp[i]->Fill(pT, dcaxy);
  }
  ;
  void FillDCAzPropCut(int i, float pT, float dcaz) {
    if (!fMinimalBooking)
      fDCAzProp[i]->Fill(pT, dcaz);
  }
  ;
  void FillTPCCrossedRowCut(int i, float Crossed) {
    if (!fMinimalBooking)
      fTPCCrossedRows[i]->Fill(Crossed);
  }
  ;
  void FillTPCRatioCut(int i, float ratio) {
    if (!fMinimalBooking)
      fTPCRatio[i]->Fill(ratio);
  }
  ;
  void FillTPCClsS(int i, float TPCClsS) {
    if (!fMinimalBooking)
      fTPCClsS[i]->Fill(TPCClsS);
  }
  ;
  void FillHasSharedClsITS(int i, int layer, int yesno) {
    if (!fMinimalBooking)
      fShrdClsITS[i]->Fill(layer, yesno);
  }
  ;
  void FillITSdedx(int i, float mom, float dedx) {
    if (!fMinimalBooking)
      fITSdedx[i]->Fill(mom, dedx);
  }
  ;
  void FillTPCdedx(int i, float mom, float dedx) {
    if (!fMinimalBooking)
      fTPCdedx[i]->Fill(mom, dedx);
  }
  ;
  void FillTOFbeta(int i, float mom, float beta) {
    if (!fMinimalBooking)
      fTOFbeta[i]->Fill(mom, beta);
  }
  ;
  void FillNSigITS(int i, float mom, float nSigITS) {
    if (!fMinimalBooking)
      fNSigITS[i]->Fill(mom, nSigITS);
  }
  ;
  void FillNSigITSMod(int i, float mom,float nSigITS){
    if (!fMinimalBooking)
      fNSigITSMod[i]->Fill(mom, TMath::Abs(nSigITS));
  }
  ;
  void FillNSigTPC(int i, float mom, float nSigTPC) {
    if (!fMinimalBooking)
      fNSigTPC[i]->Fill(mom, nSigTPC);
  }
  ;
  void FillNSigTPCMod(int i, float mom,float nSigTPC){
    if (!fMinimalBooking)
      fNSigTPCMod[i]->Fill(mom, TMath::Abs(nSigTPC));
  }
  ;
  void FillNSigTOF(int i, float mom, float nSigTOF) {
    if (!fMinimalBooking)
      fNSigTOF[i]->Fill(mom, nSigTOF);
  }
  ;
  void FillITSStatus(int i, AliPIDResponse::EDetPidStatus statusITS) {
    if (!fMinimalBooking)
      fITSStatus[i]->Fill(statusITS);
  }
  ;
  void FillTPCStatus(int i, AliPIDResponse::EDetPidStatus statusTPC) {
    if (!fMinimalBooking)
      fTPCStatus[i]->Fill(statusTPC);
  }
  ;
  void FillTOFStatus(int i, AliPIDResponse::EDetPidStatus statusTOF) {
    if (!fMinimalBooking)
      fTOFStatus[i]->Fill(statusTOF);
  }
  ;
  void FillNSigComb(float pT, float nSigTPC, float nSigTOF);

  void FillNSigComITSTPC(int i, float mom,float nSigITS, float nSigTPC){
    if (!fMinimalBooking)
      fNSigComITSTPC[i]->Fill(mom, TMath::Sqrt(nSigITS*nSigITS+nSigTPC*nSigTPC));
  }
  ;
  void FillNSigComTPCTOF(int i, float mom,float nSigTPC, float nSigTOF){
    if (!fMinimalBooking)
      fNSigComTPCTOF[i]->Fill(mom, TMath::Sqrt(nSigTPC*nSigTPC+nSigTOF*nSigTOF));
  }
  ;


  void FillDCAXYPtBins(float pT, float dcaxy, int multiplicity);

  void FillTOFMass(float mom, float beta) {
    if (fTOFMass) {
      if (beta > 0 && beta <= 1) {
        fTOFMass->Fill(mom, mom / beta * TMath::Sqrt(1 - beta * beta));
      } else {
        fTOFMass->Fill(0., -999.);
      }
    }
  }

  void FillTPCClsCPileUp(int i, int iCrit, float TPCClsC) {
    if (!fMinimalBooking)
      fTPCClsCPiluUp[i]->Fill(iCrit, TPCClsC);
  }
  void FillITSSharedPileUp(int i, int iCrit, int yesno) {
    if (!fMinimalBooking)
      fITShrdClsPileUp[i]->Fill(iCrit, yesno);
  }
  void FillTrackChiSquare(int i, float pT, float chi2) {
    if (!fMinimalBooking)
      fTrackChi2[i]->Fill(pT, chi2);
  }

  void SetMultRangeLow(int range) {
    fMultRangeLow = range;
  }
  void SetMultRangeHigh(int range) {
    fMultRangeHigh = range;
  }

  void SetName(TString name) {
    fHistList->SetName(name.Data());
  }
  ;
  TList *GetHistList() {
    return fHistList;
  }
  ;
 private:
  AliFemtoDreamTrackHist &operator=(const AliFemtoDreamTrackHist &obj);
  AliFemtoDreamTrackHist(const AliFemtoDreamTrackHist&);
  bool fMinimalBooking;     //!
  int fMultRangeLow;		//!
  int fMultRangeHigh;		//!
  float fpTmin;              //!
  float fpTmax;              //!
  TList *fHistList;         //!
  TList *fTrackCutQA[2];    //!
  TProfile *fConfig;        //!
  TH1F *fCutCounter;        //!
  TH1F *fpTDist[2];         //!
  TH1F *fpDist[2];          //!
  TH1F *fpTPCDist[2];       //!
  TH1F *fetaDist[2];        //!
  TH1F *fphiDist[2];        //!
  TH1F *fTPCCls[2];         //!
  TH2F *fShrdClsITS[2];     //!
  TH2F *fDCAxy[2];          //!
  TH2F *fDCAz[2];           //!
  TH2F *fDCAxyProp[2];      //!
  TH2F *fDCAzProp[2];       //!
  TH2F *fDCAXYPtBins;       //!
  TH2F *fTOFMass;           //!
  TH2F *fDCAXYPtBinsMult[3];//!
  TH1F *fTPCCrossedRows[2]; //!
  TH1F *fTPCRatio[2];       //!
  TH1F *fTPCClsS[2];        //!
  TH2F *fTrackChi2[2];      //!
  TH2F *fITSdedx[2];        //!
  TH2F *fTPCdedx[2];        //!
  TH2F *fTOFbeta[2];        //!
  TH2F *fNSigITS[2];        //!
  TH2F *fNSigITSMod[2];     //!
  TH2F *fNSigTPC[2];        //!
  TH2F *fNSigTPCMod[2];     //!
  TH2F *fNSigTOF[2];        //!
  TH1F *fITSStatus[2];      //!
  TH1F *fTPCStatus[2];      //!
  TH1F *fTOFStatus[2];      //!
  TH3F *fNSigCom;           //!
  TH2F  *fNSigComITSTPC[2];   //!
  TH2F  *fNSigComTPCTOF[2];   //!
  TH2F *fTPCClsCPiluUp[2];  //!
  TH2F *fITShrdClsPileUp[2];  //!
ClassDef(AliFemtoDreamTrackHist,5)
  ;
};

#endif /* ALIFEMTODREAMTRACKHIST_H_ */
