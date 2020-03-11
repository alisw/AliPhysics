/*
 * AliFemtoDreamTrackMCHist.h
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMTRACKMCHIST_H_
#define ALIFEMTODREAMTRACKMCHIST_H_
#include "AliFemtoDreamBasePart.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TString.h"

class AliFemtoDreamTrackMCHist {
 public:
  AliFemtoDreamTrackMCHist();
  AliFemtoDreamTrackMCHist(bool contribSplitting, bool DCADist,
                           bool DoMultBinning = false,
                           bool checkMother = false,
			   float pTmin = 0.5, float pTmax = 4.05);
  virtual ~AliFemtoDreamTrackMCHist();
  void FillMCDCAXYPtBins(AliFemtoDreamBasePart::PartOrigin org, int PDGCodeMoth,
                         float pT, float dcaxy, int multiplicity);
  void FillMCCorr(float pT) {
    fMCCorrPt->Fill(pT);
  }
  ;
  void FillMCIdent(float pT) {
    fMCIdentPt->Fill(pT);
  }
  ;
  void FillMCGen(float pT) {
    fMCGenPt->Fill(pT);
  }
  ;
  void FillMCCont(float pT) {
    fMCContPt->Fill(pT);
  }
  ;
  void FillMCUnkn(float pT) {
    fMCUnknownPt->Fill(pT);
  }
  ;
  void FillMCPrimary(float pT) {
    fMCPrimaryPt->Fill(pT);
  }
  ;
  void FillMCMaterial(float pT) {
    fMCMaterialPt->Fill(pT);
  }
  ;
  void FillMCFeeddown(float pT, float pdg) {
    fMCFeeddownWeakPt->Fill(pT, pdg);
  }
  ;
  void FillMCpTPCCut(int i, float pTPC) {
    fMCpTPCDist[i]->Fill(pTPC);
  }
  ;
  void FillMCetaCut(int i, float eta) {
    fMCetaDist[i]->Fill(eta);
  }
  ;
  void FillMCphiCut(int i, float phi) {
    fMCphiDist[i]->Fill(phi);
  }
  ;
  void FillMCTPCclsCut(int i, float pT, float nCls) {
    fMCTPCCls[i]->Fill(pT, nCls);
  }
  ;
  void FillMCDCAxyCut(int i, float pT, float dcaxy) {
    fMCDCAxy[i]->Fill(pT, dcaxy);
  }
  ;
  void FillMCDCAzCut(int i, float pT, float dcaz) {
    fMCDCAz[i]->Fill(pT, dcaz);
  }
  ;
  void FillMCTPCCrossedRowCut(int i, float pT, float Crossed) {
    fMCTPCCrossedRows[i]->Fill(pT, Crossed);
  }
  ;
  void FillMCTPCRatioCut(int i, float pT, float ratio) {
    fMCTPCRatio[i]->Fill(pT, ratio);
  }
  ;
  void FillMCTPCdedx(int i, float mom, float dedx) {
    fMCTPCdedx[i]->Fill(mom, dedx);
  }
  ;
  void FillMCTOFbeta(int i, float mom, float beta) {
    fMCTOFbeta[i]->Fill(mom, beta);
  }
  ;
  void FillMCNSigTPC(int i, float mom, float nSigTPC) {
    fMCNSigTPC[i]->Fill(mom, nSigTPC);
  }
  ;
  void FillMCNSigTOF(int i, float mom, float nSigTOF) {
    fMCNSigTOF[i]->Fill(mom, nSigTOF);
  }
  ;
  void FillMCPtResolution(float pTTrue, float pTReco);
  void FillMCThetaResolution(float ThetaTrue, float ThetaReco, float pTTrue);
  void FillMCPhiResolution(float PhiTrue, float PhiReco, float pTTrue);
  void FillMultiplicityHistos(int multiplicity, float pT, float dcaxy,
                              TH2F *histo1, TH2F *histo2, TH2F *histo3);
  void FillMCMother(float pT, int pdg) { fHistMCMother->Fill(pT, std::abs(pdg));
    fHistMCMotherPDG->Fill(std::abs(pdg));
  }

  void SetMultRangeLow(int range) {
    fMultRangeLow = range;
  }
  void SetMultRangeHigh(int range) {
    fMultRangeHigh = range;
  }
  void SetName(TString name) {
    fMCList->SetName(name.Data());
  }
  ;
  TList *GetHistList() const {
    return fMCList;
  }
  ;
  TString ClassName() {
    return "AliFemtoDreamTrackMCHist";
  }
  ;
 private:
  AliFemtoDreamTrackMCHist &operator=(const AliFemtoDreamTrackMCHist &obj);
  AliFemtoDreamTrackMCHist(const AliFemtoDreamTrackMCHist&);
  float fpTbins;                 //!
  float fMultRangeLow;			 //!
  float fMultRangeHigh;			 //!
  bool fDoSplitting;              //!
  bool fDoDCAPlots;               //!
  bool fDoMultiplicityBinning;    //!
  float fpTmin;                   //!
  float fpTmax;                   //!

  TList *fMCList;                 //!
  TList *fMCQAPlots[4];           //!
  TList *fDCAPlots;               //!
  TH1F *fMCCorrPt;                //!
  TH1F *fMCIdentPt;               //!
  TH1F *fMCGenPt;                 //!
  TH1F *fMCContPt;                //!
  TH1F *fMCUnknownPt;             //!
  TH1F *fMCPrimaryPt;             //!
  TH1F *fMCMaterialPt;            //!
  TH2F *fMCFeeddownWeakPt;        //!
  TH2F *fHistMCMother;            //!
  TH1I *fHistMCMotherPDG;         //!

  TH2F *fMCPrimDCAXYPtBins;           //!
  TH2F *fMCMaterialDCAXYPtBins;       //!
  TH2F *fMCSecondaryDCAXYPtBins;      //!
  TH2F *fMCSecLambdaDCAXYPtBins;      //!
  TH2F *fMCSecSigmaDCAXYPtBins;       //!
  TH2F *fMCSecSigmaPlusDCAXYPtBins;   //!
  TH2F *fMCSecSigmaMinusDCAXYPtBins;  //!
  TH2F *fMCSecXiDCAXYPtBins;  	      //!
  TH2F *fMCSecOmegaDCAXYPtBins;       //!
  TH2F *fMCSecKlongDCAXYPtBins;       //!
  TH2F *fMCSecKshortDCAXYPtBins;      //!
  TH2F *fMCSecKchDCAXYPtBins;         //!

  TH2F *fMCPrimDCAXYPtBinsMult[3];       //!
  TH2F *fMCMaterialDCAXYPtBinsMult[3];   //!
  TH2F *fMCSecondaryDCAXYPtBinsMult[3];	 //!
  TH2F *fMCSecLambdaDCAXYPtBinsMult[3];  //!
  TH2F *fMCSecSigmaDCAXYPtBinsMult[3];   //!

  TH1F *fMCpTPCDist[4];           //!
  TH1F *fMCetaDist[4];            //!
  TH1F *fMCphiDist[4];            //!
  TH2F *fMCTPCCls[4];             //!
  TH2F *fMCDCAxy[4];              //!
  TH2F *fMCDCAz[4];               //!

  TH2F *fMCTPCCrossedRows[4];	  //!
  TH2F *fMCTPCRatio[4];           //!

  TH2F *fMCTPCdedx[4];            //!
  TH2F *fMCTOFbeta[4];            //!
  TH2F *fMCNSigTPC[4];            //!
  TH2F *fMCNSigTOF[4];            //!

  TH2F *fPtResolution;            //!
  TH2F *fThetaResolution;         //!
  TH2F *fPhiResolution;           //!
ClassDef(AliFemtoDreamTrackMCHist,4)
  ;
};

#endif /* ALIFEMTODREAMTRACKMCHIST_H_ */
