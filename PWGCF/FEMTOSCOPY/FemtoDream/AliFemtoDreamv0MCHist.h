/*
 * AliFemtoDreamv0MCHist.h
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMV0MCHIST_H_
#define ALIFEMTODREAMV0MCHIST_H_
#include "Rtypes.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliFemtoDreamBasePart.h"
class AliFemtoDreamv0MCHist {
 public:
  AliFemtoDreamv0MCHist();
  AliFemtoDreamv0MCHist(
      int MassNBins,double MassMin,double MassMax,bool contribSplitting,
      bool CPADist);
  virtual ~AliFemtoDreamv0MCHist();
  void FillMCCorr(double pT){fMCCorrPt->Fill(pT);};
  void FillMCIdent(double pT){fMCIdentPt->Fill(pT);};
  void FillMCGen(double pT){fMCGenPt->Fill(pT);};
  void FillMCCont(double pT){fMCContPt->Fill(pT);};
  void FillMCPrimary(double pT){fMCPrimaryPt->Fill(pT);};
  void FillMCMaterial(double pT){fMCMaterialPt->Fill(pT);};
  void FillMCFeeddown(double pT, double pdg){
    fMCFeeddownWeakPt->Fill(pT, pdg);
  };
  void FillMCpT(int i, double val){fMCpTDist[i]->Fill(val);};
  void FillMCEta(int i, double val){fMCetaDist[i]->Fill(val);};
  void FillMCPhi(int i, double val){fMCphiDist[i]->Fill(val);};
  void FillMCDCAVtxX(int i, double pT, double val){
    fMCDecayVtxv0X[i]->Fill(pT, val);
  };
  void FillMCDCAVtxY(int i, double pT, double val){
    fMCDecayVtxv0Y[i]->Fill(pT, val);};
  void FillMCDCAVtxZ(int i, double pT, double val){
    fMCDecayVtxv0Z[i]->Fill(pT, val);};
  void FillMCTransverseRadius(int i, double pT, double val){
    fMCTransverseRadius[i]->Fill(pT, val);
  };
  void FillMCDCAPosDaugPrimVtx(int i, double pT, double val){
    fMCDCAPosDaugToPV[i]->Fill(pT, val);};
  void FillMCDCANegDaugPrimVtx(int i, double pT, double val){
    fMCDCANegDaugToPV[i]->Fill(pT, val);};
  void FillMCDCADaugVtx(int i, double pT, double val){
    fMCDCADaugToVtx[i]->Fill(pT, val);};
  void FillMCCosPoint(int i, double pT, double val){
    fMCCosPointing[i]->Fill(pT, val);};
  void FillMCInvMass(int i, double massVal){
    fMCInvMass[i]->Fill(massVal);};
  void FillMCCPAPtBins(
      AliFemtoDreamBasePart::PartOrigin org, double pT,double cpa);
  TList *GetHistList(){return fMCList;};
  TString ClassName() {return "AliFemtoDreamv0MCHist";};
 private:
  TList *fMCList;
  TList *fCPAPlots;
  TList *fMCQAPlots[4];
  TH1F *fMCCorrPt;
  TH1F *fMCIdentPt;
  TH1F *fMCGenPt;
  TH1F *fMCContPt;
  TH1F *fMCPrimaryPt;
  TH1F *fMCMaterialPt;
  TH2F *fMCFeeddownWeakPt;
  TH1F *fMCpTDist[4];
  TH1F *fMCetaDist[4];
  TH1F *fMCphiDist[4];
  TH2F *fMCDecayVtxv0X[4];
  TH2F *fMCDecayVtxv0Y[4];
  TH2F *fMCDecayVtxv0Z[4];
  TH2F *fMCTransverseRadius[4];
  TH2F *fMCDCAPosDaugToPV[4];
  TH2F *fMCDCANegDaugToPV[4];
  TH2F *fMCDCADaugToVtx[4];
  TH2F *fMCCosPointing[4];
  TH1F *fMCInvMass[4];
  TH2F *fMCPrimCPAPtBins;
  TH2F *fMCMaterialCPAPtBins;
  TH2F *fMCSecondaryCPAPtBins;
  TH2F *fMCContCPAPtBins;

  ClassDef(AliFemtoDreamv0MCHist,1)
};

#endif /* ALIFEMTODREAMV0MCHIST_H_ */
