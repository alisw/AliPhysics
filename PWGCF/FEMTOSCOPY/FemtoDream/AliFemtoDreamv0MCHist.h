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
  AliFemtoDreamv0MCHist(int MassNBins, float MassMin, float MassMax,
                        bool contribSplitting, bool CPADist,
                        bool DoMultBinning = false, bool checkMother = false);
  virtual ~AliFemtoDreamv0MCHist();
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
  void FillMCpT(int i, float val) {
    fMCpTDist[i]->Fill(val);
  }
  ;
  void FillMCEta(int i, float val) {
    fMCetaDist[i]->Fill(val);
  }
  ;
  void FillMCPhi(int i, float val) {
    fMCphiDist[i]->Fill(val);
  }
  ;
  void FillMCDCAVtxX(int i, float pT, float val) {
    fMCDecayVtxv0X[i]->Fill(pT, val);
  }
  ;
  void FillMCDCAVtxY(int i, float pT, float val) {
    fMCDecayVtxv0Y[i]->Fill(pT, val);
  }
  ;
  void FillMCDCAVtxZ(int i, float pT, float val) {
    fMCDecayVtxv0Z[i]->Fill(pT, val);
  }
  ;
  void FillMCTransverseRadius(int i, float pT, float val) {
    fMCTransverseRadius[i]->Fill(pT, val);
  }
  ;
  void FillMCDCAPosDaugPrimVtx(int i, float pT, float val) {
    fMCDCAPosDaugToPV[i]->Fill(pT, val);
  }
  ;
  void FillMCDCANegDaugPrimVtx(int i, float pT, float val) {
    fMCDCANegDaugToPV[i]->Fill(pT, val);
  }
  ;
  void FillMCDCADaugVtx(int i, float pT, float val) {
    fMCDCADaugToVtx[i]->Fill(pT, val);
  }
  ;
  void FillMCCosPoint(int i, float pT, float val) {
    fMCCosPointing[i]->Fill(pT, val);
  }
  ;
  void FillMCInvMass(int i, float massVal) {
    fMCInvMass[i]->Fill(massVal);
  }
  ;
  void FillMCCPAPtBins(AliFemtoDreamBasePart::PartOrigin org, float pT,
                       float cpa, int multiplicity);
  void FillMultiplicityHistos(int multiplicity, float pT, float cpa,
                              TH2F *histo1, TH2F *histo2, TH2F *histo3);
  void FillMCBachDCAToPV(int i, float pT, float val) {
    fMCBachDCAToPV[i]->Fill(pT, val);
  }
  ;
  void FillMCv0DecayLength(int i, float pT, float val) {
    fMCv0DecayLength[i]->Fill(pT, val);
  }
  ;
  void FillMCv0CPA(int i, float pT, float val) {
    fMCv0CPA[i]->Fill(pT, val);
  }
  ;
  void FillMCDecayLength(int i, float pT, float val) {
    fMCDecayLength[i]->Fill(pT, val);
  }
  ;
  void FillMCXiRapidity(int i, float p, float val) {
    fMCXiRapidity[i]->Fill(val, p);
  }
  ;
  void FillMCOmegaRapidity(int i, float p, float val) {
    fMCOmegaRapidity[i]->Fill(val, p);
  }
  ;
  void FillMCPodolanski(int i, float pT, float alpha) {
    fMCPodolanski[i]->Fill(alpha, pT);
  }
  ;
  void FillMCXiInvMass(int i, float pT, float mass) {
    fMCXiInvMass[i]->Fill(pT, mass);
  }
  ;
  void FillMCOmegaInvMass(int i, float pT, float mass) {
    fMCOmegaInvMass[i]->Fill(pT, mass);
  }
  ;
  void FillMCPtResolution(float pTTrue, float pTReco);
  void FillMCThetaResolution(float ThetaTrue, float ThetaReco, float pTTrue);
  void FillMCPhiResolution(float PhiTrue, float PhiReco, float pTTrue);
  void FillMCMother(float pT, int pdg) {
	fHistMCMother->Fill(pT, std::abs(pdg));
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
  TList *GetHistList() {
    return fMCList;
  }
  ;
  TString ClassName() {
    return "AliFemtoDreamv0MCHist";
  }
  ;
 private:
  AliFemtoDreamv0MCHist &operator=(const AliFemtoDreamv0MCHist &obj);
  AliFemtoDreamv0MCHist(const AliFemtoDreamv0MCHist&);
  TList *fMCList;
  TList *fCPAPlots;
  TList *fMCQAPlots[5];
  float fMultRangeLow;  //!
  float fMultRangeHigh;  //!
  bool fDoMultiplicityBinning; //!
  TH1F *fMCCorrPt;
  TH1F *fMCIdentPt;
  TH1F *fMCGenPt;
  TH1F *fMCContPt;
  TH1F *fMCPrimaryPt;
  TH1F *fMCMaterialPt;
  TH2F *fMCFeeddownWeakPt;
  TH2F *fHistMCMother;            //!
  TH1I *fHistMCMotherPDG;         //!
  TH1F *fMCpTDist[5];
  TH1F *fMCetaDist[5];
  TH1F *fMCphiDist[5];
  TH2F *fMCDecayVtxv0X[5];
  TH2F *fMCDecayVtxv0Y[5];
  TH2F *fMCDecayVtxv0Z[5];
  TH2F *fMCTransverseRadius[5];
  TH2F *fMCDCAPosDaugToPV[5];
  TH2F *fMCDCANegDaugToPV[5];
  TH2F *fMCDCADaugToVtx[5];
  TH2F *fMCCosPointing[5];
  TH1F *fMCInvMass[5];
  TH2F *fMCBachDCAToPV[5];
  TH2F *fMCv0DecayLength[5];
  TH2F *fMCv0CPA[5];
  TH2F *fMCDecayLength[5];
  TH2F *fMCXiRapidity[5];
  TH2F *fMCOmegaRapidity[5];
  TH2F *fMCPodolanski[5];
  TH2F *fMCXiInvMass[5];
  TH2F *fMCOmegaInvMass[5];
  TH2F *fMCPrimCPAPtBins;
  TH2F *fMCMaterialCPAPtBins;
  TH2F *fMCSecondaryCPAPtBins;
  TH2F *fMCContCPAPtBins;
  TH2F *fMCPrimCPAPtBinsMult[3];
  TH2F *fMCMaterialCPAPtBinsMult[3];
  TH2F *fMCSecondaryCPAPtBinsMult[3];
  TH2F *fMCContCPAPtBinsMult[3];
  TH2F *fPtResolution;            //!
  TH2F *fThetaResolution;         //!
  TH2F *fPhiResolution;           //!

ClassDef(AliFemtoDreamv0MCHist,4)
};

#endif /* ALIFEMTODREAMV0MCHIST_H_ */
