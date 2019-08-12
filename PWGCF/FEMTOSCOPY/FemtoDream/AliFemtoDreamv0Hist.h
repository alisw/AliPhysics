/*
 * AliFemtoDreamv0Hist.h
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMV0HIST_H_
#define ALIFEMTODREAMV0HIST_H_
#include "Rtypes.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

class AliFemtoDreamv0Hist {
 public:
  AliFemtoDreamv0Hist();
  AliFemtoDreamv0Hist(int MassNBins, float MassMin, float MassMax,
                      bool CPAPlots, bool perRunnumber, int iRunMin,
                      int iRunMax);
  AliFemtoDreamv0Hist(TString MinimalBooking, int MassNBins, float MassMin,
                      float MassMax);
  AliFemtoDreamv0Hist(const AliFemtoDreamv0Hist& hists);
  AliFemtoDreamv0Hist& operator=(const AliFemtoDreamv0Hist& hists);
  virtual ~AliFemtoDreamv0Hist();
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
  void FillOnFlyStatus(int i, float val) {
    if (!fMinimalBooking)
      fOnFly[i]->Fill(val);
  }
  ;
  void FillpTCut(int i, float pT) {
    if (!fMinimalBooking)
      fpTDist[i]->Fill(pT);
  }
  ;
  void FillPhi(int i, float phi) {
    if (!fMinimalBooking)
      fPhiDist[i]->Fill(phi);
  }
  ;
  void FillEtaCut(int i, float eta) {
    fetaDist[i]->Fill(eta);
  }
  ;
  void Fillv0DecayVtxXCut(int i, float vtx) {
    if (!fMinimalBooking)
      fDecayVtxv0X[i]->Fill(vtx);
  }
  ;
  void Fillv0DecayVtxYCut(int i, float vtx) {
    if (!fMinimalBooking)
      fDecayVtxv0Y[i]->Fill(vtx);
  }
  ;
  void Fillv0DecayVtxZCut(int i, float vtx) {
    if (!fMinimalBooking)
      fDecayVtxv0Z[i]->Fill(vtx);
  }
  ;
  void FillTransverRadiusCut(int i, float rad) {
    if (!fMinimalBooking)
      fTransRadius[i]->Fill(rad);
  }
  ;
  void FillDCAPosDaugToPrimVtxCut(int i, float dca) {
    if (!fMinimalBooking)
      fDCAPosDaugToPrimVtx[i]->Fill(dca);
  }
  ;
  void FillDCANegDaugToPrimVtxCut(int i, float dca) {
    if (!fMinimalBooking)
      fDCANegDaugToPrimVtx[i]->Fill(dca);
  }
  ;
  void FillDCADaugTov0VtxCut(int i, float dca) {
    if (!fMinimalBooking)
      fDCADaugToVtx[i]->Fill(dca);
  }
  ;
  void FillCPACut(int i, float cpa) {
    if (!fMinimalBooking)
      fCPA[i]->Fill(cpa);
  }
  ;
  void FillInvMass(int i, float mass) {
    if (!fMinimalBooking)
      fInvMass[i]->Fill(mass);
  }
  ;
  void FillInvMassBefKaonRej(float mass) {
    if (!fMinimalBooking)
      fInvMassBefKaonRej->Fill(mass);
  }
  ;
  void FillInvMassKaon(float mass) {
    if (!fMinimalBooking)
      fInvMassKaon->Fill(mass);
  }
  ;
  void Fillv0MassDist(float mass) {
    if (!fMinimalBooking)
      fInvMassBefSelection->Fill(mass);
  }
  ;
  void FillInvMassPtBins(float pT, float mass) {
    fInvMassPt->Fill(pT, mass);
  }
  ;
  void FillArmenterosPodolandski(int i, float qT, float alpha) {
    if (!fMinimalBooking)
      fArmenterosPodolandski[i]->Fill(qT, alpha);
  }
  void FillCPAPtBins(float pT, float cpa, int multiplicity);
  void FillInvMassPerRunNumber(int RunNumber, float mass) {
    if (!fMinimalBooking)
      fInvMassPerRunNumber->Fill(RunNumber, mass);
  }
  ;
  void SetName(TString name) {
    fHistList->SetName(name.Data());
  }
  ;
  void SetMultRangeLow(int range) {
    fMultRangeLow = range;
  }
  void SetMultRangeHigh(int range) {
    fMultRangeHigh = range;
  }
  TList *GetHistList() {
    return fHistList;
  }
  ;
 private:
  bool fMinimalBooking;
  TList *fHistList;
  TList *fv0CutQA[2];
  float fMultRangeLow;  //!
  float fMultRangeHigh;  //!
  TProfile *fConfig;
  TH1F *fCutCounter;
  TH1F *fOnFly[2];
  TH1F *fpTDist[2];
  TH1F *fPhiDist[2];
  TH1F *fetaDist[2];
  TH1F *fDecayVtxv0X[2];
  TH1F *fDecayVtxv0Y[2];
  TH1F *fDecayVtxv0Z[2];
  TH1F *fTransRadius[2];
  TH1F *fDCAPosDaugToPrimVtx[2];
  TH1F *fDCANegDaugToPrimVtx[2];
  TH1F *fDCADaugToVtx[2];
  TH1F *fCPA[2];
  TH1F *fInvMass[2];
  TH1F *fInvMassBefKaonRej;
  TH1F *fInvMassKaon;
  TH1F *fInvMassBefSelection;
  TH2F *fInvMassPt;
  TH2F *fArmenterosPodolandski[2];
  TH2F *fCPAPtBins;
  TH2F *fCPAPtBinsMult[3];  //!
  TH2F *fInvMassPerRunNumber;ClassDef(AliFemtoDreamv0Hist,3)
};

#endif /* ALIFEMTODREAMV0HIST_H_ */
