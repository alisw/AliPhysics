/*
 * AliFemtoDreamCascadeHist.h
 *
 *  Created on: Jan 11, 2018
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMCASCADEHIST_H_
#define ALIFEMTODREAMCASCADEHIST_H_
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TList.h"

class AliFemtoDreamCascadeHist {
 public:
  AliFemtoDreamCascadeHist();
  AliFemtoDreamCascadeHist(float mass, bool perRunnumber, int iRunMin,
                           int iRunMax);
  AliFemtoDreamCascadeHist(TString minimalBooking, float mass);
  virtual ~AliFemtoDreamCascadeHist();
  void FillCutCounter(int bin) {
    if (!fMinimalBooking)
      fCutCounter->Fill(bin);
  }
  ;
  void FillConfig(int iBin, float val) {
    if (!fMinimalBooking)
      fConfig->Fill(iBin, val);
  }
  ;
  void FillInvMass(int iBin, float mass) {
    if (!fMinimalBooking)
      fInvMass[iBin]->Fill(mass);
  }
  ;
  void FillInvMassPt(float Pt, float mass) {
    fInvMassPt->Fill(Pt, mass);
  }
  ;
  void FillInvMassPtv0(float Pt, float mass) {
    if (!fMinimalBooking)
      fInvMassPtv0->Fill(Pt, mass);
  }
  ;
  void FillInvMassLambda(int iBin, float mass) {
    if (!fMinimalBooking)
      fInvMassv0[iBin]->Fill(mass);
  }
  ;
  void FillXiPt(int iBin, float pT) {
    if (!fMinimalBooking)
      fXiPt[iBin]->Fill(pT);
  }
  ;
  void FillXiEta(int iBin, float Eta) {
    if (!fMinimalBooking)
      fXiEta[iBin]->Fill(Eta);
  }
  ;
  void FillXiPhi(int iBin, float phi) {
    if (!fMinimalBooking)
      fXiPhi[iBin]->Fill(phi);
  }
  ;
  void FillMomRapXi(int iBin, float rap, float pT) {
    if (!fMinimalBooking)
      fP_Y_Xi[iBin]->Fill(rap, pT);
  }
  ;
  void FillMomRapOmega(int iBin, float rap, float pT) {
    if (!fMinimalBooking)
      fP_Y_Omega[iBin]->Fill(rap, pT);
  }
  ;
  void FillDCAXiPrimVtx(int iBin, float DCA) {
    if (!fMinimalBooking)
      fDCAXi[iBin]->Fill(DCA);
  }
  ;
  void FillDCAXiDaug(int iBin, float DCA) {
    if (!fMinimalBooking)
      fDCAXiDaug[iBin]->Fill(DCA);
  }
  ;
  void FillMinDistPrimVtxBach(int iBin, float dist) {
    if (!fMinimalBooking)
      fMinDistVtxBach[iBin]->Fill(dist);
  }
  ;
  void FillCPAXi(int iBin, float cpa) {
    if (!fMinimalBooking)
      fCPAXi[iBin]->Fill(cpa);
  }
  ;
  void FillDecayLength(int iBin, float len) {
    if (!fMinimalBooking)
      fDecayLength[iBin]->Fill(len);
  }
  ;
  void Fillv0DecayLength(int iBin, float len) {
    if (!fMinimalBooking)
      fv0DecayLength[iBin]->Fill(len);
  }
  ;
  void FillTransverseRadiusXi(int iBin, float rad) {
    if (!fMinimalBooking)
      fTransRadiusXi[iBin]->Fill(rad);
  }
  ;
  void FillMaxDCAv0Daug(int iBin, float dca) {
    if (!fMinimalBooking)
      fv0MaxDCADaug[iBin]->Fill(dca);
  }
  ;
  void FillCPAv0(int iBin, float cpa) {
    if (!fMinimalBooking)
      fCPAv0[iBin]->Fill(cpa);
  }
  ;
  void FillCPAv0Xi(int iBin, float cpa) {
    if (!fMinimalBooking)
      fCPAv0Xi[iBin]->Fill(cpa);
  }
  ;
  void Fillv0Pt(int iBin, float v0pT) {
    if (!fMinimalBooking)
      fv0Pt[iBin]->Fill(v0pT);
  }
  ;
  void Fillv0Eta(int iBin, float v0eta) {
    if (!fMinimalBooking)
      fv0Eta[iBin]->Fill(v0eta);
  }
  ;
  void Fillv0Phi(int iBin, float v0phi) {
    if (!fMinimalBooking)
      fv0Phi[iBin]->Fill(v0phi);
  }
  ;
  void FillTransverseRadiusv0(int iBin, float rad) {
    if (!fMinimalBooking)
      fTransRadiusv0[iBin]->Fill(rad);
  }
  ;
  void FillMinDistPrimVtxv0(int iBin, float dist) {
    if (!fMinimalBooking)
      fMinDistVtxv0[iBin]->Fill(dist);
  }
  ;
  void FillMinDistPrimVtxv0DaugPos(int iBin, float dist) {
    if (!fMinimalBooking)
      fMinDistVtxv0DaugPos[iBin]->Fill(dist);
  }
  ;
  void FillMinDistPrimVtxv0DaugNeg(int iBin, float dist) {
    if (!fMinimalBooking)
      fMinDistVtxv0DaugNeg[iBin]->Fill(dist);
  }
  ;
  void FillPodolandski(int iBin, float alpha, float qt) {
    if (!fMinimalBooking)
      fPodolandski[iBin]->Fill(alpha, qt);
  }
  ;
  void FillInvMassPerRunNumber(int RunNumber, float mass) {
    if (!fMinimalBooking)
      fInvMassPerRunNumber->Fill(RunNumber, mass);
  }
  ;
  void SetName(TString name) {
    fHistList->SetName(name.Data());
  }
  ;
  TList *GetHistList() {
    return fHistList;
  }
  ;
 private:
  AliFemtoDreamCascadeHist &operator=(const AliFemtoDreamCascadeHist &obj);
  AliFemtoDreamCascadeHist(const AliFemtoDreamCascadeHist&);
  bool fMinimalBooking;
  TList *fHistList;
  TH1F *fCutCounter;
  TProfile *fConfig;
  TList *fCascadeQA[2];
  TH2F *fInvMassPt;
  TH2F *fInvMassPtv0;
  TH1F *fInvMass[2];
  TH1F *fInvMassv0[2];
  TH1F *fXiPt[2];
  TH1F *fXiEta[2];
  TH1F *fXiPhi[2];
  TH2F *fP_Y_Xi[2];
  TH2F *fP_Y_Omega[2];
  TH1F *fDCAXi[2];
  TH1F *fDCAXiDaug[2];
  TH1F *fMinDistVtxBach[2];
  TH1F *fCPAXi[2];
  TH1F *fDecayLength[2];
  TH1F *fv0DecayLength[2];
  TH1F *fTransRadiusXi[2];
  TH1F *fv0MaxDCADaug[2];
  TH1F *fCPAv0[2];
  TH1F *fCPAv0Xi[2];
  TH1F *fv0Pt[2];
  TH1F *fv0Eta[2];
  TH1F *fv0Phi[2];
  TH1F *fTransRadiusv0[2];
  TH1F *fMinDistVtxv0[2];
  TH1F *fMinDistVtxv0DaugPos[2];
  TH1F *fMinDistVtxv0DaugNeg[2];
  TH2F *fPodolandski[2];
  TH2F *fInvMassPerRunNumber;

ClassDef(AliFemtoDreamCascadeHist,2)
};

#endif /* ALIFEMTODREAMCASCADEHIST_H_ */
