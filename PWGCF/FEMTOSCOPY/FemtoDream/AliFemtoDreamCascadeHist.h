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
  AliFemtoDreamCascadeHist(float mass,bool perRunnumber, int iRunMin, int iRunMax);
  AliFemtoDreamCascadeHist(TString minimalBooking,float mass);
  virtual ~AliFemtoDreamCascadeHist();
  void FillCutCounter(int bin) {if(!fMinimalBooking)fCutCounter->Fill(bin);};
  void FillConfig(int iBin,float val){if(!fMinimalBooking)fConfig->Fill(iBin,val);};
  void FillInvMassXi(int iBin,float mass){if(!fMinimalBooking)fInvMass[iBin]->Fill(mass);};
  void FillInvMassPtXi(float Pt,float mass){fInvMassPtXi->Fill(Pt,mass);};
  void FillInvMassPtv0(float Pt,float mass){if(!fMinimalBooking)fInvMassPtv0->Fill(Pt,mass);};
  void FillInvMassLambda(int iBin,float mass){if(!fMinimalBooking)fInvMassv0[iBin]->Fill(mass);};
  void FillXiPt(int iBin,float pT){if(!fMinimalBooking)fXiPt[iBin]->Fill(pT);};
  void FillMomRapXi(int iBin,float rap,float pT){if(!fMinimalBooking)fP_Y_Xi[iBin]->Fill(rap,pT);};
  void FillDCAXiDaug(int iBin,float DCA){if(!fMinimalBooking)fDCAXiDaug[iBin]->Fill(DCA);};
  void FillMinDistPrimVtxBach(int iBin,float dist){if(!fMinimalBooking)fMinDistVtxBach[iBin]->Fill(dist);};
  void FillCPAXi(int iBin,float cpa){if(!fMinimalBooking)fCPAXi[iBin]->Fill(cpa);};
  void FillTransverseRadiusXi(int iBin,float rad){if(!fMinimalBooking)fTransRadiusXi[iBin]->Fill(rad);};
  void FillMaxDCAv0Daug(int iBin,float dca){if(!fMinimalBooking)fv0MaxDCADaug[iBin]->Fill(dca);};
  void FillCPAv0(int iBin,float cpa){if(!fMinimalBooking)fCPAv0[iBin]->Fill(cpa);};
  void Fillv0Pt(int iBin,float v0pT){if(!fMinimalBooking)fv0Pt[iBin]->Fill(v0pT);};
  void FillTransverseRadiusv0(int iBin,float rad){if(!fMinimalBooking)fTransRadiusv0[iBin]->Fill(rad);};
  void FillMinDistPrimVtxv0(int iBin,float dist){if(!fMinimalBooking)fMinDistVtxv0[iBin]->Fill(dist);};
  void FillMinDistPrimVtxv0DaugPos(int iBin,float dist){if(!fMinimalBooking)fMinDistVtxv0DaugPos[iBin]->Fill(dist);};
  void FillMinDistPrimVtxv0DaugNeg(int iBin,float dist){if(!fMinimalBooking)fMinDistVtxv0DaugNeg[iBin]->Fill(dist);};
  void FillPodolandski(int iBin,float alpha,float qt){if(!fMinimalBooking)fPodolandski[iBin]->Fill(alpha,qt);};
  void FillInvMassPerRunNumber(int RunNumber,float mass){
    if(!fMinimalBooking)fInvMassPerRunNumber->Fill(RunNumber,mass);};
  TList *GetHistList() {return fHistList;};
  private:
  bool fMinimalBooking;
  TList *fHistList;
  TH1F *fCutCounter;
  TProfile *fConfig;
  TList *fCascadeQA[2];
  TH2F *fInvMassPtXi;
  TH2F *fInvMassPtv0;
  TH1F *fInvMass[2];
  TH1F *fInvMassv0[2];
  TH1F *fXiPt[2];
  TH2F *fP_Y_Xi[2];
  TH1F *fDCAXiDaug[2];
  TH1F *fMinDistVtxBach[2];
  TH1F *fCPAXi[2];
  TH1F *fTransRadiusXi[2];
  TH1F *fv0MaxDCADaug[2];
  TH1F *fCPAv0[2];
  TH1F *fv0Pt[2];
  TH1F *fTransRadiusv0[2];
  TH1F *fMinDistVtxv0[2];
  TH1F *fMinDistVtxv0DaugPos[2];
  TH1F *fMinDistVtxv0DaugNeg[2];
  TH2F *fPodolandski[2];
  TH2F *fInvMassPerRunNumber;
  ClassDef(AliFemtoDreamCascadeHist,2)
};

#endif /* ALIFEMTODREAMCASCADEHIST_H_ */
