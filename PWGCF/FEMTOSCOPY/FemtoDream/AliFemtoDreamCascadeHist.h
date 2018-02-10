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
  AliFemtoDreamCascadeHist(double mass);
  virtual ~AliFemtoDreamCascadeHist();
  void FillCutCounter(int bin) {fCutCounter->Fill(bin);};
  void FillConfig(int iBin,double val){fConfig->Fill(iBin,val);};
  void FillInvMassXi(int iBin,double mass){fInvMass[iBin]->Fill(mass);};
  void FillInvMassPtXi(double Pt,double mass){fInvMassPtXi->Fill(Pt,mass);};
  void FillInvMassPtv0(double Pt,double mass){fInvMassPtv0->Fill(Pt,mass);};
  void FillInvMassLambda(int iBin,double mass){fInvMassv0[iBin]->Fill(mass);};
  void FillXiPt(int iBin,double pT){fXiPt[iBin]->Fill(pT);};
  void FillMomRapXi(int iBin,double rap,double pT){fP_Y_Xi[iBin]->Fill(rap,pT);};
  void FillDCAXiDaug(int iBin,double DCA){fDCAXiDaug[iBin]->Fill(DCA);};
  void FillMinDistPrimVtxBach(int iBin,double dist){fMinDistVtxBach[iBin]->Fill(dist);};
  void FillCPAXi(int iBin,double cpa){fCPAXi[iBin]->Fill(cpa);};
  void FillTransverseRadiusXi(int iBin,double rad){fTransRadiusXi[iBin]->Fill(rad);};
  void FillMaxDCAv0Daug(int iBin,double dca){fv0MaxDCADaug[iBin]->Fill(dca);};
  void FillCPAv0(int iBin,double cpa){fCPAv0[iBin]->Fill(cpa);};
  void FillTransverseRadiusv0(int iBin,double rad){fTransRadiusv0[iBin]->Fill(rad);};
  void FillMinDistPrimVtxv0(int iBin,double dist){fMinDistVtxv0[iBin]->Fill(dist);};
  void FillMinDistPrimVtxv0DaugPos(int iBin,double dist){fMinDistVtxv0DaugPos[iBin]->Fill(dist);};
  void FillMinDistPrimVtxv0DaugNeg(int iBin,double dist){fMinDistVtxv0DaugNeg[iBin]->Fill(dist);};
  TList *GetHistList() {return fHistList;};
  private:
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
  TH1F *fTransRadiusv0[2];
  TH1F *fMinDistVtxv0[2];
  TH1F *fMinDistVtxv0DaugPos[2];
  TH1F *fMinDistVtxv0DaugNeg[2];
  ClassDef(AliFemtoDreamCascadeHist,1)
};

#endif /* ALIFEMTODREAMCASCADEHIST_H_ */
