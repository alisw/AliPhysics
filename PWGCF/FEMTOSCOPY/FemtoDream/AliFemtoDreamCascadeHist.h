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
  AliFemtoDreamCascadeHist(float mass);
  virtual ~AliFemtoDreamCascadeHist();
  void FillCutCounter(int bin) {fCutCounter->Fill(bin);};
  void FillConfig(int iBin,float val){fConfig->Fill(iBin,val);};
  void FillInvMassXi(int iBin,float mass){fInvMass[iBin]->Fill(mass);};
  void FillInvMassPtXi(float Pt,float mass){fInvMassPtXi->Fill(Pt,mass);};
  void FillInvMassPtv0(float Pt,float mass){fInvMassPtv0->Fill(Pt,mass);};
  void FillInvMassLambda(int iBin,float mass){fInvMassv0[iBin]->Fill(mass);};
  void FillXiPt(int iBin,float pT){fXiPt[iBin]->Fill(pT);};
  void FillMomRapXi(int iBin,float rap,float pT){fP_Y_Xi[iBin]->Fill(rap,pT);};
  void FillDCAXiDaug(int iBin,float DCA){fDCAXiDaug[iBin]->Fill(DCA);};
  void FillMinDistPrimVtxBach(int iBin,float dist){fMinDistVtxBach[iBin]->Fill(dist);};
  void FillCPAXi(int iBin,float cpa){fCPAXi[iBin]->Fill(cpa);};
  void FillTransverseRadiusXi(int iBin,float rad){fTransRadiusXi[iBin]->Fill(rad);};
  void FillMaxDCAv0Daug(int iBin,float dca){fv0MaxDCADaug[iBin]->Fill(dca);};
  void FillCPAv0(int iBin,float cpa){fCPAv0[iBin]->Fill(cpa);};
  void FillTransverseRadiusv0(int iBin,float rad){fTransRadiusv0[iBin]->Fill(rad);};
  void FillMinDistPrimVtxv0(int iBin,float dist){fMinDistVtxv0[iBin]->Fill(dist);};
  void FillMinDistPrimVtxv0DaugPos(int iBin,float dist){fMinDistVtxv0DaugPos[iBin]->Fill(dist);};
  void FillMinDistPrimVtxv0DaugNeg(int iBin,float dist){fMinDistVtxv0DaugNeg[iBin]->Fill(dist);};
  void FillPodolandski(int iBin,float alpha,float qt){fPodolandski[iBin]->Fill(alpha,qt);};
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
  TH2F *fPodolandski[2];
  ClassDef(AliFemtoDreamCascadeHist,1)
};

#endif /* ALIFEMTODREAMCASCADEHIST_H_ */
