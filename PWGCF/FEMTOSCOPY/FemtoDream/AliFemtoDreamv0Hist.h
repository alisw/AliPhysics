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
  AliFemtoDreamv0Hist(
      int MassNBins,double MassMin,double MassMax,bool CPAPlots);
  virtual ~AliFemtoDreamv0Hist();
  void FillConfig(int iBin,double val){fConfig->Fill(iBin,val);};
  void FillTrackCounter(int iBin){fCutCounter->Fill(iBin);};
  void FillOnFlyStatus(int i,double val){fOnFly[i]->Fill(val);};
  void FillpTCut(int i,double pT){fpTDist[i]->Fill(pT);};
  void FillPhi(int i,double phi){fPhiDist[i]->Fill(phi);};
  void FillEtaCut(int i,double eta){fetaDist[i]->Fill(eta);};
  void Fillv0DecayVtxXCut(int i,double vtx){fDecayVtxv0X[i]->Fill(vtx);};
  void Fillv0DecayVtxYCut(int i,double vtx){fDecayVtxv0Y[i]->Fill(vtx);};
  void Fillv0DecayVtxZCut(int i,double vtx){fDecayVtxv0Z[i]->Fill(vtx);};
  void FillTransverRadiusCut(int i,double rad){fTransRadius[i]->Fill(rad);};
  void FillDCAPosDaugToPrimVtxCut(int i,double dca){
    fDCAPosDaugToPrimVtx[i]->Fill(dca);
  };
  void FillDCANegDaugToPrimVtxCut(int i,double dca){
    fDCANegDaugToPrimVtx[i]->Fill(dca);
  };
  void FillDCADaugTov0VtxCut(int i, double dca){fDCADaugToVtx[i]->Fill(dca);};
  void FillCPACut(int i,double cpa){fCPA[i]->Fill(cpa);};
  void FillInvMass(int i,double mass){fInvMass[i]->Fill(mass);};
  void FillInvMassBefKaonRej(double mass){fInvMassBefKaonRej->Fill(mass);};
  void FillInvMassKaon(double mass){fInvMassKaon->Fill(mass);};
  void Fillv0MassDist(double mass){fInvMassBefSelection->Fill(mass);};
  void FillInvMassPtBins(double pT,double mass){fInvMassPt->Fill(pT,mass);};
  void FillCPAPtBins(double pT,double cpa){fCPAPtBins->Fill(pT,cpa);};
  TList *GetHistList(){return fHistList;};
 private:
  TList *fHistList;
  TList *fv0CutQA[2];
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
  TH2F *fCPAPtBins;
  ClassDef(AliFemtoDreamv0Hist,1)
};

#endif /* ALIFEMTODREAMV0HIST_H_ */
