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
      int MassNBins,float MassMin,float MassMax,bool CPAPlots);
  virtual ~AliFemtoDreamv0Hist();
  void FillConfig(int iBin,float val){fConfig->Fill(iBin,val);};
  void FillTrackCounter(int iBin){fCutCounter->Fill(iBin);};
  void FillOnFlyStatus(int i,float val){fOnFly[i]->Fill(val);};
  void FillpTCut(int i,float pT){fpTDist[i]->Fill(pT);};
  void FillPhi(int i,float phi){fPhiDist[i]->Fill(phi);};
  void FillEtaCut(int i,float eta){fetaDist[i]->Fill(eta);};
  void Fillv0DecayVtxXCut(int i,float vtx){fDecayVtxv0X[i]->Fill(vtx);};
  void Fillv0DecayVtxYCut(int i,float vtx){fDecayVtxv0Y[i]->Fill(vtx);};
  void Fillv0DecayVtxZCut(int i,float vtx){fDecayVtxv0Z[i]->Fill(vtx);};
  void FillTransverRadiusCut(int i,float rad){fTransRadius[i]->Fill(rad);};
  void FillDCAPosDaugToPrimVtxCut(int i,float dca){
    fDCAPosDaugToPrimVtx[i]->Fill(dca);
  };
  void FillDCANegDaugToPrimVtxCut(int i,float dca){
    fDCANegDaugToPrimVtx[i]->Fill(dca);
  };
  void FillDCADaugTov0VtxCut(int i, float dca){fDCADaugToVtx[i]->Fill(dca);};
  void FillCPACut(int i,float cpa){fCPA[i]->Fill(cpa);};
  void FillInvMass(int i,float mass){fInvMass[i]->Fill(mass);};
  void FillInvMassBefKaonRej(float mass){fInvMassBefKaonRej->Fill(mass);};
  void FillInvMassKaon(float mass){fInvMassKaon->Fill(mass);};
  void Fillv0MassDist(float mass){fInvMassBefSelection->Fill(mass);};
  void FillInvMassPtBins(float pT,float mass){fInvMassPt->Fill(pT,mass);};
  void FillCPAPtBins(float pT,float cpa){fCPAPtBins->Fill(pT,cpa);};
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
