#ifndef ALIANALYSISTASKVNV0_H
#define ALIANALYSISTASKVNV0_H

// ROOT includes
#include <TObject.h>
#include <TClonesArray.h>
#include "TTree.h"
#include <TList.h>
#include <TProfile.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include <AliCFContainer.h>
#include "AliFlowBayesianPID.h"
#include "AliFlowVZEROResults.h"

class TH2F;

class AliAnalysisTaskVnV0 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskVnV0();
  AliAnalysisTaskVnV0(const char *name);

  virtual ~AliAnalysisTaskVnV0();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *); 
  
  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     
  Double_t GetMinPt() { return fMinPt; }   

  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetMinPt(Double_t value) {fMinPt = value;}   
  virtual void SetV2(Bool_t val){fV2 = val;};
  virtual void SetV3(Bool_t val){fV3 = val;};

  virtual void SetMC(Bool_t flag = kTRUE){fIsMC = flag;};
  virtual void SetQA(Bool_t flag = kTRUE){fQAsw = flag;};

  void OpenInfoCalbration(Int_t run);

 private:
  AliAnalysisTaskVnV0(const AliAnalysisTaskVnV0 &old);
  AliAnalysisTaskVnV0& operator=(const AliAnalysisTaskVnV0 &source);

  virtual Float_t GetVertex(AliAODEvent* aod) const;
  virtual void Analyze(AliAODEvent* aodEvent, Float_t v0Centr); 

  AliAODEvent* fAOD;                //! AOD object

  static const Int_t nCentrBin = 9;          // # cenrality bins

  //
  // Cuts and options
  //
  Double_t     fVtxCut;             // Vtx cut on z position in cm
  Double_t     fEtaCut;             // Eta cut used to select particles
  Double_t     fMinPt;              // Min pt - for histogram limits

  Int_t fRun;                       // current run checked to load VZERO calibrations

  TList *fList,*fList2,*fList3,*fList4;             // List for output objects
  //
  // Output objects
  TProfile *fMultV0;                // object containing VZERO calibration information
  Float_t fV0Cpol,fV0Apol;          // loaded by OADB
  Float_t fMeanQ[nCentrBin][2][2];           // and recentering
  Float_t fWidthQ[nCentrBin][2][2];          // ...
  Float_t fMeanQv3[nCentrBin][2][2];         // also for v3
  Float_t fWidthQv3[nCentrBin][2][2];        // ...

  TProfile *fHResTPCv0A2,*fHResTPCv0C2,*fHResv0Cv0A2;   // TProfile for subevent resolution (output)
  TProfile *fHResTPCv0A3,*fHResTPCv0C3,*fHResv0Cv0A3;   // also for v3

  TH2F *fPhiRPv0A,*fPhiRPv0C;          // EP distribution vs. centrality (v2)
  TH2F *fPhiRPv0Av3,*fPhiRPv0Cv3;      // EP distribution vs. centrality (v3)

  AliCFContainer *fPhiTracks;          // phi distribution of particles (if needed)


  AliCFContainer *fQA,*fQA2;           // QA container (v2)
  AliCFContainer *fQAv3,*fQA2v3;       // QA container (v3)

  AliFlowBayesianPID *fPID;            // PID class for the Bayesian probabilities
 
  TTree *fTree;                        // tree to debug EP (if needed)

  Float_t fCentrality;  // current centrality for the tree
  Float_t evPlAngV0ACor2,evPlAngV0CCor2,evPlAng2;   // subevent EPs (v2)
  Float_t evPlAngV0ACor3,evPlAngV0CCor3,evPlAng3;   // subevent EPs (v3)

  Bool_t fV2,fV3; // swith to set the armonics

  AliFlowVZEROResults *fContAllChargesV0A,*fContAllChargesV0C,*fContAllChargesV0Av3,*fContAllChargesV0Cv3,*fContAllChargesMC,*fContAllChargesMCv3;

  Bool_t fIsMC; // if MC
  Bool_t fQAsw;   // if QA

  ClassDef(AliAnalysisTaskVnV0, 2);    //Analysis task v2 and v3 analysis on AOD
};

#endif
