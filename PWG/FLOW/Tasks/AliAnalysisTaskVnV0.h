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
#include "AliFlowBayesianPID.h"
#include "AliFlowVZEROResults.h"
#include "AliFlowVZEROQA.h"

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
  
  Double_t     fVtxCut;             // Vtx cut on z position in cm
  Double_t     fEtaCut;             // Eta cut used to select particles
  Double_t     fMinPt;              // Min pt - for histogram limits

  Bool_t fV2; // switch to set the armonics
  Bool_t fV3; // switch to set the armonics
  Bool_t fIsMC; // if MC
  Bool_t fQAsw;   // if QA

  static const Int_t nCentrBin = 9;          //! # cenrality bins

  //
  // Cuts and options
  //

  Int_t fRun;                       //! current run checked to load VZERO calibrations

  TList *fList;              //! List for output objects
  TList *fList2;             //! List for output objects
  TList *fList3;             //! List for output objects
  TList *fList4;             //! List for output objects
  //
  // Output objects
  TProfile *fMultV0;                //! object containing VZERO calibration information
  Float_t fV0Cpol;          //! loaded by OADB
  Float_t fV0Apol;          //! loaded by OADB
  Float_t fMeanQ[nCentrBin][2][2];           //! and recentering
  Float_t fWidthQ[nCentrBin][2][2];          //! ...
  Float_t fMeanQv3[nCentrBin][2][2];         //! also for v3
  Float_t fWidthQv3[nCentrBin][2][2];        //! ...

  TProfile *fHResTPCv0A2;   //! TProfile for subevent resolution (output)
  TProfile *fHResTPCv0C2;   //! TProfile for subevent resolution (output)
  TProfile *fHResv0Cv0A2;   //! TProfile for subevent resolution (output)
  TProfile *fHResTPCv0A3;    //! also for v3
  TProfile *fHResTPCv0C3;   //! also for v3
  TProfile *fHResv0Cv0A3;   //! also for v3

  TH2F *fPhiRPv0A;          //! EP distribution vs. centrality (v2)
  TH2F *fPhiRPv0C;          //! EP distribution vs. centrality (v2)
  TH2F *fPhiRPv0Av3;      //! EP distribution vs. centrality (v3)
  TH2F *fPhiRPv0Cv3;      //! EP distribution vs. centrality (v3)

  AliFlowVZEROQA *fQA;            //! QA histos (v2)
  AliFlowVZEROQA *fQA2;           //! QA histos (v2)
  AliFlowVZEROQA *fQAv3;        //! QA histos (v3)
  AliFlowVZEROQA *fQA2v3;       //! QA histos (v3)

  AliFlowBayesianPID *fPID;            //! PID class for the Bayesian probabilities
 
  TTree *fTree;                        //! tree to debug EP (if needed)

  Float_t fCentrality;  //! current centrality for the tree
  Float_t evPlAngV0ACor2;   //! subevent EPs (v2)
  Float_t evPlAngV0CCor2;   //! subevent EPs (v2)
  Float_t evPlAng2;   //! subevent EPs (v2)
  Float_t evPlAngV0ACor3;   //! subevent EPs (v3)
  Float_t evPlAngV0CCor3;   //! subevent EPs (v3)
  Float_t evPlAng3;   //! subevent EPs (v3)

  AliFlowVZEROResults *fContAllChargesV0A; //! results
  AliFlowVZEROResults *fContAllChargesV0C; //! results
  AliFlowVZEROResults *fContAllChargesV0Av3; //! results
  AliFlowVZEROResults *fContAllChargesV0Cv3; //! results
  AliFlowVZEROResults *fContAllChargesMC; //! results
  AliFlowVZEROResults *fContAllChargesMCv3; //! results


  ClassDef(AliAnalysisTaskVnV0, 4);    //Analysis task v2 and v3 analysis on AOD
};

#endif
