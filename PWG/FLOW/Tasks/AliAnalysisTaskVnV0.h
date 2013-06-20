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
class AliESDtrackCuts;

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
  virtual void  SetMinDistV0(Double_t value){fMinDistV0=value;}
  virtual void  SetMaxDistV0(Double_t value){fMaxDistV0=value;}
  virtual void SetV2(Bool_t val){fV2 = val;};
  virtual void SetV3(Bool_t val){fV3 = val;};

  virtual void SetMC(Bool_t flag = kTRUE){fIsMC = flag;};
  virtual void SetQA(Bool_t flag = kTRUE){fQAsw = flag;};

  void OpenInfoCalbration(Int_t run);

  void SetFillDCAinfo(Bool_t flag=kTRUE){fFillDCA = flag;};

  void SetModulationDEDx(Bool_t flag=kTRUE){fModulationDEDx=flag;};

  void SetTPCclusterN(Int_t ncl){fNcluster=ncl;};
  static Bool_t IsPsiComputed(){return fgIsPsiComputed;};
  static Float_t GetPsi2V0A(){return fgPsi2v0a;};
  static Float_t GetPsi2V0C(){return fgPsi2v0c;};
  static Float_t GetPsi2TPC(){return fgPsi2tpc;};
  static Float_t GetPsi3V0A(){return fgPsi3v0a;};
  static Float_t GetPsi3V0C(){return fgPsi3v0c;};
  static Float_t GetPsi3TPC(){return fgPsi3tpc;};
  static Float_t GetPsi2V0AMC(){return fgPsi2v0aMC;};
  static Float_t GetPsi2V0CMC(){return fgPsi2v0cMC;};
  static Float_t GetPsi2TPCMC(){return fgPsi2tpcMC;};
  static Float_t GetPsi3V0AMC(){return fgPsi3v0aMC;};
  static Float_t GetPsi3V0CMC(){return fgPsi3v0cMC;};
  static Float_t GetPsi3TPCMC(){return fgPsi3tpcMC;};

 private:
  AliAnalysisTaskVnV0(const AliAnalysisTaskVnV0 &old); 
  AliAnalysisTaskVnV0& operator=(const AliAnalysisTaskVnV0 &source); 
 
  Int_t PassesAODCuts(AliAODv0 *myV0, AliAODEvent *tAOD,Int_t specie);

  static Bool_t fgIsPsiComputed; // flag which return if event was processed
  static Float_t fgPsi2v0a,fgPsi2v0c,fgPsi2tpc; // current Psi2
  static Float_t fgPsi2v0aMC,fgPsi2v0cMC,fgPsi2tpcMC; // current Psi2
  static Float_t fgPsi3v0aMC,fgPsi3v0cMC,fgPsi3tpcMC; // current Psi3
  static Float_t fgPsi3v0a,fgPsi3v0c,fgPsi3tpc; // current Psi3

  virtual Float_t GetVertex(AliAODEvent* aod) const;
  virtual void Analyze(AliAODEvent* aodEvent, Float_t v0Centr); 
  
  Double_t     fVtxCut;             // Vtx cut on z position in cm
  Double_t     fEtaCut;             // Eta cut used to select particles
  Double_t     fMinPt;              // Min pt - for histogram limits
  Double_t     fMinDistV0;          // Minimal distance for V0s
  Double_t     fMaxDistV0;          // Maximal distance for V0s

  Bool_t fV2; // switch to set the armonics
  Bool_t fV3; // switch to set the armonics
  Bool_t fIsMC; // if MC
  Bool_t fQAsw;   // if QA

  static const Int_t nCentrBin = 9;          //! # cenrality bins

  //
  // Cuts and options
  //

  Int_t fRun;                       //! current run checked to load VZERO calibrations

  Int_t fNcluster;           // Numer of TPC cluster required
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

  // MC
  TProfile *fHResMA2;   //! TProfile for subevent resolution (output)
  TProfile *fHResMC2;   //! TProfile for subevent resolution (output)
  TProfile *fHResAC2;   //! TProfile for subevent resolution (output)
  TProfile *fHResMA3;    //! also for v3
  TProfile *fHResMC3;   //! also for v3
  TProfile *fHResAC3;   //! also for v3

  AliFlowVZEROResults *fContAllChargesMCA; //! results
  AliFlowVZEROResults *fContAllChargesMCC; //! results
  AliFlowVZEROResults *fContAllChargesMCAv3; //! results
  AliFlowVZEROResults *fContAllChargesMCCv3; //! results

  Bool_t fFillDCA; // require to fill also DCA info
  TH2D *fHdcaPt[nCentrBin][7]; //! DCA distribution (for MC primary)
  TH2D *fHdcaPtSec[nCentrBin][7]; //! DCA distribution (for MC secondary, not used for data)
  AliFlowVZEROResults *fContQApid; //! QA pid object

  Bool_t fModulationDEDx; //add a modulation on the dE/dx response w.r.t. EP (kFALSE default)

  AliESDtrackCuts *fCutsDaughter;
  ClassDef(AliAnalysisTaskVnV0, 7);    //Analysis task v2 and v3 analysis on AOD
};

#endif
