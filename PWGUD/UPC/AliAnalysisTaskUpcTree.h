/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */

#ifndef AliAnalysisTaskUpcTree_h
#define AliAnalysisTaskUpcTree_h 1

// Task to create upc tree
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskSE.h"
#include "TBits.h"
#define NTRIGGERS 4

class TList;
class TTree;
class TH1D;
class TH2D;
class TObjString;
class TClonesArray;
class AliMuonTrackCuts;

class AliAnalysisTaskUpcTree : public AliAnalysisTaskSE {
 public: 
  AliAnalysisTaskUpcTree(const char* name = "AliAnalysisTaskUpcTree");
  virtual ~AliAnalysisTaskUpcTree(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
  void SetMC(Bool_t mc)   { fIsMC  = mc;  }

 protected:
  AliAnalysisTaskUpcTree(const  AliAnalysisTaskUpcTree &task);
  AliAnalysisTaskUpcTree& operator=(const  AliAnalysisTaskUpcTree &task);
  
  Bool_t fIsMC;
  Bool_t fIsAOD;
  AliMuonTrackCuts* fMuonTrackCuts; //
  TList* fListOfHistos;             //! list of output histograms
  TH1D*  fEventStatistics;          //!
  TH2D*  fTriggersPerRun;           //!
  TTree* fTree;                     //! analysis tree
  TObjString* fChunkFileName;       //
  Bool_t fTriggerFired[NTRIGGERS];
  Int_t fEventInFile;
  UInt_t fPeriod;
  UInt_t fOrbit;
  UShort_t fBC;
  UInt_t fL0inputs;
  UInt_t fL1inputs;
  Int_t fRunNumber; 
  Int_t fNofTracklets;
  Float_t fV0AMult[32];
  Float_t fV0CMult[32];
  Float_t fV0ATime[32];
  Float_t fV0CTime[32];
  Bool_t fBBFlag[64];
  Bool_t fBGFlag[64];
  Bool_t fBBTriggerV0A[32];
  Bool_t fBGTriggerV0A[32];
  Bool_t fBBTriggerV0C[32];
  Bool_t fBGTriggerV0C[32];
  Bool_t fBBonlineV0A;
  Bool_t fBGonlineV0A;
  Bool_t fBBonlineV0C;
  Bool_t fBGonlineV0C;
  Char_t fV0ADecision;
  Char_t fV0CDecision;
  Float_t fADATime;
  Float_t fADCTime;
  Char_t fADADecision;
  Char_t fADCDecision;
  Float_t fMTotADA;
  Float_t fMTotADC;
  UShort_t fTriggerChargeADA;
  UShort_t fTriggerChargeADC;
  Bool_t fZNAtdc;
  Bool_t fZNCtdc;
  Bool_t fZPAtdc;
  Bool_t fZPCtdc;
  Bool_t fZEM1tdc;
  Bool_t fZEM2tdc;
  Float_t fZNAenergy;
  Float_t fZNCenergy;
  Float_t fZPAenergy;
  Float_t fZPCenergy;
  Float_t fZEM1energy;
  Float_t fZEM2energy;
  Float_t fZNAtower0;
  Float_t fZNCtower0;
  Float_t fZPAtower0;
  Float_t fZPCtower0;
  Float_t fVtxX;
  Float_t fVtxY;
  Float_t fVtxZ;
  Bool_t fVtxTPC;
  UInt_t fNofITSClusters[6];
  TBits fIR1;
  TBits fIR2;
  TBits fFOmap;
  TBits fFiredChipMap;
  Int_t fNofDimuons;
  Float_t fDimuonM;
  Float_t fDimuonY;
  Float_t fDimuonPt;
  Float_t fEta1;
  Float_t fEta2;
  Float_t fPhi1;
  Float_t fPhi2;
  Float_t fPt1;
  Float_t fPt2;
  Float_t fPx1;
  Float_t fPy1;
  Float_t fPz1;
  Float_t fPx2;
  Float_t fPy2;
  Float_t fPz2;
  Char_t  fCharge1;
  Char_t  fCharge2;
  Char_t  fMatch1;
  Char_t  fMatch2;
  Float_t fRabs1;
  Float_t fRabs2;
  Bool_t fPdca1;
  Bool_t fPdca2;
  Float_t fDca1;
  Float_t fDca2;
  Float_t fChi2perNDF1;
  Float_t fChi2perNDF2;

  ClassDef(AliAnalysisTaskUpcTree,2)
};

#endif
