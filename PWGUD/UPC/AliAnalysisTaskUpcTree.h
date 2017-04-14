/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */

#ifndef AliAnalysisTaskUpcTree_h
#define AliAnalysisTaskUpcTree_h 1

// Task to create upc tree
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskSE.h"
#include "TBits.h"

class TList;
class TTree;
class TH1D;
class TH2D;
class TObjString;
class TClonesArray;
class AliMuonTrackCuts;
class AliPIDCombined;

class AliAnalysisTaskUpcTree : public AliAnalysisTaskSE {
 public: 
  AliAnalysisTaskUpcTree(const char* name = "AliAnalysisTaskUpcTree");
  virtual ~AliAnalysisTaskUpcTree(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
  void SetMC(Bool_t mc)     { fIsMC  = mc;  }
  void SetMUP(Bool_t mup=1) { fIsMUP = mup; }
 protected:
  AliAnalysisTaskUpcTree(const  AliAnalysisTaskUpcTree &task);
  AliAnalysisTaskUpcTree& operator=(const  AliAnalysisTaskUpcTree &task);
  Bool_t fIsMUP;
  Bool_t fIsMC;
  Bool_t fIsAOD;
  AliMuonTrackCuts* fMuonTrackCuts; //
  AliPIDCombined* fPIDCombined; //
  TList* fListOfHistos;             //! list of output histograms
  TH2D* fTriggersVsRun;
  TH1D* fEventStatistics;
  TTree* fTree;                     //! analysis tree
  TObjString* fChunkFileName;       //
  TClonesArray* fMcParticles;       //!
  TClonesArray* fMuons;             //!
  TClonesArray* fTracks;            //!
  TClonesArray* fTracklets;         //!
  TClonesArray* fSATracks;         //!
  TObjString fClassesFired;
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
  UInt_t fBBAFlags;
  UInt_t fBBCFlags;
  UInt_t fBGAFlags;
  UInt_t fBGCFlags;
  Bool_t fBBTriggerV0A[32];
  Bool_t fBGTriggerV0A[32];
  Bool_t fBBTriggerV0C[32];
  Bool_t fBGTriggerV0C[32];
  Char_t fV0ADecision;
  Char_t fV0CDecision;
  Float_t fMTotV0A;
  Float_t fMTotV0C;
  Float_t fADATime;
  Float_t fADCTime;
  Char_t fADADecision;
  Char_t fADCDecision;
  Float_t fMTotADA;
  Float_t fMTotADC;
  UShort_t fTriggerChargeADA;
  UShort_t fTriggerChargeADC;
  Float_t fZNATDC[4];
  Float_t fZNCTDC[4];
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
  Int_t fNofITSClusters[6];
  Bool_t fIRproblem;
  TBits fIR1;
  TBits fIR2;
  TBits fVBA;
  TBits fVBC;
  TBits fUBA;
  TBits fUBC;
  TBits fSTG;
  TBits fOM2;
  TBits fOMU;
  TBits fFOmap;
  TBits fFiredChipMap;
  UInt_t fTriggerMask[72];
  Bool_t fPFBBFlagV0[64][21];
  Bool_t fPFBGFlagV0[64][21];
  Bool_t fPFBBFlagAD[16][21];
  Bool_t fPFBGFlagAD[16][21];

  ClassDef(AliAnalysisTaskUpcTree,3)
};

#endif
