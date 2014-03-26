/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */

#ifndef AliAnalysisTaskUpcTree_h
#define AliAnalysisTaskUpcTree_h 1

// Task to create upc tree
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskSE.h"
#include "TBits.h"
#define NTRIGGERS 7

class TList;
class TTree;
class TH1I;
class TH2I;
class TObjString;
class TClonesArray;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliMuonTrackCuts;

class AliAnalysisTaskUpcTree : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcTree(const char* name = "AliAnalysisTaskUpcTree");
  virtual ~AliAnalysisTaskUpcTree(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
  void SetAOD(Bool_t aod) { fIsAOD = aod; }
  void SetMC(Bool_t mc)   { fIsMC  = mc; }
  void SetTrackFilter(AliAnalysisFilter* filter) { fTrackFilter = filter; }

 protected:
  AliAnalysisTaskUpcTree(const  AliAnalysisTaskUpcTree &task);
  AliAnalysisTaskUpcTree& operator=(const  AliAnalysisTaskUpcTree &task);
  
  Bool_t fIsMC;
  Bool_t fIsAOD;
  AliMuonTrackCuts* fMuonTrackCuts; //
  AliAnalysisFilter* fTrackFilter;  //
  TList* fListOfHistos;             //! list of output histograms
  TH1I*  fEventStatistics;          //!
  TH2I*  fTriggersPerRun;           //!
  TTree* fTree;                     //! analysis tree
  TClonesArray* fTPCtracks;         //!
  TClonesArray* fMUONtracks;        //!
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

  ClassDef(AliAnalysisTaskUpcTree,1)
};

#endif
