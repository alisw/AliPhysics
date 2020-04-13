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

class AliAnalysisTaskUpcTree : public AliAnalysisTaskSE {
 public: 
  AliAnalysisTaskUpcTree(const char* name = "AliAnalysisTaskUpcTree");
  virtual ~AliAnalysisTaskUpcTree(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
 protected:
  AliAnalysisTaskUpcTree(const  AliAnalysisTaskUpcTree &task);
  AliAnalysisTaskUpcTree& operator=(const  AliAnalysisTaskUpcTree &task);
  AliMuonTrackCuts* fMuonTrackCuts; //
  TList* fListOfHistos;             //! list of output histograms
  TH2D* fTriggersVsRun;
  TTree* fTree;                     //! analysis tree
  TObjString* fChunkFileName;       //
  TClonesArray* fMcParticles;       //!
  TClonesArray* fMuons;             //!
  TObjString fClassesFired;
  Int_t fEventInFile;
  UInt_t fPeriod;
  UInt_t fOrbit;
  UShort_t fBC;
  UInt_t fL0inputs;
  UInt_t fL1inputs;
  Int_t fRunNumber; 
  Int_t fNofTracklets;
  Bool_t fBBFlag[64];
  Bool_t fBGFlag[64];
  Char_t fV0ADecision;
  Char_t fV0CDecision;
  Char_t fADADecision;
  Char_t fADCDecision;
  Float_t fZNAtower0;
  Float_t fZNCtower0;
  Float_t fZNATDC[4];
  Float_t fZNCTDC[4];
  Float_t fVtxX;
  Float_t fVtxY;
  Float_t fVtxZ;
  Bool_t fVtxTPC;
  TBits fIR1;
  TBits fIR2;
  TBits fFOmap;
  TBits fFiredChipMap;
  
  ClassDef(AliAnalysisTaskUpcTree,3)
};

#endif
