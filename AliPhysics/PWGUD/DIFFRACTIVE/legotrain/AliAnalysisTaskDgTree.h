/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */

#ifndef AliAnalysisTaskDgTree_h
#define AliAnalysisTaskDgTree_h 1

// Task to create double gap tree
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskSE.h"
#include "TBits.h"
#include "TArrayI.h"
#include "AliPIDCombined.h"
#include "AliESDtrackCuts.h"
#include "AliTimeRangeCut.h"

class TList;
class TTree;
class TH1D;
class TH2D;
class TObjString;
class TClonesArray;


class AliAnalysisTaskDgTree : public AliAnalysisTaskSE {
 public: 
  AliAnalysisTaskDgTree(const char* name = "AliAnalysisTaskDgTree");
  virtual ~AliAnalysisTaskDgTree(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
 protected:
  AliAnalysisTaskDgTree(const  AliAnalysisTaskDgTree &task);
  AliAnalysisTaskDgTree& operator=(const  AliAnalysisTaskDgTree &task);
  AliTimeRangeCut fTimeRangeCut;
  AliPIDCombined* fPIDCombined;     //
  AliESDtrackCuts* fTrackCutsBit0;  //!
  AliESDtrackCuts* fTrackCutsBit1;  //!
  AliESDtrackCuts* fTrackCutsBit5;  //!
  TList* fListOfHistos;             //! list of output histograms
  TH2D* fTriggersVsRun;             //!
  TH1D* fEventStatistics;           //!
  TTree* fTree;                     //! analysis tree
  TObjString* fChunkFileName;       //
  TClonesArray* fTracks;            //!
  TClonesArray* fTracklets;         //!
  TClonesArray* fSATracks;          //!
  TClonesArray* fMcParticles;       //!
  Float_t fMcEventWeight;
  TObjString fClassesFired;
  Int_t fEventInFile;
  UInt_t fPeriod;
  UInt_t fOrbit;
  UShort_t fBC;
  UInt_t fL0inputs;
  UInt_t fL1inputs;
  Int_t fRunNumber; 
  Char_t fV0ADecision;
  Char_t fV0CDecision;
  Char_t fADADecision;
  Char_t fADCDecision;
  Float_t fVtxX;
  Float_t fVtxY;
  Float_t fVtxZ;
  Bool_t fVtxTPC;
  Int_t fVtxContributors;
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
  Int_t fNofTracklets;
  Int_t fNofITSClusters[6];
  TBits fFOmap;
  TBits fFiredChipMap;
  UInt_t fTriggerMask[72];
  TArrayI fTOFhits;
  TArrayF fTOFhitTimes;
  TArrayI fTrackIndices;
  Float_t fZNAtower0;
  Float_t fZNCtower0;
  Float_t fZNATDC[4];
  Float_t fZNCTDC[4];

  
  ClassDef(AliAnalysisTaskDgTree,4)
};

#endif
