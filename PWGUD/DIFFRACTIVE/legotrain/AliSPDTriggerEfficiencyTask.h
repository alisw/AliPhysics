/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */

#ifndef AliSPDTriggerEfficiencyTask_h
#define AliSPDTriggerEfficiencyTask_h 1

// Task to create a tree for SPD trigger efficiency studies
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskSE.h"
#include "TBits.h"
#include "TArrayI.h"
#include "TArrayF.h"
class TList;
class TTree;
class TH1D;
class TH2D;
class TObjString;
class TClonesArray;
class AliESDtrackCuts;

class AliSPDTriggerEfficiencyTask : public AliAnalysisTaskSE {
 public: 
  AliSPDTriggerEfficiencyTask(const char* name = "AliSPDTriggerEfficiencyTask");
  virtual ~AliSPDTriggerEfficiencyTask(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
 protected:
  AliSPDTriggerEfficiencyTask(const  AliSPDTriggerEfficiencyTask &task);
  AliSPDTriggerEfficiencyTask& operator=(const  AliSPDTriggerEfficiencyTask &task);
  AliESDtrackCuts* fTrackCutsBit0;  //!
  AliESDtrackCuts* fTrackCutsBit5;  //!
  TList* fListOfHistos;             //! list of output histograms
  TTree* fTree;                     //! analysis tree
  TH2D* fTriggersVsRun;             //!
  TClonesArray* fTracks;            //!
  TObjString fClassesFired;
  Int_t fRunNumber; 
  UInt_t fPeriod;
  UInt_t fOrbit;
  UShort_t fBC;
  UInt_t fL0inputs;
  Float_t fVtxX;
  Float_t fVtxY;
  Float_t fVtxZ;
  Bool_t fVtxTPC;
  Int_t fVtxContributors;
  Int_t fNofTracklets;
  Int_t fNofITSClusters[6];
  TBits fIR1;
  TBits fIR2;
  TBits fFOmap;
  TBits fFiredChipMap;

  
  ClassDef(AliSPDTriggerEfficiencyTask,1)
};

#endif
