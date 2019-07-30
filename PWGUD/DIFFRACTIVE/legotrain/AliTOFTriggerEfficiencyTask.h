/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */

#ifndef AliTOFTriggerEfficiencyTask_h
#define AliTOFTriggerEfficiencyTask_h 1

// Task to create a tree for TOF trigger efficiency studies
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskSE.h"
#include "TBits.h"
#include "TArrayI.h"
#include "TArrayF.h"
class TList;
class TTree;
class TH2D;
class TObjString;

class AliTOFTriggerEfficiencyTask : public AliAnalysisTaskSE {
 public: 
  AliTOFTriggerEfficiencyTask(const char* name = "AliTOFTriggerEfficiencyTask");
  virtual ~AliTOFTriggerEfficiencyTask(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
 protected:
  AliTOFTriggerEfficiencyTask(const  AliTOFTriggerEfficiencyTask &task);
  AliTOFTriggerEfficiencyTask& operator=(const  AliTOFTriggerEfficiencyTask &task);
  TList* fListOfHistos;             //! list of output histograms
  TTree* fTree;                     //! analysis tree
  TH2D* fTriggersVsRun;             //!
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
  TBits fIR1;
  TBits fIR2;
  UInt_t fTriggerMask[72];
  TArrayI fTOFhits;
  TArrayF fTOFhitTimes;
  ClassDef(AliTOFTriggerEfficiencyTask,2)
};

#endif
