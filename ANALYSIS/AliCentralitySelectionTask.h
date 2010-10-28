#ifndef ALICENTRALITYSELECTIONTASK_H
#define ALICENTRALITYSELECTIONTASK_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliCentralitySelectionTask
//   author: Alberica Toia
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TFile;
class TH2F;

class AliCentralitySelectionTask : public AliAnalysisTaskSE {

 public:

  AliCentralitySelectionTask();
  AliCentralitySelectionTask(const char *name);
  AliCentralitySelectionTask& operator= (const AliCentralitySelectionTask& ana);
  AliCentralitySelectionTask(const AliCentralitySelectionTask& c);
  virtual ~AliCentralitySelectionTask();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInput(const char* input)         {fAnalysisInput = input;}
  void SetMCInput()                        {fIsMCInput = kTRUE;}
  
  void SetPercentileFile(TString filename) {fCentfilename = filename;}
  void SetCentralityMethod(const char* x);

 private:
  Int_t    fDebug;	   	// Debug flag
  TString  fAnalysisInput; 	// "ESD", "AOD"
  Bool_t   fIsMCInput;          // true when input is MC
  TFile   *fFile;               // file that holds the centrality vs multiplicity
  TString  fCentfilename;       // name of this file
  TString  fMethod;             // method to select centrality
  Float_t  fCent;               // percentile centrality
  TH1D    *fHtemp;              // histogram with centrality vs multiplicity
  
  
  ClassDef(AliCentralitySelectionTask,1); 

};

#endif

