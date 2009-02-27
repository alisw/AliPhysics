/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */


#ifndef AliAnalysisTaskLeeYangZeros_H
#define AliAnalysisTaskLeeYangZeros_H

// AliAnalysisTaskLeeYangZeros:
// analysis task for 
// Lee Yang Zeroes method
// Author: 
// Naomi van der Kolk (kolk@nikhef.nl)             

class AliFlowEventSimple;
class AliFlowAnalysisWithLeeYangZeros;
class TFile;
class TList;

#include "TString.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskLeeYangZeros : public AliAnalysisTask {
 public:
  AliAnalysisTaskLeeYangZeros();
  AliAnalysisTaskLeeYangZeros(const char *name, Bool_t firstrun);
  virtual ~AliAnalysisTaskLeeYangZeros();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  //lyz flags
  void           SetFirstRunLYZ(Bool_t kt)     { this->fFirstRunLYZ = kt ;  }
  Bool_t         GetFirstRunLYZ() const        { return this->fFirstRunLYZ ; }
  void           SetUseSumLYZ(Bool_t kt)       { this->fUseSumLYZ = kt ;  }
  Bool_t         GetUseSumLYZ() const          { return this->fUseSumLYZ ; }
 
 private:
 
  AliAnalysisTaskLeeYangZeros(const AliAnalysisTaskLeeYangZeros& aAnalysis);
  AliAnalysisTaskLeeYangZeros& operator=(const AliAnalysisTaskLeeYangZeros& aAnalysis);

  AliFlowEventSimple* fEvent;             // input event
  AliFlowAnalysisWithLeeYangZeros* fLyz;  // LYZ analysis object
  
  TFile*           fFirstRunFile;         // file from the first loop over events
  TList*           fListHistos;           // collection of output
  
  //flags
  Bool_t fFirstRunLYZ ;    // flag for lyz analysis 
  Bool_t fUseSumLYZ ;      // flag for lyz analysis 
        
  ClassDef(AliAnalysisTaskLeeYangZeros, 1);  //AliAnalysisTaskLeeYangZeros class object
};

#endif

