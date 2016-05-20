#ifndef ALIJDIHADRONJTTASK_H
#define ALIJDIHADRONJTTASK_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: Jussi Viinikainen
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////


#include <AliAnalysisTaskSE.h>
#include "AliJJtAnalysis.h"

class AliJCORRANTask;
//==============================================================

using namespace std;

class AliJDiHadronJtTask : public AliAnalysisTaskSE {

 public:
  AliJDiHadronJtTask(); // Default constructor
  AliJDiHadronJtTask(const char *name,  TString inputformat); // Constructor
  AliJDiHadronJtTask(const AliJDiHadronJtTask& ap); // Copy constructor
  AliJDiHadronJtTask& operator = (const AliJDiHadronJtTask& ap); // Equal sign operator
  virtual ~AliJDiHadronJtTask(); // Destructor

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); // Output object creation
  virtual void Init();  // Initialization
  virtual void LocalInit() { Init(); } // Initialization
  virtual void UserExec(Option_t *option); // Functionality of the task
  virtual void Terminate(Option_t *); // Closing formalities
  virtual Bool_t UserNotify() { std::cout<<"DEBUG UserNotify"<<std::endl; return kTRUE;} // Debug function

  // Methods specific for this class
  void SetFilterTaskName(TString name){ fFilterTaskName=name; } // Setter for filter task name
  void SetJtAnalysis(AliJJtAnalysis *ana){ fJtAnalysis=ana; } // Setter for analysis

 private:
  
  AliJCORRANTask *fFilterTask;  // Filter task formatting the data for JCORRAN format
  TString         fFilterTaskName; // Name for the filter task
  AliJJtAnalysis *fJtAnalysis; // Analysis for the data
  TDirectory     *fOutput; // Output directory

  ClassDef(AliJDiHadronJtTask, 1); 
};
#endif // AliJDiHadronJtTask_H
