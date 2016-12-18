#ifndef AliJDiHadronIaaTASK_H
#define AliJDiHadronIaaTASK_H

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
#include "AliJIaaAnalysis.h"

class AliJCORRANTask;
//==============================================================

using namespace std;

class AliJDiHadronIaaTask : public AliAnalysisTaskSE {

 public:
  AliJDiHadronIaaTask(); // Default constructor
  AliJDiHadronIaaTask(const char *name,  TString inputformat); // Constructor
  AliJDiHadronIaaTask(const AliJDiHadronIaaTask& ap); // Copy constructor
  AliJDiHadronIaaTask& operator = (const AliJDiHadronIaaTask& ap); // Equal sign operator
  virtual ~AliJDiHadronIaaTask(); // Destructor

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); // Output object creation
  virtual void Init();  // Initialization
  virtual void LocalInit() { Init(); } // Initialization
  virtual void UserExec(Option_t *option); // Functionality of the task
  virtual void Terminate(Option_t *); // Closing formalities
  virtual Bool_t UserNotify() { std::cout<<"DEBUG UserNotify"<<std::endl; return kTRUE;} // Debug function

  // Methods specific for this class
  void SetFilterTaskName(TString name){ fFilterTaskName=name; } // Setter for filter task name
  void SetIaaAnalysis(AliJIaaAnalysis *ana){ fJtAnalysis=ana; } // Setter for analysis

 private:
  
  AliJCORRANTask *fFilterTask;  // Filter task formatting the data for JCORRAN format
  TString         fFilterTaskName; // Name for the filter task
  AliJIaaAnalysis *fJtAnalysis; // Analysis for the data
  TDirectory     *fOutput; // Output directory

  ClassDef(AliJDiHadronIaaTask, 1);
};
#endif // AliJDiHadronIaaTask_H
