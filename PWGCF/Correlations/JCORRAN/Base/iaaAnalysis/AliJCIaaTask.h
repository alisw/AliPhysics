#ifndef ALIJCIAATASK_H
#define ALIJCIAATASK_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: D.J~ Kim
// ALICE Group University of Jyvaskyla 
// Finland 
// Fill the analysis containers for AOD
//////////////////////////////////////////////////////////////////////////////

#include <AliAnalysisTaskSE.h>
#include "AliJIaaAna.h"

class AliJCatalystTask;
//==============================================================


class AliJCIaaTask : public AliAnalysisTaskSE {

 public:
  AliJCIaaTask(); // Default constructor
  AliJCIaaTask(const char *name,  TString inputformat); // Constructor
  AliJCIaaTask(const AliJCIaaTask& ap); // Copy constructor
  AliJCIaaTask& operator = (const AliJCIaaTask& ap); // Equal sign operator
  virtual ~AliJCIaaTask(); // Destructor

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); // Output object creation
  virtual void Init();  // Initialization
  virtual void LocalInit() { Init(); } // Initialization
  virtual void UserExec(Option_t *option); // Functionality of the task
  virtual void Terminate(Option_t *); // Closing formalities

  // Methods specific for this class
  void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name
  void SetAnalysis(AliJIaaAna *ana){ fIaaAna=ana; } // Setter for analysis

 private:
  
  AliJCatalystTask *fJCatalystTask;  // 
  TString           fJCatalystTaskName; // Name for JCatalyst task
  AliJIaaAna *fIaaAna; // Analysis for the data
  TDirectory     *fOutput; // Output directory

  ClassDef(AliJCIaaTask, 1); 
};
#endif // AliJCIaaTask_H
