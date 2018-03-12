#ifndef ALIJJTTASK_H
#define ALIJJTTASK_H

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
#include "AliJJtAna.h"

class AliJCatalystTask;
//==============================================================


class AliJJtTask : public AliAnalysisTaskSE {

 public:
  AliJJtTask(); // Default constructor
  AliJJtTask(const char *name,  TString inputformat); // Constructor
  AliJJtTask(const AliJJtTask& ap); // Copy constructor
  AliJJtTask& operator = (const AliJJtTask& ap); // Equal sign operator
  virtual ~AliJJtTask(); // Destructor

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); // Output object creation
  virtual void Init();  // Initialization
  virtual void LocalInit() { Init(); } // Initialization
  virtual void UserExec(Option_t *option); // Functionality of the task
  virtual void Terminate(Option_t *); // Closing formalities

  // Methods specific for this class
  void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name
  void SetJtAnalysis(AliJJtAna *ana){ fJtAna=ana; } // Setter for analysis

 private:
  
  AliJCatalystTask *fJCatalystTask;  // 
  TString           fJCatalystTaskName; // Name for JCatalyst task
  AliJJtAna *fJtAna; // Analysis for the data
  TDirectory     *fOutput; // Output directory

  ClassDef(AliJJtTask, 1); 
};
#endif // AliJJtTask_H
