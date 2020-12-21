#ifndef ALIJCIAAEPTASK_H
#define ALIJCIAAEPTASK_H

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
#include <AliJIaaAna.h>
#include "AliJFlowBaseTask.h"

//==============================================================


class AliJCIaaEPTask : public AliAnalysisTaskSE {

 public:
  AliJCIaaEPTask(); // Default constructor
  AliJCIaaEPTask(const char *name,  TString inputformat); // Constructor
  AliJCIaaEPTask(const AliJCIaaEPTask& ap); // Copy constructor
  AliJCIaaEPTask& operator = (const AliJCIaaEPTask& ap); // Equal sign operator
  virtual ~AliJCIaaEPTask(); // Destructor

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); // Output object creation
  virtual void Init();  // Initialization
  virtual void LocalInit() { Init(); } // Initialization
  virtual void UserExec(Option_t *option); // Functionality of the task
  virtual void Terminate(Option_t *); // Closing formalities

  // Methods specific for this class
  void SetJFlowBaseTaskName(TString name){ fJFlowBaseTaskName=name; } // Setter for filter task name
  void SetEPDector(int i) { fEPDetID = i;}
  void SetAnalysis(AliJIaaAna *ana){ fIaaAna=ana; } // Setter for analysis
  void SetEPmin(double min) { fEPmin = min;}
  void SetEPmax(double max) { fEPmax = max;}


 private:

  AliJFlowBaseTask *fJFlowBaseTask;  // 
  TString           fJFlowBaseTaskName; // Name for JCatalyst task
  int fEPDetID; 
  double fEPmin;
  double fEPmax;
  AliJIaaAna *fIaaAna; // Analysis for the data
  TDirectory     *fOutput; // Output directory

  ClassDef(AliJCIaaEPTask, 1); 
};
#endif // AliJCIaaEPTask_H
