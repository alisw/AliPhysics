#ifndef ALIJEFFICIENCYTASK_H
#define ALIJEFFICIENCYTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////


#include "AliAnalysisTaskSE.h"
#include "AliJEfficiencyScanner.h"

class AliJCORRANTask;
//==============================================================

using namespace std;

class AliJEfficiencyTask : public AliAnalysisTaskSE {

 public:
  AliJEfficiencyTask();
  AliJEfficiencyTask(const char *name,  TString inputformat);
  AliJEfficiencyTask(const AliJEfficiencyTask& ap);   
  AliJEfficiencyTask& operator = (const AliJEfficiencyTask& ap);
  virtual ~AliJEfficiencyTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t * opt = "");

  void SetFilterTaskName(TString name){ fFilterTaskName=name; }
  void SetJEfficiencyScanner(AliJEfficiencyScanner * jco){ fEfficiencyScanner=jco; }
  
  AliJEfficiencyScanner * GetEfficiencyScanner(){ return fEfficiencyScanner; }
  
 private:

  AliJCORRANTask * fFilterTask;
  TString          fFilterTaskName;
  AliJEfficiencyScanner * fEfficiencyScanner;
  TDirectory * fEffHistDir;


  ClassDef(AliJEfficiencyTask, 1); 
};
#endif // AliJEfficiencyTask_H
