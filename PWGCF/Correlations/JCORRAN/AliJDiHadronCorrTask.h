#ifndef ALIJDIHADRONCORRTASK_H
#define ALIJDIHADRONCORRTASK_H

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
#include "AliJCORRAN.h"

class AliJCORRANTask;
//==============================================================

using namespace std;

class AliJDiHadronCorrTask : public AliAnalysisTaskSE {

 public:
  AliJDiHadronCorrTask();
  AliJDiHadronCorrTask(const char *name,  TString inputformat);
  AliJDiHadronCorrTask(const AliJDiHadronCorrTask& ap);   
  AliJDiHadronCorrTask& operator = (const AliJDiHadronCorrTask& ap);
  virtual ~AliJDiHadronCorrTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual Bool_t UserNotify() { std::cout<<"DEBUG UserNotify"<<std::endl; return kTRUE;}

  void SetFilterTaskName(TString name){ fFilterTaskName=name; }
  void SetJCORRAN(AliJCORRAN * jco){ fJCORRAN=jco; }

 private:
  
  AliJCORRANTask * fFilterTask;
  TString          fFilterTaskName;
  AliJCORRAN     * fJCORRAN;
  TDirectory     * fOutput;

  ClassDef(AliJDiHadronCorrTask, 1); 
};
#endif // AliJDiHadronCorrTask_H
