#ifndef ALIJHSCTASK_H
#define ALIJHSCTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for providing various flow informations 
// author: D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <TDirectory.h>
#include <TComplex.h>
#include <AliLog.h>
#include <AliAnalysisTaskSE.h>
#include <AliJCatalystTask.h>
#include "AliAnalysisAnaTwoMultiCorrelations.h" // Cindy's ana code


using namespace std;
//==============================================================
class TClonesArray;
class AliJFlowHistos;

class AliJHSCTask : public AliAnalysisTaskSE {

 public:
  AliJHSCTask();
  AliJHSCTask(const char *name);
  AliJHSCTask(const AliJHSCTask& ap);   
  AliJHSCTask& operator = (const AliJHSCTask& ap);
  virtual ~AliJHSCTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t* );
  AliJCatalystTask *GetJCatalystTask() {return fJCatalystTask;}
  void BookHistos(TClonesArray *inList);
  Bool_t  IsMC()const{ return fIsMC; }
  void    SetIsMC(Bool_t b){ fIsMC=b; }

  // Methods specific for this class
  void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name
  TString GetJCatalystTaskName(){ return fJCatalystTaskName; } // Setter for filter task name

 private:

  AliJCatalystTask *fJCatalystTask;  //
  TString           fJCatalystTaskName; // Name for JCatalyst task
  Bool_t fFirstEvent; //
  Bool_t      fIsMC;       // MC data or real data
  AliAnalysisAnaTwoMultiCorrelations *fTwoMultiAna;
  

  ClassDef(AliJHSCTask, 1); 
};
#endif // AliJHSCTask_H
