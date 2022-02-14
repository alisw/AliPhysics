#ifndef ALIJSOUNDSTASK_H
#define ALIJSOUNDSTASK_H

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
#include "AliJFlowHistos.h"
#include "AliJEfficiency.h"


using namespace std;
//==============================================================
class TClonesArray;
class AliJFlowHistos;

class AliJSoundsTask : public AliAnalysisTaskSE {

 public:
  AliJSoundsTask();
  AliJSoundsTask(const char *name,  TString inputformat);
  AliJSoundsTask(const AliJSoundsTask& ap);   
  AliJSoundsTask& operator = (const AliJSoundsTask& ap);
  virtual ~AliJSoundsTask();

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

 private:

  AliJCatalystTask *fJCatalystTask;  //
  TString           fJCatalystTaskName; // Name for JCatalyst task
  AliJEfficiency *fEfficiency; //
  Bool_t fFirstEvent; //
  Bool_t      fIsMC;       // MC data or real data
  AliJFlowHistos *fhistos;
  int fCBin;
  TDirectory     *fOutput; // Output directory

  ClassDef(AliJSoundsTask, 1); 
};
#endif // AliJSoundsTask_H
