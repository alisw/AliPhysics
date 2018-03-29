#ifndef ALIJFLOWBASETASK_H
#define ALIJFLOWBASETASK_H

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
#include "AliJCatalystTask.h"

using namespace std;
//==============================================================
class TClonesArray;
//class AliJCatalystTask;

class AliJFlowBaseTask : public AliAnalysisTaskSE {

 public:
  AliJFlowBaseTask();
  AliJFlowBaseTask(const char *name,  TString inputformat);
  AliJFlowBaseTask(const AliJFlowBaseTask& ap);   
  AliJFlowBaseTask& operator = (const AliJFlowBaseTask& ap);
  virtual ~AliJFlowBaseTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t* );
  void CalculateEventPlane(TClonesArray *inList);
  TComplex GetQvectorsEP(int id, int ih) { return QvectorsEP[id][ih];}
  AliJCatalystTask *GetJCatalystTask() {return fJCatalystTask;}

  // Methods specific for this class
  void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name
  enum DETECTOR{
	  D_TPC, //TPC full coverage
	  D_TPC_ETAA, //TPC with eta gap
	  D_TPC_ETAC,
	  D_V0A,
	  D_V0C,
	  D_V0P, //V0+
	  D_COUNT
  };

 private:

  AliJCatalystTask *fJCatalystTask;  //
  TString           fJCatalystTaskName; // Name for JCatalyst task
  TComplex QvectorsEP[D_COUNT][2]; // 0->2th 1->3rd
  TDirectory     *fOutput; // Output directory

  ClassDef(AliJFlowBaseTask, 1); 
};
#endif // AliJFlowBaseTask_H
