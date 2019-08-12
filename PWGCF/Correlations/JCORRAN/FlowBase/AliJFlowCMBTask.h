#ifndef ALIJFLOWCMBTASK_H
#define ALIJFLOWCMBTASK_H

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
#include <AliAnalysisTaskFlowVectorCorrections.h>
#include "AliJFlowHistos.h"


using namespace std;
//==============================================================
class TClonesArray;
class AliJFlowHistos;

class AliJFlowCMBTask : public AliAnalysisTaskSE {

 public:
  AliJFlowCMBTask();
  AliJFlowCMBTask(const char *name,  TString inputformat);
  AliJFlowCMBTask(const AliJFlowCMBTask& ap);   
  AliJFlowCMBTask& operator = (const AliJFlowCMBTask& ap);
  virtual ~AliJFlowCMBTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t* );
  void CalculateEventPlane(TClonesArray *inList);
  TComplex GetQvectorsEP(int id, int ih) { return QvectorsEP[id][ih];}
  double GetEventPlaneALICE(int id, int ih) { return fEventPlaneALICE[id][ih];}
  AliJCatalystTask *GetJCatalystTask() {return fJCatalystTask;}
  Bool_t  IsMC()const{ return fIsMC; }
  void    SetIsMC(Bool_t b){ fIsMC=b; }

  // Methods specific for this class
  void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name
  enum DETECTOR{
	  D_TPC, //TPC full coverage
	  D_TPC_ETAA, //TPC with eta gap
	  D_TPC_ETAC,
	  D_V0A,
	  D_V0C,
	  D_V0P, // V0+
	  D_VIRT, //FWD virtual
	  D_COUNT
  };

 private:

  AliJCatalystTask *fJCatalystTask;  //
  TString           fJCatalystTaskName; // Name for JCatalyst task
  AliAnalysisTaskFlowVectorCorrections *fFlowVectorTask; //
  TComplex QvectorsEP[D_COUNT][2]; // 0->2th 1->3rd
  double  fEventPlaneALICE[D_COUNT][2]; // 0->2th 1->3rd
  Bool_t      fIsMC;       // MC data or real data
  AliJFlowHistos *fhistos;
  int fCBin;
  TDirectory     *fOutput; // Output directory

  ClassDef(AliJFlowCMBTask, 1); 
};
#endif // AliJFlowCMBTask_H
