#ifndef ALIEBYECFANALYSISTASK_H
#define ALIEBYECFANALYSISTASK_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*-------------------------------------------------------------------------
 *                     AliEbyECFAnalysisTask Class  
 *       This class deals with runing  Charge Fluctuation Task 
 *                 origin: Satyajit Jena <sjena@cern.ch>
 * 
 *------------------------------------------------------------------------*/

class TString;
class TH1F;
class TH2F;

class TList;

#include "AliAnalysisTaskSE.h"

class AliEbyEEventBase;
class AliEbyEChargeFluctuationAnalysis;

class  AliEbyECFAnalysisTask: public AliAnalysisTaskSE {
 public:
  AliEbyECFAnalysisTask(const char *name = "AliEbyECFAnalysisTask");
  virtual ~AliEbyECFAnalysisTask() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
  void SetAnalysisObject(AliEbyEChargeFluctuationAnalysis *const analysis) { fEbyECFBase = analysis;}

 private:

  TList *fListPhy;

  AliEbyEChargeFluctuationAnalysis *fEbyECFBase;
  TH1F *fEvtCounter;
  Int_t fECnt;
    
  AliEbyECFAnalysisTask(const AliEbyECFAnalysisTask&); 
  AliEbyECFAnalysisTask& operator=(const AliEbyECFAnalysisTask&); 
  
  ClassDef(AliEbyECFAnalysisTask, 1); 
};

#endif

