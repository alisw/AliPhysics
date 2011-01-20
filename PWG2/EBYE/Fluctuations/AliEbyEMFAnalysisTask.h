#ifndef ALIEBYEMFANALYSISTASK_H
#define ALIEBYEMFANALYSISTASK_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*-------------------------------------------------------------------------
 *                     AliEbyEMFAnalysisTask Class  
 *       This class deals with runing  Multiplicity Fluctuation Task 
 *                 origin: Satyajit Jena <sjena@cern.ch>
 * 
 *------------------------------------------------------------------------*/

class TString;
class TH1F;
class TH2F;

class TList;

#include "AliAnalysisTaskSE.h"

//class AliEbyEEventBase;
class AliEbyEMultiplicityFluctuationAnalysis;

class  AliEbyEMFAnalysisTask: public AliAnalysisTaskSE {
 public:
  AliEbyEMFAnalysisTask(const char *name = "AliEbyEMFAnalysisTask");
  virtual ~AliEbyEMFAnalysisTask() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
  void SetAnalysisObject(AliEbyEMultiplicityFluctuationAnalysis *const analysis) { fEbyEMFBase = analysis;}

 

 private:

  TList *fListPhy;

  AliEbyEMultiplicityFluctuationAnalysis *fEbyEMFBase;
  TH1F *fEvtCounter;
  Int_t fECnt;
 
  AliEbyEMFAnalysisTask(const AliEbyEMFAnalysisTask&); 
  AliEbyEMFAnalysisTask& operator=(const AliEbyEMFAnalysisTask&); 
  
  ClassDef(AliEbyEMFAnalysisTask, 1); 
};

#endif

