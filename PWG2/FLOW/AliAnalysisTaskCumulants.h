/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */


#ifndef AliAnalysisTaskCumulants_H
#define AliAnalysisTaskCumulants_H

// AliAnalysisTaskCumulants:
// analysis task for 
// Cumulant method
// with many authors (N.K. R.S. A.B.)
// who do something

class AliESDEvent;
class AliAODEvent;
class AliCFManager;
class AliFlowAnalysisWithCumulants;
class AliFlowEventSimpleMaker;
class TFile;

#include "AliAnalysisTask.h"

class AliAnalysisTaskCumulants : public AliAnalysisTask {
 public:
  //AliAnalysisTaskCumulants(const char *name = "AliAnalysisTaskCumulants");
  AliAnalysisTaskCumulants();
  AliAnalysisTaskCumulants(const char *name);
  virtual ~AliAnalysisTaskCumulants() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void           SetAnalysisType(TString type) {this->fAnalysisType = type;}
  TString GetAnalysisType() const    { return this->fAnalysisType; }

  void SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr; } 
  AliCFManager* GetCFManager1()           {return this->fCFManager1; }
  void SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr; } 
  AliCFManager* GetCFManager2()           {return this->fCFManager2; } 
 
 private:
 
  AliAnalysisTaskCumulants(const AliAnalysisTaskCumulants& aAnalysis);
  AliAnalysisTaskCumulants& operator=(const AliAnalysisTaskCumulants& aAnalysis);

  AliESDEvent *fESD;                      //ESD object
  AliAODEvent* fAOD;                      //AOD object
  AliFlowAnalysisWithCumulants* fMyCumuAnalysis;  //Cumulant analysis object
  AliFlowEventSimpleMaker* fEventMaker;   //FlowEventSimple maker object
  TString fAnalysisType;                  //string to select which kind of input to analyse: ESD, AOD or MC
  AliCFManager* fCFManager1;              // correction framework manager
  AliCFManager* fCFManager2;              // correction framework manager
  TList  *fListHistos;                    //collection of output //NEW
   
  ClassDef(AliAnalysisTaskCumulants, 1); // example of analysis
};

#endif











