/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */


#ifndef AliAnalysisTaskCumulants_H
#define AliAnalysisTaskCumulants_H

// AliAnalysisTaskCumulants:
// analysis task for 
// Cumulant method
// with many authors
// who do something

class AliESDEvent;
class AliAODEvent;
class AliFlowAnalysisWithCumulants;
class AliFlowEventSimpleMaker;
class TFile;

#include "AliAnalysisTask.h"

class AliAnalysisTaskCumulants : public AliAnalysisTask {
 public:
  AliAnalysisTaskCumulants(const char *name = "AliAnalysisTaskCumulants");
  virtual ~AliAnalysisTaskCumulants() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void           SetAnalysisType(TString type) {this->fAnalysisType = type;}
 
 
 private:
 
  AliAnalysisTaskCumulants(const AliAnalysisTaskCumulants& aAnalysis);
  AliAnalysisTaskCumulants& operator=(const AliAnalysisTaskCumulants& aAnalysis);

  AliESDEvent *fESD;                      //ESD object
  AliAODEvent* fAOD;                      //AOD object
  TString fAnalysisType;                  //string to select which kind of input to analyse: ESD, AOD or MC
  AliFlowAnalysisWithCumulants* fMyCumuAnalysis;  //Cumulant analysis object
  AliFlowEventSimpleMaker* fEventMaker;   //FlowEventSimple maker object

  ClassDef(AliAnalysisTaskCumulants, 1); // example of analysis
};

#endif

