/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskMCEventPlane_H
#define AliAnalysisTaskMCEventPlane_H

// AliAnalysisTaskMCEventPlane:
// analysis task for 
// Monte Carlo Event Plane
// Author: 
//        Naomi van der Kolk (kolk@nikhef.nl)

class AliESDEvent;
class AliAODEvent;
class AliFlowAnalysisWithMCEventPlane;
class AliFlowEventSimpleMaker;

#include "TString.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskMCEventPlane : public AliAnalysisTask {
 public:
  AliAnalysisTaskMCEventPlane(const char *name = "AliAnalysisTaskMCEventPlane");
  virtual ~AliAnalysisTaskMCEventPlane() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisType(TString type) { this->fAnalysisType = type; }
  TString GetAnalysisType() const    { return this->fAnalysisType; }
  
 private:
 
  AliAnalysisTaskMCEventPlane(const AliAnalysisTaskMCEventPlane& aAnalysis);
  AliAnalysisTaskMCEventPlane& operator=(const AliAnalysisTaskMCEventPlane& aAnalysis);
  
  AliESDEvent *fESD;                      // ESD object
  AliAODEvent *fAOD;                      // AOD object
  AliFlowAnalysisWithMCEventPlane* fMc;   // MC EP analysis object
  AliFlowEventSimpleMaker* fEventMaker;   // FlowEventSimple maker object
  TString fAnalysisType;                  // can be MC, ESD or AOD

  ClassDef(AliAnalysisTaskMCEventPlane, 1); // example of analysis
};

#endif

