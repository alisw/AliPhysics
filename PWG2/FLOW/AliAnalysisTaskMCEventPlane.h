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
class AliCFManager;
class AliFlowAnalysisWithMCEventPlane;
class AliFlowEventSimpleMaker;
class TList;

#include "TString.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskMCEventPlane : public AliAnalysisTask {
 public:

  AliAnalysisTaskMCEventPlane();
  AliAnalysisTaskMCEventPlane(const char *name);
  virtual ~AliAnalysisTaskMCEventPlane();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisType(TString type) { this->fAnalysisType = type; }
  TString GetAnalysisType() const    { return this->fAnalysisType; }
  
  void           SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr; } 
  AliCFManager*  GetCFManager1()               {return this->fCFManager1; }
  void           SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr; } 
  AliCFManager*  GetCFManager2()               {return this->fCFManager2; }

 private:
 
  AliAnalysisTaskMCEventPlane(const AliAnalysisTaskMCEventPlane& aAnalysis);
  AliAnalysisTaskMCEventPlane& operator=(const AliAnalysisTaskMCEventPlane& aAnalysis);
  
  AliESDEvent *fESD;                      // ESD object
  AliAODEvent *fAOD;                      // AOD object
  TString fAnalysisType;                  // can be MC, ESD or AOD
  AliCFManager*    fCFManager1;           // correction framework manager
  AliCFManager*    fCFManager2;           // correction framework manager
  AliFlowAnalysisWithMCEventPlane* fMc;   // MC EP analysis object
  AliFlowEventSimpleMaker* fEventMaker;   // FlowEventSimple maker object
  
  TList*           fListHistos;           // collection of output

  ClassDef(AliAnalysisTaskMCEventPlane, 1); // example of analysis
};

#endif

