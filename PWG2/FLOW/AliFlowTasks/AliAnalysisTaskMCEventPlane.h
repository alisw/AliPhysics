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

class AliFlowEventSimple;
class AliFlowAnalysisWithMCEventPlane;
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

 private:
 
  AliAnalysisTaskMCEventPlane(const AliAnalysisTaskMCEventPlane& aAnalysis);
  AliAnalysisTaskMCEventPlane& operator=(const AliAnalysisTaskMCEventPlane& aAnalysis);
  
  AliFlowEventSimple*               fEvent;        //input event
  AliFlowAnalysisWithMCEventPlane*  fMc;           // MC EP analysis object
  TList*                            fListHistos;   // collection of output

  ClassDef(AliAnalysisTaskMCEventPlane, 0); // AliAnalysisTaskMCEventPlane class object
};

#endif

