/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKMCEVENTPLANE_H
#define ALIANALYSISTASKMCEVENTPLANE_H

// AliAnalysisTaskMCEventPlane:
// analysis task for 
// Monte Carlo Event Plane
// Author: 
//        Naomi van der Kolk (kolk@nikhef.nl)

class AliFlowEventSimple;
class AliFlowAnalysisWithMCEventPlane;
class TList;

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMCEventPlane : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskMCEventPlane();
  AliAnalysisTaskMCEventPlane(const char *name);
  virtual ~AliAnalysisTaskMCEventPlane();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

 private:
 
  AliAnalysisTaskMCEventPlane(const AliAnalysisTaskMCEventPlane& aAnalysis);
  AliAnalysisTaskMCEventPlane& operator=(const AliAnalysisTaskMCEventPlane& aAnalysis);
  
  AliFlowEventSimple*               fEvent;        //input event
  AliFlowAnalysisWithMCEventPlane*  fMc;           // MC EP analysis object
  TList*                            fListHistos;   // collection of output

  ClassDef(AliAnalysisTaskMCEventPlane, 1); // AliAnalysisTaskMCEventPlane class object
};

#endif

