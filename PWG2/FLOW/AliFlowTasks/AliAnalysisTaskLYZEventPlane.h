/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskLYZEventPlane_H
#define AliAnalysisTaskLYZEventPlane_H

// AliAnalysisTaskLYZEventPlane:
// analysis task for 
// Lee Yang Zeros Event Plane
// Author: 
//        Naomi van der Kolk (kolk@nikhef.nl)

class AliFlowEventSimple;
class AliFlowLYZEventPlane;
class AliFlowAnalysisWithLYZEventPlane;
class TFile;
class TList;

#include "TString.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskLYZEventPlane : public AliAnalysisTask {
 public:
  AliAnalysisTaskLYZEventPlane();
  AliAnalysisTaskLYZEventPlane(const char *name);
  virtual ~AliAnalysisTaskLYZEventPlane();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
 
  AliAnalysisTaskLYZEventPlane(const AliAnalysisTaskLYZEventPlane& aAnalysis);
  AliAnalysisTaskLYZEventPlane& operator=(const AliAnalysisTaskLYZEventPlane& aAnalysis);

  AliFlowEventSimple*               fEvent;         // input event
  AliFlowLYZEventPlane*             fLyzEp;         //LYZ EP object
  AliFlowAnalysisWithLYZEventPlane* fLyz;           //LYZ EP analysis object
  TList*                            fListHistos;    //collection of output hists
  TFile*                            fSecondRunFile; //output from the second LYZ loop
      
  ClassDef(AliAnalysisTaskLYZEventPlane, 1); // example of analysis
};

#endif

