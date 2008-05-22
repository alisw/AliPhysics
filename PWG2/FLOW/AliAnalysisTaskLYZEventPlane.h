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

class AliESDEvent;
class AliAODEvent;
class AliFlowLYZEventPlane;
class AliFlowAnalysisWithLYZEventPlane;
class AliFlowEventSimpleMaker;
class TFile;

#include "AliAnalysisTask.h"

class AliAnalysisTaskLYZEventPlane : public AliAnalysisTask {
 public:
  AliAnalysisTaskLYZEventPlane(const char *name = "AliAnalysisTaskLYZEventPlane");
  virtual ~AliAnalysisTaskLYZEventPlane() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void           SetAnalysisType(TString type) {this->fAnalysisType = type ; }

 private:
 
  AliAnalysisTaskLYZEventPlane(const AliAnalysisTaskLYZEventPlane& aAnalysis);
  AliAnalysisTaskLYZEventPlane& operator=(const AliAnalysisTaskLYZEventPlane& aAnalysis);

  AliESDEvent *fESD;                      //ESD object
  AliAODEvent *fAOD;                      //AOD object
  TString fAnalysisType;                  //string to set the kind of input for the analysis: ESD, AOD or MC
  AliFlowLYZEventPlane* fLyzEp;           //LYZ EP object
  AliFlowAnalysisWithLYZEventPlane* fLyz; //LYZ EP analysis object
  AliFlowEventSimpleMaker* fEventMaker;   //FlowEventSimple maker object

  TFile* fFirstRunFile;                   //! output from the first LYZ loop
  TFile* fSecondRunFile;                  //! output from the second LYZ loop
    
  ClassDef(AliAnalysisTaskLYZEventPlane, 1); // example of analysis
};

#endif

