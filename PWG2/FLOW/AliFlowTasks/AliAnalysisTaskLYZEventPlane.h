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
class AliCFManager;
class AliFlowLYZEventPlane;
class AliFlowAnalysisWithLYZEventPlane;
class AliFlowEventSimpleMaker;
class TFile;
class TList;

#include "TString.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskLYZEventPlane : public AliAnalysisTask {
 public:
  AliAnalysisTaskLYZEventPlane();
  AliAnalysisTaskLYZEventPlane(const char *name, Bool_t QAon = kFALSE);
  virtual ~AliAnalysisTaskLYZEventPlane();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void           SetAnalysisType(TString type) {this->fAnalysisType = type ; }
  TString        GetAnalysisType() const       {return this->fAnalysisType; }

  void           SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr; }
  AliCFManager*  GetCFManager1() const         {return this->fCFManager1; }
  void           SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr; }
  AliCFManager*  GetCFManager2() const         {return this->fCFManager2; }
  void           SetQAList1(TList* list)       {this->fQAInt = list; }
  TList*         GetQAList1()                  {return this->fQAInt; }
  void           SetQAList2(TList* list)       {this->fQADiff = list; }
  TList*         GetQAList2()                  {return this->fQADiff; }
  void           SetQAOn(Bool_t kt)            {this->fQA = kt; }
  Bool_t         GetQAOn()                     {return this->fQA; }

 private:
 
  AliAnalysisTaskLYZEventPlane(const AliAnalysisTaskLYZEventPlane& aAnalysis);
  AliAnalysisTaskLYZEventPlane& operator=(const AliAnalysisTaskLYZEventPlane& aAnalysis);

  AliESDEvent*  fESD;                     //ESD object
  AliAODEvent*  fAOD;                     //AOD object
  TString fAnalysisType;                  //string to set the kind of input for the analysis: ESD, AOD or MC
  AliFlowLYZEventPlane* fLyzEp;           //LYZ EP object
  AliFlowAnalysisWithLYZEventPlane* fLyz; //LYZ EP analysis object
  AliFlowEventSimpleMaker* fEventMaker;   //FlowEventSimple maker object
  AliCFManager* fCFManager1;              //Correction framework manager
  AliCFManager* fCFManager2;              //Correction framework manager
  TList*        fListHistos;              //collection of output hists
  TFile*        fSecondRunFile;           //output from the second LYZ loop
  TList*        fQAInt;                   // QA histogram list
  TList*        fQADiff;                  // QA histogram list

  Bool_t        fQA;                      // flag to set the filling of the QA hostograms
    
  ClassDef(AliAnalysisTaskLYZEventPlane, 1); // example of analysis
};

#endif

