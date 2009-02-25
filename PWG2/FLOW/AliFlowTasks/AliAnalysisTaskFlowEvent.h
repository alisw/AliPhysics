/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskFlowEvent_H
#define AliAnalysisTaskFlowEvent_H

// AliAnalysisTaskFlowEvent:
// analysis task to fill the flow event and make it available to the methods


class AliESDEvent;
class AliAODEvent;
class AliCFManager;
class AliFlowEventSimpleMaker;
class TList;

#include "TString.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskFlowEvent : public AliAnalysisTask {
 public:
  AliAnalysisTaskFlowEvent();
  AliAnalysisTaskFlowEvent(const char *name,Bool_t QAon = kFALSE);
  virtual ~AliAnalysisTaskFlowEvent();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisType(TString type) { this->fAnalysisType = type; }
  TString GetAnalysisType() const    { return this->fAnalysisType; }

  void SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr; } 
  AliCFManager* GetCFManager1()           {return this->fCFManager1; }
  void SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr; } 
  AliCFManager* GetCFManager2()           {return this->fCFManager2; }
  void          SetQAList1(TList* list)   {this->fQAInt = list; }
  TList*        GetQAList1()              {return this->fQAInt; }
  void          SetQAList2(TList* list)   {this->fQADiff = list; }
  TList*        GetQAList2()              {return this->fQADiff; }
  void          SetQAOn(Bool_t kt)        {this->fQA = kt; }
  Bool_t        GetQAOn()                 {return this->fQA; }

 private:

  AliAnalysisTaskFlowEvent(const AliAnalysisTaskFlowEvent& aAnalysisTask);
  AliAnalysisTaskFlowEvent& operator=(const AliAnalysisTaskFlowEvent& aAnalysisTask); 

  //  TFile*        fOutputFile;              // temporary output file for testing
  AliESDEvent*  fESD;                     // ESD object
  AliAODEvent*  fAOD;                     // AOD object
  AliFlowEventSimpleMaker* fEventMaker;   // FlowEventSimple maker object
  TString       fAnalysisType;            // can be MC, ESD or AOD
  AliCFManager* fCFManager1;              // correction framework manager
  AliCFManager* fCFManager2;              // correction framework manager
  TList*        fQAInt;                   // QA histogram list
  TList*        fQADiff;                  // QA histogram list

  Bool_t fQA;                             // flag to set the filling of the QA hostograms

  ClassDef(AliAnalysisTaskFlowEvent, 0); // example of analysis
};

#endif
