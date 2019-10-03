/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUpcTriggerCounter_H
#define AliAnalysisTaskUpcTriggerCounter_H

class TH1;
class TH2;
class TList;
class AliAODEvent;

#define ntrig 13
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpcTriggerCounter : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcTriggerCounter();
  AliAnalysisTaskUpcTriggerCounter(const char *name);
  virtual ~AliAnalysisTaskUpcTriggerCounter();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
 
  virtual void Terminate(Option_t *);

 private:
  TList *fListHist;		//<
  
  TH1I *fHistTriggerCounter;    //!
  TH2I *fHistTriggerCounterIR1; //!
  TH2I *fHistTriggerCounterIR2; //!

  AliAnalysisTaskUpcTriggerCounter(const AliAnalysisTaskUpcTriggerCounter&); //not implemented
  AliAnalysisTaskUpcTriggerCounter& operator =(const AliAnalysisTaskUpcTriggerCounter&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcTriggerCounter, 2); 
};

#endif
