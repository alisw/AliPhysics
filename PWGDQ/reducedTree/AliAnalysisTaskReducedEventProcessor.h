/*
***********************************************************
  Wrapper class for AliReducedAnalysisTaskSE to be used in the AliAnalysisTask framework
  contact: jaap onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  2015/10/01
  *********************************************************
*/

//#include "AliSysInfo.h"

#ifndef ALIANALYSISTASKREDUCEDEVENTPROCESSOR_H
#define ALIANALYSISTASKREDUCEDEVENTPROCESSOR_H

#include "AliAnalysisTaskSE.h"

class TObject;
class AliAnalysis;
class AliReducedAnalysisTaskSE;

//_________________________________________________________
class AliAnalysisTaskReducedEventProcessor : public AliAnalysisTaskSE {

  enum Constants {
    kMaxOutputs=100
  };

 public:
  AliAnalysisTaskReducedEventProcessor();
  AliAnalysisTaskReducedEventProcessor(const char *name);
  virtual ~AliAnalysisTaskReducedEventProcessor(){}

  void AddTask(AliReducedAnalysisTaskSE* task) {fReducedTask=task;}

  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();

 private:
  AliReducedAnalysisTaskSE* fReducedTask;
  TObject* fOutputSlot[kMaxOutputs];
  Int_t fContainerType[kMaxOutputs]; // 0: output   1: exchange
  Int_t fNoutputSlots;
  
  AliAnalysisTaskReducedEventProcessor(const AliAnalysisTaskReducedEventProcessor &c);
  AliAnalysisTaskReducedEventProcessor& operator= (const AliAnalysisTaskReducedEventProcessor &c);

  ClassDef(AliAnalysisTaskReducedEventProcessor, 1);
};

#endif


