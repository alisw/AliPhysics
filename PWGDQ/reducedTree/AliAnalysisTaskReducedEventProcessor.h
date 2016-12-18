/*
***********************************************************
  Wrapper class for AliReducedAnalysisTaskSE to be used in the AliAnalysisTask framework
  contact: 
  Ionut-Cristian Arsene, iarsene@cern.ch
  Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  2015/10/01
  *********************************************************
*/

#ifndef ALIANALYSISTASKREDUCEDEVENTPROCESSOR_H
#define ALIANALYSISTASKREDUCEDEVENTPROCESSOR_H

#include "AliAnalysisTaskSE.h"
#include "AliReducedBaseEvent.h"

class TObject;
class AliAnalysis;
class AliReducedAnalysisTaskSE;

//_________________________________________________________
class AliAnalysisTaskReducedEventProcessor : public AliAnalysisTaskSE {

public:
  enum Constants {
    //kMaxOutputs=100,
    kUseOnTheFlyReducedEvents=1,      // use events from exchange container published by the tree maker task
    kUseEventsFromTree=2                     // run directly over trees of reduced events
  };
  
 public:
  AliAnalysisTaskReducedEventProcessor();
  AliAnalysisTaskReducedEventProcessor(const char *name, Int_t runningMode=kUseEventsFromTree);
  virtual ~AliAnalysisTaskReducedEventProcessor(){}

  void AddTask(AliReducedAnalysisTaskSE* task) {fReducedTask=task;}

  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  
  virtual void ConnectInputData(Option_t *option = "");

  Int_t GetRunningMode() const {return fRunningMode;}  
  AliReducedAnalysisTaskSE* GetReducedTask() const {return fReducedTask;}
  
 protected:
  AliReducedAnalysisTaskSE* fReducedTask;      // Pointer to the analysis task which will process the reduced events
 // TObject* fOutputSlot[kMaxOutputs];
 // Int_t fContainerType[kMaxOutputs]; // 0: output   1: exchange
 // Int_t fNoutputSlots;
  
  Int_t fRunningMode;                               // Running mode, as specified in options 1 and 2 from Constants
  Long_t fEventNumber;                  // event number
  
  AliReducedBaseEvent* fReducedEvent;   //! reduced event
  
  AliAnalysisTaskReducedEventProcessor(const AliAnalysisTaskReducedEventProcessor &c);
  AliAnalysisTaskReducedEventProcessor& operator= (const AliAnalysisTaskReducedEventProcessor &c);

  ClassDef(AliAnalysisTaskReducedEventProcessor, 3);
};

#endif
