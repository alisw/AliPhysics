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
    kUseOnTheFlyReducedEvents=1,      // use events from exchange container published by the tree maker task
    kUseEventsFromTree=2                     // run directly over trees of reduced events
  };
  
 public:
  AliAnalysisTaskReducedEventProcessor();
  AliAnalysisTaskReducedEventProcessor(const char *name, Int_t runningMode=kUseEventsFromTree, Bool_t writeFilteredTree=kFALSE);
  virtual ~AliAnalysisTaskReducedEventProcessor(){}

  void AddTask(AliReducedAnalysisTaskSE* task) {fReducedTask=task;}

  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  
  virtual void ConnectInputData(Option_t *option = "");

  Int_t GetRunningMode() const {return fRunningMode;}  
  AliReducedAnalysisTaskSE* GetReducedTask() const {return fReducedTask;}
  
  Bool_t GetWriteFilteredTree() const {return fWriteFilteredTree;}
  
 protected:
  AliReducedAnalysisTaskSE* fReducedTask;      // Pointer to the analysis task which will process the reduced events
  
  Int_t fRunningMode;                               // Running mode, as specified in options 1 and 2 from Constants
  
  AliReducedBaseEvent* fReducedEvent;   //! reduced event
  
  Bool_t fWriteFilteredTree;                   // if kTRUE, the reduced task will produce filtered reduced trees
  
  AliAnalysisTaskReducedEventProcessor(const AliAnalysisTaskReducedEventProcessor &c);
  AliAnalysisTaskReducedEventProcessor& operator= (const AliAnalysisTaskReducedEventProcessor &c);

  ClassDef(AliAnalysisTaskReducedEventProcessor, 4);
};

#endif
