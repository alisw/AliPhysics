//
// Creation date: 2015/10/01
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#ifndef ALIREDUCEDANALYSISTEST_H
#define ALIREDUCEDANALYSISTEST_H

#include "AliReducedAnalysisTaskSE.h"


//________________________________________________________________
class AliReducedAnalysisTest : public AliReducedAnalysisTaskSE {
  
public:
  AliReducedAnalysisTest();
  AliReducedAnalysisTest(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisTest();
  
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();
  
  // setters
  
  // getters
  
  
protected:
  
  virtual Bool_t IsEventSelected(AliReducedBaseEvent* event);
  virtual Bool_t IsTrackSelected(AliReducedBaseTrack* track);
  virtual Bool_t IsPairSelected(AliReducedBaseTrack* pair);
  
  ClassDef(AliReducedAnalysisTest,1);
};

#endif
