//
// Creation date: 2015/10/01
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#ifndef ALIREDUCEDANALYSISTEST_H
#define ALIREDUCEDANALYSISTEST_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliHistogramManager.h"

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
  void AddEventCut(AliReducedInfoCut* cut) {fEventCuts.Add(cut);}
  void AddTrackCut(AliReducedInfoCut* cut) {fTrackCuts.Add(cut);}
  void AddPairCut(AliReducedInfoCut* cut) {fPairCuts.Add(cut);}
  
  // getters
  virtual AliHistogramManager* GetHistogramManager() const {return fHistosManager;}
  
protected:
   AliHistogramManager* fHistosManager;   // Histogram manager
   
   TList fEventCuts;     // array of event cuts
   TList fTrackCuts;     // array of track cuts
   TList fPairCuts;      // array of pair cuts
   
  Bool_t IsEventSelected(AliReducedBaseEvent* event);
  Bool_t IsTrackSelected(AliReducedBaseTrack* track);
  Bool_t IsPairSelected(AliReducedBaseTrack* pair);
  
  ClassDef(AliReducedAnalysisTest,2);
};

#endif
