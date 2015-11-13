/*
Creation date: 2015/10/01
Author: Jacobus Onderwaater, Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
 */


#ifndef ALIREDUCEDANALYSISQNCORRECTIONS_H
#define ALIREDUCEDANALYSISQNCORRECTIONS_H

#include "AliReducedAnalysisTaskSE.h"


class AliQnCorrectionsManager;
class AliReducedEventInfo;
class AliReducedQnFillEvent;


//________________________________________________________________
class AliReducedAnalysisQnCorrections : public AliReducedAnalysisTaskSE {
  
public:
  AliReducedAnalysisQnCorrections();
  AliReducedAnalysisQnCorrections(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisQnCorrections();
  
  // Virtual functions, to be implemented in the inheriting analysis classes
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();
  
  void SetQnCorrectionsManager( AliQnCorrectionsManager* man) {fQnManager=man;}


  // setters
  
  // getters
  
  
protected:
  
  virtual Bool_t IsEventSelected(AliReducedBaseEvent* event);
  virtual Bool_t IsTrackSelected(AliReducedBaseTrack* track);
  virtual Bool_t IsPairSelected(AliReducedBaseTrack* pair);

  AliReducedQnFillEvent* fQnFillEvent;
  AliQnCorrectionsManager* fQnManager;

  ClassDef(AliReducedAnalysisQnCorrections,1);
};

#endif
