//
// Creation date: 2016/09/06
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#ifndef ALIREDUCEDANALYSISJPSI2EE_H
#define ALIREDUCEDANALYSISJPSI2EE_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliHistogramManager.h"
#include "AliMixingHandler.h"

//________________________________________________________________
class AliReducedAnalysisJpsi2ee : public AliReducedAnalysisTaskSE {
  
public:
  AliReducedAnalysisJpsi2ee();
  AliReducedAnalysisJpsi2ee(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisJpsi2ee();
  
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
  virtual AliMixingHandler* GetMixingHandler() const {return fMixingHandler;}
  
protected:
   AliHistogramManager* fHistosManager;   // Histogram manager
   AliMixingHandler*         fMixingHandler;    // mixing handler
   
   TList fEventCuts;               // array of event cuts
   TList fTrackCuts;               // array of track cuts
   TList fPreFilterTrackCuts;  // track cuts to be used at the prefilter stage
   TList fPairCuts;                  // array of pair cuts
   TList fPreFilterPairCuts;     // pair cuts to be used at the prefilter stage
   
  Bool_t IsEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
  Bool_t IsTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
  Bool_t IsTrackPrefilterSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
  Bool_t IsPairSelected(AliReducedBaseTrack* pair, Float_t* values=0x0);
  Bool_t IsPairPreFilterSelected(AliReducedBaseTrack* pair, Float_t* values=0x0);
  
  ClassDef(AliReducedAnalysisJpsi2ee,1);
};

#endif
