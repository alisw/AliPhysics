///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Analysis task for single tracks                                       //
//                                                                       //
// creation date: 10/10/2018                                             //
// author: Lucas Altenkamper, lucas.altenkamper@cern.ch                  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIREDUCEDANALYSISSINGLETRACK_H
#define ALIREDUCEDANALYSISSINGLETRACK_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliHistogramManager.h"
#include "AliMixingHandler.h"

//________________________________________________________________
class AliReducedAnalysisSingleTrack : public AliReducedAnalysisTaskSE {
  
public:
  AliReducedAnalysisSingleTrack();
  AliReducedAnalysisSingleTrack(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisSingleTrack();
  
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();

  // setters
  void AddEventCut(AliReducedInfoCut* cut) {fEventCuts.Add(cut);}
  void AddTrackCut(AliReducedInfoCut* cut) {fTrackCuts.Add(cut);}
  void AddMCSignalCut(AliReducedInfoCut* cut) {if (fMCSignalCuts.GetEntries()>=32) return; fMCSignalCuts.Add(cut);}
  void SetRunOverMC(Bool_t option) {fOptionRunOverMC = option;};

  //getters
  virtual AliHistogramManager*  GetHistogramManager() const {return fHistosManager;}
  Int_t                         GetNTrackCuts() const {return fTrackCuts.GetEntries();}
  const Char_t*                 GetTrackCutName(Int_t i) const {return (i<fTrackCuts.GetEntries() ? fTrackCuts.At(i)->GetName() : "");}
  Bool_t                        GetRunOverMC() const {return fOptionRunOverMC;};
  Int_t                         GetNMCSignalCuts() const {return fMCSignalCuts.GetEntries();}
  const Char_t*                 GetMCSignalCutName(Int_t i) const {return (i<fMCSignalCuts.GetEntries() ? fMCSignalCuts.At(i)->GetName() : "");}
  
protected:
  AliHistogramManager* fHistosManager;   // Histogram manager
  
  Bool_t fOptionRunOverMC;  // true: trees contain MC info -> fill histos to compute efficiencies, false: run normally as on data
  
  TList fEventCuts;     // array of event cuts
  TList fTrackCuts;     // array of track cuts
  TList fMCSignalCuts;  // array of MC signal cuts

  TList fTracks; // list of selected tracks

  Bool_t  IsEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
  Bool_t  IsTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
  UInt_t  CheckTrackMCTruth(AliReducedBaseTrack* track);
  void    FillMCTruthHistograms();
  void    RunTrackSelection();
  void    FillTrackHistograms(TString trackClass="Track");
  void    FillTrackHistograms(AliReducedBaseTrack* track, TString trackClass="Track");

  ClassDef(AliReducedAnalysisSingleTrack,1);
};

#endif

