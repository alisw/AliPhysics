//
// Creation date: 2016/09/06
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#ifndef ALIREDUCEDANALYSISJPSI2EEMULT_H
#define ALIREDUCEDANALYSISJPSI2EEMULT_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliHistogramManager.h"
#include "AliMixingHandler.h"


//________________________________________________________________
class AliReducedAnalysisJpsi2eeMult : public AliReducedAnalysisTaskSE {
  
public:
  AliReducedAnalysisJpsi2eeMult();
  AliReducedAnalysisJpsi2eeMult(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisJpsi2eeMult();
  
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();
  
  // setters
  void AddEventCut(AliReducedInfoCut* cut) {fEventCuts.Add(cut);}
  void AddTrackCut(AliReducedInfoCut* cut);
  void AddPrefilterTrackCut(AliReducedInfoCut* cut) {fPreFilterTrackCuts.Add(cut);}
  void AddPairCut(AliReducedInfoCut* cut) {fPairCuts.Add(cut);}
  void AddPrefilterPairCut(AliReducedInfoCut* cut) {fPreFilterPairCuts.Add(cut);}
  void SetRunEventMixing(Bool_t option) {fOptionRunMixing = option;};
  void SetRunPairing(Bool_t option) {fOptionRunPairing = option;};
  void SetRunOverMC(Bool_t option) {fOptionRunOverMC = option;};
  void SetRunLikeSignPairing(Bool_t option) {fOptionRunLikeSignPairing = option;}
  void SetLoopOverTracks(Bool_t option) {
     fOptionLoopOverTracks = option; 
     if(!fOptionLoopOverTracks) {fOptionRunPairing = kFALSE; fOptionRunMixing = kFALSE; fOptionRunLikeSignPairing = kFALSE;}     
  }
  
  // getters
  virtual AliHistogramManager* GetHistogramManager() const {return fHistosManager;}
  virtual AliMixingHandler* GetMixingHandler() const {return fMixingHandler;}
  Int_t GetNTrackCuts() const {return fTrackCuts.GetEntries();}
  const Char_t* GetTrackCutName(Int_t i) const {return (i<fTrackCuts.GetEntries() ? fTrackCuts.At(i)->GetName() : "");} 
  Bool_t GetRunOverMC() {return fOptionRunOverMC;};
  Bool_t GetRunLikeSignPairing() {return fOptionRunLikeSignPairing;}
  Bool_t GetRunEventMixing() {return fOptionRunMixing;}
  Bool_t GetRunPairing() {return fOptionRunPairing;}
  Bool_t GetLoopOverTracks() {return fOptionLoopOverTracks;}
  
protected:
   AliHistogramManager* fHistosManager;   // Histogram manager
   AliMixingHandler*         fMixingHandler;    // mixing handler
   
   Bool_t fOptionRunMixing;    // true: run event mixing, false: no event mixing
   Bool_t fOptionRunPairing;    // true: run pairing, false: only apply the track cuts
   Bool_t fOptionRunOverMC;  // true: trees contain MC info -> fill histos to compute efficiencies, false: run normally as on data
   Bool_t fOptionRunLikeSignPairing;   // true (default): performs the like sign pairing in addition to the opposite pairing
   Bool_t fOptionLoopOverTracks;       // true (default); if false do not loop over tracks and consequently no pairing
  
   TList fEventCuts;               // array of event cuts
   TList fTrackCuts;               // array of track cuts
   TList fPreFilterTrackCuts;  // track cuts to be used at the prefilter stage
   TList fPairCuts;                  // array of pair cuts
   TList fPreFilterPairCuts;     // pair cuts to be used at the prefilter stage
   
   TList fPosTracks;               // list of selected positive tracks in the current event
   TList fNegTracks;              // list of selected negative tracks in the current event
   TList fPrefilterPosTracks;  // list of prefilter selected positive tracks in the current event
   TList fPrefilterNegTracks; // list of prefilter selected negative tracks in the current event
   
   ULong_t fEventCounter;   // event counter
   
  Bool_t IsEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
  Bool_t IsTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
  Bool_t IsTrackPrefilterSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
  Bool_t IsPairSelected(Float_t* values);
  Bool_t IsPairPreFilterSelected(Float_t* values);
  Bool_t IsMCTruth(AliReducedTrackInfo* ptrack, AliReducedTrackInfo* ntrack);
  Bool_t IsMCTruth(AliReducedTrackInfo* track);
  void    FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1, Int_t& leg2);
  
  void RunPrefilter();
  void RunSameEventPairing(TString pairClass = "PairSE");
  void RunTrackSelection();
  void FillTrackHistograms(TString trackClass = "Track");
  void FillTrackHistograms(AliReducedTrackInfo* track, TString trackClass = "Track");
  void FillPairHistograms(ULong_t mask, Int_t pairType, TString pairClass = "PairSE", Bool_t isMCTruth = kFALSE);
  void FillMCTruthHistograms();
  
  ClassDef(AliReducedAnalysisJpsi2eeMult,3);
};

#endif
