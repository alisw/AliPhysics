//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Analysis task for J/psi - hadron correlations                        //
//                                                                      //
// creation date: 07/12/2017                                            //
// authors:       I. Arsene, ionut.arsene@cern.ch                       //
//                L. Altenkämper, lucas.altenkamper@cern.ch             //
//                A. Lardeux, antoine.lardeux@cern.ch                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ALIREDUCEDANALYSISJPSI2EECORRELATIONS_H
#define ALIREDUCEDANALYSISJPSI2EECORRELATIONS_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedAnalysisJpsi2ee.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"
#include "AliMixingHandler.h"

//________________________________________________________________
class AliReducedAnalysisJpsi2eeCorrelations : public AliReducedAnalysisJpsi2ee {

public:
  AliReducedAnalysisJpsi2eeCorrelations();
  AliReducedAnalysisJpsi2eeCorrelations(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisJpsi2eeCorrelations();
  
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();

  // setters
  void AddMBEventCut(AliReducedInfoCut* cut) {fMBEventCuts.Add(cut);}
  void AddAssociatedTrackCut(AliReducedInfoCut* cut);
  void AddJpsiElectronCut(AliReducedInfoCut* cut) {AddTrackCut(cut);}      // synonim function to AddTrackCut
  void SetUseLikeSignPairs(Bool_t option) {fOptionUseLikeSignPairs = option;}
  void SetRunCorrelation(Bool_t option) {fOptionRunCorrelation = option;}
  void SetLoopOverTracks(Bool_t option) {
    fOptionLoopOverTracks = option;
    if(!fOptionLoopOverTracks) {fOptionRunPairing = kFALSE; fOptionRunMixing = kFALSE; fOptionUseLikeSignPairs = kFALSE; fOptionRunCorrelation = kFALSE;}
  }
  void SetRunCorrelationMixing(Bool_t option) {fOptionRunCorrelationMixing = option;};

  // getters
  Int_t GetNAssociatedTrackCuts() const {return fAssociatedTrackCuts.GetEntries();}
  const Char_t* GetAssociatedTrackCutName(Int_t i) const {return (i<fAssociatedTrackCuts.GetEntries() ? fAssociatedTrackCuts.At(i)->GetName() : "");}
  Bool_t GetRunCorrelation() {return fOptionRunCorrelation;}
  AliMixingHandler* GetCorrelationMixingHandler() const {return fCorrelationsMixingHandler;};
  Bool_t GetUseLikeSignPairs() const {return fOptionUseLikeSignPairs;}
  Bool_t GetRunCorrelationMixing() const {return fOptionRunCorrelationMixing;}
  
protected:
  AliMixingHandler* fCorrelationsMixingHandler;

  Bool_t fOptionUseLikeSignPairs;     // true: use like sign pairs for correlation
  Bool_t fOptionRunCorrelation;       // true: run correlation of candidates and associated tracks, false: apply only associated track cuts
  Bool_t fOptionRunCorrelationMixing; // true: run the event mixing for the correlation function

  TList fMBEventCuts;           // array of MB event cuts (for mixing with EMCal trigger)
  TList fAssociatedTrackCuts;   // array of associated track cuts
  TList fAssociatedTracks;      // list of selected associated tracks
  TList fAssociatedTracksMB;    // list of selected associated tracks in MB events (for mixing with EMCal trigger)

  Bool_t IsMBEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
  Bool_t IsAssociatedTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);

  void RunAssociatedTrackSelection(Bool_t fillHistograms=kTRUE, Bool_t fillMBTracks=kFALSE);
  void LoopOverAssociatedTracks(Int_t arrayOption=1, Bool_t fillHistograms=kTRUE, Bool_t fillMBTracks=kFALSE);
  void RunSameEventCorrelation(TString pairClass="PairSE");
  void FillAssociatedTrackHistograms(TString trackClass = "AssociatedTrack");
  void FillAssociatedTrackHistograms(AliReducedTrackInfo* track, TString trackClass = "AssociatedTrack");
  void FillCorrelationHistograms(AliReducedPairInfo* jpsi, AliReducedBaseTrack* assoc, TString corrClass="CorrSE", Bool_t isMCTruth=kFALSE);

  ClassDef(AliReducedAnalysisJpsi2eeCorrelations, 4);
};

#endif
