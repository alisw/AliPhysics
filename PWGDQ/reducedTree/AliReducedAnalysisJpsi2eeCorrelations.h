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
  void AddAssociatedTrackCut(AliReducedInfoCut* cut);
  void SetElectronBit(Int_t bit) {fElectronBit = bit;}
  void SetRunCorrelation(Bool_t option) {fOptionRunCorrelation = option;}
  void SetAssociatedTracksSecondArray(Bool_t option) {fOptionAssociatedTracks = option;}
  void SetLoopOverTracks(Bool_t option) {
    fOptionLoopOverTracks = option;
    if(!fOptionLoopOverTracks) {fOptionRunPairing = kFALSE; fOptionRunMixing = kFALSE; fOptionRunLikeSignPairing = kFALSE; fOptionRunCorrelation = kFALSE;}
  }

  // getters
  Int_t GetNAssociatedTrackCuts() const {return fAssociatedTrackCuts.GetEntries();}
  const Char_t* GetAssociatedTrackCutName(Int_t i) const {return (i<fAssociatedTrackCuts.GetEntries() ? fAssociatedTrackCuts.At(i)->GetName() : "");}
  Int_t GetElectronBit() {return fElectronBit;}
  Bool_t GetRunCorrelation() {return fOptionRunCorrelation;}

protected:
  Int_t   fElectronBit;                           // toggled BIT in fQualityFlags for electron/positron selection
  Float_t fValues2[AliReducedVarManager::kNVars]; // array of values to hold information for associated tracks for histograms

  Bool_t fOptionRunCorrelation;   // true: run correlation of candidates and associated tracks, false: apply only associated track cuts
  Bool_t fOptionAssociatedTracks; // true: associated tracks are to be take from second track array, false: first track array

  TList fAssociatedTrackCuts;   // array of associated track cuts
  TList fAssociatedTracks;      // list of selected associated tracks

  Bool_t IsAssociatedTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);

  void RunAssociatedTrackSelection();
  void RunSameEventCorrelation(TString pairClass="PairSE");
  void FillAssociatedTrackHistograms(TString trackClass = "AssociatedTrack");
  void FillAssociatedTrackHistograms(AliReducedTrackInfo* track, TString trackClass = "AssociatedTrack");
  void FillCorrelationHistograms(ULong_t maskTrack, ULong_t maskAssocTrack, TString corrClass="CorrSE", Bool_t isMCTruth=kFALSE);

  ClassDef(AliReducedAnalysisJpsi2eeCorrelations, 1);
};

#endif
