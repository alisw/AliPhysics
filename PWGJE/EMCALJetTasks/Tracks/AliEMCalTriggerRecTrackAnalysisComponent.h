#ifndef ALIEMCALTRIGGERRECTRACKANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERRECTRACKANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <Tracks/AliEMCalTriggerTracksAnalysisComponent.h>

class TString;
class AliVParticle;
class AliVTrack;
class AliMCEvent;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData;
class AliEMCalPtTaskVTrackSelection;

class AliEMCalTriggerRecTrackAnalysisComponent : public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerRecTrackAnalysisComponent();
  AliEMCalTriggerRecTrackAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerRecTrackAnalysisComponent();

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  void SetSwapEta(Bool_t doSwap = kTRUE) { fSwapEta = doSwap; }
  void SetUsePatches(Bool_t doUse = kTRUE) { fUsePatches = doUse; }
  void SetRequestMCtrueTracks(Bool_t doRequest = kTRUE) { fRequestMCtrue = doRequest; }
  void SetTrackSelection(AliEMCalPtTaskVTrackSelection *trackSel) { fTrackSelection = trackSel; }

protected:
  const AliVParticle *IsMCTrueTrack(const AliVTrack *const trk, const AliMCEvent *evnt) const;
  void FillHistogram(const TString &histname, const AliVTrack *const trk, const AliVParticle *assocMC, const AliVEvent * const recev, Bool_t useMCkine, Double_t weight);
  void FillCorrelation(const AliVParticle *const genparticle, const AliVParticle * const recparticle, double weight = 1.);

  AliEMCalPtTaskVTrackSelection *   fTrackSelection;          // Track selection cuts used in the analysis
  Bool_t                            fSwapEta;                 // Swap eta sign
  Bool_t                            fUsePatches;              // Use patches for trigger decision
  Bool_t                            fRequestMCtrue;           // Request MC true track

  ClassDef(AliEMCalTriggerRecTrackAnalysisComponent, 1);      // Analysis component for charged tracks
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERRECTRACKANALYSISCOMPONENT_H */
