#ifndef ALIEMCALTRIGGERRECJETANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERRECJETANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliEMCalTriggerTracksAnalysisComponent.h"

class TString;
class AliEmcalJet;
class AliMCEvnet;
class AliVParticle;

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

class AliEMCalPtTaskVTrackSelection;
class AliEMCalTriggerEventData;

class AliEMCalTriggerRecJetAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerRecJetAnalysisComponent();
  AliEMCalTriggerRecJetAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerRecJetAnalysisComponent();

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  void SetUsePatches(Bool_t doUse = kTRUE) { fUsePatches = doUse; }
  void SetMinimumJetPt(Double_t minpt) { fMinimumJetPt = minpt; }
  void SetSingleTrackCuts(AliEMCalPtTaskVTrackSelection * trackcuts) { fTrackSelection = trackcuts; }
  void SetSwapEta(Bool_t doSwap = kTRUE) { fSwapEta = doSwap; }

protected:
  const AliVParticle * IsMCTrueTrack(const AliVTrack* const trk, const AliMCEvent* evnt) const;
  void FillHistogram(const TString &histname, const AliVParticle *track, const AliEmcalJet *jet, double vz, double weight);
  void FillJetHistogram(const TString &histname, const AliEmcalJet *recjet, double vz, double weight);
  void FillTrackHistogramCentrality(const TString &histname, const AliVTrack * const trk, const AliEmcalJet *jet, double centpercent, double weight);
  AliEMCalPtTaskVTrackSelection     *fTrackSelection;         // Track selection cuts used in the analysis
  Double_t                          fMinimumJetPt;            // Minimum jet pt
  Bool_t                            fRequestMCtrue;           // Request MC true track
  Bool_t                            fSwapEta;                 // Swap eta sign on request
  Bool_t                            fUsePatches;              // Use patches for trigger decision

  ClassDef(AliEMCalTriggerRecJetAnalysisComponent, 1);        // Analysis component for reconstructed Jets
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERRECJETANALYSISCOMPONENT_H */
