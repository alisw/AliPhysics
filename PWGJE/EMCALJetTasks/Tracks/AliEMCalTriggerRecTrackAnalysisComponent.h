/**
 * \file AliEMCalTriggerRecTrackAnalysisComponent.h
 * \brief Analysis component for reconstructed tracks
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALTRIGGERRECTRACKANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERRECTRACKANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEMCalTriggerTracksAnalysisComponent.h"

class TList;
class TString;
class AliVParticle;
class AliVTrack;
class AliMCEvent;
class AliEmcalTrackSelection;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerEventData;

/**
 * \class AliEMCalTriggerRecTrackAnalysisComponent
 * \brief Analysis component for reconstructed tracks
 *
 * Track analysis component: Loops over tracks from the EMCal track container and
 * counts the tracks in histograms. Separate histograms are filled for tracks which
 * have an EMCAL cluster assigned, and the Monte-Carlo information in case the analysis
 * is performed on Monte-Carlo events. In this case also the correlation matrix between
 * reconstructed \f$ p_{t} \f$ and generated \f$ p_{t} \f$ is filled.
 */
class AliEMCalTriggerRecTrackAnalysisComponent : public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerRecTrackAnalysisComponent();
  AliEMCalTriggerRecTrackAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerRecTrackAnalysisComponent();

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  /**
   * Defines swapping of the \f$ \eta \f$ sign in the analysis.
   *
   * \param doSwap If true, the \f$ \eta \f$ sign will be swapped
   */
  void SetSwapEta(Bool_t doSwap = kTRUE) { fSwapEta = doSwap; }

  /**
   * Defines whether tracks are required to be MC-true tracks, defined as track with
   * associated MC particle which fulfills the physical primary condition. Only relevant
   * for analysis of Monte-Carlo data sets
   *
   * \param doRequest If true, only MC-true tracks are selected
   */
  void SetRequestMCtrueTracks(Bool_t doRequest = kTRUE) { fRequestMCtrue = doRequest; }

  /**
   * Defines whether track patch matching is performed and correponding histograms are created and filled.
   *
   * \param doRequest
   */
  void SetRequestMatchTriggerPatches(Bool_t doRequest = kTRUE) { fDoMatchPatches = doRequest; }

  /**
   * Set the track selection applied in the analysis. Track cuts are always of type
   * AliEMCalPtTaskVTrackSelection. Transparent to ESD and AOD analysis.
   *
   * \param trackSel The (virtual) track selection applied in the analysis
   */
  void SetTrackSelection(AliEmcalTrackSelection *trackSel) { fTrackSelection = trackSel; }

protected:
  const AliVParticle *IsMCTrueTrack(const AliVTrack *const trk, const AliMCEvent *evnt) const;
  void FillHistogram(const TString &histname, const AliVTrack *const trk, const AliVParticle *assocMC, const AliVEvent * const recev, Bool_t useMCkine, Double_t weight);
  void FillCorrelation(const AliVParticle *const genparticle, const AliVParticle * const recparticle, double weight = 1.);
  void MatchTriggerPatches(const AliVTrack *rectrack, const TClonesArray *inputpatches, TList &outputpatches) const;
  Bool_t HasMatchedPatchOfType(TString triggertype, const TList & patches) const;

  AliEmcalTrackSelection *fTrackSelection;         ///< Track selection cuts used in the analysis
  Bool_t                            fSwapEta;                 ///< Swap eta sign
  Bool_t                            fRequestMCtrue;           ///< Request MC true track
  Bool_t                            fDoMatchPatches;          ///< Request matching with trigger patches

  ClassDef(AliEMCalTriggerRecTrackAnalysisComponent, 1);      // Analysis component for charged tracks
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERRECTRACKANALYSISCOMPONENT_H */
