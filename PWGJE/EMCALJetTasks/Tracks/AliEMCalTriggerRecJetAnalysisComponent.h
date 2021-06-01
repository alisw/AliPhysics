/**
 * \file AliEMCalTriggerRecJetAnalysisComponent.h
 * \brief Declaration of the analysis component on reconstructed jets
 *
 * This class defines the analysis components of reconstructed particles
 * in reconstructed jets.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALTRIGGERRECJETANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERRECJETANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEMCalTriggerTracksAnalysisComponent.h"

class TString;
class AliEmcalJet;
class AliMCEvnet;
class AliVParticle;

class AliEmcalTrackSelection;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerEventData;

/**
 * \class AliEMCalTriggerRecJetAnalysisComponent
 * \brief Analysis component for tracks in reconstructed jets
 *
 * Analysis component for tracks in reconstructed jets: Connects reconstructed
 * particles to jets, and fills track based histograms for tracks in jet with
 * a minimum jet \f$ p_{t} \f$.
 */
class AliEMCalTriggerRecJetAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerRecJetAnalysisComponent();
  AliEMCalTriggerRecJetAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerRecJetAnalysisComponent();

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  /**
   * Set the minimum \f$ p_{t} \f$ allowed to select reconstructed jets.
   * \param minpt The minimum jet \f$ p_{t} \f$
   */
  void SetMinimumJetPt(Double_t minpt) { fMinimumJetPt = minpt; }

  /**
   * Set quality cuts used to select reconstructed tracks.
   * \param trackcuts Track selection cuts used in this analysis component.
   */
  void SetSingleTrackCuts(AliEmcalTrackSelection * trackcuts) { fTrackSelection = trackcuts; }

  /**
   * Defines whether we swap the sign of \f$ eta\f$.
   * \param doSwap If true we swap the sign of \f$ eta\f$
   */
  void SetSwapEta(Bool_t doSwap = kTRUE) { fSwapEta = doSwap; }

protected:
  const AliVParticle * IsMCTrueTrack(const AliVTrack* const trk, const AliMCEvent* evnt) const;
  void FillHistogram(const TString &histname, const AliVParticle *track, const AliEmcalJet *jet, double vz, double weight);
  void FillJetHistogram(const TString &histname, const AliEmcalJet *recjet, double vz, double weight);
  void FillTrackHistogramCentrality(const TString &histname, const AliVTrack * const trk, const AliEmcalJet *jet, double centpercent, double weight);
  AliEmcalTrackSelection  *fTrackSelection;         ///< Track selection cuts used in the analysis
  Double_t                          fMinimumJetPt;            ///< Minimum jet \f$ p_{t} \f$
  Bool_t                            fRequestMCtrue;           ///< Request MC true track
  Bool_t                            fSwapEta;                 ///< Swap eta sign on request

  ClassDef(AliEMCalTriggerRecJetAnalysisComponent, 1);        // Analysis component for reconstructed Jets
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERRECJETANALYSISCOMPONENT_H */
