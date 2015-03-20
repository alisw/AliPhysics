/**
 * \file AliEMCalTriggerPatchAnalysisComponent.h
 * \brief Analysis component for EMCAL trigger patches
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALTRIGGERPATCHANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERPATCHANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliEMCalTriggerTracksAnalysisComponent.h"

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData;

/**
 * \class AliEMCalTriggerPatchAnalysisComponent
 * \brief Analysis component for EMCAL trigger patches
 *
 * Analysis components for trigger patches. Fills THnSparses with the different energy definitions
 * (amplitude, estimated (rough) patch energy, calibrated (offline) energy and the patch position on
 * the EMCAL surface for EMCAL trigger patches of different categories created by the EMCAL trigger
 * patch maker. Main patches defined as the patches of a given type with the highest energy.
 */
class AliEMCalTriggerPatchAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerPatchAnalysisComponent();
  AliEMCalTriggerPatchAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerPatchAnalysisComponent() { }

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  /**
   * Check whether thresholds for online patches are swapped
   *
   * \return true if thresholds are swapped, fales otherwise
   */
  Bool_t IsSwapOnlineThresholds() const { return fSwapOnlineThresholds; }

  /**
   * Check whether thresholds for offline patches are swapped
   *
   * \return true if thresholds are swapped, fales otherwise
   */
  Bool_t IsSwapOfflineThresholds() const { return fSwapOfflineThresholds; }

  /**
   * Swap online thresholds (i.e. low threshold becomes high threshold and vice versa)
   *
   * \param doSwap Do swap the thresholds
   */
  void SetSwapOnlineThresholds(Bool_t doSwap = kTRUE) { fSwapOnlineThresholds = doSwap; }

  /**
   * Swap offline thresholds (i.e. low threshold becomes high threshold and vice versa)
   *
   * \param doSwap Do swap the thresholds
   */
  void SetSwapOfflineThresholds(Bool_t doSwap = kTRUE) { fSwapOfflineThresholds = doSwap; }

protected:

  Bool_t                        fSwapOnlineThresholds;          ///< Swap trigger thresholds for online patches
  Bool_t                        fSwapOfflineThresholds;         ///< Swap trigger thresholds for offline patches

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerPatchAnalysisComponent, 1);     // Component for trigger patch analysis
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERPATCHANALYSISCOMPONENT_H */
