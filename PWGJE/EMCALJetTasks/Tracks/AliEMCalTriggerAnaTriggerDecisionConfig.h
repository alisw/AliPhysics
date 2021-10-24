#ifndef ALIEMCALTRIGGERANATRIGGERDECISIONCONFIG_H
#define ALIEMCALTRIGGERANATRIGGERDECISIONCONFIG_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include "AliEMCalTriggerAnaHelper.h"

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerAnaTriggerDecisionConfig : public TObject {
public:
  AliEMCalTriggerAnaTriggerDecisionConfig();

  /**
   * Destructor, nothing to do.
   */
  virtual ~AliEMCalTriggerAnaTriggerDecisionConfig() {}

  /**
   * Define whether we swap the low and high energy thresholds (necessary for Monte-Carlo patches)
   * \param doSwap If true we swap the thresholds
   */
  void SetSwapThresholds(Bool_t doSwap = kTRUE) { fSwapThresholds = doSwap;}

  /**
   * Define whether we use online or offline patches for the trigger decision
   * \param useOffline If true we use offline patches instead of online patches.
   */
  void SetUseOfflinePatches(Bool_t doUse = kTRUE ) { fUseOfflinePatches = doUse; }

  /**
   * Set the energy threshold for the trigger selection using patches for a given trigger class
   * \param triggerClass Trigger class for which to set the threshold
   * \param threshold New trigger threshold
   */
  void SetEnergyThreshold(ETATriggerType trigger, double threshold){
    fEnergyThresholds[trigger] = threshold;
  }

  /**
   * Specify which type of energy is used for the patch energy cut
   * \param energyType Energy type used in the selection
   */
  void SetPatchEnergyType(EPatchEnergyType_t energyType) { fEnergyType = energyType; }

  /**
   * Check whether energy threshols are swapped (only relevant for online patches in MC)
   * \return True if thresholds are required to be swapped
   */
  Bool_t IsSwapThresholds() const { return fSwapThresholds; }

  /**
   * Check whether offline patches are used in the analysis.
   * \return True if we use offline patches
   */
  Bool_t IsUsingOfflinePatches() const { return fUseOfflinePatches; }

  /**
   * Get the energy threshold for a trigger class for further patch selection
   * \param trigger Trigger class
   * \return The energy threshold
   */
  Double_t GetEnergyThreshold(ETATriggerType trigger) const {
    return fEnergyThresholds[trigger];
  }

  /**
   * Check if an energy threshold is defined for a given trigger class type (> 0)
   * \param trigger Trigger class type
   * \return True if a threshold is set
   */
  Bool_t HasEnergyThreshold(ETATriggerType trigger) const {
    return fEnergyThresholds[trigger] > 0;
  }

  /**
   * Get the type of energy used in the trigger selection.
   * \return Type of the energy used
   */
  EPatchEnergyType_t GetPatchEnergyType() const { return fEnergyType; }

private:
  Bool_t                        fSwapThresholds;                            ///< Flag for swapping high and low energy threshold
  Bool_t                        fUseOfflinePatches;                         ///< Switch for using offline patches for event selection
  Double_t                      fEnergyThresholds[4];                       ///< Energy thresholds applied in the analysis
  EPatchEnergyType_t            fEnergyType;                                ///< Energy type from patch used for the patch energy selection

  ClassDef(AliEMCalTriggerAnaTriggerDecisionConfig, 1);
};

}

}
#endif /* */
