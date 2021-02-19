#ifndef ALIEMCALTRIGGEROFFLINESELECTION_H
#define ALIEMCALTRIGGEROFFLINESELECTION_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TString.h>

class TClonesArray;
class TH2;
class AliVEvent;

namespace PWGJE {
  
namespace EMCALJetTasks {

/**
 * @class AliEmcalTriggerOfflineSelection
 * @brief Helper class selecting events on the presence of a trigger patch for the given type above threshold
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since Feb 23, 2016
 *
 * AliEmcalTriggerOffline selection provides the functionality to select triggered events
 * based on the presence of a trigger patch for the given type above energy threshold. As convention
 * the calibrated energy from cells is used to select patches above threshold. Energy thresholds
 * can be set via
 *
 * ~~~{.cxx}
 * AliEmcalTriggerOfflineSelection sel
 * // Setting 10 GeV threshold for the L1 gamma trigger
 * sel.SetOfflineEnergyThreshold(AliEmcalTriggerOfflineSelection::kTrgEG1, 10.);
 * ~~~
 *
 * In order to mimic non-uniform trigger acceptance, acceptance maps from Data can be
 * used in the trigger patch selection. When available, trigger patches are selected
 * on a statistical basis using the map, which is expected to be normalized to the
 * position with the maximum trigger counts, in order to provide a probability value.
 * Acceptance maps are set via
 *
 * ~~~{.cxx}
 * TH2 *accEG1; // User is responsible to load the histogram
 * sel.SetAcceptanceMap(AliEmcalTriggerOfflineSelection::kTrgEG1, accEG1);
 * ~~~
 *
 * Attention: This class takes over ownership and expects the histogram to not belong to a
 * directory. In case no acceptance map is providede, no acceptance mimicing is applied.
 */
class AliEmcalTriggerOfflineSelection : public TObject {
public:
  /**
   * @enum EmcalTriggerClass
   * @brief Definition of the various supported trigger types
   */
  enum EmcalTriggerClass{
    kTrgEL0 = 0,      ///< EMCAL L0 trigger
    kTrgEG1,          ///< EMCAL L1 Gamma trigger, high threshold
    kTrgEG2,          ///< EMCAL L1 Gamma trigger, low threshold
    kTrgEJ1,          ///< EMCAL L1 Jet trigger, high threshold
    kTrgEJ2,          ///< EMCAL L1 Jet trigger, low threshold
    kTrgDL0,          ///< DCAL L0 trigger
    kTrgDG1,          ///< DCAL L1 Gamma trigger, high threshold
    kTrgDG2,          ///< DCAL L1 Gamma trigger, low threshold
    kTrgDJ1,          ///< DCAL L1 Jet trigger, high threshold
    kTrgDJ2,          ///< DCAL L1 Jet trigger, low threshold
    kTrgn             ///< Number of supported triggers
  };

  /**
   * @enum EmcalEnergyDefinition_t
   * @brief Definition of EMCAL patch energy measurements
   */
  enum EmcalEnergyDefinition_t {
     kFEEEnergy,                ///< FEE energy
     kFEETransverseEnergy,      ///< FEE transverse energy
     kFEEADC,                   ///< FEE energy converted to L1 ADC
     kFEETransverseADC,         ///< FEE transverse converted to L1 transverse ADC
     kClusterEnergy,            ///< Cluster energy
     kClusterTransverseEnergy   ///< Cluster transverse energy
  };

  /**
   * @brief Default constructor
   */
  AliEmcalTriggerOfflineSelection();

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalTriggerOfflineSelection();

  /**
   * @brief Specify threshold for a given offline trigger class.
   *
   * Convention is a threshold on the patch energy (from cells) in GeV
   *
   * @param[in] trgcls Trigger class for which to set the threshold
   * @param[in] threshold Threshold value for the trigger class
   */
  void SetOfflineEnergyThreshold(EmcalTriggerClass trgcls, double threshold) { fOfflineEnergyThreshold[trgcls] = threshold; }

  /**
   * @brief Define according to which energy measurement a patch is selected as trigger patch.
   *
   * See @ref EmcalEnergyDefinition_t for possible settings.
   *
   * @param[in] endef Type of the energy measurement.
   */
  void SetEnergyDefinition(EmcalEnergyDefinition_t endef) { fEnergyDefinition = endef; }

  /**
   * @brief Set acceptance map for a given trigger class.
   *
   * The acceptance map is expected to be normalized to 1 for the position with the largest
   * trigger patch.
   *
   * @param[in] trgcls Trigger class for which to set then acceptance map
   * @param[in] accmap Acceptance map as 2D histogram
   */
  void SetAcceptanceMap(EmcalTriggerClass trgcls, const TH2 *accmap) { fAcceptanceMaps[trgcls] = accmap; }

  /**
   * Define name of the cluster container (for cluster trigger)
   * @param[in] clustercont Name of the cluster container
   */
  void SetClusterContainer(const TString &clustercont) { fNameClusterContainer = clustercont; }

  /**
   * @brief Select event as triggered event.
   *
   * Apply additional cut requiring at least one offline patch above a given energy (not fake ADC!)
   * Attention: This task groups into single shower triggers (L0, EG1, EG2) and jet triggers (EJ1 and EJ2).
   * Per convention the low threshold patch is selected. No energy cut should be applied in the trigger maker.
   *
   * Patches are supposed to contain masking done in the trigger maker. If available, an acceptance map can be used
   * to mimic patch azimuthal acceptance in simulation.
   *
   * Attention: The component take over ownership of the histogram. Need to be set to nullptr once
   * the histogram is read in.
   *
   * @param[in] trgcls Trigger class for which to apply additional offline patch selection
   * @param[in] triggerpatches Array of trigger patches
   * @return True if at least on patch above threshold is found or no cut is applied
   */
  Bool_t IsOfflineSelected(EmcalTriggerClass trgcls, const AliVEvent * const data) const;

  /**
   * @brief Get the offline trigger threshold (on energy) for a given trigger class.
   * @param[in] trgcls Trigger class for which to check the threshold
   * @return Threshold for the given trigger class (0 if not set)
   */
  Double_t GetThresholdForTrigger(EmcalTriggerClass trgcls) const {return fOfflineEnergyThreshold[trgcls]; }

  /**
   * @brief Get the name of the cluster container
   * @return Name of the cluster container
   */
  const TString &GetNameClusterContainer() const  { return fNameClusterContainer; }

  /**
   * @brief Checks if the trigger class is a single shower patch trigger class
   * @param[in] cls Type of the trigger class to check
   * @return True if the trigger class is a single shower patch trigger class
   */
  static Bool_t IsSingleShower(EmcalTriggerClass cls);

  /**
   * @brief Checks if the trigger class is a jet patch trigger class
   * @param[in] cls Type of the trigger class to check
   * @return True if the trigger class is a single shower patch trigger class
   */
  static Bool_t IsDCAL(EmcalTriggerClass cls);

  /**
   * @brief Get the name of the trigger class as string representation
   * @param[in] cls EMCAL/DCAL trigger class
   * @return Name of the trigger class
   */
  static const TString &GetTriggerName(EmcalTriggerClass cls) { return fgkTriggerNames[cls]; }

  /**
   * @brief Set the energy resolution used to smear the threshold.
   *
   * Note that an absolute value for the energy is used. In case the
   * resolution is set to 0, no smearing is applied.
   *
   * @param[in] resolution Energy resolution used for smearing
   */
  void SetEnergyResolution(Double_t resolution) { fResolution = resolution; }

  /**
   * @brief Switch whether to use original or smeared energy.
   * @param[in] doUse If true then the smeared energy is used
   */
  void SetUseSmearedEnergy(Bool_t doUse) { fUseSmearedEnergy = doUse; }


protected:

  /**
   * @brief Run event selection using trigger patches.
   *
   * Events are selected in case at least one trigger patch is found
   * with an energy above threshold. The energy can be of different
   * types (energy, transverse energy, ADC, transverse ADC). Note
   * that triggers with transverse energy definitions are hypothetical
   * triggers, used for principal studies, and are not implemented in
   * the hardware.
   *
   * Trigger patches must match the correct type. This means that for
   * L0 or gamma triggers only gamma patches are used, while for
   * jet triggers only jet patches are used.
   *
   * Effects of the energy resolution are studied using randomized
   * thresholds.
   *
   * @param[in] trgcls EMCAL trigger class to be triggered
   * @param[in] triggerpatches Container of trigger patches used to trigger the event
   * @return True if the event is selected, false otherwise
   */
  bool ApplyPatchTrigger(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const;

  /**
   * @brief Run event selection using a EMCAL clusters
   *
   * Events are selected requiring at least one EMCAL cluster with
   * energy above threshold. The energy can be either the total
   * cluster energy ot the cluster transverse energy. Note
   * that this trigger is a hypothetical trigger, used for
   * principal studies, and is not implemented as such in the
   * hardware.
   *
   * @param[in] trgcls EMCAL trigger classs
   * @param[in] event Input event, containing the cluster container
   * @return True if the event is selected, false otherwise
   */
  bool ApplyClusterTrigger(EmcalTriggerClass trgcls, const AliVEvent * const event) const;

  /**
   * @brief Check whether the trigger observable is based on clusters
   * @return True if the trigger is based on clusters, false otherwise
   */
  bool UseClusters() const;

  /**
   * @brief Check whether the trigger observable is based on trigger patches
   * @return True if the observable is based on trigger patches
   */
  bool UsePatches() const;

  static const TString        fgkTriggerNames[kTrgn];                     ///< Names of the various trigger classes
  Double_t                    fOfflineEnergyThreshold[kTrgn];             ///< Thresholds applied on offline energy
  const TH2                   *fAcceptanceMaps[kTrgn];                    //!<! Online acceptance distribution
  EmcalEnergyDefinition_t     fEnergyDefinition;                          ///< Define type of energy to be use for the patch selection
  TString                     fNameClusterContainer;                      ///< Name of the cluster container
  Double_t                    fResolution;                                ///< Resolution for threshold smearing
  Bool_t                      fUseSmearedEnergy;                          ///< Switch whether to use smeared or original energy

  ClassDef(AliEmcalTriggerOfflineSelection, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGEROFFLINESELECTION_H */
