/**
 * \file AliEMCalTriggerAnaTriggerDecision.h
 * \brief Declaration of class AliEMCalTriggerAnaTriggerDecision, a container for trigger decision
 * in EMCAL-triggered events using several methods. The configuration of the trigger decision handler
 * is outsourced to AliEMCalTriggerAnaTriggerDecisionConfig, declared as well in this header in order to
 * avoid inclusion problem.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef ALIEMCALTRIGGERANATRIGGERDECISION_H
#define ALIEMCALTRIGGERANATRIGGERDECISION_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TClonesArray;
class TString;
class AliEMCALTriggerPatchInfo;

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-p_{t} tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-p_{t} tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData;
class AliEMCalTriggerDecisionConfig;

/**
 * \enum ETATriggerType
 * \brief Trigger types defined for this analysis
 *
 * This enumeration declared switches for then trigger types defined for
 * and used in this analysis.
 */
enum ETATriggerType{
  kTAEMCJHigh       = 0,            ///< Jet trigger, high threshold
  kTAEMCJLow        = 1,            ///< Jet trigger, low threshold
  kTAEMCGHigh       = 2,            ///< Gamma trigger, high threshold
  kTAEMCGLow        = 3,            ///< Gamma trigger, low threshold
  kTAUndef          = -1            ///< Default, undefined
};

/**
 * \enum EPatchEnergyType
 * Definition of the method to obtain the energy for an additional energy cut in the
 * selection of trigger patches.
 */
enum EPatchEnergyType_t{
  kAmplitudeOnline,                 ///< Online amplitude of the patch, from L0 time sums
  kAmplitudeOffline,                ///< Offline amplitude, estimated from EMCAL cells
  kEnergyOnline,                    ///< Online energy, estimated from L0 time sums
  kEnergyOffline                    ///< Offline energy, from cells, calibrated, exluding hot towers
};

/**
 * \enum ETriggerMethod_t
 * \brief Methods available to select event as triggered events
 *
 * This enumeration defines the possible methods to select events
 * as triggered events using the information stored in the trigger
 * decision object.
 */
enum ETriggerMethod_t{
  kTriggerString,                   ///< kTriggerString
  kTriggerPatches,                  ///< kPatches
  kTriggerMixed                     ///< kMixed
};

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

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerAnaTriggerDecisionConfig, 1);
  /// \endcond
};

/**
 * \class AliEMCalTriggerAnaTriggerDecision
 * \brief Class performing the selection of triggered events
 *
 * This class generates a trigger decision for an event for the given EMCAL trigger class.
 * The decision can come either from trigger patches or trigger strings. For the trigger decision
 * from patches also an energy threshold can be applied. The trigger decision is done in the
 * function Create, which builds the decision for both patches and trigger string. The event selection,
 * done using the function IsTriggerd, can come from the patches, from the string, or from both at the
 * same time. Note that different selection modes need to run in data and Monte-Carlo since in Monte-
 * Carlo trigger strings do not exists, and in data the trigger is already applied at hardware level.
 */
class AliEMCalTriggerAnaTriggerDecision : public TObject {
public:


  AliEMCalTriggerAnaTriggerDecision();

  /**
   * Destructor, nothing to do
   */
  virtual ~AliEMCalTriggerAnaTriggerDecision(){}

  void Create(const AliEMCalTriggerEventData * const data);

  /**
   * Apply event selection using a given trigger class. Three methods are avialable:
   *  -# Using the trigger string (kTriggerString)
   *  -# Using the reconstructed patches, online or offline (kTriggerPatches)
   *  -# Using patches and trigger string (kTriggerMixed): In this case both conditions need to be fulfilled
   *
   * \param trigger Type of the trigger requested in order to select the event
   * \param method Method (see above) used for the event selection
   * \return True if event is selected as triggered event, false otherwise
   */
  Bool_t IsTriggered(ETATriggerType trigger, ETriggerMethod_t method = kTriggerString) const {
    Bool_t result = kFALSE;
    switch(method){
    case kTriggerString: result = fDecisionFromString[trigger]; break;
    case kTriggerPatches: result = fDecisionFromPatches[trigger]; break;
    case kTriggerMixed: result = fDecisionFromPatches[trigger] && fDecisionFromString[trigger]; break;
    };
    return result;
  }

  /**
   * Get the trigger decision configuration
   * \return the trigger decision configuration
   */
  const AliEMCalTriggerAnaTriggerDecisionConfig * GetConfiguration() const { return &fConfiguration; }

  /**
   * Set the configuration used for the trigger selection.
   * \param conf Configuration used for the trigger selection
   */
  void ConfigureTriggerDecision(const AliEMCalTriggerAnaTriggerDecisionConfig &conf) { fConfiguration = conf; }

  void Reset();

  void Print(Option_t * opt = NULL) const;

  /**
   * Set the analysis into debug mode (enabling debug printouts)
   * \param doDebug If true we enable debug printouts.
   */
  void SetDebugMode(Bool_t doDebug = true) { fDoDebug = doDebug; }

  bool CheckConsistency() const;

protected:
  void MakeDecisionFromString(const TString &triggerstring);
  void MakeDecisionFromPatches(const TClonesArray &listOfPatches);

  Bool_t SelectTriggerPatch(ETATriggerType trigger, const AliEMCALTriggerPatchInfo * const recpatch) const;
  Double_t GetPatchEnergy(EPatchEnergyType_t energytype,  const AliEMCALTriggerPatchInfo *const patch) const;

  AliEMCalTriggerAnaTriggerDecisionConfig     fConfiguration;                      ///< Configuration for the trigger decision handler
  Bool_t                                      fDecisionFromPatches[4];             ///< Storage for result from trigger string
  Bool_t                                      fDecisionFromString[4];              ///< Storage for result from trigger patches

  Bool_t fDoDebug;                                                                 ///< Switch for debug mode

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerAnaTriggerDecision, 1);     // EMCal trigger decision
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERANATRIGGERDECISION_H */
