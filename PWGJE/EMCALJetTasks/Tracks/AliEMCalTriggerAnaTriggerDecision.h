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

#include "AliEMCalTriggerAnaHelper.h"
#include "AliEMCalTriggerAnaTriggerDecisionConfig.h"

class TClonesArray;
class TString;
class AliEMCALTriggerPatchInfo;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerEventData;

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

  ClassDef(AliEMCalTriggerAnaTriggerDecision, 1);     // EMCal trigger decision
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERANATRIGGERDECISION_H */
