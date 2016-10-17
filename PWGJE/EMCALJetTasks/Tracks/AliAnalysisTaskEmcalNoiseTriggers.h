#ifndef ALIANALYSISTASKEMCALNOISETRIGGERS_H
#define ALIANALYSISTASKEMCALNOISETRIGGERS_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliEMCALTriggerDataGrid.h"
#include <TString.h>

class AliEMCALTriggerPatchADCInfoAP;
class AliEMCALTriggerPatchInfo;

namespace EMCalTriggerPtAnalysis {

/**
 * @class AliAnalysisTaskEmcalNoiseTriggers
 * @brief Analysis of trigger quantities in rejected (noise) events
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Oct 17, 2016
 *
 * Tool to study trigger related quatities (FastORs and trigger patches) in events
 * rejected by the online rejection method (require at least one patch after FastOR
 * masking with energy above threshold). This allows a more detailed investigation
 * how single noisy FastORs contribute to a trigger decision both for the jet and
 * the gamma trigger.
 *
 * Important: Trigger maker required, and bad fastor map needs to applied in the
 * trigger maker
 */
class AliAnalysisTaskEmcalNoiseTriggers : public AliAnalysisTaskEmcalTriggerBase {
public:
  enum SelectPatchType_t{
      kOnline = 0,
      kRecalc = 1
  };

  /**
   * Dummy I/O constructor
   */
  AliAnalysisTaskEmcalNoiseTriggers();

  /**
   * Named constructor, initializing also output list
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalNoiseTriggers(const char *name);

  /**
   * Destructor, nothing dynamically allocated within the task
   */
  virtual ~AliAnalysisTaskEmcalNoiseTriggers() {}

  /**
   * Further event selection according to trigger bits and trigger string
   * @param[in] triggerbits Trigger bit selection
   * @param[in] triggerstring Trigger string selection (optional)
   */
  void SetSelectTrigger(ULong_t triggerbits, const TString &triggerstring = "") {
    fTriggerBits = triggerbits;
    fTriggerString = triggerstring;
  }

protected:

  /**
   * Implementation of framework function IsUserEventSelected: Allows only one trigger
   * class to be selected.
   * @return True if the event is in the required trigger class (and event is selected)
   */
  virtual bool IsUserEventSelected();

  /**
   * Actual analysis: Selecting patches connected to the given trigger class (only gamma
   * patches for the gamma trigger and jet patches for the jet trigger), and filling
   * histograms at FastOR and patch level. Histograms are described in the relevant
   * functions.
   * @return Always true
   */
  virtual bool Run();

  /**
   * Creating user histograms, on fastor and on patch level
   */
  virtual void CreateUserHistos();

  /**
   * Creating additional user objects, optional (not implemented for this task)
   */
  virtual void CreateUserObjects() {}

  /**
   * Implementation of framework function: Fill simple histograms (event counter,
   * z-vertex distribution)
   */
  virtual void UserFillHistosAfterEventSelection();

  /**
   * Fill L1 Fastor ADC data grid with ADC values. To be called at the beginning
   * of a new event, once it is selected.
   */
  void PrepareL1FastorADC();

  /**
   * Create FastOR energy grid corresponding to the reconstructed trigger patch.
   * Allocates memory which needs to be handled by the user.
   * @param[in] patch Input trigger patch (for which fastor ADC values are calculated)
   * @return Data grid with single fastor ADC values contributing to the trigger patch
   */
  AliEMCALTriggerPatchADCInfoAP *MakeFastorADCValuesForPatch(const AliEMCALTriggerPatchInfo &patch) const;

  /**
   * Select trigger patch according to the sum of the fastor ADC values. Patches are
   * selected without masking. Energy cuts are for the moment hard coded to the 2013
   * pPb thresholds.
   * @param[in] sumADC Sum of the fastor ADC values within a trigger patch
   * @return True if the patch is selected as firing the trigger
   */
  Bool_t SelectFiredPatch(Int_t sumADC) const;

  /**
   * Filling histograms on FastOR Level
   * - ADC spectrum of each FastOR
   */
  void AnalyseFastors();

  /**
   * Filling historgrams on trigger patch level
   * - Sum good vs sum all
   * - Max fastor ADC vs. sum all
   * - Number of contributing (non-zero) fastors vs. fractional ADC value of the highest FastOR
   * @param[in] recpatch Patch to be analysed
   * @param[in] pt Type of the trigger patch (online or recalc)
   */
  void AnalyseTriggerPatch(const AliEMCALTriggerPatchInfo &recpatch, SelectPatchType_t pt);

  static const TString                             fgkPatchNames[2];  ///< Names of the trigger patch types (for histograms)
  AliEMCALTriggerDataGrid<Int_t>                   fL1ADC;            ///< Level1 fastor ADCs
  ULong_t                                          fTriggerBits;      ///< Trigger bits for the event selection
  TString                                          fTriggerString;    ///< Trigger string for event selection

private:
  AliAnalysisTaskEmcalNoiseTriggers(const AliAnalysisTaskEmcalNoiseTriggers &);
  AliAnalysisTaskEmcalNoiseTriggers &operator=(const AliAnalysisTaskEmcalNoiseTriggers &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalNoiseTriggers, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALNOISETRIGGERS_H */
