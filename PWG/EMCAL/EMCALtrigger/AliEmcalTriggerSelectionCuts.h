/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIEMCALTRIGGERSELECTIONCUTS_H
#define ALIEMCALTRIGGERSELECTIONCUTS_H

#include <TObject.h>

class AliEMCALTriggerPatchInfo;

namespace PWG {
namespace EMCAL {

/**
 * @class AliEmcalTriggerSelectionCuts
 * @brief Class for the selection of trigger patches in the EMCAL triggered event selection
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch> Oak Ridge National Laboratory
 * @since Dec 17, 2014
 */
class AliEmcalTriggerSelectionCuts: public TObject {
public:
  enum SelectionMethod_t {
    kADC = 0,
    kEnergyRough = 1,
    kEnergyOffline = 2,
    kEnergyOfflineSmeared = 3
  };
  enum PatchType_t {
    kAnyPatch = 0,
    kL1JetPatch = 1,
    kL1GammaPatch = 2,
    kL0Patch = 3,
    kL1JetLowPatch = 4,
    kL1JetHighPatch = 5,
    kL1GammaLowPatch = 6,
    kL1GammaHighPatch = 7
  };
  enum AcceptanceType_t {
    kEMCALAcceptance = 0,
    kDCALAcceptance = 1,
    kEMCALDCALAcceptance = 2
  };

  /**
   * @brief Dummy constructor
   */
  AliEmcalTriggerSelectionCuts();
  virtual ~AliEmcalTriggerSelectionCuts() {}

  PatchType_t GetPatchType() const { return fPatchType; }
  SelectionMethod_t GetSelectionMethod() const { return fSelectionMethod; }
  Double_t GetThreshold() const { return fThreshold; }
  Bool_t IsRequestingSimpleOfflinePatches() const { return fUseSimpleOffline; }

  void SetPatchType(PatchType_t patchType) { fPatchType = patchType; }
  void SetAcceptanceType(AcceptanceType_t acceptance) { fAcceptanceType = acceptance; }
  void SetSelectionMethod(SelectionMethod_t selectionMethod) { fSelectionMethod = selectionMethod; }
  void SetThreshold(Double_t threshold) { fThreshold = threshold; }
  void SetUseSimpleOfflinePatches(Bool_t doUse = kTRUE) { fUseSimpleOffline = doUse; }
  void SetUseRecalcPatches(Bool_t doUse = kTRUE) { fUseRecalc = doUse; }

  /**
   * @brief Apply selection of the given trigger patch according to the selections described in the object
   *
   * @param[in] patch the trigger patch to check
   * @return the decision (true if selected, false otherwise)
   */
  Bool_t IsSelected(const AliEMCALTriggerPatchInfo * const patch) const;

  /**
   * @brief Compare two patches according to the energy measure specified in the cut object
   *
   * @param[in] first the first patch
   * @param[in] second the second patch
   * @return the result of the comparison (0 if equal, 1 if the first patch has a larger primitive,
   *          -1 if the second patch has a larger primitive)
   */
  Int_t CompareTriggerPatches(const AliEMCALTriggerPatchInfo *first, const AliEMCALTriggerPatchInfo *second) const;

protected:

  /**
   * @brief Return (energy) measure we cut on, depending on the selection method specified
   *
   * @param[in] patch The patch from which to obtain the value
   * @return The energy measure of the patch
   */
  Double_t GetCutPrimitive(const AliEMCALTriggerPatchInfo * const patch) const;

  /**
   * @brief Select type of the patch according the definitions in the header file
   *
   * @param[in] patch the patch to be checked
   * @return selection result (true ig the patch is selected)
   */
  Bool_t SelectPatchType(const AliEMCALTriggerPatchInfo * const patch) const;

  /**
   * @brief Select detector acceptance
   *
   * Detector acceptance can be either of
   * - EMCAL
   * - DCAL
   * - EMCAL and DCAL
   * @param patch Trigger patch to be checked
   * @return True if the patch is within the acceptance, false otherwise
   */
  Bool_t SelectAcceptance(const AliEMCALTriggerPatchInfo * const patch) const;

  SelectionMethod_t     fSelectionMethod;           ///< Variable to cut on
  PatchType_t           fPatchType;                 ///< Type of the patch to be selected
  AcceptanceType_t      fAcceptanceType;            ///< Acceptance type (EMCAL or DCAL acceptance)
  Double_t              fThreshold;                 ///< Threshold used
  Bool_t                fUseSimpleOffline;          ///< Request simple offline patches
  Bool_t                fUseRecalc;                 ///< Request recalc patch

  /// \cond CLASSIMP
  ClassDef(AliEmcalTriggerSelectionCuts, 1);         // Cuts for the EMCAL Trigger selection
  /// \endcond
};

}
}

#endif /* ALIEMCALTRIGGERSELECTIONCUTS_H */
