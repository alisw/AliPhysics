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

#include <iosfwd>
#include <TObject.h>

class AliEMCALTriggerPatchInfo;

// operator<< has to be forward declared carefully to stay in the global namespace so that it works with CINT.
// For generally how to keep the operator in the global namespace, See: https://stackoverflow.com/a/38801633
namespace PWG { namespace EMCAL { class AliEmcalTriggerSelectionCuts; } }
std::ostream & operator<< (std::ostream &in, const PWG::EMCAL::AliEmcalTriggerSelectionCuts &cuts);

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
   * @enum RhoMethod_t
   * @brief Method based on which Rho is calculated
   */
  enum RhoMethod_t {
    kOnlineRho = 0,         ///< Online Rho (from STU, in ESD/AOD), in ADC counts
    kRecalcRho = 1,         ///< Recalc Rho (recalculated from background recalc patches), in ADC counts
    kOfflineRho = 2,        ///< Offline Rho (recalculated from simple offline patches), in Energy
    kNoRho = 3              ///< No rho
  };

  /**
   * @struct RhoForTrigger
   * @brief Rho estimate, subtracted from trigger patches
   * 
   * Rho values are mentioned for the detector it is subtracted from.
   * The values are estimated with the opposite side detector.
   */
  struct RhoForTrigger {
    double fRhoForEmcalOnline;      ///< Rho estimated with DCAL (from online data), subtracted from EMCAL patches
    double fRhoForEmcalRecalc;      ///< Rho estimated with DCAL (recalculated from background patches), subtracted from EMCAL patches
    double fRhoForEmcalOffline;     ///< Rho estimated with DCAL (recalculated from offline simple patches), subtracted from EMCAL patches
    double fRhoForDCALOnline;       ///< Rho estimated with EMCAL (from online data), subtracted from DCAL patches
    double fRhoForDCALRecalc;       ///< Rho estimated with EMCAL (recalculated from background patches), subtracted from DCAL patches
    double fRhoForDCALOffline;      ///< Rho estimated with EMCAL (recalculated from offline patches), subtracted from DCAL patches
  };

  /**
   * @brief Dummy constructor
   */
  AliEmcalTriggerSelectionCuts();

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalTriggerSelectionCuts() {}

  PatchType_t GetPatchType() const { return fPatchType; }
  SelectionMethod_t GetSelectionMethod() const { return fSelectionMethod; }
  Double_t GetThreshold() const { return fThreshold; }
  Bool_t IsRequestingSimpleOfflinePatches() const { return fUseSimpleOffline; }

  /**
   * @brief Check whether rho calculation is requested by the selection cuts
   * @return True if patch selection uses rho subtraction, false otherwise 
   */
  Bool_t IsSubtractRho() const { return fSubtractRho; }

  /**
   * @brief Get method for rho subtraction
   * @return Method for rho subtraction 
   */
  RhoMethod_t GetRhoMethod() const { return fRhoMethod; }

  void SetPatchType(PatchType_t patchType) { fPatchType = patchType; }
  void SetAcceptanceType(AcceptanceType_t acceptance) { fAcceptanceType = acceptance; }
  void SetSelectionMethod(SelectionMethod_t selectionMethod) { fSelectionMethod = selectionMethod; }
  void SetThreshold(Double_t threshold) { fThreshold = threshold; }
  void SetUseSimpleOfflinePatches(Bool_t doUse = kTRUE) { fUseSimpleOffline = doUse; }
  void SetUseRecalcPatches(Bool_t doUse = kTRUE) { fUseRecalc = doUse; }

  /**
   * @brief Enable/Disable background subtraction from the patch energy / ADC
   * @param doSubtract If true rho * A is subtracted from the patch energy ADC
   * @param method Method used to estimate rho
   * 
   * Only to be used for PbPb data / simulation
   */
  void SetSubtractRho(Bool_t doSubtract = kTRUE, RhoMethod_t method = kOnlineRho) { fSubtractRho = doSubtract; fRhoMethod = method; }

  /**
   * @brief Apply selection of the given trigger patch according to the selections described in the object
   *
   * @param[in] patch The trigger patch to check
   * @param[in] rhocontainer Container with rho values
   * @return the decision (true if selected, false otherwise)
   */
  Bool_t IsSelected(const AliEMCALTriggerPatchInfo * const patch, const RhoForTrigger &rhocontainer) const;

  /**
   * @brief Compare two patches according to the energy measure specified in the cut object
   *
   * @param[in] first the first patch
   * @param[in] second the second patch
   * @param[in] rhocontainer Container with rho values
   * @return the result of the comparison (0 if equal, 1 if the first patch has a larger primitive,
   *          -1 if the second patch has a larger primitive)
   */
  Int_t CompareTriggerPatches(const AliEMCALTriggerPatchInfo *first, const AliEMCALTriggerPatchInfo *second, const RhoForTrigger &rhocontainer) const;

  /**
   * @brief Output stream operator
   * 
   * Print cut settings to the stream.
   * 
   * @param stream Output stream used for logging
   * @param cuts Selection cuts object to be put on the stream
   * @return Stream after logging
   */
  friend std::ostream& ::operator<<(std::ostream &stream, const AliEmcalTriggerSelectionCuts &cuts);

protected:

  /**
   * @brief Return (energy) measure we cut on, depending on the selection method specified
   *
   * @param[in] patch The patch from which to obtain the value
   * @param[in] rho Event rho
   * @return The energy measure of the patch
   */
  Double_t GetCutPrimitive(const AliEMCALTriggerPatchInfo * const patch, double rho) const;

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

  /**
   * @brief Get the rho value for the detector according to the selection criterion
   * @param rhocontaienr Container with rho values
   * @param isEMCAL If true the rho is obtained for EMCAL, otherwise for DCAL
   * @return Rho value based on the selection criterion
   */
  Double_t  GetRho(const RhoForTrigger &rhocontaienr, Bool_t isEMCAL) const; 

  SelectionMethod_t     fSelectionMethod;           ///< Variable to cut on
  PatchType_t           fPatchType;                 ///< Type of the patch to be selected
  AcceptanceType_t      fAcceptanceType;            ///< Acceptance type (EMCAL or DCAL acceptance)
  Double_t              fThreshold;                 ///< Threshold used
  Bool_t                fUseSimpleOffline;          ///< Request simple offline patches
  Bool_t                fUseRecalc;                 ///< Request recalc patch
  Bool_t                fSubtractRho;               ///< Subtract rho from trigger patches
  RhoMethod_t           fRhoMethod;                 ///< Method to calculate rho

private:

  /**
   * @brief Helper function for output stream operator
   * 
   * Putting all cut settings to the output stream
   * 
   * @param stream Output stream used for logging
   */
  void PrintStream(std::ostream &stream) const;

  ClassDef(AliEmcalTriggerSelectionCuts, 1);         // Cuts for the EMCAL Trigger selection
};

}
}

#endif /* ALIEMCALTRIGGERSELECTIONCUTS_H */
