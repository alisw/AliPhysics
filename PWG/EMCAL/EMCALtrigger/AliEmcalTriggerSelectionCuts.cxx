/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/*
 * Class for the selection of trigger patches in the EMCAL triggered event selection
 *
 * Author: Markus Fasel
 */
#include "AliEmcalTriggerSelectionCuts.h"
#include "AliEMCALTriggerPatchInfo.h"

ClassImp(AliEmcalTriggerSelectionCuts)

//______________________________________________________________________________
AliEmcalTriggerSelectionCuts::AliEmcalTriggerSelectionCuts() :
  TObject(),
  fSelectionMethod(kADC),
  fPatchType(kAnyPatch),
  fThreshold(0),
  fUseSimpleOffline(kFALSE)
{
  /*
   * Dummy constructor
   */
}

//______________________________________________________________________________
Bool_t AliEmcalTriggerSelectionCuts::IsSelected(const AliEMCALTriggerPatchInfo * const patch) const {
  /*
   * Apply selection of the given trigger patch according to the selections described in the object
   *
   * @param patch: the trigger patch to check
   * @return" the decision (true if selected, false otherwise)
   */
  if(fUseSimpleOffline && !patch->IsOfflineSimple()) return kFALSE;
  else if(!fUseSimpleOffline && patch->IsOfflineSimple()) return kFALSE;
  if(!SelectPatchType(patch)) return kFALSE;
  if(GetCutPrimitive(patch) <= fThreshold) return kFALSE;
  return kTRUE;
}

//______________________________________________________________________________
Int_t AliEmcalTriggerSelectionCuts::CompareTriggerPatches(const AliEMCALTriggerPatchInfo *first, const AliEMCALTriggerPatchInfo *second) const {
  /*
   * Compare two patches according to the energy measure specified in the cut object
   *
   * @param first: the first patch
   * @param second: the second patch
   * @return: the result of the comparison (0 if equal, 1 if the first patch has a larger primitive,
   *          -1 if the second patch has a larger primitive)
   */
  Double_t valfirst = GetCutPrimitive(first), valsecond = GetCutPrimitive(second);
  if(valfirst == valsecond) return 0;
  if(valfirst > valsecond) return 1;
  return -1;
}

//______________________________________________________________________________
Double_t AliEmcalTriggerSelectionCuts::GetCutPrimitive(const AliEMCALTriggerPatchInfo * const patch) const{
  /*
   * Return (energy) measure we cut on, depending on the selection method specified
   *
   * @param patch: The patch from which to obtain the value
   * @return: The energy measure of the patch
   */
  if(fSelectionMethod == kADC) return static_cast<Double_t>(patch->GetADCAmp());
  else if(fSelectionMethod == kEnergyRough) return patch->GetADCAmpGeVRough();
  return patch->GetPatchE();
}

//______________________________________________________________________________
Bool_t AliEmcalTriggerSelectionCuts::SelectPatchType(const AliEMCALTriggerPatchInfo * const patch) const{
  /*
   * Select type of the patch according the definitions in the header file
   *
   * @param patch: the patch to be checked
   * @return: selection result (true ig the patch is selected)
   */
  if(fPatchType == kAnyPatch) return kTRUE;
  if(fUseSimpleOffline){
    if(patch->IsJetLowSimple() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetLowPatch))) return kTRUE;
    if(patch->IsJetHighSimple() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetHighPatch))) return kTRUE;
    if(patch->IsGammaLowSimple() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaLowPatch))) return kTRUE;
    if(patch->IsGammaHighSimple() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaHighPatch))) return kTRUE;
  } else {
    if(patch->IsJetLow() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetLowPatch))) return kTRUE;
    if(patch->IsJetHigh() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetHighPatch))) return kTRUE;
    if(patch->IsGammaLow() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaLowPatch))) return kTRUE;
    if(patch->IsGammaHigh() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaHighPatch))) return kTRUE;
    if(patch->IsLevel0() && fPatchType == kL0Patch) return kTRUE;
  }
  return kFALSE;
}
