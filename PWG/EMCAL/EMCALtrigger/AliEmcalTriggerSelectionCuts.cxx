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
#include "AliEmcalTriggerSelectionCuts.h"
#include "AliEMCALTriggerPatchInfo.h"

ClassImp(PWG::EMCAL::AliEmcalTriggerSelectionCuts)

namespace PWG {
namespace EMCAL {


AliEmcalTriggerSelectionCuts::AliEmcalTriggerSelectionCuts() :
  TObject(),
  fSelectionMethod(kADC),
  fPatchType(kAnyPatch),
  fAcceptanceType(kEMCALDCALAcceptance),
  fThreshold(0),
  fUseSimpleOffline(kFALSE),
  fUseRecalc(kFALSE)
{
}

Bool_t AliEmcalTriggerSelectionCuts::IsSelected(const AliEMCALTriggerPatchInfo * const patch) const {
  //if(fUseSimpleOffline && !patch->IsOfflineSimple()) return kFALSE;
  //else if(!fUseSimpleOffline && patch->IsOfflineSimple()) return kFALSE;
  if(!SelectAcceptance(patch)) return kFALSE;
  if(!SelectPatchType(patch)) return kFALSE;
  if(GetCutPrimitive(patch) < fThreshold) return kFALSE;
  return kTRUE;
}

Int_t AliEmcalTriggerSelectionCuts::CompareTriggerPatches(const AliEMCALTriggerPatchInfo *first, const AliEMCALTriggerPatchInfo *second) const {
  Double_t valfirst = GetCutPrimitive(first), valsecond = GetCutPrimitive(second);
  if(valfirst == valsecond) return 0;
  if(valfirst > valsecond) return 1;
  return -1;
}

Double_t AliEmcalTriggerSelectionCuts::GetCutPrimitive(const AliEMCALTriggerPatchInfo * const patch) const{
  double energy(0);
  switch(fSelectionMethod){
  case kADC: energy = static_cast<Double_t>(patch->GetADCAmp()); break;
  case kEnergyRough: energy = patch->GetADCAmpGeVRough(); break;
  case kEnergyOfflineSmeared: energy = patch->GetSmearedEnergy(); break;
  case kEnergyOffline: energy = patch->GetPatchE(); break;
  default: energy = -1.;
  };
  return energy;
}

Bool_t AliEmcalTriggerSelectionCuts::SelectPatchType(const AliEMCALTriggerPatchInfo * const patch) const{
  if(fPatchType == kAnyPatch) return kTRUE;
  if(fUseSimpleOffline){
    if(patch->IsJetLowSimple() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetLowPatch))) return kTRUE;
    if(patch->IsJetHighSimple() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetHighPatch))) return kTRUE;
    if(patch->IsGammaLowSimple() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaLowPatch))) return kTRUE;
    if(patch->IsGammaHighSimple() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaHighPatch))) return kTRUE;
  } else if(fUseRecalc) {
    if(patch->IsJetLowRecalc() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetLowPatch))) return kTRUE;
    if(patch->IsJetHighRecalc() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetHighPatch))) return kTRUE;
    if(patch->IsGammaLowRecalc() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaLowPatch))) return kTRUE;
    if(patch->IsGammaHighRecalc() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaHighPatch))) return kTRUE;
  } else {
    if(patch->IsJetLow() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetLowPatch))) return kTRUE;
    if(patch->IsJetHigh() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetHighPatch))) return kTRUE;
    if(patch->IsGammaLow() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaLowPatch))) return kTRUE;
    if(patch->IsGammaHigh() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaHighPatch))) return kTRUE;
    if(patch->IsLevel0() && fPatchType == kL0Patch) return kTRUE;
  }
  return kFALSE;
}

Bool_t AliEmcalTriggerSelectionCuts::SelectAcceptance(const AliEMCALTriggerPatchInfo * const patch) const {
  Bool_t selected(false);
  switch(fAcceptanceType){
  case kEMCALAcceptance: selected = patch->IsEMCal(); break;
  case kDCALAcceptance: selected = patch->IsDCalPHOS(); break;
  case kEMCALDCALAcceptance: selected = patch->IsEMCal() || patch->IsDCalPHOS(); break;
  default: selected = false;
  };
  return selected;
}

}
}
