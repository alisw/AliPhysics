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
#include <iostream>
#include <string>
#include <map>
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
  fUseRecalc(kFALSE),
  fSubtractRho(kFALSE),
  fRhoMethod(kNoRho)
{
}

Bool_t AliEmcalTriggerSelectionCuts::IsSelected(const AliEMCALTriggerPatchInfo * const patch, const RhoForTrigger &rhocontainer) const {
  //if(fUseSimpleOffline && !patch->IsOfflineSimple()) return kFALSE;
  //else if(!fUseSimpleOffline && patch->IsOfflineSimple()) return kFALSE;
  if(!SelectAcceptance(patch)) return kFALSE;
  if(!SelectPatchType(patch)) return kFALSE;
  double rho = 0;
  if(fSubtractRho) rho = GetRho(rhocontainer, patch->IsEMCal());
  if(GetCutPrimitive(patch, rho) < fThreshold) return kFALSE;
  return kTRUE;
}

Int_t AliEmcalTriggerSelectionCuts::CompareTriggerPatches(const AliEMCALTriggerPatchInfo *first, const AliEMCALTriggerPatchInfo *second, const RhoForTrigger &rhocontainer) const {
  double rho = 0.;
  if(fSubtractRho) rho = GetRho(rhocontainer, first->IsEMCal());
  Double_t valfirst = GetCutPrimitive(first, rho), valsecond = GetCutPrimitive(second, rho);
  if(valfirst == valsecond) return 0;
  if(valfirst > valsecond) return 1;
  return -1;
}

Double_t AliEmcalTriggerSelectionCuts::GetCutPrimitive(const AliEMCALTriggerPatchInfo * const patch, double rho) const{
  const double kPatchSizeBackground = 8;     // background patches for rho estimate have a fixed patch size of 8
  double energy(0);
  switch(fSelectionMethod){
  case kADC: energy = static_cast<Double_t>(patch->GetADCAmp()); break;
  case kEnergyRough: energy = patch->GetADCAmpGeVRough(); break;
  case kEnergyOfflineSmeared: energy = patch->GetSmearedEnergy(); break;
  case kEnergyOffline: energy = patch->GetPatchE(); break;
  default: energy = -1.;
  };
  if(fSubtractRho) {
    // scale rho by the ratio of the patch sizes
    double scale1D = double(patch->GetPatchSize())/kPatchSizeBackground;
    energy -= rho * scale1D * scale1D; 
  }
  return energy;
}

Double_t AliEmcalTriggerSelectionCuts::GetRho(const AliEmcalTriggerSelectionCuts::RhoForTrigger &rhocontainer, Bool_t isEMCAL) const {
  double rho = 0;
  switch(fRhoMethod) {
    case kOnlineRho: rho = isEMCAL ? rhocontainer.fRhoForEmcalOnline : rhocontainer.fRhoForDCALOnline; break;
    case kRecalcRho: rho = isEMCAL ? rhocontainer.fRhoForEmcalRecalc : rhocontainer.fRhoForDCALRecalc; break;
    case kOfflineRho: rho = isEMCAL ? rhocontainer.fRhoForEmcalOffline : rhocontainer.fRhoForDCALOffline; break;
    case kNoRho:
    default: break;
  };
  return rho;
}

Bool_t AliEmcalTriggerSelectionCuts::SelectPatchType(const AliEMCALTriggerPatchInfo * const patch) const{
  if(fPatchType == kAnyPatch) return kTRUE;
  if(fUseSimpleOffline){
    if(patch->IsJetLowSimple() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetLowPatch))) return kTRUE;
    if(patch->IsJetHighSimple() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetHighPatch))) return kTRUE;
    if(patch->IsGammaLowSimple() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaLowPatch))) return kTRUE;
    if(patch->IsGammaHighSimple() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaHighPatch))) return kTRUE;
    if(patch->IsLevel0Simple() && fPatchType == kL0Patch) return kTRUE;
  } else if(fUseRecalc) {
    if(patch->IsJetLowRecalc() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetLowPatch))) return kTRUE;
    if(patch->IsJetHighRecalc() && ((fPatchType == kL1JetPatch) || (fPatchType == kL1JetHighPatch))) return kTRUE;
    if(patch->IsGammaLowRecalc() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaLowPatch))) return kTRUE;
    if(patch->IsGammaHighRecalc() && ((fPatchType == kL1GammaPatch) || (fPatchType == kL1GammaHighPatch))) return kTRUE;
    if(patch->IsLevel0Recalc() && fPatchType == kL0Patch) return kTRUE;
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

void AliEmcalTriggerSelectionCuts::PrintStream(std::ostream &stream) const {
  std::map<AcceptanceType_t, std::string> acceptancetext = {{kEMCALAcceptance, "EMCAL"},
                                                                      {kDCALAcceptance, "DCAL"},
                                                                      {kEMCALDCALAcceptance, "EMCAL-DCAL"}};
  std::map<PatchType_t, std::string> patchtypetext = {{kL0Patch, "Level0"},
                                                                {kL1GammaPatch, "L1-gamma"},
                                                                {kL1GammaHighPatch, "L1-gamma, high threshold"},
                                                                {kL1GammaLowPatch, "L1-gamma, low threshold"},
                                                                {kL1JetPatch, "L1-jet"},
                                                                {kL1JetHighPatch, "L1-jet, high threshold"},
                                                                {kL1JetLowPatch, "L1-jet, low threshold"}};
  std::map<SelectionMethod_t, std::string> selmodetext = {{kADC, "FastOR ADC"},
                                                                    {kEnergyRough, "FastOR Energy"},
                                                                    {kEnergyOffline, "FEE Energy"},
                                                                    {kEnergyOfflineSmeared, "FEE Energy, decalibrated"}};
  stream << "  Cut settings:" << std::endl;
  stream << "    acceptance:      " << acceptancetext.find(fAcceptanceType)->second << std::endl;
  stream << "    patchtype:       " << patchtypetext.find(fPatchType)->second << std::endl;
  stream << "    sel mode:        " << selmodetext.find(fSelectionMethod)->second << std::endl;
  stream << "    Offline Patches: " << (fUseSimpleOffline ? "yes" : "no") << std::endl;
  stream << "    Recalc Patches:  " << (fUseRecalc ? "yes" : "no") << std::endl;
  stream << "    Threshold:       " << fThreshold << std::endl;
}

}
}

std::ostream &operator<<(std::ostream &stream, const PWG::EMCAL::AliEmcalTriggerSelectionCuts &cuts){
  cuts.PrintStream(stream);
  return stream;
}
