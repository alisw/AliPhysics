/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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
#include <iostream>
#include <TClonesArray.h>
#include <TString.h>
#include "AliVEvent.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerEventData.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerDecision)

using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy (I/O) and main constructor
 */
AliEMCalTriggerAnaTriggerDecision::AliEMCalTriggerAnaTriggerDecision() :
    fDoDebug(kFALSE)
{
  Reset();
}

/**
 * Steer creation of the trigger decision. Delegates the creation of the trigger
 * decision from trigger strings to the function MakeDecisionFromString and of
 * the decision from patches to the function MakeDecisionFromPatches.
 *
 * \param data all event information
 */
void AliEMCalTriggerAnaTriggerDecision::Create(const AliEMCalTriggerEventData* const data) {
  Reset();
  MakeDecisionFromPatches(*(data->GetTriggerPatchContainer()));
  MakeDecisionFromString(data->GetRecEvent()->GetFiredTriggerClasses());
}

/**
 * Reset trigger decisions stored in this object
 */
void AliEMCalTriggerAnaTriggerDecision::Reset() {
  for(int itrg = 0; itrg < 4; itrg++){
    fDecisionFromPatches[itrg] = kFALSE;
    fDecisionFromString[itrg] = kFALSE;
  }
}

/**
 * Create trigger decision from trigger string. For each trigger class we only check if the
 * name of the trigger class appears in the trigger string.
 *
 * \param triggerstring the trigger string stored in the reconstructed event.
 */
void AliEMCalTriggerAnaTriggerDecision::MakeDecisionFromString(const TString& triggerstring) {
  if(triggerstring.Contains("EJ1") || triggerstring.Contains("EJE")) fDecisionFromString[kTAEMCJHigh]  = kTRUE;
  if(triggerstring.Contains("EJ2")) fDecisionFromString[kTAEMCJLow]   = kTRUE;
  if(triggerstring.Contains("EG1") || triggerstring.Contains("EGA")) fDecisionFromString[kTAEMCGHigh]  = kTRUE;
  if(triggerstring.Contains("EG2")) fDecisionFromString[kTAEMCGLow]   = kTRUE;
}

/**
 * Check whether the trigger decision from the trigger strings and the trigger patches are the same
 *
 * \return result of the comparison
 */
bool AliEMCalTriggerAnaTriggerDecision::CheckConsistency() const{
  bool result = true;
  for(int icls = 0; icls < 4; icls++){
    if(fDecisionFromString[icls] != fDecisionFromPatches[icls]) result = false;
  }
  return result;
}

/**
 * Create trigger decision from trigger patches. In case swap thresholds is requested, the low threshold
 * triggers are replaced by the high threshold triggers and vice versa
 *
 * \param listOfPatches the TClonesArray of the trigger patches, created by the trigger patch maker
 */
void AliEMCalTriggerAnaTriggerDecision::MakeDecisionFromPatches(const TClonesArray& listOfPatches) {
  TIter patchIter(&listOfPatches);
  AliEMCALTriggerPatchInfo *mypatch(NULL);
  if(fDoDebug){
    std::cout << "Generating trigger decision from found patches: " << listOfPatches.GetEntries() << std::endl;
    if(fConfiguration.IsUsingOfflinePatches()) std::cout << "Using offline patches\n";
    else std::cout << "Using online patches\n";
  }
  int foundpatches[4] = {0,0,0,0};
  while((mypatch = dynamic_cast<AliEMCALTriggerPatchInfo *>(patchIter()))){
    if(fDoDebug) std::cout << "Next patch: " << (mypatch->IsOfflineSimple() ? "offline" : "online:") << std::endl;
	  for(int icase = 0; icase < 4; icase++){
	    if(SelectTriggerPatch(ETATriggerType(icase), mypatch)){
	      fDecisionFromPatches[icase] = kTRUE;
        foundpatches[icase]++;
	    }
	  }
  }
  if(fDoDebug){
    std::cout << "Found patches:" << std::endl;
    std::cout << "Jet high:    " << foundpatches[kTAEMCJHigh] << std::endl;
    std::cout << "Jet low:     " << foundpatches[kTAEMCJLow] << std::endl;
    std::cout << "Gamma high:  " << foundpatches[kTAEMCGHigh] << std::endl;
    std::cout << "Gamma low:   " << foundpatches[kTAEMCGLow] << std::endl;
  }
}

/**
 * Select trigger patch for a given trigger type:
 *  -# Check whether the patch is of patch type (online/offline) we require
 *  -# Check if the patch belongs to the given trigger class.
 *  -# If required perform additional selection of patches avbove energy.
 *
 * \param trigger Type of the trigger class
 * \param recpatch Patch to accept
 * \return True if the patch is selected, false otherwise
 */
Bool_t AliEMCalTriggerAnaTriggerDecision::SelectTriggerPatch(ETATriggerType trigger, const AliEMCALTriggerPatchInfo* const recpatch) const {
  bool swapThresholds = fConfiguration.IsSwapThresholds();
	bool selectPatchType = kFALSE;
	if(fConfiguration.IsUsingOfflinePatches()){
	  if(!recpatch->IsOfflineSimple()) return kFALSE;
		switch(trigger){
		case kTAEMCJHigh: selectPatchType = swapThresholds ? recpatch->IsJetLowSimple() :  recpatch->IsJetHighSimple(); break;
		case kTAEMCJLow: selectPatchType = swapThresholds ? recpatch->IsJetHighSimple() :  recpatch->IsJetLowSimple(); break;
		case kTAEMCGHigh: selectPatchType = swapThresholds ? recpatch->IsGammaLowSimple() :  recpatch->IsGammaHighSimple(); break;
		case kTAEMCGLow: selectPatchType = swapThresholds ? recpatch->IsGammaHighSimple() :  recpatch->IsGammaLowSimple(); break;
		case kTAUndef: break;
		};
	} else {
	  if(recpatch->IsOfflineSimple()) return kFALSE;
		switch(trigger){
		case kTAEMCJHigh: selectPatchType = swapThresholds ? recpatch->IsJetLow() :  recpatch->IsJetHigh(); break;
		case kTAEMCJLow: selectPatchType = swapThresholds ? recpatch->IsJetHigh() :  recpatch->IsJetLow(); break;
		case kTAEMCGHigh: selectPatchType = swapThresholds ? recpatch->IsGammaLow() :  recpatch->IsGammaHigh(); break;
		case kTAEMCGLow: selectPatchType = swapThresholds ? recpatch->IsGammaHigh() :  recpatch->IsGammaLow(); break;
    case kTAUndef: break;
		};
	}

	if(!selectPatchType) return kFALSE;

	if(fConfiguration.HasEnergyThreshold(trigger)){
		// Additional threshold on energy requested to select the patch
		return GetPatchEnergy(fConfiguration.GetPatchEnergyType(), recpatch) > fConfiguration.GetEnergyThreshold(trigger);
	}
	return kTRUE;
}
/**
 * Retrieve patch energy from the reconstruced patch according to the energy definition specified.
 * \param energytype Energy type
 * \param patch The reconstructed patch
 * \return The patch energy
 */
Double_t AliEMCalTriggerAnaTriggerDecision::GetPatchEnergy(EPatchEnergyType_t energytype, const AliEMCALTriggerPatchInfo *const patch) const {
  Double_t energy = 0.;
  switch(energytype){
  case kAmplitudeOnline:      energy = patch->GetADCAmp(); break;
  case kAmplitudeOffline:     energy = patch->GetADCOfflineAmp(); break;
  case kEnergyOnline:         energy = patch->GetADCAmpGeVRough(); break;
  case kEnergyOffline:        energy = patch->GetPatchE(); break;
  };
  return energy;
}

/**
 * Print status of the trigger decision
 *
 * \param Parameter required by the interface, not used here
 */
void AliEMCalTriggerAnaTriggerDecision::Print(Option_t*) const {
  std::cout << "Trigger decision" << std::endl;
  std::cout << "===============================" << std::endl;
  std::string triggertitles[4] = {"Jet High", "Jet Low", "Gamma High", "Gamma Low"};
  for(int icase = 0; icase < 4; icase++){
    std::cout << triggertitles[icase] << ": String[" << (fDecisionFromString[icase] ? "yes" : "no")
        << "], Patches[" << (fDecisionFromPatches[icase] ? "yes" : "no") << "]" << std::endl;
  }
}
