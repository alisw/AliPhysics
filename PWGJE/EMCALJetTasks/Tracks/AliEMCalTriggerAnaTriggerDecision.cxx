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
/*
 * Class performing the selection of triggered events
 *
 * Author:
 *    Markus Fasel
 */
#include <iostream>
#include <TClonesArray.h>
#include <TString.h>
#include "AliVEvent.h"
#include "AliEmcalTriggerPatchInfo.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerEventData.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerAnaTriggerDecision)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerAnaTriggerDecision::AliEMCalTriggerAnaTriggerDecision() :
    fSwapThresholds(kFALSE),
    fIsMinBias(kFALSE),
	fUseOfflinePatches(kFALSE),
    fDoDebug(kFALSE)
{
  /*
   * Main constructor
   */
  Reset();
  memset(fEnergyThresholds, 0, sizeof(double) * 4);
}

//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecision::Create(const AliEMCalTriggerEventData* const data) {
  /*
   * Steer creation of the trigger decision
   *
   * @param data: all event information
   */
  Reset();
  MakeDecisionFromPatches(*(data->GetTriggerPatchContainer()));
  MakeDecisionFromString(data->GetRecEvent()->GetFiredTriggerClasses());
}

//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecision::Reset() {
  for(int itrg = 0; itrg < 4; itrg++){
    fDecisionFromPatches[itrg] = kFALSE;
    fDecisionFromString[itrg] = kFALSE;
  }
}

//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecision::MakeDecisionFromString(const TString& triggerstring) {
  /*
   * Create trigger decision from trigger string
   *
   * @param triggerstring: the trigger string
   */
  if(triggerstring.Contains("EJ1")) fDecisionFromString[kTAEMCJHigh]  = kTRUE;
  if(triggerstring.Contains("EJ2")) fDecisionFromString[kTAEMCJLow]   = kTRUE;
  if(triggerstring.Contains("EG1")) fDecisionFromString[kTAEMCGHigh]  = kTRUE;
  if(triggerstring.Contains("EG2")) fDecisionFromString[kTAEMCGLow]   = kTRUE;
}

//______________________________________________________________________________
bool AliEMCalTriggerAnaTriggerDecision::CheckConsistency() const{
  /*
   * Check whether the trigger decision from the trigger strings and the trigger patches are the same
   *
   * @return: result of the comparison
   */
  bool result = true;
  for(int icls = 0; icls < 4; icls++){
    if(fDecisionFromString[icls] != fDecisionFromPatches[icls]) result = false;
  }
  return result;
}

//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecision::MakeDecisionFromPatches(const TClonesArray& listOfPatches) {
  /*
   * Create trigger decision from trigger patches. In case swap thresholds is requested, the low threshold
   * triggers are replaced by the high threshold triggers and vice versa
   *
   * @param triggerstring: the TClonesArray of the trigger patches, created by the trigger patch maker
   */
  TIter patchIter(&listOfPatches);
  AliEmcalTriggerPatchInfo *mypatch(NULL);
  if(fDoDebug)
    std::cout << "Generating trigger decision from found patches: " << listOfPatches.GetEntries() << std::endl;
  int foundpatches[4] = {0,0,0,0};
  int index = -1;
  while((mypatch = dynamic_cast<AliEmcalTriggerPatchInfo *>(patchIter()))){
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

//______________________________________________________________________________
Bool_t AliEMCalTriggerAnaTriggerDecision::SelectTriggerPatch(ETATriggerType trigger, const AliEmcalTriggerPatchInfo* const recpatch) const {
	/*
	 *
	 */
	bool selectPatchType = kFALSE;
	if(fUseOfflinePatches){
		switch(trigger){
		case kTAEMCJHigh: selectPatchType = fSwapThresholds ? recpatch->IsJetLowSimple() :  recpatch->IsJetHighSimple(); break;
		case kTAEMCJLow: selectPatchType = fSwapThresholds ? recpatch->IsJetHighSimple() :  recpatch->IsJetLowSimple(); break;
		case kTAEMCGHigh: selectPatchType = fSwapThresholds ? recpatch->IsGammaLowSimple() :  recpatch->IsGammaHighSimple(); break;
		case kTAEMCGLow: selectPatchType = fSwapThresholds ? recpatch->IsGammaHighSimple() :  recpatch->IsGammaLowSimple(); break;
		};
	} else {
		switch(trigger){
		case kTAEMCJHigh: selectPatchType = fSwapThresholds ? recpatch->IsJetLow() :  recpatch->IsJetHigh(); break;
		case kTAEMCJLow: selectPatchType = fSwapThresholds ? recpatch->IsJetHigh() :  recpatch->IsJetLow(); break;
		case kTAEMCGHigh: selectPatchType = fSwapThresholds ? recpatch->IsGammaLow() :  recpatch->IsGammaHigh(); break;
		case kTAEMCGLow: selectPatchType = fSwapThresholds ? recpatch->IsGammaHigh() :  recpatch->IsGammaLow(); break;
		};
	}

	if(!selectPatchType) return kFALSE;

	if(fEnergyThresholds[trigger]){
		// Additional threshold on energy requested to select the patch
		return recpatch->GetPatchE() > fEnergyThresholds[trigger];
	}
	return kTRUE;
}


//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecision::Print(Option_t*) const {
  /*
   * Print status of the trigger decision
   */
  std::cout << "Trigger decision" << std::endl;
  std::cout << "===============================" << std::endl;
  std::cout << "MinBias:                   " << (fIsMinBias ? "yes" : "no") << std::endl;
  std::string triggertitles[4] = {"Jet High", "Jet Low", "Gamma High", "Gamma Low"};
  for(int icase = 0; icase < 4; icase++){
    std::cout << triggertitles[icase] << ": String[" << (fDecisionFromString[icase] ? "yes" : "no")
        << "], Patches[" << (fDecisionFromPatches[icase] ? "yes" : "no") << "]" << std::endl;
  }
}

} /* namespace EMCalTriggerPtAnalysis */

