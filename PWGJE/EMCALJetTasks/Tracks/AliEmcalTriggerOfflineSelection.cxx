/*
 * AliEmcalTriggerOfflineSelection.cxx
 *
 *  Created on: Feb 23, 2016
 *      Author: markus
 */
#include <TClonesArray.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"

ClassImp(EMCalTriggerPtAnalysis::AliEmcalTriggerOfflineSelection)

namespace EMCalTriggerPtAnalysis {

AliEmcalTriggerOfflineSelection::AliEmcalTriggerOfflineSelection() {
  enum EmcalTriggerClass{
    kCPREL0 = 0,
    kCPREG1,
    kCPREG2,
    kCPREJ1,
    kCPREJ2,
    kCPRntrig
  };
}

/**
 * Apply additional cut requiring at least one offline patch above a given energy (not fake ADC!)
 * Attention: This task groups into single shower triggers (L0, EG1, EG2) and jet triggers (EJ1 and EJ2).
 * Per convention the low threshold patch is selected. No energy cut should be applied in the trigger maker
 * @param trgcls Trigger class for which to apply additional offline patch selection
 * @param triggerpatches Array of trigger patches
 * @return True if at least on patch above threshold is found or no cut is applied
 */
Bool_t AliEmcalTriggerOfflineSelection::IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const {
  if(fOfflineEnergyThreshold[trgcls] < 0) return true;
  bool isSingleShower = ((trgcls == kTrgEL0) || (trgcls == kTrgEG1) || (trgcls == kTrgEG2));
  int nfound = 0;
  AliEMCALTriggerPatchInfo *patch = NULL;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    patch = static_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(isSingleShower){
     if(!patch->IsGammaLowSimple()) continue;
    } else {
      if(!patch->IsJetLowSimple()) continue;
    }
    if(patch->GetPatchE() > fOfflineEnergyThreshold[trgcls]) nfound++;
  }
  return nfound > 0;
}

Bool_t AliEmcalTriggerOfflineSelection::IsSingleShower(EmcalTriggerClass cls){
 return ((cls == kTrgEG1) || (cls == kTrgEG2) || (cls == kTrgEL0) || (cls == kTrgDG1) || (cls == kTrgDG2) || (cls == kTrgDL0));

}

Bool_t AliEmcalTriggerOfflineSelection::IsDCAL(EmcalTriggerClass cls){
  return ((cls == kTrgDL0) || (cls == kTrgDG1) || (cls == kTrgDG2) || (cls == kTrgDJ1) || (cls == kTrgDJ2));
}

} /* namespace EMCalTriggerPtAnalysis */
