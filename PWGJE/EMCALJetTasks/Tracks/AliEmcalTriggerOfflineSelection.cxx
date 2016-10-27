/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <TH2.h>
#include <TRandom.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliEmcalTriggerOfflineSelection)
/// \endcond

namespace EMCalTriggerPtAnalysis {

const TString AliEmcalTriggerOfflineSelection::fgkTriggerNames[AliEmcalTriggerOfflineSelection::kTrgn] = {
    "EMC7", "EG1", "EG2", "EJ1", "EJ2", "DMC7", "DG1", "DG2", "DJ1", "DJ2"
};

AliEmcalTriggerOfflineSelection::AliEmcalTriggerOfflineSelection() {
  for(int itrg = 0; itrg < kTrgn; itrg++) fOfflineEnergyThreshold[itrg] = 100000.;  // unimplemented triggers get very high threshold assinged, so that the result is automatically false
  memset(fAcceptanceMaps, 0, sizeof(TH2 *) * kTrgn);
}

AliEmcalTriggerOfflineSelection::~AliEmcalTriggerOfflineSelection(){
  for(int itrg = 0; itrg < kTrgn; itrg++){
    if(fAcceptanceMaps[itrg]) delete fAcceptanceMaps[itrg];
  }
}

Bool_t AliEmcalTriggerOfflineSelection::IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const {
  if(fOfflineEnergyThreshold[trgcls] < 0) return true;
  AliDebugStream(1) << "Applying offline threshold " << fOfflineEnergyThreshold[trgcls] << " for trigger class " << GetTriggerName(trgcls) << std::endl;
  bool isSingleShower = IsSingleShower(trgcls), isDCAL = IsDCAL(trgcls);
  int nfound = 0;
  AliEMCALTriggerPatchInfo *patch = NULL;
  for(auto patchIter : *triggerpatches){
    patch = static_cast<AliEMCALTriggerPatchInfo *>(patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if((isDCAL && !patch->IsDCalPHOS()) || (!isDCAL && patch->IsDCalPHOS())) continue;      // reject patches in opposite detector
    if(isSingleShower){
     if(!patch->IsGammaLowSimple()) continue;
    } else {
      if(!patch->IsJetLowSimple()) continue;
    }
    if(patch->GetPatchE() > fOfflineEnergyThreshold[trgcls]){
      AliDebugStream(2) << GetTriggerName(trgcls) << " patch above threshold (" << patch->GetPatchE() << " | " << fOfflineEnergyThreshold[trgcls] << ")" <<  std::endl;
      // Handle offline acceptance maps (if providede)
      // sampling means probe (random value) has to be
      if(fAcceptanceMaps[trgcls]){
        double acceptanceprob = fAcceptanceMaps[trgcls]->GetBinContent(patch->GetColStart(), patch->GetRowStart()),
                patchprob = gRandom->Uniform(0., 1.);
        AliDebugStream(2) << "Sampling trigger " << GetTriggerName(trgcls) << ": Efficiency(" << acceptanceprob << "), sampled (" << patchprob << ")" << std::endl;
        if(patchprob < acceptanceprob){
          AliDebugStream(2) << "Patch selected" << std::endl;
          nfound++;
        }
      } else{
        nfound++;
      }
    }
  }
  if(nfound){
    AliDebugStream(1) << "Event selected for trigger class " << GetTriggerName(trgcls) << ", " << nfound << " good patch(es) found" << std::endl;
    return true;
  }
  return false;
}

Bool_t AliEmcalTriggerOfflineSelection::IsSingleShower(EmcalTriggerClass cls){
 return ((cls == kTrgEG1) || (cls == kTrgEG2) || (cls == kTrgEL0) || (cls == kTrgDG1) || (cls == kTrgDG2) || (cls == kTrgDL0));

}

Bool_t AliEmcalTriggerOfflineSelection::IsDCAL(EmcalTriggerClass cls){
  return ((cls == kTrgDL0) || (cls == kTrgDG1) || (cls == kTrgDG2) || (cls == kTrgDJ1) || (cls == kTrgDJ2));
}

} /* namespace EMCalTriggerPtAnalysis */
