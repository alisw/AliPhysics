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
#include <algorithm>
#include <iostream>
#include <vector>

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
  std::vector<double> patchefficiencies;
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
      if(fAcceptanceMaps[trgcls]){
        // Handle azimuthal inefficiencies of the trigger observed online:
        // For each patch provide an efficiency. Once all patches in the event
        // is determined, only the patch with the maximum efficiency is chosen.
        // The event is selected in this case if the sample value is below the
        // efficiency for the chosen patch.
        double peff = fAcceptanceMaps[trgcls]->GetBinContent(patch->GetColStart(), patch->GetRowStart());
        patchefficiencies.push_back(peff);
        AliDebugStream(2) << "Spatial Efficiency " << peff
            << " for trigger patch at position (" << patch->GetColStart()
            << "," << patch->GetRowStart() << ")" << std::endl;
      } else{
        patchefficiencies.push_back(1.);
      }
    }
  }
  if(patchefficiencies.size()){
    std::sort(patchefficiencies.begin(), patchefficiencies.end(), std::greater<double>());
    double sample = gRandom->Uniform(0., 1.);
    if(sample < patchefficiencies[0]){
      AliDebugStream(1) << "Event selected for trigger class " << GetTriggerName(trgcls) << ", " << nfound << " good patch(es) found" << std::endl;
      return true;
    }
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
