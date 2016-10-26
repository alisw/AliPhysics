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
#include <TClonesArray.h>
#include <TH2.h>
#include <TRandom.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliEmcalTriggerOfflineSelection)
/// \endcond

namespace EMCalTriggerPtAnalysis {

AliEmcalTriggerOfflineSelection::AliEmcalTriggerOfflineSelection() {
  memset(fOfflineEnergyThreshold, 0, sizeof(Double_t) * kTrgn);
  memset(fAcceptanceMaps, 0, sizeof(TH2 *) * kTrgn);
}

AliEmcalTriggerOfflineSelection::~AliEmcalTriggerOfflineSelection(){
  for(int itrg = 0; itrg < kTrgn; itrg++){
    if(fAcceptanceMaps[itrg]) delete fAcceptanceMaps[itrg];
  }
}

Bool_t AliEmcalTriggerOfflineSelection::IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const {
  if(fOfflineEnergyThreshold[trgcls] < 0) return true;
  bool isSingleShower = IsSingleShower(trgcls);
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
    if(patch->GetPatchE() > fOfflineEnergyThreshold[trgcls]){
      // Handle offline acceptance maps (if providede)
      if(fAcceptanceMaps[trgcls]){
        double acceptanceprob = fAcceptanceMaps[trgcls]->GetBinContent(patch->GetColStart(), patch->GetRowStart()),
                patchprob = gRandom->Uniform(0., 1.);
        if(patchprob < acceptanceprob) nfound++;
      } else{
        nfound++;
      }
    }
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
