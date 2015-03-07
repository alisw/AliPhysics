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
 * Analysis component for different trigger patches
 *
 *   Author: Markus Fasel
 */
#include <TClonesArray.h>

#include "AliEmcalTriggerPatchInfo.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerPatchAnalysisComponent.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerPatchAnalysisComponent)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerPatchAnalysisComponent::AliEMCalTriggerPatchAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent()
{
  /*
   * Dummy (I/O) constructor, not to be used
   */
}

//______________________________________________________________________________
AliEMCalTriggerPatchAnalysisComponent::AliEMCalTriggerPatchAnalysisComponent(const char *name) :
  AliEMCalTriggerTracksAnalysisComponent(name)
{
  /*
   * Main constructor, to be used by the users
   */
}

//______________________________________________________________________________
void AliEMCalTriggerPatchAnalysisComponent::CreateHistos() {
  /*
   * Create histograms for the trigger patch analysis
   */
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();

  AliEMCalTriggerBinningDimension *etabinning = fBinning->GetBinning("eta"),
      *phibinning = fBinning->GetBinning("phi");
  const TAxis *patchaxes[6] = {
    DefineAxis("energy", 100, 0., 100),
    DefineAxis("energyRough", 100, 0., 100),
    DefineAxis("amplitude", 5000, 0., 5000.),       // limit for the moment
    DefineAxis("eta", etabinning),
    DefineAxis("phi", phibinning),
    DefineAxis("isMain", 2, -0.5, 1.5)
  };

  std::string patchnames[] = {"Level0", "JetHigh", "JetLow", "GammaHigh", "GammaLow"};
  std::string triggermodes[] = {"Online", "Offline"};
  for(std::string * triggerpatch = patchnames; triggerpatch < patchnames + sizeof(patchnames)/sizeof(std::string); ++triggerpatch){
    for(std::string *triggermode = triggermodes; triggermode < triggermodes + sizeof(triggermodes)/sizeof(std::string); ++triggermode){
      if((!strcmp(triggermode->c_str(), "Offline")) && (!strcmp(triggerpatch->c_str(), "Level0"))) continue; // Don't process L0 in case of offline
      printf("Adding patch for trigger %s in case of %s\n", triggerpatch->c_str(), triggermode->c_str());
      fHistos->CreateTHnSparse(Form("PatchInfo%s%s", triggerpatch->c_str(), triggermode->c_str()), Form("Patch energy for %s %s trigger patches", triggerpatch->c_str(), triggermode->c_str()), 6, patchaxes, "s");
    }
  }

}

//______________________________________________________________________________
void AliEMCalTriggerPatchAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  /*
   * Run trigger patch analysis
   */
  AliEmcalTriggerPatchInfo *triggerpatch(NULL);
  TIter patchIter(data->GetTriggerPatchContainer());
  while((triggerpatch = dynamic_cast<AliEmcalTriggerPatchInfo *>(patchIter()))){
    double triggerpatchinfo[6] = {triggerpatch->GetPatchE(),triggerpatch->GetADCAmpGeVRough(),
        static_cast<double>(triggerpatch->GetADCAmp()), triggerpatch->GetEtaGeo(),
        triggerpatch->GetPhiGeo(), triggerpatch->IsMainTrigger() ? 1. : 0.};
    if(triggerpatch->IsOfflineSimple()){
      if(triggerpatch->IsJetHighSimple()){
        fHistos->FillTHnSparse("PatchInfoJetHighOffline", triggerpatchinfo);
      }
      if(triggerpatch->IsJetLowSimple()){
    	  fHistos->FillTHnSparse("PatchInfoJetLowOffline", triggerpatchinfo);
      }
      if(triggerpatch->IsGammaHighSimple()){
        fHistos->FillTHnSparse("PatchInfoGammaHighOffline", triggerpatchinfo);
      }
      if(triggerpatch->IsGammaLowSimple()){
        fHistos->FillTHnSparse("PatchInfoGammaLowOffline", triggerpatchinfo);
      }
    } else{
      if(triggerpatch->IsJetHigh()){
        fHistos->FillTHnSparse("PatchInfoJetHighOnline", triggerpatchinfo);
      }
      if(triggerpatch->IsJetLow()){
        fHistos->FillTHnSparse("PatchInfoJetLowOnline", triggerpatchinfo);
      }
      if(triggerpatch->IsGammaHigh()){
        fHistos->FillTHnSparse("PatchInfoGammaHighOnline", triggerpatchinfo);
      }
      if(triggerpatch->IsGammaLow()){
        fHistos->FillTHnSparse("PatchInfoGammaLowOnline", triggerpatchinfo);
      }
      if(triggerpatch->IsLevel0()){
        fHistos->FillTHnSparse("PatchInfoLevel0Online", triggerpatchinfo);
      }
    }
  }
}

} /* namespace EMCalTriggerPtAnalysis */
