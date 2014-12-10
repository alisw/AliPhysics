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
  AliEMCalTriggerBinningDimension *etabinning = fBinning->GetBinning("eta"),
      *phibinning = fBinning->GetBinning("phi");
  const TAxis *patchenergyaxes[4] = {
    DefineAxis("energy", 100, 0., 100),
    DefineAxis("eta", etabinning),
    DefineAxis("phi", phibinning),
    DefineAxis("isMain", 2, -0.5, 1.5)
  };
  const TAxis *patchampaxes[4] = {
    DefineAxis("amplitude", 10000, 0., 10000.),
    DefineAxis("eta", etabinning),
    DefineAxis("phi",  phibinning),
    DefineAxis("isMain", 2, -0.5, 1.5)
  };

  std::string patchnames[] = {"Level0", "JetHigh", "JetLow", "GammaHigh", "GammaLow"};
  for(std::string * triggerpatch = patchnames; triggerpatch < patchnames + sizeof(patchnames)/sizeof(std::string); ++triggerpatch){
    fHistos->CreateTHnSparse(Form("Energy%s", triggerpatch->c_str()), Form("Patch energy for %s trigger patches", triggerpatch->c_str()), 4, patchenergyaxes, "s");
    fHistos->CreateTHnSparse(Form("EnergyRough%s", triggerpatch->c_str()), Form("Rough patch energy for %s trigger patches", triggerpatch->c_str()), 4, patchenergyaxes, "s");
    fHistos->CreateTHnSparse(Form("Amplitude%s", triggerpatch->c_str()), Form("Patch amplitude for %s trigger patches", triggerpatch->c_str()), 4, patchampaxes, "s");
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
    double triggerpatchinfo[4] = {triggerpatch->GetPatchE(), triggerpatch->GetEtaGeo(), triggerpatch->GetPhiGeo(), triggerpatch->IsMainTrigger() ? 1. : 0.};
    double triggerpatchinfoamp[4] = {static_cast<double>(triggerpatch->GetADCAmp()), triggerpatch->GetEtaGeo(), triggerpatch->GetPhiGeo(), triggerpatch->IsMainTrigger() ? 1. : 0.};
    double triggerpatchinfoer[4] = {triggerpatch->GetADCAmpGeVRough(), triggerpatch->GetEtaGeo(), triggerpatch->GetPhiGeo(), triggerpatch->IsMainTrigger() ? 1. : 0.};
    if(triggerpatch->IsJetHigh()){
      fHistos->FillTHnSparse("EnergyJetHigh", triggerpatchinfo);
      fHistos->FillTHnSparse("AmplitudeJetHigh", triggerpatchinfoamp);
      fHistos->FillTHnSparse("EnergyRoughJetHigh", triggerpatchinfoer);
    }
    if(triggerpatch->IsJetLow()){
      fHistos->FillTHnSparse("EnergyJetLow", triggerpatchinfo);
      fHistos->FillTHnSparse("AmplitudeJetLow", triggerpatchinfoamp);
      fHistos->FillTHnSparse("EnergyRoughJetLow", triggerpatchinfoer);
    }
    if(triggerpatch->IsGammaHigh()){
      fHistos->FillTHnSparse("EnergyGammaHigh", triggerpatchinfo);
      fHistos->FillTHnSparse("AmplitudeGammaHigh", triggerpatchinfoamp);
      fHistos->FillTHnSparse("EnergyRoughGammaHigh", triggerpatchinfoer);
    }
    if(triggerpatch->IsGammaLow()){
      fHistos->FillTHnSparse("EnergyGammaLow", triggerpatchinfo);
      fHistos->FillTHnSparse("AmplitudeGammaLow", triggerpatchinfoamp);
      fHistos->FillTHnSparse("EnergyRoughGammaLow", triggerpatchinfoer);
    }
    if(triggerpatch->IsLevel0()){
      fHistos->FillTHnSparse("EnergyLevel0", triggerpatchinfo);
      fHistos->FillTHnSparse("AmplitudeLevel0", triggerpatchinfoamp);
      fHistos->FillTHnSparse("EnergyRoughLevel0", triggerpatchinfoer);
    }
  }
}

} /* namespace EMCalTriggerPtAnalysis */
