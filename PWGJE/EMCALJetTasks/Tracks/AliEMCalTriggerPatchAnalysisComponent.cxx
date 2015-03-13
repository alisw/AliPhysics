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
#include <TClonesArray.h>

#include "AliEmcalTriggerPatchInfo.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerPatchAnalysisComponent.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerPatchAnalysisComponent)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * \brief Default constructor
 *
 * Dummy constructor, only for ROOT I/O, not to be used by the users. Sets as default swapping of
 * the trigger thresholds to false for both types of patches.
 */
AliEMCalTriggerPatchAnalysisComponent::AliEMCalTriggerPatchAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent(),
  fSwapOnlineThresholds(kFALSE),
  fSwapOfflineThresholds(kFALSE)
{
}

/**
 * \brief Main constructor
 *
 * Main constructor, to be used to create the analysis component. Sets as default swapping of
 * the trigger thresholds to false for both types of patches.
 *
 * \param name Name of the component
 */
AliEMCalTriggerPatchAnalysisComponent::AliEMCalTriggerPatchAnalysisComponent(const char *name) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fSwapOnlineThresholds(kFALSE),
  fSwapOfflineThresholds(kFALSE)
{
}

/**
 * \brief Create histograms for the trigger patch analysis
 *
 * Create histograms for the trigger patch analysis. Two types of histograms are currently defined.
 *  -# A histogram correlating different energy values (amplitude, estimated energy, calibrated energy)
 *  -# A histogram correlating online amplitude and offline amplitude
 * All histograms contain the trigger patch position in \f$\eta\f$ and \f$\phi\f$ as well. In case the
 * patch is a main patch, defined as the patch of a given trigger type with the highest energy, this
 * patch is marked as well.
 *
 * Each trigger type has its own histogram. Currently the following types are supported:
 *  -# EMCal Level1 Jet, high threshold (EMCJHigh)
 *  -# EMCal Level1 Jet, low threshold (EMCJLow)
 *  -# EMCal Level1 Gamma, high threshold (EMCGHigh)
 *  -# EMCal Level1 Gamma, low threshold (EMCGLow)
 *  -# EMCal Level0
 * For all Level1 triggers, histograms are separated for online and offline patches. For Level0 only offline
 * patches are available.
 *
 * This function is the implementation of the abstract method CreateHistos declared in AliEMCalTriggerTracksAnalysisComponent.
 */
void AliEMCalTriggerPatchAnalysisComponent::CreateHistos() {
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

  const TAxis *ampaxes[5] = {
     DefineAxis("amplitudeOnline", 5000, 0., 5000.),
     DefineAxis("amplitudeOffline", 5000, 0., 5000.),
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
      if(strcmp(triggerpatch->c_str(), "Level0")){
        // Add histogram for online-offline of amplitudes
        fHistos->CreateTHnSparse(Form("PatchAmplitudes%s%s", triggerpatch->c_str(), triggermode->c_str()), Form("Patch amplitudes for %s %s trigger patches", triggerpatch->c_str(), triggermode->c_str()), 5, ampaxes, "s");
      }
    }
  }
}

/**
 * \brief Performs analysis on trigger patches found in the event
 *
 * Perform analysis on the event data set. Loops over all trigger patches found by the trigger patch
 * maker and fill relevant histograms for the different trigger types associated with the trigger patch. In
 * case the thresholds are requested to be swapped (separately for online and offline patches), then low
 * threshold patches will be used for high threshold histograms and vice versa.
 *
 * This function is the implementation of the abstract method Process declared in AliEMCalTriggerTracksAnalysisComponent.
 *
 * \param data Event information
 */
void AliEMCalTriggerPatchAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  AliEmcalTriggerPatchInfo *triggerpatch(NULL);
  TIter patchIter(data->GetTriggerPatchContainer());
  while((triggerpatch = dynamic_cast<AliEmcalTriggerPatchInfo *>(patchIter()))){
	bool isMain = triggerpatch->IsOfflineSimple() ? triggerpatch->IsMainTriggerSimple() : triggerpatch->IsMainTrigger();
    double triggerpatchinfo[6] = {triggerpatch->GetPatchE(),triggerpatch->GetADCAmpGeVRough(),
        static_cast<double>(triggerpatch->IsOfflineSimple() ? triggerpatch->GetADCOfflineAmp() : triggerpatch->GetADCAmp()), triggerpatch->GetEtaGeo(),
        triggerpatch->GetPhiGeo(), isMain ? 1. : 0.};
    double amplitudeinfo[5] = {triggerpatch->GetADCAmp(), triggerpatch->GetADCOfflineAmp(), triggerpatch->GetEtaGeo(), triggerpatch->GetPhiGeo(), isMain ? 1. : 0.};
    if(triggerpatch->IsOfflineSimple()){
      if((!fSwapOfflineThresholds && triggerpatch->IsJetHighSimple()) || (fSwapOfflineThresholds && triggerpatch->IsJetLowSimple())){
        fHistos->FillTHnSparse("PatchInfoJetHighOffline", triggerpatchinfo);
        fHistos->FillTHnSparse("PatchAmplitudesJetHighOffline", amplitudeinfo);
      }
      if((!fSwapOfflineThresholds && triggerpatch->IsJetLowSimple()) || (fSwapOfflineThresholds && triggerpatch->IsJetHighSimple())){
    	  fHistos->FillTHnSparse("PatchInfoJetLowOffline", triggerpatchinfo);
    	  fHistos->FillTHnSparse("PatchAmplitudesJetLowOffline", amplitudeinfo);
      }
      if((!fSwapOfflineThresholds && triggerpatch->IsGammaHighSimple()) || (fSwapOfflineThresholds && triggerpatch->IsGammaLowSimple())){
        fHistos->FillTHnSparse("PatchInfoGammaHighOffline", triggerpatchinfo);
        fHistos->FillTHnSparse("PatchAmplitudeGammaHighOffline", amplitudeinfo);
      }
      if((!fSwapOfflineThresholds && triggerpatch->IsGammaLowSimple()) || (fSwapOfflineThresholds && triggerpatch->IsGammaHighSimple())){
        fHistos->FillTHnSparse("PatchInfoGammaLowOffline", triggerpatchinfo);
        fHistos->FillTHnSparse("PatchAmplitudeGammaLowOffline", amplitudeinfo);
      }
    } else{
      if((!fSwapOnlineThresholds && triggerpatch->IsJetHigh()) || (fSwapOnlineThresholds && triggerpatch->IsJetLow())){
        fHistos->FillTHnSparse("PatchInfoJetHighOnline", triggerpatchinfo);
        fHistos->FillTHnSparse("PatchAmplitudeJetHighOnline", amplitudeinfo);
      }
      if((!fSwapOnlineThresholds && triggerpatch->IsJetLow()) || (fSwapOnlineThresholds && triggerpatch->IsJetHigh())){
        fHistos->FillTHnSparse("PatchInfoJetLowOnline", triggerpatchinfo);
        fHistos->FillTHnSparse("PatchAmplitudeJetLowOnline", amplitudeinfo);
      }
      if((!fSwapOnlineThresholds && triggerpatch->IsGammaHigh()) || (fSwapOnlineThresholds && triggerpatch->IsGammaLow())){
        fHistos->FillTHnSparse("PatchInfoGammaHighOnline", triggerpatchinfo);
        fHistos->FillTHnSparse("PatchAmplitudeGammaHighOnline", amplitudeinfo);
      }
      if((!fSwapOnlineThresholds && triggerpatch->IsGammaLow()) || (fSwapOnlineThresholds && triggerpatch->IsGammaHigh())){
        fHistos->FillTHnSparse("PatchInfoGammaLowOnline", triggerpatchinfo);
        fHistos->FillTHnSparse("PatchAmplitudeGammaLowOnline", amplitudeinfo);
      }
      if(triggerpatch->IsLevel0()){
        fHistos->FillTHnSparse("PatchInfoLevel0Online", triggerpatchinfo);
      }
    }
  }
}

} /* namespace EMCalTriggerPtAnalysis */
