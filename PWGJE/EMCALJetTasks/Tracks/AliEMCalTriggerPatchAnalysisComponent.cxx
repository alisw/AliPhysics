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
#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <TClonesArray.h>
#include <THistManager.h>
#include <TLinearBinning.h>

#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerPatchAnalysisComponent.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerPatchAnalysisComponent)

using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy constructor, only for ROOT I/O, not to be used by the users. Sets as default swapping of
 * the trigger thresholds to false for both types of patches.
 */
AliEMCalTriggerPatchAnalysisComponent::AliEMCalTriggerPatchAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent(),
  fSwapOnlineThresholds(kFALSE),
  fSwapOfflineThresholds(kFALSE),
  fWithEventSelection(kFALSE)
{
}

/**
 * Main constructor, to be used to create the analysis component. Sets as default swapping of
 * the trigger thresholds to false for both types of patches.
 *
 * \param name Name of the component
 * \param withEventSelection In case of true, histograms are created for different event selection classes
 */
AliEMCalTriggerPatchAnalysisComponent::AliEMCalTriggerPatchAnalysisComponent(const char *name, Bool_t withEventSelection) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fSwapOnlineThresholds(kFALSE),
  fSwapOfflineThresholds(kFALSE),
  fWithEventSelection(withEventSelection)
{
}

/**
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

  const TBinning *etabinning = fBinning->GetBinning("eta"),
      *phibinning = fBinning->GetBinning("phi");
  const TAxis *patchaxes[6] = {
    DefineAxis("energy", TLinearBinning(100, 0., 100)),
    DefineAxis("energyRough", TLinearBinning(100, 0., 100)),
    DefineAxis("amplitude", TLinearBinning(2100, 0., 2100.)),       // limit for the moment
    DefineAxis("eta", *etabinning),
    DefineAxis("phi", *phibinning),
    DefineAxis("isMain", TLinearBinning(2, -0.5, 1.5))
  };

  const TAxis *ampaxes[5] = {
     DefineAxis("amplitudeOnline", TLinearBinning(2000, 0., 2100.)),
     DefineAxis("amplitudeOffline", TLinearBinning(2000, 0., 2100.)),
     DefineAxis("eta", *etabinning),
     DefineAxis("phi", *phibinning),
     DefineAxis("isMain", TLinearBinning(2, -0.5, 1.5))
  };

  // Create trigger definitions
  std::map<std::string, std::string> triggerCombinations;
  GetAllTriggerNamesAndTitles(triggerCombinations);

  std::string patchnames[] = {"Level0", "JetHigh", "JetLow", "GammaHigh", "GammaLow"};
  std::string triggermodes[] = {"Online", "Offline"};
  for(std::string * triggerpatch = patchnames; triggerpatch < patchnames + sizeof(patchnames)/sizeof(std::string); ++triggerpatch){
    for(std::string *triggermode = triggermodes; triggermode < triggermodes + sizeof(triggermodes)/sizeof(std::string); ++triggermode){
      if(fWithEventSelection){
        for(std::map<std::string, std::string>::iterator trgiter = triggerCombinations.begin(); trgiter != triggerCombinations.end(); ++trgiter){
          printf("Adding patch for trigger %s in case of %s for event selection %s\n", triggerpatch->c_str(), triggermode->c_str(), trgiter->second.c_str());
          fHistos->CreateTHnSparse(Form("PatchInfo%s%s%s", triggerpatch->c_str(), triggermode->c_str(), trgiter->first.c_str()), Form("Patch energy for %s %s trigger patches", triggerpatch->c_str(), triggermode->c_str()), 6, patchaxes, "s");
        }
      } else {
        if((!strcmp(triggermode->c_str(), "Offline")) && (!strcmp(triggerpatch->c_str(), "Level0"))) continue; // Don't process L0 in case of offline
        printf("Adding patch for trigger %s in case of %s\n", triggerpatch->c_str(), triggermode->c_str());
        fHistos->CreateTHnSparse(Form("PatchInfo%s%s", triggerpatch->c_str(), triggermode->c_str()), Form("Patch energy for %s %s trigger patches", triggerpatch->c_str(), triggermode->c_str()), 6, patchaxes, "s");
        if(strcmp(triggerpatch->c_str(), "Level0")){
          // Add histogram for online-offline of amplitudes
          fHistos->CreateTHnSparse(Form("PatchAmplitude%s%s", triggerpatch->c_str(), triggermode->c_str()), Form("Patch amplitudes for %s %s trigger patches", triggerpatch->c_str(), triggermode->c_str()), 5, ampaxes, "s");
        }
      }
    }
  }
}

/**
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
  AliEMCALTriggerPatchInfo *triggerpatch(NULL);
  TIter patchIter(data->GetTriggerPatchContainer());
  while((triggerpatch = dynamic_cast<AliEMCALTriggerPatchInfo *>(patchIter()))){
    if(fWithEventSelection){
      std::vector<std::string> triggernames;
      GetMachingTriggerNames(triggernames);
      for(std::vector<std::string>::iterator trgclassit = triggernames.begin(); trgclassit != triggernames.end(); ++trgclassit){
        FillStandardMonitoring(triggerpatch, trgclassit->c_str());
      }
    } else {
      FillStandardMonitoring(triggerpatch);
    }
  }
}

/**
 * Fill monitoring of the trigger patches (energies and amplitudes). If event type is specified, only the general histogram
 * will be filled. Otherwise all histograms of classes which are fulfilled by the histogram are filled
 * \param patch Trigger patch to be processed
 * \param eventType Trigger class the event was selected by (optional)
 */
void AliEMCalTriggerPatchAnalysisComponent::FillStandardMonitoring(const AliEMCALTriggerPatchInfo* const patch, TString eventType) {
  std::vector<TString> triggerclasses;
  triggerclasses.push_back("JetHigh");
  triggerclasses.push_back("JetLow");
  triggerclasses.push_back("GammaHigh");
  triggerclasses.push_back("GammaLow");
  triggerclasses.push_back("Level0");
  AliEmcalTriggerPatchHandlerFactory mypatchhandler(fSwapOnlineThresholds, fSwapOfflineThresholds);
  for(std::vector<TString>::iterator trgclassit = triggerclasses.begin(); trgclassit != triggerclasses.end(); ++trgclassit){
    if(mypatchhandler.IsPatchOfType(patch, *trgclassit)){
      std::stringstream infohistname;
      infohistname << "PatchInfo" << trgclassit->Data() << (patch->IsOfflineSimple() ? "Offline" : "Online");
      if(eventType.Length()) infohistname << eventType.Data();
      FillTriggerInfoHistogram(infohistname.str().c_str(), patch);
      if((!eventType.Length()) &&(*trgclassit != "Level0")){
        std::stringstream amphistname;
        amphistname << "PatchAmplitude" << trgclassit->Data() << (patch->IsOfflineSimple() ? "Offline" : "Online");
        FillAmplitudeHistogram(amphistname.str().c_str(), patch);
      }
    }
  }
}

/**
 * Fill standard trigger patch info histogram
 * \param histo Name of the histogram to fill
 * \param patch Patch with information
 */
void AliEMCalTriggerPatchAnalysisComponent::FillTriggerInfoHistogram(TString histo, const AliEMCALTriggerPatchInfo* const patch) {
	bool isMain = patch->IsOfflineSimple() ? patch->IsMainTriggerSimple() : patch->IsMainTrigger();
	double triggerpatchinfo[6] = {patch->GetPatchE(),patch->GetADCAmpGeVRough(),
	    static_cast<double>(patch->IsOfflineSimple() ? patch->GetADCOfflineAmp() : patch->GetADCAmp()), patch->GetEtaGeo(),
	    patch->GetPhiGeo(), isMain ? 1. : 0.};
	fHistos->FillTHnSparse(histo.Data(), triggerpatchinfo);
}


/**
 * Fill histogram for patch amplitude
 * \param histo Name of the histogram to fill
 * \param patch Patch with information
 */
void AliEMCalTriggerPatchAnalysisComponent::FillAmplitudeHistogram(TString histo, const AliEMCALTriggerPatchInfo* const patch) {
	bool isMain = patch->IsOfflineSimple() ? patch->IsMainTriggerSimple() : patch->IsMainTrigger();
	double amplitudeinfo[5] = {(double)patch->GetADCAmp(), (double)patch->GetADCOfflineAmp(), (double)patch->GetEtaGeo(), (double)patch->GetPhiGeo(), isMain ? 1. : 0.};
  fHistos->FillTHnSparse(histo.Data(), amplitudeinfo);
}

//
// Implementation of patch handlers and patch handler factory
//

/**
 * Patch handler for trigger class Jet Low: checks if the patch is a jet low patch. Handles online and offline patches and swapping of the
 * threshold.
 * \param patch The patch to check
 * \return True if it is a jet low patch, false otherwise
 */
Bool_t AliEMCalTriggerPatchAnalysisComponent::AliEmcalTriggerPatchHandlerFactory::AliEmcalTriggerPatchHandlerJetLow::IsOfType(const AliEMCALTriggerPatchInfo* const patch) const {
  Bool_t result = false;
  if(patch->IsOfflineSimple()){
    result = ((!fPatchSwapThresholdsOffline && patch->IsJetLowSimple()) || (fPatchSwapThresholdsOffline && patch->IsJetHighSimple()));
  } else {
    result = ((!fPatchSwapThresholdsOnline && patch->IsJetLow()) || (fPatchSwapThresholdsOnline && patch->IsJetHigh()));
  }
  return result;
}

/**
 * Patch handler for trigger class Jet High: checks if the patch is a jet low patch. Handles online and offline patches and swapping of the
 * threshold.
 * \param patch The patch to check
 * \return True if it is a jet high patch, false otherwise
 */
Bool_t AliEMCalTriggerPatchAnalysisComponent::AliEmcalTriggerPatchHandlerFactory::AliEmcalTriggerPatchHandlerJetHigh::IsOfType(const AliEMCALTriggerPatchInfo* const patch) const {
  Bool_t result = false;
  if(patch->IsOfflineSimple()){
    result = ((!fPatchSwapThresholdsOffline && patch->IsJetHighSimple()) || (fPatchSwapThresholdsOffline && patch->IsJetLowSimple()));
  } else {
    result = ((!fPatchSwapThresholdsOnline && patch->IsJetHigh()) || (fPatchSwapThresholdsOnline && patch->IsJetLow()));
  }
  return result;
}

/**
 * Patch handler for trigger class Gamma Low: checks if the patch is a jet low patch. Handles online and offline patches and swapping of the
 * threshold.
 * \param patch The patch to check
 * \return True if it is a gamma low patch, false otherwise
 */
Bool_t AliEMCalTriggerPatchAnalysisComponent::AliEmcalTriggerPatchHandlerFactory::AliEmcalTriggerPatchHandlerGammaLow::IsOfType(const AliEMCALTriggerPatchInfo* const patch) const {
  Bool_t result = false;
  if(patch->IsOfflineSimple()){
    result = ((!fPatchSwapThresholdsOffline && patch->IsGammaLowSimple()) || (fPatchSwapThresholdsOffline && patch->IsGammaHighSimple()));
  } else {
    result = ((!fPatchSwapThresholdsOnline && patch->IsGammaLow()) || (fPatchSwapThresholdsOnline && patch->IsGammaHigh()));
  }
  return result;
}

/**
 * Patch handler for trigger class Gamma High: checks if the patch is a jet low patch. Handles online and offline patches and swapping of the
 * threshold.
 * \param patch The patch to check
 * \return True if it is a gamma high patch, false otherwise
 */
Bool_t AliEMCalTriggerPatchAnalysisComponent::AliEmcalTriggerPatchHandlerFactory::AliEmcalTriggerPatchHandlerGammaHigh::IsOfType(const AliEMCALTriggerPatchInfo* const patch) const {
  Bool_t result = false;
  if(patch->IsOfflineSimple()){
    result = ((!fPatchSwapThresholdsOffline && patch->IsGammaHighSimple()) || (fPatchSwapThresholdsOffline && patch->IsGammaLowSimple()));
  } else {
    result = ((!fPatchSwapThresholdsOnline && patch->IsGammaHigh()) || (fPatchSwapThresholdsOnline && patch->IsGammaLow()));
  }
  return result;
}

/**
 * Patch handler for trigger class Gamma High: checks if the patch is a jet low patch. Handles online and offline patches and swapping of the
 * threshold.
 * \param patch The patch to check
 * \return True if it is a gamma high patch, false otherwise
 */
Bool_t AliEMCalTriggerPatchAnalysisComponent::AliEmcalTriggerPatchHandlerFactory::AliEmcalTriggerPatchHandlerLevel0::IsOfType(const AliEMCALTriggerPatchInfo* const patch) const {
  if(patch->IsLevel0()) return true;
  return false;
}

/**
 *
 * \param patch The patch to check
 * \param patchtype The trigger class to check
 * \return True if the trigger class is
 */
Bool_t AliEMCalTriggerPatchAnalysisComponent::AliEmcalTriggerPatchHandlerFactory::IsPatchOfType(
  const AliEMCALTriggerPatchInfo* const patch, TString patchtype) const {
  std::unique_ptr<AliEmcalTriggerPatchHandler> myhandler;
  if (patchtype == "JetHigh") myhandler = std::unique_ptr<AliEmcalTriggerPatchHandler>(new AliEmcalTriggerPatchHandlerGammaHigh(fSwapThresholdsOnline, fSwapThresholdsOffline));
  else if (patchtype == "JetLow") myhandler = std::unique_ptr<AliEmcalTriggerPatchHandler>(new AliEmcalTriggerPatchHandlerGammaLow(fSwapThresholdsOnline, fSwapThresholdsOffline));
  else if (patchtype == "GammaHigh") myhandler = std::unique_ptr<AliEmcalTriggerPatchHandler>(new AliEmcalTriggerPatchHandlerJetHigh(fSwapThresholdsOnline, fSwapThresholdsOffline));
  else if (patchtype == "GammaLow") myhandler = std::unique_ptr<AliEmcalTriggerPatchHandler>(new AliEmcalTriggerPatchHandlerJetLow(fSwapThresholdsOnline, fSwapThresholdsOffline));
  else if (patchtype == "Level0") myhandler = std::unique_ptr<AliEmcalTriggerPatchHandler>(new AliEmcalTriggerPatchHandlerLevel0(fSwapThresholdsOnline, fSwapThresholdsOffline));
  if(!myhandler.get()) return false;
  return myhandler->IsOfType(patch);
}
