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
#include <functional>
#include <iostream>

#include <TClonesArray.h>
#include <TGrid.h>
#include <THistManager.h>
#include <TParameter.h>

#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliAnalysisUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerMapping.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliOADBContainer.h"
#include "AliVVertex.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalTriggerBase)
/// \endcond

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskEmcalTriggerBase::AliAnalysisTaskEmcalTriggerBase():
  AliAnalysisTaskEmcal(),
  fTriggerSelection(nullptr),
  fHistos(nullptr),
  fTriggerStringFromPatches(kFALSE),
  fSelectedTriggers(),
  fVertexCut(-10., 10.),
  fNameDownscaleOADB(""),
  fDownscaleOADB(nullptr),
  fDownscaleFactors(nullptr),
  fNameMaskedFastorOADB(),
  fMaskedFastorOADB(nullptr),
  fMaskedFastors(),
  fSelectNoiseEvents(false),
  fRejectNoiseEvents(false),
  fUseOnlineTriggerSelectorV1(false)
{
  SetNeedEmcalGeom(true);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalTriggerBase::AliAnalysisTaskEmcalTriggerBase(const char *name):
  AliAnalysisTaskEmcal(name, true),
  fTriggerSelection(nullptr),
  fHistos(nullptr),
  fTriggerStringFromPatches(kFALSE),
  fSelectedTriggers(),
  fVertexCut(-10., 10.),
  fNameDownscaleOADB(""),
  fDownscaleOADB(nullptr),
  fDownscaleFactors(nullptr),
  fNameMaskedFastorOADB(),
  fMaskedFastorOADB(nullptr),
  fMaskedFastors(),
  fSelectNoiseEvents(false),
  fRejectNoiseEvents(false),
  fUseOnlineTriggerSelectorV1(false)
{
  SetNeedEmcalGeom(true);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalTriggerBase::~AliAnalysisTaskEmcalTriggerBase() {
  if(fTriggerSelection) delete fTriggerSelection;
  if(fHistos) delete fHistos;
  if(fDownscaleOADB) delete fDownscaleOADB;
}

/**
 * Create the output histograms
 */
void AliAnalysisTaskEmcalTriggerBase::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  if(!fAliAnalysisUtils) fAliAnalysisUtils = new AliAnalysisUtils;

  fHistos = new THistManager(Form("Histos_%s", GetName()));

  CreateUserObjects();
  CreateUserHistos();

  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalTriggerBase::IsEventSelected(){
  // Apply trigger selection
  TriggerSelection();
  if(!fSelectedTriggers.size()) return false;

  UserFillHistosBeforeEventSelection();

  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return false;
  if(fInputEvent->IsA() == AliESDEvent::Class() && fAliAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return false;
  if(fAliAnalysisUtils->IsPileUpEvent(fInputEvent)) return false;       // Apply new vertex cut

  if(!fAliAnalysisUtils->IsVertexSelected2013pA(fInputEvent))return false;       // Apply new vertex cut
  if(!fVertexCut.IsInRange(fVertex[2])) return false;
  if(!IsUserEventSelected()) return false;

  // Do MC outlier cut
  if(fIsPythia){
    if(!CheckMCOutliers()){
      AliDebugStream(1) << GetName() << ": Reject MC outliers" << std::endl;
      return false;
    }
  }

  UserFillHistosAfterEventSelection();
  return true;
}


void AliAnalysisTaskEmcalTriggerBase::TriggerSelection(){
  fSelectedTriggers.clear();

  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(fTriggerPatchInfo);
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1"),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2"),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1"),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2"),
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("CEMC7"),
      isDJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ1"),
      isDJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ2"),
      isDG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG1"),
      isDG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG2"),
      isDMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("DMC7");
  if(fTriggerPatchInfo && fTriggerSelection){
      isEJ1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ1, fTriggerPatchInfo);
      isEJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ2, fTriggerPatchInfo);
      isEG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG1, fTriggerPatchInfo);
      isEG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG2, fTriggerPatchInfo);
      isEMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEL0, fTriggerPatchInfo);
      isDJ1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ1, fTriggerPatchInfo);
      isDJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ2, fTriggerPatchInfo);
      isDG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG1, fTriggerPatchInfo);
      isDG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG2, fTriggerPatchInfo);
      isDMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDL0, fTriggerPatchInfo);
  }
  // Online selection / rejection
  if(fRejectNoiseEvents || fSelectNoiseEvents){
    std::function<bool (OnlineTrigger_t)> selector;
    if(fUseOnlineTriggerSelectorV1){
        selector = [this](OnlineTrigger_t t) -> bool { return this->SelectOnlineTriggerV1(t); } ;
    } else{
        selector = [this](OnlineTrigger_t t) -> bool { return this->SelectOnlineTrigger(t); };
    }
    if(fRejectNoiseEvents){
      if(isEJ1) isEJ1 &= selector(kCPREJ1);
      if(isEJ2) isEJ2 &= selector(kCPREJ2);
      if(isEG1) isEG1 &= selector(kCPREG1);
      if(isEG2) isEG2 &= selector(kCPREG2);
      if(isEMC7) isEMC7 &= selector(kCPREL0);
    } else {
      if(isEJ1) isEJ1 &= !selector(kCPREJ1);
      if(isEJ2) isEJ2 &= !selector(kCPREJ2);
      if(isEG1) isEG1 &= !selector(kCPREG1);
      if(isEG2) isEG2 &= !selector(kCPREG2);
      if(isEMC7) isEMC7 &= !selector(kCPREL0);
    }
  }
  if(isMinBias) fSelectedTriggers.push_back("MB");
  if(isEMC7){
    fSelectedTriggers.push_back("EMC7");
    if(!isMinBias) fSelectedTriggers.push_back("EMC7excl");
  }
  if(isDMC7){
    fSelectedTriggers.push_back("DMC7");
    if(!isMinBias) fSelectedTriggers.push_back("DMC7excl");
  }
  if(isEJ2){
    fSelectedTriggers.push_back("EJ2");
    if(!(isMinBias || isEMC7)) fSelectedTriggers.push_back("EJ2excl");
  }
  if(isEJ1){
    fSelectedTriggers.push_back("EJ1");
    if(!(isMinBias || isEMC7 || isEJ2)) fSelectedTriggers.push_back("EJ1excl");
  }
  if(isDJ2){
    fSelectedTriggers.push_back("DJ2");
    if(!(isMinBias || isDMC7)) fSelectedTriggers.push_back("DJ2excl");
  }
  if(isDJ1){
    fSelectedTriggers.push_back("DJ1");
    if(!(isMinBias || isDMC7 || isDJ2)) fSelectedTriggers.push_back("DJ1excl");
  }
  if(isEG2){
    fSelectedTriggers.push_back("EG2");
    if(!(isMinBias || isEMC7)) fSelectedTriggers.push_back("EG2excl");
  }
  if(isEG1){
    fSelectedTriggers.push_back("EG1");
    if(!(isMinBias || isEMC7 || isEG2)) fSelectedTriggers.push_back("EG1excl");
  }
  if(isDG2){
    fSelectedTriggers.push_back("DG2");
    if(!(isMinBias || isDMC7)) fSelectedTriggers.push_back("DG2excl");
  }
  if(isDG1){
    fSelectedTriggers.push_back("DG1");
    if(!(isMinBias || isDMC7 || isDG2)) fSelectedTriggers.push_back("DG1excl");
  }
}

/**
 * Perform gloabl initializations
 * Used for the moment for
 * - Initialize OADB container with downscaling factors
 */
void AliAnalysisTaskEmcalTriggerBase::ExecOnce(){
  AliAnalysisTaskEmcal::ExecOnce();

  if(!fLocalInitialized){
    return;
  }

  // Handle OADB container with downscaling factors
  if(fNameDownscaleOADB.Length()){
    if(fNameDownscaleOADB.Contains("alien://") && ! gGrid) TGrid::Connect("alien://");
    fDownscaleOADB = new AliOADBContainer("AliEmcalDownscaleFactors");
    fDownscaleOADB->InitFromFile(fNameDownscaleOADB.Data(), "AliEmcalDownscaleFactors");
  }
  // Load OADB container with masked fastors (in case fastor masking is switched on)
  if(fNameMaskedFastorOADB.Length() && (fRejectNoiseEvents || fSelectNoiseEvents)){
    if(fNameMaskedFastorOADB.Contains("alien://") && ! gGrid) TGrid::Connect("alien://");
    fMaskedFastorOADB = new AliOADBContainer("AliEmcalMaskedFastors");
    fMaskedFastorOADB->InitFromFile(fNameMaskedFastorOADB.Data(), "AliEmcalMaskedFastors");
  }
}

/**
 * Run change method. Called when the run number of the new event
 * is different compared to the run number of the previous event.
 * Used for loading of the downscale factor for a given
 * run from the downscale OADB.
 * @param[in] runnumber Number of the new run.
 */
void AliAnalysisTaskEmcalTriggerBase::RunChanged(Int_t runnumber){
  if(fDownscaleOADB){
    fDownscaleFactors = static_cast<TObjArray *>(fDownscaleOADB->GetObject(runnumber));
  }
  if(fMaskedFastorOADB){
   fMaskedFastors.clear();
   TObjArray *ids = static_cast<TObjArray *>(fMaskedFastorOADB->GetObject(runnumber));
   for(auto m : *ids){
     TParameter<int> *id = static_cast<TParameter<int> *>(m);
     fMaskedFastors.push_back(id->GetVal());
   }
 }
}

std::vector<TString> AliAnalysisTaskEmcalTriggerBase::GetSupportedTriggers(){
  // Exclusive means classes without lower trigger classes (which are downscaled) -
  // in order to make samples statistically independent: MBExcl means MinBias && !EMCAL trigger
  std::vector<TString> triggers = {
    "MB", "EMC7", "DMC7",
    "EJ1", "EJ2", "EG1", "EG2", "DJ1", "DJ2", "DG1", "DG2",
    "EMC7excl", "DMC7excl", "EG2excl", "EJ2excl", "DG2excl", "DJ2excl",
    "EJ1excl", "DJ1excl", "EG1excl", "DG1excl"
  };
  return triggers;
}


/**
 * Select events which contain good patches (at least one patch without noisy fastor).
 * @param trg L0/L1 trigger class to be checked
 * @return True if the event has at least one good patch, false otherwise
 */
bool AliAnalysisTaskEmcalTriggerBase::SelectOnlineTrigger(OnlineTrigger_t trg){
  AliDebugStream(1) << "Using nominal online trigger selector" << std::endl;
  int ngood(0);
  const int kEJ1threshold = 223, kEJ2threshold = 140;
  std::function<bool(const AliEMCALTriggerPatchInfo *)> PatchSelector[5] = {
      [](const AliEMCALTriggerPatchInfo * patch) -> bool {
          return patch->IsGammaHigh();
      },
      [](const AliEMCALTriggerPatchInfo * patch) -> bool {
          return patch->IsGammaLow();
      },
      [kEJ1threshold](const AliEMCALTriggerPatchInfo * patch) -> bool {
          return patch->IsJetLowRecalc() && patch->GetADCAmp() > kEJ1threshold;
      },
      [kEJ2threshold](const AliEMCALTriggerPatchInfo * patch) -> bool {
          return patch->IsJetLowRecalc() && patch->GetADCAmp() > kEJ2threshold;
      },
      [](const AliEMCALTriggerPatchInfo * patch) -> bool {
          return patch->IsLevel0();
      }
  };
  for(auto p : *fTriggerPatchInfo){
    AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(p);
    if(PatchSelector[trg](patch)){
      bool patchMasked(false);
      for(int icol = patch->GetColStart(); icol < patch->GetColStart() + patch->GetPatchSize(); icol++){
        for(int irow = patch->GetRowStart(); irow < patch->GetRowStart() + patch->GetPatchSize(); irow++){
          int fastorAbs(-1);
          fGeom->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, fastorAbs);
          if(std::find(fMaskedFastors.begin(), fMaskedFastors.end(), fastorAbs) != fMaskedFastors.end()){
            patchMasked = true;
            break;
          }
        }
      }
      if(!patchMasked) ngood++;
    }
  }
  return ngood > 0;
}

/**
 * Second approach: We assume masked fastors are already handled in the
 * trigger maker. In this case we treat masked fastors similarly to online
 * masked fastors and apply online cuts on the recalc ADC value. Implemented
 * for the moment only for L1 triggers.
 * @return True if the event has at least 1 recalc patch above threshold
 */
bool AliAnalysisTaskEmcalTriggerBase::SelectOnlineTriggerV1(OnlineTrigger_t trigger){
  AliDebugStream(1) << "Using V1 online trigger selector" << std::endl;
  if(trigger == kCPREL0) return true;
  double onlinethresholds[5] = {140., 89., 260., 127., 0.};
  int ngood(0);
  for(auto p : *fTriggerPatchInfo){
    AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(p);
    if(((trigger == kCPREG1 || trigger == kCPREG2) && patch->IsGammaLowRecalc()) ||
        ((trigger == kCPREJ1 || trigger == kCPREJ2) && patch->IsJetLowRecalc())){
      if(patch->GetADCAmp() > onlinethresholds[trigger]) ngood++;
    }
  }
  return ngood > 0;
}

/**
 * Get a trigger class dependent event weight. The weight
 * is defined as 1/downscalefactor. The downscale factor
 * is taken from the OADB. For triggers which are not downscaled
 * the weight is always 1.
 * @param[in] triggerclass Class for which to obtain the trigger.
 * @return Downscale facror for the trigger class (1 if trigger is not downscaled or no OADB container is available)
 */
Double_t AliAnalysisTaskEmcalTriggerBase::GetTriggerWeight(const TString &triggerclass) const {
  if(fDownscaleFactors){
    TParameter<double> *result(nullptr);
    // Downscaling only done on MB, L0 and the low threshold triggers
    if(triggerclass.Contains("MB")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("INT7"));
    else if(triggerclass.Contains("EMC7")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EMC7"));
    else if(triggerclass.Contains("EJ2")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ2"));
    else if(triggerclass.Contains("EG2")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG2"));
    if(result) return 1./result->GetVal();
  }
  return 1.;
}

/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param triggerpatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskEmcalTriggerBase::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
  TString triggerstring = "";
  Int_t nEJ1 = 0, nEJ2 = 0, nEG1 = 0, nEG2 = 0;
  double  minADC_EJ1 = 260.,
          minADC_EJ2 = 127.,
          minADC_EG1 = 140.,
          minADC_EG2 = 89.;
  for(auto patchIter : *triggerpatches){
    AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(patch->IsJetHighSimple() && patch->GetADCOfflineAmp() > minADC_EJ1) nEJ1++;
    if(patch->IsJetLowSimple() && patch->GetADCOfflineAmp() > minADC_EJ2) nEJ2++;
    if(patch->IsGammaHighSimple() && patch->GetADCOfflineAmp() > minADC_EG1) nEG1++;
    if(patch->IsGammaLowSimple() && patch->GetADCOfflineAmp() > minADC_EG2) nEG2++;
  }
  if(nEJ1) triggerstring += "EJ1";
  if(nEJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EJ2";
  }
  if(nEG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG1";
  }
  if(nEG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG2";
  }
  return triggerstring;
}

} /* namespace EMCalTriggerPtAnalysis */
