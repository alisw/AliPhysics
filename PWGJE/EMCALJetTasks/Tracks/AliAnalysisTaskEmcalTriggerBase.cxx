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
#include <array>
#include <functional>
#include <iostream>

#include <TClonesArray.h>
#include <TGrid.h>
#include <THistManager.h>
#include <TParameter.h>

#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerMapping.h"
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
  fUseTriggerBits(kTRUE),
  fRequireBunchCrossing(kTRUE),
  fHistos(nullptr),
  fTriggerStringFromPatches(kFALSE),
  fSelectedTriggers(),
  fNameClusterContainer(""),
  fRequireAnalysisUtils(kTRUE),
  fVertexCut(-10., 10.),
  fNameDownscaleOADB(""),
  fDownscaleOADB(nullptr),
  fDownscaleFactors(nullptr),
  fNameMaskedFastorOADB(),
  fMaskedFastorOADB(nullptr),
  fMaskedFastors(),
  fOnlineTriggerThresholds(),
  fNameAcceptanceOADB(),
  fSelectNoiseEvents(false),
  fRejectNoiseEvents(false),
  fEnableDCALTriggers(true)
{
  SetNeedEmcalGeom(true);
  SetMakeGeneralHistograms(kTRUE);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalTriggerBase::AliAnalysisTaskEmcalTriggerBase(const char *name):
  AliAnalysisTaskEmcal(name, true),
  fTriggerSelection(nullptr),
  fUseTriggerBits(kTRUE),
  fRequireBunchCrossing(kTRUE),
  fHistos(nullptr),
  fTriggerStringFromPatches(kFALSE),
  fSelectedTriggers(),
  fNameClusterContainer(""),
  fRequireAnalysisUtils(kTRUE),
  fVertexCut(-10., 10.),
  fNameDownscaleOADB(""),
  fDownscaleOADB(nullptr),
  fDownscaleFactors(nullptr),
  fNameMaskedFastorOADB(),
  fMaskedFastorOADB(nullptr),
  fMaskedFastors(),
  fOnlineTriggerThresholds(),
  fNameAcceptanceOADB(),
  fSelectNoiseEvents(false),
  fRejectNoiseEvents(false),
  fEnableDCALTriggers(true)
{
  SetNeedEmcalGeom(true);
  SetMakeGeneralHistograms(kTRUE);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalTriggerBase::~AliAnalysisTaskEmcalTriggerBase() {
  if(fTriggerSelection) delete fTriggerSelection;
  if(fHistos) delete fHistos;
  if(fDownscaleOADB) delete fDownscaleOADB;
}

void AliAnalysisTaskEmcalTriggerBase::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  if(fRequireAnalysisUtils && !fAliAnalysisUtils) fAliAnalysisUtils = new AliAnalysisUtils;

  if(!fNameClusterContainer.Length()) fNameClusterContainer = AliEmcalAnalysisFactory::ClusterContainerNameFactory(fInputHandler->IsA() == AliAODInputHandler::Class());
  if(fTriggerSelection && !fTriggerSelection->GetNameClusterContainer().Length()){
    fTriggerSelection->SetClusterContainer(fNameClusterContainer);
  }

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

  if(fRequireAnalysisUtils){
    if(fInputEvent->IsA() == AliESDEvent::Class() && fAliAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return false;
    if(fAliAnalysisUtils->IsPileUpEvent(fInputEvent)) return false;       // Apply new vertex cut
    if(!fAliAnalysisUtils->IsVertexSelected2013pA(fInputEvent))return false;       // Apply new vertex cut
  }

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
  Bool_t isMC = MCEvent() != nullptr;

  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(fTriggerPatchInfo);
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }


  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgn];
  for(int itrg = 0; itrg < AliEmcalTriggerOfflineSelection::kTrgn; itrg++) emcalTriggers[itrg] = true;
  if(!isMC){
    // In case of data select events as bunch-bunch (-B-) events.
    // Cut not applied in simulations
    if(fRequireBunchCrossing && ! triggerstring.Contains("-B-")) return;

    // In case of data use information from the physics selection and event record
    // Further cleanup of trigger events can be performed depending on the presence
    // of recalc patches (after masking hot fastors in the trigger maker) above
    // threshold
    if(fUseTriggerBits){
      const std::array<ULong_t, AliEmcalTriggerOfflineSelection::kTrgn> kSelectTriggerBits = {
    	  AliVEvent::kEMC7, AliVEvent::kEMCEGA, AliVEvent::kEMCEGA, AliVEvent::kEMCEJE, AliVEvent::kEMCEJE,
		  AliVEvent::kEMC7, AliVEvent::kEMCEGA, AliVEvent::kEMCEGA, AliVEvent::kEMCEJE, AliVEvent::kEMCEJE,
      };
      for(int iclass = 0; iclass < AliEmcalTriggerOfflineSelection::kTrgn; iclass++)
        emcalTriggers[iclass] &= bool(selectionstatus & kSelectTriggerBits[iclass]);
    }

    // Apply cut on the trigger string - this basically discriminates high- and low-threshold
    // triggers
    const std::array<TString, AliEmcalTriggerOfflineSelection::kTrgn> kSelectTriggerStrings = {
    		"CEMC7", "EG1", "EG2", "EJ1", "EJ2", "CDMC7", "DG1", "DG2", "DJ1", "DJ2"
    };
    if(triggerstring.Contains("EMC")) AliDebugStream(1) << GetName() << ": Trigger string " << triggerstring << std::endl;
    for(int iclass = 0; iclass < AliEmcalTriggerOfflineSelection::kTrgn; iclass++){
      emcalTriggers[iclass] &= triggerstring.Contains(kSelectTriggerStrings[iclass]);
      if(emcalTriggers[iclass])
        AliDebugStream(1) << GetName() << ": Event selected as trigger " << kSelectTriggerStrings[iclass] << std::endl;
    }

    // Online selection / rejection
    if(fRejectNoiseEvents || fSelectNoiseEvents){
      if(fRejectNoiseEvents){
        for(int itrg = 0; itrg < AliEmcalTriggerOfflineSelection::kTrgn; itrg++){
          if(emcalTriggers[itrg])
            emcalTriggers[itrg] &= SelectOnlineTrigger(AliEmcalTriggerOfflineSelection::EmcalTriggerClass(itrg));
        }
      } else {
        for(int itrg = 0; itrg < AliEmcalTriggerOfflineSelection::kTrgn; itrg++){
          if(emcalTriggers[itrg])
            emcalTriggers[itrg] &= !SelectOnlineTrigger(AliEmcalTriggerOfflineSelection::EmcalTriggerClass(itrg));
        }
      }
    }
  }
  // Apply offline trigger selection: In this case cuts are performed on the
  // patch energy from EMCAL cells after calibration. This method is most relevant
  // for simulations. It can have a special use case in data in case a stronger
  // offline selection is applied in addition to the online selection.
  if(fTriggerPatchInfo && fTriggerSelection){
    for(int itrg = 0; itrg < AliEmcalTriggerOfflineSelection::kTrgn; itrg++)
      emcalTriggers[itrg] &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::EmcalTriggerClass(itrg), fInputEvent);
  }
  if(isMinBias) fSelectedTriggers.push_back("MB");
  if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0]){
    fSelectedTriggers.push_back("EMC7");
    if(!isMinBias) fSelectedTriggers.push_back("EMC7excl");
  }
  if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2]){
    fSelectedTriggers.push_back("EJ2");
    if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EJ2excl");
  }
  if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ1]){
    fSelectedTriggers.push_back("EJ1");
    if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2])) fSelectedTriggers.push_back("EJ1excl");
  }
  if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2]){
    fSelectedTriggers.push_back("EG2");
    if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EG2excl");
  }
  if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEG1]){
    fSelectedTriggers.push_back("EG1");
    if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2])) fSelectedTriggers.push_back("EG1excl");
  }

  if(fEnableDCALTriggers){
    // Handle DCAL triggers only in case DCAL triggers are enabled,
    // otherwise ignore results of the online/offline trigger selection
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0]){
      fSelectedTriggers.push_back("DMC7");
      if(!isMinBias) fSelectedTriggers.push_back("DMC7excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2]){
      fSelectedTriggers.push_back("DJ2");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DJ2excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ1]){
      fSelectedTriggers.push_back("DJ1");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2])) fSelectedTriggers.push_back("DJ1excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDG2]){
      fSelectedTriggers.push_back("DG2");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DG2excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDG1]){
      fSelectedTriggers.push_back("DG1");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDG2])) fSelectedTriggers.push_back("DG1excl");
    }
  }
}

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

  // Setting online threshold for trigger
  if(!OnlineThresholdsInitialized()){
    if(fInputEvent->GetRunNumber() >= 15344 && fInputEvent->GetRunNumber() <= 197388){
      if(!GetOnlineTriggerThresholdByName("EG1")) SetOnlineTriggerThreshold("EG1", 140);
      if(!GetOnlineTriggerThresholdByName("EG2")) SetOnlineTriggerThreshold("EG2", 89);
      if(!GetOnlineTriggerThresholdByName("EJ1")) SetOnlineTriggerThreshold("EJ1", 260);
      if(!GetOnlineTriggerThresholdByName("EJ2")) SetOnlineTriggerThreshold("EJ2", 127);
      SetOnlineTriggerThreshold("DG1", 0);
      SetOnlineTriggerThreshold("DG2", 0);
      SetOnlineTriggerThreshold("DJ1", 0);
      SetOnlineTriggerThreshold("DJ2", 0);
      SetOnlineTriggerThreshold("EMC7", 0);
      SetOnlineTriggerThreshold("DMC7", 0);
    }
  }

  // Load acceptance OADB
  if(fNameAcceptanceOADB.Length() && fTriggerSelection){
    AliDebugStream(1) << GetName() << ": Loading acceptance map from OADB file " <<  fNameAcceptanceOADB << std::endl;
    AliOADBContainer acceptanceCont("AliEmcalTriggerAcceptance");
    acceptanceCont.InitFromFile(fNameAcceptanceOADB.Data(), "AliEmcalTriggerAcceptance");
    TObjArray *acceptanceMaps = dynamic_cast<TObjArray *>(acceptanceCont.GetObject(fInputEvent->GetRunNumber()));
    TH2 *map(nullptr);
    if((map = dynamic_cast<TH2 *>(acceptanceMaps->FindObject("EG1")))){
      AliDebugStream(1) << GetName() << ": Found acceptance map for trigger EG1" << std::endl;
      map->SetDirectory(nullptr);
      fTriggerSelection->SetAcceptanceMap(AliEmcalTriggerOfflineSelection::kTrgEG1, map);
    }
    if((map = dynamic_cast<TH2 *>(acceptanceMaps->FindObject("EG2")))){
      AliDebugStream(1) << GetName() << ": Found acceptance map for trigger EG2" << std::endl;
      map->SetDirectory(nullptr);
      fTriggerSelection->SetAcceptanceMap(AliEmcalTriggerOfflineSelection::kTrgEG2, map);
    }
    if((map = dynamic_cast<TH2 *>(acceptanceMaps->FindObject("DG1")))){
      AliDebugStream(1) << GetName() << ": Found acceptance map for trigger DG1" << std::endl;
      map->SetDirectory(nullptr);
      fTriggerSelection->SetAcceptanceMap(AliEmcalTriggerOfflineSelection::kTrgDG1, map);
    }
    if((map = dynamic_cast<TH2 *>(acceptanceMaps->FindObject("DG2")))){
      AliDebugStream(1) << GetName() << ": Found acceptance map for trigger DG2" << std::endl;
      map->SetDirectory(nullptr);
      fTriggerSelection->SetAcceptanceMap(AliEmcalTriggerOfflineSelection::kTrgDG1, map);
    }
    if((map = dynamic_cast<TH2 *>(acceptanceMaps->FindObject("EJ1")))){
      AliDebugStream(1) << GetName() << ": Found acceptance map for trigger EJ1" << std::endl;
      map->SetDirectory(nullptr);
      fTriggerSelection->SetAcceptanceMap(AliEmcalTriggerOfflineSelection::kTrgEJ1, map);
    }
    if((map = dynamic_cast<TH2 *>(acceptanceMaps->FindObject("EJ2")))){
      AliDebugStream(1) << GetName() << ": Found acceptance map for trigger EJ2" << std::endl;
      map->SetDirectory(nullptr);
      fTriggerSelection->SetAcceptanceMap(AliEmcalTriggerOfflineSelection::kTrgEJ2, map);
    }
    if((map = dynamic_cast<TH2 *>(acceptanceMaps->FindObject("DJ1")))){
      AliDebugStream(1) << GetName() << ": Found acceptance map for trigger DJ1" << std::endl;
      map->SetDirectory(nullptr);
      fTriggerSelection->SetAcceptanceMap(AliEmcalTriggerOfflineSelection::kTrgDJ1, map);
    }
    if((map = dynamic_cast<TH2 *>(acceptanceMaps->FindObject("DJ2")))){
      AliDebugStream(1) << GetName() << ": Found acceptance map for trigger DJ2" << std::endl;
      map->SetDirectory(nullptr);
      fTriggerSelection->SetAcceptanceMap(AliEmcalTriggerOfflineSelection::kTrgDJ1, map);
    }
  }
}

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
    "MB", "EMC7", "EJ1", "EJ2", "EG1", "EG2", "EMC7excl", "EG2excl", "EJ2excl", "EJ1excl", "EG1excl"
  };
  std::vector<TString> dcaltriggers = {
      "DMC7", "DJ1", "DJ2", "DG1", "DG2", "DMC7excl", "DG2excl", "DJ2excl", "DJ1excl", "DG1excl"
  };
  if(fEnableDCALTriggers){
    for(const auto &t : dcaltriggers) triggers.push_back(t);
  }
  return triggers;
}

bool AliAnalysisTaskEmcalTriggerBase::SelectOnlineTrigger(AliEmcalTriggerOfflineSelection::EmcalTriggerClass trigger) const{
  AliDebugStream(1) << "Using V1 online trigger selector" << std::endl;
  if(trigger == AliEmcalTriggerOfflineSelection::kTrgEL0 || trigger == AliEmcalTriggerOfflineSelection::kTrgDL0) return true;
  int ngood(0);
  for(auto p : *fTriggerPatchInfo){
    AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(p);
    if((AliEmcalTriggerOfflineSelection::IsSingleShower(trigger) && patch->IsGammaLowRecalc()) ||
        (!AliEmcalTriggerOfflineSelection::IsSingleShower(trigger) && patch->IsJetLowRecalc())){
      if(patch->GetADCAmp() >= GetOnlineTriggerThresholdByIndex(trigger)) ngood++;
    }
  }
  return ngood > 0;
}

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

void AliAnalysisTaskEmcalTriggerBase::SetOnlineTriggerThreshold(const TString &triggerclass, Int_t threshold){
  TParameter<int> *threshpar(nullptr);
  if((threshpar = static_cast<TParameter<int> *>(fOnlineTriggerThresholds.FindObject(triggerclass.Data())))){
    threshpar->SetVal(threshold);
  } else {
    fOnlineTriggerThresholds.Add(new TParameter<int>(triggerclass.Data(), threshold));
  }
}

Int_t AliAnalysisTaskEmcalTriggerBase::GetOnlineTriggerThresholdByName(const TString &name) const {
  Int_t threshold(0);
  TParameter<int> *val(nullptr);
  if((val = static_cast<TParameter<int> *>(fOnlineTriggerThresholds.FindObject(name))))
    threshold = val->GetVal();
  return threshold;
}

Int_t AliAnalysisTaskEmcalTriggerBase::GetOnlineTriggerThresholdByIndex(AliEmcalTriggerOfflineSelection::EmcalTriggerClass trigger) const {
  const TString triggernames[AliEmcalTriggerOfflineSelection::kTrgn] = {"EMC7", "EG1", "EG2", "EJ1", "EJ2", "DMC7",
      "DG1", "DG2", "DJ1", "DJ2"};
  return GetOnlineTriggerThresholdByName(triggernames[trigger]);
}

Bool_t AliAnalysisTaskEmcalTriggerBase::OnlineThresholdsInitialized() const {
  const TString triggernames[AliEmcalTriggerOfflineSelection::kTrgn] = {"EMC7", "EG1", "EG2", "EJ1", "EJ2", "DMC7",
      "DG1", "DG2", "DJ1", "DJ2"};
  bool isInitialized = true;
  for(int itrg = 0; itrg < AliEmcalTriggerOfflineSelection::kTrgn; itrg++){
    if(!fOnlineTriggerThresholds.FindObject(triggernames[itrg].Data())) {
      isInitialized = false;
      break;
    }
  }
  return isInitialized;
}

Bool_t AliAnalysisTaskEmcalTriggerBase::SelectFiredPatch(const TString &triggerclass, Int_t adc) const {
  return adc > GetOnlineTriggerThresholdByName(triggerclass);
}

TString AliAnalysisTaskEmcalTriggerBase::GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const {
  TString triggerstring = "";
  Int_t nEJ1 = 0, nEJ2 = 0, nEG1 = 0, nEG2 = 0, nDJ1 = 0, nDJ2 = 0, nDG1 = 0, nDG2 = 0;
  for(auto patchIter : *triggerpatches){
    AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(patchIter);
    if(!patch->IsRecalc()) continue;
    if(patch->IsEMCal()){
      if(patch->IsGammaLowRecalc() && SelectFiredPatch("EG1", patch->GetADCAmp())) nEG1++;
      if(patch->IsGammaLowRecalc() && SelectFiredPatch("EG2", patch->GetADCAmp())) nEG2++;
      if(patch->IsJetLowRecalc() && SelectFiredPatch("EJ1", patch->GetADCAmp())) nEJ1++;
      if(patch->IsJetLowRecalc() && SelectFiredPatch("EJ2", patch->GetADCAmp())) nEJ2++;
    } else {
      if(patch->IsGammaLowRecalc() && SelectFiredPatch("DG1", patch->GetADCAmp())) nDG1++;
      if(patch->IsGammaLowRecalc() && SelectFiredPatch("DG2", patch->GetADCAmp())) nDG2++;
      if(patch->IsJetLowRecalc() && SelectFiredPatch("DJ1", patch->GetADCAmp())) nDJ1++;
      if(patch->IsJetLowRecalc() && SelectFiredPatch("DJ2", patch->GetADCAmp())) nDJ2++;
    }
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
  if(nDJ1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DJ2";
  }
  if(nDJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DJ2";
  }
  if(nDG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DG1";
  }
  if(nDG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DG2";
  }
  return triggerstring;
}

} /* namespace EMCalTriggerPtAnalysis */
