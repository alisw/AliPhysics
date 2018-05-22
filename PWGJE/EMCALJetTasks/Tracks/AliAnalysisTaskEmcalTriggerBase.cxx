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
#include <array>
#include <cfloat>
#include <functional>
#include <iostream>
#include <memory>

#include <TClonesArray.h>
#include <TGrid.h>
#include <THistManager.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParameter.h>

#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerMapping.h"
#include "AliEmcalTriggerDecisionContainer.h"
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
  fUseDownscaleCorrectionFormOCDB(kFALSE),
  fHistos(nullptr),
  fTriggerStringFromPatches(kFALSE),
  fSelectedTriggers(),
  fNameClusterContainer(""),
  fRequireAnalysisUtils(kTRUE),
  fUseSPDVertex(false),
  fApplyVertexCuts(true),
  fVertexCut(-10., 10.),
  fNameDownscaleOADB(""),
  fDownscaleOADB(nullptr),
  fDownscaleFactors(nullptr),
  fNameMaskedFastorOADB(),
  fMaskedFastorOADB(nullptr),
  fMaskedFastors(),
  fOnlineTriggerThresholds(),
  fNameAcceptanceOADB(),
  fNameTriggerSelectionContainer("EmcalTriggerDecision"),
  fSelectNoiseEvents(false),
  fRejectNoiseEvents(false),
  fEnableDCALTriggers(true),
  fEnableCentralityTriggers(false),
  fEnableV0Triggers(true),
  fEnableT0Triggers(false),
  fEnableNoINTTriggers(false),
  fRequireL0forL1(false),
  fExclusiveMinBias(false),
  fUseTriggerSelectionContainer(false)
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
  fUseDownscaleCorrectionFormOCDB(kFALSE),
  fHistos(nullptr),
  fTriggerStringFromPatches(kFALSE),
  fSelectedTriggers(),
  fNameClusterContainer(""),
  fRequireAnalysisUtils(kTRUE),
  fUseSPDVertex(false),
  fApplyVertexCuts(true),
  fVertexCut(-10., 10.),
  fNameDownscaleOADB(""),
  fDownscaleOADB(nullptr),
  fDownscaleFactors(nullptr),
  fNameMaskedFastorOADB(),
  fMaskedFastorOADB(nullptr),
  fMaskedFastors(),
  fOnlineTriggerThresholds(),
  fNameAcceptanceOADB(),
  fNameTriggerSelectionContainer("EmcalTriggerDecision"),
  fSelectNoiseEvents(false),
  fRejectNoiseEvents(false),
  fEnableDCALTriggers(true),
  fEnableCentralityTriggers(false),
  fEnableV0Triggers(true),
  fEnableT0Triggers(false),
  fEnableNoINTTriggers(false),
  fRequireL0forL1(false),
  fExclusiveMinBias(false),
  fUseTriggerSelectionContainer(false)
{
  SetNeedEmcalGeom(true);
  SetMakeGeneralHistograms(kTRUE);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalTriggerBase::~AliAnalysisTaskEmcalTriggerBase() {
  if(fTriggerSelection) delete fTriggerSelection;
  if(fHistos) delete fHistos;
  if(fDownscaleOADB) delete fDownscaleOADB;
  if(fMaskedFastorOADB) delete fMaskedFastorOADB;
  if(fDownscaleFactors) delete fDownscaleFactors;
}

void AliAnalysisTaskEmcalTriggerBase::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  if(fRequireAnalysisUtils && !fAliAnalysisUtils) fAliAnalysisUtils = new AliAnalysisUtils;

  if((!fNameClusterContainer.Length()) || fNameClusterContainer == "usedefault") fNameClusterContainer = AliEmcalAnalysisFactory::ClusterContainerNameFactory(fInputHandler->IsA() == AliAODInputHandler::Class());
  if(fTriggerSelection && !fTriggerSelection->GetNameClusterContainer().Length()){
    fTriggerSelection->SetClusterContainer(fNameClusterContainer);
  }

  fHistos = new THistManager(Form("Histos_%s", GetName()));

  // Create trigger correlation histogram
  const std::array<const TString, 11> binlabels = {"MB", "EMC7", "EG1", "EG2", "EJ1", "EJ2", "DMC7", "DG1", "DG2", "DJ1", "DJ2"};
  fHistos->CreateTH2("hTriggerCorrelation", "Correlation selected trigger classes", binlabels.size(), -0.5, binlabels.size() - 0.5, binlabels.size(), -0.5, binlabels.size() - 0.5);
  TH1 *correlationHist = static_cast<TH1 *>(fHistos->FindObject("hTriggerCorrelation"));
  for(decltype(binlabels.size()) ib = 0; ib < binlabels.size(); ib++){
    correlationHist->GetXaxis()->SetBinLabel(ib+1, binlabels[ib]);
    correlationHist->GetYaxis()->SetBinLabel(ib+1, binlabels[ib]);
  }

  CreateUserObjects();
  CreateUserHistos();

  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);
  fHistos->GetListOfHistograms()->SetOwner(false);

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalTriggerBase::IsEventSelected(){
  // Apply trigger selection
  TriggerSelection();
  if(!fSelectedTriggers.size()) {
    AliDebugStream(1) << "Failed trigger selection" << std::endl;
    return false;
  }

  UserFillHistosBeforeEventSelection();

  if(fApplyVertexCuts){
    const AliVVertex *vtx = fUseSPDVertex ? fInputEvent->GetPrimaryVertexSPD() : fInputEvent->GetPrimaryVertex();
    //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
    if(vtx->GetNContributors() < 1) {
      AliDebugStream(1) << "Failed Vertex Selection" << std::endl;
      return false;
    }
    if(!fVertexCut.IsInRange(fUseSPDVertex ? fVertexSPD[2] : fVertex[2])) {
      AliDebugStream(1) << "Failed vertex-z cut" << std::endl;
      return false;
    }
  }

  if(fRequireAnalysisUtils){
    AliDebugStream(1) << "Checking cuts in AliAnalysisUtils" << std::endl;
    if(fInputEvent->IsA() == AliESDEvent::Class() && fAliAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return false;
    if(fAliAnalysisUtils->IsPileUpEvent(fInputEvent)) return false;       // Apply new vertex cut
    if(!fAliAnalysisUtils->IsVertexSelected2013pA(fInputEvent))return false;       // Apply new vertex cut
  }

  if(!IsUserEventSelected()) {
    AliDebugStream(1) << "Failed user extra cuts" << std::endl;
    return false;
  }

  // Do MC outlier cut
  if(fIsPythia){
    if(!CheckMCOutliers()){
      AliDebugStream(1) << GetName() << ": Reject MC outliers" << std::endl;
      return false;
    }
  }

  // Fill histogram with trigger correlation
  // self-correlations included
  const std::array<const TString, 11> kAbsTriggers = {"MB", "EMC7", "EG1", "EG2", "EJ1", "EJ2", "DMC7", "DG1", "DG2", "DJ1", "DJ2"};
  for(int itrg = 0; itrg < kAbsTriggers.size(); itrg++){
    bool hasTriggerA = (std::find(fSelectedTriggers.begin(), fSelectedTriggers.end(), kAbsTriggers[itrg]) != fSelectedTriggers.end());
    if(hasTriggerA) {
      for(int jtrg = 0; jtrg < kAbsTriggers.size(); jtrg++){
        bool hasTriggerB = (std::find(fSelectedTriggers.begin(), fSelectedTriggers.end(), kAbsTriggers[jtrg]) != fSelectedTriggers.end());
        if(hasTriggerB)
          fHistos->FillTH2("hTriggerCorrelation", kAbsTriggers[itrg], kAbsTriggers[jtrg]);
      }
    }
  }

  UserFillHistosAfterEventSelection();
  AliDebugStream(1) << "Event is selected" << std::endl;
  return true;
}


void AliAnalysisTaskEmcalTriggerBase::TriggerSelection(){
  AliDebugStream(1) << "Entering trigger selection\n";
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
         isMinBiasT0 = selectionstatus & AliVEvent::kINT8,
         isCENT = selectionstatus & AliVEvent::kCentral,
         isSemiCENT = selectionstatus & AliVEvent::kSemiCentral,
      emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgn],
      emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgn],
      emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgn];

  if(fExclusiveMinBias){
    AliDebugStream(1) << "Min bias mode\n";
    // do not perform EMCAL trigger selection in case only
    // min. bias trigger is requested:w
    if(isMinBias) fSelectedTriggers.push_back("MB");
    if(isMinBiasT0 && fEnableT0Triggers) fSelectedTriggers.push_back("MBT0");
    if(fEnableCentralityTriggers){
      if(isCENT) fSelectedTriggers.push_back("CENT");
      if(isSemiCENT) fSelectedTriggers.push_back("SEMICENT");
    }
    return;
  }

  PWG::EMCAL::AliEmcalTriggerDecisionContainer *triggersel(nullptr);
  if(fUseTriggerSelectionContainer) {
    triggersel = dynamic_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->GetList()->FindObject(fNameTriggerSelectionContainer.Data()));
    if(!triggersel) {
      AliErrorStream() << "Trigger selection container requested but not found - not possible to select EMCAL triggers" << std::endl;
      return;
    }
  }

  AliDebugStream(1) << "Found triggers " << fInputEvent->GetFiredTriggerClasses() << std::endl;

  for(int itrg = 0; itrg < AliEmcalTriggerOfflineSelection::kTrgn; itrg++) emcalTriggers[itrg] = true;
  if(fEnableT0Triggers) for(int itrg = 0; itrg < AliEmcalTriggerOfflineSelection::kTrgn; itrg++) emc8Triggers[itrg] = true;
  if(fEnableNoINTTriggers) for(int itrg = 0; itrg < AliEmcalTriggerOfflineSelection::kTrgn; itrg++) emcNoIntTriggers[itrg] = true;
  if(!isMC){
    // In case of data select events as bunch-bunch (-B-) events.
    // Cut not applied in simulations
    if(fRequireBunchCrossing && ! (triggerstring.Contains("-B-") || triggerstring.Contains("-S-"))) return;

    // In case of data use information from the physics selection and event record
    // Further cleanup of trigger events can be performed depending on the presence
    // of recalc patches (after masking hot fastors in the trigger maker) above
    // threshold
    if(fUseTriggerBits){
      AliDebugStream(1) << "Require trigger bits" << std::endl;
      const std::array<ULong_t, AliEmcalTriggerOfflineSelection::kTrgn> kSelectTriggerBits = {
    	  AliVEvent::kEMC7|AliVEvent::kEMC8, AliVEvent::kEMCEGA, AliVEvent::kEMCEGA, AliVEvent::kEMCEJE, AliVEvent::kEMCEJE,
		  AliVEvent::kEMC7|AliVEvent::kEMC8, AliVEvent::kEMCEGA, AliVEvent::kEMCEGA, AliVEvent::kEMCEJE, AliVEvent::kEMCEJE
      };
      for(int iclass = 0; iclass < AliEmcalTriggerOfflineSelection::kTrgn; iclass++){
        if(!(selectionstatus & kSelectTriggerBits[iclass])) {
          emcNoIntTriggers[iclass] = emc8Triggers[iclass] = emcalTriggers[iclass] = false;
        }
        if(fRequireL0forL1 && !bool(selectionstatus & (AliVEvent::kEMC7|AliVEvent::kEMC8))) {
          emcNoIntTriggers[iclass] = emc8Triggers[iclass] = emcalTriggers[iclass] = false;
        }
      }
    }
    
    // Apply cut on the trigger string - this basically discriminates high- and low-threshold
    // triggers
    const std::array<TString, AliEmcalTriggerOfflineSelection::kTrgn> kSelectTriggerStrings = {
    		"CEMC7-|CEMC8-|C0EMC-", "EG1|EGA", "EG2", "EJ1|EJE", "EJ2", "CDMC7-|CDMC8-|C0DMC-", "DG1", "DG2", "DJ1", "DJ2"
    };
    if(triggerstring.Contains("EMC")) AliDebugStream(1) << GetName() << ": Trigger string " << triggerstring << std::endl;
    bool isT0trigger = triggerstring.Contains("INT8") || triggerstring.Contains("TVX") || triggerstring.Contains("EMC8") || triggerstring.Contains("DMC8"),
         isVZEROtrigger = triggerstring.Contains("INT7") || triggerstring.Contains("EMC7") || triggerstring.Contains("DMC7");
    for(int iclass = 0; iclass < AliEmcalTriggerOfflineSelection::kTrgn; iclass++){
      AliDebugStream(1) << "Next trigger: " << kSelectTriggerStrings[iclass] << std::endl;
      bool selectionStatus = false;
      if(kSelectTriggerStrings[iclass].Contains("|")){
        std::unique_ptr<TObjArray> options(kSelectTriggerStrings[iclass].Tokenize("|"));
        for(auto o : *options){
          TObjString *optstring = static_cast<TObjString *>(o);
          AliDebugStream(1) << "Checking trigger " << optstring->String() << std::endl;
          if(triggerstring.Contains(optstring->String())){
            AliDebugStream(1) << "Found " << optstring->String() << " ... " << std::endl;
            selectionStatus = true;
            if(fUseTriggerSelectionContainer){
              AliDebugStream(1) << "Checking trigger patch container for additional online patch" << std::endl;
              selectionStatus = selectionStatus && triggersel->IsEventSelected(optstring->String().Data());
            }
          }
        }
      } else {
        selectionStatus = triggerstring.Contains(kSelectTriggerStrings[iclass]);
        if(selectionStatus){
          AliDebugStream(1) << "Found " << kSelectTriggerStrings[iclass] << " ... " << std::endl;
          
        }
        if(fUseTriggerSelectionContainer){
          AliDebugStream(1) << "Checking trigger patch container for additional online patch" << std::endl;
          selectionStatus = selectionStatus && triggersel->IsEventSelected(kSelectTriggerStrings[iclass]);
        }
      }
      if(isT0trigger) {
        emc8Triggers[iclass] &= selectionStatus;
        emcalTriggers[iclass] = false;
        emcNoIntTriggers[iclass] = false; 
      } else if(isVZEROtrigger){
        emcalTriggers[iclass] &= selectionStatus;
        emc8Triggers[iclass] = false;
        emcNoIntTriggers[iclass] = false; 
      } else {
        // No coincidence with interaction trigger
        emcNoIntTriggers[iclass] &= selectionStatus; 
        emcalTriggers[iclass] = false;
        emc8Triggers[iclass] = false;
      }
      if(emcalTriggers[iclass])
        AliDebugStream(1) << GetName() << ": Event selected as trigger " << kSelectTriggerStrings[iclass] << " (INT7 suite)" << std::endl;
      if(emc8Triggers[iclass])
        AliDebugStream(1) << GetName() << ": Event selected as trigger " << kSelectTriggerStrings[iclass] << " (INT8 suite)" << std::endl;
      if(emcNoIntTriggers[iclass])
        AliDebugStream(1) << GetName() << ": Event selected as trigger " << kSelectTriggerStrings[iclass] << " (No INT coincidence)" << std::endl;
    }

    // Online selection / rejection
    if(fRejectNoiseEvents || fSelectNoiseEvents){
      for(int itrg = 0; itrg < AliEmcalTriggerOfflineSelection::kTrgn; itrg++){
        if(emcalTriggers[itrg] || emc8Triggers[itrg]){
          Bool_t onlinestatus = SelectOnlineTrigger(AliEmcalTriggerOfflineSelection::EmcalTriggerClass(itrg));;
          if(fRejectNoiseEvents){
            if(emcalTriggers[itrg]) emcalTriggers[itrg] &= onlinestatus;
            if(emc8Triggers[itrg]) emc8Triggers[itrg] &= onlinestatus;
          } else {
            if(emcalTriggers[itrg]) emcalTriggers[itrg] &= !onlinestatus;
            if(emc8Triggers[itrg]) emc8Triggers[itrg] &= !onlinestatus;
          }
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
  if(fEnableV0Triggers){
    if(isMinBias) fSelectedTriggers.push_back("MB");
    if(fEnableCentralityTriggers){
      if(isCENT) fSelectedTriggers.push_back("CENT");
      if(isSemiCENT) fSelectedTriggers.push_back("SEMICENT");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0]){
      AliDebugStream(1) << "Event selected as EMC7" << std::endl;
      fSelectedTriggers.push_back("EMC7");
      if(!isMinBias) fSelectedTriggers.push_back("EMC7excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2]){
      AliDebugStream(1) << "Event selected as EJ2" << std::endl;
      fSelectedTriggers.push_back("EJ2");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EJ2excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ1]){
      AliDebugStream(1) << "Event selected as EJ1" << std::endl;
      fSelectedTriggers.push_back("EJ1");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2])) fSelectedTriggers.push_back("EJ1excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2]){
      AliDebugStream(1) << "Event selected as EG2" << std::endl;
      fSelectedTriggers.push_back("EG2");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EG2excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEG1]){
      AliDebugStream(1) << "Event selected as EG1" << std::endl;
      fSelectedTriggers.push_back("EG1");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2])) fSelectedTriggers.push_back("EG1excl");
    }

    if(fEnableDCALTriggers){
      // Handle DCAL triggers only in case DCAL triggers are enabled,
      // otherwise ignore results of the online/offline trigger selection
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0]){
        AliDebugStream(1) << "Event selected as DMC7" << std::endl;
        fSelectedTriggers.push_back("DMC7");
        if(!isMinBias) fSelectedTriggers.push_back("DMC7excl");
      }
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2]){
        AliDebugStream(1) << "Event selected as DJ2" << std::endl;
        fSelectedTriggers.push_back("DJ2");
        if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DJ2excl");
      }
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ1]){
        AliDebugStream(1) << "Event selected as DJ1" << std::endl;
        fSelectedTriggers.push_back("DJ1");
        if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2])) fSelectedTriggers.push_back("DJ1excl");
      }
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDG2]){
        AliDebugStream(1) << "Event selected as DG2" << std::endl;
        fSelectedTriggers.push_back("DG2");
        if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DG2excl");
      }
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDG1]){
        AliDebugStream(1) << "Event selected as DG1" << std::endl;
        fSelectedTriggers.push_back("DG1");
        if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDG2])) fSelectedTriggers.push_back("DG1excl");
      }
    }
  }

  if(fEnableT0Triggers) {
    if(isMinBiasT0) fSelectedTriggers.push_back("MBT0");
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEL0]) {
      // EMC8 trigger
      AliDebugStream(1) << "Event selected as EMC8" << std::endl;
      fSelectedTriggers.push_back("EMC8");
      if(!isMinBiasT0) fSelectedTriggers.push_back("EMC8excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEJ2]){
      AliDebugStream(1) << "Event selected as EJ2 (EMC8)" << std::endl;
      fSelectedTriggers.push_back("EMC8EJ2");
      if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EMC8EJ2excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEJ1]){
      AliDebugStream(1) << "Event selected as EJ1 (EMC8)" << std::endl;
      fSelectedTriggers.push_back("EMC8EJ1");
      if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEJ2])) fSelectedTriggers.push_back("EMC8EJ1excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEG2]){
      AliDebugStream(1) << "Event selected as EG2 (EMC8)" << std::endl;
      fSelectedTriggers.push_back("EMC8EG2");
      if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EMC8EG2excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEG1]){
      AliDebugStream(1) << "Event selected as EG1 (EMC8)" << std::endl;
      fSelectedTriggers.push_back("EMC8EG1");
      if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEG2])) fSelectedTriggers.push_back("EMC8EG1excl");
    }

    if(fEnableDCALTriggers){
      // Handle DCAL triggers only in case DCAL triggers are enabled,
      // otherwise ignore results of the online/offline trigger selection
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0]){
        AliDebugStream(1) << "Event selected as DMC8" << std::endl;
        fSelectedTriggers.push_back("DMC8");
        if(!isMinBiasT0) fSelectedTriggers.push_back("DMC8excl");
      }
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDJ2]){
        AliDebugStream(1) << "Event selected as DJ2 (DMC8)" << std::endl;
        fSelectedTriggers.push_back("DMC8DJ2");
        if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DMC8DJ2excl");
      }
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDJ1]){
        AliDebugStream(1) << "Event selected as DJ1 (DMC8)" << std::endl;
        fSelectedTriggers.push_back("DMC8DJ1");
        if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDJ2])) fSelectedTriggers.push_back("DMC8DJ1excl");
      }
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDG2]){
        AliDebugStream(1) << "Event selected as DG2 (DMC8)" << std::endl;
        fSelectedTriggers.push_back("DMC8DG2");
        if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DMC8DG2excl");
      }
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDG1]){
        AliDebugStream(1) << "Event selected as DG1 (DMC8)" << std::endl;
        fSelectedTriggers.push_back("DMC8DG1");
        if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDG2])) fSelectedTriggers.push_back("DMC8DG1excl");
      }
    }
  }

  if(fEnableNoINTTriggers) {
    if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0]) {
      // EMC8 trigger
      AliDebugStream(1) << "Event selected as 0EMC" << std::endl;
      fSelectedTriggers.push_back("0EMC");
    }
    if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2]){
      AliDebugStream(1) << "Event selected as EJ2 (E0MC)" << std::endl;
      fSelectedTriggers.push_back("0EMCEJ2");
      if(!emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0]) fSelectedTriggers.push_back("0EMCEJ2excl");
    }
    if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ1]){
      AliDebugStream(1) << "Event selected as EJ1 (0EMC)" << std::endl;
      fSelectedTriggers.push_back("0EMCEJ1");
      if(!(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2])) fSelectedTriggers.push_back("0EMCEJ1excl");
    }
    if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2]){
      AliDebugStream(1) << "Event selected as EG2 (0EMC)" << std::endl;
      fSelectedTriggers.push_back("0EMCEG2");
      if(!emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0]) fSelectedTriggers.push_back("0EMCEG2excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEG1]){
      AliDebugStream(1) << "Event selected as EG1 (EMC8)" << std::endl;
      fSelectedTriggers.push_back("0EMCEG1");
      if(!(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2])) fSelectedTriggers.push_back("0EMCEG1excl");
    }

    if(fEnableDCALTriggers){
      // Handle DCAL triggers only in case DCAL triggers are enabled,
      // otherwise ignore results of the online/offline trigger selection
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0]){
        AliDebugStream(1) << "Event selected as 0DMC" << std::endl;
        fSelectedTriggers.push_back("0DMC");
      }
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2]){
        AliDebugStream(1) << "Event selected as DJ2 (0DMC)" << std::endl;
        fSelectedTriggers.push_back("0DMCDJ2");
        if(!emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0]) fSelectedTriggers.push_back("0DMCDJ2excl");
      }
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ1]){
        AliDebugStream(1) << "Event selected as DJ1 (0DMC)" << std::endl;
        fSelectedTriggers.push_back("0DMCDJ1");
        if(!(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2])) fSelectedTriggers.push_back("0DMCDJ1excl");
      }
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDG2]){
        AliDebugStream(1) << "Event selected as DG2 (0DMC)" << std::endl;
        fSelectedTriggers.push_back("0DMCDG2");
        if(!emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0]) fSelectedTriggers.push_back("0DMCDG2excl");
      }
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDG1]){
        AliDebugStream(1) << "Event selected as DG1 (0DMC)" << std::endl;
        fSelectedTriggers.push_back("0DMCDG1");
        if(!(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDG2])) fSelectedTriggers.push_back("0DMCDG1excl");
      }
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

  if(!fExclusiveMinBias){
    // Load EMCAL trigger OADB in case EMCAL triggers
    // are enabled

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
}

void AliAnalysisTaskEmcalTriggerBase::RunChanged(Int_t runnumber){
  if(fDownscaleOADB){
    fDownscaleFactors = static_cast<TObjArray *>(fDownscaleOADB->GetObject(runnumber));
  }
  if(fUseDownscaleCorrectionFormOCDB) PrepareDownscaleFactorsFormOCDB();

  // Log info for downscale factors
  if(!fDownscaleFactors || !fDownscaleFactors->GetEntries()){
    AliInfoStream() << GetName() << ": No downscale factors provided for run " << runnumber << std::endl;
  } else {
    AliInfoStream() << GetName() << ": Downscale factors used for run " << runnumber << std::endl;
    for(auto e : *fDownscaleFactors){
      TParameter<double> *dfactor = static_cast<TParameter<double> *>(e);
      AliInfoStream() << GetName() << ": Trigger " << dfactor->GetName() << ", downscale factor " << dfactor->GetVal() << std::endl;
    }
  }

  if(!fExclusiveMinBias){
    if(fMaskedFastorOADB){
      fMaskedFastors.clear();
      TObjArray *ids = static_cast<TObjArray *>(fMaskedFastorOADB->GetObject(runnumber));
      for(auto m : *ids){
        TParameter<int> *id = static_cast<TParameter<int> *>(m);
        fMaskedFastors.push_back(id->GetVal());
      }
    }
  }
}

std::vector<TString> AliAnalysisTaskEmcalTriggerBase::GetSupportedTriggers(Bool_t useExclusiveTriggers) const {
  // Exclusive means classes without lower trigger classes (which are downscaled) -
  // in order to make samples statistically independent: MBExcl means MinBias && !EMCAL trigger
  std::vector<TString> triggers; 
  const std::array<TString, 5>  emcaltriggers = {{"EMC7", "EJ1", "EJ2", "EG1", "EG2"}},
                                dcaltriggers = {{"DMC7", "DJ1", "DJ2", "DG1", "DG2"}},
                                emcalexclusive = {{"EMC7excl", "EG2excl", "EJ2excl", "EJ1excl", "EG1excl"}},
                                dcalexclusive = {{"DMC7excl", "DG2excl", "DJ2excl", "DJ1excl", "DG1excl"}},
                                t0triggers = {{"EMC8",  "EMC8EJ1", "EMC8EJ2", "EMC8EG1", "EMC8EG2"}},
                                t0exclusive = {{"EMC8excl", "EMC8EG2excl", "EMC8EJ2excl", "EMC8EJ1excl", "EMC8EG1excl"}},
                                t0dcaltriggers = {{"DMC8", "DMC8DJ1", "DMC8DJ2", "DMC8DG1", "DMC8DG2"}},
                                t0dcalexclusive = {{"DMC8excl", "DMC8DG2excl", "DMC8DJ2excl", "DMC8DJ1excl", "DMC8DG1excl"}},
                                nointEMCAL = {{"0EMC", "0EMCEJ1", "0EMCEJ2", "0EMCEG1", "0EMCEG2"}},
                                nointDCAL = {{"0DMC", "0DMCDJ1", "0DMCDJ2", "0DMCDG1", "0DMCDG2"}};
  const std::array<TString, 4>  nointemcalexclusive = {{"0EMCEG2excl", "0EMCEJ2excl", "0EMCEJ1excl", "0EMCEG1excl"}}, 
                                nointdcalexclusive = {{"0DMCDG2excl", "0DMCDJ2excl", "0DMCDJ1excl", "0DMCDG1excl"}};
  const std::array<TString, 2> centralitytriggers = {{"CENT", "SEMICENT"}};
  if(fEnableV0Triggers){
    triggers.push_back("MB"); // Min. Bias always enabled
    if(!fExclusiveMinBias){
      for(const auto &t : emcaltriggers) triggers.push_back(t);
      if(useExclusiveTriggers)
        for(const auto &t : emcalexclusive) triggers.push_back(t);
    }
    if(fEnableDCALTriggers){
      for(const auto &t : dcaltriggers) triggers.push_back(t);
      if(useExclusiveTriggers)
        for(const auto &t : dcalexclusive) triggers.push_back(t);
    }
    if(fEnableCentralityTriggers){
      for(const auto &t : centralitytriggers) triggers.push_back(t);
    }
  }
  if(fEnableT0Triggers){
    triggers.push_back("MBT0");
    if(!fExclusiveMinBias){
      for(const auto &t: t0triggers) triggers.push_back(t);
      if(useExclusiveTriggers)
        for(const auto &t : t0exclusive) triggers.push_back(t);
    }
    if(fEnableDCALTriggers){
      for(const auto &t: t0dcaltriggers) triggers.push_back(t);
      if(useExclusiveTriggers)
        for(const auto &t : t0dcalexclusive) triggers.push_back(t);
    }
  }
  if(fEnableNoINTTriggers) { 
    // No MB trigger since no interaction trigger
    if(!fExclusiveMinBias){
      for(const auto &t: nointEMCAL) triggers.push_back(t);
      if(useExclusiveTriggers)
        for(const auto &t : nointemcalexclusive) triggers.push_back(t);
    }
    if(fEnableDCALTriggers){
      for(const auto &t: nointDCAL) triggers.push_back(t);
      if(useExclusiveTriggers)
        for(const auto &t : nointdcalexclusive) triggers.push_back(t);
    }
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
    else if(triggerclass.Contains("DMC7")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DMC7"));
    else if(triggerclass.Contains("EJ2")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ2"));
    else if(triggerclass.Contains("EJ1")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ1"));
    else if(triggerclass.Contains("EG2")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG2"));
    else if(triggerclass.Contains("EG1")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG1"));
    else if(triggerclass.Contains("DG2")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DG2"));
    else if(triggerclass.Contains("DG1")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DG1"));
    double triggerweight = 1.;
    if(result) triggerweight = 1./result->GetVal();
    AliDebugStream(1) << "Using trigger weight " << triggerweight << " for trigger " << triggerclass << std::endl;
    return triggerweight;
  } else {
    AliDebugStream(1) << "No downscale factors loaded - using trigger weight 1" << std::endl;
  }
  return 1.;
}

void AliAnalysisTaskEmcalTriggerBase::PrepareDownscaleFactorsFormOCDB(){
  AliInfoStream() << "Reading downscale factors from OCDB for run " << fInputEvent->GetRunNumber() << std::endl;
  if(!fDownscaleFactors){
    fDownscaleFactors = new TObjArray;
    fDownscaleFactors->SetOwner(true);
  }
  fDownscaleFactors->Clear();
  auto downscaleOCDB = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
  if(downscaleOCDB->GetCurrentRun() != fInputEvent->GetRunNumber()) downscaleOCDB->SetRun(fInputEvent->GetRunNumber());
  const std::array<TString, 11> khwtriggers = {"INT7", "EMC7", "DMC7", "EJ1", "EJ2", "DJ1", "DJ2", "EG1", "EG2", "DG1", "DG2"};
  std::vector<TString> runtriggers = downscaleOCDB->GetTriggerClasses();
  for(const auto &t : khwtriggers){
    std::function<bool (TString)> triggerfinder = [t](const TString &test) -> bool {
      if(!test.Contains(t + "-B-")) return false;
      return true;
    };
    auto entry = std::find_if(runtriggers.begin(), runtriggers.end(), triggerfinder);
    if(entry  != runtriggers.end()){
      TString triggername = *entry;
      double downscalefactor = downscaleOCDB->GetDownscaleFactorForTriggerClass(triggername);
      AliInfoStream() << "Applying downscale factor " << downscalefactor << " for trigger " << t << " (" << triggername << ") for run " << fInputEvent->GetRunNumber() << std::endl;
      fDownscaleFactors->Add(new TParameter<double>(t, TMath::Abs(downscalefactor) < DBL_EPSILON ? 1 : downscalefactor));
    } else {
      AliErrorStream() << "No downscale factor found for trigger " << t << " for run " << fInputEvent->GetRunNumber() << std::endl;
    }
  }
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
