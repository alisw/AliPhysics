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
#include <set>
#include <sstream>

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
#include "AliEmcalTriggerStringDecoder.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliOADBContainer.h"
#include "AliVVertex.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerBase)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalTriggerBase::AliAnalysisTaskEmcalTriggerBase():
  AliAnalysisTaskEmcal(),
  fHistos(nullptr),
  fUseTriggerBits(kTRUE),
  fRequireBunchCrossing(kTRUE),
  fUseDownscaleCorrectionFormOCDB(kFALSE),
  fTriggerSelection(nullptr),
  fSelectedTriggers(),
  fNameClusterContainer(""),
  fRequireAnalysisUtils(kTRUE),
  fUseSPDVertex(false),
  fApplyVertexCuts(true),
  fVertexCut(-10., 10.),
  fNameDownscaleOADB(""),
  fDownscaleOADB(nullptr),
  fDownscaleFactors(nullptr),
  fNameTriggerSelectionContainer("EmcalTriggerDecision"),
  fEnableDCALTriggers(true),
  fEnableEDCombinedTriggers(false),
  fEnableV0Triggers(true),
  fEnableT0Triggers(false),
  fEnableNoINTTriggers(false),
  fEnableCentralityTriggers(false),
  fExclusiveMinBias(false),
  fUseTriggerSelectionContainer(false),
  fSelectCentralityTriggers2018(false)
{
  SetNeedEmcalGeom(true);
  SetMakeGeneralHistograms(kTRUE);
}

AliAnalysisTaskEmcalTriggerBase::AliAnalysisTaskEmcalTriggerBase(const char *name):
  AliAnalysisTaskEmcal(name, true),
  fHistos(nullptr),
  fUseTriggerBits(kTRUE),
  fRequireBunchCrossing(kTRUE),
  fUseDownscaleCorrectionFormOCDB(kFALSE),
  fTriggerSelection(nullptr),
  fSelectedTriggers(),
  fNameClusterContainer(""),
  fRequireAnalysisUtils(kTRUE),
  fUseSPDVertex(false),
  fApplyVertexCuts(true),
  fVertexCut(-10., 10.),
  fNameDownscaleOADB(""),
  fDownscaleOADB(nullptr),
  fDownscaleFactors(nullptr),
  fNameTriggerSelectionContainer("EmcalTriggerDecision"),
  fEnableDCALTriggers(true),
  fEnableEDCombinedTriggers(false),
  fEnableV0Triggers(true),
  fEnableT0Triggers(false),
  fEnableNoINTTriggers(false),
  fEnableCentralityTriggers(false),
  fExclusiveMinBias(false),
  fUseTriggerSelectionContainer(false),
  fSelectCentralityTriggers2018(false)
{
  SetNeedEmcalGeom(true);
  SetMakeGeneralHistograms(kTRUE);
}

AliAnalysisTaskEmcalTriggerBase::~AliAnalysisTaskEmcalTriggerBase() {
  if(fTriggerSelection) delete fTriggerSelection;
  if(fHistos) delete fHistos;
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
  std::vector<std::string> binlabels = {"MB"};
  if(fEnableCentralityTriggers) {
    binlabels.emplace_back("CENT");
    binlabels.emplace_back("SEMICENT");
  }
  if(fEnableV0Triggers){
    const std::array<const std::string, 5> vzlabels = {{"EMC7", "EG1", "EG2", "EJ1", "EJ2"}};
    for(const auto & vlab : vzlabels) binlabels.emplace_back(vlab);
    if(fEnableDCALTriggers){
      const std::array<const std::string, 5> dclabels = {{"DMC7", "DG1", "DG2", "DJ1", "DJ2"}};
      for(const auto & dlab : dclabels) binlabels.emplace_back(dlab);
    }
  }
  if(fEnableT0Triggers) {
    binlabels.emplace_back("MBT0");
    const std::array<const std::string, 5> t0labels = {{"EMC8", "EMC8EG1", "EMC8EG2", "EMC8EJ1", "EMC8EJ2"}};
    for(const auto & tlab : t0labels) binlabels.emplace_back(tlab);
    if(fEnableDCALTriggers){
      const std::array<const std::string, 5> dtclabels = {{"DMC8", "DMC8DG1", "DMC8DG2", "DMC8DJ1", "DMC8DJ2"}};
      for(const auto & dtlab : dtclabels) binlabels.emplace_back(dtlab);
    }
  }
  fHistos->CreateTH2("hTriggerCorrelation", "Correlation selected trigger classes", binlabels.size(), -0.5, binlabels.size() - 0.5, binlabels.size(), -0.5, binlabels.size() - 0.5);
  TH1 *correlationHist = static_cast<TH1 *>(fHistos->FindObject("hTriggerCorrelation"));
  for(decltype(binlabels.size()) ib = 0; ib < binlabels.size(); ib++){
    correlationHist->GetXaxis()->SetBinLabel(ib+1, binlabels[ib].data());
    correlationHist->GetYaxis()->SetBinLabel(ib+1, binlabels[ib].data());
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

  if(!AliAnalysisTaskEmcal::IsEventSelected()) return false;

  if(!IsUserEventSelected()) {
    AliDebugStream(1) << "Failed user extra cuts" << std::endl;
    return false;
  }

  // Fill histogram with trigger correlation
  // self-correlations included
  auto *corrhist = static_cast<TH2 *>(fHistos->GetListOfHistograms()->FindObject("hTriggerCorrelation"));
  for(int itrg = 0; itrg < corrhist->GetXaxis()->GetNbins(); itrg++){
    const char *xlabel = corrhist->GetXaxis()->GetBinLabel(itrg+1);
    bool hasTriggerA = (std::find(fSelectedTriggers.begin(), fSelectedTriggers.end(), xlabel) != fSelectedTriggers.end());
    if(hasTriggerA) {
      for(int jtrg = 0; jtrg < corrhist->GetYaxis()->GetNbins(); jtrg++){
        const char *ylabel = corrhist->GetYaxis()->GetBinLabel(jtrg+1);
        bool hasTriggerB = (std::find(fSelectedTriggers.begin(), fSelectedTriggers.end(), ylabel) != fSelectedTriggers.end());
        if(hasTriggerB) fHistos->FillTH2("hTriggerCorrelation", xlabel, ylabel);
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

  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
         isMinBiasT0 = selectionstatus & AliVEvent::kINT8,
         isCENT = selectionstatus & AliVEvent::kCentral,
         isSemiCENT = selectionstatus & AliVEvent::kSemiCentral,
         emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgn],
         emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgn],
         emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgn];
  
  // check the centrality triggers
  // temp hack to overcome missing support for 2018 by the physics selection 
  if(fEnableCentralityTriggers) {
    if(fSelectCentralityTriggers2018) {
      auto triggers = PWG::EMCAL::Triggerinfo::DecodeTriggerString(fInputEvent->GetFiredTriggerClasses().Data());
      for(auto t : triggers) {
        if(t.Triggercluster() != "CENT") continue;
        if(t.Triggerclass() == "CV0H7") isCENT = true;
        else if(t.Triggerclass() == "CMID7") isSemiCENT = true;
      }
    }
  }

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
  const std::array<std::string, AliEmcalTriggerOfflineSelection::kTrgn> kEmcalSelectTriggerStrings = {
    		"=CEMC7|CEMC8|C0EMC|EMCL0", "EG1|EGA", "EG2", "EJ1|EJE", "EJ2", "=CDMC7|CDMC8|C0DMC|DMCL0", "DG1", "DG2", "DJ1", "DJ2"
  };
  if(!isMC){
    // In case of data select events as bunch-bunch (-B-) events.
    // Cut not applied in simulations
    if(fRequireBunchCrossing && ! (fInputEvent->GetFiredTriggerClasses().Contains("-B-") || fInputEvent->GetFiredTriggerClasses().Contains("-S-"))) return;

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
      }
    }
    auto triggerstring =  fInputEvent->GetFiredTriggerClasses();
    if(triggerstring.Contains("EMC") || triggerstring.Contains("DMC") || 
       triggerstring.Contains("INT7E") || triggerstring.Contains("INT7D")){   // special conditions for 2015 PbPb
      // Apply cut on the trigger string - this basically discriminates high- and low-threshold
      // triggers
      auto triggers = PWG::EMCAL::Triggerinfo::DecodeTriggerString(fInputEvent->GetFiredTriggerClasses().Data());
      std::map<int, std::array<bool, 3>> matchedTriggers;
      for(auto t : triggers) {  
        const auto &triggerclass = t.Triggerclass();
        if((triggerclass.find("EMC") != std::string::npos) || (triggerclass.find("DMC") != std::string::npos) || 
           (triggerclass.find("INT7E") != std::string::npos) || (triggerclass.find("INT7D") != std::string::npos)) 
          AliDebugStream(1) << GetName() << ": Trigger string " << t.ExpandClassName() << std::endl;
        else continue;    // No EMC / DMC trigger - not to be checked
        bool isT0trigger = (triggerclass.find("EMC8") != std::string::npos) || (triggerclass.find("DMC8") != std::string::npos) || (triggerclass.find("INT8") != std::string::npos),
             isVZEROtrigger = (triggerclass.find("EMC7") != std::string::npos) || (triggerclass.find("DMC7") != std::string::npos) || (triggerclass.find("INT7") != std::string::npos);
        for(auto iclass = 0; iclass < AliEmcalTriggerOfflineSelection::kTrgn; iclass++){
          AliDebugStream(1) << "Next trigger: " << kEmcalSelectTriggerStrings[iclass] << std::endl;
          bool emcalSelectionStatus = MatchTriggerFromPattern(kEmcalSelectTriggerStrings[iclass], triggerclass);
          if(fUseTriggerSelectionContainer) emcalSelectionStatus = emcalSelectionStatus & MatchTriggerFromContainer(kEmcalSelectTriggerStrings[iclass], triggersel);
          if(emcalSelectionStatus) {
            auto entry = matchedTriggers.find(iclass);
            if(entry == matchedTriggers.end()) {
              std::array<bool, 3> interactions = {{isVZEROtrigger, isT0trigger, !(isVZEROtrigger || isT0trigger)}};
              matchedTriggers.insert(std::pair<int, std::array<bool, 3>>(iclass, interactions));
            } else {
              auto &interactions = entry->second;
              if(isVZEROtrigger) interactions[0] = true;
              if(isT0trigger) interactions[1] = true;
              if(!(isVZEROtrigger || isT0trigger)) interactions[2] = true;
            }
          }
        }
      }

      for(auto iclass = 0; iclass < AliEmcalTriggerOfflineSelection::kTrgn; iclass++){
        auto entry = matchedTriggers.find(iclass);
        if(entry != matchedTriggers.end()){
          // trigger selected among at least one L0 class
          const auto &interactions = entry->second;
          emcalTriggers[iclass] &= interactions[0];
          emc8Triggers[iclass] &= interactions[1];
          emcNoIntTriggers[iclass] &= interactions[2]; 
        } else {
          // trigger not selected - mark all cases as false
          emcalTriggers[iclass] = false;
          emc8Triggers[iclass] = false;
          emcNoIntTriggers[iclass] = false; 
        }

        if(emcalTriggers[iclass])
          AliDebugStream(1) << GetName() << ": Event selected as trigger " << kEmcalSelectTriggerStrings[iclass] << " (INT7 suite)" << std::endl;
        if(emc8Triggers[iclass])
          AliDebugStream(1) << GetName() << ": Event selected as trigger " << kEmcalSelectTriggerStrings[iclass] << " (INT8 suite)" << std::endl;
        if(emcNoIntTriggers[iclass])
          AliDebugStream(1) << GetName() << ": Event selected as trigger " << kEmcalSelectTriggerStrings[iclass] << " (No INT coincidence)" << std::endl;
      }
    }
  } else {
    // MC: Use INT7/INT08 for VZERO/TZERO triggers, for EMCAL trigger use trigger selection container
    bool isVZEROtrigger = selectionstatus & AliVEvent::kINT7, isT0trigger = selectionstatus & AliVEvent::kINT8;
    if(fUseTriggerSelectionContainer){
      for(int iclass = 0; iclass < AliEmcalTriggerOfflineSelection::kTrgn; iclass++){
        auto emcalSelectionStatus = MatchTriggerFromContainer(kEmcalSelectTriggerStrings[iclass], triggersel);
        if(emcalSelectionStatus) {
          // trigger selected - correlate with interaction trigger status
          emcalTriggers[iclass] &= isVZEROtrigger;
          emc8Triggers[iclass] &= isT0trigger;
          emcNoIntTriggers[iclass] &= !(isT0trigger || isVZEROtrigger);
        } else {
          // trigger not selected, set trigger bit to false
          emcalTriggers[iclass] = false;
          emc8Triggers[iclass] = false;
          emcNoIntTriggers[iclass] = false;
        }
        if(emcalTriggers[iclass])
          AliDebugStream(1) << GetName() << ": Event selected as trigger " << kEmcalSelectTriggerStrings[iclass] << " (INT7 suite)" << std::endl;
        if(emc8Triggers[iclass])
          AliDebugStream(1) << GetName() << ": Event selected as trigger " << kEmcalSelectTriggerStrings[iclass] << " (INT8 suite)" << std::endl;
        if(emcNoIntTriggers[iclass])
          AliDebugStream(1) << GetName() << ": Event selected as trigger " << kEmcalSelectTriggerStrings[iclass] << " (No INT coincidence)" << std::endl;
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
  std::set<std::string> combinedtriggers;
  if(fEnableV0Triggers){
    if(isMinBias) fSelectedTriggers.push_back("MB");
    if(fEnableCentralityTriggers){
      if(isCENT) fSelectedTriggers.push_back("CENT");
      if(isSemiCENT) fSelectedTriggers.push_back("SEMICENT");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0]){
      AliDebugStream(1) << "Event selected as EMC7" << std::endl;
      fSelectedTriggers.push_back("EMC7");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC7") == combinedtriggers.end())) combinedtriggers.insert("EDMC7");
      if(!isMinBias) fSelectedTriggers.push_back("EMC7excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2]){
      AliDebugStream(1) << "Event selected as EJ2" << std::endl;
      fSelectedTriggers.push_back("EJ2");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDJ2") == combinedtriggers.end())) combinedtriggers.insert("EDJ2");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EJ2excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ1]){
      AliDebugStream(1) << "Event selected as EJ1" << std::endl;
      fSelectedTriggers.push_back("EJ1");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDJ1") == combinedtriggers.end())) combinedtriggers.insert("EDJ1");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2])) fSelectedTriggers.push_back("EJ1excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2]){
      AliDebugStream(1) << "Event selected as EG2" << std::endl;
      fSelectedTriggers.push_back("EG2");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDG2") == combinedtriggers.end())) combinedtriggers.insert("EDG2");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EG2excl");
    }
    if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEG1]){
      AliDebugStream(1) << "Event selected as EG1" << std::endl;
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDG1") == combinedtriggers.end()))  combinedtriggers.insert("EDG1");
      fSelectedTriggers.push_back("EG1");
      if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2])) fSelectedTriggers.push_back("EG1excl");
    }

    if(fEnableDCALTriggers){
      // Handle DCAL triggers only in case DCAL triggers are enabled,
      // otherwise ignore results of the online/offline trigger selection
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0]){
        AliDebugStream(1) << "Event selected as DMC7" << std::endl;
        fSelectedTriggers.push_back("DMC7");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC7") == combinedtriggers.end())) combinedtriggers.insert("EDMC7");
        if(!isMinBias) fSelectedTriggers.push_back("DMC7excl");
      }
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2]){
        AliDebugStream(1) << "Event selected as DJ2" << std::endl;
        fSelectedTriggers.push_back("DJ2");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDJ2") == combinedtriggers.end())) combinedtriggers.insert("EDJ2");
        if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DJ2excl");
      }
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ1]){
        AliDebugStream(1) << "Event selected as DJ1" << std::endl;
        fSelectedTriggers.push_back("DJ1");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDJ1") == combinedtriggers.end())) combinedtriggers.insert("EDJ1");
        if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2])) fSelectedTriggers.push_back("DJ1excl");
      }
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDG2]){
        AliDebugStream(1) << "Event selected as DG2" << std::endl;
        fSelectedTriggers.push_back("DG2");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDG2") == combinedtriggers.end())) combinedtriggers.insert("EDG2");
        if(!(isMinBias || emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DG2excl");
      }
      if(emcalTriggers[AliEmcalTriggerOfflineSelection::kTrgDG1]){
        AliDebugStream(1) << "Event selected as DG1" << std::endl;
        fSelectedTriggers.push_back("DG1");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDG1") == combinedtriggers.end()))  combinedtriggers.insert("EDG1");
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
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8") == combinedtriggers.end())) combinedtriggers.insert("EDMC8");
      if(!isMinBiasT0) fSelectedTriggers.push_back("EMC8excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEJ2]){
      AliDebugStream(1) << "Event selected as EJ2 (EMC8)" << std::endl;
      fSelectedTriggers.push_back("EMC8EJ2");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8J2") == combinedtriggers.end())) combinedtriggers.insert("EDMC8J2");
      if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EMC8EJ2excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEJ1]){
      AliDebugStream(1) << "Event selected as EJ1 (EMC8)" << std::endl;
      fSelectedTriggers.push_back("EMC8EJ1");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8J1") == combinedtriggers.end())) combinedtriggers.insert("EDMC8J1");
      if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEJ2])) fSelectedTriggers.push_back("EMC8EJ1excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEG2]){
      AliDebugStream(1) << "Event selected as EG2 (EMC8)" << std::endl;
      fSelectedTriggers.push_back("EMC8EG2");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8G2") == combinedtriggers.end())) combinedtriggers.insert("EDMC8G2");
      if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEL0])) fSelectedTriggers.push_back("EMC8EG2excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEG1]){
      AliDebugStream(1) << "Event selected as EG1 (EMC8)" << std::endl;
      fSelectedTriggers.push_back("EMC8EG1");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8G1") == combinedtriggers.end())) combinedtriggers.insert("EDMC8G1");
      if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEG2])) fSelectedTriggers.push_back("EMC8EG1excl");
    }

    if(fEnableDCALTriggers){
      // Handle DCAL triggers only in case DCAL triggers are enabled,
      // otherwise ignore results of the online/offline trigger selection
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0]){
        AliDebugStream(1) << "Event selected as DMC8" << std::endl;
        fSelectedTriggers.push_back("DMC8");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8") == combinedtriggers.end())) combinedtriggers.insert("EDMC8");
        if(!isMinBiasT0) fSelectedTriggers.push_back("DMC8excl");
      }
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDJ2]){
        AliDebugStream(1) << "Event selected as DJ2 (DMC8)" << std::endl;
        fSelectedTriggers.push_back("DMC8DJ2");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8J1") == combinedtriggers.end())) combinedtriggers.insert("EDMC8J1");
        if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DMC8DJ2excl");
      }
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDJ1]){
        AliDebugStream(1) << "Event selected as DJ1 (DMC8)" << std::endl;
        fSelectedTriggers.push_back("DMC8DJ1");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8J1") == combinedtriggers.end())) combinedtriggers.insert("EDMC8J1");
        if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDJ2])) fSelectedTriggers.push_back("DMC8DJ1excl");
      }
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDG2]){
        AliDebugStream(1) << "Event selected as DG2 (DMC8)" << std::endl;
        fSelectedTriggers.push_back("DMC8DG2");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8G2") == combinedtriggers.end())) combinedtriggers.insert("EDMC8G2");
        if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0])) fSelectedTriggers.push_back("DMC8DG2excl");
      }
      if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDG1]){
        AliDebugStream(1) << "Event selected as DG1 (DMC8)" << std::endl;
        fSelectedTriggers.push_back("DMC8DG1");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("EDMC8G1") == combinedtriggers.end())) combinedtriggers.insert("EDMC8G1");
        if(!(isMinBiasT0 || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgDG2])) fSelectedTriggers.push_back("DMC8DG1excl");
      }
    }
  }

  if(fEnableNoINTTriggers) {
    if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0]) {
      // EMC8 trigger
      AliDebugStream(1) << "Event selected as 0EMC" << std::endl;
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMC") == combinedtriggers.end())) combinedtriggers.insert("0EDMC");
      fSelectedTriggers.push_back("0EMC");
    }
    if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2]){
      AliDebugStream(1) << "Event selected as EJ2 (E0MC)" << std::endl;
      fSelectedTriggers.push_back("0EMCEJ2");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMCJ2") == combinedtriggers.end())) combinedtriggers.insert("0EDMCJ2");
      if(!emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0]) fSelectedTriggers.push_back("0EMCEJ2excl");
    }
    if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ1]){
      AliDebugStream(1) << "Event selected as EJ1 (0EMC)" << std::endl;
      fSelectedTriggers.push_back("0EMCEJ1");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMCJ1") == combinedtriggers.end())) combinedtriggers.insert("0EDMCJ1");
      if(!(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEJ2])) fSelectedTriggers.push_back("0EMCEJ1excl");
    }
    if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2]){
      AliDebugStream(1) << "Event selected as EG2 (0EMC)" << std::endl;
      fSelectedTriggers.push_back("0EMCEG2");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMCG2") == combinedtriggers.end())) combinedtriggers.insert("0EDMCG2");
      if(!emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0]) fSelectedTriggers.push_back("0EMCEG2excl");
    }
    if(emc8Triggers[AliEmcalTriggerOfflineSelection::kTrgEG1]){
      AliDebugStream(1) << "Event selected as EG1 (0EMC)" << std::endl;
      fSelectedTriggers.push_back("0EMCEG1");
      if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMCG1") == combinedtriggers.end())) combinedtriggers.insert("0EDMCG1");
      if(!(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEL0] || emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgEG2])) fSelectedTriggers.push_back("0EMCEG1excl");
    }

    if(fEnableDCALTriggers){
      // Handle DCAL triggers only in case DCAL triggers are enabled,
      // otherwise ignore results of the online/offline trigger selection
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0]){
        AliDebugStream(1) << "Event selected as 0DMC" << std::endl;
        fSelectedTriggers.push_back("0DMC");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMC") == combinedtriggers.end())) combinedtriggers.insert("0EDMC");
      }
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2]){
        AliDebugStream(1) << "Event selected as DJ2 (0DMC)" << std::endl;
        fSelectedTriggers.push_back("0DMCDJ2");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMCJ2") == combinedtriggers.end())) combinedtriggers.insert("0EDMCJ2");
        if(!emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0]) fSelectedTriggers.push_back("0DMCDJ2excl");
      }
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ1]){
        AliDebugStream(1) << "Event selected as DJ1 (0DMC)" << std::endl;
        fSelectedTriggers.push_back("0DMCDJ1");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMCJ1") == combinedtriggers.end())) combinedtriggers.insert("0EDMCJ1");
        if(!(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDJ2])) fSelectedTriggers.push_back("0DMCDJ1excl");
      }
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDG2]){
        AliDebugStream(1) << "Event selected as DG2 (0DMC)" << std::endl;
        fSelectedTriggers.push_back("0DMCDG2");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMCG2") == combinedtriggers.end())) combinedtriggers.insert("0EDMCG2");
        if(!emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0]) fSelectedTriggers.push_back("0DMCDG2excl");
      }
      if(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDG1]){
        AliDebugStream(1) << "Event selected as DG1 (0DMC)" << std::endl;
        fSelectedTriggers.push_back("0DMCDG1");
        if(fEnableEDCombinedTriggers && (combinedtriggers.find("0EDMCG1") == combinedtriggers.end())) combinedtriggers.insert("0EDMCG1");
        if(!(emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDL0] || emcNoIntTriggers[AliEmcalTriggerOfflineSelection::kTrgDG2])) fSelectedTriggers.push_back("0DMCDG1excl");
      }
    }
  }

  for(const auto &ct : combinedtriggers) fSelectedTriggers.push_back(ct.data());
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
  const std::array<TString, 5>  edccombinedV0 = {{"EDMC7", "EDG1", "EDG2", "EDJ1", "EDJ2"}},
                                edccombinedT0 = {{"EDMC8", "EDMC8G1", "EDMC8G2", "EDMC8J1", "EDMC8J2"}},
                                edccombinedNoINT = {{"0EDMC", "0EDG1", "0EDG2", "0EDJ1", "EDJ2"}};
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
    if(fEnableEDCombinedTriggers) {
      for(const auto &t : edccombinedV0) triggers.push_back(t);
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
    if(fEnableEDCombinedTriggers) {
      for(const auto &t : edccombinedT0) triggers.push_back(t);
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
    if(fEnableEDCombinedTriggers) {
      for(const auto &t : edccombinedNoINT) triggers.push_back(t);
    }
  }
  return triggers;
}

bool AliAnalysisTaskEmcalTriggerBase::MatchTriggerFromPattern(EMCAL_STRINGVIEW pattern, EMCAL_STRINGVIEW trigger) const {
  std::string patternstring = pattern.data();
  bool isEqual = pattern[0] == '=';
  if(pattern[0] == '=') patternstring = pattern.substr(1).data();
  std::vector<std::string> classes;
  if(patternstring.find("|") != std::string::npos){
    std::stringstream decoder(patternstring);
    std::string tmp;
    while(std::getline(decoder, tmp, '|')) classes.emplace_back(tmp);
  } else classes.emplace_back(patternstring);
  bool found(false);
  for(const auto &t : classes){
    if(isEqual && (trigger == t)) {
      found = true;
      break;
    }
    if(!isEqual && (trigger.find(t) != std::string::npos)){
      found = true;
      break;
    }
  }
  return found;
}

bool AliAnalysisTaskEmcalTriggerBase::MatchTriggerFromContainer(EMCAL_STRINGVIEW pattern, const PWG::EMCAL::AliEmcalTriggerDecisionContainer *trgsel) const {
  if(!trgsel) return false;
  std::string patternstring = pattern.data();
  if(pattern[0] == '=') patternstring = pattern.substr(1).data();
  std::vector<std::string> classes;
  if(patternstring.find("|") != std::string::npos){
    std::stringstream decoder(patternstring);
    std::string tmp;
    while(std::getline(decoder, tmp, '|')) classes.emplace_back(tmp);
  } else classes.emplace_back(patternstring);
  bool found(false);
  for(const auto &t : classes) {
    if(trgsel->IsEventSelected(t.data())) {
      found = true;
      break;
    }
  }
  return found;
}

Double_t AliAnalysisTaskEmcalTriggerBase::GetTriggerWeight(EMCAL_STRINGVIEW triggerclass) const {
  if(fDownscaleFactors){
    TParameter<double> *result(nullptr);
    // Downscaling only done on MB, L0 and the low threshold triggers
    if(triggerclass.find("MB") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("INT7"));
    else if(triggerclass.find("EMC7") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EMC7"));
    else if(triggerclass.find("DMC7") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DMC7"));
    else if(triggerclass.find("EJ2") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ2"));
    else if(triggerclass.find("EJ1") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ1"));
    else if(triggerclass.find("DJ2") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DJ2"));
    else if(triggerclass.find("DJ1") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DJ1"));
    else if(triggerclass.find("EG2") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG2"));
    else if(triggerclass.find("EG1") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG1"));
    else if(triggerclass.find("DG2") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DG2"));
    else if(triggerclass.find("DG1") != std::string::npos) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DG1"));
    // for the combined classes check the downscaling of EMCAL and DCAL, usually it was the same
    // in case they are the same return downscaling of the EMCAL, otherwise return the larger of 
    // the two corresponding to the larger luminosity of the trigger
    if(triggerclass.find("ED") != std::string::npos) {
      TParameter<double> *resultEMCAL(nullptr), *resultDCAL(nullptr);
      if(triggerclass.find("EDMC7") != std::string::npos) {
        resultEMCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EMC7"));
        resultDCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DMC7"));
      } 
      else if(triggerclass.find("EDJ2") != std::string::npos) {
        resultEMCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ2"));
        resultDCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DJ2"));
      }
      else if(triggerclass.find("EDJ1") != std::string::npos) {
        resultEMCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ1"));
        resultDCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DJ1"));
      }
      else if(triggerclass.find("EDG2") != std::string::npos) {
        resultEMCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG2"));
        resultDCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DG2"));
      }
      else if(triggerclass.find("EDG1") != std::string::npos) {
        resultEMCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG1"));
        resultDCAL = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("DG1"));
      }
      if(resultEMCAL && resultDCAL) {
        auto valEMCAL = resultEMCAL->GetVal(), valDCAL = resultDCAL->GetVal();
        if(TMath::Abs(valEMCAL - valDCAL) < DBL_EPSILON) result = resultEMCAL;
        else if(valEMCAL > valDCAL) result = resultEMCAL;
        else result = resultDCAL;
      } else if(resultEMCAL) result = resultEMCAL;
      else if(resultDCAL) result = resultDCAL;
    }
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
      if(test.Contains("WU") || test.Contains("H")) return false; // Run1: Reject TRD triggers
      if(!(test.Contains(t + "-B-") || test.Contains(t + "-S-"))) return false;
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
