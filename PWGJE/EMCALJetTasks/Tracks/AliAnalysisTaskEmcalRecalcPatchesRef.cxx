/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <algorithm>
#include <array>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <string>
#include <TClonesArray.h>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TList.h>
#include <TObjString.h>
#include <TString.h>
#include "AliAnalysisManager.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerStringDecoder.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTaskEmcalRecalcPatchesRef.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalRecalcPatchesRef);

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalRecalcPatchesRef::AliAnalysisTaskEmcalRecalcPatchesRef():
  AliAnalysisTaskEmcalTriggerBase(),
  fEnableSumw2(false),
  fOnlineThresholds(),
  fSwapPatches(false),
  fRequiredOverlaps(),
  fExcludedOverlaps(),
  fCentralityRange(-999., 999.),
  fUseRecalcPatches(false),
  fRequestCentrality(false),
  fFillTHnSparse(true),
  fEventCentrality(0)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalRecalcPatchesRef::AliAnalysisTaskEmcalRecalcPatchesRef(const char * name):
  AliAnalysisTaskEmcalTriggerBase(name),
  fEnableSumw2(false),
  fOnlineThresholds(),
  fSwapPatches(false),
  fRequiredOverlaps(),
  fExcludedOverlaps(),
  fCentralityRange(-999., 999.),
  fUseRecalcPatches(false),
  fRequestCentrality(false),
  fFillTHnSparse(true),
  fEventCentrality(0)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  fOnlineThresholds.Set(8);
}

void AliAnalysisTaskEmcalRecalcPatchesRef::CreateUserHistos(){
  const std::array<std::string, 9> kNamesTriggerClasses = {{"MB", "EG1", "EG2", "DG1", "DG2", "EJ1", "EJ2", "DJ1", "DJ2"}};
  const std::array<std::string, 6> kNamesTriggerClusters = {{"ANY", "CENT", "CENTNOTRD", "BOTH", "ONLYCENT", "ONLYCENTNOTRD"}};
  const std::array<std::string, 4> kNamesPatchTypes = {{"EGA", "DGA", "EJE", "DJE"}};

  TLinearBinning ADCBinning(2000, 0., 2000.), colbinning(48, -0.5, 47.5), rowbinning(104, -0.5, 103.5), npatchbinning(51, -0.5, 50.5), noverlapbinning(21, -0.5, 20.5);
  const TBinning *firedpatchbinning[5] = {&ADCBinning, &colbinning, &rowbinning, &npatchbinning, &noverlapbinning},
                 *allpatchbinning[3] = {&ADCBinning, &npatchbinning, &noverlapbinning};
  // Create event counters
  for(const auto &kt : kNamesTriggerClasses){
    // Min. Bias: Only CENT cluster - no distinction between trigger clusters
    // EMCAL triggers: Distinction between CENT and CENTNOTRD necessary
    if(kt.find("MB") != std::string::npos){
      fHistos->CreateTH1(Form("hEventCounter%s", kt.data()), Form("Event counter for %s", kt.data()), 1, 0.5, 1.5);
      fHistos->CreateTH1(Form("hEventCounterWeighted%s", kt.data()), Form("Event counter for %s", kt.data()), 1, 0.5, 1.5);
    } else {
      for(const auto &kc : kNamesTriggerClusters) {
        fHistos->CreateTH1(Form("hEventCounter%s%s", kt.data(), kc.data()), Form("Event counter for %s in cluster %s", kt.data(), kc.data()), 1, 0.5, 1.5);
        fHistos->CreateTH1(Form("hEventCounterWeighted%s%s", kt.data(), kc.data()), Form("Event counter for %s in cluster %s", kt.data(), kc.data()), 1, 0.5, 1.5);
     }
    }
  } 

  // Min. Bias: Create patch energy spectra (all patches) for EMCAL and DCAL
  // Min. Bias trigger only in CENT cluster
  for(const auto &kp : kNamesPatchTypes) {
    fHistos->CreateTH1(Form("hPatchADC%sMB", kp.data()), Form("Patch ADC spectra for %s patches in MB events", kp.data()), 2000, 0., 2000., fEnableSumw2 ? "s" : "");
    fHistos->CreateTH1(Form("hPatchADCWeighted%sMB", kp.data()), Form("Patch ADC spectra for %s patches in MB events", kp.data()), 2000, 0., 2000., fEnableSumw2 ? "s" : "");
  }

  // Triggers: Create trigger spectra and THnSparse of firing patches
  for(const auto &kt : kNamesTriggerClasses) {
    if(kt == "MB") continue;
    const char detector = kt[0];
    const char *patchtype = ((kt[1] == 'G') ? "GA" : "JE");
    // distinction between trigger clusters
    for(const auto &kc : kNamesTriggerClusters){
      fHistos->CreateTH1(Form("hPatchADC%c%s%s%s", detector, patchtype, kt.data(), kc.data()), Form("Patch ADC spectra for %c%s patches in %s events (cluster %s)", detector, patchtype, kt.data(), kc.data()), 2000, 0., 2000., fEnableSumw2 ? "s" : "");
      fHistos->CreateTH1(Form("hPatchADCWeighted%c%s%s%s", detector, patchtype, kt.data(), kc.data()), Form("Patch ADC spectra for %c%s patches in %s events (cluster %s)", detector, patchtype, kt.data(), kc.data()), 2000, 0., 2000., fEnableSumw2 ? "s" : "");
      if(fFillTHnSparse) {
        fHistos->CreateTHnSparse(Form("hFiredPatches%c%s%s%s", detector, patchtype, kt.data(), kc.data()), Form("Fired %c%s patches for trigger %s (cluster %s)", detector, patchtype, kt.data(), kc.data()), 5, firedpatchbinning, fEnableSumw2 ? "s" : "");
        fHistos->CreateTHnSparse(Form("hAllPatches%c%s%s%s", detector, patchtype, kt.data(), kc.data()), Form("Fired %c%s patches for trigger %s (cluster %s)", detector, patchtype, kt.data(), kc.data()), 3, allpatchbinning, fEnableSumw2 ? "s" : "");
        fHistos->CreateTHnSparse(Form("hFiredPatchesWeighted%c%s%s%s", detector, patchtype, kt.data(), kc.data()), Form("Fired %c%s patches for trigger %s (cluster %s)", detector, patchtype, kt.data(), kc.data()), 5, firedpatchbinning, fEnableSumw2 ? "s" : "");
        fHistos->CreateTHnSparse(Form("hAllPatchesWeighted%c%s%s%s", detector, patchtype, kt.data(), kc.data()), Form("Fired %c%s patches for trigger %s (cluster %s)", detector, patchtype, kt.data(), kc.data()), 3, allpatchbinning, fEnableSumw2 ? "s" : "");
      }
    } 
  }
}

bool AliAnalysisTaskEmcalRecalcPatchesRef::IsUserEventSelected(){
  
  fEventCentrality = -1;
  if(fRequestCentrality){
    AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if(!mult){
      AliErrorStream() << GetName() << ": Centrality selection enabled but no centrality estimator found" << std::endl;
      return false;
    }
    if(mult->IsEventSelected()) return false;
    fEventCentrality = mult->GetEstimator("V0M")->GetPercentile();
    AliDebugStream(1) << GetName() << ": Centrality " <<  fEventCentrality << std::endl;
    if(!fCentralityRange.IsInRange(fEventCentrality)){
      AliDebugStream(1) << GetName() << ": reject centrality: " << fEventCentrality << std::endl;
      return false;
    } else {
      AliDebugStream(1) << GetName() << ": select centrality " << fEventCentrality << std::endl;
    }
  } else {
    AliDebugStream(1) << GetName() << ": No centrality selection applied" << std::endl;
  }

  // handle overlaps
  if(fRequiredOverlaps.GetEntries() || fExcludedOverlaps.GetEntries()){
    if(fRequiredOverlaps.GetEntries()){
      Bool_t allFound(true);
      for(auto t : fRequiredOverlaps){
        auto trgstr = static_cast<TObjString *>(t)->String();
        if(std::find_if(fSelectedTriggers.begin(), fSelectedTriggers.end(), [&trgstr](const TString &seltrigger) -> bool { return seltrigger.Contains(trgstr); }) == fSelectedTriggers.end()) {
          allFound = false;
          break;
        }
      }
      if(!allFound) return false;
    }

    if(fExcludedOverlaps.GetEntries()){
      Bool_t oneFound(false);
      for(auto t : fExcludedOverlaps){
        auto trgstr = static_cast<TObjString *>(t)->String();
        if(std::find_if(fSelectedTriggers.begin(), fSelectedTriggers.end(), [&trgstr](const TString &seltrigger) -> bool { return seltrigger.Contains(trgstr); }) != fSelectedTriggers.end()) {
          oneFound = true;
          break;
        }
      }
      if(oneFound) return false;
    }
  }
  return true;
}

void AliAnalysisTaskEmcalRecalcPatchesRef::UserFillHistosAfterEventSelection(){
  const std::array<std::string, 9> kNamesTriggerClasses = {{"MB", "EG1", "EG2", "DG1", "DG2", "EJ1", "EJ2", "DJ1", "DJ2"}};
  const auto selclusters = GetAcceptedTriggerClusters(fInputEvent->GetFiredTriggerClasses().Data());

  auto findTriggerType = [](const std::vector<TString> &triggers, TString type) -> bool {
    bool found = false;
    for(const auto &t : triggers) {
      if(t.Contains(type)) {
        found = true;
        break;
      }
    }
    return found;
  };

  std::vector<TString> handledtriggers;
  for(auto t  : kNamesTriggerClasses) {
    if(findTriggerType(fSelectedTriggers, t.data())){
      handledtriggers.emplace_back(t);
    } 
  }

  for(const auto &kt : handledtriggers) {
    if(kt == "MB") {
      fHistos->FillTH1(Form("hEventCounter%s", kt.Data()), 1.);
      fHistos->FillTH1(Form("hEventCounterWeighted%s", kt.Data()), GetTriggerWeight("MB"));
    } else {
      auto weight = GetTriggerWeight(kt.Data());
      for(const auto &kc : selclusters) {
        fHistos->FillTH1(Form("hEventCounter%s%s", kt.Data(), kc.data()), 1.);
        fHistos->FillTH1(Form("hEventCounterWeighted%s%s", kt.Data(), kc.data()), weight);
      }
    }
  }
}

void AliAnalysisTaskEmcalRecalcPatchesRef::RunChanged(Int_t newrun){
  LoadTriggerThresholdsFromCDB();
}

bool AliAnalysisTaskEmcalRecalcPatchesRef::Run(){
  const std::array<std::string, 9> kNamesTriggerClasses = {{"MB", "EG1", "EG2", "DG1", "DG2", "EJ1", "EJ2", "DJ1", "DJ2"}};
  const std::unordered_map<std::string, ETriggerThreshold_t> kPatchIndex = {{"EG1", kThresholdEG1},{"EG2", kThresholdEG2},{"DG1", kThresholdDG1},{"DG2", kThresholdDG2},
                                                                            {"EJ1", kThresholdEJ1},{"EJ2", kThresholdEJ2},{"DJ1", kThresholdDJ1},{"DJ2", kThresholdDJ2}};
  if(!fSelectedTriggers.size()) return false;       // no trigger selected
  AliDebugStream(1) << "Found triggers" << std::endl;
  for(auto t : fSelectedTriggers) AliDebugStream(1) << t << std::endl;
  AliDebugStream(1) << "Trigger patch container has " << fTriggerPatchInfo->GetEntries() << " patches" << std::endl;

  // Decode trigger clusters
  const auto selclusters = GetAcceptedTriggerClusters(fInputEvent->GetFiredTriggerClasses().Data());

  auto findTriggerType = [](const std::vector<TString> &triggers, TString type) -> bool {
    bool found = false;
    for(const auto &t : triggers) {
      if(t.Contains(type)) {
        found = true;
        break;
      }
    }
    return found;
  };
  bool isMB = std::find(fSelectedTriggers.begin(), fSelectedTriggers.end(), "MB") != fSelectedTriggers.end(),
       isEG = findTriggerType(fSelectedTriggers, "EG"),
       isDG = findTriggerType(fSelectedTriggers, "DG"),
       isEJ = findTriggerType(fSelectedTriggers, "EJ"),
       isDJ = findTriggerType(fSelectedTriggers, "DJ");

  std::vector<TString> handledtriggers;
  for(auto t  : kNamesTriggerClasses) {
    if(findTriggerType(fSelectedTriggers, t.data())){
      handledtriggers.emplace_back(t);
    } 
  }

  AliDebugStream(1) << "Processing triggers" << std::endl;
  for(auto t : handledtriggers) AliDebugStream(1) << t << std::endl;
  if(!handledtriggers.size()){
    AliDebugStream(1) << "No supported trigger class found " << std::endl;
    return false;
  }

  std::vector<const AliEMCALTriggerPatchInfo *> EGApatches, DGApatches, EJEpatches, DJEpatches;
  if(isMB || isEG) EGApatches = SelectAllPatchesByType(*fTriggerPatchInfo, fSwapPatches ? kEJEpatches : kEGApatches);
  if(isMB || isDG) DGApatches = SelectAllPatchesByType(*fTriggerPatchInfo, fSwapPatches ? kDJEpatches : kDGApatches);
  if(isMB || isEJ) EJEpatches = SelectAllPatchesByType(*fTriggerPatchInfo, fSwapPatches ? kEGApatches : kEJEpatches);
  if(isMB || isDJ) DJEpatches = SelectAllPatchesByType(*fTriggerPatchInfo, fSwapPatches ? kDGApatches : kDJEpatches);
  
  for(const auto &t : handledtriggers) {
    if(t == "MB") {
      // Min bias: Only fill patch ADC spectra all patches
      auto mbweight = GetTriggerWeight("MB");
      for(auto patch :  EGApatches) {
        fHistos->FillTH1("hPatchADCEGAMB", patch->GetADCAmp()); 
        fHistos->FillTH1("hPatchADCWeightedEGAMB", patch->GetADCAmp(), mbweight); 
      }
      for(auto patch :  DGApatches) {
        fHistos->FillTH1("hPatchADCDGAMB", patch->GetADCAmp()); 
        fHistos->FillTH1("hPatchADCWeightedDGAMB", patch->GetADCAmp(), mbweight); 
      }
      for(auto patch :  EJEpatches) {
        fHistos->FillTH1("hPatchADCEJEMB", patch->GetADCAmp());
        fHistos->FillTH1("hPatchADCWeightedEJEMB", patch->GetADCAmp(), mbweight);
      }
      for(auto patch :  DJEpatches) {
        fHistos->FillTH1("hPatchADCDJEMB", patch->GetADCAmp()); 
        fHistos->FillTH1("hPatchADCWeightedDJEMB", patch->GetADCAmp(), mbweight); 
      }
    } else {
      const char detector = t[0];
      const char *patchtype = ((t[1] == 'G') ? "GA" : "JE");
      auto triggerweight = GetTriggerWeight(t.Data());
      std::vector<const AliEMCALTriggerPatchInfo *> &patchhandler = (detector == 'E' ? (t[1] == 'G' ? EGApatches : EJEpatches) : (t[1] == 'G' ? DGApatches : DJEpatches)); 
      auto firedpatches = SelectFiredPatchesByTrigger(*fTriggerPatchInfo, kPatchIndex.find(t.Data())->second);
      auto patchareas = GetNumberNonOverlappingPatchAreas(firedpatches);
      AliDebugStream(3) << "Trigger " << t << ", patches " << patchhandler.size() << ", firing " << firedpatches.size() << std::endl;
      for(auto p : patchhandler){
        double point[3] = {static_cast<double>(p->GetADCAmp()), static_cast<double>(firedpatches.size()), static_cast<double>(patchareas)};
        for(const auto &kc : selclusters) {
          fHistos->FillTH1(Form("hPatchADC%c%s%s%s", detector, patchtype, t.Data(), kc.data()), p->GetADCAmp());
          fHistos->FillTH1(Form("hPatchADCWeighted%c%s%s%s", detector, patchtype, t.Data(), kc.data()), p->GetADCAmp(), triggerweight);
          if(fFillTHnSparse){
            fHistos->FillTHnSparse(Form("hAllPatches%c%s%s%s", detector, patchtype, t.Data(), kc.data()), point);
            fHistos->FillTHnSparse(Form("hAllPatchesWeighted%c%s%s%s", detector, patchtype, t.Data(), kc.data()), point, triggerweight);
          }
        }
      }
      if(fFillTHnSparse) {
        for(auto p : firedpatches) {
          double point[5] = {static_cast<double>(p->GetADCAmp()), static_cast<double>(p->GetColStart()), static_cast<double>(p->GetRowStart()), static_cast<double>(firedpatches.size()), static_cast<double>(patchareas)};
          for(const auto &kc : selclusters) {
            fHistos->FillTHnSparse(Form("hFiredPatches%c%s%s%s", detector, patchtype, t.Data(), kc.data()), point);
            fHistos->FillTHnSparse(Form("hFiredPatchesWeighted%c%s%s%s", detector, patchtype, t.Data(), kc.data()), point, triggerweight);
          }
        }
      }
    }
  }
  return true;
}

void AliAnalysisTaskEmcalRecalcPatchesRef::LoadTriggerThresholdsFromCDB(){
  AliCDBManager *mgr = AliCDBManager::Instance();   // Expect to be configured elsewhere
  AliCDBEntry * en = mgr->Get("EMCAL/Calib/Trigger");
  AliEMCALTriggerDCSConfig *trgconf = static_cast<AliEMCALTriggerDCSConfig *>(en->GetObject());

  auto emcalstu = trgconf->GetSTUDCSConfig(false),
       dcalstu = trgconf->GetSTUDCSConfig(true);

  if(emcalstu) {
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdEG1, emcalstu->GetG(2,0));
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdEG2, emcalstu->GetG(2,1));
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdEJ1, emcalstu->GetJ(2,0));
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdEJ2, emcalstu->GetJ(2,1));
  } else {
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdEG1, 0);
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdEG2, 0);
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdEJ1, 0);
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdEJ2, 0);
  }
  if(dcalstu) {
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdDG1, dcalstu->GetG(2,0));
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdDG2, dcalstu->GetG(2,1));
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdDJ1, dcalstu->GetJ(2,0));
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdDJ2, dcalstu->GetJ(2,1));  
  } else {
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdDG1, 0);
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdDG2, 0);
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdDJ1, 0);
    SetOnlineThreshold(AliAnalysisTaskEmcalRecalcPatchesRef::kThresholdDJ2, 0);  
  }
}

std::vector<const AliEMCALTriggerPatchInfo *> AliAnalysisTaskEmcalRecalcPatchesRef::SelectAllPatchesByType(const TClonesArray &list, EPatchType_t patchtype) const {
  AliDebugStream(2) << "Selecting all patches for trigger " << static_cast<int>(patchtype) << std::endl;
  std::vector<const AliEMCALTriggerPatchInfo *> result;
  for(auto p : list){
    AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(p);
    if(!patch->IsRecalc()) continue;
    bool selected(true);
    switch(patchtype){
    case kEGApatches: if(patch->IsDCalPHOS() || !patch->IsGammaHighRecalc()) selected = false; break;
    case kDGApatches: if(patch->IsEMCal() || !patch->IsGammaHighRecalc()) selected = false; break;
    case kEJEpatches: if(patch->IsDCalPHOS() || !patch->IsJetHighRecalc()) selected = false; break;
    case kDJEpatches: if(patch->IsEMCal() || !patch->IsJetHighRecalc()) selected = false; break;
    };
    if(selected) result.emplace_back(patch);
  }
  AliDebugStream(2) << "In: " << list.GetEntries() << ", out: " << result.size() << std::endl;
  return result;
}

std::vector<const AliEMCALTriggerPatchInfo *> AliAnalysisTaskEmcalRecalcPatchesRef::SelectFiredPatchesByTrigger(const TClonesArray &list, ETriggerThreshold_t trigger) const {
  std::vector<const AliEMCALTriggerPatchInfo *> result;
  EPatchType_t patchtype;
  switch(trigger) {
  case kThresholdEG1: patchtype = kEGApatches; break;
  case kThresholdEG2: patchtype = kEGApatches; break;
  case kThresholdDG1: patchtype = kDGApatches; break;
  case kThresholdDG2: patchtype = kDGApatches; break;
  case kThresholdEJ1: patchtype = kEJEpatches; break;
  case kThresholdEJ2: patchtype = kEJEpatches; break;
  case kThresholdDJ1: patchtype = kDJEpatches; break;
  case kThresholdDJ2: patchtype = kDJEpatches; break;
  default: return result;     // unsupported patch type - return empty list
  };
  for(auto patch : SelectAllPatchesByType(list, patchtype)) {
    if(patch->GetADCAmp() > fOnlineThresholds[trigger]) result.emplace_back(patch);
  }
  return result;
}

std::vector<std::string> AliAnalysisTaskEmcalRecalcPatchesRef::GetAcceptedTriggerClusters(const char *triggerstring) const {
  auto clusters = PWG::EMCAL::Triggerinfo::DecodeTriggerString(triggerstring);
  std::vector<std::string> selclusters;
  selclusters.push_back("ANY");
  bool isCENT(false), isCENTNOTRD(false);
  for(const auto &c : clusters){
    if((c.Triggercluster() == "CENT") && !isCENT) { // don't count clusters multiple times
      selclusters.push_back("CENT");
      isCENT = true;
    } else if((c.Triggercluster() == "CENTNOTRD") && !isCENTNOTRD) { // don't count clusters multiple times
      selclusters.push_back("CENTNOTRD");
      isCENTNOTRD = true;
    }
  }
  if(isCENT && isCENTNOTRD) selclusters.push_back("BOTH");
  else {
    if(isCENT) selclusters.push_back("ONLYCENT");
    if(isCENTNOTRD) selclusters.push_back("ONLYCENTNOTRD");
  }
  return selclusters;
}

int AliAnalysisTaskEmcalRecalcPatchesRef::GetNumberNonOverlappingPatchAreas(const std::vector<const AliEMCALTriggerPatchInfo *> &firedpatches) const {
  std::vector<const AliEMCALTriggerPatchInfo *> patchareas;
  for(const auto patch : firedpatches) {
    if(!patchareas.size()) patchareas.push_back(patch); // first patch always considered as new area
    else {
      bool overlapFound = false;
      for(const auto refpatch : patchareas) {
        if(!refpatch) {
          AliErrorStream() << "Ref patch null" << std::endl;
          AliErrorStream() << "Patchareas has size " << patchareas.size() << std::endl;
          AliErrorStream() << "Firedpatches has size " << firedpatches.size() << std::endl;
        } 
        if(!patch){
          AliErrorStream() << "Test patch null" << std::endl;
          AliErrorStream() << "Patchareas has size " << patchareas.size() << std::endl;
          AliErrorStream() << "Firedpatches has size " << firedpatches.size() << std::endl;
        }
        if(HasOverlap(*refpatch, *patch)) {
          // patch has overlap with allready accepted patch - discard
          overlapFound = true;
          break;
        }
      }
      if(!overlapFound) patchareas.emplace_back(patch); // New non-overlapping patch found
    }
  }
  return patchareas.size();
}

bool AliAnalysisTaskEmcalRecalcPatchesRef::HasOverlap(const AliEMCALTriggerPatchInfo &ref, const AliEMCALTriggerPatchInfo &test) const {
  int testcolmin = test.GetColStart(), testcolmax = test.GetColStart()+test.GetPatchSize()-1,
      testrowmin = test.GetRowStart(), testrowmax = test.GetRowStart()+test.GetPatchSize()-1,
      refcolmin = ref.GetColStart(), refcolmax = ref.GetColStart()+ref.GetPatchSize()-1,
      refrowmin = ref.GetRowStart(), refrowmax = ref.GetRowStart()+ref.GetPatchSize()-1;
  if((InRange(testcolmin, refcolmin, refcolmax) && InRange(testrowmin, refrowmin, refrowmax)) ||
     (InRange(testcolmax, refcolmin, refcolmax) && InRange(testrowmax, refrowmin, refrowmax))) return true;
  return false;
}

void AliAnalysisTaskEmcalRecalcPatchesRef::AddRequiredTriggerOverlap(const char *trigger){
  if(fRequiredOverlaps.FindObject(trigger)) return;
  fRequiredOverlaps.Add(new TObjString(trigger));
}

void AliAnalysisTaskEmcalRecalcPatchesRef::AddExcludedTriggerOverlap(const char *trigger){
  if(fExcludedOverlaps.FindObject(trigger)) return;
  fExcludedOverlaps.Add(new TObjString(trigger));
}

AliAnalysisTaskEmcalRecalcPatchesRef *AliAnalysisTaskEmcalRecalcPatchesRef::AddTaskEmcalRecalcPatches(const char *suffix) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    std::cerr << "No analysis manager defined" << std::endl;
    return nullptr;
  }

  std::stringstream taskname;
  taskname << "EmcalRecalcPatches_" << suffix;
  AliAnalysisTaskEmcalRecalcPatchesRef *task = new AliAnalysisTaskEmcalRecalcPatchesRef(taskname.str().data());
  mgr->AddTask(task);
  task->SetEnableSumw2(true);

  std::stringstream outfilename, outlistname;
  outfilename << mgr->GetCommonFileName() << ":" << taskname.str();
  outlistname << "Histos_" << taskname.str();
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(outlistname.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));

  return task;
}
