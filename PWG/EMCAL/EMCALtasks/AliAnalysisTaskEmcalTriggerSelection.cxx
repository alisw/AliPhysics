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
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <TH1.h>
#include "AliEmcalTriggerAlias.h"
#include "AliEmcalTriggerDecision.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalTriggerSelection.h"
#include "AliAnalysisTaskEmcalTriggerSelection.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliYAMLConfiguration.h"

ClassImp(PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection)

namespace PWG {
namespace EMCAL {

AliAnalysisTaskEmcalTriggerSelection::AliAnalysisTaskEmcalTriggerSelection():
  AliAnalysisTaskEmcal(),
  fTriggerDecisionContainer(nullptr),
  fGlobalDecisionContainerName("EmcalTriggerDecision"),
  fTriggerSelections(),
  fSelectionQA(),
  fUseRho(false)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetCaloTriggersName("usedefault");
  fTriggerSelections.SetOwner(kTRUE);
}

AliAnalysisTaskEmcalTriggerSelection::AliAnalysisTaskEmcalTriggerSelection(const char* name):
  AliAnalysisTaskEmcal(name, kTRUE),
  fTriggerDecisionContainer(nullptr),
  fGlobalDecisionContainerName("EmcalTriggerDecision"),
  fTriggerSelections(),
  fSelectionQA()
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetCaloTriggersName("usedefault");
  SetMakeGeneralHistograms(true);
  fTriggerSelections.SetOwner(kTRUE);
}

void AliAnalysisTaskEmcalTriggerSelection::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  for(auto s : fTriggerSelections) InitQA(static_cast<AliEmcalTriggerSelection *>(s));
  for(auto q : fSelectionQA) static_cast<AliEmcalTriggerSelectionQA *>(q)->GetHistos(fOutput);
}

void AliAnalysisTaskEmcalTriggerSelection::UserExecOnce(){
  if(!fTriggerDecisionContainer) fTriggerDecisionContainer = new AliEmcalTriggerDecisionContainer(fGlobalDecisionContainerName.Data());
  fInputEvent->AddObject(fTriggerDecisionContainer);

  // check whether any of the cuts is using rho subtraction
  fUseRho = false;
  for(auto selectionMethod : fTriggerSelections) {
    AliEmcalTriggerSelection *selection = static_cast<AliEmcalTriggerSelection *>(selectionMethod);
    if(selection->GetSelectionCuts()->IsSubtractRho()){
      fUseRho = true;
      break;
    }
  }
}

void AliAnalysisTaskEmcalTriggerSelection::AddTriggerSelection(AliEmcalTriggerSelection * const selection){
  fTriggerSelections.Add(selection);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Run(){
  AliDebugStream(1) << GetName() <<  ": Using trigger patch container " << fCaloTriggerPatchInfoName << std::endl;
  AliDebugStream(1) << GetName() <<  ": Writing trigger selection to " << fGlobalDecisionContainerName << std::endl;
  fTriggerDecisionContainer->Reset();
  AliEmcalTriggerSelection *selection(NULL);
  TIter selectionIter(&fTriggerSelections);
  // Build rho object
  AliEmcalTriggerSelectionCuts::RhoForTrigger rhocontainer {0., 0.,0.,0.,0.,0.};
  if(fUseRho) {
    // Assumpution: Rho in data stored for detector:
    rhocontainer.fRhoForEmcalOnline = fCaloTriggers->GetMedian(0);
    rhocontainer.fRhoForDCALOnline = fCaloTriggers->GetMedian(1);
    rhocontainer.fRhoForEmcalRecalc = CalculateRho(fTriggerPatchInfo, true, false);
    rhocontainer.fRhoForDCALRecalc = CalculateRho(fTriggerPatchInfo, true, true);
    rhocontainer.fRhoForEmcalRecalc = CalculateRho(fTriggerPatchInfo, false, false);
    rhocontainer.fRhoForDCALRecalc = CalculateRho(fTriggerPatchInfo, false, true);
  }
  while((selection = dynamic_cast<AliEmcalTriggerSelection *>(selectionIter()))){
    fTriggerDecisionContainer->AddTriggerDecision(selection->MakeDecison(fTriggerPatchInfo, rhocontainer));
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::FillHistograms() {
  MakeQA(GetGlobalTriggerDecisionContainer());
  return kTRUE;
}

void AliAnalysisTaskEmcalTriggerSelection::InitQA(const AliEmcalTriggerSelection *sel){
  AliEmcalTriggerSelectionQA *qa = new AliEmcalTriggerSelectionQA(sel);
  fSelectionQA.Add(qa);
}

void AliAnalysisTaskEmcalTriggerSelection::MakeQA(const AliEmcalTriggerDecisionContainer *cont) {
  for(auto d : *(cont->GetListOfTriggerDecisions())) {
    AliEmcalTriggerDecision *myd = static_cast<AliEmcalTriggerDecision *>(d);
    static_cast<AliEmcalTriggerSelectionQA *>(fSelectionQA.FindObject(myd->GetName()))->Fill(myd);
  }
}

double AliAnalysisTaskEmcalTriggerSelection::CalculateRho(const TClonesArray *const patches, Bool_t recalc, Bool_t isEmcal) {
  std::vector<double> patchenergies;
  for(auto patchobject : *patches) {
    AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(patchobject);
    if(!(patch->IsBkgSimple() || patch->IsBkgRecalc())) continue;
    if((isEmcal && !patch->IsEMCal()) || (!isEmcal && patch->IsEMCal())) continue;
    if(recalc && patch->IsBkgRecalc()) patchenergies.push_back(patch->GetADCAmp());
    else patchenergies.push_back(patch->GetPatchE());
  }
  std::sort(patchenergies.begin(), patchenergies.end(), std::less<double>());
  return TMath::Median(patchenergies.size(), patchenergies.data());
}

void AliAnalysisTaskEmcalTriggerSelection::AutoConfigure(const char *period) {
  // data
  if(Is2011PP7TeV(period)) ConfigurePP7TeV2011();
  if(Is2012PP(period)) ConfigurePP2012();
  if(Is2013PPB(period)) ConfigurePPB5TeV2013();
  if(Is2015PP5TeV(period)) ConfigurePP5TeV2015();
  if(Is2016PP(period)) ConfigurePP2016();
  if(Is2016PPB5TeV(period)) ConfigurePPB5TeV2016();
  if(Is2016PPB8TeV(period)) ConfigurePPB8TeV2016();
  if(Is2017PP5TeV(period)) ConfigurePP5TeV2017();
  // MC
  if(Is2011MCPP7TeV(period)) ConfigureMCPP7TeV2011();
  if(Is2012MCPP(period)) ConfigureMCPP2012();
  if(Is2013MCPPB(period)) ConfigureMCPPB5TeV2013();
  if(Is2015MCPP5TeV(period)) ConfigureMCPP5TeV2015();
  if(Is2016MCPPB5TeV(period)) ConfigureMCPPB5TeV2016();
  if(Is2016MCPPB8TeV(period)) ConfigureMCPPB8TeV2016();
  if(Is2016MCPP(period)) ConfigureMCPP2016();
  if(Is2017MCPP5TeV(period)) ConfigureMCPP5TeV2017();
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2011PP7TeV(const char *dataset) const {
  TString datasetstring(dataset);
  datasetstring.ToLower();
  if(datasetstring.Length() != 6) return false;     // not data period
  if(datasetstring.Contains("lhc11")){
    auto subperiod = datasetstring[5];
    if(subperiod == 'c' || subperiod == 'd') return true;
  }
  return false;
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2012PP(const char *dataset) const {
  TString datasetstring(dataset);
  datasetstring.ToLower();
  if(datasetstring.Length() != 6) return false;     // not data period
  if(datasetstring.Contains("lhc12")){
    auto subperiod = datasetstring[5];
    if(subperiod > 'b' && subperiod < 'j') return true;
  }
  return false;
}
Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2013PPB(const char *dataset) const {
  TString datasetstring(dataset);
  datasetstring.ToLower();
  if(datasetstring.Length() != 6) return false;     // not data period
  if(datasetstring.Contains("lhc13")){
    auto subperiod = datasetstring[5];
    if(subperiod > 'a' && subperiod < 'g') return true;
  }
  return false;
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2015PP5TeV(const char *dataset) const {
  TString datasetstring(dataset);
  datasetstring.ToLower();
  if(datasetstring.Length() != 6) return false;     // not data period
  if(datasetstring.Contains("lhc15n")) return true;
  return false;
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2016PP(const char *dataset) const {
  TString datasetstring(dataset);
  datasetstring.ToLower();
  if(datasetstring.Length() != 6) return false;     // not data period
  if(datasetstring.Contains("lhc16") || datasetstring.Contains("lhc17") || datasetstring.Contains("lhc18")){
    auto subperiod = datasetstring[5];
    if(datasetstring.Contains("lhc16")){
      if(subperiod > 'g' && subperiod < 'q') return true;
    }
    if(datasetstring.Contains("lhc17")) {
      if((subperiod > 'c' && subperiod < 'n') || (subperiod == 'o') || (subperiod < 'r')) return true;
    }
    if(datasetstring.Contains("lhc18")) {
      // 2018 runs will follow when taken
      return true;
    }
  }
  return false;
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2016PPB5TeV(const char *dataset) const {
  TString datasetstring(dataset);
  datasetstring.ToLower();
  if(datasetstring.Length() != 6) return false;     // not data period
  if(datasetstring.Contains("lhc16")){
    auto subperiod = datasetstring[5];
    if(subperiod == 'q' || subperiod == 't') return true;
  }
  return false;
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2016PPB8TeV(const char *dataset) const {
  TString datasetstring(dataset);
  datasetstring.ToLower();
  if(datasetstring.Length() != 6) return false;     // not data period
  if(datasetstring.Contains("lhc16")){
    auto subperiod = datasetstring[5];
    if(subperiod == 'r' || subperiod == 's') return true;
  }
  return false;
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2017PP5TeV(const char *dataset) const {
  TString datasetstring(dataset);
  datasetstring.ToLower();
  if(datasetstring.Length() != 6) return false;     // not data period
  if(datasetstring.Contains("lhc17")){
    auto subperiod = datasetstring[5];
    if(subperiod == 'p' || subperiod == 'q') return true;
  }
  return false;
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2011MCPP7TeV(const char *dataset) const {
  std::vector<TString> supportedProductions = {"lhc14k1a", "lhc14b7", "lhc14k1b"};
  return IsSupportedMCSample(dataset, supportedProductions);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2012MCPP(const char *dataset) const {
  std::vector<TString> supportedProductions = {"lhc15h1", "lhc15h2", "lhc16a1", "lhc16c2", "lhc17g5a", "lhc17g5"};
  return IsSupportedMCSample(dataset, supportedProductions);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2013MCPPB(const char *dataset) const {
  std::vector<TString> supportedProductions = {"lhc17g6a", "lhc16c3a", "lhc16c3b", "lhc18j5", "lhc19a4"};
  return IsSupportedMCSample(dataset, supportedProductions);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2015MCPP5TeV(const char *dataset) const {
  std::vector<TString> supportedProductions = {"lhc17e2", "lhc18j3"};
  return IsSupportedMCSample(dataset, supportedProductions);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2016MCPP(const char *dataset) const {
  std::vector<TString> supportedProductions = {"lhc17f8", "lhc18f5", "lhc18g2", "lhc19a1", "lhc19d3", "lhc20k1"};
  return IsSupportedMCSample(dataset, supportedProductions);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2016MCPPB5TeV(const char *dataset) const {
  std::vector<TString> supportedProductions = {"lhc18f3", "lhc17g8a"};
  return IsSupportedMCSample(dataset, supportedProductions);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2016MCPPB8TeV(const char *dataset) const {
  std::vector<TString> supportedProductions = {"lhc18f3b", "lhc18f3c", "lhc18b9b", "lhc18b9c"};
  return IsSupportedMCSample(dataset, supportedProductions);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Is2017MCPP5TeV(const char *dataset) const {
  std::vector<TString> supportedProductions = {"lhc17l3b", "lhc18j2"};
  return IsSupportedMCSample(dataset, supportedProductions);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::IsSupportedMCSample(const char *dataset, std::vector<TString> &supportedProductions) const{
  TString datasetstring(dataset);
  datasetstring.ToLower();
  bool found(false);
  for(const auto & prod : supportedProductions) {
    if(datasetstring.Contains(prod)) {
      found = true;
      break;
    }
  }
  return found;
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigurePP7TeV2011(){
  std::cout << "AutoConfigure: Configuring for data pp 7 TeV" << std::endl;
  AliEmcalTriggerSelectionCuts *emcl0cuts = new AliEmcalTriggerSelectionCuts;
  emcl0cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  emcl0cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  emcl0cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  emcl0cuts->SetUseRecalcPatches(true);
  emcl0cuts->SetThreshold(292);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EMCL0", emcl0cuts));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigureMCPP7TeV2011() {
  std::cout << "AutoConfigure: Configuring for MC pp 7 TeV" << std::endl;
  AliEmcalTriggerSelectionCuts *emcl0cuts = new AliEmcalTriggerSelectionCuts;
  emcl0cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  emcl0cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  emcl0cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  emcl0cuts->SetUseSimpleOfflinePatches(true);
  emcl0cuts->SetThreshold(5.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EMCL0", emcl0cuts));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigurePP2012(){
  std::cout << "AutoConfigure: Configuring for data pp 8 TeV" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg1cuts->SetUseRecalcPatches(true);
  eg1cuts->SetThreshold(130);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EGA", eg1cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej1cuts->SetUseRecalcPatches(true);
  ej1cuts->SetThreshold(200);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJE", ej1cuts));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigureMCPP2012() {
  std::cout << "AutoConfigure: Configuring for MC pp 8 TeV" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg1cuts->SetUseSimpleOfflinePatches(true);
  eg1cuts->SetThreshold(10.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EGA", eg1cuts, new AliEmcalTriggerAlias("EGA;EG1")));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej1cuts->SetUseSimpleOfflinePatches(true);
  ej1cuts->SetThreshold(15.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJE", ej1cuts, new AliEmcalTriggerAlias("EJE;EJ1")));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigurePPB5TeV2013(){
  std::cout << "AutoConfigure: Configuring for data p-Pb 5.02 TeV (2013)" << std::endl;
  AliEmcalTriggerSelectionCuts *emcl0cuts = new AliEmcalTriggerSelectionCuts;
  emcl0cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  emcl0cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  emcl0cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  emcl0cuts->SetUseRecalcPatches(true);
  emcl0cuts->SetThreshold(158);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EMCL0", emcl0cuts));

  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg1cuts->SetUseRecalcPatches(true);
  eg1cuts->SetThreshold(140);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg2cuts->SetUseRecalcPatches(true);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(89);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej1cuts->SetUseRecalcPatches(true);
  ej1cuts->SetThreshold(260);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej2cuts->SetUseRecalcPatches(true);
  ej2cuts->SetThreshold(127);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigureMCPPB5TeV2013() {
  std::cout << "AutoConfigure: Configuring for MC p-Pb 5.02 TeV (2013)" << std::endl;
  AliEmcalTriggerSelectionCuts *emcl0cuts = new AliEmcalTriggerSelectionCuts;
  emcl0cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  emcl0cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  emcl0cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  emcl0cuts->SetUseSimpleOfflinePatches(true);
  emcl0cuts->SetThreshold(3.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EMCL0", emcl0cuts));

  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg1cuts->SetUseSimpleOfflinePatches(true);
  eg1cuts->SetThreshold(10.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(6.7);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej1cuts->SetUseSimpleOfflinePatches(true);
  ej1cuts->SetThreshold(20.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej2cuts->SetUseSimpleOfflinePatches(true);
  ej2cuts->SetThreshold(10.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));

}

void AliAnalysisTaskEmcalTriggerSelection::ConfigurePP5TeV2015(){
  std::cout << "AutoConfigure: Configuring for data pp 5.02 TeV (2015)" << std::endl;
  AliEmcalTriggerSelectionCuts *emcl0cuts = new AliEmcalTriggerSelectionCuts;
  emcl0cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  emcl0cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  emcl0cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  emcl0cuts->SetUseRecalcPatches(true);
  emcl0cuts->SetThreshold(263);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EMCL0", emcl0cuts));

  AliEmcalTriggerSelectionCuts *dmcl0cuts = new AliEmcalTriggerSelectionCuts;
  dmcl0cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  dmcl0cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  dmcl0cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dmcl0cuts->SetUseRecalcPatches(true);
  dmcl0cuts->SetThreshold(263);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DMCL0", dmcl0cuts));

}

void AliAnalysisTaskEmcalTriggerSelection::ConfigureMCPP5TeV2015() {
  std::cout << "AutoConfigure: Configuring for MC pp 5.02 TeV (2015)" << std::endl;
  AliEmcalTriggerSelectionCuts *emcl0cuts = new AliEmcalTriggerSelectionCuts;
  emcl0cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  emcl0cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  emcl0cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  emcl0cuts->SetUseSimpleOfflinePatches(true);
  emcl0cuts->SetThreshold(5.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EMCL0", emcl0cuts));

  AliEmcalTriggerSelectionCuts *dmcl0cuts = new AliEmcalTriggerSelectionCuts;
  dmcl0cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dmcl0cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  dmcl0cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dmcl0cuts->SetUseSimpleOfflinePatches(true);
  dmcl0cuts->SetThreshold(5.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DMCL0", dmcl0cuts));

}

void AliAnalysisTaskEmcalTriggerSelection::ConfigurePP2016(){
  std::cout << "AutoConfigure: Configuring for data pp 13 TeV" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg1cuts->SetUseRecalcPatches(true);
  eg1cuts->SetThreshold(115);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg2cuts->SetUseRecalcPatches(true);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(51);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *dg1cuts = new AliEmcalTriggerSelectionCuts;
  dg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  dg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dg1cuts->SetUseRecalcPatches(true);
  dg1cuts->SetThreshold(115);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG1", dg1cuts));

  AliEmcalTriggerSelectionCuts *dg2cuts = new AliEmcalTriggerSelectionCuts;
  dg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  dg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dg2cuts->SetUseRecalcPatches(true);
  dg2cuts->SetThreshold(51);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG2", dg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej1cuts->SetUseRecalcPatches(true);
  ej1cuts->SetThreshold(255);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej2cuts->SetUseRecalcPatches(true);
  ej2cuts->SetThreshold(204);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));

  AliEmcalTriggerSelectionCuts *dj1cuts = new AliEmcalTriggerSelectionCuts;
  dj1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  dj1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dj1cuts->SetUseRecalcPatches(true);
  dj1cuts->SetThreshold(255);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ1", dj1cuts));

  AliEmcalTriggerSelectionCuts *dj2cuts = new AliEmcalTriggerSelectionCuts;
  dj2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  dj2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dj2cuts->SetUseRecalcPatches(true);
  dj2cuts->SetThreshold(204);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ2", dj2cuts));

}

void AliAnalysisTaskEmcalTriggerSelection::ConfigureMCPP2016() {
  std::cout << "AutoConfigure: Configuring for MC pp 13 TeV" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg1cuts->SetUseSimpleOfflinePatches(true);
  eg1cuts->SetThreshold(9.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(4.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *dg1cuts = new AliEmcalTriggerSelectionCuts;
  dg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  dg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dg1cuts->SetUseSimpleOfflinePatches(true);
  dg1cuts->SetThreshold(9.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG1", dg1cuts));

  AliEmcalTriggerSelectionCuts *dg2cuts = new AliEmcalTriggerSelectionCuts;
  dg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  dg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dg2cuts->SetUseSimpleOfflinePatches(true);
  dg2cuts->SetThreshold(4.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG2", dg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej1cuts->SetUseSimpleOfflinePatches(true);
  ej1cuts->SetThreshold(20.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej2cuts->SetUseSimpleOfflinePatches(true);
  ej2cuts->SetThreshold(16.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));

  AliEmcalTriggerSelectionCuts *dj1cuts = new AliEmcalTriggerSelectionCuts;
  dj1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  dj1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dj1cuts->SetUseSimpleOfflinePatches(true);
  dj1cuts->SetThreshold(20.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ1", dj1cuts));

  AliEmcalTriggerSelectionCuts *dj2cuts = new AliEmcalTriggerSelectionCuts;
  dj2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  dj2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dj2cuts->SetUseSimpleOfflinePatches(true);
  dj2cuts->SetThreshold(16.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ2", dj2cuts));

  AliEmcalTriggerSelectionCuts *emc7cuts = new AliEmcalTriggerSelectionCuts;
  emc7cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  emc7cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  emc7cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  emc7cuts->SetUseSimpleOfflinePatches(true);
  emc7cuts->SetThreshold(2.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EMCL0", emc7cuts));

  AliEmcalTriggerSelectionCuts *dmc7cuts = new AliEmcalTriggerSelectionCuts;
  dmc7cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dmc7cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  dmc7cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dmc7cuts->SetUseSimpleOfflinePatches(true);
  dmc7cuts->SetThreshold(2.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DMCL0", dmc7cuts));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigurePPB5TeV2016(){
  std::cout << "AutoConfigure: Configuring for data p-Pb 5.02 TeV (2016)" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg1cuts->SetUseRecalcPatches(true);
  eg1cuts->SetThreshold(140);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg2cuts->SetUseRecalcPatches(true);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(83);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *dg1cuts = new AliEmcalTriggerSelectionCuts;
  dg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  dg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dg1cuts->SetUseRecalcPatches(true);
  dg1cuts->SetThreshold(140);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG1", dg1cuts));

  AliEmcalTriggerSelectionCuts *dg2cuts = new AliEmcalTriggerSelectionCuts;
  dg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  dg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dg2cuts->SetUseRecalcPatches(true);
  dg2cuts->SetThreshold(83);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG2", dg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej1cuts->SetUseRecalcPatches(true);
  ej1cuts->SetThreshold(318);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej2cuts->SetUseRecalcPatches(true);
  ej2cuts->SetThreshold(255);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));

  AliEmcalTriggerSelectionCuts *dj1cuts = new AliEmcalTriggerSelectionCuts;
  dj1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  dj1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dj1cuts->SetUseRecalcPatches(true);
  dj1cuts->SetThreshold(318);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ1", dj1cuts));

  AliEmcalTriggerSelectionCuts *dj2cuts = new AliEmcalTriggerSelectionCuts;
  dj2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  dj2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dj2cuts->SetUseRecalcPatches(true);
  dj2cuts->SetThreshold(255);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ2", dj2cuts));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigureMCPPB5TeV2016() {
  std::cout << "AutoConfigure: Configuring for MC p-Pb 5.02 TeV (2016)" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg1cuts->SetUseSimpleOfflinePatches(true);
  eg1cuts->SetThreshold(11.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(6.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *dg1cuts = new AliEmcalTriggerSelectionCuts;
  dg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  dg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dg1cuts->SetUseSimpleOfflinePatches(true);
  dg1cuts->SetThreshold(11.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG1", dg1cuts));

  AliEmcalTriggerSelectionCuts *dg2cuts = new AliEmcalTriggerSelectionCuts;
  dg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  dg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dg2cuts->SetUseSimpleOfflinePatches(true);
  dg2cuts->SetThreshold(6.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG2", dg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej1cuts->SetUseSimpleOfflinePatches(true);
  ej1cuts->SetThreshold(25.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej2cuts->SetUseSimpleOfflinePatches(true);
  ej2cuts->SetThreshold(20.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));

  AliEmcalTriggerSelectionCuts *dj1cuts = new AliEmcalTriggerSelectionCuts;
  dj1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  dj1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dj1cuts->SetUseSimpleOfflinePatches(true);
  dj1cuts->SetThreshold(25.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ1", dj1cuts));

  AliEmcalTriggerSelectionCuts *dj2cuts = new AliEmcalTriggerSelectionCuts;
  dj2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  dj2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dj2cuts->SetUseSimpleOfflinePatches(true);
  dj2cuts->SetThreshold(20.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ2", dj2cuts));

  AliEmcalTriggerSelectionCuts *emc7cuts = new AliEmcalTriggerSelectionCuts;
  emc7cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  emc7cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  emc7cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  emc7cuts->SetUseSimpleOfflinePatches(true);
  emc7cuts->SetThreshold(2.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EMCL0", emc7cuts));

  AliEmcalTriggerSelectionCuts *dmc7cuts = new AliEmcalTriggerSelectionCuts;
  dmc7cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dmc7cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  dmc7cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dmc7cuts->SetUseSimpleOfflinePatches(true);
  dmc7cuts->SetThreshold(2.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DMCL0", dmc7cuts));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigurePPB8TeV2016(){
  std::cout << "AutoConfigure: Configuring for data p-Pb 8.16 TeV (2016)" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg1cuts->SetUseRecalcPatches(true);
  eg1cuts->SetThreshold(102);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg2cuts->SetUseRecalcPatches(true);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(70);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *dg1cuts = new AliEmcalTriggerSelectionCuts;
  dg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  dg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dg1cuts->SetUseRecalcPatches(true);
  dg1cuts->SetThreshold(102);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG1", dg1cuts));

  AliEmcalTriggerSelectionCuts *dg2cuts = new AliEmcalTriggerSelectionCuts;
  dg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  dg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dg2cuts->SetUseRecalcPatches(true);
  dg2cuts->SetThreshold(70);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG2", dg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej1cuts->SetUseRecalcPatches(true);
  ej1cuts->SetThreshold(293);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej2cuts->SetUseRecalcPatches(true);
  ej2cuts->SetThreshold(229);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));

  AliEmcalTriggerSelectionCuts *dj1cuts = new AliEmcalTriggerSelectionCuts;
  dj1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  dj1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dj1cuts->SetUseRecalcPatches(true);
  dj1cuts->SetThreshold(293);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ1", dj1cuts));

  AliEmcalTriggerSelectionCuts *dj2cuts = new AliEmcalTriggerSelectionCuts;
  dj2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  dj2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dj2cuts->SetUseRecalcPatches(true);
  dj2cuts->SetThreshold(229);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ2", dj2cuts));

}

void AliAnalysisTaskEmcalTriggerSelection::ConfigureMCPPB8TeV2016() {
  std::cout << "AutoConfigure: Configuring for MC p-Pb 8.16 TeV (2016)" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg1cuts->SetUseSimpleOfflinePatches(true);
  eg1cuts->SetThreshold(7.8);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(5.2);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *dg1cuts = new AliEmcalTriggerSelectionCuts;
  dg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  dg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dg1cuts->SetUseSimpleOfflinePatches(true);
  dg1cuts->SetThreshold(7.8);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG1", dg1cuts));

  AliEmcalTriggerSelectionCuts *dg2cuts = new AliEmcalTriggerSelectionCuts;
  dg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  dg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dg2cuts->SetUseSimpleOfflinePatches(true);
  dg2cuts->SetThreshold(5.2);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG2", dg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej1cuts->SetUseSimpleOfflinePatches(true);
  ej1cuts->SetThreshold(23.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej2cuts->SetUseSimpleOfflinePatches(true);
  ej2cuts->SetThreshold(18.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));

  AliEmcalTriggerSelectionCuts *dj1cuts = new AliEmcalTriggerSelectionCuts;
  dj1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  dj1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dj1cuts->SetUseSimpleOfflinePatches(true);
  dj1cuts->SetThreshold(23.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ1", dj1cuts));

  AliEmcalTriggerSelectionCuts *dj2cuts = new AliEmcalTriggerSelectionCuts;
  dj2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  dj2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dj2cuts->SetUseSimpleOfflinePatches(true);
  dj2cuts->SetThreshold(18.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ2", dj2cuts));

  AliEmcalTriggerSelectionCuts *emc7cuts = new AliEmcalTriggerSelectionCuts;
  emc7cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  emc7cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  emc7cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  emc7cuts->SetUseSimpleOfflinePatches(true);
  emc7cuts->SetThreshold(2.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EMCL0", emc7cuts));

  AliEmcalTriggerSelectionCuts *dmc7cuts = new AliEmcalTriggerSelectionCuts;
  dmc7cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dmc7cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL0Patch);
  dmc7cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dmc7cuts->SetUseSimpleOfflinePatches(true);
  dmc7cuts->SetThreshold(2.5);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DMCL0", dmc7cuts));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigurePP5TeV2017(){
  std::cout << "AutoConfigure: Configuring for data pp 5.02 TeV (2017)" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg1cuts->SetUseRecalcPatches(true);
  eg1cuts->SetThreshold(115);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  eg2cuts->SetUseRecalcPatches(true);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(51);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *dg1cuts = new AliEmcalTriggerSelectionCuts;
  dg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  dg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dg1cuts->SetUseRecalcPatches(true);
  dg1cuts->SetThreshold(115);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG1", dg1cuts));

  AliEmcalTriggerSelectionCuts *dg2cuts = new AliEmcalTriggerSelectionCuts;
  dg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  dg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dg2cuts->SetUseRecalcPatches(true);
  dg2cuts->SetThreshold(51);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG2", dg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej1cuts->SetUseRecalcPatches(true);
  ej1cuts->SetThreshold(255);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  ej2cuts->SetUseRecalcPatches(true);
  ej2cuts->SetThreshold(204);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));

  AliEmcalTriggerSelectionCuts *dj1cuts = new AliEmcalTriggerSelectionCuts;
  dj1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  dj1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dj1cuts->SetUseRecalcPatches(true);
  dj1cuts->SetThreshold(255);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ1", dj1cuts));

  AliEmcalTriggerSelectionCuts *dj2cuts = new AliEmcalTriggerSelectionCuts;
  dj2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  dj2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kADC);
  dj2cuts->SetUseRecalcPatches(true);
  dj2cuts->SetThreshold(204);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ2", dj2cuts));

}

void AliAnalysisTaskEmcalTriggerSelection::ConfigureMCPP5TeV2017() {
  std::cout << "AutoConfigure: Configuring for MC pp 5.02 TeV (2017)" << std::endl;
  AliEmcalTriggerSelectionCuts *eg1cuts = new AliEmcalTriggerSelectionCuts;
  eg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  eg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg1cuts->SetUseSimpleOfflinePatches(true);
  eg1cuts->SetThreshold(9.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG1", eg1cuts));

  AliEmcalTriggerSelectionCuts *eg2cuts = new AliEmcalTriggerSelectionCuts;
  eg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  eg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  eg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  eg2cuts->SetUseSimpleOfflinePatches(true);
  eg2cuts->SetThreshold(4.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EG2", eg2cuts));

  AliEmcalTriggerSelectionCuts *dg1cuts = new AliEmcalTriggerSelectionCuts;
  dg1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaHighPatch);
  dg1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dg1cuts->SetUseSimpleOfflinePatches(true);
  dg1cuts->SetThreshold(9.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG1", dg1cuts));

  AliEmcalTriggerSelectionCuts *dg2cuts = new AliEmcalTriggerSelectionCuts;
  dg2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dg2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1GammaLowPatch);
  dg2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dg2cuts->SetUseSimpleOfflinePatches(true);
  dg2cuts->SetThreshold(4.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DG2", dg2cuts));

  AliEmcalTriggerSelectionCuts *ej1cuts = new AliEmcalTriggerSelectionCuts;
  ej1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  ej1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej1cuts->SetUseSimpleOfflinePatches(true);
  ej1cuts->SetThreshold(19.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej2cuts->SetUseSimpleOfflinePatches(true);
  ej2cuts->SetThreshold(14.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ2", ej2cuts));

  AliEmcalTriggerSelectionCuts *dj1cuts = new AliEmcalTriggerSelectionCuts;
  dj1cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj1cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetHighPatch);
  dj1cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dj1cuts->SetUseSimpleOfflinePatches(true);
  dj1cuts->SetThreshold(19.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ1", dj1cuts));

  AliEmcalTriggerSelectionCuts *dj2cuts = new AliEmcalTriggerSelectionCuts;
  dj2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kDCALAcceptance);
  dj2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  dj2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  dj2cuts->SetUseSimpleOfflinePatches(true);
  dj2cuts->SetThreshold(14.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ2", dj2cuts));
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigureFromYAML(const char *configfile) {
  using YAMLhandler = PWG::Tools::AliYAMLConfiguration;
  YAMLhandler configuration;
  configuration.AddConfiguration(configfile, "user");
  configuration.Initialize();
  std::string namecontainer, acceptance, patchtype, energydef, energysource;
  std::vector<std::string> triggerclasses;
  configuration.GetProperty("containername", namecontainer);
  configuration.GetProperty("energydef", energydef);
  configuration.GetProperty("energysource", energysource);
  configuration.GetProperty("triggerclasses", triggerclasses);
  bool isOfflineSimple = energysource.find("Offline") != std::string::npos,
       isRecalc = energysource.find("Recalc") != std::string::npos;

  SetGlobalDecisionContainerName(namecontainer.data());

  AliEmcalTriggerSelectionCuts::SelectionMethod_t selectionmethod;
  try {
    selectionmethod = DecodeEnergyDefinition(energydef);
  } catch(ConfigValueException &e) {
    AliErrorStream() << e.what() << " - not processing trigger classes" << std::endl;
    return;
  }
  for(auto t : triggerclasses) {
    double threshold;
    configuration.GetProperty(Form("%s:acceptance", t.data()), acceptance);
    configuration.GetProperty(Form("%s:patchtype", t.data()), patchtype);
    configuration.GetProperty(Form("%s:threshold", t.data()), threshold);

    AliEmcalTriggerSelectionCuts *cuts = new AliEmcalTriggerSelectionCuts;
    try {
      cuts->SetAcceptanceType(DecodeAcceptanceString(acceptance));
      cuts->SetPatchType(DecodePatchTypeString(patchtype));
    } catch(ConfigValueException &e){
      AliErrorStream() << e.what() << " - not adding trigger class " << t << std::endl;
      delete cuts;
      continue;
    }

    cuts->SetSelectionMethod(selectionmethod);
    if(isOfflineSimple) cuts->SetUseSimpleOfflinePatches();
    if(isRecalc) cuts->SetUseRecalcPatches();
    cuts->SetThreshold(threshold);
    this->AddTriggerSelection(new AliEmcalTriggerSelection(t.data(), cuts));
  }
}

AliEmcalTriggerSelectionCuts::AcceptanceType_t AliAnalysisTaskEmcalTriggerSelection::DecodeAcceptanceString(const std::string &acceptancestring){
  std::unordered_map<std::string, AliEmcalTriggerSelectionCuts::AcceptanceType_t> mapacceptance = {
    {"EMCAL", AliEmcalTriggerSelectionCuts::kEMCALAcceptance},
    {"DCAL", AliEmcalTriggerSelectionCuts::kDCALAcceptance}
  };
  auto result = mapacceptance.find(acceptancestring);
  if(result == mapacceptance.end()) throw ConfigValueException("accpetance", acceptancestring.data());
  return result->second;
}

AliEmcalTriggerSelectionCuts::PatchType_t AliAnalysisTaskEmcalTriggerSelection::DecodePatchTypeString(const std::string &patchtypestring) {
  std::unordered_map<std::string, AliEmcalTriggerSelectionCuts::PatchType_t> mappatchtype = {
    {"L1Gamma", AliEmcalTriggerSelectionCuts::kL1GammaPatch},
    {"L1GammaHigh", AliEmcalTriggerSelectionCuts::kL1GammaHighPatch},
    {"L1GammaLow", AliEmcalTriggerSelectionCuts::kL1GammaLowPatch},
    {"L1Jet", AliEmcalTriggerSelectionCuts::kL1JetPatch},
    {"L1JetHigh", AliEmcalTriggerSelectionCuts::kL1JetHighPatch},
    {"L1JetLow", AliEmcalTriggerSelectionCuts::kL1JetLowPatch}
  };
  auto result = mappatchtype.find(patchtypestring);
  if(result == mappatchtype.end()) throw ConfigValueException("accpetance", patchtypestring.data());
  return result->second;
}

AliEmcalTriggerSelectionCuts::SelectionMethod_t AliAnalysisTaskEmcalTriggerSelection::DecodeEnergyDefinition(const std::string &energydefstring){
  std::unordered_map<std::string, AliEmcalTriggerSelectionCuts::SelectionMethod_t> mapenergydef = {
    {"ADC", AliEmcalTriggerSelectionCuts::kADC},
    {"Energy", AliEmcalTriggerSelectionCuts::kEnergyOffline},
    {"EnergyRough", AliEmcalTriggerSelectionCuts::kEnergyRough},
    {"EnergySmeared", AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared}
  };
  auto result = mapenergydef.find(energydefstring);
  if(result == mapenergydef.end()) throw ConfigValueException("accpetance", energydefstring.data());
  return result->second;

}

void AliAnalysisTaskEmcalTriggerSelection::PrintStream(std::ostream &stream) const {
    stream << "Task: " << GetName() << ", name of the output container: " << fGlobalDecisionContainerName << std::endl << std::endl;
    stream << "Trigger classes: " << std::endl;
    for(const auto c : this->fTriggerSelections){
      PWG::EMCAL::AliEmcalTriggerSelection *sel = static_cast<PWG::EMCAL::AliEmcalTriggerSelection *>(c);
      stream << *sel << std::endl;
    }
}

AliAnalysisTaskEmcalTriggerSelection::AliEmcalTriggerSelectionQA::AliEmcalTriggerSelectionQA():
    TNamed(),
    fMaxPatchADC(nullptr),
    fMaxPatchEnergy(nullptr),
    fMaxPatchEnergySmeared(nullptr)
{
}

AliAnalysisTaskEmcalTriggerSelection::AliEmcalTriggerSelectionQA::AliEmcalTriggerSelectionQA(const AliEmcalTriggerSelection * const sel):
    TNamed(sel->GetName(), ""),
    fMaxPatchADC(nullptr),
    fMaxPatchEnergy(nullptr),
    fMaxPatchEnergySmeared(nullptr)
{
  fMaxPatchADC = new TH1D(Form("hMaxPatchADC%s", GetName()), "Max. patch ADC", 1000, 0., 1000);
  fMaxPatchEnergy = new TH1D(Form("hMaxPatchEnergy%s", GetName()), "Max. patch energy", 1000, 0., 100);
  fMaxPatchEnergySmeared = new TH1D(Form("hMaxPatchEnergySmeared%s", GetName()), "Max. patch smeared energy", 1000, 0., 100);
}

AliAnalysisTaskEmcalTriggerSelection::AliEmcalTriggerSelectionQA::AliEmcalTriggerSelectionQA(const AliEmcalTriggerSelectionQA &ref):
    TNamed(ref),
    fMaxPatchADC(ref.fMaxPatchADC),
    fMaxPatchEnergy(ref.fMaxPatchEnergy),
    fMaxPatchEnergySmeared(ref.fMaxPatchEnergySmeared)
{
}

AliAnalysisTaskEmcalTriggerSelection::AliEmcalTriggerSelectionQA &AliAnalysisTaskEmcalTriggerSelection::AliEmcalTriggerSelectionQA::operator=(const AliAnalysisTaskEmcalTriggerSelection::AliEmcalTriggerSelectionQA &ref) {
  TNamed::operator=(ref);
  if(this != &ref) {
    fMaxPatchADC = ref.fMaxPatchADC;
    fMaxPatchEnergy = ref.fMaxPatchEnergy;
    fMaxPatchEnergySmeared = ref.fMaxPatchEnergySmeared;
  }
  return *this;
}

void AliAnalysisTaskEmcalTriggerSelection::AliEmcalTriggerSelectionQA::Fill(const AliEmcalTriggerDecision * const decision){
  if(decision->GetMainPatch()){
    fMaxPatchADC->Fill(decision->GetMainPatch()->GetADCAmp());
    fMaxPatchEnergy->Fill(decision->GetMainPatch()->GetPatchE());
    fMaxPatchEnergySmeared->Fill(decision->GetMainPatch()->GetSmearedEnergy());
  }
}

void AliAnalysisTaskEmcalTriggerSelection::AliEmcalTriggerSelectionQA::GetHistos(TList *targetlist) const {
  targetlist->Add(fMaxPatchADC);
  targetlist->Add(fMaxPatchEnergy);
  targetlist->Add(fMaxPatchEnergySmeared);

}

AliAnalysisTaskEmcalTriggerSelection::ConfigValueException::ConfigValueException(const char *key, const char *value):
  fKey(key),
  fValue(value),
  fMessage()
{
  std::stringstream msgbuilder;
  msgbuilder << "Improper value for key " << fKey << ": " << fValue;
  fMessage = msgbuilder.str();
}

}
}

std::ostream &operator<<(std::ostream &stream, const PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection &task) {
  task.PrintStream(stream);
  return stream;
}
