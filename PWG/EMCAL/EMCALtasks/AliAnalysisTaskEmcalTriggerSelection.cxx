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
#include <vector>
#include <TH1.h>
#include "AliEmcalTriggerDecision.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalTriggerSelection.h"
#include "AliEmcalTriggerSelectionCuts.h"
#include "AliAnalysisTaskEmcalTriggerSelection.h"
#include "AliEMCALTriggerPatchInfo.h"

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection)
/// \endcond

namespace PWG {
namespace EMCAL {

AliAnalysisTaskEmcalTriggerSelection::AliAnalysisTaskEmcalTriggerSelection():
  AliAnalysisTaskEmcal(),
  fTriggerDecisionContainer(nullptr),
  fGlobalDecisionContainerName("EmcalTriggerDecision"),
  fTriggerSelections(),
  fSelectionQA()
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
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
}

void AliAnalysisTaskEmcalTriggerSelection::AddTriggerSelection(AliEmcalTriggerSelection * const selection){
  fTriggerSelections.Add(selection);
}

Bool_t AliAnalysisTaskEmcalTriggerSelection::Run(){
  fTriggerDecisionContainer->Reset();
  AliEmcalTriggerSelection *selection(NULL);
  TIter selectionIter(&fTriggerSelections);
  while((selection = dynamic_cast<AliEmcalTriggerSelection *>(selectionIter()))){
    fTriggerDecisionContainer->AddTriggerDecision(selection->MakeDecison(fTriggerPatchInfo));
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

void AliAnalysisTaskEmcalTriggerSelection::AutoConfigure(const char *period) {
  std::vector<TString> pp2016periods = {"LHC16h", "LHC16i", "LHC16j", "LHC16k", "LHC16l", "LHC16o", "LHC16p"};
  std::vector<TString> mcpp2016periods = {"LHC17f8", "LHC17f8a", "LHC17f8b", "LHC178c", "LHC17f8d", "LHC17f8e",
                                          "LHC17f8f", "LHC17f8g", "LHC17f8h", "LHC17f8i", "LHC17f8j", "LHC17f8k"};
  TString periodstring(period);
  if(std::find(pp2016periods.begin(), pp2016periods.end(), periodstring) != pp2016periods.end()) ConfigurePP2016();
  if(std::find(mcpp2016periods.begin(), mcpp2016periods.end(), periodstring) != mcpp2016periods.end()) ConfigureMCPP2016();
}

void AliAnalysisTaskEmcalTriggerSelection::ConfigurePP2016(){
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
  ej1cuts->SetThreshold(9.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("EJ1", ej1cuts));

  AliEmcalTriggerSelectionCuts *ej2cuts = new AliEmcalTriggerSelectionCuts;
  ej2cuts->SetAcceptanceType(AliEmcalTriggerSelectionCuts::kEMCALAcceptance);
  ej2cuts->SetPatchType(AliEmcalTriggerSelectionCuts::kL1JetLowPatch);
  ej2cuts->SetSelectionMethod(AliEmcalTriggerSelectionCuts::kEnergyOfflineSmeared);
  ej2cuts->SetUseSimpleOfflinePatches(true);
  ej2cuts->SetThreshold(4.);
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
  dj2cuts->SetThreshold(4.);
  this->AddTriggerSelection(new AliEmcalTriggerSelection("DJ2", dj2cuts));
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

}
}
