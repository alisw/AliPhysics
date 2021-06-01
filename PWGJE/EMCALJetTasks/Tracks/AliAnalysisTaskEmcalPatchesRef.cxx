/************************************************************************************
 * Copyright (C) 2015, Copyright Holders of the ALICE Collaboration                 *
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

#include <array>
#include <iostream>
#include <map>
#include <vector>

#include <TClonesArray.h>
#include <TGrid.h>
#include <THistManager.h>
#include <THashList.h>
#include <TLinearBinning.h>
#include <TObjArray.h>
#include <TParameter.h>

#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliOADBContainer.h"

#include "AliAnalysisTaskEmcalPatchesRef.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalPatchesRef)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef() :
    AliAnalysisTaskEmcalTriggerBase(),
    fCentralityRange(-999., 999.),
    fCellTimeRange(-1., 1.),
    fOfflineTriggerData(),
    fMinNumberFastors(0),
    fEnableSumw2(false),
    fUseRecalcPatches(false),
    fRequestCentrality(false),
    fEventCentrality(0)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetNeedEmcalGeom(true);
}

AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef(const char *name):
    AliAnalysisTaskEmcalTriggerBase(name),
    fCentralityRange(-999., 999.),
    fCellTimeRange(-1., 1.),
    fOfflineTriggerData(),
    fMinNumberFastors(0),
    fEnableSumw2(false),
    fUseRecalcPatches(false),
    fRequestCentrality(false),
    fEventCentrality(0)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetNeedEmcalGeom(true);
}

void AliAnalysisTaskEmcalPatchesRef::CreateUserHistos(){
  AliInfoStream() <<  "Creating histograms for task " << GetName() << std::endl;

  EnergyBinning energybinning;
  TLinearBinning etabinning(100, -0.7, 0.7);
  const std::array<const TString, 10> patchtypes = {"EG1", "EG2", "EJ1", "EJ2", "EMC7", "DG1", "DG2", "DJ1", "DJ2", "DMC7"};
  const std::array<double, 5> encuts = {1., 2., 5., 10., 20.};
  TString optionstring = fEnableSumw2 ? "s" : "";
  for(const auto &trg : GetSupportedTriggers()){
    fHistos->CreateTH1("EventCount" + trg, "Event count for trigger class " + trg, 1, 0.5, 1.5, optionstring);
    fHistos->CreateTH1("EventCentrality" + trg, "Event centrality for trigger class " + trg, 103, -2., 101., optionstring);
    fHistos->CreateTH1("VertexZ" + trg, "z-position of the primary vertex for trigger class " + trg, 200, -40., 40., optionstring);
    for(const auto &patch : patchtypes){
      fHistos->CreateTH1(patch + "PatchEnergy" + trg,  patch + "-patch energy for trigger class " + trg, energybinning, optionstring);
      fHistos->CreateTH1(patch + "PatchET" + trg,  patch +"-patch transverse energy for trigger class "+ trg, energybinning, optionstring);
      fHistos->CreateTH2(patch + "PatchEnergyEsmear" + trg, patch + "-patch energy vs. smeared energy for trigger class " + trg, energybinning, energybinning);
      fHistos->CreateTH2(patch + "PatchEnergyEta" + trg, patch + "%s-patch energy for trigger class " + trg, energybinning, etabinning, optionstring);
      fHistos->CreateTH2(patch + "PatchETEta" +trg, patch + "-patch transverse energy for trigger class " + trg, energybinning, etabinning, optionstring);
      for(auto energy : encuts){
        fHistos->CreateTH2(Form("%sEtaPhi%dG%s", patch.Data(), static_cast<int>(energy), trg.Data()), Form("%s-patch #eta-#phi map for patches with energy larger than %f GeV/c for trigger class %s", patch.Data(), energy, trg.Data()), 100, -0.7, 0.7, 200, 0, TMath::TwoPi(), optionstring);
        fHistos->CreateTH2(Form("%sColRow%dG%s", patch.Data(), static_cast<int>(energy), trg.Data()), Form("%s-patch col-row map for patches with energy larger than %f GeV/c for trigger class %s", patch.Data(), energy, trg.Data()), 48, -0.5, 47.5, 104, -0.5, 103.5, optionstring);
      }
    }
  }

  // creating trigger data grid
  fOfflineTriggerData.Allocate(48, 104);

  AliDebugStream(1) << "Histograms done" << std::endl;
}

bool AliAnalysisTaskEmcalPatchesRef::IsUserEventSelected(){
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
  return true;
}

bool AliAnalysisTaskEmcalPatchesRef::Run(){
  AliDebugStream(1) << GetName() << ": Start function" << std::endl;
  FillDataGrid();

  AliDebugStream(1) << GetName() << ": Number of trigger patches " << fTriggerPatchInfo->GetEntries() << std::endl;

  Double_t energy, eta, phi, et, smearedenergy;
  Int_t col, row;
  for(auto patchIter : *fTriggerPatchInfo){
    AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(GetNumberOfFastORs(patch) < fMinNumberFastors) continue;

    bool isDCAL         = patch->IsDCalPHOS(),
        isSingleShower  = SelectSingleShowerPatch(patch),
        isJetPatch      = SelectJetPatch(patch);

    std::vector<TString> patchnames;
    if(isJetPatch){
      if(isDCAL){
        patchnames.push_back("DJ1");
        patchnames.push_back("DJ2");
      } else {
        patchnames.push_back("EJ1");
        patchnames.push_back("EJ2");
      }
    }
    if(isSingleShower){
      if(isDCAL){
        patchnames.push_back("DMC7");
        patchnames.push_back("DG1");
        patchnames.push_back("DG2");
      } else {
        patchnames.push_back("EMC7");
        patchnames.push_back("EG1");
        patchnames.push_back("EG2");
      }
    }
    if(!patchnames.size()){
      // Undefined patch type - ignore
      continue;
    }

    TLorentzVector posvec;
    energy = fUseRecalcPatches ? patch->GetADCAmpGeVRough() : patch->GetPatchE();
    smearedenergy = patch->GetSmearedEnergy();
    eta = patch->GetEtaGeo();
    phi = patch->GetPhiGeo();
    col = patch->GetColStart();
    row = patch->GetRowStart();
    et = patch->GetLorentzVectorCenterGeo().Et();

    // fill histograms allEta
    for(const auto &nameit : patchnames){
      for(const auto &trg : fSelectedTriggers){
        FillPatchHistograms(trg.Data(), nameit, energy, et, smearedenergy, eta, phi, col, row);
      }
    }
  }
  return true;
}

void AliAnalysisTaskEmcalPatchesRef::FillPatchHistograms(TString triggerclass, TString patchname, double energy, double transverseenergy, double smearedenergy, double eta, double phi, int col, int row){
  Double_t weight = GetTriggerWeight(triggerclass.Data());
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << triggerclass << " in patch histograms." << std::endl;
  fHistos->FillTH1(patchname + "PatchEnergy" + triggerclass, energy, weight);
  fHistos->FillTH1(patchname + "PatchET" + triggerclass, transverseenergy, weight);
  fHistos->FillTH2(patchname + "PatchEnergyEta" + triggerclass, energy, eta, weight);
  fHistos->FillTH2(patchname + "PatchETEta" + triggerclass, transverseenergy, eta, weight);
  fHistos->FillTH2(patchname + "PatchEnergyEsmear" + triggerclass, energy, smearedenergy, weight);
  const std::array<double, 5> encuts = {1., 2., 5., 10., 20.};
  for(auto etest : encuts){
    if(energy > etest){
      fHistos->FillTH2(Form("%sEtaPhi%dG%s", patchname.Data(), static_cast<int>(etest), triggerclass.Data()), eta, phi, weight);
      fHistos->FillTH2(Form("%sColRow%dG%s", patchname.Data(), static_cast<int>(etest), triggerclass.Data()), col, row, weight);
    }
  }
}

void AliAnalysisTaskEmcalPatchesRef::UserFillHistosAfterEventSelection(){
  // Fill Event counter and reference vertex distributions for the different trigger classes
  for(const auto &trg : fSelectedTriggers){
    Double_t weight = GetTriggerWeight(trg.Data());
    AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << trg << " in event histograms." << std::endl;
    fHistos->FillTH1("EventCount" + trg, 1, weight);
    fHistos->FillTH1("EventCentrality" + trg, fEventCentrality, weight);
    fHistos->FillTH1("VertexZ" + trg, fVertex[2], weight);
  }

}

void AliAnalysisTaskEmcalPatchesRef::GetPatchBoundaries(const AliEMCALTriggerPatchInfo *patch, Double_t *boundaries) const {
  boundaries[0] = patch->GetEtaMin();
  boundaries[1] = patch->GetEtaMax();
  boundaries[2] = patch->GetPhiMin();
  boundaries[3] = patch->GetPhiMax();
}

bool AliAnalysisTaskEmcalPatchesRef::SelectSingleShowerPatch(const AliEMCALTriggerPatchInfo *patch) const{
  if(fUseRecalcPatches){
    if(!patch->IsRecalc()) return false;
    return patch->IsGammaLowRecalc();
  } else {
    if(!patch->IsOfflineSimple()) return false;
    return patch->IsGammaLowSimple();
  }
}

bool AliAnalysisTaskEmcalPatchesRef::SelectJetPatch(const AliEMCALTriggerPatchInfo *patch) const{
  if(fUseRecalcPatches){
    if(!patch->IsRecalc()) return false;
    return patch->IsJetLowRecalc();
  } else {
    if(!patch->IsOfflineSimple()) return false;
    return patch->IsJetLowSimple();
  }
}

void AliAnalysisTaskEmcalPatchesRef::FillDataGrid(){
  fOfflineTriggerData.Reset();
  Short_t cellID;
  Double_t cellamplitude, celltime, efrac;
  Int_t mclabel, fastorID, row, col;
  for(int icell = 0; icell < fCaloCells->GetNumberOfCells(); icell++) {
    fCaloCells->GetCell(icell, cellID, cellamplitude, celltime, mclabel, efrac);
    if(!fCellTimeRange.IsInRange(celltime)) continue;
    if(cellamplitude > 0.) {
      fGeom->GetFastORIndexFromCellIndex(cellID, fastorID);
      fGeom->GetPositionInEMCALFromAbsFastORIndex(fastorID, col, row);
      fOfflineTriggerData(col, row) += cellamplitude;
    } 
  }
}

Int_t AliAnalysisTaskEmcalPatchesRef::GetNumberOfFastORs(const AliEMCALTriggerPatchInfo *patch) const {
  Int_t nfastor = 0;
  for(int icol = patch->GetColStart(); icol < patch->GetColStart() + patch->GetPatchSize(); icol++) {
    for(int irow = patch->GetRowStart(); irow < patch->GetRowStart() + patch->GetPatchSize(); irow++) {
      if(fOfflineTriggerData(icol, irow) > 0) nfastor++;
    }
  }
  return nfastor;
}
AliAnalysisTaskEmcalPatchesRef::EnergyBinning::EnergyBinning():
    TCustomBinning()
{
  this->SetMinimum(0.);
  this->AddStep(1., 0.05);
  this->AddStep(2., 0.1);
  this->AddStep(4, 0.2);
  this->AddStep(7, 0.5);
  this->AddStep(16, 1);
  this->AddStep(32, 2);
  this->AddStep(40, 4);
  this->AddStep(50, 5);
  this->AddStep(100, 10);
  this->AddStep(200, 20);
}
