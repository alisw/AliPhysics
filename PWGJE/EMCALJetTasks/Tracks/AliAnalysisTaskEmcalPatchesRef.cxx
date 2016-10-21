/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef)
/// \endcond

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef() :
    AliAnalysisTaskEmcalTriggerBase(),
    fCentralityRange(-999., 999.),
    fRequestCentrality(false),
    fEventCentrality(0)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalPatchesRef::AliAnalysisTaskEmcalPatchesRef(const char *name):
    AliAnalysisTaskEmcalTriggerBase(name),
    fCentralityRange(-999., 999.),
    fRequestCentrality(false),
    fEventCentrality(0)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

void AliAnalysisTaskEmcalPatchesRef::CreateUserHistos(){
  AliInfoStream() <<  "Creating histograms for task " << GetName() << std::endl;

  EnergyBinning energybinning;
  TLinearBinning etabinning(100, -0.7, 0.7);
  std::array<TString, 10> patchtypes = {"EG1", "EG2", "EJ1", "EJ2", "EMC7", "DG1", "DG2", "DJ1", "DJ2", "DMC7"};
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(auto trg : GetSupportedTriggers()){
    fHistos->CreateTH1(Form("hEventCount%s", trg.Data()), Form("Event count for trigger class %s", trg.Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hEventCentrality%s", trg.Data()), Form("Event centrality for trigger class %s", trg.Data()), 103, -2., 101.);
    fHistos->CreateTH1(Form("hVertexZ%s", trg.Data()), Form("z-position of the primary vertex for trigger class %s", trg.Data()), 200, -40., 40.);
    for(auto patch : patchtypes){
      fHistos->CreateTH1(Form("h%sPatchEnergy%s", patch.Data(), trg.Data()), Form("%s-patch energy for trigger class %s", patch.Data(), trg.Data()), energybinning);
      fHistos->CreateTH1(Form("h%sPatchET%s", patch.Data(), trg.Data()), Form("%s-patch transverse energy for trigger class %s", patch.Data(), trg.Data()), energybinning);
      fHistos->CreateTH2(Form("h%sPatchEnergyEta%s", patch.Data(), trg.Data()), Form("%s-patch energy for trigger class %s", patch.Data(), trg.Data()), energybinning, etabinning);
      fHistos->CreateTH2(Form("h%sPatchETEta%s", patch.Data(), trg.Data()), Form("%s-patch transverse energy for trigger class %s", patch.Data(), trg.Data()), energybinning, etabinning);
      for(int ien = 0; ien < 5; ien++){
        fHistos->CreateTH2(Form("h%sEtaPhi%dG%s", patch.Data(), static_cast<int>(encuts[ien]), trg.Data()), Form("%s-patch #eta-#phi map for patches with energy larger than %f GeV/c for trigger class %s", patch.Data(), encuts[ien], trg.Data()), 100, -0.7, 0.7, 200, 0, TMath::TwoPi());
        fHistos->CreateTH2(Form("h%sColRow%dG%s", patch.Data(), static_cast<int>(encuts[ien]), trg.Data()), Form("%s-patch col-row map for patches with energy larger than %f GeV/c for trigger class %s", patch.Data(), encuts[ien], trg.Data()), 48, -0.5, 47.5, 104, -0.5, 103.5);
      }
    }
  }
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

  AliDebugStream(1) << GetName() << ": Number of trigger patches " << fTriggerPatchInfo->GetEntries() << std::endl;

  Double_t energy, eta, phi, et;
  Int_t col, row;
  for(auto patchIter : *fTriggerPatchInfo){
    AliEMCALTriggerPatchInfo *patch = static_cast<AliEMCALTriggerPatchInfo *>(patchIter);
    if(!patch->IsOfflineSimple()) continue;

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
    energy = patch->GetPatchE();
    eta = patch->GetEtaGeo();
    phi = patch->GetPhiGeo();
    col = patch->GetColStart();
    row = patch->GetRowStart();
    et = patch->GetLorentzVectorCenterGeo().Et();

    // fill histograms allEta
    for(const auto &nameit : patchnames){
      for(const auto &trg : fSelectedTriggers){
        FillPatchHistograms(trg.Data(), nameit, energy, et, eta, phi, col, row);
      }
    }
  }
}

void AliAnalysisTaskEmcalPatchesRef::FillPatchHistograms(TString triggerclass, TString patchname, double energy, double transverseenergy, double eta, double phi, int col, int row){
  Double_t weight = GetTriggerWeight(triggerclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << triggerclass << " in patch histograms." << std::endl;
  fHistos->FillTH1(Form("h%sPatchEnergy%s", patchname.Data(), triggerclass.Data()), energy, weight);
  fHistos->FillTH1(Form("h%sPatchET%s", patchname.Data(), triggerclass.Data()), transverseenergy, weight);
  fHistos->FillTH2(Form("h%sPatchEnergyEta%s", patchname.Data(), triggerclass.Data()), energy, eta, weight);
  fHistos->FillTH2(Form("h%sPatchETEta%s", patchname.Data(), triggerclass.Data()), transverseenergy, eta, weight);
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(int ien = 0; ien < 5; ien++){
    if(energy > encuts[ien]){
      fHistos->FillTH2(Form("h%sEtaPhi%dG%s", patchname.Data(), static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi, weight);
      fHistos->FillTH2(Form("h%sColRow%dG%s", patchname.Data(), static_cast<int>(encuts[ien]), triggerclass.Data()), col, row, weight);
    }
  }
}

void AliAnalysisTaskEmcalPatchesRef::UserFillHistosAfterEventSelection(){
  // Fill Event counter and reference vertex distributions for the different trigger classes
  for(const auto &trg : fSelectedTriggers){
    Double_t weight = GetTriggerWeight(trg);
    AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << trg << " in event histograms." << std::endl;
    fHistos->FillTH1(Form("hEventCount%s", trg.Data()), 1, weight);
    fHistos->FillTH1(Form("hEventCentrality%s", trg.Data()), fEventCentrality, weight);
    fHistos->FillTH1(Form("hVertexZ%s", trg.Data()), fVertex[2], weight);
  }

}

void AliAnalysisTaskEmcalPatchesRef::GetPatchBoundaries(const AliEMCALTriggerPatchInfo *patch, Double_t *boundaries) const {
  boundaries[0] = patch->GetEtaMin();
  boundaries[1] = patch->GetEtaMax();
  boundaries[2] = patch->GetPhiMin();
  boundaries[3] = patch->GetPhiMax();
}

bool AliAnalysisTaskEmcalPatchesRef::SelectSingleShowerPatch(const AliEMCALTriggerPatchInfo *patch) const{
  if(!patch->IsOfflineSimple()) return false;
  return patch->IsGammaLowSimple();
}

bool AliAnalysisTaskEmcalPatchesRef::SelectJetPatch(const AliEMCALTriggerPatchInfo *patch) const{
  if(!patch->IsOfflineSimple()) return false;
  return patch->IsJetLowSimple();
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


} /* namespace EMCalTriggerPtAnalysis */
