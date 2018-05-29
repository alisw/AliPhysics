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
#include <string>
#include <sstream>
#include <vector>
#include <THistManager.h>
#include <TCustomBinning.h>
#include <TLinearBinning.h>

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetEnergyScale.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliLog.h"
#include "AliVEventHandler.h"

ClassImp(EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergyScale)

using namespace EmcalTriggerJets;

AliAnalysisTaskEmcalJetEnergyScale::AliAnalysisTaskEmcalJetEnergyScale():
  AliAnalysisTaskEmcalJet(),
  fHistos(nullptr),
  fNameDetectorJets(),
  fNameParticleJets(),
  fTriggerSelectionString(),
  fNameTriggerDecisionContainer("EmcalTriggerDecision")
{
}

AliAnalysisTaskEmcalJetEnergyScale::AliAnalysisTaskEmcalJetEnergyScale(const char *name):
  AliAnalysisTaskEmcalJet(name, true),
  fHistos(nullptr),
  fNameDetectorJets(),
  fNameParticleJets(),
  fTriggerSelectionString(),
  fNameTriggerDecisionContainer("EmcalTriggerDecision")
{
  SetUseAliAnaUtils(true);
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEmcalJetEnergyScale::~AliAnalysisTaskEmcalJetEnergyScale() {
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalJetEnergyScale::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  TLinearBinning jetPtBinning(20, 0., 200.), nefbinning(100, 0., 1.), ptdiffbinning(200, -1., 1.), jetEtaBinning(100, -0.9, 0.9), jetPhiBinning(100, 0., TMath::TwoPi());

  const TBinning *diffbinning[3] = {&jetPtBinning, &nefbinning, &ptdiffbinning},
                 *corrbinning[3] = {&jetPtBinning, &jetPtBinning, &nefbinning},
                 *effbinning[3] = {&jetPtBinning, &jetEtaBinning, &jetPhiBinning};

  fHistos = new THistManager("energyScaleHistos");
  fHistos->CreateTH1("hEventCounter", "Event counter", 1, 0.5, 1.5);
  fHistos->CreateTHnSparse("hPtDiff", "pt diff det/part", 3., diffbinning, "s");
  fHistos->CreateTHnSparse("hPtCorr", "Correlation det pt / part pt", 3., corrbinning, "s");
  fHistos->CreateTHnSparse("hPartJetsAccepted", "Accepted particle level jets", 3, effbinning, "s");
  fHistos->CreateTHnSparse("hPartJetsReconstructed", "Accepted and reconstructed particle level jets", 3, effbinning, "s");
  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalJetEnergyScale::Run(){
  AliDebugStream(1) << "Next event" << std::endl;
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;
  if(IsSelectEmcalTriggers(fTriggerSelectionString.Data())){
    auto mctrigger = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer));
    AliDebugStream(1) << "Found trigger decision object: " << (mctrigger ? "yes" : "no") << std::endl;
    if(!mctrigger){
      AliErrorStream() <<  "Trigger decision container with name " << fNameTriggerDecisionContainer << " not found in event - not possible to select EMCAL triggers" << std::endl;
      return false;
    }
    if(!mctrigger->IsEventSelected(fTriggerSelectionString)) return false;
  }
  AliDebugStream(1) << "event selected" << std::endl;
  fHistos->FillTH1("hEventCounter", 1);

  auto detjets = GetJetContainer(fNameDetectorJets),
       partjets = GetJetContainer(fNameParticleJets);
  if(!detjets || !partjets) {
    AliErrorStream() << "At least one jet container missing, exiting ..." << std::endl;
    return false;
  }
  AliDebugStream(1) << "Have both jet containers: part(" << partjets->GetNAcceptedJets() << "|" << partjets->GetNJets() << "), det(" << detjets->GetNAcceptedJets() << "|" << detjets->GetNJets() << ")" << std::endl;

  std::vector<AliEmcalJet *> taggedjets;
  for(auto detjet : detjets->accepted()){
    AliDebugStream(2) << "Next jet" << std::endl;
    auto partjet = detjet->ClosestJet();
    if(!partjet) {
      AliDebugStream(2) << "No tagged jet" << std::endl;
      continue;
    }
    taggedjets.emplace_back(partjet);
    double pointCorr[3] = {partjet->Pt(), detjet->Pt(), detjet->NEF()},
           pointDiff[3] = {partjet->Pt(), (detjet->Pt()-partjet->Pt())/partjet->Pt(), detjet->NEF()};
    fHistos->FillTHnSparse("hPtDiff", pointDiff);
    fHistos->FillTHnSparse("hPtCorr", pointCorr);
  }

  // efficiency x acceptance: Add histos for all accepted and reconstucted accepted jets
  for(auto partjet : partjets->accepted()){
    double pvect[3] = {partjet->Pt(), partjet->Eta(), partjet->Phi()};
    if(pvect[2] < 0) pvect[2] += TMath::TwoPi();
    fHistos->FillTHnSparse("hPartJetsAccepted", pvect);
    if(std::find(taggedjets.begin(), taggedjets.end(), partjet) != taggedjets.end()) fHistos->FillTHnSparse("hPartJetsReconstructed", pvect);
  }
  return true;
}

bool AliAnalysisTaskEmcalJetEnergyScale::IsSelectEmcalTriggers(const TString &triggerstring) const {
  const std::array<TString, 8> kEMCALTriggers = {
    "EJ1", "EJ2", "DJ1", "DJ2", "EG1", "EG2", "DG1", "DG2"
  };
  bool isEMCAL = false;
  for(auto emcaltrg : kEMCALTriggers) {
    if(triggerstring.Contains(emcaltrg)) {
      isEMCAL = true;
      break;
    }
  }
  return isEMCAL;
}

AliAnalysisTaskEmcalJetEnergyScale *AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale(AliJetContainer::EJetType_t jettype, Double_t jetradius, Bool_t useDCAL, const char *trigger) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale", "No analysis manager available");
    return nullptr;
  } 

  auto inputhandler = mgr->GetInputEventHandler();
  auto isAOD = inputhandler->IsA() == AliAODInputHandler::Class();

  std::string jettypename;
  AliJetContainer::JetAcceptanceType acceptance(AliJetContainer::kTPCfid);
  bool addClusterContainer(false), addTrackContainer(false);
  switch(jettype){
    case AliJetContainer::kFullJet:
        jettypename = "FullJet";
        acceptance = useDCAL ? AliJetContainer::kDCALfid : AliJetContainer::kEMCALfid;
        addClusterContainer = addTrackContainer = true;
        break;
    case AliJetContainer::kChargedJet:
        jettypename = "ChargedJet";
        acceptance = AliJetContainer::kTPCfid;
        addTrackContainer = true;
        break;
    case AliJetContainer::kNeutralJet:
        jettypename = "NeutralJet";
        acceptance = useDCAL ? AliJetContainer::kDCALfid : AliJetContainer::kEMCALfid;
        addClusterContainer = true;
        break;
  };

  std::stringstream taskname, tag;
  tag << jettypename << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  taskname << "EnergyScaleTask_" << tag.str();
  AliAnalysisTaskEmcalJetEnergyScale *energyscaletask = new AliAnalysisTaskEmcalJetEnergyScale(taskname.str().data());
  mgr->AddTask(energyscaletask);
  energyscaletask->SetTriggerName(trigger);

  auto partcont = energyscaletask->AddMCParticleContainer("mcparticles");
  partcont->SetMinPt(0.);

  AliClusterContainer *clusters(nullptr);
  if(addClusterContainer) {
    clusters = energyscaletask->AddClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    switch(jettype){
      case AliJetContainer::kChargedJet: break;     // Silence compiler
      case AliJetContainer::kFullJet:
        clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
        clusters->SetClusHadCorrEnergyCut(0.3);
        break;
      case AliJetContainer::kNeutralJet:
        clusters->SetDefaultClusterEnergy(AliVCluster::kNonLinCorr);
        clusters->SetClusNonLinCorrEnergyCut(0.3);
        break;
    };
  }
  AliTrackContainer *tracks(nullptr);
  if(addTrackContainer) {
    tracks = energyscaletask->AddTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
  }

  auto contpartjet = energyscaletask->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, jetradius,
                                                      AliJetContainer::kTPCfid, partcont, nullptr);
  contpartjet->SetName("particleLevelJets");
  energyscaletask->SetNamePartJetContainer("particleLevelJets");
  std::cout << "Adding particle-level jet container with underling array: " << contpartjet->GetArrayName() << std::endl;

  auto contdetjet = energyscaletask->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, jetradius,
                                                     acceptance, tracks, clusters);
  contdetjet->SetName("detectorLevelJets");
  energyscaletask->SetNameDetJetContainer("detectorLevelJets");
  std::cout << "Adding detector-level jet container with underling array: " << contdetjet->GetArrayName() << std::endl;

  std::stringstream outnamebuilder, listnamebuilder;
  listnamebuilder << "EnergyScaleHists_" << tag.str();
  outnamebuilder << mgr->GetCommonFileName() << ":EnergyScaleResults_" << tag.str();

  mgr->ConnectInput(energyscaletask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(energyscaletask, 1, mgr->CreateContainer(listnamebuilder.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, outnamebuilder.str().data()));
  return energyscaletask;
} 