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
#include <sstream>
#include <vector>

#include <TLorentzVector.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskForwardJets.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliLog.h"
#include "AliVVZERO.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskForwardJets)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskForwardJets::AliAnalysisTaskForwardJets() : AliAnalysisTaskEmcalJet() {

}

AliAnalysisTaskForwardJets::AliAnalysisTaskForwardJets(const char *name):
  AliAnalysisTaskEmcalJet(name, true)
{
  this->SetNeedEmcalGeom(true);
}

AliAnalysisTaskForwardJets::~AliAnalysisTaskForwardJets() {
  if(fHistos) {
    delete fHistos;
  }
}

void AliAnalysisTaskForwardJets::UserCreateOutputObjects(){
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fHistos = new THistManager(Form("Histos_%s", GetName()));
  fHistos->CreateTH1("hJetEtaDistAll", "Jet eta distribution (all jets)", 1000, -7, 7);
  fHistos->CreateTH1("hJetEtaDistSel", "Jet eta distribution (sel jets)", 1000, -7, 7);
  fHistos->CreateTH1("hJetEtaDistSelEvents", "Jet eta distribution (all jets sel events)", 1000, -7, 7);
  fHistos->CreateTH1("hSumMultV0A", "Integrated V0A multiplicity", 1000, 0., 10000.);
  fHistos->CreateTH1("hSumMultV0C", "Integrated V0C multiplicity", 1000, 0., 10000.);
  fHistos->CreateTH1("hMaxMultV0A", "Max V0A multiplicity", 1000, 0., 10000.);
  fHistos->CreateTH1("hMaxMultV0C", "Max V0C multiplicity", 1000, 0., 10000.);
  fHistos->CreateTH1("hMeanMultV0A", "Mean V0A multiplicity", 1000, 0., 10000.);
  fHistos->CreateTH1("hMeanMultV0C", "Mean V0C multiplicity", 1000, 0., 10000.);
  fHistos->CreateTH2("hCorrJetESumMultV0A", "Jet energy vs. integrated V0A multiplicity", 200, 0., 200, 1000, 0., 10000.);
  fHistos->CreateTH2("hCorrJetESumMultV0C", "Jet energy vs. integrated V0C multiplicity", 200, 0., 200, 1000, 0., 10000.);
  fHistos->CreateTH2("hCorrJetEMaxMultV0A", "Jet energy vs. max V0A channel multiplicity", 200, 0., 200, 1000, 0., 10000.);
  fHistos->CreateTH2("hCorrJetEMaxMultV0C", "Jet energy vs. max V0C channel multiplicity", 200, 0., 200, 1000, 0., 10000.);
  fHistos->CreateTH2("hCorrJetEMeanMultV0A", "Jet energy vs. mean V0A channel multiplicity", 200, 0., 200, 1000, 0., 10000.);
  fHistos->CreateTH2("hCorrJetEMeanMultV0C", "Jet energy vs. mean V0C channel multiplicity", 200, 0., 200, 1000, 0., 10000.);
  fHistos->CreateTH2("hChannelMultV0A", "Multiplcity per V0A channel", 32, -0.5, 31.5, 1000, 0., 10000.);
  fHistos->CreateTH2("hChannelMultV0C", "Multiplcity per V0C channel", 32, -0.5, 31.5, 1000, 0., 10000.);
  for(int ichan = 0; ichan < 32; ichan++) {
    fHistos->CreateTH2(Form("hCorrJetEMultV0A%d", ichan), Form("Jet energy vs. V0A multiplcity channel %d", ichan), 200, 0., 200, 1000, 0., 10000.);
    fHistos->CreateTH2(Form("hCorrJetEMultV0C%d", ichan), Form("Jet energy vs. V0C multiplcity channel %d", ichan), 200, 0., 200, 1000, 0., 10000.);
  }
  std::vector<std::pair<int, int>> eranges = {{5, 10}, {10, 20}, {20, 50}, {50, 100}, {100, 150}, {150, 200}};
  for(auto erange : eranges) {
    fHistos->CreateTH2(Form("hChannelMultV0A_%d_%d", erange.first, erange.second), Form("Multiplcity per V0A channel for %d GeV < Ejet < %d GeV", erange.first, erange.second), 32, -0.5, 31.5, 1000, 0., 10000.);
    fHistos->CreateTH2(Form("hChannelMultV0C_%d_%d", erange.first, erange.second), Form("Multiplcity per V0C channel for %d GeV < Ejet < %d GeV", erange.first, erange.second), 32, -0.5, 31.5, 1000, 0., 10000.);
  }
  fHistos->CreateTH1("hClusterMultiplicityEMCAL", "EMCAL cluster multiplicity", 200, 0., 200);
  fHistos->CreateTH1("hClusterMultiplicityDCAL", "DCAL cluster multiplicity", 200, 0., 200);
  fHistos->CreateTH1("hClusterMultiplicityEDCAL", "EMCAL+DCAL cluster multiplicity", 200, 0., 200);
  fHistos->CreateTH1("hTotalEnergyEMCAL", "Total energy in EMCAL", 1000, 0., 1000.);
  fHistos->CreateTH1("hTotalEnergyDCAL", "Total energy in DCAL", 1000, 0., 1000.);
  fHistos->CreateTH1("hTotalEnergyEDCAL", "Total energy in EMCAL+DCAL", 1000, 0., 1000.);
  fHistos->CreateTH1("hLeadingClusterEMCAL", "Leading cluster energy in EMCAL", 200, 0., 200);
  fHistos->CreateTH1("hLeadingClusterDCAL", "Leading cluster energy in DCAL", 200, 0., 200);
  fHistos->CreateTH1("hLeadingClusterEDCAL", "Leading cluster energy in EDCAL", 200, 0., 200);

  fHistos->CreateTH2("hCorrJetEClusterMultiplicityEMCAL", "FOCAL jet energy vs. EMCAL cluster multiplicity", 200, 0., 200., 200, 0., 200);
  fHistos->CreateTH2("hCorrJetEClusterMultiplicityDCAL", "FOCAL jet energy vs. DCAL cluster multiplicity", 200, 0., 200., 200, 0., 200);
  fHistos->CreateTH2("hCorrJetEClusterMultiplicityEDCAL", "FOCAL jet energy vs. EMCAL+DCAL cluster multiplicity", 200, 0., 200., 200, 0., 200);
  fHistos->CreateTH2("hCorrJetETotalEnergyEMCAL", "FOCAL jet energy vs. Total energy in EMCAL", 200, 0., 200., 1000, 0., 1000.);
  fHistos->CreateTH2("hCorrJetETotalEnergyDCAL", "FOCAL jet energy vs. Total energy in DCAL", 200, 0., 200., 1000, 0., 1000.);
  fHistos->CreateTH2("hCorrJetETotalEnergyEDCAL", "FOCAL jet energy vs. Total energy in EMCAL+DCAL", 200, 0., 200., 1000, 0., 1000.);
  fHistos->CreateTH2("hCorrJetELeadingClusterEMCAL", "FOCAL jet energy vs. Leading cluster energy in EMCAL", 200, 0., 200., 200, 0., 200);
  fHistos->CreateTH2("hCorrJetELeadingClusterDCAL", "FOCAL jet energy vs. Leading cluster energy in DCAL", 200, 0., 200., 200, 0., 200);
  fHistos->CreateTH2("hCorrJetELeadingClusterEDCAL", "FOCAL jet energy vs. Leading cluster energy in EDCAL", 200, 0., 200., 200, 0., 200);
  if(fDebugMaxJetOutliers) {
    fHistos->CreateTH2("hDebugMaxJetPt", "Debug Detection of max. part jet; p_{t,stl}; p_{t,man}", 350, 0., 350, 350., 0., 350.);
    fHistos->CreateTH1("hDebugMaxJetError", "Number of cases where the max. jet pt differs between methods", 1., 0.5, 1.5);
  }
  for(auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskForwardJets::Run() {
  double elead = 0;
  if(fRequireFOCALJet) {
    const double kMinEtaFOCAL = 3.2, kMaxEtaFOCAL = 5.8;
    auto partjets = this->GetJetContainer(fNamePartJetCont.Data());
    bool acceptedEvent = false;
    for(auto jet : partjets->accepted()) {
      fHistos->FillTH1("hJetEtaDistAll", jet->Eta());
      if(jet->Eta() > kMinEtaFOCAL + partjets->GetJetRadius() && jet->Eta() < kMaxEtaFOCAL - partjets->GetJetRadius()) {
        fHistos->FillTH1("hJetEtaDistSel", jet->Eta());
        acceptedEvent = true;
      }
    }
    if(!acceptedEvent) {
      return false;
    }

    for(auto jet : partjets->accepted()) {
      fHistos->FillTH1("hJetEtaDistSelEvents", jet->Eta());
      if(jet->Eta() > kMinEtaFOCAL + partjets->GetJetRadius() && jet->Eta() < kMaxEtaFOCAL - partjets->GetJetRadius()) {
        if(jet->E() > elead) {
          elead = jet->E();
        }
      }
    }
  }

  std::vector<std::pair<int, int>> eranges = {{5, 10}, {10, 20}, {20, 50}, {50, 100}, {100, 150}, {150, 200}};
  auto foundrange = eranges.end();
  if(fRequireFOCALJet) {
    for(auto rangeit = eranges.begin(); rangeit != eranges.end(); ++rangeit) {
      if(elead >= rangeit->first && elead < rangeit->second) {
        foundrange = rangeit;
      }
    }
  }

  // Look at VZERO properties
  // VZERO-A: Same side as FOCAL
  // VZERO-C: Opposite site
  double sumV0a = 0, sumV0c = 0;
  std::array<double, 32> v0amults, v0cmults;
  auto vzerodata = fInputEvent->GetVZEROData();
  for(int ichan = 0; ichan < 32; ichan++) {
    auto v0amult = vzerodata->GetMultiplicityV0A(ichan), v0cmult = vzerodata->GetMultiplicityV0C(ichan);
    v0amults[ichan] = v0amult;
    v0cmults[ichan] = v0cmult;
    sumV0a += v0amult;
    sumV0c += v0cmult;
    fHistos->FillTH2("hChannelMultV0A", ichan, v0amult);
    fHistos->FillTH2("hChannelMultV0C", ichan, v0cmult);
    fHistos->FillTH2(Form("hCorrJetEMultV0A%d", ichan), ichan, v0amult);
    fHistos->FillTH2(Form("hCorrJetEMultV0C%d", ichan), ichan, v0cmult);
    if(foundrange != eranges.end()){
      fHistos->FillTH2(Form("hChannelMultV0A_%d_%d", foundrange->first, foundrange->second), ichan, v0amult);
      fHistos->FillTH2(Form("hChannelMultV0C_%d_%d", foundrange->first, foundrange->second), ichan, v0cmult);
    }
  }
  fHistos->FillTH1("hSumMultV0A", sumV0a);
  fHistos->FillTH1("hSumMultV0C", sumV0c);
  fHistos->FillTH1("hMaxMultV0A", *std::max_element(v0amults.begin(), v0amults.end()));
  fHistos->FillTH1("hMaxMultV0C", *std::max_element(v0cmults.begin(), v0cmults.end()));
  fHistos->FillTH1("hMeanMultV0A", TMath::Mean(v0amults.begin(), v0amults.end()));
  fHistos->FillTH1("hMeanMultV0C", TMath::Mean(v0cmults.begin(), v0cmults.end()));
  if(fRequireFOCALJet) {
    fHistos->FillTH2("hCorrJetESumMultV0A", elead, sumV0a);
    fHistos->FillTH2("hCorrJetESumMultV0C", elead, sumV0c);
    fHistos->FillTH2("hCorrJetEMaxMultV0A", elead, *std::max_element(v0amults.begin(), v0amults.end()));
    fHistos->FillTH2("hCorrJetEMaxMultV0C", elead, *std::max_element(v0cmults.begin(), v0cmults.end()));
    fHistos->FillTH2("hCorrJetEMeanMultV0A", elead, TMath::Mean(v0amults.begin(), v0amults.end()));
    fHistos->FillTH2("hCorrJetEMeanMultV0C", elead, TMath::Mean(v0cmults.begin(), v0cmults.end()));
  }

  // Look at cluster Spectrum in EMCAL
  double totalEnergyEMCAL = 0, totalEnergyDCAL = 0, leadingEMCAL = 0, leadingDCAL = 0, multiplicityEMCAL = 0, multiplicityDCAL = 0;
  auto clusters = this->GetClusterContainer(0);
  for(auto clust : clusters->accepted()) {
    auto clusterE = clust->GetNonLinCorrEnergy();
    TLorentzVector clustervec;
    clust->GetMomentum(clustervec, fVertex, AliVCluster::VCluUserDefEnergy_t::kNonLinCorr);
    int absID = -1;
    fGeom->GetAbsCellIdFromEtaPhi(clustervec.Eta(), clustervec.Phi(),absID);
    int supermoduleID, moduleID, phiModule, etaModule;
    fGeom->GetCellIndex(absID, supermoduleID, moduleID, phiModule, etaModule);
    if(supermoduleID > 11) {
      multiplicityDCAL += 1.;
      if(clusterE > leadingDCAL) {
        leadingDCAL = clusterE;
      }
    } else {
      multiplicityEMCAL += 1.;
      if(clusterE > leadingEMCAL) {
        leadingEMCAL = clusterE;
      }
    }
  }
  auto cells = fInputEvent->GetEMCALCells();
  for(int icell = 0; icell < cells->GetNumberOfCells(); icell++) {
    auto towerID = cells->GetCellPosition(icell);
    auto energy = cells->GetAmplitude(icell);
    // MC - don't cut on time
    int supermoduleID, moduleID, phiModule, etaModule;
    fGeom->GetCellIndex(towerID, supermoduleID, moduleID, phiModule, etaModule);
    if(supermoduleID > 11) {
      totalEnergyDCAL += energy;
    } else {
      totalEnergyEMCAL += energy;
    }
  }
  fHistos->FillTH1("hClusterMultiplicityEMCAL", multiplicityEMCAL);
  fHistos->FillTH1("hClusterMultiplicityDCAL", multiplicityDCAL);
  fHistos->FillTH1("hClusterMultiplicityEDCAL", multiplicityEMCAL+multiplicityDCAL);
  fHistos->FillTH1("hTotalEnergyEMCAL", totalEnergyEMCAL);
  fHistos->FillTH1("hTotalEnergyDCAL", totalEnergyDCAL);
  fHistos->FillTH1("hTotalEnergyEDCAL", totalEnergyEMCAL+totalEnergyDCAL);
  fHistos->FillTH1("hLeadingClusterEMCAL", leadingEMCAL);
  fHistos->FillTH1("hLeadingClusterDCAL", leadingDCAL);
  fHistos->FillTH1("hLeadingClusterEDCAL", std::max(leadingEMCAL, leadingDCAL));

  if(fRequireFOCALJet) {
    fHistos->FillTH2("hCorrJetEClusterMultiplicityEMCAL", elead, multiplicityEMCAL);
    fHistos->FillTH2("hCorrJetEClusterMultiplicityDCAL", elead, multiplicityDCAL);
    fHistos->FillTH2("hCorrJetEClusterMultiplicityEDCAL", elead, multiplicityEMCAL+multiplicityDCAL);
    fHistos->FillTH2("hCorrJetETotalEnergyEMCAL", elead, totalEnergyEMCAL);
    fHistos->FillTH2("hCorrJetETotalEnergyDCAL", elead, totalEnergyDCAL);
    fHistos->FillTH2("hCorrJetETotalEnergyEDCAL", elead, totalEnergyEMCAL+totalEnergyDCAL);
    fHistos->FillTH2("hCorrJetELeadingClusterEMCAL", elead, leadingEMCAL);
    fHistos->FillTH2("hCorrJetELeadingClusterDCAL", elead, leadingDCAL);
    fHistos->FillTH2("hCorrJetELeadingClusterEDCAL", elead, std::max(leadingEMCAL, leadingDCAL));
  }

  return true;
}

Bool_t AliAnalysisTaskForwardJets::CheckMCOutliers(){
  auto outlierjets = this->GetJetContainer(fNamePartJetCont.Data());
  // Check whether there is at least one particle level jet with pt above n * event pt-hard
  auto jetiter = outlierjets->accepted();
  auto max = std::max_element(jetiter.begin(), jetiter.end(), [](const AliEmcalJet *lhs, const AliEmcalJet *rhs ) { return lhs->Pt() < rhs->Pt(); });
  if(max != jetiter.end())  {
    // At least one jet found with pt > n * pt-hard
    AliDebugStream(2) << "Found max jet with pt " << (*max)->Pt() << " GeV/c" << std::endl;
    if(fDebugMaxJetOutliers) {
      // cross check whether implemenation using stl gives the same result as a trivial manual iteration
      // over all jets
      // for debug purposes
      auto jetiter1 = outlierjets->accepted();
      AliEmcalJet *debugmax(nullptr);
      for(auto testjet : jetiter1) {
        if(!debugmax) debugmax = testjet;
        else {
          if(testjet->Pt() > debugmax->Pt()) debugmax = testjet;
        }
      }
      fHistos->FillTH2("hDebugMaxJetPt", (*max)->Pt(), debugmax->Pt());
      if(*max != debugmax) fHistos->FillTH1("hDebugMaxJetError", 1.);
    }
    if((*max)->Pt() > fPtHardAndJetPtFactor * fPtHard) return false;
  }
  return true;
}

AliAnalysisTaskForwardJets *AliAnalysisTaskForwardJets::AddTaskForwardJets(const char *suffix){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale", "No analysis manager available");
    return nullptr;
  }

  auto inputhandler = mgr->GetInputEventHandler();
  auto isAOD = inputhandler->IsA() == AliAODInputHandler::Class();

  std::stringstream tasknamebuilder, outnamebuilder, listnamebuilder;
  tasknamebuilder << "FOCALJetTask";
  listnamebuilder << "FOCALJetHists";
  outnamebuilder << mgr->GetCommonFileName() << ":FOCALJetResults";
  if(strlen(suffix)) {
    tasknamebuilder << "_" << suffix;
    listnamebuilder << "_" << suffix;
    outnamebuilder << "_" << suffix;
  }

  auto task = new AliAnalysisTaskForwardJets(tasknamebuilder.str().data());
  mgr->AddTask(task);

  auto partcont = task->AddMCParticleContainer("mcparticles");
  partcont->SetMinPt(0.);
  const std::string kNameMCParticles = "MCParticles";
  partcont->SetName(kNameMCParticles.data());

  const std::string kNameJetsPart = "particleLevelJets";
  auto contpartjet = task->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, 0.4,
                                                      AliJetContainer::kUser, partcont, nullptr);
  contpartjet->SetName(kNameJetsPart.data());
  task->setNameParticleJetContainer(kNameJetsPart.data());
  std::cout << "Adding particle-level jet container with underling array: " << contpartjet->GetArrayName() << std::endl;

  const std::string kNameClusterContainer = "EMCALClusters";
  auto clusters = task->AddClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
  clusters->SetDefaultClusterEnergy(AliVCluster::VCluUserDefEnergy_t::kNonLinCorr);
  clusters->SetClusUserDefEnergyCut(AliVCluster::VCluUserDefEnergy_t::kNonLinCorr, 0.3);
  clusters->SetName(kNameClusterContainer.data());


  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(listnamebuilder.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, outnamebuilder.str().data()));

  return task;
}
