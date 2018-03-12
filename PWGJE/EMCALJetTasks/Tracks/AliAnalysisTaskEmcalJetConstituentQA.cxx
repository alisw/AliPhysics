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
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetConstituentQA.h"
#include "AliAODInputHandler.h"
#include "AliClusterContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliTrackContainer.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVTrack.h"

/// \cond CLASSIMP
ClassImp(EmcalTriggerJets::AliAnalysisTaskEmcalJetConstituentQA);
/// \endcond

using namespace EmcalTriggerJets;

AliAnalysisTaskEmcalJetConstituentQA::AliAnalysisTaskEmcalJetConstituentQA():
  AliAnalysisTaskEmcalJet(),
  fHistos(nullptr),
  fNameTrackContainer(""),
  fNameClusterContainer(""),
  fTriggerSelectionString(""),
  fUseTriggerSelection(kFALSE),
  fNameTriggerDecisionContainer("EmcalTriggerDecision")
{
  this->SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalJetConstituentQA::AliAnalysisTaskEmcalJetConstituentQA(const char *name):
  AliAnalysisTaskEmcalJet(name, true),
  fHistos(nullptr),
  fNameTrackContainer(""),
  fNameClusterContainer(""),
  fTriggerSelectionString(""),
  fUseTriggerSelection(kFALSE),
  fNameTriggerDecisionContainer("EmcalTriggerDecision")
{
  this->SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalJetConstituentQA::~AliAnalysisTaskEmcalJetConstituentQA(){
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalJetConstituentQA::UserCreateOutputObjects(){
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  TLinearBinning binningz(50, 0., 1), multbinning(51, -0.5, 50.5), binningnef(50, 0., 1.), binningR(50, 0., 0.5), binningptconst(200, 0., 200.), binningptjet(20, 0., 200.),
                 binningNCell(101, -0.5, 100.5), binningFracCellLeading(100, 0., 1.), binningM02(100, 0., 1.);

  const TBinning *jetbinning[4] = {&binningptjet, &binningnef, &multbinning, &multbinning},
                 *chargedbinning[7] = {&binningptjet, &binningnef, &multbinning, &multbinning, &binningptconst, &binningz, &binningR},
                 *neutralbinning[8] = {&binningptjet, &binningnef, &multbinning, &multbinning, &binningptconst, &binningptconst, &binningz, &binningR},
                 *binningHighZClusters[7] = {&binningptjet, &binningnef, &binningptconst, &binningz, &binningNCell, &binningFracCellLeading, &binningM02};

  fHistos = new THistManager(Form("histos_%s", GetName()));
  for(auto c : fNamesJetContainers){
    auto contname = dynamic_cast<TObjString *>(c);
    if(!contname) continue;
    fHistos->CreateTHnSparse(Form("hJetCounter%s", contname->String().Data()), Form("jet counter for jets %s", contname->String().Data()), 4, jetbinning);
    fHistos->CreateTHnSparse(Form("hChargedConstituents%s", contname->String().Data()), Form("charged constituents in jets %s", contname->String().Data()), 7, chargedbinning);
    fHistos->CreateTHnSparse(Form("hNeutralConstituents%s", contname->String().Data()), Form("neutral constituents in jets %s", contname->String().Data()), 8, neutralbinning);
    fHistos->CreateTHnSparse(Form("hHighZClusters"), "Properties of high-z clusters", 7, binningHighZClusters);
  }

  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);  
  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalJetConstituentQA::Run(){
  AliParticleContainer * tracks = GetTrackContainer(fNameTrackContainer);
  if(!tracks) tracks = GetParticleContainer(fNameTrackContainer);
  const auto clusters = GetClusterContainer(fNameClusterContainer);
  if(fNameTrackContainer.Length() && !tracks){
    AliErrorStream() << "Track container " << fNameTrackContainer << " required but missing ..." << std::endl;
    return kFALSE;
  }
  if(fNameClusterContainer.Length() && !clusters){
    AliErrorStream() << "Cluster container " << fNameClusterContainer << " required but missing ..." << std::endl;
    return kFALSE;
  }

  // Event selection
  AliDebugStream(1) << "Trigger selection string: " << fTriggerSelectionString << ", fired trigger classes: " << fInputEvent->GetFiredTriggerClasses() << std::endl;
  if(fTriggerSelectionString.Contains("INT7")){
    // INT7 trigger
    if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;
  } else if(fTriggerSelectionString.Contains("EJ")){
    auto triggerclass = fTriggerSelectionString(fTriggerSelectionString.Index("EJ"),3); // cleanup trigger string from possible tags
    AliDebugStream(1) << "Inspecting trigger class " << triggerclass << std::endl;
    // EMCAL JET trigger
    if(!fMCEvent){ // Running on Data
      if(!(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE)) return false;
      if(!fInputEvent->GetFiredTriggerClasses().Contains(triggerclass)) return false;
    }

    if(fUseTriggerSelection) {
      auto triggerdecisions = dynamic_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer.Data()));
      if(!triggerdecisions) {
        AliErrorStream() << "No offline trigger selection available" << std::endl;
        return false;
      }
      else if(!triggerdecisions->IsEventSelected(triggerclass.Data())) return false;
    }
  } else return false;

  AliDebugStream(1) << "Event is selected" << std::endl;

  for(auto jc : fNamesJetContainers){
    auto contname = dynamic_cast<TObjString *>(jc);
    if(!contname) {
      AliErrorStream() << "Non-string object in the list of jet container names" << std::endl;
      continue;
    } 
    const auto jetcont = GetJetContainer(contname->String().Data());
    if(!jetcont){
      AliErrorStream() << "Jet container with name " << contname->String() << " not found in the list of jet containers" << std::endl;
      continue;
    } 
    AliDebugStream(2) << "Reading " << jetcont->GetArray()->GetName() << std::endl;

    for(auto jet : jetcont->accepted()){
      AliDebugStream(3) << "Next accepted jet, found " << jet->GetNumberOfTracks() << " tracks and " << jet->GetNumberOfClusters() << " clusters." << std::endl;
      Double_t pointjet[4] = {std::abs(jet->Pt()), jet->NEF(), static_cast<double>(jet->GetNumberOfTracks()), static_cast<double>(jet->GetNumberOfClusters())}, 
               pointcharged[7] = {std::abs(jet->Pt()), jet->NEF(), static_cast<double>(jet->GetNumberOfTracks()), static_cast<double>(jet->GetNumberOfClusters()), -1., 1., -1.}, 
               pointneutral[8] = {std::abs(jet->Pt()), jet->NEF(), static_cast<double>(jet->GetNumberOfTracks()), static_cast<double>(jet->GetNumberOfClusters()), -1., 1., -1., -1.},
               pointHighZCluster[7] = {std::abs(jet->Pt()), jet->NEF(), -1., -1., -1., -1., -1.};
      fHistos->FillTHnSparse(Form("hJetCounter%s", contname->String().Data()), pointjet);
      TVector3 jetvec{jet->Px(), jet->Py(), jet->Pz()};
      if(tracks){
        for(decltype(jet->GetNumberOfTracks()) itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++){
          const auto trk = jet->TrackAt(itrk, tracks->GetArray());
          if(!trk) continue;
          if(trk->Charge()){
            pointcharged[4] = std::abs(trk->Pt());
            pointcharged[5] = std::abs(jet->GetZ(trk));
            pointcharged[6] = jet->DeltaR(trk);
            fHistos->FillTHnSparse(Form("hChargedConstituents%s", contname->String().Data()), pointcharged);
          } else {
            // particle level jets
            pointneutral[4] = pointneutral[5] = std::abs(trk->E());
            pointneutral[6] = std::abs(jet->GetZ(trk));
            pointneutral[7] = jet->DeltaR(trk);
            fHistos->FillTHnSparse(Form("hNeutralConstituents%s", contname->String().Data()), pointneutral);
          }
        }
      }
      if(clusters){
        for(decltype(jet->GetNumberOfClusters()) icl = 0; icl < jet->GetNumberOfClusters(); icl++){
          const auto clust = jet->ClusterAt(icl, clusters->GetArray());
          if(!clust) continue; 
          TLorentzVector ptvec;
          clust->GetMomentum(ptvec, this->fVertex, AliVCluster::kHadCorr);
          pointneutral[4] = std::abs(clust->GetHadCorrEnergy());
          pointneutral[5] = std::abs(clust->GetNonLinCorrEnergy());
          pointneutral[6] = jet->GetZ(ptvec.Px(), ptvec.Py(), ptvec.Pz());
          pointneutral[7] = jetvec.DeltaR(ptvec.Vect());
          fHistos->FillTHnSparse(Form("hNeutralConstituents%s", contname->String().Data()), pointneutral);

          if(pointneutral[6] > 0.95) {
            pointHighZCluster[2] = pointneutral[4];
            pointHighZCluster[3] = pointneutral[6];
            pointHighZCluster[4] = clust->GetNCells();
            pointHighZCluster[5] = *std::max_element(clust->GetCellsAmplitudeFraction(), clust->GetCellsAmplitudeFraction()+clust->GetNCells());
            pointHighZCluster[6] = clust->GetM02();
            fHistos->FillTHnSparse("hHighZClusters", pointHighZCluster);
          }
        }
      }
    }
  }
  
  return kTRUE;
}

AliAnalysisTaskEmcalJetConstituentQA *AliAnalysisTaskEmcalJetConstituentQA::AddTaskEmcalJetConstituentQA(const char *trigger, bool partmode){
  using AnalysisHelpers = EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    std::cerr << "[AliAnalysisTaskJetConstituentQA::AddTaskEmcalJetConstituentQA(EE)] No analysis manager provided ..." << std::endl;
    return nullptr;
  }
  
  std::stringstream taskname;
  taskname << "constituentQA_" << trigger;
  auto task = new AliAnalysisTaskEmcalJetConstituentQA(taskname.str().data());
  task->SetTriggerSelection(trigger);
  mgr->AddTask(task);

  auto inputhandler = mgr->GetInputEventHandler();
  auto isAOD = false;
  if(inputhandler->IsA() == AliAODInputHandler::Class()) isAOD = true;

  TString tracksname, clustername;
  AliParticleContainer *tracks(nullptr);
  AliClusterContainer *clusters(nullptr);
  if(!partmode) {
    tracksname = AnalysisHelpers::TrackContainerNameFactory(isAOD);
    tracks = task->AddTrackContainer(tracksname);
    task->SetNameTrackContainer(tracksname);
    tracks->SetMinPt(0.15);

    clustername = AnalysisHelpers::ClusterContainerNameFactory(isAOD);
    clusters = task->AddClusterContainer(clustername);
    task->SetNameClusterContainer(clustername);
    clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    clusters->SetClusHadCorrEnergyCut(0.3);
  } else {
    tracksname = "mcparticles";
    tracks = task->AddParticleContainer(tracksname);
    task->SetNameTrackContainer(tracksname);
    tracks->SetMinPt(0.);
  }

  // create jet containers for R02 and R04 jets
  std::array<double, 2> jetradii = {{0.2, 0.4}};
  for(auto r : jetradii) {
    std::stringstream contname;
    contname << "fulljets_R" << std::setw(2) << std::setfill('0') << int(r*10.);
    auto jcont = task->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, r, AliJetContainer::kEMCALfid, tracks, clusters, "Jet");
    jcont->SetName(contname.str().data());
    task->AddNameJetContainer(contname.str().data());
    jcont->SetMinPt(20.);
  }

  std::stringstream contname, outfilename;
  contname << "JetConstituentQA_" << trigger;
  outfilename << mgr->GetCommonFileName() << ":JetConstituentQA_" << trigger;
  if(partmode) {
    contname << "_part";
    outfilename << "_part";
  }
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname.str().data(), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));

  return task;
}
