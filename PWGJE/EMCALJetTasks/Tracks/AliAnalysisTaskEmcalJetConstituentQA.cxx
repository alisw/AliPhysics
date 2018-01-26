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
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <THistManager.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetConstituentQA.h"
#include "AliAODInputHandler.h"
#include "AliClusterContainer.h"
#include "AliEmcalAnalysisFactory.h"
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
  fTriggerSelectionString("")
{
  this->SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalJetConstituentQA::AliAnalysisTaskEmcalJetConstituentQA(const char *name):
  AliAnalysisTaskEmcalJet(name, true),
  fHistos(nullptr),
  fNameTrackContainer(""),
  fNameClusterContainer(""),
  fTriggerSelectionString("")
{
  this->SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalJetConstituentQA::~AliAnalysisTaskEmcalJetConstituentQA(){
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalJetConstituentQA::UserCreateOutputObjects(){
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fHistos = new THistManager(Form("histos_%s", GetName()));
  for(auto c : fNamesJetContainers){
    auto contname = dynamic_cast<TObjString *>(c);
    if(!contname) continue;
    fHistos->CreateTH2(Form("hPtChargedConstituents%s", contname->String().Data()), "Charged jet constituents; p_{t, jet} (GeV/c); p_{t,ch}", 200, 0., 200., 200., 0., 200.);
    fHistos->CreateTH2(Form("hPtNeutralConstituents%s", contname->String().Data()), "Neutral jet constituents; p_{t, jet} (GeV/c); p_{t,ne}", 200, 0., 200., 200., 0., 200.);
    fHistos->CreateTH2(Form("hENeutralConstituents%s", contname->String().Data()), "Neutral jet constituents; p_{t, jet} (GeV/c); E_{ne}", 200, 0., 200., 200., 0., 200.);
    fHistos->CreateTH2(Form("hENeutralConstituentsNoSub%s", contname->String().Data()), "Neutral jet constituent without hadronic correction; p_{t, jet} (GeV/c); E_{ne}", 200, 0., 200., 200., 0., 200.);
    fHistos->CreateTH2(Form("hZChargedConstituents%s", contname->String().Data()), "Charged jet constituents; p_{t, jet} (GeV/c); z", 200, 0., 200., 100., 0., 1.);
    fHistos->CreateTH2(Form("hZNeutralConstituents%s", contname->String().Data()), "Neutral jet constituent without hadronic correction; p_{t, jet} (GeV/c); z", 200, 0., 200., 100., 0., 1.);
  }

  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);  
  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalJetConstituentQA::Run(){
  auto tracks = GetTrackContainer(fNameTrackContainer);
  const auto clusters = GetClusterContainer(fNameClusterContainer);
  if(!(tracks && clusters)){
    AliErrorStream() << "Either track or cluster container missing, aborting ..." << std::endl;
  }

  // Event selection
  if(fTriggerSelectionString.Contains("INT7")){
    // INT7 trigger
    if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;
  } else if(fTriggerSelectionString.Contains("EJ")){
     // EMCAL JET trigger
     if(!(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE)) return false;
     if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerSelectionString)) return false;
  } else return false;

  for(auto jc : fNamesJetContainers){
    auto contname = dynamic_cast<TObjString *>(jc);
    if(!contname) continue;
    const auto jetcont = GetJetContainer(contname->String().Data());
    if(!jetcont) continue;

    for(auto jet : jetcont->accepted()){
      for(decltype(jet->GetNumberOfTracks()) itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++){
        const auto trk = jet->TrackAt(itrk, tracks->GetArray());
        if(!trk) continue;
        fHistos->FillTH2(Form("hPtChargedConstituents%s", contname->String().Data()), std::abs(jet->Pt()), std::abs(trk->Pt()));
        fHistos->FillTH2(Form("hZChargedConstituents%s", contname->String().Data()), std::abs(jet->Pt()), jet->GetZ(trk));
      }
      for(decltype(jet->GetNumberOfClusters()) icl = 0; icl < jet->GetNumberOfClusters(); icl++){
        const auto clust = jet->ClusterAt(icl, clusters->GetArray());
        if(!clust) continue; 
        TLorentzVector ptvec;
        clust->GetMomentum(ptvec, this->fVertex, AliVCluster::kHadCorr);
        fHistos->FillTH2(Form("hPtNeutralConstituents%s", contname->String().Data()), std::abs(jet->Pt()), std::abs(ptvec.Pt()));
        fHistos->FillTH2(Form("hENeutralConstituents%s",  contname->String().Data()), std::abs(jet->Pt()), std::abs(clust->GetHadCorrEnergy()));
        fHistos->FillTH2(Form("hENeutralConstituentsNoSub%s",  contname->String().Data()), std::abs(jet->Pt()), std::abs(clust->GetNonLinCorrEnergy()));
        fHistos->FillTH2(Form("hZNeutralConstituents%s", contname->String().Data()), std::abs(jet->Pt()), jet->GetZ(ptvec.Px(), ptvec.Py(), ptvec.Pz()));
      }
    }
  }
  
  return kTRUE;
}

AliAnalysisTaskEmcalJetConstituentQA *AliAnalysisTaskEmcalJetConstituentQA::AddTaskEmcalJetConstituentQA(const char *trigger){
  using AnalysisHelpers = EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    std::cerr << "[AliAnalysisTaskJetConstituentQA::AddTaskEmcalJetConstituentQA(EE)] No analysis manager provided ..." << std::endl;
    return nullptr;
  }
  
  std::stringstream taskname;
  taskname << "constituentQA_" << trigger;
  auto task = new AliAnalysisTaskEmcalJetConstituentQA(taskname.str().data());
  mgr->AddTask(task);

  auto inputhandler = mgr->GetInputEventHandler();
  auto isAOD = false;
  if(inputhandler->IsA() == AliAODInputHandler::Class()) isAOD = true;

  auto tracksname = AnalysisHelpers::TrackContainerNameFactory(isAOD);
  auto tracks = task->AddTrackContainer(tracksname);
  task->SetNameTrackContainer(tracksname);
  tracks->SetMinPt(0.15);

  auto clustername = AnalysisHelpers::ClusterContainerNameFactory(isAOD);
  auto clusters = task->AddClusterContainer(clustername);
  task->SetNameClusterContainer(clustername);
  clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  clusters->SetClusHadCorrEnergyCut(0.3);

  // create jet containers for R02 and R04 jets
  std::array<double, 2> jetradii = {{0.2, 0.4}};
  for(auto r : jetradii) {
    std::stringstream contname;
    contname << "fulljets_R" << std::setw(2) << std::setfill('0') << r;
    auto jcont = task->AddJetContainer(AliJetContainer::kFullJet, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, r, AliJetContainer::kEMCALfid, tracks, clusters, "Jet");
    jcont->SetName(contname.str().data());
    task->AddNameJetContainer(contname.str().data());
    jcont->SetMinPt(20.);
  }

  std::stringstream contname, outfilename;
  contname << "JetConstituentQA_" << trigger;
  outfilename << mgr->GetCommonInputContainer() << ":JetConstituentQA_" << trigger;
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname.str().data(), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));

  return task;
}