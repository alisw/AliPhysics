/************************************************************************************
 * Copyright (C) 2018, Copyright Holders of the ALICE Collaboration                 *
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
#include <sstream>
#include <string>
#include <vector>

#include <THistManager.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjString.h>

#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskEmcalClustersInJets.h"
#include "AliClusterContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"
#include "AliVEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalClustersInJets)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalClustersInJets::AliAnalysisTaskEmcalClustersInJets():
  AliAnalysisTaskEmcalJet(),
  fHistos(nullptr),
  fNamesJetContainers(),
  fNameClusterContainer(),
  fNameTriggerClass()
{
  this->SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalClustersInJets::AliAnalysisTaskEmcalClustersInJets(const char *name):
  AliAnalysisTaskEmcalJet(name, true),
  fHistos(nullptr),
  fNamesJetContainers(),
  fNameClusterContainer(),
  fNameTriggerClass()
{
  this->SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalClustersInJets::~AliAnalysisTaskEmcalClustersInJets(){  
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalClustersInJets::UserCreateOutputObjects(){
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  fHistos = new THistManager(Form("histos_%s", GetName()));

  TLinearBinning energybinning(200, 0., 100.), etabinning(100, -0.8, 0.8), phibinning(100, 0., TMath::TwoPi());
  fHistos->CreateTH1("hEventCounter", "Event counter", 1, 0.5, 1.5);
  for(auto c : fNamesJetContainers){
    auto contname = static_cast<TObjString *>(c)->String();
    fHistos->CreateTH3(Form("hClustersAttached%s", contname.Data()), Form("Clusters attached to jets %s", contname.Data()), energybinning, etabinning, phibinning);
    fHistos->CreateTH3(Form("hClustersNotAttached%s", contname.Data()), Form("Clusters attached to jets %s", contname.Data()), energybinning, etabinning, phibinning);
  }
  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);

  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalClustersInJets::Run(){
  UInt_t triggersel = AliVEvent::kAny;
  if(fNameTriggerClass == "INT7") triggersel = AliVEvent::kINT7;
  else if(fNameTriggerClass.Contains("EG") || fNameTriggerClass.Contains("DG")) triggersel = AliVEvent::kEMCEGA;
  else if(fNameTriggerClass.Contains("EJ") || fNameTriggerClass.Contains("DJ")) triggersel = AliVEvent::kEMCEJE;

  bool isEMCALAcceptance = fNameTriggerClass == "INT7" || fNameTriggerClass.Contains("E");

  if(!(fInputHandler->IsEventSelected() && triggersel)) return false;
  if(fNameTriggerClass.Length() && !fInputEvent->GetFiredTriggerClasses().Contains(fNameTriggerClass)) return false; 

  auto clusters = this->GetClusterContainer();

  for(auto c : fNamesJetContainers){
    auto contname = static_cast<TObjString *>(c)->String();
    auto jetcontainer = this->GetJetContainer(contname.Data());
    std::vector<AliVCluster *> attachedclusters;
    for(auto j : jetcontainer->accepted()){
      for(decltype(j->GetNumberOfClusters()) iclust = 0; iclust < j->GetNumberOfClusters(); iclust++){
        auto clust = j->ClusterAt(iclust, clusters->GetArray());
        attachedclusters.emplace_back(clust);
        TLorentzVector pvec;
        clust->GetMomentum(pvec, fVertex);
        fHistos->FillTH3(Form("hClustersAttached%s", contname.Data()), clust->GetNonLinCorrEnergy(), pvec.Eta(), pvec.Phi());
      }
    }

    // Find all clusters not associtated to a jet
    for(auto c : clusters->accepted()){
      if(c->GetIsExotic()) continue;
      if(std::find(attachedclusters.begin(), attachedclusters.end(), c) != attachedclusters.end()) continue;
      TLorentzVector pvec;
      c->GetMomentum(pvec, fVertex);
      bool isEMCAL = (pvec.Phi() > 0. && pvec.Phi() < 3.5);
      if((isEMCALAcceptance && !isEMCAL) || (!isEMCAL && isEMCALAcceptance)) continue;
      fHistos->FillTH3(Form("hClustersNotAttached%s", contname.Data()), c->GetNonLinCorrEnergy(), pvec.Eta(), pvec.Phi());
    }
  }

  return true;
}

void AliAnalysisTaskEmcalClustersInJets::AddNameJetContainer(const char *name){
  if(!fNamesJetContainers.FindObject(name)) fNamesJetContainers.Add(new TObjString(name));
}

AliAnalysisTaskEmcalClustersInJets *AliAnalysisTaskEmcalClustersInJets::AddTaskEmcalClustersInJets(AliJetContainer::EJetType_t jettype, const char * trigger){
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    std::cerr << "No analysis manager available. Exiting ..." << std::endl;
    return nullptr;
  }
  bool isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();

  std::string jetname;
  switch(jettype){
    case AliJetContainer::kFullJet: jetname = "fulljets"; break;
    case AliJetContainer::kNeutralJet: jetname = "neutraljets"; break;
    case AliJetContainer::kChargedJet:
    case AliJetContainer::kUndefinedJetType: break;
  };

  if(!jetname.length()){
    std::cerr << "Unsuported jet type. Exiting ..." << std::endl;
    return nullptr;
  }

  std::stringstream taskname;
  taskname << jetname << "_" << trigger;

  auto task = new AliAnalysisTaskEmcalClustersInJets(taskname.str().data());
  mgr->AddTask(task);


  // adding cluster container
  auto nameclusters = AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD);
  auto clusters = task->AddClusterContainer(nameclusters);
  task->SetNameClusterContainer(nameclusters);
  clusters->SetClusTimeCut(-20e-9, 15e-9);
  if(jettype == AliJetContainer::kNeutralJet){
    clusters->SetDefaultClusterEnergy(AliVCluster:: kNonLinCorr);
    clusters->SetClusNonLinCorrEnergyCut(0.3);
    clusters->SetClusHadCorrEnergyCut(0.);
  } else {
    clusters->SetClusHadCorrEnergyCut(0.3);
  }

  AliTrackContainer *tracks = nullptr;
  if(jettype == AliJetContainer::kFullJet){
    tracks = task->AddTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
    tracks->SetMinPt(0.15);
  }

  // Adding jet containers
  const std::array<double, 4> kJetRadii = {{0.2, 0.3, 0.4, 0.5}};
  for(auto r : kJetRadii){
    std::stringstream namejetcont;
    namejetcont << jetname << "_R" << std::setw(2) << std::setfill('0') << int(r*10.);
    auto jcont = task->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, r, AliJetContainer::kEMCALfid, tracks, clusters);
    jcont->SetName(namejetcont.str().data());
    task->AddNameJetContainer(namejetcont.str().data());
    jcont->SetMinPt(5.);
  }

  // Connect output containers
  std::stringstream outfilename;
  outfilename << mgr->GetCommonFileName() << ":ClustersInJets_" << taskname.str();
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1 , mgr->CreateContainer(Form("Histos_%s", taskname.str().data()), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));

  return task;
}
