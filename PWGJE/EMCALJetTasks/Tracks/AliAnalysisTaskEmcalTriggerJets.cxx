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
#include <THistManager.h>
#include <TMath.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliClusterContainer.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliTrackContainer.h"
#include "AliAnalysisTaskEmcalTriggerJets.h"
#include "AliVEvent.h"

#include <iostream>

/// \cond CLASSIMP
ClassImp(EmcalTriggerJets::AliAnalysisTaskEmcalTriggerJets)
/// \endcond

namespace EmcalTriggerJets {

AliAnalysisTaskEmcalTriggerJets::AliAnalysisTaskEmcalTriggerJets():
    AliAnalysisTaskEmcalJet(),
    fHistos(nullptr)
{

}

AliAnalysisTaskEmcalTriggerJets::AliAnalysisTaskEmcalTriggerJets(const char *name):
    AliAnalysisTaskEmcalJet(name, true),
    fHistos(nullptr)
{

}

AliAnalysisTaskEmcalTriggerJets::~AliAnalysisTaskEmcalTriggerJets() {
  delete fHistos;
}

void AliAnalysisTaskEmcalTriggerJets::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  std::vector<TString> kEmcalTriggers = {"INT7", "EJ1", "EJ2", "DJ1", "DJ2"};
  fHistos = new THistManager("EmcalJetHistos");
  for(auto t : kEmcalTriggers){
    fHistos->CreateTH1("hEventCount" + t, "Event counter for trigger " + t, 1., 0.5, 1.5);
    fHistos->CreateTH1("hPtRawFullJetR02EMCAL" + t, "Raw pt spectrum for full jets with R=0.2 in EMCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawFullJetR04EMCAL" + t, "Raw pt spectrum for full jets with R=0.4 in EMCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawChargedJetR02EMCAL" + t, "Raw pt spectrum for charged jets with R=0.2 in EMCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawChargedJetR04EMCAL" + t, "Raw pt spectrum for charged jets with R=0.4 in EMCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawNeutralJetR02EMCAL" + t, "Raw pt spectrum for neutral jets with R=0.2 in EMCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawNeutralJetR04EMCAL" + t, "Raw pt spectrum for neutral jets with R=0.4 in EMCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawFullJetR02DCAL" + t, "Raw pt spectrum for full jets with R=0.2 in DCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawFullJetR04DCAL" + t, "Raw pt spectrum for full jets with R=0.4 in DCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawChargedJetR02DCAL" + t, "Raw pt spectrum for charged jets with R=0.2 in DCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawChargedJetR04DCAL" + t, "Raw pt spectrum for charged jets with R=0.4 in DCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawNeutralJetR02DCAL" + t, "Raw pt spectrum for neutral jets with R=0.2 in DCAL for trigger " + t, 200, 0., 200);
    fHistos->CreateTH1("hPtRawNeutralJetR04DCAL" + t, "Raw pt spectrum for neutral jets with R=0.4 in DCAL for trigger " + t, 200, 0., 200);
  }
  for(auto h : *(fHistos->GetListOfHistograms())){
    fOutput->Add(h);
  }
  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalTriggerJets::Run(){
  std::vector<TString> triggers, kEmcalTriggers = {"EJ1", "EJ2", "DJ1", "DJ2"};
  if(fInputHandler->IsEventSelected() & AliVEvent::kINT7) triggers.push_back("INT7");
  if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE){
    TString fired = fInputEvent->GetFiredTriggerClasses();
    for(auto e : kEmcalTriggers){
      if(fired.Contains(e)) triggers.push_back(e);
    }
  }
  if(!triggers.size()) return false;

  for(auto t : triggers) fHistos->FillTH1("hEventCount" + t, 1);
  std::vector<TString> jettypes =  {"Full", "Charged", "Neutral"}, detectors = {"EMCAL", "DCAL"}, radii = {"R02", "R04"};
  for(auto jt : jettypes) {
    for(auto det : detectors){
      for(auto r : radii) {
        TString namejcont = jt + "Jets" + r + det,
                histnamebase = "hPtRaw" + jt + "Jet" + r + det;
        AliJetContainer *c = this->GetJetContainer(namejcont);
        if(!c) AliErrorStream() << "Not found jet container " << namejcont << std::endl;
        for(auto j : c->accepted()){
          for(auto t : triggers) {
            fHistos->FillTH1(histnamebase + t, TMath::Abs(j->Pt()));
          }
        }
      }
    }
  }
  return true;
}

AliAnalysisTaskEmcalTriggerJets *AliAnalysisTaskEmcalTriggerJets::AddTaskEmcalTriggerJets(const char *name){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisTaskEmcalTriggerJets *task = new AliAnalysisTaskEmcalTriggerJets(name);
  mgr->AddTask(task);

  // Adding containers for clusters and tracks
  AliClusterContainer *clustercont = task->AddClusterContainer("caloClusters");
  clustercont->SetMinE(0.3);
  AliTrackContainer *trackcont = task->AddTrackContainer("tracks");
  trackcont->SetMinPt(0.15);

  // Adding Jet containers
  // - Using jet radii 0.2 and 0.4
  // - Splitting between EMCAL and DCAL side
  // - processing full, charged and neutral jets

  // Full jets, R=0.2, EMCAL
  AliJetContainer *cont = task->AddJetContainer(
                              AliJetContainer::kFullJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.2,
                              AliEmcalJet::kEMCALfid,
                              trackcont, clustercont
                              );
  cont->SetName("FullJetsR02EMCAL");

  // Full jets, R=0.4, EMCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kFullJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.4,
                              AliEmcalJet::kEMCALfid,
                              trackcont, clustercont
                              );
  cont->SetName("FullJetsR04EMCAL");

  // Full jets, R=0.2, DCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kFullJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.2,
                              AliEmcalJet::kDCALfid,
                              trackcont, clustercont
                              );
  cont->SetName("FullJetsR02DCAL");

  // Full jets, R=0.4, EMCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kFullJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.4,
                              AliEmcalJet::kDCALfid,
                              trackcont, clustercont
                              );
  cont->SetName("FullJetsR04DCAL");


  // Charged jets, R=0.2, EMCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kChargedJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.2,
                              AliEmcalJet::kEMCALfid,
                              trackcont, nullptr
                              );
  cont->SetName("ChargedJetsR02EMCAL");

  // Charged jets, R=0.4, EMCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kChargedJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.4,
                              AliEmcalJet::kEMCALfid,
                              trackcont, nullptr
                              );
  cont->SetName("ChargedJetsR04EMCAL");

  // Charged jets, R=0.2, DCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kChargedJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.2,
                              AliEmcalJet::kDCALfid,
                              trackcont, nullptr
                              );
  cont->SetName("ChargedJetsR02DCAL");

  // Charged jets, R=0.4, DCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kChargedJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.4,
                              AliEmcalJet::kDCALfid,
                              trackcont, nullptr
                              );
  cont->SetName("ChargedJetsR04DCAL");


  // Neutral jets, R=0.2, EMCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kNeutralJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.2,
                              AliEmcalJet::kEMCALfid,
                              nullptr, clustercont
                              );
  cont->SetName("NeutralJetsR02EMCAL");

  // Neutral jets, R=0.4, EMCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kNeutralJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.4,
                              AliEmcalJet::kEMCALfid,
                              nullptr, clustercont
                              );
  cont->SetName("NeutralJetsR04EMCAL");

  // Neutral jets, R=0.2, DCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kNeutralJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.2,
                              AliEmcalJet::kDCALfid,
                              nullptr, clustercont
                              );
  cont->SetName("NeutralJetsR02DCAL");

  // Neutral jets, R=0.4, DCAL
  cont = task->AddJetContainer(
                              AliJetContainer::kNeutralJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.4,
                              AliEmcalJet::kDCALfid,
                              nullptr, clustercont
                              );
  cont->SetName("NeutralJetsR04DCAL");

  // Connect Input / Output containers
  TString outfilename = mgr->GetCommonFileName();
  outfilename += ":EmcalTriggerJets";
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("HistsEmcalTriggerJets", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename));

  return task;
}

} /* namespace EmcalTriggerJets */
