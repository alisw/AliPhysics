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
#include <string>
#include <vector>

#include <THistManager.h>
#include <TLinearBinning.h>
#include <TVariableBinning.h>

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetPbPbTriggerTests.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliVEvent.h"

ClassImp(EmcalTriggerJets::AliAnalysisTaskEmcalJetPbPbTriggerTests);

using namespace EmcalTriggerJets;

AliAnalysisTaskEmcalJetPbPbTriggerTests::AliAnalysisTaskEmcalJetPbPbTriggerTests():
  AliAnalysisTaskEmcalJet(),
  fHistos(nullptr),
  fIsMC(false),
  fNameJetContainer("datajets"),
  fTriggerSelectionBits(AliVEvent::kAny),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0)
{
  SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalJetPbPbTriggerTests::AliAnalysisTaskEmcalJetPbPbTriggerTests(const char *name):
  AliAnalysisTaskEmcalJet(name, true),
  fHistos(nullptr),
  fIsMC(false),
  fNameJetContainer("datajets"),
  fTriggerSelectionBits(AliVEvent::kAny),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0)
{
  SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalJetPbPbTriggerTests::~AliAnalysisTaskEmcalJetPbPbTriggerTests(){
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalJetPbPbTriggerTests::UserCreateOutputObjects(){
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  if (fUseAliEventCuts) fEventCuts.OverrideAutomaticTriggerSelection(fOffTrigger);

  if(!fUserPtBinning.GetSize()) {
    // binning not set. apply default binning
    AliInfoStream() << "Using default pt binning";
    fUserPtBinning.Set(201);
    double current(0.);
    for(int istep = 0; istep < 201; istep++) {
      fUserPtBinning[istep] = current;
      current += 1; 
    }
  }

  TLinearBinning etabinning(100, -1., 1.), phibinning(100., 0., 7.), nefbinning(100, 0., 1.), centralitybinning(10, 0, 100);
  TVariableBinning jetptbinning(fUserPtBinning);
  const TBinning *binnings[5] = {&jetptbinning, &etabinning, &phibinning, &nefbinning, &centralitybinning};
  fHistos = new THistManager(Form("Histos_%s", GetName()));
  fHistos->CreateTH1("hEventCounter", "Event counter histogram", 1, 0.5, 1.5);
  fHistos->CreateTHnSparse("hJetTHnSparse", "jet thnsparse", 6, binnings, "s");
  fHistos->CreateTHnSparse("hMaxJetTHnSparse", "jet thnsparse", 6, binnings, "s");

  for(auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalJetPbPbTriggerTests::Run(){
  auto datajets = this->GetJetContainer(fNameJetContainer);
  if(!datajets) {
    AliErrorStream() << "Jet container " << fNameJetContainer << " not found" << std::endl;
    return false;
  }

  fHistos->FillTH1("hEventCounter", 1);
  AliEmcalJet *maxjet(nullptr);
  for(auto j : datajets->accepted()){
    if(!maxjet || (j->E() > maxjet->E())) maxjet = j;
    double datapoint[5] = {j->Pt(), j->Eta(), j->Phi(), j->NEF(), 0.};
      fHistos->FillTHnSparse("hJetTHnSparse", datapoint);
  }

  double maxdata[6];
  memset(maxdata, 0., sizeof(double) * 6);
  if(maxjet){
    maxdata[0] = maxjet->Pt();
    maxdata[1] = maxjet->Eta();
    maxdata[2] = maxjet->Phi();
    maxdata[3] = maxjet->NEF();
    maxdata[4] = fCent;
      fHistos->FillTHnSparse("hMaxJetTHnSparse", maxdata);
  }
  return true;
}

bool AliAnalysisTaskEmcalJetPbPbTriggerTests::IsTriggerSelected() {
    if(!(fInputHandler->IsEventSelected() & fTriggerSelectionBits)) return false;
    if (fUseAliEventCuts && !fEventCuts.AcceptEvent(InputEvent())) return false;
  return true;
}

AliAnalysisTaskEmcalJetPbPbTriggerTests *AliAnalysisTaskEmcalJetPbPbTriggerTests::AddTaskJetPbPbTriggerTests(Bool_t isMC, AliJetContainer::EJetType_t jettype, double radius, const char *suffix){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    std::cerr << "Analysis manager not initialized" << std::endl;
    return nullptr;
  }

  Bool_t isAOD(kFALSE);
  AliInputEventHandler *inputhandler = static_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
  if(inputhandler) {
    if(inputhandler->IsA() == AliAODInputHandler::Class()){
      std::cout << "Analysing AOD events\n";
      isAOD = kTRUE;
    } else {
      std::cout << "Analysing ESD events\n";
    }
  }

  std::string jettypestring;
  UInt_t acctype(AliJetContainer::kTPCfid);
  switch(jettype){
    case AliJetContainer::kChargedJet:  jettypestring = "ChargedJets"; acctype = AliJetContainer::kTPCfid; break;
    case AliJetContainer::kFullJet:     jettypestring = "FullJets";    acctype = AliJetContainer::kEMCALfid; break;
    case AliJetContainer::kNeutralJet:  jettypestring = "NeutralJets"; acctype = AliJetContainer::kEMCALfid; break;
  };

  std::stringstream tag, outfilename;
  tag << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(radius * 10.) << "_";
  if(strlen(suffix)) {
    tag << "_" << suffix;
  }
  auto task = new AliAnalysisTaskEmcalJetPbPbTriggerTests(Form("JetEnergySpectrum_%s", tag.str().data()));
  task->SetIsMC(isMC);
  mgr->AddTask(task);

  // Connect particle and cluster container
  AliTrackContainer *tracks(nullptr);
  AliClusterContainer *clusters(nullptr);
  if(jettype == AliJetContainer::kChargedJet || jettype == AliJetContainer::kFullJet) {
    tracks = task->AddTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
    tracks->SetMinPt(0.15);
  }
  if(jettype == AliJetContainer::kNeutralJet || jettype == AliJetContainer::kFullJet){
    clusters = task->AddClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    clusters->SetClusHadCorrEnergyCut(0.3);
  }


  // Create proper jet container
  auto jetcont = task->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, radius, acctype, tracks, clusters);
  jetcont->SetName("datajets");
  task->SetNameJetContainer("datajets");
  std::cout << "Adding jet container with underlying array:" << jetcont->GetArrayName() << std::endl;

  // Link input and output container
  outfilename << mgr->GetCommonFileName() << ":JetSpectrum_" << tag.str().data();
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("JetSpectrum_%s", tag.str().data()), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));

  return task;
}
