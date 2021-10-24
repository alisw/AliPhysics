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
#include <iostream>

#include <THistManager.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskEmcalTriggerSelectionTest.h"
#include "AliClusterContainer.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalTriggerDecision.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVEvent.h"

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::Test::AliAnalysisTaskEmcalTriggerSelectionTest)
/// \endcond

namespace PWGJE {
namespace EMCALJetTasks {
namespace Test {

AliAnalysisTaskEmcalTriggerSelectionTest::AliAnalysisTaskEmcalTriggerSelectionTest():
    AliAnalysisTaskEmcal(),
    fTriggerDecision(nullptr),
    fHistos(nullptr),
    fSelectedTriggers()
{
}

AliAnalysisTaskEmcalTriggerSelectionTest::AliAnalysisTaskEmcalTriggerSelectionTest(const char *name):
    AliAnalysisTaskEmcal(name, true),
    fTriggerDecision(nullptr),
    fHistos(nullptr),
    fSelectedTriggers()
{
  this->SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalTriggerSelectionTest::~AliAnalysisTaskEmcalTriggerSelectionTest() {
}

void AliAnalysisTaskEmcalTriggerSelectionTest::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  const std::array<TString, 9> EMCALL1TRIGGERS = {"INT7", "EG1", "EG2", "EJ1", "EJ2", "DG1", "DG2", "DJ1", "DJ2"};
  fHistos = new THistManager("triggerselectiontest");

  for(auto t : EMCALL1TRIGGERS){
    fHistos->CreateTH1(Form("hEventCount%s", t.Data()), Form("Event counter for trigger %s", t.Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexDist%s", t.Data()), Form("Vertex distribution for trigger %s", t.Data()), 160, -40., 40);
    fHistos->CreateTH1(Form("hClusterEnergyNonLinCorr%s", t.Data()), Form("Cluster energy corrected for non-linearity for trigger %s", t.Data()), 1000, 0., 100);
  }

  for(auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
}

void AliAnalysisTaskEmcalTriggerSelectionTest::UserExecOnce(){
  fTriggerDecision = dynamic_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject("EmcalTriggerDecision"));
  if(!fTriggerDecision)
    AliFatal("Trigger decision object not found in the event - task cannot run\n");
}

bool AliAnalysisTaskEmcalTriggerSelectionTest::Run(){
  // check if we have at least one trigger selected
  fSelectedTriggers.clear();
  if(fInputHandler->IsEventSelected() & AliVEvent::kINT7) fSelectedTriggers.push_back("INT7");
  for(auto t : *fTriggerDecision->GetListOfTriggerDecisions()) {
    PWG::EMCAL::AliEmcalTriggerDecision *sel = static_cast<PWG::EMCAL::AliEmcalTriggerDecision *>(t);
    AliDebugStream(2) << "Processing decision " << sel->GetName() << std::endl;
    if(sel->IsSelected()) {
      AliDebugStream(2) << "Selected" << std::endl;
      fSelectedTriggers.push_back(sel->GetName());
    }
  }
  AliDebugStream(1) << "Found " << fSelectedTriggers.size() << " triggers\n";
  return fSelectedTriggers.size() > 0;
}

bool AliAnalysisTaskEmcalTriggerSelectionTest::FillHistograms(){
  for(auto t : fSelectedTriggers) this->FillHistosForTrigger(t);
  return true;
}

void AliAnalysisTaskEmcalTriggerSelectionTest::FillHistosForTrigger(const char *trigger){
  fHistos->FillTH1(Form("hEventCount%s", trigger), 1);
  fHistos->FillTH1(Form("hVertexDist%s", trigger), fVertex[2]);
  for(auto c : this->GetClusterContainer(0)->accepted()){
    fHistos->FillTH1(Form("hClusterEnergyNonLinCorr%s", trigger), c->GetNonLinCorrEnergy());
  }
}

AliAnalysisTaskEmcalTriggerSelectionTest *AliAnalysisTaskEmcalTriggerSelectionTest::AddTaskEmcalTriggerSelectionTest(TString suffix){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    std::cerr << "AddTaskEmcalTriggerSelectionTest: Analysis manager not available" << std::endl;
    return nullptr;
  }

  TString taskname("TriggerSelectionTest"),outfilename(TString::Format("%s:TriggerSelectionTest", mgr->GetCommonFileName()));
  if(suffix.Length()){
    taskname += "_" + suffix;
    outfilename += "_" + suffix;
  }

  AliAnalysisTaskEmcalTriggerSelectionTest *testtask = new AliAnalysisTaskEmcalTriggerSelectionTest(taskname);
  mgr->AddTask(testtask);

  AliClusterContainer *clusters = testtask->AddClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class()));
  clusters->SetClusNonLinCorrEnergyCut(0.5);

  mgr->ConnectInput(testtask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(testtask, 1, mgr->CreateContainer(Form("Histos%s", taskname.Data()), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outfilename));

  return testtask;
}


} /* namespace Test */
} /* namespace EMCALJetTasks */
} /* namespace PWGJE */
