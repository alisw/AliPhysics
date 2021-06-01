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
#include <array>
#include <bitset>
#include <iostream>
#include <sstream>
#include <string>

#include <TH2.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalTriggerCorrelationMC.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliVEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerCorrelationMC)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalTriggerCorrelationMC::AliAnalysisTaskEmcalTriggerCorrelationMC():
  AliAnalysisTaskEmcal(),
  fTriggerCorrelationHist(nullptr),
  fNameTriggerDecisionContainer("EmcalTriggerDecision")
{

}

AliAnalysisTaskEmcalTriggerCorrelationMC::AliAnalysisTaskEmcalTriggerCorrelationMC(const char *name):
  AliAnalysisTaskEmcal(name, true),
  fTriggerCorrelationHist(nullptr),
  fNameTriggerDecisionContainer("EmcalTriggerDecision")
{
  SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalTriggerCorrelationMC::~AliAnalysisTaskEmcalTriggerCorrelationMC(){

}

void AliAnalysisTaskEmcalTriggerCorrelationMC::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  const std::array<const std::string, 9> kTriggers = {{"INT7", "EG1", "EG2", "EJ1", "EJ2", "DG1", "DG2", "DJ1", "DJ2"}};
  fTriggerCorrelationHist = new TH2F("hTriggerCorrelationMC", "Trigger correlation hist", kTriggers.size(), -0.5, kTriggers.size() - 0.5, kTriggers.size(), -0.5, kTriggers.size() - 0.5);

  for(decltype(kTriggers.size()) itrg = 0; itrg < kTriggers.size(); itrg++){
    fTriggerCorrelationHist->GetXaxis()->SetBinLabel(itrg+1, kTriggers[itrg].data());
    fTriggerCorrelationHist->GetYaxis()->SetBinLabel(itrg+1, kTriggers[itrg].data());
  }
  fOutput->Add(fTriggerCorrelationHist);  
  
  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalTriggerCorrelationMC::Run(){
  const std::array<const std::string, 9> kTriggers = {{"INT7", "EG1", "EG2", "EJ1", "EJ2", "DG1", "DG2", "DJ1", "DJ2"}};
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;

  auto emcaltriggers = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer.Data()));
  if(!emcaltriggers){
    AliErrorStream() << "EMCAL trigger decision container with name " << fNameTriggerDecisionContainer << " not found in the event. Exiting ..." << std::endl;
    return false;
  }

  // Get all triggers selected in the event
  std::bitset<sizeof(Short_t) * 8> triggersSelected = 0;
  for(decltype(kTriggers.size()) itrg = 0; itrg < kTriggers.size(); itrg++){
    const auto &triggername = kTriggers[itrg];
    if(triggername == "INT7") triggersSelected.set(itrg, true);
    else if(emcaltriggers->IsEventSelected(triggername.data())) triggersSelected.set(itrg, true);
  }

  for(decltype(kTriggers.size()) itrg = 0; itrg < kTriggers.size(); itrg++){
    if(!triggersSelected.test(itrg)) continue;
    const auto &xlabel = kTriggers[itrg]; 
    for(decltype(kTriggers.size()) jtrg = 0; jtrg < kTriggers.size(); jtrg++){
      if(!triggersSelected.test(jtrg)) continue;
      const auto &ylabel = kTriggers[jtrg];
      fTriggerCorrelationHist->Fill(xlabel.data(), ylabel.data(), 1.);
    }
  }
  return true;
}

AliAnalysisTaskEmcalTriggerCorrelationMC *AliAnalysisTaskEmcalTriggerCorrelationMC::AddTaskEmcalTriggerCorrelationMC(const char *suffix) {
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    std::cerr << "No Analysis Manager found. Exiting ..." << std::endl;
    return nullptr;
  }

  std::stringstream tasknamebuilder;
  tasknamebuilder << "EmcalTriggerCorrelationMC";
  if(strlen(suffix)) {
    tasknamebuilder << "_" << suffix;
  }
  auto correlationtask = new AliAnalysisTaskEmcalTriggerCorrelationMC(tasknamebuilder.str().data());
  mgr->AddTask(correlationtask);

  std::stringstream outcontname, outfilename;
  outcontname << "HistosTriggerCorrelationMC";
  outfilename << mgr->GetCommonFileName() << ":EmcalTriggerCorrelationMC";
  if(strlen(suffix)){
    outcontname << "_" << suffix;
    outfilename << "_" << suffix;
  }
  mgr->ConnectInput(correlationtask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(correlationtask, 1, mgr->CreateContainer(outcontname.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));

  return correlationtask;
}