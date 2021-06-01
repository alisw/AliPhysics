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
#include <bitset>
#include <iostream>
#include <sstream>

#include <TH2.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalEG1Correlation.h"
#include "AliEmcalTriggerStringDecoder.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalEG1Correlation)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalEG1Correlation::AliAnalysisTaskEmcalEG1Correlation():
  AliAnalysisTaskSE(),
  fOutput(nullptr),
  fCorrelationHistNoPhysSel(nullptr),
  fCorrelationHistWithPhysSel(nullptr)
{
}

AliAnalysisTaskEmcalEG1Correlation::AliAnalysisTaskEmcalEG1Correlation(const char *taskname):
  AliAnalysisTaskSE(taskname),
  fOutput(nullptr),
  fCorrelationHistNoPhysSel(nullptr),
  fCorrelationHistWithPhysSel(nullptr)
{
  DefineOutput(1, TList::Class()); 
}

AliAnalysisTaskEmcalEG1Correlation::~AliAnalysisTaskEmcalEG1Correlation(){
  if(fOutput) delete fOutput;
}

void AliAnalysisTaskEmcalEG1Correlation::UserCreateOutputObjects(){
  fOutput = new TList;
  fOutput->SetOwner(true);

  auto triggers = GetSupportedTriggers();
  fCorrelationHistNoPhysSel = new TH2D("fCorrelationHistNoPhysSel", "Trigger Correlation without physics selection", triggers.size(), -0.5, triggers.size() - 0.5, triggers.size(), -0.5, triggers.size() - 0.5);
  fCorrelationHistWithPhysSel = new TH2D("fCorrelationHistWithPhysSel", "Trigger Correlation without physics selection", triggers.size(), -0.5, triggers.size() - 0.5, triggers.size(), -0.5, triggers.size() - 0.5);

  for(decltype(triggers.size()) itrg = 0; itrg < triggers.size(); itrg++){
    fCorrelationHistNoPhysSel->GetXaxis()->SetBinLabel(itrg+1, triggers[itrg].data());
    fCorrelationHistNoPhysSel->GetYaxis()->SetBinLabel(itrg+1, triggers[itrg].data());
    fCorrelationHistWithPhysSel->GetXaxis()->SetBinLabel(itrg+1, triggers[itrg].data());
    fCorrelationHistWithPhysSel->GetYaxis()->SetBinLabel(itrg+1, triggers[itrg].data());
  }

  fOutput->Add(fCorrelationHistNoPhysSel);
  fOutput->Add(fCorrelationHistWithPhysSel);
 
  PostData(1, fOutput);
}



void AliAnalysisTaskEmcalEG1Correlation::UserExec(Option_t *){
  auto triggers = PWG::EMCAL::Triggerinfo::DecodeTriggerString(fInputEvent->GetFiredTriggerClasses().Data());
  auto supported = GetSupportedTriggers();
  std::bitset<8> firedTS(0), firedPS(0);
  for(decltype(supported.size()) itrg = 0; itrg < supported.size(); itrg++){
    if(IsTriggerSelectedTS("C" + supported[itrg], triggers)) firedTS.set(itrg);
    if(IsTriggerSelectedPS(supported[itrg], triggers)) firedPS.set(itrg);
  }

  // From trigger string
  for(decltype(supported.size()) itrg = 0; itrg < supported.size(); itrg++){
    if(!firedTS.test(itrg)) continue;
    const auto &xlabel = supported[itrg];
    for(decltype(supported.size()) jtrg = 0; jtrg < supported.size(); jtrg++){
      if(!firedTS.test(jtrg)) continue;
      const auto &ylabel = supported[jtrg];
      fCorrelationHistNoPhysSel->Fill(xlabel.data(), ylabel.data(), 1.);
    }
  }

  // From physics selection
  for(decltype(supported.size()) itrg = 0; itrg < supported.size(); itrg++){
    if(!firedPS.test(itrg)) continue;
    const auto &xlabel = supported[itrg];
    for(decltype(supported.size()) jtrg = 0; jtrg < supported.size(); jtrg++){
      if(!firedPS.test(jtrg)) continue;
      const auto &ylabel = supported[jtrg];
      fCorrelationHistWithPhysSel->Fill(xlabel.data(), ylabel.data(), 1.);
    }
  }
} 

bool AliAnalysisTaskEmcalEG1Correlation::IsTriggerSelectedTS(const std::string &trigger, const std::vector<PWG::EMCAL::Triggerinfo> &decoded) const {
  bool found(false);  
  for(auto d : decoded){
    if(d.Triggerclass() == trigger) {
      found = true;
      break;
    }
  }
  return found;
}

bool AliAnalysisTaskEmcalEG1Correlation::IsTriggerSelectedPS(const std::string &trigger, const std::vector<PWG::EMCAL::Triggerinfo> &decoded) const {
  bool isPhysSel = false;
  if(trigger.find("EG") != std::string::npos){
    isPhysSel = fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA;
  } else if(trigger.find("EJ") != std::string::npos){
    isPhysSel = fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE;
  } else if(trigger == "EMC7") {
    isPhysSel = fInputHandler->IsEventSelected() & AliVEvent::kEMC7;
  } else if(trigger == "EMC8") {
    isPhysSel = fInputHandler->IsEventSelected() & AliVEvent::kEMC8;
  }
  if(!isPhysSel) return false;
  return IsTriggerSelectedTS("C"+trigger, decoded);
}

std::vector<std::string> AliAnalysisTaskEmcalEG1Correlation::GetSupportedTriggers() const {
  return {
    "EMC7EG1", "EMC8EG1", "EMC7EG2", "EMC7EJ1", "EMC7EJ2", "EMC7", "EMC8"
  };
}

AliAnalysisTaskEmcalEG1Correlation *AliAnalysisTaskEmcalEG1Correlation::AddTaskEmcalEG1Correlation(const char *taskname){
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    std::cerr << "No analysis manager found ..." << std::endl;
    return nullptr;
  }

  auto task = new AliAnalysisTaskEmcalEG1Correlation(taskname);
  task->SelectCollisionCandidates(AliVEvent::kAny);
  mgr->AddTask(task);

  std::stringstream outfilename, outcontname;
  outfilename << mgr->GetCommonFileName() << ":EG1Correlation" << taskname;
  outcontname << "HistosEG1Correlation" << taskname;

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(outcontname.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));

  return task;
}
