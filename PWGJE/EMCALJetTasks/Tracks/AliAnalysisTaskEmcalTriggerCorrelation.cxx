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
#include <iostream>
#include <TList.h>
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalTriggerCorrelation.h"
#include "AliMultSelection.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerCorrelation)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalTriggerCorrelation::AliAnalysisTaskEmcalTriggerCorrelation():
  AliAnalysisTaskEmcalTriggerBase(),
  fRequestCentrality(false),
  fEventCentrality(99.),
  fCentralityEstimator("V0M"),
  fCentralityRange(0.,100.)
{
  SetRequireAnalysisUtils(true);
}

AliAnalysisTaskEmcalTriggerCorrelation::AliAnalysisTaskEmcalTriggerCorrelation(const char *name):
  AliAnalysisTaskEmcalTriggerBase(name),
  fRequestCentrality(false),
  fEventCentrality(99.),
  fCentralityEstimator("V0M"),
  fCentralityRange(0.,100.)
{
  SetRequireAnalysisUtils(true);
}

AliAnalysisTaskEmcalTriggerCorrelation::~AliAnalysisTaskEmcalTriggerCorrelation(){

}

bool AliAnalysisTaskEmcalTriggerCorrelation::IsUserEventSelected(){
  fEventCentrality = 99;   // without centrality put everything in the peripheral bin
  if(fRequestCentrality){
    AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if(!mult){
      AliErrorStream() << GetName() << ": Centrality selection enabled but no centrality estimator found" << std::endl;
      return false;
    }
    //if(mult->IsEventSelected()) return false;
    fEventCentrality = mult->GetEstimator(fCentralityEstimator)->GetPercentile();
    AliDebugStream(1) << GetName() << ": Centrality " <<  fEventCentrality << std::endl;
    if(!fCentralityRange.IsInRange(fEventCentrality)){
      AliDebugStream(1) << GetName() << ": reject centrality: " << fEventCentrality << std::endl;
      return false;
    } else {
      AliDebugStream(1) << GetName() << ": select centrality " << fEventCentrality << std::endl;
    }
  } else {
    AliDebugStream(1) << GetName() << ": No centrality selection applied" << std::endl;
  }

  return true;
}

AliAnalysisTaskEmcalTriggerCorrelation *AliAnalysisTaskEmcalTriggerCorrelation::AddTaskTriggerCorrelation(const char *name){
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    std::cerr << "No analysis manager available. Exiting ..." << std::endl;
    return nullptr;
  }

  auto task = new AliAnalysisTaskEmcalTriggerCorrelation(name);
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("histosTriggerCorrelation_%s", name), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TriggerCorrelationQA_%s", mgr->GetCommonFileName(), name)));
  return task;
}