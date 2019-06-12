/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
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

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalTriggerNormalization.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEmcalTriggerStringDecoder.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliVEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerNormalization)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalTriggerNormalization::AliAnalysisTaskEmcalTriggerNormalization():
    AliAnalysisTaskEmcal(),
    fHistos(nullptr),
    fTriggerCluster(),
    fMBTriggerClasses()
{
}

AliAnalysisTaskEmcalTriggerNormalization::AliAnalysisTaskEmcalTriggerNormalization(const char *name):
    AliAnalysisTaskEmcal(name, kTRUE),
    fHistos(nullptr),
    fTriggerCluster(),
    fMBTriggerClasses()
{
}

void AliAnalysisTaskEmcalTriggerNormalization::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fHistos = new THistManager("HistosTriggerNorm");

  fHistos->CreateTH2("hTriggerNorm", "Histogram for the trigger normalization", 11, -0.5, 10.5, 100, 0., 100.);
  auto normhist = static_cast<TH1 *>(fHistos->GetListOfHistograms()->FindObject("hTriggerNorm"));
  std::array<std::string, 9> triggers = {"INT7", "EG", "EG2", "EJ1", "EJ2", "DG1", "DG2", "DJ1", "DJ2"};
  for(auto ib = 0; ib < triggers.size(); ib++){
    normhist->GetXaxis()->SetBinLabel(ib+1, triggers[ib].data());
  }

  for(auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalTriggerNormalization::Run(){
  if(!fTriggerCluster.length()) throw TriggerClusterNotSetException(); 
  if(!fMBTriggerClasses.size()) throw MBTriggerNotSetException();

  double centralitypercentile = 99.;
  if(this->GetBeamType() == AliAnalysisTaskEmcal::kAA) {
    AliMultSelection *mult = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
    if(mult){
      centralitypercentile = mult->GetMultiplicityPercentile(fCentEst.Data());
    } else {
      throw CentralityNotSetException();
    }
  }

  EMCAL_STRINGVIEW triggerstring(fInputEvent->GetFiredTriggerClasses().Data());

  // Min. bias trigger (reference trigger)
  if(fInputHandler->IsEventSelected() & AliVEvent::kINT7) {
    auto match_triggerclass = MatchTrigger(triggerstring, fMBTriggerClasses, fTriggerCluster);
    if(match_triggerclass.length()){
      double weight = 1./PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(match_triggerclass.data());
      auto normhist = static_cast<TH2 *>(fHistos->FindObject("hTriggerNorm"));
      auto triggerbin = normhist->GetXaxis()->FindBin("INT7");
      normhist->Fill(normhist->GetXaxis()->GetBinCenter(triggerbin), centralitypercentile, weight);
    }
  }

  // EMCAL triggers
  const UInt_t EGABIT = AliVEvent::kEMCEGA,
               EJEBIT = AliVEvent::kEMCEJE;
  const Int_t NEMCAL_TRIGGERS = 8;
  const std::array<std::string, NEMCAL_TRIGGERS> EMCAL_TRIGGERS = {{"EG1", "EG2", "EJ1", "EJ2", "DG1", "DG2", "DJ1", "DJ2"}};
  const std::array<UInt_t, NEMCAL_TRIGGERS> EMCAL_TRIGGERBITS = {{EGABIT, EGABIT, EJEBIT, EJEBIT, EGABIT, EGABIT, EJEBIT, EJEBIT}};
  for(Int_t i = 0; i < NEMCAL_TRIGGERS; i++) {
    std::vector<std::string> match_emctriggers = {EMCAL_TRIGGERS[i]};
    if(fInputHandler->IsEventSelected() & EMCAL_TRIGGERBITS[i]){
      auto match_triggerclass = MatchTrigger(triggerstring, match_emctriggers, fTriggerCluster);
      if(!match_triggerclass.length()) continue;      // match trigger class not found
      double weight = 1./PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(match_triggerclass.data());
      auto normhist = static_cast<TH2 *>(fHistos->FindObject("hTriggerNorm"));
      auto triggerbin = normhist->GetXaxis()->FindBin(EMCAL_TRIGGERS[i].data());
      normhist->Fill(normhist->GetXaxis()->GetBinCenter(triggerbin), centralitypercentile, weight);
    }
  }
  return true;
}

std::string AliAnalysisTaskEmcalTriggerNormalization::MatchTrigger(EMCAL_STRINGVIEW triggerstring, const std::vector<std::string> &triggerclasses, EMCAL_STRINGVIEW triggercluster) const {
  std::string result;
  auto triggers = PWG::EMCAL::Triggerinfo::DecodeTriggerString(triggerstring);
  for(const auto &t : triggers) {
    if(t.Triggercluster() != triggercluster) continue;
    for(const auto &c : triggerclasses) {
      if(t.IsTriggerClass(c)) {
        result = t.ExpandClassName();
        break;
      }
    }
  }
  return result;
}

void AliAnalysisTaskEmcalTriggerNormalization::RunChanged(int newrun){
  PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->SetRun(newrun);
}

AliAnalysisTaskEmcalTriggerNormalization *AliAnalysisTaskEmcalTriggerNormalization::AddTaskTriggerNormalization(const char *name) {
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    AliErrorGeneralStream("AliAnalysisTaskEmcalTriggerNormalization::AddTaskTriggerNormalization") << "No analysis manager defined" << std::endl;
    return nullptr;
  }

  auto task = new AliAnalysisTaskEmcalTriggerNormalization(name);
  mgr->AddTask(task);

  std::string outputfile = mgr->GetCommonFileName();
  outputfile += Form(":EmcalTriggerNorm%s", name);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("EmcalNormalizationHistograms", TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.data()));

  return task;
}