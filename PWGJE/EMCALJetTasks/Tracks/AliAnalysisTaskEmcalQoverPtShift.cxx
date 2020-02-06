/************************************************************************************
 * Copyright (C) 2020, Copyright Holders of the ALICE Collaboration                 *
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
#include <string>
#include <sstream>
#include "THashList.h"
#include "THistManager.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalQoverPtShift.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliLog.h"
#include "AliTrackContainer.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliVTrack.h"
 
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalQoverPtShift)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalQoverPtShift::AliAnalysisTaskEmcalQoverPtShift() :
  AliAnalysisTaskEmcalLight(),
  fHistos(nullptr),
  fQOverPtShift(0.),
  fTriggerBits(0),
  fTriggerString()
{
}

AliAnalysisTaskEmcalQoverPtShift::AliAnalysisTaskEmcalQoverPtShift(const char *name):
  AliAnalysisTaskEmcalLight(name, kTRUE),
  fHistos(nullptr),
  fQOverPtShift(0.),
  fTriggerBits(0),
  fTriggerString()
{
}

AliAnalysisTaskEmcalQoverPtShift::~AliAnalysisTaskEmcalQoverPtShift()
{
}

void AliAnalysisTaskEmcalQoverPtShift::UserCreateOutputObjects() {
  AliAnalysisTaskEmcalLight::UserCreateOutputObjects();

  fHistos = new THistManager("fHistos");

  std::array<std::string, 2> chargetypes = {{"pos", "neg"}};
  for(auto charge : chargetypes) {
    fHistos->CreateTH2(Form("fHistShift%s", charge.data()), Form("Pt-shift for charge %s; p_{t}^{orig} (GeV/c); p_{t}^{shift}", charge.data()), 300, 0., 300, 300, 0., 300.);  
  }
  
  for(auto hist : *fHistos->GetListOfHistograms()) fOutput->Add(hist);

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalQoverPtShift::Run() {
  auto tracks = GetTrackContainer("detector");
  if(!tracks) {
    AliErrorStream() << "No track container attached, not entering track loop" << std::endl;
    return kFALSE;
  }

  for(auto trk : tracks->accepted()) {
    Double_t chargeval = trk->Charge() > 0 ? 1 : -1;
    std::string chargestring = chargeval > 0 ? "pos" : "neg";
    Double_t ptorig = abs(trk->Pt());
    Double_t ptshift = 1./(fQOverPtShift*chargeval  + 1./ptorig);
    fHistos->FillTH2(Form("fHistShift%s", chargestring.data()), ptorig, ptshift);
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskEmcalQoverPtShift::IsTriggerSelected() {
  if(!(fInputHandler->IsEventSelected() & fTriggerBits)) return kFALSE;  
  if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerString)) return kFALSE;
  return kTRUE;
}

AliAnalysisTaskEmcalQoverPtShift *AliAnalysisTaskEmcalQoverPtShift::AddTaskQOverPtShift(const char *trigger, double shift) {
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    std::cerr << "Analysis Manager not initialized" << std::endl;
    return nullptr;
  }

  std::string track_name;
  auto inputhandler = mgr->GetInputEventHandler();
  if(inputhandler->InheritsFrom(AliAODInputHandler::Class())) {
    track_name = "tracks"; 
    std::cout << "Setting track container \"" << track_name << "\" for AODs" << std::endl;
  } else if(inputhandler->InheritsFrom(AliESDInputHandler::Class())) {
    track_name = "Tracks"; 
    std::cout << "Setting track container \"" << track_name << "\" for ESDs" << std::endl;
  } else {
    std::cerr << "Unknown input format or input handler not set" << std::endl;
    return nullptr;
  }

  UInt_t triggerbits(0);
  EMCAL_STRINGVIEW triggerstring(trigger);
  if(triggerstring == "INT7") triggerbits = AliVEvent::kINT7;
  else if(triggerstring == "EJ1" || triggerstring == "EJ2") triggerbits = AliVEvent::kEMCEJE;
  else if(triggerstring == "EG1" || triggerstring == "EG2") triggerbits = AliVEvent::kEMCEGA;
  if(!triggerbits) {
    std::cerr << "Unknown trigger class " << triggerstring << ", cannot configure task" << std::endl;
    return nullptr;
  }

  std::stringstream taskname;
  std::string qoverptstring(Form("%c%05d", shift > 0 ? 'p' : 'm', int(shift * 1e5)));
  taskname << "QOverPtTask_" <<  qoverptstring << "_" << trigger;
  auto task = new AliAnalysisTaskEmcalQoverPtShift(taskname.str().data());
  task->AddParticleContainer(track_name.data(), "detector");
  task->SetQOverPtShift(shift);
  task->SetTriggerSelection(triggerbits, trigger);
  mgr->AddTask(task);

  std::stringstream outputfile, listname;
  outputfile << mgr->GetCommonFileName() << ":QOverPtShift_" << qoverptstring << "_" << trigger;
  listname << "QOverPtShiftHistos_" << qoverptstring << "_" << trigger;

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(listname.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.str().data()));

  return task;
}