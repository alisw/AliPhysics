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
#include <iostream>
#include <THistManager.h>
#include <TList.h>
#include <TMath.h>
#include <TParticle.h>
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliAnalysisTaskGammaConvDtrue.h"

ClassImp(AliAnalysisTaskGammaConvDtrue)

AliAnalysisTaskGammaConvDtrue::AliAnalysisTaskGammaConvDtrue():
  AliAnalysisTaskEmcal(),
  fHistos(nullptr)
{

}

AliAnalysisTaskGammaConvDtrue::AliAnalysisTaskGammaConvDtrue(const char *name):
  AliAnalysisTaskEmcal(name, true),
  fHistos(nullptr)
{
  this->SetMakeGeneralHistograms(true);
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaConvDtrue::~AliAnalysisTaskGammaConvDtrue(){
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskGammaConvDtrue::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fHistos = new THistManager(Form("Histos_%s", GetName()));

  fHistos->CreateTH1("hTrueD0pt", "True D0 pt spectrum", 200, 0., 200.);
  fHistos->CreateTH1("hTrueD0pi0piplpimi", "True D0 pt spectrum (3 pion channel)", 200, 0., 200.);
  fHistos->CreateTH1("hTrueD0pi0piplpimiAcc",  "True D0 pt spectrum (3 pion channel, in acceptance)", 200, 0., 200.);
  fHistos->CreateTH2("hPi0FromD0pt", "pt of the pi0 from D0 decay", 200, 0., 200., 200, 0., 200.);
  fHistos->CreateTH2("hPiPlFromD0pt", "pt of the pi+ from D0 decay", 200, 0., 200., 200, 0., 200.);
  fHistos->CreateTH2("hPiMiFromD0pt", "pt of the pi- from D0 decay", 200, 0., 200., 200, 0., 200.);

  for(auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
  PostData(1, fOutput);
}

bool AliAnalysisTaskGammaConvDtrue::Run(){
  if(!fMCEvent) {
    AliErrorStream() << "No MC event available" << std::endl;
    return false;
  }
  for(int ipart = 0; ipart < fMCEvent->GetNumberOfTracks(); ipart++){
    auto part = fMCEvent->GetTrack(ipart);
    if(TMath::Abs(part->PdgCode()) != 421) continue;
    if(TMath::Abs(part->Eta()) < 0.5) continue;
    fHistos->FillTH1("hTrueD0pt", part->Pt());
    // check for 3 pion decay
    auto decay = CheckFor3PionDecay(part);
    if(!(decay.fPi0daughter && decay.fPiMidaughter && decay.fPiPldaughter)) continue;
    fHistos->FillTH1("hTrueD0pi0piplpimi", part->Pt());
    // check whether all daughters are in the acceptance;
    if(TMath::Abs(decay.fPiPldaughter->Eta()) > 0.8) continue;
    if(TMath::Abs(decay.fPiMidaughter->Eta()) > 0.8) continue;
    if(TMath::Abs(decay.fPi0daughter->Eta()) > 0.7) continue;   // EMCAL
    fHistos->FillTH1("hTrueD0pi0piplpimiAcc", part->Pt());

    fHistos->FillTH2("hPi0FromD0pt", part->Pt(), decay.fPi0daughter->Pt());
    fHistos->FillTH2("hPiPlFromD0pt", part->Pt(), decay.fPiPldaughter->Pt());
    fHistos->FillTH2("hPiMiFromD0pt", part->Pt(), decay.fPiMidaughter->Pt());
  }
  return true;
}

AliAnalysisTaskGammaConvDtrue::D03PionDecay AliAnalysisTaskGammaConvDtrue::CheckFor3PionDecay(AliVParticle *d0mother) const {
  D03PionDecay decay;
  decay.fPi0daughter = nullptr;
  decay.fPiPldaughter = nullptr;
  decay.fPiMidaughter = nullptr;
  AliDebug(1, "Found D0\n");
  int did = 0;
  for(auto idaughter = d0mother->GetDaughterFirst(); idaughter <= d0mother->GetDaughterLast(); idaughter++){
    auto daughterpart = fMCEvent->GetTrack(idaughter);
    AliDebug(2, Form("Daughter %d: %d\n", did, daughterpart->PdgCode()));
    if(TMath::Abs(daughterpart->PdgCode()) == 111)  decay.fPi0daughter = daughterpart;  
    if(daughterpart->PdgCode() == 211) decay.fPiPldaughter = daughterpart;
    if(daughterpart->PdgCode() == -211) decay.fPiMidaughter = daughterpart;
  }

  return decay;
}

AliAnalysisTaskGammaConvDtrue *AliAnalysisTaskGammaConvDtrue::AddTaskGammaConvDtrue(){
  auto mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    std::cerr << "AliAnalysisTaskGammaConvDtrue: No AnalysisManager defined" << std::endl;
    return nullptr;
  }

  auto d0task = new AliAnalysisTaskGammaConvDtrue("d0task");
  mgr->AddTask(d0task);

  TString outputcont = mgr->GetCommonFileName();
  outputcont += ":D0MCtrue";
  
  mgr->ConnectInput(d0task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(d0task, 1, mgr->CreateContainer("D0true", TList::Class(), AliAnalysisManager::kOutputContainer, outputcont));
  return d0task;
}
