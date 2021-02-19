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
#include <TLinearBinning.h>
#include <TMath.h>
#include <TString.h>
#include <TVector3.h>

#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliClusterContainer.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliPIDResponse.h"
#include "AliTrackContainer.h"
#include "AliAnalysisTaskEmcalTriggerJetsIDcorr.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVTrack.h"

#include <array>
#include <iostream>

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerJetsIDcorr)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalTriggerJetsIDcorr::AliAnalysisTaskEmcalTriggerJetsIDcorr():
    AliAnalysisTaskEmcalJet(),
    fJetCont(nullptr),
    fPIDResponse(nullptr),
    fHistos(nullptr)
{

}

AliAnalysisTaskEmcalTriggerJetsIDcorr::AliAnalysisTaskEmcalTriggerJetsIDcorr(const char *name):
    AliAnalysisTaskEmcalJet(name, true),
    fJetCont(nullptr),
    fPIDResponse(nullptr),
    fHistos(nullptr)
{

}

AliAnalysisTaskEmcalTriggerJetsIDcorr::~AliAnalysisTaskEmcalTriggerJetsIDcorr() {
  delete fHistos;
}

void AliAnalysisTaskEmcalTriggerJetsIDcorr::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  const std::array<TString, 3> kEmcalTriggers = {"INT7", "EJ1", "EJ2"};

  TLinearBinning jetptbinning(9, 20, 200), pbinning(300, 0., 30.), massbinning(600., 0., 6.), radiusBinning(10, 0., 1.);
  const TBinning *binningPID[4] = {&jetptbinning, &pbinning, &radiusBinning, &massbinning};
  const std::array<TString, 2> kSpecies = {"Proton", "Deuteron"};
  fHistos = new THistManager("EmcalJetHistos");
  fHistos->CreateTH1("hTPCdEdxErrors", "Error counter for TPC dE/dx", AliPID::kSPECIESC, -0.5, AliPID::kSPECIESC - 0.5);
  for(auto t : kEmcalTriggers){
    fHistos->CreateTH1("hEventCount" + t, "Number of events for trigger " + t, 1, 0.5, 1.5);
    fHistos->CreateTH1("hPtRawAllJet" + t, " Raw jet pt spectrum for trigger " + t, 1000, 0., 1000.);
    fHistos->CreateTH1("hPtRawSelJet" + t, " Raw jet pt spectrum for trigger " + t, 1000, 0., 1000.);
    fHistos->CreateTH1("hNJetsAll" + t, "All reconstructed jets for trigger " + t, 101, -0.5, 100.5);
    fHistos->CreateTH1("hNJetsSelected" + t, "Selected reconstructed jets for trigger " + t, 101, -0.5, 101.5);
    for(auto s: kSpecies){
      fHistos->CreateTH1("hNCandidatesPerEvent" + s + t, "Number of " + s + " candidates per events ", 101, -0.5, 100.5);
      fHistos->CreateTHnSparse("hPIDAssociateFullJet" + s + t, TString::Format("PID (TOF masss) for particles associated to full jets for R=0.4f in EMCAL for %s candidates for trigger %s", s.Data(), t.Data()), 4, binningPID);
    }
  }
  for(auto h : *(fHistos->GetListOfHistograms())){
    fOutput->Add(h);
  }
  PostData(1, fOutput);
}

void AliAnalysisTaskEmcalTriggerJetsIDcorr::UserExecOnce() {
  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse){
    AliErrorStream() << "PID Response not available - PID plots will not be filled" << std::endl;
  }
  fJetCont = this->GetJetContainer("FullJetsR04EMCAL");
}

bool AliAnalysisTaskEmcalTriggerJetsIDcorr::Run(){
  std::vector<TString> triggers, kEmcalTriggers = {"EJ1", "EJ2"};
  if(fInputHandler->IsEventSelected() & AliVEvent::kINT7) triggers.push_back("INT7");
  if(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE){
    TString fired = fInputEvent->GetFiredTriggerClasses();
    for(auto e : kEmcalTriggers){
      if(fired.Contains(e))
        triggers.push_back(e);
    }
  }
  if(!triggers.size()) return false;

  for(const auto &t : triggers) fHistos->FillTH1("hEventCount" + t, 1);

  // Get (rough) TPC candidates for protons and deuterons
  std::vector<AliVTrack *> protonCandidates , deuteronCandidates;
  try{
    protonCandidates = this->GetTPCPIDCandidates(AliPID::kProton);
  } catch(TPCdEdxException &e) {
    fHistos->FillTH1("hTPCdEdxErrors", e.GetParticleType());
  }
  try{
    deuteronCandidates = this->GetTPCPIDCandidates(AliPID::kDeuteron);
  } catch(TPCdEdxException &e) {
    fHistos->FillTH1("hTPCdEdxErrors", e.GetParticleType());
  }
  for(const auto &t : triggers){
    fHistos->FillTH1("hNCandidatesPerEventProton" + t, protonCandidates.size());
    fHistos->FillTH1("hNCandidatesPerEventDeuteron" + t, deuteronCandidates.size());
  }

  int njetsAll(0), njetsSel(0);
  for(auto j : fJetCont->all()) {
    njetsAll++;
    for(const auto &t : triggers) fHistos->FillTH1("hPtRawAllJet" + t, j->Pt());
  }
  for(auto j : fJetCont->accepted()){
    njetsSel++;
    for(const auto &t : triggers) fHistos->FillTH1("hPtRawSelJet" + t, j->Pt());
    TVector3 jetmomentum(j->Px(), j->Py(), j->Pz());
    // Proton-jet correlations
    std::vector<CorrParticleInfo> protonsCorrelated = CorrelateCandidatesToJet(jetmomentum, protonCandidates);
    for(auto c : protonsCorrelated) {
      Double_t point[4] = {TMath::Abs(j->Pt()), c.fPt, c.fDR, c.fMass*c.fMass};
      for(auto t : triggers) fHistos->FillTHnSparse("hPIDAssociateFullJetProton" +t, point);
    }

    // deuteron-jet correlations
    std::vector<CorrParticleInfo> deuteronsCorrelated = CorrelateCandidatesToJet(jetmomentum, deuteronCandidates);
    for(auto c : deuteronsCorrelated) {
      Double_t point[4] = {TMath::Abs(j->Pt()), c.fPt, c.fDR, c.fMass*c.fMass};
      for(auto t : triggers) fHistos->FillTHnSparse("hPIDAssociateFullJetDeuteron" +t, point);
    }
  }
  for(const auto &t : triggers){
    fHistos->FillTH1("hNJetsAll" + t, njetsAll);
    fHistos->FillTH1("hNJetsSelected" + t, njetsSel);
  }

  return true;
}

std::vector<CorrParticleInfo> AliAnalysisTaskEmcalTriggerJetsIDcorr::CorrelateCandidatesToJet(const TVector3 &jet, std::vector<AliVTrack *> candidates) const {
  std::vector<CorrParticleInfo> correlatedParticles;
  for(auto c : candidates) {
    TVector3 particleMom(c->Px(), c->Py(), c->Pz());
    double dR = particleMom.DeltaR(jet);
    if(dR > 1.) continue;
    double mTOF = -1;
    try {
      mTOF = GetTOFMass(c);
    } catch(TOFMassException &e) {
      continue;
    }
    correlatedParticles.push_back({c->Pt(), dR, mTOF});
  }
  return correlatedParticles;
}

double AliAnalysisTaskEmcalTriggerJetsIDcorr::GetTOFMass(const AliVTrack *const track) const {
  if(!((track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME))) throw TOFMassException();
  Double_t trtime = (track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetTimeZero()) * 1e-12;
  Double_t v = track->GetIntegratedLength()/(100. * trtime);
  Double_t beta =  v / TMath::C(), gamma = 1 / TMath::Sqrt(1-beta*beta);
  return track->P() / (beta*gamma);
}

std::vector<AliVTrack *>  AliAnalysisTaskEmcalTriggerJetsIDcorr::GetTPCPIDCandidates(AliPID::EParticleType type) const {
  std::vector<AliVTrack *> result;
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); itrk++){
    AliVTrack *trk = static_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
    if(!trk) continue;
    if(trk->IsA() == AliAODTrack::Class()){
      AliAODTrack *aodtrk = static_cast<AliAODTrack *>(trk);
      if(!(aodtrk->IsHybridGlobalConstrainedGlobal() || aodtrk->IsHybridTPCConstrainedGlobal())) continue;
      if(aodtrk->GetTPCsignalN() < 30) continue;
      if(aodtrk->GetTPCSharedMapPtr()->CountBits(0) > 70) continue; // exclude tracks with more than 70 shared clusters ->
      Double_t nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(trk, type);
      if(nSigmaTPC <= -999.) continue; //throw TPCdEdxException(type);
      if(TMath::Abs(nSigmaTPC) > 3) continue;
      result.push_back(trk);
    }
  }
  return result;
}

AliAnalysisTaskEmcalTriggerJetsIDcorr *AliAnalysisTaskEmcalTriggerJetsIDcorr::AddTaskEmcalTriggerJetsIDcorr(const char *name){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisTaskEmcalTriggerJetsIDcorr *task = new AliAnalysisTaskEmcalTriggerJetsIDcorr(name);
  mgr->AddTask(task);

  // Adding containers for clusters and tracks
  AliClusterContainer *clustercont = task->AddClusterContainer("caloClusters");
  clustercont->SetMinE(0.3);
  AliTrackContainer *trackcont = task->AddTrackContainer("tracks");
  trackcont->SetMinPt(0.15);

  // Adding Jet containers
  // - Using jet radius 0.4
  // - Only EMCAL side
  // - processing only full jets


  // Full jets, R=0.4, EMCAL
  AliJetContainer *cont = task->AddJetContainer(
                              AliJetContainer::kFullJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              0.4,
                              AliEmcalJet::kEMCALfid,
                              trackcont, clustercont
                              );
  cont->SetName("FullJetsR04EMCAL");
  cont->SetJetPtCut(20.);
  cont->SetJetPtCutMax(200.);


  // Connect Input / Output containers
  TString outfilename = mgr->GetCommonFileName();
  outfilename += ":EmcalTriggerJetsIDcorr";
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("HistsEmcalTriggerJetsIDcorr", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename));

  return task;
}
