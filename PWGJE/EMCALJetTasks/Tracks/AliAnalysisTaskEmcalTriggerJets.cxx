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

#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliClusterContainer.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliPIDResponse.h"
#include "AliTrackContainer.h"
#include "AliAnalysisTaskEmcalTriggerJets.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVTrack.h"

#include <array>
#include <iostream>

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerJets)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalTriggerJets::AliAnalysisTaskEmcalTriggerJets():
    AliAnalysisTaskEmcalJet(),
    fPIDResponse(nullptr),
    fHistos(nullptr)
{

}

AliAnalysisTaskEmcalTriggerJets::AliAnalysisTaskEmcalTriggerJets(const char *name):
    AliAnalysisTaskEmcalJet(name, true),
    fPIDResponse(nullptr),
    fHistos(nullptr)
{

}

AliAnalysisTaskEmcalTriggerJets::~AliAnalysisTaskEmcalTriggerJets() {
  delete fHistos;
}

void AliAnalysisTaskEmcalTriggerJets::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  const std::array<TString, 5> kEmcalTriggers = {"INT7", "EJ1", "EJ2", "DJ1", "DJ2"};

  TLinearBinning jetptbinning(9, 20, 200), pbinning(300, 0., 30.), dEdxbinning(600, 0., 600.), massbinning(400., 0., 4.), eopbinning(150, 0., 1.5), radiusBinning(10, 0., 1.);
  const TBinning *binningPID[5] {&jetptbinning, &pbinning, &dEdxbinning, &massbinning, &eopbinning};
  const std::array<int, 2> kJetRadii = {2, 4};
  const std::array<TString, 3> kJetTypes = {"Charged", "Full", "Neutral"};
  const std::array<TString, 2> kDetectors = {"EMCAL", "DCAL"};
  const std::array<TString, 2> kConstituentType = {"Charged", "Neutral"};
  fHistos = new THistManager("EmcalJetHistos");
  for(auto t : kEmcalTriggers){
    fHistos->CreateTH1("hEventCount" + t, "Event counter for trigger " + t, 1., 0.5, 1.5);
    for(auto det: kDetectors){
      for(auto jt : kJetTypes) {
        for(auto radius : kJetRadii){
          fHistos->CreateTH1("hPtRaw" + jt + "JetR" + Form("%02d", radius) + det + t,
              "Raw pt spectrum for " + jt + " jets with R=" + Form("%.1f", double(radius)/10.) + " in " + det +" for trigger " + t,
              200, 0., 200);
          if(jt == "Full") {
            // Neutral energy fraction vs. jet pt
            fHistos->CreateTH2("hNefFullJetR" + TString::Format("%02d", radius) + det + t,
                "Neutral energy fraction vs. jet pt for R=" + TString::Format("%.1f", double(radius)/10.) + " full jets in " + det + " for trigger " + t,
                200, 0., 200., 100, 0., 1.);
            // pt, z and pt_rel of the leading particle (cluster and track) vs jet pt
            for(auto constituent : kConstituentType){
              fHistos->CreateTH2("hLeading" + constituent + "PtFullJetR" + TString::Format("%02d", radius) + det + t,
                  "p{t_const} vs. p_{t, jet} for leading " + constituent + " constituents in full jets with R=" + TString::Format("%.1f", double(radius)/10.) + " in " + det + "for trigger " + t,
                  200, 0., 200., 200, 0., 200.);
              fHistos->CreateTH2("hLeading" + constituent + "ZFullJetR" + TString::Format("%02d", radius) + det + t,
                  "z vs. p_{t, jet} for leading " + constituent + " constituents in full jets with R=" + TString::Format("%.1f", double(radius)/10.) + " in " + det + " for trigger " + t,
                  200, 0., 200., 100, 0., 1.);
              fHistos->CreateTH2("hLeading" + constituent + "PtrelFullJetR" + TString::Format("%02d", radius) + det + t,
                  "p_{t,rel} vs. p_{t, jet} for leading " + constituent + " constituents in full jets with R=" + TString::Format("%.1f", double(radius)/10.) + " in " + det + " for trigger " + t,
                  200., 0., 200., 100, 0., 1.);
            }
            fHistos->CreateTHnSparse("hPIDConstituentFullJet" + TString::Format("R%02d", radius) + det + t, TString::Format("PID for full jet constituents for R=%0.2f in %s for trigger %s", double(radius)/10., det.Data(), t.Data()), 5, binningPID);
            fHistos->CreateTHnSparse("hPIDLeadingFullJet" + TString::Format("R%02d", radius) + det + t, TString::Format("PID for full jet leading constituents for R=%0.2f in %s for trigger %s", double(radius)/10., det.Data(), t.Data()), 5, binningPID);
          }
        }
      }
    }
  }
  for(auto h : *(fHistos->GetListOfHistograms())){
    fOutput->Add(h);
  }
  PostData(1, fOutput);
}

void AliAnalysisTaskEmcalTriggerJets::UserExecOnce() {
  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse){
    AliErrorStream() << "PID Response not available - PID plots will not be filled" << std::endl;
  }
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
        double radius = double(TString(r).ReplaceAll("R0","").Atoi())/10.;
        TString namejcont = jt + "Jets" + r + det,
                rawhistnamebase = "hPtRaw" + jt + "Jet" + r + det,
                tagForLeading = "FullJet" + r + det;
        AliJetContainer *c = this->GetJetContainer(namejcont);
        if(!c) AliErrorStream() << "Not found jet container " << namejcont << std::endl;
        for(auto j : c->accepted()){
          TLorentzVector jetvec(j->Px(), j->Py(), j->Pt(), j->E());
          // get the leading particle
          AliVCluster *leadingClust = (jt == "Full") ? j->GetLeadingCluster(this->GetClusterContainer("caloClusters")->GetArray()) : nullptr;
          AliVTrack *leadingTrack = (jt == "Full") ? static_cast<AliVTrack *>(j->GetLeadingTrack(this->GetTrackContainer("tracks")->GetArray())) : nullptr;
          TLorentzVector clustervec;
          TVector3 trackvec;
          if(leadingClust){
            leadingClust->GetMomentum(clustervec, fVertex);
          }
          if(leadingTrack) {
            trackvec.SetXYZ(leadingTrack->Px(), leadingTrack->Py(), leadingTrack->Pz());
          }
          double absJetPt = TMath::Abs(j->Pt());
          for(auto t : triggers) {
            fHistos->FillTH1(rawhistnamebase + t, absJetPt);
            if(jt == "Full") {
              AliDebugStream(1) << "Filling full jet leading histograms, constituents c[" << j->GetNumberOfTracks() << "], n[" << j->GetNumberOfClusters() << "]" << std::endl;
              FillJetPIDPlots(j, radius,  t, det);
              fHistos->FillTH2("hNefFullJet" + r + det + t, absJetPt, j->NEF());
              // fill also histograms for the neutral energy and the leading particles
              if(leadingClust){
                fHistos->FillTH2("hLeadingNeutralPt" + tagForLeading + t, absJetPt, TMath::Abs(clustervec.Pt()));
                fHistos->FillTH2("hLeadingNeutralZ" + tagForLeading + t, absJetPt, j->GetZ(clustervec.Px(), clustervec.Py(), clustervec.Pz()));
                fHistos->FillTH2("hLeadingNeutralPtrel" + tagForLeading + t, absJetPt, clustervec.Pt(jetvec.Vect()));
              } else {
                AliDebugStream(2) << "No leading cluster found" << std::endl;
              }
              if(leadingTrack){
                fHistos->FillTH2("hLeadingChargedPt" + tagForLeading + t, absJetPt, TMath::Abs(leadingTrack->Pt()));
                fHistos->FillTH2("hLeadingChargedZ" + tagForLeading + t, absJetPt, j->GetZ(leadingTrack));
                fHistos->FillTH2("hLeadingChargedPtrel" + tagForLeading + t, absJetPt, trackvec.Pt(jetvec.Vect()));
                AliDebugStream(2) << "Filling full jet leading PID histograms" << std::endl;
                this->FillJetPIDPlotsLeading(leadingTrack, j->Pt(), radius, t, det);
              } else {
                AliDebugStream(2) << "No leading track found" << std::endl;
              }
              AliDebugStream(1) << "Filling full jet leading constituent histograms done" << std::endl;
            }
          }
        }
      }
    }
  }
  return true;
}

void AliAnalysisTaskEmcalTriggerJets::FillJetPIDPlots(const AliEmcalJet *jet, double radius, const char *trigger, const char *detector){
  if(TMath::Abs(jet->Pt()) < 20 || TMath::Abs(jet->Pt()) > 200) return;
  TString histname = TString::Format("hPIDConstituentFullJetR%02d%s%s", int(radius * 10), detector, trigger);
  AliTrackContainer *tc = GetTrackContainer("tracks");

  for(int icharged = 0; icharged < jet->GetNumberOfTracks(); icharged++){
    AliVTrack *constituent = static_cast<AliVTrack *>(jet->TrackAt(icharged, tc->GetArray()));
    // Select only constituents with sufficient PID information in both TPC and TOF
    if(constituent->GetTPCsignalN() < 30) continue;
    if(constituent->GetTPCSharedMapPtr()->CountBits(0) > 70) continue;
    if(!((constituent->GetStatus() & AliVTrack::kTOFout) && (constituent->GetStatus() & AliVTrack::kTIME))) continue;
    Double_t trtime = (constituent->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetTimeZero()) * 1e-12;
    Double_t v = constituent->GetIntegratedLength()/(100. * trtime);
    Double_t beta =  v / TMath::C(), gamma = 1 / TMath::Sqrt(1-beta*beta), mtof = constituent->P() / (beta*gamma);
    Double_t eop = -1.;
    if(constituent->GetEMCALcluster() >= 0) {
      AliVCluster *matched = static_cast<AliVCluster *>((*GetClusterContainer("caloClusters"))[constituent->GetEMCALcluster()]);
      if(matched) eop = TMath::Abs(matched->GetNonLinCorrEnergy()/constituent->P());
    }
    Double_t datapoint[5] = {TMath::Abs(jet->Pt()), TMath::Abs(constituent->P()), constituent->GetTPCsignal(), mtof, eop};
    fHistos->FillTHnSparse(histname, datapoint);
  }
}

void AliAnalysisTaskEmcalTriggerJets::FillJetPIDPlotsLeading(const AliVTrack *leading, double ptjet, double radius, const char *trigger, const char *detector){
  if(ptjet < 20 || ptjet > 200) return;

  TString histname = TString::Format("hPIDLeadingFullJetR%02d%s%s", int(radius*10.), detector, trigger);
  if(leading->GetTPCsignalN() < 30) return;
  if(!((leading->GetStatus() & AliVTrack::kTOFout) && (leading->GetStatus() & AliVTrack::kTIME))) return;
  if(leading->GetTPCSharedMapPtr()->CountBits(0) > 70) return;
  Double_t trtime = (leading->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetTimeZero()) * 1e-12;
  Double_t v = leading->GetIntegratedLength()/(100. * trtime);
  Double_t beta =  v / TMath::C(), gamma = 1 / TMath::Sqrt(1-beta*beta), mtof = leading->P() / (beta*gamma);
  Double_t eop = -1.;
  if(leading->GetEMCALcluster() >= 0){
      AliVCluster *matched = static_cast<AliVCluster *>((*GetClusterContainer("caloClusters"))[leading->GetEMCALcluster()]);
      if(matched) eop = TMath::Abs(matched->GetNonLinCorrEnergy()/leading->P());
  }
  Double_t datapoint[5] = {ptjet, TMath::Abs(leading->P()), leading->GetTPCsignal(), mtof, eop};
  fHistos->FillTHnSparse(histname, datapoint);
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
