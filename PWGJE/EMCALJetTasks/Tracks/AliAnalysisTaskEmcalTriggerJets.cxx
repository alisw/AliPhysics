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

/// \cond CLASSIMP
ClassImp(EmcalTriggerJets::AliAnalysisTaskEmcalTriggerJets)
/// \endcond

namespace EmcalTriggerJets {

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
  const int kNJetPtBins = 9;
  const int kNJetRadiusBins = 7;
  const std::array<int, kNJetPtBins+1> kJetPtBins = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200};
  const std::array<int, 2> kJetRadii = {2, 4};
  const std::array<TString, 3> kJetTypes = {"Charged", "Full", "Neutral"};
  const std::array<TString, 2> kDetectors = {"EMCAL", "DCAL"};
  const std::array<TString, 2> kConstituentType = {"Charged", "Neutral"};
  const std::array<double, kNJetRadiusBins+1> kRadiusBins = {0., 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.};
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
            // PID histograms for full jets (all and leading charged)
            for(int ib  = 0; ib < kNJetPtBins; ib++){
              fHistos->CreateTH2("hTPCdEdxFullJet" + TString::Format("Min%dMax%dR%02d", kJetPtBins[ib], kJetPtBins[ib+1], radius) + det + t,
                                 "TPC dE/dx vs. p for jet constituents for full jets with " + TString::Format("%d < p_{t} < %d and R=%.1f", kJetPtBins[ib], kJetPtBins[ib+1], float(radius)/10.) + "in " + det + " for trigger " + t,
                                 300, 0., 30., 1000., 0., 500.);
              fHistos->CreateTH2("hTOFBetaFullJet" + TString::Format("Min%dMax%dR%02d", kJetPtBins[ib], kJetPtBins[ib+1], radius) + det + t,
                                 "TOF #beta vs. p for jet constituents for full jets with " + TString::Format("%d < p_{t} < %d and R=%.1f", kJetPtBins[ib], kJetPtBins[ib+1], float(radius)/10.) + " in " + det + " for trigger " + t,
                                 300., 0., 30., 120., 0., 1.2);
              fHistos->CreateTH2("hEMCALEoPFullJet" + TString::Format("Min%dMax%dR%02d", kJetPtBins[ib], kJetPtBins[ib+1], radius) + det + t,
                                 "TOF #beta vs. p for jet constituents for full jets with " + TString::Format("%d < p_{t} < %d and R=%.1f", kJetPtBins[ib], kJetPtBins[ib+1], float(radius)/10.) + " in " + det + " for trigger " + t,
                                 300., 0., 30., 150., 0., 1.5);
              fHistos->CreateTH2("hTPCdEdxLeadingFullJet" + TString::Format("Min%dMax%dR%02d", kJetPtBins[ib], kJetPtBins[ib+1], radius) + det + t,
                                 "TPC dE/dx vs. p for leading jet constituents for full jets with " + TString::Format("%d < p_{t} < %d and R=%.1f", kJetPtBins[ib], kJetPtBins[ib+1], float(radius)/10.) + "in " + det + " for trigger " + t,
                                 300, 0., 30., 1000., 0., 500.);
              fHistos->CreateTH2("hTOFBetaLeadingFullJet" + TString::Format("Min%dMax%dR%02d", kJetPtBins[ib], kJetPtBins[ib+1], radius) + det + t,
                                 "TOF #beta vs. p for leading jet constituents for full jets with " + TString::Format("%d < p_{t} < %d and R=%.1f", kJetPtBins[ib], kJetPtBins[ib+1], float(radius)/10.) + " in " + det + " for trigger " + t,
                                 300., 0., 30., 120., 0., 1.2);
              fHistos->CreateTH2("hEMCALEoPLeadingFullJet" + TString::Format("Min%dMax%dR%02d", kJetPtBins[ib], kJetPtBins[ib+1], radius) + det + t,
                                 "TOF #beta vs. p for jet constituents for full jets with " + TString::Format("%d < p_{t} < %d and R=%.1f", kJetPtBins[ib], kJetPtBins[ib+1], float(radius)/10.) + " in " + det + " for trigger " + t,
                                 300., 0., 30., 150., 0., 1.5);

              // PID plots in bins of the jet radius
              for(int ir = 0; ir < kNJetRadiusBins; ir++){
                fHistos->CreateTH2("hTPCdEdxFullJet" + TString::Format("DMin%dDMax%dPtMin%dPtMax%d", int(kRadiusBins[ir]*10.), int(kRadiusBins[ir+1]*10.), kJetPtBins[ib], kJetPtBins[ib+1]) + det + t,
                    "TPC dE/dx vs. p for particles with " + TString::Format("%.1f < d < %.1f and %d < p_{t} < %d", kRadiusBins[ir], kRadiusBins[ir+1], kJetPtBins[ib], kJetPtBins[ib+1]) + " in " + det + " for trigger " + t,
                    300, 0., 30., 1000., 0., 500.);
                fHistos->CreateTH2("hTOFBetaFullJet" + TString::Format("DMin%dDMax%dPtMin%dPtMax%d", int(kRadiusBins[ir]*10.), int(kRadiusBins[ir+1]*10.), kJetPtBins[ib], kJetPtBins[ib+1]) + det + t,
                    "TOF #beta vs. p for particles with " + TString::Format("%.1f < d < %.1f and %d < p_{t} < %d", kRadiusBins[ir], kRadiusBins[ir+1], kJetPtBins[ib], kJetPtBins[ib+1]) + " in " + det + " for trigger " + t,
                    300, 0., 30., 150., 0., 1.5);
                fHistos->CreateTH2("hEMCALEoPFullJet" + TString::Format("DMin%dDMax%dPtMin%dPtMax%d", int(kRadiusBins[ir]*10.), int(kRadiusBins[ir+1]*10.), kJetPtBins[ib], kJetPtBins[ib+1]) + det + t,
                    "TPC dE/dx vs. p for particles with " + TString::Format("%.1f < d < %.1f and %d < p_{t} < %d", kRadiusBins[ir], kRadiusBins[ir+1], kJetPtBins[ib], kJetPtBins[ib+1]) + " in " + det + " for trigger " + t,
                    300, 0., 30., 150., 0., 1.5);
              }
            }
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
        bool doPID = fPIDResponse && (jt == "Full") && (r == "R04");
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
              // Fill PID plots for particles around the main jet axis
              if(radius == 0.4) {
                FillPIDCorrelationPlot(j, this->GetTrackContainer("tracks"), t, det);
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
  const int kNJetPtBins = 9;
  const std::array<int, kNJetPtBins+1> kJetPtBins = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200};
  int jetptbin = -1;
  for(int ib = 0; ib < kNJetPtBins; ib++){
    if(TMath::Abs(jet->Pt()) >= kJetPtBins[ib] && TMath::Abs(jet->Pt()) < kJetPtBins[ib+1]) {
      jetptbin = ib;
      break;
    }
  }
  if(jetptbin < 0) return;
  TString histnameTPC = TString::Format("hTPCdEdxFullJetMin%dMax%dR%02d%s%s", kJetPtBins[jetptbin], kJetPtBins[jetptbin+1], int(radius*10.), detector, trigger),
          histnameTOF = TString::Format("hTOFBetaFullJetMin%dMax%dR%02d%s%s", kJetPtBins[jetptbin], kJetPtBins[jetptbin+1], int(radius*10.), detector, trigger),
          histnameEMCAL = TString::Format("hEMCALEoPFullJetMin%dMax%dR%02d%s%s", kJetPtBins[jetptbin], kJetPtBins[jetptbin+1], int(radius*10.), detector, trigger);
  AliTrackContainer *tc = GetTrackContainer("tracks");

  for(int icharged = 0; icharged < jet->GetNumberOfTracks(); icharged++){
    AliVTrack *constituent = static_cast<AliVTrack *>(jet->TrackAt(icharged, tc->GetArray()));
    // Select only constituents with sufficient PID information in both TPC and TOF
    if(constituent->GetTPCsignalN() < 30) continue;
    if(constituent->GetTPCSharedMapPtr()->CountBits(0) > 70) continue;
    if(!((constituent->GetStatus() & AliVTrack::kTOFout) && (constituent->GetStatus() & AliVTrack::kTIME))) continue;
    Double_t trtime = (constituent->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetTimeZero()) * 1e-12;
    Double_t v = constituent->GetIntegratedLength()/(100. * trtime);
    Double_t beta =  v / TMath::C();
    fHistos->FillTH2(histnameTPC, TMath::Abs(constituent->P()), constituent->GetTPCsignal());
    fHistos->FillTH2(histnameTOF, TMath::Abs(constituent->P()), beta);
    if(constituent->GetEMCALcluster() >= 0) {
      AliVCluster *matched = static_cast<AliVCluster *>((*GetClusterContainer("caloClusters"))[constituent->GetEMCALcluster()]);
      fHistos->FillTH2(histnameEMCAL, TMath::Abs(constituent->P()), TMath::Abs(matched->GetNonLinCorrEnergy()/constituent->P()));
    }
  }
}

void AliAnalysisTaskEmcalTriggerJets::FillJetPIDPlotsLeading(const AliVTrack *leading, double ptjet, double radius, const char *trigger, const char *detector){
  const int kNJetPtBins = 9;
  const std::array<int, kNJetPtBins+1> kJetPtBins = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200};
  int jetptbin = -1;
  for(int ib = 0; ib < kNJetPtBins; ib++){
    if(TMath::Abs(ptjet) >= kJetPtBins[ib] && TMath::Abs(ptjet) < kJetPtBins[ib+1]) {
      jetptbin = ib;
      break;
    }
  }
  if(jetptbin < 0) return;

  TString histnameTPC = TString::Format("hTPCdEdxLeadingFullJetMin%dMax%dR%02d%s%s", kJetPtBins[jetptbin], kJetPtBins[jetptbin+1], int(radius*10.), detector, trigger),
          histnameTOF = TString::Format("hTOFBetaLeadingFullJetMin%dMax%dR%02d%s%s", kJetPtBins[jetptbin], kJetPtBins[jetptbin+1], int(radius*10.), detector, trigger),
          histnameEMCAL = TString::Format("hEMCALEoPLeadingFullJetMin%dMax%dR%02d%s%s", kJetPtBins[jetptbin], kJetPtBins[jetptbin+1], int(radius*10.), detector, trigger);
  if(leading->GetTPCsignalN() < 30) return;
  if(!((leading->GetStatus() & AliVTrack::kTOFout) && (leading->GetStatus() & AliVTrack::kTIME))) return;
  if(leading->GetTPCSharedMapPtr()->CountBits(0) > 70) return;
  Double_t trtime = (leading->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetTimeZero()) * 1e-12;
  Double_t v = leading->GetIntegratedLength()/(100. * trtime);
  Double_t beta =  v / TMath::C();
  fHistos->FillTH2(histnameTPC, TMath::Abs(leading->P()), leading->GetTPCsignal());
  fHistos->FillTH2(histnameTOF, TMath::Abs(leading->P()), beta);
  if(leading->GetEMCALcluster() >= 0){
      AliVCluster *matched = static_cast<AliVCluster *>((*GetClusterContainer("caloClusters"))[leading->GetEMCALcluster()]);
      fHistos->FillTH2(histnameEMCAL, TMath::Abs(leading->P()), TMath::Abs(matched->GetNonLinCorrEnergy()/leading->P()));
  }
}

void AliAnalysisTaskEmcalTriggerJets::FillPIDCorrelationPlot(const AliEmcalJet *jet, const AliTrackContainer *particles, const char *trigger, const char *detector) {
  const int kNJetPtBins = 9;
  const int kNRadiusBins = 8;
  const std::array<int, kNJetPtBins+1> kJetPtBins = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200};
  const std::array<double, kNRadiusBins+1> kRadiusBins = {0., 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.};
  int jetptbin = -1;
  for(int ib = 0; ib < kNJetPtBins; ib++){
    if(TMath::Abs(jet->Pt()) >= kJetPtBins[ib] && TMath::Abs(jet->Pt()) < kJetPtBins[ib+1]) {
      jetptbin = ib;
      break;
    }
  }
  if(jetptbin < 0) return;

  for(auto track : particles->all()) {
   if(track->IsA() == AliAODTrack::Class()) {
      AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
      if(!(aodtrack->IsHybridGlobalConstrainedGlobal() || aodtrack->IsHybridTPCConstrainedGlobal())) continue;
      if(track->GetTPCsignalN() < 30) continue;
      if(track->GetTPCSharedMapPtr()->CountBits(0) > 70) continue; // exclude tracks with more than 70 shared clusters ->
      if(!((track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME))) continue;
      TVector3 partvector(track->Px(), track->Py(), track->Pz()), jetvector(jet->Px(), jet->Py(), jet->Pz());
      double dr = jetvector.DeltaR(partvector);
      if(dr >= 1.) continue;

      // track is accepted, find delta_r bin
      int drbin = 0;
      for(int ib = 0; ib < kNRadiusBins; ib++){
        if(dr >= kRadiusBins[ib] && dr < kRadiusBins[ib+1]) {
          drbin = ib;
          break;
        }
      }
      TString histnameTPC = TString::Format("hTPCdEdxFullJetDMin%dDMax%dPtMin%dPtMax%d%s%s", static_cast<int>(kRadiusBins[drbin] * 10.), static_cast<int>(kRadiusBins[drbin+1] * 10.), kJetPtBins[jetptbin], kJetPtBins[jetptbin+1], detector, trigger),
              histnameTOF = TString::Format("hTOFBetaFullJetDMin%dDMax%dPtMin%dPtMax%d%s%s", static_cast<int>(kRadiusBins[drbin] * 10.), static_cast<int>(kRadiusBins[drbin+1] * 10.), kJetPtBins[jetptbin], kJetPtBins[jetptbin+1], detector, trigger),
              histnameEMCAL = TString::Format("hEMCALEoPFullJetDMin%dDMax%dPtMin%dPtMax%d%s%s", static_cast<int>(kRadiusBins[drbin] * 10.), static_cast<int>(kRadiusBins[drbin+1] * 10.), kJetPtBins[jetptbin], kJetPtBins[jetptbin+1], detector, trigger);
      Double_t trtime = (track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetTimeZero()) * 1e-12;
      Double_t v = track->GetIntegratedLength()/(100. * trtime);
      Double_t beta =  v / TMath::C();
      fHistos->FillTH2(histnameTPC, TMath::Abs(track->P()), track->GetTPCsignal());
      fHistos->FillTH2(histnameTOF, TMath::Abs(track->P()), beta);
      if(track->GetEMCALcluster() >= 0){
          AliVCluster *matched = static_cast<AliVCluster *>((*GetClusterContainer("caloClusters"))[track->GetEMCALcluster()]);
          fHistos->FillTH2(histnameEMCAL, TMath::Abs(track->P()), TMath::Abs(matched->GetNonLinCorrEnergy()/track->P()));
      }
    }
  }
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

} /* namespace EmcalTriggerJets */
