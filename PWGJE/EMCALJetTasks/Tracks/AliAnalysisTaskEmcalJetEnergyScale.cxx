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
#include <algorithm>
#include <array>
#include <string>
#include <sstream>
#include <vector>
#include <THistManager.h>
#include <TCustomBinning.h>
#include <TLinearBinning.h>
#include <TCustomBinning.h>
#include <TRandom.h>

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetEnergyScale.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliLog.h"
#include "AliVEventHandler.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalJetEnergyScale::AliAnalysisTaskEmcalJetEnergyScale():
  AliAnalysisTaskEmcalJet(),
  fHistos(nullptr),
  fNameDetectorJets(),
  fNameParticleJets(),
  fTriggerSelectionString(),
  fNameTriggerDecisionContainer("EmcalTriggerDecision"),
  fFractionResponseClosure(0.8),
  fFillHSparse(false),
  fScaleShift(0.),
  fSampleSplitter(nullptr)
{
}

AliAnalysisTaskEmcalJetEnergyScale::AliAnalysisTaskEmcalJetEnergyScale(const char *name):
  AliAnalysisTaskEmcalJet(name, true),
  fHistos(nullptr),
  fNameDetectorJets(),
  fNameParticleJets(),
  fTriggerSelectionString(),
  fNameTriggerDecisionContainer("EmcalTriggerDecision"),
  fFractionResponseClosure(0.8),
  fFillHSparse(false),
  fScaleShift(0.),
  fSampleSplitter(nullptr)
{
  SetUseAliAnaUtils(true);
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEmcalJetEnergyScale::~AliAnalysisTaskEmcalJetEnergyScale() {
  if(fHistos) delete fHistos;
  if(fSampleSplitter) delete fSampleSplitter;
}

void AliAnalysisTaskEmcalJetEnergyScale::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  TCustomBinning jetPtBinningCoarseDet, jetPtBinningCoarsePart;
  jetPtBinningCoarseDet.SetMinimum(20.);
  jetPtBinningCoarseDet.AddStep(40., 2.);
  jetPtBinningCoarseDet.AddStep(60., 5.);
  jetPtBinningCoarseDet.AddStep(120., 10.);
  jetPtBinningCoarseDet.AddStep(200., 20.);
  jetPtBinningCoarsePart.SetMinimum(0);
  jetPtBinningCoarsePart.AddStep(20., 20.);
  jetPtBinningCoarsePart.AddStep(80., 10.);
  jetPtBinningCoarsePart.AddStep(200., 20.);
  jetPtBinningCoarsePart.AddStep(280., 40.);
  jetPtBinningCoarsePart.AddStep(500., 220.);

  fHistos = new THistManager("energyScaleHistos");
  fHistos->CreateTH1("hEventCounter", "Event counter", 1, 0.5, 1.5);
  fHistos->CreateTH2("hJetEnergyScale", "Jet Energy scale; p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleDet", "Jet Energy scale (det); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetResponseFine", "Response matrix, fine binning", 350, 0., 350., 800, 0., 800.);
  fHistos->CreateTH2("hJetResponseFineClosure", "Response matrix, fine binning, for closure test", 350, 0., 350., 800, 0., 800.);
  fHistos->CreateTH2("hJetResponseFineNoClosure", "Response matrix, fine binning, for closure test", 350, 0., 350., 800, 0., 800.);
  fHistos->CreateTH1("hJetSpectrumPartAll", "Part level jet pt spectrum ", 800, 0., 800.);
  if(fFillHSparse){
    TLinearBinning jetPtBinningDet(350, 0., 350.), jetPtBinningPart(600, 0., 600), nefbinning(100, 0., 1.), ptdiffbinning(200, -1., 1.), jetEtaBinning(100, -0.9, 0.9), jetPhiBinning(100, 0., TMath::TwoPi()),
                   subsampleBinning(2, -0.5, 1.5), deltaRbinning(20, 0., 1.), statusbinningEff(3, -0.5, 2.5);

    const TBinning *diffbinning[3] = {&jetPtBinningPart, &nefbinning, &ptdiffbinning},
                   *corrbinning[6] = {&jetPtBinningPart, &jetPtBinningDet, &nefbinning, &deltaRbinning,&subsampleBinning,&subsampleBinning},
                   *effbinning[3] = {&jetPtBinningPart, &jetPtBinningDet, &statusbinningEff};

    fHistos->CreateTHnSparse("hPtDiff", "pt diff det/part", 3, diffbinning, "s");
    fHistos->CreateTHnSparse("hPtCorr", "Correlation det pt / part pt", 6, corrbinning, "s");
    fHistos->CreateTHnSparse("hJetfindingEfficiency", "Jet finding efficiency", 3, effbinning, "s");
  }

  // A bit of QA stuff
  fHistos->CreateTH2("hQANEFPtPart", "Neutral energy fraction at part. level; p_{t} (GeV/c); NEF", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQANEFPtDet", "Neutral energy fraction at det. level; p_{t} (GeV/c); NEF", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZchPtPart", "z_{ch,max} at part. level; p_{t} (GeV/c); z_{ch,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZchPtDet", "z_{ch,max} at det. level; p_{t} (GeV/c); z_{ch,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZnePtPart", "z_{ne,max} at part. level; p_{t} (GeV/c); z_{ne,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZnePtDet", "z_{ne,max} at det. level; p_{t} (GeV/c); z_{ne,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQANChPtPart", "Number of charged constituents at part. level; p_{t} (GeV/c); N_{ch}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQANChPtDet", "Number of charged constituents at det. level; p_{t} (GeV/c); N_{ch}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQANnePtPart", "Number of neutral constituents at part. level; p_{t} (GeV/c); N_{ne}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQANnePtDet", "Number of neutral constituents at det. level; p_{t} (GeV/c); N_{ne}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQAConstPtChPart", "p_{t} of charged constituents (part. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtChDet", "p_{t} of charged constituents (det. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtNePart", "p_{t} of neutral constituents (part. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtNeDet", "p_{t} of neutral constituents (det. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtChMaxPart", "p_{t} of max charged constituents (part. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtChMaxDet", "p_{t} of max charged constituents (det. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtNeMaxPart", "p_{t} of max neutral constituents (part. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtNeMaxDet", "p_{t} of max neutral constituents (det. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAEtaPhiPart", "#eta vs. #phi for selected part. level jets; #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiDet", "#eta vs. #phi for selected det. level jets; #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstChPart", "#eta vs. #phi for charged constituents (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstChDet", "#eta vs. #phi for charged constituents (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstNePart", "#eta vs. #phi for neutral constituents (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstNeDet", "#eta vs. #phi for neutral constituents (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstMaxChPart", "#eta vs. #phi for max charged constituents (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstMaxChDet", "#eta vs. #phi for max charged constituents (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstMaxNePart", "#eta vs. #phi for max neutral constituents (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstMaxNeDet", "#eta vs. #phi for max neutral constituents (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQADeltaRChargedPart", "#DeltaR vs. p_{t,jet} of charged constituents (part. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1.);
  fHistos->CreateTH2("hQADeltaRChargedDet", "#DeltaR vs. p_{t,jet} of charged constituents (det. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1.);
  fHistos->CreateTH2("hQADeltaRNeutralPart", "#DeltaR vs. p_{t,jet} of neutral constituents (part. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1);
  fHistos->CreateTH2("hQADeltaRNeutralDet", "#DeltaR vs. p_{t,jet} of neutral constituents (det. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1);
  fHistos->CreateTH2("hQADeltaRMaxChargedPart", "#DeltaR vs. p_{t,jet} of charged constituents (part. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1.);
  fHistos->CreateTH2("hQADeltaRMaxChargedDet", "#DeltaR vs. p_{t,jet} of charged constituents (det. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1.);
  fHistos->CreateTH2("hQADeltaRMaxNeutralPart", "#DeltaR vs. p_{t,jet} of neutral constituents (part. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1);
  fHistos->CreateTH2("hQADeltaRMaxNeutralDet", "#DeltaR vs. p_{t,jet} of neutral constituents (det. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1);
  fHistos->CreateTH2("hQAJetAreaVsJetPtPart", "Jet area vs. jet pt at particle level; p_{t} (GeV/c); Area", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAJetAreaVsJetPtDet", "Jet area vs. jet pt at detector level; p_{t} (GeV/c); Area", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAJetAreaVsNEFPart", "Jet area vs. NEF at particle level; NEF; Area", 100, 0., 1., 100, 0., 1.);
  fHistos->CreateTH2("hQAJetAreaVsNEFDet", "Jet area vs. NEF at detector level; NEF; Area", 100, 0., 1., 100, 0., 1.);
  fHistos->CreateTH2("hQAJetAreaVsNConstPart", "Jet area vs. number of consituents at particle level; Number of constituents; Area", 101, -0.5, 100.5, 100, 0., 1.);
  fHistos->CreateTH2("hQAJetAreaVsNConstDet", "Jet area vs. number of consituents at detector level; Number of constituents; Area", 101, -0.5, 100.5, 100, 0., 1.);
  fHistos->CreateTH1("hQAMatchingDRAbs", "Distance between part. level jet and  det. level jet", 100, 0., 1.);
  fHistos->CreateTH1("hQAMatchingDRel", "Distance between part. level jet and  det. level jet", 100, 0., 1.);
  fHistos->CreateTH1("hFracPtHardPart", "Part. level jet Pt relative to the Pt-hard of the event", 100, 0., 10.);
  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);

  fSampleSplitter = new TRandom;

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalJetEnergyScale::CheckMCOutliers() {
  if(!fMCRejectFilter) return true;
  if(!(fIsPythia || fIsHerwig)) return true;    // Only relevant for pt-hard production
  AliDebugStream(1) << "Using custom MC outlier rejection" << std::endl;
  auto partjets = GetJetContainer(fNameParticleJets);
  if(!partjets) return true;

  // Check whether there is at least one particle level jet with pt above n * event pt-hard
  auto jetiter = partjets->accepted();
  auto max = std::max_element(jetiter.begin(), jetiter.end(), [](const AliEmcalJet *lhs, const AliEmcalJet *rhs ) { return lhs->Pt() < rhs->Pt(); });
  if(max != jetiter.end())  {
    // At least one jet found with pt > n * pt-hard
    AliDebugStream(1) << "Found max jet with pt " << (*max)->Pt() << " GeV/c" << std::endl;
    if((*max)->Pt() > fPtHardAndJetPtFactor * fPtHard) return false;
  }
  return true;
}

Bool_t AliAnalysisTaskEmcalJetEnergyScale::Run(){
  AliDebugStream(1) << "Next event" << std::endl;
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;
  if(IsSelectEmcalTriggers(fTriggerSelectionString.Data())){
    auto mctrigger = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer));
    AliDebugStream(1) << "Found trigger decision object: " << (mctrigger ? "yes" : "no") << std::endl;
    if(!mctrigger){
      AliErrorStream() <<  "Trigger decision container with name " << fNameTriggerDecisionContainer << " not found in event - not possible to select EMCAL triggers" << std::endl;
      return false;
    }
    if(!mctrigger->IsEventSelected(fTriggerSelectionString)) return false;
  }
  AliDebugStream(1) << "event selected" << std::endl;
  fHistos->FillTH1("hEventCounter", 1);

  auto detjets = GetJetContainer(fNameDetectorJets),
       partjets = GetJetContainer(fNameParticleJets);
  if(!detjets || !partjets) {
    AliErrorStream() << "At least one jet container missing, exiting ..." << std::endl;
    return false;
  }
  AliClusterContainer *clusters(detjets->GetClusterContainer());
  AliTrackContainer *tracks(static_cast<AliTrackContainer *>(detjets->GetParticleContainer()));
  AliParticleContainer *particles(partjets->GetParticleContainer());
  AliDebugStream(1) << "Have both jet containers: part(" << partjets->GetNAcceptedJets() << "|" << partjets->GetNJets() << "), det(" << detjets->GetNAcceptedJets() << "|" << detjets->GetNJets() << ")" << std::endl;

  std::vector<AliEmcalJet *> acceptedjets;
  for(auto detjet : detjets->accepted()){
    AliDebugStream(2) << "Next jet" << std::endl;
    acceptedjets.push_back(detjet);
    auto partjet = detjet->ClosestJet();
    if(!partjet) {
      AliDebugStream(2) << "No tagged jet" << std::endl;
      continue;
    }
    bool isClosure = fSampleSplitter->Uniform() < fFractionResponseClosure;
    Double_t detpt = detjet->Pt();
    if(TMath::Abs(fScaleShift) > DBL_EPSILON){
      detpt += fScaleShift * detpt;
    }
    if(fFillHSparse) {
      Bool_t acceptancematch = false;
      if (partjet->GetJetAcceptanceType() & detjets->GetAcceptanceType()) acceptancematch = true;
      TVector3 basevec, tagvec;
      basevec.SetPtEtaPhi(detjet->Pt(), detjet->Eta(), detjet->Phi());
      tagvec.SetPtEtaPhi(partjet->Pt(), partjet->Eta(), partjet->Phi());
      double pointCorr[6] = {partjet->Pt(), detpt, detjet->NEF(), basevec.DeltaR(tagvec), acceptancematch ? 1. : 0.,  isClosure ? 0. : 1.},
             pointDiff[3] = {partjet->Pt(), detjet->NEF(), (detpt-partjet->Pt())/partjet->Pt()};
      fHistos->FillTHnSparse("hPtDiff", pointDiff);
      fHistos->FillTHnSparse("hPtCorr", pointCorr);
    }
    fHistos->FillTH2("hJetResponseFine", detpt, partjet->Pt());
    fHistos->FillTH1("hJetEnergyScale", partjet->Pt(), (detpt - partjet->Pt())/partjet->Pt());
    fHistos->FillTH1("hJetEnergyScaleDet", detpt, (detpt - partjet->Pt())/partjet->Pt());
    // splitting for closure test
    if(isClosure) {
      fHistos->FillTH2("hJetResponseFineClosure", detpt, partjet->Pt());
    } else {
      fHistos->FillTH2("hJetResponseFineNoClosure", detpt, partjet->Pt());
    }

    // Fill QA histograms
    fHistos->FillTH2("hQANEFPtPart", partjet->Pt(), partjet->NEF());
    fHistos->FillTH2("hQANEFPtDet", detjet->Pt(), detjet->NEF());
    fHistos->FillTH2("hQAEtaPhiPart", partjet->Eta(), TVector2::Phi_0_2pi(partjet->Phi()));
    fHistos->FillTH2("hQAEtaPhiDet", detjet->Eta(), TVector2::Phi_0_2pi(detjet->Phi()));
    auto deltaR = TMath::Abs(partjet->DeltaR(detjet));
    fHistos->FillTH1("hQAMatchingDRAbs", deltaR);
    fHistos->FillTH1("hQAMatchingDRel", deltaR/partjets->GetJetRadius());
    TVector3 jetvecDet(detjet->Px(), detjet->Py(), detjet->Px());
    if(clusters){
      auto leadcluster = detjet->GetLeadingCluster(clusters->GetArray());
      fHistos->FillTH2("hQANnePtDet", detjet->Pt(), detjet->GetNumberOfClusters());

      for(auto iclust = 0; iclust < detjet->GetNumberOfClusters(); iclust++) {
        auto cluster = detjet->Cluster(iclust);
        TLorentzVector clustervec;
        cluster->GetMomentum(clustervec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
        fHistos->FillTH2("hQAConstPtNeDet", detjet->Pt(), clustervec.Pt());
        fHistos->FillTH2("hQAEtaPhiConstNeDet", clustervec.Eta(), TVector2::Phi_0_2pi(clustervec.Phi()));
        fHistos->FillTH2("hQADeltaRNeutralDet", detjet->Pt(), jetvecDet.DeltaR(clustervec.Vect()));
      }

      if(leadcluster){
        TLorentzVector ptvec;
        leadcluster->GetMomentum(ptvec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
        fHistos->FillTH2("hQAZnePtDet", detjet->Pt(), detjet->GetZ(ptvec.Px(), ptvec.Py(), ptvec.Pz()));
        fHistos->FillTH2("hQAConstPtNeMaxDet", detjet->Pt(), ptvec.Pt());
        fHistos->FillTH2("hQAEtaPhiConstMaxNeDet", ptvec.Eta(), TVector2::Phi_0_2pi(ptvec.Phi()));
        fHistos->FillTH2("hQADeltaRMaxNeutralDet", detjet->Pt(), jetvecDet.DeltaR(ptvec.Vect()));
      }
    }
    if(tracks){
      fHistos->FillTH2("hQANChPtDet", detjet->Pt(),  detjet->GetNumberOfTracks());
      auto leadingtrack = detjet->GetLeadingTrack(tracks->GetArray());

      for(int itrk = 0; itrk < detjet->GetNumberOfTracks(); itrk++) {
        auto trk = detjet->Track(itrk);
        fHistos->FillTH2("hQAConstPtChPart", detjet->Pt(), trk->Pt());
        fHistos->FillTH2("hQAEtaPhiConstChDet", trk->Eta(), TVector2::Phi_0_2pi(trk->Phi()));
        fHistos->FillTH2("hQADeltaRChargedDet", detjet->Pt(), detjet->DeltaR(trk));
      }
      
      if(leadingtrack){
        fHistos->FillTH2("hQAZchPtDet", detjet->Pt(), detjet->GetZ(leadingtrack->Px(), leadingtrack->Py(), leadingtrack->Pz()));
        fHistos->FillTH2("hQAConstPtChMaxDet", detjet->Pt(), leadingtrack->Pt());
        fHistos->FillTH2("hQAEtaPhiConstMaxChDet", leadingtrack->Eta(), leadingtrack->Phi());
        fHistos->FillTH2("hQADeltaRMaxChargedDet", detjet->Pt(), detjet->DeltaR(leadingtrack));
      }
    }
    if(particles){
      AliVParticle *leadingcharged(nullptr), *leadingneutral(nullptr);
      int ncharged(0), nneutral(0);
      for(int ipart = 0; ipart < partjet->GetNumberOfTracks(); ipart++) {
        auto particle = partjet->Track(ipart);
        if(particle->Charge()) {
          ncharged++;
          fHistos->FillTH2("hQAConstPtChPart", partjet->Pt(), particle->Pt());
          fHistos->FillTH2("hQAEtaPhiConstChPart", particle->Eta(), TVector2::Phi_0_2pi(particle->Phi()));
          fHistos->FillTH2("hQADeltaRChargedPart", partjet->Pt(), partjet->DeltaR(particle));
          if(!leadingcharged) leadingcharged = particle;
          else {
            if(particle->E() > leadingcharged->E()) leadingcharged = particle;
          }
        } else {
          nneutral++;
          fHistos->FillTH2("hQAConstPtNePart", partjet->Pt(), particle->Pt());
          fHistos->FillTH2("hQAEtaPhiConstNePart", particle->Eta(), TVector2::Phi_0_2pi(particle->Phi()));
          fHistos->FillTH2("hQADeltaRNeutralPart", partjet->Pt(), partjet->DeltaR(particle));
          if(!leadingneutral) leadingneutral = particle;
          else {
            if(particle->E() > leadingneutral->E()) leadingneutral = particle;
          }
        }
      }
      if(leadingcharged) {
        fHistos->FillTH2("hQAConstPtChMaxPart", partjet->Pt(), leadingcharged->Pt());
        fHistos->FillTH2("hQAEtaPhiConstMaxChPart", leadingcharged->Eta(), TVector2::Phi_0_2pi(leadingcharged->Phi()));
        fHistos->FillTH2("hQAZchPtPart", partjet->Pt(), partjet->GetZ(leadingcharged));
        fHistos->FillTH2("hQADeltaRMaxChargedPart", partjet->Pt(), partjet->DeltaR(leadingcharged));
      }
      if(leadingneutral) {
        fHistos->FillTH2("hQAConstPtNeMaxPart", partjet->Pt(), leadingneutral->Pt());
        fHistos->FillTH2("hQAEtaPhiConstMaxNePart", leadingneutral->Eta(), TVector2::Phi_0_2pi(leadingneutral->Phi()));
        fHistos->FillTH2("hQAZnePtPart", partjet->Pt(), partjet->GetZ(leadingneutral));
        fHistos->FillTH2("hQADeltaRMaxNeutralPart", partjet->Pt(), partjet->DeltaR(leadingneutral));
      }
      fHistos->FillTH2("hQANChPtPart", partjet->Pt(), ncharged);
      fHistos->FillTH2("hQANnePtPart", partjet->Pt(), nneutral);
    }
    fHistos->FillTH2("hQAJetAreaVsJetPtPart", partjet->Pt(), partjet->Area());
    fHistos->FillTH2("hQAJetAreaVsJetPtDet", detjet->Pt(), detjet->Area());
    fHistos->FillTH2("hQAJetAreaVsNEFPart", partjet->NEF(), partjet->Area());
    fHistos->FillTH2("hQAJetAreaVsNEFDet", detjet->NEF(), detjet->Area());
    fHistos->FillTH2("hQAJetAreaVsNConstPart", partjet->GetNumberOfTracks(), partjet->Area());
    fHistos->FillTH2("hQAJetAreaVsNConstDet", detjet->GetNumberOfClusters() + detjet->GetNumberOfTracks(), detjet->Area());
    fHistos->FillTH1("hFracPtHardPart", partjet->Pt()/fPtHard);
  }

  // efficiency x acceptance: Add histos for all accepted and reconstucted accepted jets
  for(auto partjet : partjets->accepted()){
    if(fFillHSparse){
      auto detjet = partjet->ClosestJet();
      double effvec[3] = {partjet->Pt(), 0., 0.};
      if(detjet) {
        // Found a match
        effvec[1] = detjet->Pt();
        if(TMath::Abs(fScaleShift) > DBL_EPSILON){
          effvec[1] += fScaleShift * effvec[1];
        }
        effvec[2] = 1;    // Tagged
        if(std::find(acceptedjets.begin(), acceptedjets.end(), detjet) != acceptedjets.end()) effvec[2] = 2;
      }
      fHistos->FillTHnSparse("hJetfindingEfficiency", effvec);
    }
    fHistos->FillTH1("hJetSpectrumPartAll", partjet->Pt());
  }

  return true;
}

bool AliAnalysisTaskEmcalJetEnergyScale::IsSelectEmcalTriggers(const TString &triggerstring) const {
  const std::array<TString, 8> kEMCALTriggers = {
    "EJ1", "EJ2", "DJ1", "DJ2", "EG1", "EG2", "DG1", "DG2"
  };
  bool isEMCAL = false;
  for(auto emcaltrg : kEMCALTriggers) {
    if(triggerstring.Contains(emcaltrg)) {
      isEMCAL = true;
      break;
    }
  }
  return isEMCAL;
}

AliAnalysisTaskEmcalJetEnergyScale *AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t jetradius, Bool_t useDCAL, const char *namepartcont, const char *trigger, const char *suffix) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("EmcalTriggerJets::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale", "No analysis manager available");
    return nullptr;
  } 

  auto inputhandler = mgr->GetInputEventHandler();
  auto isAOD = inputhandler->IsA() == AliAODInputHandler::Class();

  std::string jettypename;
  AliJetContainer::JetAcceptanceType acceptance(AliJetContainer::kTPCfid);
  AliJetContainer::EJetType_t mcjettype(jettype);
  bool addClusterContainer(false), addTrackContainer(false);
  switch(jettype){
    case AliJetContainer::kFullJet:
        jettypename = "FullJet";
        acceptance = useDCAL ? AliJetContainer::kDCALfid : AliJetContainer::kEMCALfid;
        addClusterContainer = addTrackContainer = true;
        break;
    case AliJetContainer::kChargedJet:
        jettypename = "ChargedJet";
        acceptance = AliJetContainer::kTPCfid;
        addTrackContainer = true;
        mcjettype = AliJetContainer::kChargedJet; 
        break;
    case AliJetContainer::kNeutralJet:
        jettypename = "NeutralJet";
        acceptance = useDCAL ? AliJetContainer::kDCALfid : AliJetContainer::kEMCALfid;
        addClusterContainer = true;
        mcjettype = AliJetContainer::kFullJet;    // Correct back neutral detector-level jets to full particle level jets
        break;
    case AliJetContainer::kUndefinedJetType:
        break;
  };

  std::stringstream taskname, tag;
  tag << jettypename << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  if(strlen(suffix)) tag << "_" << suffix;
  taskname << "EnergyScaleTask_" << tag.str();
  AliAnalysisTaskEmcalJetEnergyScale *energyscaletask = new AliAnalysisTaskEmcalJetEnergyScale(taskname.str().data());
  mgr->AddTask(energyscaletask);
  energyscaletask->SetTriggerName(trigger);

  TString partcontname(namepartcont);
  if(partcontname == "usedefault") partcontname = "mcparticles";
  auto partcont = energyscaletask->AddMCParticleContainer(partcontname.Data());
  partcont->SetMinPt(0.);

  AliClusterContainer *clusters(nullptr);
  if(addClusterContainer) {
    clusters = energyscaletask->AddClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    clusters->SetDefaultClusterEnergy(energydef);
    clusters->SetClusUserDefEnergyCut(energydef, 0.3);
  }
  AliTrackContainer *tracks(nullptr);
  if(addTrackContainer) {
    tracks = energyscaletask->AddTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
  }

  auto contpartjet = energyscaletask->AddJetContainer(mcjettype, AliJetContainer::antikt_algorithm, recoscheme, jetradius,
                                                      acceptance, partcont, nullptr);
  contpartjet->SetName("particleLevelJets");
  energyscaletask->SetNamePartJetContainer("particleLevelJets");
  std::cout << "Adding particle-level jet container with underling array: " << contpartjet->GetArrayName() << std::endl;

  auto contdetjet = energyscaletask->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, recoscheme, jetradius,
                                                     acceptance, tracks, clusters);
  contdetjet->SetName("detectorLevelJets");
  energyscaletask->SetNameDetJetContainer("detectorLevelJets");
  std::cout << "Adding detector-level jet container with underling array: " << contdetjet->GetArrayName() << std::endl;

  std::stringstream outnamebuilder, listnamebuilder;
  listnamebuilder << "EnergyScaleHists_" << tag.str();
  outnamebuilder << mgr->GetCommonFileName() << ":EnergyScaleResults_" << tag.str();

  mgr->ConnectInput(energyscaletask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(energyscaletask, 1, mgr->CreateContainer(listnamebuilder.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, outnamebuilder.str().data()));
  return energyscaletask;
} 