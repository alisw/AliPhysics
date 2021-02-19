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
#include <memory>
#include <sstream>

#include <TArray.h>
#include <TCustomBinning.h>
#include <THistManager.h>
#include <TLinearBinning.h>

#include "AliAnalysisTaskEmcalSoftDropData.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliClusterContainer.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliTrackContainer.h"
#include "AliVCluster.h"
#include "AliVTrack.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalSoftDropData)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalSoftDropData::AliAnalysisTaskEmcalSoftDropData() : 
  AliAnalysisTaskEmcalJet(),
  AliAnalysisEmcalSoftdropHelperImpl(),
  AliAnalysisEmcalTriggerSelectionHelperImpl(),
  fTriggerBits(AliVEvent::kAny),
  fTriggerString(""),
  fUseDownscaleWeight(false),
  fBeta(0.),
  fZcut(0.1),
  fReclusterizer(kCAAlgo),
  fUseChargedConstituents(kTRUE),
  fUseNeutralConstituents(kTRUE), 
  fDropMass0Jets(false),
  fHistos(nullptr),
  fPtBinning(nullptr)
{

}

AliAnalysisTaskEmcalSoftDropData::AliAnalysisTaskEmcalSoftDropData(EMCAL_STRINGVIEW name) : 
  AliAnalysisTaskEmcalJet(name.data(), kTRUE),
  AliAnalysisEmcalSoftdropHelperImpl(),
  AliAnalysisEmcalTriggerSelectionHelperImpl(),
  fTriggerBits(AliVEvent::kAny),
  fTriggerString(""),
  fUseDownscaleWeight(false),
  fBeta(0.),
  fZcut(0.1),
  fReclusterizer(kCAAlgo),
  fUseChargedConstituents(kTRUE),
  fUseNeutralConstituents(kTRUE),
  fDropMass0Jets(false),
  fHistos(nullptr),
  fPtBinning(nullptr)
{

}

AliAnalysisTaskEmcalSoftDropData::~AliAnalysisTaskEmcalSoftDropData() {
  if(fPtBinning) delete fPtBinning;
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalSoftDropData::UserCreateOutputObjects() {
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  double R = double(int(GetJetContainer("datajets")->GetJetRadius() * 1000.))/1000.;   // Save cast from float to double truncating after 3rd decimal digit
  if(!fPtBinning) fPtBinning = new TLinearBinning(300, 0., 300.);  // Use fine binning for data, rebin offline
  std::unique_ptr<TBinning> zgBinning(GetZgBinning(fZcut)),
                            rgBinning(GetRgBinning(R)),
                            nsdBinning(new TLinearBinning(22, -1.5, 20.5)),     // Negative bins for untagged jets
                            thetagBinning(new TLinearBinning(11, -0.1, 1.)),
                            triggerClusterBinning(new TLinearBinning(kTrgClusterN, -0.5, kTrgClusterN - 0.5));
  TArrayD edgesPt;
  fPtBinning->CreateBinEdges(edgesPt);

  // Define binnings for substructure THnSparses
  const TBinning *binningSparseZg[3] = {zgBinning.get(), fPtBinning, triggerClusterBinning.get()},
                 *binningSparseRg[3] = {rgBinning.get(), fPtBinning, triggerClusterBinning.get()},
                 *binningSparseThethag[3] = {thetagBinning.get(), fPtBinning, triggerClusterBinning.get()},
                 *binningSparseNsd[3] = {nsdBinning.get(), fPtBinning, triggerClusterBinning.get()};

  fHistos = new THistManager("histosSoftdrop");
  fHistos->CreateTH1("hEventCounter", "EventCounter", 1, 0.5, 1.5);
  fHistos->CreateTH1("hEventCounterRun", "Runwise event counter", 100000, 200000, 300000);
  fHistos->CreateTH1("hClusterCounterAbs", "Event counter histogram absolute", kTrgClusterN, -0.5, kTrgClusterN - 0.5);
  fHistos->CreateTH2("hJetPtRaw", "raw jet pt", kTrgClusterN, -0.5, kTrgClusterN - 0.5, 300, 0., 300.);
  fHistos->CreateTHnSparse("hZgVsPt", "zg vs pt", 3, binningSparseZg, "s");
  fHistos->CreateTHnSparse("hRgVsPt", "rg vs pt", 3, binningSparseRg, "s");
  fHistos->CreateTHnSparse("hNsdVsPt", "nsd vs pt", 3, binningSparseNsd, "s");
  fHistos->CreateTHnSparse("hThetagVsPt", "thetag vs pt", 3, binningSparseThethag, "s");
  fHistos->CreateTH2("hSkippedJets", "Number of skipped jets", *triggerClusterBinning, *fPtBinning);
  if(fUseDownscaleWeight){
    fHistos->CreateTH1("hClusterCounter", "Event counter histogram", kTrgClusterN, -0.5, kTrgClusterN - 0.5);
    fHistos->CreateTHnSparse("hZgVsPtWeighted", "zg vs pt (weighted)", 3, binningSparseZg, "s");
    fHistos->CreateTHnSparse("hRgVsPtWeighted", "rg vs pt (weighted)", 3, binningSparseRg, "s");
    fHistos->CreateTHnSparse("hNsdVsPtWeighted", "nsd vs pt (weighted)", 3, binningSparseNsd, "s");
    fHistos->CreateTHnSparse("hThetagVsPtWeighted", "thetag vs pt (weighted)", 3, binningSparseThethag, "s");
    fHistos->CreateTH1("hEventCounterWeighted", "Event counter, weighted", 1., 0.5, 1.5);
    fHistos->CreateTH1("hEventCounterRunWeighted", "Runwise event counter (weighted)", 100000, 200000, 300000);
    fHistos->CreateTProfile("hDownscaleFactorsRunwise", "Runwise downscale factors", 100000, 200000, 300000);
    fHistos->CreateTH2("hJetPtRawWeighted", "raw jet pt", kTrgClusterN, -0.5, kTrgClusterN - 0.5, 300, 0., 300., "s");
    fHistos->CreateTH2("hSkippedJetsWeighted", "Number of skipped jets (weighted)", *triggerClusterBinning, *fPtBinning);
  }

  // A bit of QA stuff
  fHistos->CreateTH2("hQANEFPt", "Neutral energy fraction; p_{t} (GeV/c); NEF", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAEtaPhi", "#eta vs. #phi for selected jets; #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAZchPt", "z_{ch,max}; p_{t} (GeV/c); z_{ch,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZnePt", "z_{ne,max}; p_{t} (GeV/c); z_{ne,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQANChPt", "Number of charged constituents; p_{t} (GeV/c); N_{ch}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQANnePt", "Number of neutral constituents; p_{t} (GeV/c); N_{ne}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQAJetAreaVsJetPt", "Jet area vs. jet pt at detector level; p_{t} (GeV/c); Area", 350, 0., 350., 200, 0., 2.);
  fHistos->CreateTH2("hQAJetAreaVsNEF", "Jet area vs. NEF at detector level; NEF; Area", 100, 0., 1., 200, 0., 2.);
  fHistos->CreateTH2("hQAJetAreaVsNConst", "Jet area vs. number of consituents at detector level; Number of constituents; Area", 101, -0.5, 100.5, 200, 0., 2.);
  fHistos->CreateTH2("hSDUsedChargedPtjvPtc", "p_{t,j} vs. p_{t,const} for tracks used in SD; p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350., 0., 350., 350, 0., 350.);
  fHistos->CreateTH2("hSDUsedNeutralPtjvPtc", "p_{t,j} vs. p_{t,const} for clusters used in SD; p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350., 0., 350., 350, 0., 350.);
  fHistos->CreateTH2("hSDUsedChargedPtjvPtcMax", "p_{t,j} vs. p_{t,const} for max tracks used in SD; p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hSDUsedNeutralPtjvPcMax", "p_{t,j} vs. p_{t,const} for max clusters used in SD; p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hSDUsedChargedEtaPhi", "#eta-phi for tracks used in SD; #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hSDUsedNeutralEtaPhi", "#eta vs. #phi for clusters used in SD; #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hSDUsedChargedEtaPhiMax", "#eta-phi for tracks used in SD; #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hSDUsedNeutralEtaPhiMax", "#eta vs. #phi for clusters used in SD; #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hSDUsedChargedDR", "#DeltaR vs. p_{t,jet} for tracks used in SD; p_{t,jet}; #DeltaR", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hSDUsedNeutralDR", "#DeltaR vs. p_{t,jet} for clusters used in SD; p_{t,jet}; #DeltaR", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hSDUsedChargedDRMax", "#DeltaR vs. p_{t,jet} for tracks used in SD; p_{t,jet}; #DeltaR", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hSDUsedNeutralDRMax", "#DeltaR vs. p_{t,jet} for clusters used in SD; p_{t,jet}; #DeltaR", 100, -1., 1., 100, 0., 7.);
  // Cluster constituent QA
  fHistos->CreateTH2("hSDUsedClusterTimeVsE", "Cluster time vs. energy; time (ns); E (GeV)", 1200, -600, 600, 200, 0., 200);
  fHistos->CreateTH2("hSDUsedClusterTimeVsEFine", "Cluster time vs. energy (main region); time (ns); E (GeV)", 1000, -100, 100, 200, 0., 200);
  fHistos->CreateTH2("hSDUsedClusterNCellVsE", "Cluster number of cells vs. energy; Number of cells; E (GeV)", 201, -0.5, 200.5, 200, 0., 200.);
  fHistos->CreateTH2("hSDUsedlusterM02VsE", "Cluster M02 vs energy; M02; E (GeV)", 150, 0., 1.5, 200, 0., 200.);
  fHistos->CreateTH2("hSDUsedClusterFracLeadingVsE", "Cluster frac leading cell vs energy; E (GeV); Frac. leading cell", 200, 0., 200., 110, 0., 1.1);
  fHistos->CreateTH2("hSDUsedClusterFracLeadingVsNcell", "Cluster frac leading cell vs number of cells; Number of cells; Frac. leading cell", 201, -0.5, 200.5, 110, 0., 1.1);


  for(auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalSoftDropData::IsTriggerSelected(){
  if(!(fInputHandler->IsEventSelected() & fTriggerBits)) return false;
  if(fTriggerString.length()) {
    if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerString)) return false;
  }
  return true;
}

void AliAnalysisTaskEmcalSoftDropData::RunChanged(Int_t newrun){
  if(fUseDownscaleWeight) 
    PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->SetRun(newrun);
}

Bool_t AliAnalysisTaskEmcalSoftDropData::Run() {
  auto jets = GetDetLevelJetContainer();
  if(!jets) {
    AliErrorStream() << "Jet container not found" << std::endl;
    return false;
  } 
  auto clusters = GetClusterContainer(0);
  if(fUseNeutralConstituents && !clusters) {
    AliErrorStream() << "Cluster container not found, but neutral constituents requested" << std::endl; 
  }
  auto tracks  = GetTrackContainer(0);
  if(fUseChargedConstituents &&  !tracks) {
    AliErrorStream() << "Track container not found, but charged constituent requested." << std::endl;
    return false;
  }

  auto trgclusters = GetTriggerClusterIndices(fInputEvent->GetFiredTriggerClasses().Data());
  Double_t weight = fUseDownscaleWeight ? 1./GetDownscaleWeight() : 1.;
  Double_t Rjet = jets->GetJetRadius();
  fHistos->FillTH1("hEventCounter", 1.);
  fHistos->FillTH1("hEventCounterRun", fRunNumber);
  for(auto icl : trgclusters) fHistos->FillTH1("hClusterCounterAbs", icl);
  if(fUseDownscaleWeight) {
    fHistos->FillTH1("hEventCounterWeighted", 1., weight);
    fHistos->FillTH1("hEventCounterRunWeighted", fRunNumber, weight);
    fHistos->FillProfile("hDownscaleFactorsRunwise", fRunNumber, 1./weight);
    for(auto icl : trgclusters) fHistos->FillTH1("hClusterCounter", icl, weight);
  }

  for(auto jet : jets->accepted()){
    AliDebugStream(2) << "Next accepted jet with pt " << jet->Pt() << std::endl;
    for(auto icl : trgclusters) {
      fHistos->FillTH2("hJetPtRaw", icl, jet->Pt());
      if(fUseDownscaleWeight) fHistos->FillTH2("hJetPtRawWeighted", icl, jet->Pt(), weight);
    }
    try {
      FillJetQA(*jet, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      auto sdparams = MakeSoftdrop(*jet, jets->GetJetRadius(), false, {(AliAnalysisEmcalSoftdropHelperImpl::EReclusterizer_t)fReclusterizer, fBeta, fZcut, fUseChargedConstituents, fUseNeutralConstituents}, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets);
      auto splittings = IterativeDecluster(*jet, jets->GetJetRadius(), false, {(AliAnalysisEmcalSoftdropHelperImpl::EReclusterizer_t)fReclusterizer, fBeta, fZcut, fUseChargedConstituents, fUseNeutralConstituents}, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets);
      bool untagged = sdparams.fZg < fZcut;
      AliDebugStream(2) << "Found jet with pt " << jet->Pt() << " and zg " << sdparams.fZg << std::endl;
      Double_t pointZg[3] = {sdparams.fZg, jet->Pt(), -1},
               pointRg[3] = {untagged ? -0.01 : sdparams.fRg, jet->Pt(), -1},
               pointThetaG[3] = {untagged ? -0.05 : sdparams.fRg/Rjet, jet->Pt(), -1},
               pointNSD[3] = {untagged ? -1. : double(splittings.size()), jet->Pt(), -1};
      for(auto icl : trgclusters) {
        pointZg[2] = pointRg[2] = pointNSD[2] = pointThetaG[2] = icl;
        fHistos->FillTHnSparse("hZgVsPt", pointZg);
        fHistos->FillTHnSparse("hRgVsPt", pointRg);
        fHistos->FillTHnSparse("hNsdVsPt", pointNSD);
        fHistos->FillTHnSparse("hThetagVsPt", pointThetaG);
        if(fUseDownscaleWeight) {
          fHistos->FillTHnSparse("hZgVsPtWeighted", pointZg, weight);
          fHistos->FillTHnSparse("hRgVsPtWeighted", pointRg, weight);
          fHistos->FillTHnSparse("hNsdVsPtWeighted", pointNSD, weight);
          fHistos->FillTHnSparse("hThetagVsPtWeighted", pointThetaG, weight);
        } 
      }
    } catch (int e) {
      if(fUseChargedConstituents && fUseNeutralConstituents) AliErrorStream() << "Softdrop error " << e << ": Having 0 constituents for reclustering" << std::endl;
      for(auto icl : trgclusters) {
        fHistos->FillTH2("hSkippedJets", icl, jet->Pt());
        if(fUseDownscaleWeight) 
          fHistos->FillTH2("hSkippedJetsWeighted", icl, jet->Pt(), weight);
      }
    }

    // Fill QA plots - trigger cluster independent
    // Those plots have been in before (as part of the Tree) but were 
    // removed in order to reduce the disk space consumption.
    fHistos->FillTH2("hQANEFPt", jet->Pt(), jet->NEF(), weight);
    fHistos->FillTH2("hQAEtaPhi", jet->Eta(), jet->Phi(), weight);
    if(clusters){
      auto leadcluster = jet->GetLeadingCluster(clusters->GetArray());
      if(leadcluster){
        TLorentzVector ptvec;
        leadcluster->GetMomentum(ptvec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
        fHistos->FillTH2("hQAZnePt", jet->Pt(), jet->GetZ(ptvec.Px(), ptvec.Py(), ptvec.Pz()), weight);
      }
    }
    if(tracks){
      auto leadingtrack = jet->GetLeadingTrack(tracks->GetArray());
      if(leadingtrack) fHistos->FillTH2("hQAZchPt", jet->Pt(), jet->GetZ(leadingtrack->Px(), leadingtrack->Py(), leadingtrack->Pz()), weight);
    }
    fHistos->FillTH2("hQANChPt", jet->Pt(), jet->GetNumberOfTracks(), weight);
    fHistos->FillTH2("hQANnePt", jet->Pt(), jet->GetNumberOfClusters(), weight);
    fHistos->FillTH2("hQAJetAreaVsJetPt", jet->Pt(), jet->Area(), weight);
    fHistos->FillTH2("hQAJetAreaVsNEF", jet->NEF(), jet->Area(), weight);
    fHistos->FillTH2("hQAJetAreaVsNConst", jet->GetNumberOfClusters() + jet->GetNumberOfTracks(), jet->Area(), weight);
  }
  return true;
}

Double_t AliAnalysisTaskEmcalSoftDropData::GetDownscaleWeight() const {
  Double_t weight = 1.;
  TString triggerclass;
  if(fTriggerString == "INT7") triggerclass = "CINT7-B-NOPF-CENT";
  else if(fTriggerString == "EJ1") triggerclass = "CEMC7EJ1-B-NOPF-CENTNOTRD";
  else if(fTriggerString == "EJ2") triggerclass = "CEMC7EJ2-B-NOPF-CENT";
  if(triggerclass.Length()) weight = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(triggerclass);
  return weight;
}

void AliAnalysisTaskEmcalSoftDropData::FillJetQA(const AliEmcalJet &jet, AliVCluster::VCluUserDefEnergy_t energydef) {
  TVector3 maxcharged, maxneutral, jetvec(jet.Px(), jet.Py(), jet.Pz());
  bool hasMaxCharged = false, hasMaxNeutral = false;
  if(fUseChargedConstituents){                    // Neutral particles part of particle container in case of MC
    AliDebugStream(1) << "Jet substructure: Using charged constituents" << std::endl;
    for(int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++){
      auto track = jet.Track(itrk);
      if(!track->Charge() && !fUseNeutralConstituents) continue;      // Reject neutral constituents in case of using only charged consituents
      if(track->Charge() && !fUseChargedConstituents) continue;       // Reject charged constituents in case of using only neutral consituents
      TVector3 trackvec(track->Px(), track->Py(), track->Pz());
      fHistos->FillTH2("hSDUsedChargedPtjvPtc", jet.Pt(), TMath::Abs(track->Pt()));
      fHistos->FillTH2("hSDUsedChargedEtaPhi", track->Eta(), TVector2::Phi_0_2pi(track->Phi()));
      fHistos->FillTH2("hSDUsedChargedDR", jet.Pt(), jetvec.DeltaR(trackvec));
      if(!hasMaxCharged) {
        maxcharged = trackvec;
        hasMaxCharged = true;
      } else {
        if(trackvec.Pt() > maxcharged.Pt())
          maxcharged = trackvec;
      }
    }
  }
  if(hasMaxCharged) {
    fHistos->FillTH2("hSDUsedChargedPtjvPtcMax", jet.Pt(), TMath::Abs(maxcharged.Pt()));
    fHistos->FillTH2("hSDUsedChargedEtaPhiMax", maxcharged.Eta(), TVector2::Phi_0_2pi(maxcharged.Phi()));
    fHistos->FillTH2("hSDUsedChargedDRMax", jet.Pt(), jetvec.DeltaR(maxcharged));
  }

  if(fUseNeutralConstituents){
    AliDebugStream(1) << "Jet substructure: Using neutral constituents" << std::endl;
    for(int icl = 0; icl < jet.GetNumberOfClusters(); icl++) {
      auto cluster = jet.Cluster(icl);
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, fVertex, energydef);
      TVector3 clustervec3(clustervec.Px(), clustervec.Py(), clustervec.Pz());
      fHistos->FillTH2("hSDUsedNeutralPtjvPtc", jet.Pt(), clustervec.Pt());
      fHistos->FillTH2("hSDUsedNeutralEtaPhi", clustervec.Eta(), TVector2::Phi_0_2pi(clustervec.Phi()));
      fHistos->FillTH2("hSDUsedNeutralDR", jet.Pt(), jetvec.DeltaR(clustervec3));
      fHistos->FillTH2("hSDUsedClusterTimeVsE", cluster->GetTOF() * 1e9, clustervec.E());
      fHistos->FillTH2("hSDUsedClusterTimeVsEFine", cluster->GetTOF() * 1e9, clustervec.E());
      fHistos->FillTH2("hSDUsedClusterNCellVsE", cluster->GetNCells(), clustervec.E());
      fHistos->FillTH2("hSDUsedlusterM02VsE", cluster->GetM02(), clustervec.E());
      double maxamplitude = 0.;
      for(int icell = 0; icell < cluster->GetNCells(); icell++) {
        double amplitude = fInputEvent->GetEMCALCells()->GetAmplitude(fInputEvent->GetEMCALCells()->GetCellPosition(cluster->GetCellAbsId(icell)));
        if(amplitude > maxamplitude) maxamplitude = amplitude;
      }
      fHistos->FillTH2("hSDUsedClusterFracLeadingVsE", clustervec.E(), maxamplitude/cluster->E());
      fHistos->FillTH2("hSDUsedClusterFracLeadingVsNcell", cluster->GetNCells(), maxamplitude/cluster->E());
      if(!hasMaxNeutral) {
        maxneutral = clustervec3;
        hasMaxNeutral = true;
      } else {
        if(clustervec3.Pt() > maxneutral.Pt())
          maxneutral = clustervec3;
      }
    }
  }
  if(hasMaxNeutral){
    fHistos->FillTH2("hSDUsedNeutralPtjvPcMax", jet.Pt(), maxneutral.Pt());
    fHistos->FillTH2("hSDUsedNeutralEtaPhiMax", maxneutral.Eta(), TVector2::Phi_0_2pi(maxneutral.Phi()));
    fHistos->FillTH2("hSDUsedNeutralDRMax", jet.Pt(), jetvec.DeltaR(maxneutral));
  }
}

void AliAnalysisTaskEmcalSoftDropData::ConfigureDetJetSelection(Double_t minJetPt, Double_t maxTrackPt, Double_t maxClusterPt, Double_t minAreaPerc) {
  auto detjets = GetDetLevelJetContainer();
  
  detjets->SetJetPtCut(minJetPt);
  if(detjets->GetJetType() == AliJetContainer::kFullJet || detjets->GetJetType() == AliJetContainer::kChargedJet) {
    detjets->SetMaxTrackPt(maxTrackPt);
  }
  if(detjets->GetJetType() == AliJetContainer::kFullJet || detjets->GetJetType() == AliJetContainer::kNeutralJet) {
    detjets->SetMaxClusterPt(maxClusterPt);
  }
  if(minAreaPerc >= 0.) {
    detjets->SetPercAreaCut(minAreaPerc);
  }
}

AliAnalysisTaskEmcalSoftDropData *AliAnalysisTaskEmcalSoftDropData::AddTaskEmcalSoftDropData(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, AliVCluster::VCluUserDefEnergy_t energydef, EMCAL_STRINGVIEW trigger) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  Bool_t isAOD(kFALSE);
  AliInputEventHandler *inputhandler = static_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
  if(inputhandler) {
    if(inputhandler->IsA() == AliAODInputHandler::Class()){
      std::cout << "Analysing AOD events\n";
      isAOD = kTRUE;
    } else {
      std::cout << "Analysing ESD events\n";
    }
  }

  std::stringstream taskname;
  taskname << "SoftdropDataMaker_R" << std::setw(2) << std::setfill('0') << int(jetradius*10) << trigger;  
  AliAnalysisTaskEmcalSoftDropData *datamaker = new AliAnalysisTaskEmcalSoftDropData(taskname.str().data());
  datamaker->SelectCollisionCandidates(AliVEvent::kAny);
  mgr->AddTask(datamaker);

  AliTrackContainer *tracks(nullptr);
  if((jettype == AliJetContainer::kChargedJet) || (jettype == AliJetContainer::kFullJet)){
      tracks = datamaker->AddTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
      std::cout << "Track container name: " << tracks->GetName() << std::endl;
      tracks->SetMinPt(0.15);
  }
  AliClusterContainer *clusters(nullptr);
  if((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)){
    std::cout << "Using full or neutral jets ..." << std::endl;
    clusters = datamaker->AddClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    std::cout << "Cluster container name: " << clusters->GetName() << std::endl;
    switch (energydef)
    {
    case AliVCluster::VCluUserDefEnergy_t::kHadCorr:
      clusters->SetClusHadCorrEnergyCut(0.3); // 300 MeV E-cut
      clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      break;
    case AliVCluster::VCluUserDefEnergy_t::kNonLinCorr:
      clusters->SetClusNonLinCorrEnergyCut(0.3); // 300 MeV E-cut
      clusters->SetDefaultClusterEnergy(AliVCluster::kNonLinCorr); 
      break;
    case AliVCluster::VCluUserDefEnergy_t::kUserDefEnergy1:
    case AliVCluster::VCluUserDefEnergy_t::kUserDefEnergy2:
    case AliVCluster::VCluUserDefEnergy_t::kLastUserDefEnergy:
    default:
      AliErrorGeneralStream("AliAnalysisTaskEmcalSoftDropData::AddTaskEmcalSoftDropData") << "Requested energy type not supported" << std::endl;
      break;
    }

  } else {
    std::cout << "Using charged jets ... " << std::endl;
  }

  AliJetContainer *datajets = datamaker->AddJetContainer(
                              jettype,
                              AliJetContainer::antikt_algorithm,
                              recombinationScheme,
                              jetradius,
                              ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)) ? AliEmcalJet::kEMCALfid : AliEmcalJet::kTPCfid,
                              tracks, clusters);
  datajets->SetName("datajets");
  datajets->SetJetPtCut(0.);
  datajets->SetMaxTrackPt(1000.);

  std::string jettypestring;
  switch(jettype) {
    case AliJetContainer::kFullJet: jettypestring = "FullJets"; break;
    case AliJetContainer::kChargedJet: jettypestring = "ChargedJets"; break;
    case AliJetContainer::kNeutralJet: jettypestring = "NeutralJets"; break;
    default: jettypestring = "Undef";
  };

  ULong_t triggerbits(AliVEvent::kINT7);
  std::string triggerstring(trigger);
  std::cout << "Found trigger " << triggerstring << std::endl;
  if(triggerstring == "EJ1") {
    std::cout << "Setting binning mode for EJ1" << std::endl;
    triggerbits = AliVEvent::kEMCEJE;
  } else if(triggerstring == "EJ2") {
    std::cout << "Setting binning mode for EJ2" << std::endl;
    triggerbits = AliVEvent::kEMCEJE;
  }
  datamaker->SetSelectTrigger(triggerbits, trigger.data());

  // Connecting containers
  std::stringstream outputfile, histname;
  outputfile << mgr->GetCommonFileName() << ":SoftDropResponse_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  histname << "SoftDropResponseHistos_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  mgr->ConnectInput(datamaker, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(datamaker, 1, mgr->CreateContainer(histname.str().data(), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.str().data()));

  return datamaker;
}