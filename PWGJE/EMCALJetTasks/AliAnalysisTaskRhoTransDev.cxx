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
#include "AliAnalysisTaskRhoTransDev.h"

#include <TClonesArray.h>
#include <TMath.h>

#include <AliLog.h>
#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"

ClassImp(AliAnalysisTaskRhoTransDev);

AliAnalysisTaskRhoTransDev::AliAnalysisTaskRhoTransDev() :
  AliAnalysisTaskRhoBaseDev(),
  fHistB2BRhoVsCent(0),
  fHistB2BRhoVsLeadJetPt(),
  fHistB2BRhoVsLeadTrackPt(0),
  fHistB2BRhoVsNtrack(0),
  fHistB2BRhoVsLeadClusterE(0),
  fHistB2BRhoVsNcluster(0),
  fHistB2BRhoScaledVsCent(0),
  fHistB2BRhoScaledVsNtrack(0),
  fHistB2BRhoScaledVsNcluster(0)
{
}

AliAnalysisTaskRhoTransDev::AliAnalysisTaskRhoTransDev(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBaseDev(name, histo),
  fHistB2BRhoVsCent(0),
  fHistB2BRhoVsLeadJetPt(),
  fHistB2BRhoVsLeadTrackPt(0),
  fHistB2BRhoVsNtrack(0),
  fHistB2BRhoVsLeadClusterE(0),
  fHistB2BRhoVsNcluster(0),
  fHistB2BRhoScaledVsCent(0),
  fHistB2BRhoScaledVsNtrack(0),
  fHistB2BRhoScaledVsNcluster(0)
{
}

void AliAnalysisTaskRhoTransDev::UserCreateOutputObjects()
{
  if (!fCreateHisto) return;

  AliAnalysisTaskRhoBaseDev::UserCreateOutputObjects();

  TString name;

  Int_t maxTracks = 6000;
  Double_t maxRho = 500;
  Int_t nRhoBins = 500;

  if (fForceBeamType == kpp) {
    maxRho = 50;
    maxTracks = 200;
  }
  else if (fForceBeamType == kpA) {
    maxRho = 200;
    maxTracks = 500;
  }

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);

  fHistB2BRhoVsCent = new TH2F("fHistB2BRhoVsCent", "fHistB2BRhoVsCent", 100, 0,  100, nRhoBins, 0, maxRho);
  fHistB2BRhoVsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistB2BRhoVsCent->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
  fOutput->Add(fHistB2BRhoVsCent);

  if (fParticleCollArray.size() > 0) {
    fHistB2BRhoVsNtrack = new TH2F("fHistB2BRhoVsNtrack", "fHistB2BRhoVsNtrack", 200, 0, maxTracks, nRhoBins, 0, maxRho);
    fHistB2BRhoVsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistB2BRhoVsNtrack->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistB2BRhoVsNtrack);

    fHistB2BRhoVsLeadTrackPt = new TH2F("fHistB2BRhoVsLeadTrackPt", "fHistB2BRhoVsLeadTrackPt", nPtBins, 0, fMaxPt, nRhoBins, 0, maxRho);
    fHistB2BRhoVsLeadTrackPt->GetXaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fHistB2BRhoVsLeadTrackPt->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistB2BRhoVsLeadTrackPt);
  }

  if (fClusterCollArray.size()>0) {
    fHistB2BRhoVsNcluster = new TH2F("fHistB2BRhoVsNcluster", "fHistB2BRhoVsNcluster", 50, 0, maxTracks / 4, nRhoBins, 0, maxRho);
    fHistB2BRhoVsNcluster->GetXaxis()->SetTitle("No. of clusters");
    fHistB2BRhoVsNcluster->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistB2BRhoVsNcluster);

    fHistB2BRhoVsLeadClusterE = new TH2F("fHistB2BRhoVsLeadClusterE", "fHistB2BRhoVsLeadClusterE", nPtBins, 0, fMaxPt, nRhoBins, 0, maxRho);
    fHistB2BRhoVsLeadClusterE->GetXaxis()->SetTitle("#it{p}_{T,track} (GeV/c)");
    fHistB2BRhoVsLeadClusterE->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistB2BRhoVsLeadClusterE);
  }

  for (auto jetCont : fJetCollArray) {
    name = TString::Format("%s_fHistB2BRhoVsLeadJetPt", jetCont.first.c_str());
    fHistB2BRhoVsLeadJetPt[jetCont.first] = new TH2F(name, name, nPtBins, 0, fMaxPt, nRhoBins, 0, maxRho);
    fHistB2BRhoVsLeadJetPt[jetCont.first]->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
    fHistB2BRhoVsLeadJetPt[jetCont.first]->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistB2BRhoVsLeadJetPt[jetCont.first]);
  }

  if (fScaleFunction) {
    fHistB2BRhoScaledVsCent = new TH2F("fHistB2BRhoScaledVsCent", "fHistB2BRhoScaledVsCent", 100, 0, 100, nRhoBins, 0, maxRho);
    fHistB2BRhoScaledVsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistB2BRhoScaledVsCent->GetYaxis()->SetTitle("#rho_{scaled} (GeV/#it{c} #times rad^{-1})");
    fOutput->Add(fHistB2BRhoScaledVsCent);

    if (fParticleCollArray.size() > 0) {
      fHistB2BRhoScaledVsNtrack = new TH2F("fHistB2BRhoScaledVsNtrack", "fHistB2BRhoScaledVsNtrack", 200, 0, maxTracks, nRhoBins, 0, maxRho);
      fHistB2BRhoScaledVsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistB2BRhoScaledVsNtrack->GetYaxis()->SetTitle("#rho (GeV/#it{c} #times rad^{-1})");
      fOutput->Add(fHistB2BRhoScaledVsNtrack);
    }

    if (fClusterCollArray.size() > 0) {
      fHistB2BRhoScaledVsNcluster = new TH2F("fHistB2BRhoScaledVsNcluster", "fHistB2BRhoScaledVsNcluster", 50, 0, maxTracks / 4, nRhoBins, 0, maxRho);
      fHistB2BRhoScaledVsNcluster->GetXaxis()->SetTitle("No. of clusters");
      fHistB2BRhoScaledVsNcluster->GetYaxis()->SetTitle("#rho_{scaled} (GeV/#it{c} #times rad^{-1})");
      fOutput->Add(fHistB2BRhoScaledVsNcluster);
    }
  }
}

Double_t AliAnalysisTaskRhoTransDev::GetPerpPtDensity(AliEmcalContainer* cont, AliVParticle* leadingJet)
{
  static Float_t minPhi = (3.0/8.0) * TMath::Pi();
  static Float_t maxPhi = (5.0/8.0) * TMath::Pi();

  Double_t perpPt = 0;

  for (auto mom : cont->accepted_momentum()) {
    Double_t phi_diff = TMath::Abs(AliEmcalContainer::RelativePhi(mom.first.Phi(), leadingJet->Phi()));
    if (phi_diff >=  minPhi && phi_diff <= maxPhi) perpPt += mom.first.Pt();
  }

  Double_t acc = cont->GetEtaSwing() * (maxPhi - minPhi) * 2;

  return acc > 0 ? perpPt / acc : 0;
}

void AliAnalysisTaskRhoTransDev::CalculateRho()
{
  AliEmcalJet* leadingJet = fLeadingJet["Signal"];
  if (!leadingJet) return;

  Double_t perpPtDensity = 0;

  for (auto partCont : fParticleCollArray) {
    perpPtDensity += GetPerpPtDensity(partCont.second, leadingJet);
  }

  for (auto clusCont : fClusterCollArray) {
    perpPtDensity += GetPerpPtDensity(clusCont.second, leadingJet);
  }

  fOutRho->SetVal(perpPtDensity);
}

Bool_t AliAnalysisTaskRhoTransDev::FillHistograms()
{
  Bool_t r = AliAnalysisTaskRhoBaseDev::FillHistograms();
  if (!r) return kFALSE;

  if (IsB2BEvent()) {
    fHistB2BRhoVsCent->Fill(fCent, fOutRho->GetVal());

    if (fLeadingParticle) {
      fHistB2BRhoVsLeadTrackPt->Fill(fLeadingParticle->Pt(), fOutRho->GetVal());
    }

    if (fLeadingCluster) {
      fHistB2BRhoVsLeadClusterE->Fill(fLeadingCluster->E(), fOutRho->GetVal());
    }

    if (fHistB2BRhoVsNtrack) fHistB2BRhoVsNtrack->Fill(fNtracks, fOutRho->GetVal());
    if (fHistB2BRhoVsNcluster) fHistRhoVsNcluster->Fill(fNclusters, fOutRho->GetVal());

    for (auto jetCont : fJetCollArray) {
      if (fLeadingJet[jetCont.first]) fHistB2BRhoVsLeadJetPt[jetCont.first]->Fill(fLeadingJet[jetCont.first]->Pt(), fOutRho->GetVal());
    }

    if (fOutRhoScaled) {
      fHistB2BRhoScaledVsCent->Fill(fCent, fOutRhoScaled->GetVal());
      if (fHistB2BRhoScaledVsNtrack) fHistRhoScaledVsNtrack->Fill(fNtracks, fOutRhoScaled->GetVal());
      if (fHistB2BRhoScaledVsNcluster) fHistRhoScaledVsNcluster->Fill(fNclusters,  fOutRhoScaled->GetVal());
    }
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskRhoTransDev::VerifyContainers()
{
  if (fJetCollArray.count("Signal") == 0) {
    AliError("No signal jet collection found. Task will not run!");
    return kFALSE;
  }

  if (fParticleCollArray.size() + fClusterCollArray.size() == 0) {
    AliError("No particle or cluster array was provided. Task will not run!");
    return kFALSE;
  }

  return kTRUE;
}

AliAnalysisTaskRhoTransDev* AliAnalysisTaskRhoTransDev::AddTaskRhoTransDev(TString trackName, Double_t trackPtCut, TString clusName, Double_t clusECut, TString nRho, Double_t jetradius, UInt_t acceptance, AliJetContainer::EJetType_t jetType , AliJetContainer::ERecoScheme_t rscheme, Bool_t histo, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTaskRhoTransDev::AddTaskRhoTransDev", "No analysis manager to connect to.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AliAnalysisTaskRhoTransDev::AddTaskRhoTransDev", "This task requires an input event handler");
    return nullptr;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  // Init the task and do settings
  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  TString name(TString::Format("AliAnalysisTaskRhoTransDev_%s", nRho.Data()));
  if (!suffix.IsNull()) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskRhoTransDev* mgrTask = dynamic_cast<AliAnalysisTaskRhoTransDev*>(mgr->GetTask(name.Data()));
  if (mgrTask) {
    ::Warning("AliAnalysisTaskRhoTransDev::AddTaskRhoTransDev", "Not adding the task again, since a task with the same name '%s' already exists", name.Data());
    return mgrTask;
  }

  AliAnalysisTaskRhoTransDev* rhotask = new AliAnalysisTaskRhoTransDev(name, histo);
  rhotask->SetOutRhoName(nRho);

  AliParticleContainer* partCont = rhotask->AddParticleContainer(trackName.Data());
  partCont->SetMinPt(trackPtCut);
  AliClusterContainer *clusterCont = rhotask->AddClusterContainer(clusName.Data());
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(clusECut);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  AliJetContainer *jetCont = new AliJetContainer(jetType, AliJetContainer::antikt_algorithm, rscheme, jetradius, partCont, clusterCont);
  if (jetCont) {
    jetCont->SetJetPtCut(1);
    jetCont->SetJetAcceptanceType(acceptance);
    jetCont->SetName("Signal");
    rhotask->AdoptJetContainer(jetCont);
  }

  // Final settings, pass to manager and set the containers
  mgr->AddTask(rhotask);

  // Create containers for input/output
  mgr->ConnectInput(rhotask, 0, mgr->GetCommonInputContainer());
  if (histo) {
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
        TList::Class(), AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(rhotask, 1, coutput1);
  }

  return rhotask;
}
