/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliAnalysisTaskRhoTransDev.h"

#include <TClonesArray.h>
#include <TMath.h>

#include <AliLog.h>
#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskRhoTransDev);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
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

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name  Name of the task
 * @param[in] histo If kTRUE, the task will also produce QA histograms
 */
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

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
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

/**
 * Calculates the transervse momentum density in a region that
 * is perpedincular to the leading jet. The perpedincular region
 * spans over an angle pi/4 around pi/4 and -pi/4.
 * This implements the techniques described in the analysis note
 * of the ALICE inclusive jet spectrum in pp collisions at 2.76 TeV,
 * section 7.3: https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/rma/2013-Mar-29-analysis_note-ppJet_note.pdf
 * @param cont Collection of objects (tracks, clusters or particles) used to calculate the momentum density
 * @param leadingJet Leading jet used to define the perpendicular region
 * @return The perpendicular momentum density
 */
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

/**
 * Calculates the average background by summing the transverse momenta
 * of the particles perpendicular to the leading jet axis.
 * This implements the techniques described in the analysis note
 * of the ALICE inclusive jet spectrum in pp collisions at 2.76 TeV,
 * section 7.3: https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/rma/2013-Mar-29-analysis_note-ppJet_note.pdf
 * The background densisty is stored in fOutRho.
 */
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

/**
 * Fill histograms.
 */
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

/**
 * Verify that the required particle, cluster and jet containers were provided.
 * @return kTRUE if all requirements are satisfied, kFALSE otherwise
 */
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

/**
 * Create an instance of this class and add it to the analysis manager
 * @param trackName name of the track collection
 * @param trackPtCut minimum pt of the tracks
 * @param clusName name of the calorimeter cluster collection
 * @param clusECut minimum energy of the calorimeter clustuers
 * @param nRho name of the output rho object
 * @param jetradius Radius of the kt jets used to calculate the background
 * @param acceptance Fiducial acceptance of the kt jets
 * @param jetType Jet type (full/charged)
 * @param rscheme Recombination scheme
 * @param histo If kTRUE the task will also produce QA histograms
 * @param suffix additional suffix that can be added at the end of the task name
 * @return pointer to the new AliAnalysisTaskRhoDev task
 */
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
