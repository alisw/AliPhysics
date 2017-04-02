/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskEmcalJetSample.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetSample);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalJetSample::AliAnalysisTaskEmcalJetSample() : 
  AliAnalysisTaskEmcalJet(),
  fHistManager()
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetSample::AliAnalysisTaskEmcalJetSample(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistManager(name)
{
  SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetSample::~AliAnalysisTaskEmcalJetSample()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetSample::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AllocateClusterHistograms();
  AllocateTrackHistograms();
  AllocateJetHistograms();
  AllocateCellHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for basic EMCal cluster QA.
 * A set of histograms (energy, eta, phi, number of cluster) is allocated
 * per each cluster container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetSample::AllocateClusterHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliClusterContainer* clusCont = 0;
  TIter next(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
    groupname = clusCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The cluster containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{exotic} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{non-lin.corr.} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{had.corr.} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{custer};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{custer};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      histname = TString::Format("%s/histNClusters_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of clusters;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    }
  }

  histname = "fHistSumNClusters";
  histtitle = TString::Format("%s;Sum of n clusters;events", histname.Data());
  if (fForceBeamType != kpp) {
    fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000);
  }
  else {
    fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
  }
}

/*
 * This function allocates the histograms for basic EMCal QA.
 * One 2D histogram with the cell energy spectra and the number of cells
 * per event is allocated per each centrality bin.
 */
void AliAnalysisTaskEmcalJetSample::AllocateCellHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname(fCaloCellsName);

  fHistManager.CreateHistoGroup(groupname);
  for (Int_t cent = 0; cent < fNcentBins; cent++) {
    histname = TString::Format("%s/histCellEnergy_%d", groupname.Data(), cent);
    histtitle = TString::Format("%s;#it{E}_{cell} (GeV);counts", histname.Data());
    fHistManager.CreateTH1(histname, histtitle, 300, 0, 150);

    histname = TString::Format("%s/histNCells_%d", groupname.Data(), cent);
    histtitle = TString::Format("%s;number of cells;events", histname.Data());
    if (fForceBeamType != kpp) {
      fHistManager.CreateTH1(histname, histtitle, 500, 0, 6000);
    }
    else {
      fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
    }
  }
}

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetSample::AllocateTrackHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The track containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
        histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#eta}_{track}^{vertex} - #it{#eta}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 50, -0.5, 0.5);

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#phi}_{track}^{vertex} - #it{#phi}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 200, -2, 2);

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{p}_{T,track}^{vertex} - #it{p}_{T,track}^{EMCal} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, -fMaxBinPt/2, fMaxBinPt/2);

        histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{P}_{track} (GeV/#it{c});#it{E}_{cluster} / #it{P}_{track} #it{c};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, 0, 4);
      }

      histname = TString::Format("%s/histNTracks_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of tracks;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    }
  }

  histname = "fHistSumNTracks";
  histtitle = TString::Format("%s;Sum of n tracks;events", histname.Data());
  if (fForceBeamType != kpp) {
    fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
  }
  else {
    fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
  }
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetSample::AllocateJetHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, 3);

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      histname = TString::Format("%s/histNJets_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of jets;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 500);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 100, 0, 100);
      }

      if (!jetCont->GetRhoName().IsNull()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
      }
    }
  }
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetSample::FillHistograms()
{
  DoJetLoop();
  DoTrackLoop();
  DoClusterLoop();
  DoCellLoop();

  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetSample::DoJetLoop()
{
  TString histname;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    UInt_t count = 0;
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;

      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Pt());

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Area());

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Phi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Eta());

      if (jetCont->GetRhoParameter()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());
      }
    }
    histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetSample::DoTrackLoop()
{
  AliClusterContainer* clusCont = GetClusterContainer(0);

  TString histname;
  TString groupname;
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    UInt_t count = 0;
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      count++;

      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Pt());

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Phi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Eta());

      if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
        const AliVTrack* track = static_cast<const AliVTrack*>(part);

        histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Eta() - track->GetTrackEtaOnEMCal());

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Phi() - track->GetTrackPhiOnEMCal());

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Pt() - track->GetTrackPtOnEMCal());

        if (clusCont) {
          Int_t iCluster = track->GetEMCALcluster();
          if (iCluster >= 0) {
            AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
            if (cluster) {
              histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), fCentBin);
              fHistManager.FillTH2(histname, track->P(), cluster->GetNonLinCorrEnergy() / track->P());
            }
          }
        }
      }
    }
    sumAcceptedTracks += count;

    histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }

  histname = "fHistSumNTracks";
  fHistManager.FillTH1(histname, sumAcceptedTracks);
}

/**
 * This function performs a loop over the reconstructed EMCal clusters
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetSample::DoClusterLoop()
{
  TString histname;
  TString groupname;
  UInt_t sumAcceptedClusters = 0;
  AliClusterContainer* clusCont = 0;
  TIter next(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
    groupname = clusCont->GetName();

    for(auto cluster : clusCont->all()) {
      if (!cluster) continue;

      if (cluster->GetIsExotic()) {
        histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, cluster->E());
      }
    }

    UInt_t count = 0;
    for(auto cluster : clusCont->accepted()) {
      if (!cluster) continue;
      count++;

      AliTLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);

      histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->E());

      histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->GetNonLinCorrEnergy());

      histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->GetHadCorrEnergy());

      histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, nPart.Phi_0_2pi());

      histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, nPart.Eta());
    }
    sumAcceptedClusters += count;

    histname = TString::Format("%s/histNClusters_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }

  histname = "fHistSumNClusters";
  fHistManager.FillTH1(histname, sumAcceptedClusters);
}

/**
 * This function performs a loop over the reconstructed EMCal cells
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetSample::DoCellLoop()
{
  if (!fCaloCells) return;

  TString histname;

  const Short_t ncells = fCaloCells->GetNumberOfCells();

  histname = TString::Format("%s/histNCells_%d", fCaloCellsName.Data(), fCentBin);
  fHistManager.FillTH1(histname, ncells);

  histname = TString::Format("%s/histCellEnergy_%d", fCaloCellsName.Data(), fCentBin);
  for (Short_t pos = 0; pos < ncells; pos++) {
    Double_t amp   = fCaloCells->GetAmplitude(pos);
    
    fHistManager.FillTH1(histname, amp);
  }
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetSample::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetSample::Run()
{
  return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetSample::Terminate(Option_t *) 
{
}
