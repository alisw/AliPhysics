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
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <THashList.h>

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>

#include "AliEmcalJet.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliClusterContainer.h"
#include "AliParticleContainer.h"

#include "AliAnalysisTaskEmcalJetSpectraQA.h"

// Definitions of class AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfo

/// Default constructor
AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfo::AliEmcalJetInfo() :
  AliTLorentzVector(),
  fArea(0),
  fMCPt(0),
  fNConstituents(0),
  fNEF(0),
  fCent(0),
  fEP(0),
  fCorrPt(0),
  fZ(0),
  fLeadingPt(0)
{
}

/// Constructor that uses information from an AliEmcalJet object
///
/// \param jet A const reference to an AliEmcalJet object
AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfo::AliEmcalJetInfo(const AliEmcalJet& jet) :
  AliTLorentzVector(),
  fArea(jet.Area()),
  fMCPt(jet.MCPt()),
  fNConstituents(jet.GetNumberOfConstituents()),
  fNEF(jet.NEF()),
  fCent(0),
  fEP(0),
  fCorrPt(0),
  fZ(0),
  fLeadingPt(0)
{
  jet.GetMomentum(*this);
}

// Definitions of class AliAnalysisTaskEmcalJetSpectraQA

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetSpectraQA);
/// \endcond

/// Default constructor for ROOT I/O purposes
AliAnalysisTaskEmcalJetSpectraQA::AliAnalysisTaskEmcalJetSpectraQA() :
  AliAnalysisTaskEmcalJetLight("AliAnalysisTaskEmcalJetSpectraQA", kTRUE),
  fHistoType(kTHnSparse),
  fJetEPaxis(kFALSE),
  fAreaAxis(kTRUE),
  fPtBinWidth(0.5),
  fMaxPt(250),
  fIsEmbedded(kFALSE),
  fHistManager()

{
  SetMakeGeneralHistograms(kTRUE);
}

/// Standard named constructor
///
/// \param name Name of the task
AliAnalysisTaskEmcalJetSpectraQA::AliAnalysisTaskEmcalJetSpectraQA(const char *name) :
  AliAnalysisTaskEmcalJetLight(name, kTRUE),
  fHistoType(kTHnSparse),
  fJetEPaxis(kFALSE),
  fAreaAxis(kTRUE),
  fPtBinWidth(0.5),
  fMaxPt(250),
  fIsEmbedded(kFALSE),
  fHistManager(name)
{
  SetMakeGeneralHistograms(kTRUE);
}

/// This method is not implemented. It can be overloaded in derived classes to have a tree output.
///
/// \param jets Valid pointer to an AliJetContainer object
void AliAnalysisTaskEmcalJetSpectraQA::AllocateTTree(const AliJetContainer* jets)
{
  AliError("Tree output not implemented. Falling back to THnSparse output");
  AllocateTHnSparse(jets);
}

/// Allocate output THnSparse for a jet container
///
/// \param jets Valid pointer to an AliJetContainer object
void AliAnalysisTaskEmcalJetSpectraQA::AllocateTHnSparse(const AliJetContainer* jets)
{
  Double_t jetRadius = jets->GetJetRadius();

  TString title[30]= {""};
  Int_t nbins[30]  = {0};
  Double_t min[30] = {0.};
  Double_t max[30] = {0.};
  Int_t dim = 0;

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);

  if (fForceBeamType != kpp) {
    title[dim] = "Centrality (%)";
    nbins[dim] = 20;
    min[dim] = 0;
    max[dim] = 100;
    dim++;

    if (fJetEPaxis) {
      title[dim] = "#phi_{jet} - #psi_{EP}";
      nbins[dim] = nPtBins/5;
      min[dim] = 0;
      max[dim] = TMath::Pi();
      dim++;
    }
  }

  title[dim] = "#eta_{jet}";
  nbins[dim] = nPtBins/10;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "#phi_{jet} (rad)";
  nbins[dim] = nPtBins/10*3;
  min[dim] = 0;
  max[dim] = 2*TMath::Pi();
  dim++;

  title[dim] = "#it{p}_{T} (GeV/#it{c})";
  nbins[dim] = nPtBins;
  min[dim] = 0;
  max[dim] = fMaxPt;
  dim++;

  if (fIsEmbedded) {
    title[dim] = "#it{p}_{T}^{MC} (GeV/#it{c})";
    nbins[dim] = nPtBins;
    min[dim] = 0;
    max[dim] = fMaxPt;
    dim++;
  }

  if (fForceBeamType != kpp) {
    title[dim] = "#it{p}_{T}^{corr} (GeV/#it{c})";
    nbins[dim] = nPtBins;
    min[dim] = -fMaxPt/2;
    max[dim] = fMaxPt/2;
    dim++;
  }

  if (fAreaAxis) {
    // area resolution is about 0.01 (w/ ghost area 0.005)
    // for nPtBins = 250 use bin width 0.01
    title[dim] = "#it{A}_{jet}";
    nbins[dim] = TMath::CeilNint(2.0*jetRadius*jetRadius*TMath::Pi() / 0.01 * nPtBins / 250);
    min[dim] = 0;
    max[dim] = 2.0*jetRadius*jetRadius*TMath::Pi();
    dim++;
  }

  if (fClusterCollArray.size() > 0 && fParticleCollArray.size() > 0) {
    title[dim] = "NEF";
    nbins[dim] = nPtBins/5;
    min[dim] = 0;
    max[dim] = 1.0;
    dim++;
  }

  title[dim] = "#it{z}_{leading}";
  nbins[dim] = nPtBins/5;
  min[dim] = 0;
  max[dim] = 1.0;
  dim++;

  if (fForceBeamType != kpp) {
    title[dim] = "No. of constituents";
    nbins[dim] = 125;
    min[dim] = 0;
    max[dim] = 250;
    dim++;
  }
  else {
    title[dim] = "No. of constituents";
    nbins[dim] = 50;
    min[dim] = -0.5;
    max[dim] = 49.5;
    dim++;
  }

  title[dim] = "#it{p}_{T,particle}^{leading} (GeV/#it{c})";
  nbins[dim] = nPtBins/10*3;
  min[dim] = 0;
  max[dim] = 150;
  dim++;

  TString histname = TString::Format("%s/fHistJetObservables", jets->GetArrayName().Data());
  THnSparse* hn = fHistManager.CreateTHnSparse(histname.Data(), histname.Data(), dim, nbins, min, max);
  for (Int_t i = 0; i < dim; i++) {
    hn->GetAxis(i)->SetTitle(title[i]);
  }
}

/// Allocate output TH1/TH2/TH3 for a jet container
///
/// \param jets Valid pointer to an AliJetContainer object
void AliAnalysisTaskEmcalJetSpectraQA::AllocateTHX(const AliJetContainer* jets)
{
  TString histname;
  TString title;

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);

  for (Int_t i = 0; i < GetNCentBins(); i++) {
    histname = TString::Format("%s/fHistJetPtEtaPhi_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#eta;#phi (rad)";
    fHistManager.CreateTH3(histname.Data(), title.Data(), 20, -1, 1, 41, 0, 2*TMath::Pi()*41/40, nPtBins, 0, fMaxPt);

    histname = TString::Format("%s/fHistJetPtArea_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#it{A}_{jet};counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 150, 0, 1.5);

    histname = TString::Format("%s/fHistJetPtEP_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#phi_{jet} - #psi_{EP};counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 100, 0, TMath::Pi());;

    histname = TString::Format("%s/fHistJetPtNEF_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});NEF;counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 102, 0, 1.02);

    histname = TString::Format("%s/fHistJetPtZ_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#it{z}_{leading};counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 102, 0, 1.02);

    histname = TString::Format("%s/fHistJetPtLeadingPartPt_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c});counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 120, 0, 120);

    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/fHistJetCorrPtEtaPhi_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#eta;#phi (rad)";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 20, -1, 1, 41, 0, 2*TMath::Pi()*201/200, nPtBins*2, -fMaxPt, fMaxPt);

      histname = TString::Format("%s/fHistJetCorrPtArea_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#it{A}_{jet};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 150, 0, 1.5);

      histname = TString::Format("%s/fHistJetCorrPtEP_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#phi_{jet} - #psi_{EP};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 100, 0, TMath::Pi());;

      histname = TString::Format("%s/fHistJetCorrPtNEF_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});NEF;counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 102, 0, 1.02);

      histname = TString::Format("%s/fHistJetCorrPtZ_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#it{z}_{leading};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 102, 0, 1.02);

      histname = TString::Format("%s/fHistJetCorrPtLeadingPartPt_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 120, 0, 120);

      histname = TString::Format("%s/fHistJetPtCorrPt_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T} (GeV/#it{c});#it{p}_{T,corr} (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nPtBins*2, -fMaxPt, fMaxPt);

      if (fIsEmbedded) {
        histname = TString::Format("%s/fHistJetMCPtCorrPt_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,MC} (GeV/#it{c});#it{p}_{T,corr} (GeV/#it{c});counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nPtBins*2, -fMaxPt, fMaxPt);
      }
    }

    if (fIsEmbedded) {
      histname = TString::Format("%s/fHistJetPtMCPt_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T} (GeV/#it{c});#it{p}_{T,MC} (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nPtBins, 0, fMaxPt);
    }
  }
}

/// Overloads base class method. Creates output objects
void AliAnalysisTaskEmcalJetSpectraQA::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJetLight::UserCreateOutputObjects();

  Int_t maxTracks = 6000;
  Int_t constituentsNbins = 250;
  Double_t constituentsMax = 249.5;
  Double_t maxRho = 500;

  if (fForceBeamType == kpp) {
    constituentsNbins = 50;
    constituentsMax = 49.5;
    maxRho = 50;
    maxTracks = 200;
  }
  else if (fForceBeamType == kpA) {
    constituentsNbins = 100;
    constituentsMax = 99.5;
    maxRho = 200;
    maxTracks = 500;
  }

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);

  TString histname;
  TString title;

  for (auto cont_it : fJetCollArray) {
    AliJetContainer* jets = cont_it.second;
    fHistManager.CreateHistoGroup(jets->GetArrayName());

    switch (fHistoType) {
    case kTH2:
      AllocateTHX(jets);
      break;
    case kTHnSparse:
      AllocateTHnSparse(jets);
      break;
    case kTTree:
      AllocateTTree(jets);
      break;
    }

    TString histname;

    for (Int_t i = 0; i < GetNCentBins(); i++) {
      if (jets->GetParticleContainer()) {
        histname = TString::Format("%s/fHistTracksJetPt_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,track} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins / 2, 0, fMaxPt / 2, nPtBins, 0, fMaxPt);

        histname = TString::Format("%s/fHistTracksPtDist_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,track} (GeV/#it{c});#it{d};counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins / 2, 0, fMaxPt / 2, 100, 0, 5);

        histname = TString::Format("%s/fHistTracksZJetPtJetConst_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{z}_{track} (GeV/#it{c});#it{d};No. of constituents";
        fHistManager.CreateTH3(histname.Data(), title.Data(), 120, 0.0, 1.2, nPtBins, 0, fMaxPt, constituentsNbins, -0.5, constituentsMax);
      }

      if (jets->GetClusterContainer()) {
        histname = TString::Format("%s/fHistClustersJetPt_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,cluster} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins / 2, 0, fMaxPt / 2, nPtBins, 0, fMaxPt);

        histname = TString::Format("%s/fHistClustersPtDist_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,cluster} (GeV/#it{c});#it{d};counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins / 2, 0, fMaxPt / 2, 100, 0, 5);

        histname = TString::Format("%s/fHistClustersZJetPtJetConst_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{z}_{cluster} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});No. of constituents";
        fHistManager.CreateTH3(histname.Data(), title.Data(), 120, 0.0, 1.2, nPtBins, 0, fMaxPt, constituentsNbins, -0.5, constituentsMax);
      }

      histname = TString::Format("%s/fHistRejectionReason_%d", jets->GetArrayName().Data(), i);
      title = histname + ";Rejection reason;#it{p}_{T,jet} (GeV/#it{c});counts";
      TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 100, 0, 250);
      SetRejectionReasonLabels(hist->GetXaxis());
    }

    histname = TString::Format("%s/fHistLeadJetPtVsCent", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T,jet} (GeV/#it{c});counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, nPtBins, 0, fMaxPt);

    histname = TString::Format("%s/fHistLeadJetPtVsNTracks", jets->GetArrayName().Data());
    title = histname + ";no. of tracks;#it{p}_{T,jet} (GeV/#it{c});counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 200, 0, maxTracks, nPtBins, 0, fMaxPt);

    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/fHistRhoVsCent", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#rho (GeV/#it{c} rad^{-1});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 100, 1000, 0, maxRho);

      histname = TString::Format("%s/fHistRhoVsNTracks", jets->GetArrayName().Data());
      title = histname + ";no. of tracks;#rho (GeV/#it{c} rad^{-1});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 200, 0, maxTracks, 1000, 0, maxRho);

      histname = TString::Format("%s/fHistRhoVsLeadJetPt", jets->GetArrayName().Data());
      title = histname + ";#it{p}_{T,jet} (GeV/#it{c});#rho (GeV/#it{c} rad^{-1});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, 1000, 0, maxRho);
    }
  }

  TIter nextElement(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = nextElement())) fOutput->Add(obj);
  PostData(1, fOutput);
}

/// Overloads base class method. Fills the output histograms
///
/// \return kTRUE if successful
Bool_t AliAnalysisTaskEmcalJetSpectraQA::FillHistograms()
{
  TString histname;

  if (fCentBin < 0) {
    AliError(Form("fCentBin is %d! fCent = %.3f. Fix the centrality bins to include all possible values of centrality.", fCentBin, fCent));
    return kFALSE;
  }

  for (auto cont_it : fJetCollArray) {
    AliJetContainer* jets = cont_it.second;
    Double_t rhoVal = 0;
    if (jets->GetRhoParameter()) rhoVal = jets->GetRhoVal();

    Double_t leadJetPt = 0;
    for (auto jet : jets->accepted()) {

      UInt_t rejectionReason = 0;
      if (!jets->AcceptJet(jet, rejectionReason)) {
        histname = TString::Format("%s/fHistRejectionReason_%d", jets->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH2(histname.Data(), jets->GetRejectionReasonBitPosition(rejectionReason), jet->Pt());
        continue;
      }

      Float_t ptLeading = jets->GetLeadingHadronPt(jet);
      Float_t corrPt = jet->Pt() - rhoVal * jet->Area();

      TLorentzVector leadPart;

      jets->GetLeadingHadronMomentum(leadPart, jet);

      // Fill THnSparse
      Double_t ep = jet->Phi() - fEPV0;
      while (ep < 0) ep += TMath::Pi();
      while (ep >= TMath::Pi()) ep -= TMath::Pi();

      Double_t z = GetParallelFraction(leadPart.Vect(), jet);
      if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1

      AliEmcalJetInfo jetInfo(*jet);
      jetInfo.fCent = fCent;
      jetInfo.fEP = ep;
      jetInfo.fCorrPt = corrPt;
      jetInfo.fZ = z;
      jetInfo.fLeadingPt = ptLeading;

      if (jet->Pt() > leadJetPt) leadJetPt = jet->Pt();

      FillJetHisto(jetInfo, jets);

      AliParticleContainer* tracks = jets->GetParticleContainer();
      if (tracks) {
        for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
          AliVParticle *track = jet->TrackAt(it, tracks->GetArray());
          if (track) {
            histname = TString::Format("%s/fHistTracksJetPt_%d", jets->GetArrayName().Data(), fCentBin);
            fHistManager.FillTH2(histname.Data(), track->Pt(), jet->Pt());

            Double_t dphi = TVector2::Phi_0_2pi(track->Phi() - jet->Phi());
            Double_t deta = track->Eta() - jet->Eta();
            Double_t dist = TMath::Sqrt(deta * deta + dphi * dphi);

            histname = TString::Format("%s/fHistTracksPtDist_%d", jets->GetArrayName().Data(), fCentBin);
            fHistManager.FillTH2(histname.Data(), track->Pt(), dist);

            histname = TString::Format("%s/fHistTracksZJetPtJetConst_%d", jets->GetArrayName().Data(), fCentBin);
            fHistManager.FillTH3(histname.Data(), GetParallelFraction(track, jet), jet->Pt(), jet->GetNumberOfConstituents());
          }
        }
      }

      AliClusterContainer* clusters = jets->GetClusterContainer();
      if (clusters) {
        for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
          AliVCluster *cluster = jet->ClusterAt(ic, clusters->GetArray());

          if (cluster) {
            TLorentzVector nPart;
            if (clusters->GetDefaultClusterEnergy() >=0 && clusters->GetDefaultClusterEnergy() < AliVCluster::kLastUserDefEnergy) {
              cluster->GetMomentum(nPart, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
            }
            else {
              cluster->GetMomentum(nPart, fVertex);
            }

            histname = TString::Format("%s/fHistClustersJetPt_%d", jets->GetArrayName().Data(), fCentBin);
            fHistManager.FillTH2(histname.Data(), nPart.Pt(), jet->Pt());

            Double_t dphi = TVector2::Phi_0_2pi(nPart.Phi() - jet->Phi());
            Double_t deta = nPart.Eta() - jet->Eta();
            Double_t dist = TMath::Sqrt(deta * deta + dphi * dphi);

            histname = TString::Format("%s/fHistClustersPtDist_%d", jets->GetArrayName().Data(), fCentBin);
            fHistManager.FillTH2(histname.Data(), nPart.Pt(), dist);

            histname = TString::Format("%s/fHistClustersZJetPtJetConst_%d", jets->GetArrayName().Data(), fCentBin);
            fHistManager.FillTH3(histname.Data(), GetParallelFraction(nPart.Vect(), jet), jet->Pt(), jet->GetNumberOfConstituents());
          }
        }
      }
    } //jet loop

    Int_t ntracks = 0;
    for (auto cont : this->fParticleCollArray) ntracks += cont.second->GetNAcceptEntries();

    histname = TString::Format("%s/fHistLeadJetPtVsCent", jets->GetArrayName().Data());
    fHistManager.FillTH2(histname.Data(), fCent, leadJetPt);

    histname = TString::Format("%s/fHistLeadJetPtVsNTracks", jets->GetArrayName().Data());
    fHistManager.FillTH2(histname.Data(), ntracks, leadJetPt);

    if (jets->GetRhoParameter()) {
      histname = TString::Format("%s/fHistRhoVsCent", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), fCent, rhoVal);

      histname = TString::Format("%s/fHistRhoVsNTracks", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), ntracks, rhoVal);

      histname = TString::Format("%s/fHistRhoVsLeadJetPt", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), leadJetPt, rhoVal);
    }
  }
  return kTRUE;
}

/// Fill histograms with jet
///
/// \param jet  Jet containing the information to be sent to the tree/histograms
/// \param jets Jet container
void AliAnalysisTaskEmcalJetSpectraQA::FillTHX(const AliEmcalJetInfo& jet, const AliJetContainer* jets)
{
  TString histname;

  histname = TString::Format("%s/fHistJetPtEtaPhi_%d", jets->GetArrayName().Data(), fCentBin);
  fHistManager.FillTH3(histname.Data(), jet.Eta(), jet.Phi_0_2pi(), jet.Pt());

  histname = TString::Format("%s/fHistJetPtArea_%d", jets->GetArrayName().Data(), fCentBin);
  fHistManager.FillTH2(histname.Data(), jet.Pt(), jet.fArea);

  histname = TString::Format("%s/fHistJetPtEP_%d", jets->GetArrayName().Data(), fCentBin);
  fHistManager.FillTH2(histname.Data(), jet.Pt(), jet.fEP);

  histname = TString::Format("%s/fHistJetPtNEF_%d", jets->GetArrayName().Data(), fCentBin);
  fHistManager.FillTH2(histname.Data(), jet.Pt(), jet.fNEF);

  histname = TString::Format("%s/fHistJetPtZ_%d", jets->GetArrayName().Data(), fCentBin);
  fHistManager.FillTH2(histname.Data(), jet.Pt(), jet.fZ);

  histname = TString::Format("%s/fHistJetPtLeadingPartPt_%d", jets->GetArrayName().Data(), fCentBin);
  fHistManager.FillTH2(histname.Data(), jet.Pt(), jet.fLeadingPt);

  if (fIsEmbedded) {
    histname = TString::Format("%s/fHistJetPtMCPt_%d", jets->GetArrayName().Data(), fCentBin);
    fHistManager.FillTH2(histname.Data(), jet.Pt(), jet.fMCPt);
  }

  if (!jets->GetRhoName().IsNull()) {
    histname = TString::Format("%s/fHistJetCorrPtEtaPhi_%d", jets->GetArrayName().Data(), fCentBin);
    fHistManager.FillTH3(histname.Data(), jet.Eta(), jet.Phi_0_2pi(), jet.fCorrPt);

    histname = TString::Format("%s/fHistJetCorrPtArea_%d", jets->GetArrayName().Data(), fCentBin);
    fHistManager.FillTH2(histname.Data(), jet.fCorrPt, jet.fArea);

    histname = TString::Format("%s/fHistJetCorrPtEP_%d", jets->GetArrayName().Data(), fCentBin);
    fHistManager.FillTH2(histname.Data(), jet.fCorrPt, jet.fEP);

    histname = TString::Format("%s/fHistJetCorrPtNEF_%d", jets->GetArrayName().Data(), fCentBin);
    fHistManager.FillTH2(histname.Data(), jet.fCorrPt, jet.fNEF);

    histname = TString::Format("%s/fHistJetCorrPtZ_%d", jets->GetArrayName().Data(), fCentBin);
    fHistManager.FillTH2(histname.Data(), jet.fCorrPt, jet.fZ);

    histname = TString::Format("%s/fHistJetCorrPtLeadingPartPt_%d", jets->GetArrayName().Data(), fCentBin);
    fHistManager.FillTH2(histname.Data(), jet.fCorrPt, jet.fLeadingPt);

    histname = TString::Format("%s/fHistJetPtCorrPt_%d", jets->GetArrayName().Data(), fCentBin);
    fHistManager.FillTH2(histname.Data(), jet.Pt(), jet.fCorrPt);

    if (fIsEmbedded) {
      histname = TString::Format("%s/fHistJetMCPtCorrPt_%d", jets->GetArrayName().Data(), fCentBin);
      fHistManager.FillTH2(histname.Data(), jet.fMCPt, jet.fCorrPt);
    }
  }
}

/// This method is not implemented. It can be overloaded in derived classes to have a tree output.
///
/// \param jet  Jet containing the information to be sent to the tree/histograms
/// \param jets Jet container
void AliAnalysisTaskEmcalJetSpectraQA::FillTTree(const AliEmcalJetInfo& jet, const AliJetContainer* jets)
{
  AliError("Tree output not implemented. Falling back to THnSparse output");
  FillTHnSparse(jet, jets);
}

/// Fill THnSparse histogram with jet
///
/// \param jet  Jet containing the information to be sent to the tree/histograms
/// \param jets Jet container
void AliAnalysisTaskEmcalJetSpectraQA::FillTHnSparse(const AliEmcalJetInfo& jet, const AliJetContainer* jets)
{
  TString histname;
  Double_t contents[30]={0};

  histname = TString::Format("%s/fHistJetObservables", jets->GetArrayName().Data());
  THnSparse* histJetObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));

  if (!histJetObservables) return;

  for (Int_t i = 0; i < histJetObservables->GetNdimensions(); i++) {
    TString title(histJetObservables->GetAxis(i)->GetTitle());
    if (title=="Centrality (%)")
      contents[i] = jet.fCent;
    else if (title=="#phi_{jet} - #psi_{EP}")
      contents[i] = jet.fEP;
    else if (title=="#eta_{jet}")
      contents[i] = jet.Eta();
    else if (title=="#phi_{jet} (rad)")
      contents[i] = jet.Phi_0_2pi();
    else if (title=="#it{p}_{T} (GeV/#it{c})")
      contents[i] = jet.Pt();
    else if (title=="#it{p}_{T}^{MC} (GeV/#it{c})")
      contents[i] = jet.fMCPt;
    else if (title=="#it{p}_{T}^{corr} (GeV/#it{c})")
      contents[i] = jet.fCorrPt;
    else if (title=="#it{A}_{jet}")
      contents[i] = jet.fArea;
    else if (title=="NEF")
      contents[i] = jet.fNEF;
    else if (title=="#it{z}_{leading}")
      contents[i] = jet.fZ;
    else if (title=="No. of constituents")
      contents[i] = jet.fNConstituents;
    else if (title=="#it{p}_{T,particle}^{leading} (GeV/#it{c})")
      contents[i] = jet.fLeadingPt;
    else
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }

  histJetObservables->Fill(contents);
}

/// Fill histogram with jet
///
/// \param jet  Jet containing the information to be sent to the tree/histograms
/// \param jets Jet container
void AliAnalysisTaskEmcalJetSpectraQA::FillJetHisto(const AliEmcalJetInfo& jet, const AliJetContainer* jets)
{
  switch (fHistoType) {
  case kTH2:
    FillTHX(jet, jets);
    break;

  case kTHnSparse:
    FillTHnSparse(jet, jets);
    break;

  case kTTree:
    FillTTree(jet, jets);
    break;
  }
}

/// Create an instance of this class and add it to the analysis manager
///
/// \param ntracks name of the track collection
/// \param nclusters name of the calorimeter cluster collection
/// \param trackPtCut minimum transverse momentum of tracks
/// \param clusECut minimum energy of calorimeter clusters
/// \param suffix additional suffix that can be added at the end of the task name
/// \return pointer to the new AddTaskEmcalJetSpectraQA task
AliAnalysisTaskEmcalJetSpectraQA* AliAnalysisTaskEmcalJetSpectraQA::AddTaskEmcalJetSpectraQA(TString trackName, TString clusName,
    Double_t trackPtCut, Double_t clusECut, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalJetSpectraQA", "No analysis manager to connect to.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AddTaskEmcalJetSpectraQA", "This task requires an input event handler");
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

  TString name("AliAnalysisTaskEmcalJetSpectraQA");
  if (strcmp(suffix,"")) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskEmcalJetSpectraQA* jetTask = new AliAnalysisTaskEmcalJetSpectraQA(name);
  jetTask->SetVzRange(-10,10);
  jetTask->SetNeedEmcalGeom(kFALSE);
  AliParticleContainer *partCont = jetTask->AddParticleContainer(trackName.Data());
  if (partCont) partCont->SetParticlePtCut(trackPtCut);

  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(clusName.Data());
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(clusECut);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  // Final settings, pass to manager and set the containers
  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                  TList::Class(),AliAnalysisManager::kOutputContainer,
                  Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );

  return jetTask;
}
