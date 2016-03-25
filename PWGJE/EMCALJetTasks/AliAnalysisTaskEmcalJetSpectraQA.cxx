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

// Definitions of class AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryBase

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryBase)
/// \endcond

/// Constructor that sets the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryBase::AliEmcalJetInfoSummaryBase(const AliEmcalJetInfo& source) :
  fPt(0),
  fEta(0),
  fPhi(0),
  fNEF(0),
  fZLeading(0)
{
  Set(source);
}

/// Reset the object
void AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryBase::Reset()
{
  fPt = 0;
  fEta = 0;
  fPhi = 0;
  fNEF = 0;
  fZLeading = 0;
}

/// Set the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
void AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryBase::Set(const AliEmcalJetInfo& source)
{
  fPt = source.Pt();
  fEta = source.Eta();
  fPhi = source.Phi_0_2pi();
  fNEF = source.fNEF;
  fZLeading = source.fZ;
}

// Definitions of class AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPP

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPP)
/// \endcond

/// Constructor that sets the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPP::AliEmcalJetInfoSummaryPP(const AliEmcalJetInfo& source) :
  AliEmcalJetInfoSummaryBase(),
  fNConstituents(0)
{
  Set(source);
}

/// Reset the object
void AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPP::Reset()
{
  AliEmcalJetInfoSummaryBase::Reset();
  fNConstituents = 0;
}

/// Set the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
void AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPP::Set(const AliEmcalJetInfo& source)
{
  AliEmcalJetInfoSummaryBase::Set(source);
  fNConstituents = Char_t(source.fNConstituents);
}

// Definitions of class AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPbPb

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPbPb)
/// \endcond

/// Constructor that sets the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPbPb::AliEmcalJetInfoSummaryPbPb(const AliEmcalJetInfo& source) :
  AliEmcalJetInfoSummaryBase(),
  fCent(0),
  fEP(0),
  fArea(0),
  fNConstituents(0),
  fCorrPt(0)
{
  Set(source);
}

/// Reset the object
void AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPbPb::Reset()
{
  AliEmcalJetInfoSummaryBase::Reset();
  fCent = 0;
  fEP = 0;
  fArea = 0;
  fNConstituents = 0;
  fCorrPt = 0;
}

/// Set the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
void AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPbPb::Set(const AliEmcalJetInfo& source)
{
  AliEmcalJetInfoSummaryBase::Set(source);
  fCent = Char_t(source.fCent);
  fEP = source.fEP;
  fArea = source.fArea;
  fNConstituents = Short_t(source.fNConstituents);
  fCorrPt = source.fCorrPt;
}

// Definitions of class AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPbPb

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryEmbedding)
/// \endcond

/// Constructor that sets the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryEmbedding::AliEmcalJetInfoSummaryEmbedding(const AliEmcalJetInfo& source) :
  AliEmcalJetInfoSummaryPbPb(),
  fMCPt(0)
{
  Set(source);
}

/// Reset the object
void AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryEmbedding::Reset()
{
  AliEmcalJetInfoSummaryPbPb::Reset();
  fMCPt = 0;
}

/// Set the object copying information from an AliEmcalJetInfo object
///
/// \param source Const reference to an AliEmcalJetInfo object to copy from
void AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryEmbedding::Set(const AliEmcalJetInfo& source)
{
  AliEmcalJetInfoSummaryPbPb::Set(source);
  fMCPt = source.fMCPt;
}

// Definitions of class AliAnalysisTaskEmcalJetSpectraQA

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetSpectraQA)
/// \endcond

/// Default constractor for ROOT I/O purposes
AliAnalysisTaskEmcalJetSpectraQA::AliAnalysisTaskEmcalJetSpectraQA() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetSpectraQA", kTRUE),
  fHistoType(kTHnSparse),
  fJetEPaxis(kFALSE),
  fAreaAxis(kTRUE),
  fHistManager(),
  fCurrentJetInfo(0),
  fTreeList(0)

{
  SetMakeGeneralHistograms(kTRUE);
}

/// Standard named constructor
///
/// \param name Name of the task
AliAnalysisTaskEmcalJetSpectraQA::AliAnalysisTaskEmcalJetSpectraQA(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistoType(kTHnSparse),
  fJetEPaxis(kFALSE),
  fAreaAxis(kTRUE),
  fHistManager(name),
  fCurrentJetInfo(0),
  fTreeList(0)
{
  SetMakeGeneralHistograms(kTRUE);
}

/// Allocate output TTree for a jet container
///
/// \param jets Valid pointer to an AliJetContainer object
void AliAnalysisTaskEmcalJetSpectraQA::AllocateTTree(const AliJetContainer* jets)
{
  TString classname;

  if (!fTreeList) {
    fTreeList = new THashList();
    fTreeList->SetName("Jets");
  }

  if (fForceBeamType == kpp) {
    classname = "AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPP";
    fCurrentJetInfo = new AliEmcalJetInfoSummaryPP();
  }
  else {
    if (fIsEmbedded) {
      classname = "AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryEmbedding";
      fCurrentJetInfo = new AliEmcalJetInfoSummaryEmbedding();
    }
    else {
      classname = "AliAnalysisTaskEmcalJetSpectraQA::AliEmcalJetInfoSummaryPbPb";
      fCurrentJetInfo = new AliEmcalJetInfoSummaryPbPb();
    }
  }

  TTree* tree = new TTree(jets->GetName(), jets->GetName());
  tree->Branch("Jets", classname, &fCurrentJetInfo, 32000, 0);
  fTreeList->Add(tree);
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

  if (fForceBeamType != kpp) {
    title[dim] = "Centrality (%)";
    nbins[dim] = 20;
    min[dim] = 0;
    max[dim] = 100;
    dim++;

    if (fJetEPaxis) {
      title[dim] = "#phi_{jet} - #psi_{EP}";
      nbins[dim] = fNbins/5;
      min[dim] = 0;
      max[dim] = TMath::Pi();
      dim++;
    }
  }

  title[dim] = "#eta_{jet}";
  nbins[dim] = fNbins/10;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "#phi_{jet} (rad)";
  nbins[dim] = fNbins/10*3;
  min[dim] = 0;
  max[dim] = 2*TMath::Pi();
  dim++;

  title[dim] = "#it{p}_{T} (GeV/#it{c})";
  nbins[dim] = fNbins;
  min[dim] = fMinBinPt;
  max[dim] = fMaxBinPt;
  dim++;

  if (fIsEmbedded) {
    title[dim] = "#it{p}_{T}^{MC} (GeV/#it{c})";
    nbins[dim] = fNbins;
    min[dim] = fMinBinPt;
    max[dim] = fMaxBinPt;
    dim++;
  }

  if (!GetRhoName().IsNull()) {
    title[dim] = "#it{p}_{T}^{corr} (GeV/#it{c})";
    nbins[dim] = fNbins;
    min[dim] = -fMaxBinPt/2;
    max[dim] = fMaxBinPt/2;
    dim++;
  }

  if (fAreaAxis) {
    // area resolution is about 0.01 (w/ ghost area 0.005)
    // for fNbins = 250 use bin width 0.01
    title[dim] = "#it{A}_{jet}";
    nbins[dim] = TMath::CeilNint(2.0*jetRadius*jetRadius*TMath::Pi() / 0.01 * fNbins / 250);
    min[dim] = 0;
    max[dim] = 2.0*jetRadius*jetRadius*TMath::Pi();
    dim++;
  }

  if (fClusterCollArray.GetEntriesFast() > 0 && fParticleCollArray.GetEntriesFast() > 0) {
    title[dim] = "NEF";
    nbins[dim] = fNbins/5;
    min[dim] = 0;
    max[dim] = 1.0;
    dim++;
  }

  title[dim] = "#it{z}_{leading}";
  nbins[dim] = fNbins/5;
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
  nbins[dim] = fNbins/10*3;
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

  for (Int_t i = 0; i < fNcentBins; i++) {
    histname = TString::Format("%s/fHistJetPtEtaPhi_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#eta;#phi (rad)";
    fHistManager.CreateTH3(histname.Data(), title.Data(), 20, -1, 1, 41, 0, 2*TMath::Pi()*41/40, fNbins, fMinBinPt, fMaxBinPt);

    histname = TString::Format("%s/fHistJetPtArea_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#it{A}_{jet};counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 150, 0, 1.5);

    histname = TString::Format("%s/fHistJetPtEP_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#phi_{jet} - #psi_{EP};counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 100, 0, TMath::Pi());;

    histname = TString::Format("%s/fHistJetPtNEF_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});NEF;counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 102, 0, 1.02);

    histname = TString::Format("%s/fHistJetPtZ_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#it{z}_{leading};counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 102, 0, 1.02);

    histname = TString::Format("%s/fHistJetPtLeadingPartPt_%d", jets->GetArrayName().Data(), i);
    title = histname + ";#it{p}_{T} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c});counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 120, 0, 120);

    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/fHistJetCorrPtEtaPhi_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#eta;#phi (rad)";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 20, -1, 1, 41, 0, 2*TMath::Pi()*201/200, fNbins*2, -fMaxBinPt, fMaxBinPt);

      histname = TString::Format("%s/fHistJetCorrPtArea_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#it{A}_{jet};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 150, 0, 1.5);

      histname = TString::Format("%s/fHistJetCorrPtEP_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#phi_{jet} - #psi_{EP};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 100, 0, TMath::Pi());;

      histname = TString::Format("%s/fHistJetCorrPtNEF_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});NEF;counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 102, 0, 1.02);

      histname = TString::Format("%s/fHistJetCorrPtZ_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#it{z}_{leading};counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 102, 0, 1.02);

      histname = TString::Format("%s/fHistJetCorrPtLeadingPartPt_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T,corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, 120, 0, 120);

      histname = TString::Format("%s/fHistJetPtCorrPt_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T} (GeV/#it{c});#it{p}_{T,corr} (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);

      if (fIsEmbedded) {
        histname = TString::Format("%s/fHistJetMCPtCorrPt_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,MC} (GeV/#it{c});#it{p}_{T,corr} (GeV/#it{c});counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);
      }
    }

    if (fIsEmbedded) {
      histname = TString::Format("%s/fHistJetPtMCPt_%d", jets->GetArrayName().Data(), i);
      title = histname + ";#it{p}_{T} (GeV/#it{c});#it{p}_{T,MC} (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    }
  }
}

/// Overloads base class method. Creates output objects
void AliAnalysisTaskEmcalJetSpectraQA::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Int_t constituentsNbins = 250;
  Double_t constituentsMax = 249.5;

  if (fForceBeamType == kpp) {
    constituentsNbins = 50;
    constituentsMax = 49.5;
  }

  TString histname;
  TString title;

  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
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

    for (Int_t i = 0; i < fNcentBins; i++) {
      if (jets->GetParticleContainer()) {
        histname = TString::Format("%s/fHistTracksJetPt_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,track} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2, fNbins, fMinBinPt, fMaxBinPt);

        histname = TString::Format("%s/fHistTracksPtDist_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,track} (GeV/#it{c});#it{d};counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2, 100, 0, 5);

        histname = TString::Format("%s/fHistTracksZJetPtJetConst_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{z}_{track} (GeV/#it{c});#it{d};No. of constituents";
        fHistManager.CreateTH3(histname.Data(), title.Data(), 120, 0.0, 1.2, fNbins, fMinBinPt, fMaxBinPt, constituentsNbins, -0.5, constituentsMax);
      }

      if (jets->GetClusterContainer()) {
        histname = TString::Format("%s/fHistClustersJetPt_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,cluster} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2, fNbins, fMinBinPt, fMaxBinPt);

        histname = TString::Format("%s/fHistClustersPtDist_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{p}_{T,cluster} (GeV/#it{c});#it{d};counts";
        fHistManager.CreateTH2(histname.Data(), title.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2, 100, 0, 5);

        histname = TString::Format("%s/fHistClustersZJetPtJetConst_%d", jets->GetArrayName().Data(), i);
        title = histname + ";#it{z}_{cluster} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});No. of constituents";
        fHistManager.CreateTH3(histname.Data(), title.Data(), 120, 0.0, 1.2, fNbins, fMinBinPt, fMaxBinPt, constituentsNbins, -0.5, constituentsMax);
      }

      histname = TString::Format("%s/fHistRejectionReason_%d", jets->GetArrayName().Data(), i);
      title = histname + ";Rejection reason;#it{p}_{T,jet} (GeV/#it{c});counts";
      TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 100, 0, 250);
      SetRejectionReasonLabels(hist->GetXaxis());
    }

    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/fHistRhoVsCent", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 101, 0, 101, 100, 0, 500);
    }
  }

  fOutput->Add(fHistManager.GetListOfHistograms());

  if (fHistoType == kTTree) fOutput->Add(fTreeList);

  PostData(1, fOutput);
}

/// Overloads base class method. Fills the output histograms
///
/// \return kTRUE if successful
Bool_t AliAnalysisTaskEmcalJetSpectraQA::FillHistograms()
{
  // Fill histograms.

  TString histname;

  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    if (jets->GetRhoParameter()) {
      histname = TString::Format("%s/fHistRhoVsCent", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), fCent, jets->GetRhoVal());
    }

    AliEmcalJet* jet = 0;
    jets->ResetCurrentID();
    while ((jet = jets->GetNextJet())) {

      UInt_t rejectionReason = 0;
      if (!jets->AcceptJet(jet, rejectionReason)) {
        histname = TString::Format("%s/fHistRejectionReason_%d", jets->GetArrayName().Data(), fCentBin);
        fHistManager.FillTH2(histname.Data(), jets->GetRejectionReasonBitPosition(rejectionReason), jet->Pt());
        continue;
      }

      Float_t ptLeading = jets->GetLeadingHadronPt(jet);
      Float_t corrPt = jet->Pt() - fRhoVal * jet->Area();

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

/// Fill tree with jet info
///
/// \param jet  Jet containing the information to be sent to the tree/histograms
void AliAnalysisTaskEmcalJetSpectraQA::FillTTree(const AliEmcalJetInfo& jet, const AliJetContainer* jets)
{
  static TTree* tree = 0;

  if (!tree || TString(tree->GetName()) != jets->GetName()) tree = static_cast<TTree*>(fTreeList->FindObject(jets->GetName()));
  if (!tree) return;

  fCurrentJetInfo->Set(jet);

  tree->Fill();
}

/// Fill histogram or tree with jet
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
