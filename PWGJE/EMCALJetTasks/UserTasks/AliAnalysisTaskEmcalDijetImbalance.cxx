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
#include <THnSparse.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskEmcalDijetImbalance.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalDijetImbalance);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalDijetImbalance::AliAnalysisTaskEmcalDijetImbalance() : 
  AliAnalysisTaskEmcalJet(),
  fHistManager(),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fDeltaPhiMin(0),
  fNDijetPtThresholds(1),
  fMinTrigJetPt(0),
  fMinAssJetPt(0),
  fDijetLeadingHadronPt(0),
  fMaxPt(250),
  fNCentHistBins(0),
  fCentHistBins(0),
  fPlotJetHistograms(kFALSE),
  fPlotDijetJetHistograms(kFALSE),
  fPlotDijetImbalanceHistograms(kFALSE),
  fDoMomentumBalance(kFALSE),
  fDoGeometricalMatching(kFALSE),
  fMatchingJetR(0.2),
  fTrackConstituentThreshold(0),
  fClusterConstituentThreshold(0)
{
  GenerateHistoBins();
  Dijet_t fDijet;
  Dijet_t fMatchingDijet;
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalDijetImbalance::AliAnalysisTaskEmcalDijetImbalance(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistManager(name),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fDeltaPhiMin(0),
  fNDijetPtThresholds(1),
  fMinTrigJetPt(0),
  fMinAssJetPt(0),
  fDijetLeadingHadronPt(0),
  fMaxPt(250),
  fNCentHistBins(0),
  fCentHistBins(0),
  fPlotJetHistograms(kFALSE),
  fPlotDijetJetHistograms(kFALSE),
  fPlotDijetImbalanceHistograms(kFALSE),
  fDoMomentumBalance(kFALSE),
  fDoGeometricalMatching(kFALSE),
  fMatchingJetR(0.2),
  fTrackConstituentThreshold(0),
  fClusterConstituentThreshold(0)
{
  GenerateHistoBins();
  Dijet_t fDijet;
  Dijet_t fMatchingDijet;
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalDijetImbalance::~AliAnalysisTaskEmcalDijetImbalance()
{
}

/**
 * Generate histogram binning arrays
 */
void AliAnalysisTaskEmcalDijetImbalance::GenerateHistoBins()
{
  fNCentHistBins = 4;
  fCentHistBins = new Double_t[fNCentHistBins+1];
  fCentHistBins[0] = 0;
  fCentHistBins[1] = 10;
  fCentHistBins[2] = 30;
  fCentHistBins[3] = 50;
  fCentHistBins[4] = 90;
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalDijetImbalance::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  if (fPlotJetHistograms) AllocateJetHistograms();
  if (fPlotDijetJetHistograms) AllocateDijetJetHistograms();
  if (fPlotDijetImbalanceHistograms) AllocateDijetImbalanceHistograms();
  if (fDoMomentumBalance) AllocateMomentumBalanceHistograms();
  if (fDoGeometricalMatching) AllocateGeometricalMatchingHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  // Intialize AliEventCuts
  fEventCutList = new TList();
  fEventCutList ->SetOwner();
  fEventCutList ->SetName("EventCutOutput");
  
  fEventCuts.OverrideAutomaticTriggerSelection(fOffTrigger);
  if(fUseManualEventCuts==1)
  {
    fEventCuts.SetManualMode();
    // Configure manual settings here
    // ...
  }
  fEventCuts.AddQAplotsToList(fEventCutList);
  fOutput->Add(fEventCutList);
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for single jets.
 * A set of histograms is allocated per each jet container and per each centrality bin.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateJetHistograms()
{
  
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    // Allocate THnSparse
    TString axisTitle[30]= {""};
    Int_t nbins[30]  = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    Double_t *binEdges[20] = {0};
    Int_t dim = 0;
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "Centrality (%)";
      nbins[dim] = fNCentHistBins;
      binEdges[dim] = fCentHistBins;
      min[dim] = fCentHistBins[0];
      max[dim] = fCentHistBins[fNCentHistBins];
      dim++;
    }
    
    axisTitle[dim] = "#eta_{jet}";
    nbins[dim] = 100;
    min[dim] = -1;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#phi_{jet} (rad)";
    nbins[dim] = 101;
    min[dim] = 0;
    max[dim] = TMath::Pi() * 2.02;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;

    axisTitle[dim] = "#it{p}_{T} (GeV/#it{c})";
    nbins[dim] = TMath::CeilNint(fMaxPt/2);
    min[dim] = 0;
    max[dim] = fMaxPt;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "#it{p}_{T}^{corr} (GeV/#it{c})";
      nbins[dim] = TMath::CeilNint(fMaxPt/2);
      min[dim] = -fMaxPt/2 + 25;
      max[dim] = fMaxPt/2 + 25;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
    }
    
    axisTitle[dim] = "#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    nbins[dim] = 56;
    Double_t* pTLeadingBins = new Double_t[nbins[dim]+1];
    GenerateFixedBinArray(20, 0, 10, pTLeadingBins);
    GenerateFixedBinArray(10, 10, 20, pTLeadingBins+20);
    GenerateFixedBinArray(26, 20, 150, pTLeadingBins+30);
    min[dim] = 0;
    max[dim] = pTLeadingBins[nbins[dim]];
    binEdges[dim] = pTLeadingBins;
    dim++;
    
    axisTitle[dim] = "#it{A}_{jet}";
    nbins[dim] = fMaxPt/3;
    min[dim] = 0;
    max[dim] = 1.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    if (fClusterCollArray.GetEntriesFast() > 0 && fParticleCollArray.GetEntriesFast() > 0) {
      axisTitle[dim] = "NEF";
      nbins[dim] = fMaxPt/5;
      min[dim] = 0;
      max[dim] = 1.0;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
    }
    
    axisTitle[dim] = "#it{z}_{leading}";
    nbins[dim] = fMaxPt/5;
    min[dim] = 0;
    max[dim] = 1.0;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "No. of constituents";
      nbins[dim] = fMaxPt/5;
      min[dim] = 0;
      max[dim] = fMaxPt;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
    }
    else {
      axisTitle[dim] = "No. of constituents";
      nbins[dim] = 50;
      min[dim] = -0.5;
      max[dim] = 49.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
    }
    
    TString thnname = TString::Format("%s/JetObservables", jets->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(axisTitle[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
    
    // Allocate other jet histograms
    TString histname;
    TString title;
    histname = TString::Format("%s/fHistJetRejectionReason", jets->GetArrayName().Data());
    title = histname + ";Rejection reason;#it{p}_{T,jet} (GeV/#it{c});counts";
    TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 50, 0, 250);
    SetRejectionReasonLabels(hist->GetXaxis());
    
    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/fHistRhoVsCent", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 101, 0, 101, 100, 0, 500);
    }
    
  }
}

/*
 * This function allocates the histograms for single jets that comprise dijet pairs.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateDijetJetHistograms()
{
  // Allocate dijet THnSparse
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    TString axisTitle[30]= {""};
    Int_t nbins[30]  = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    Double_t *binEdges[20] = {0};
    Int_t dim = 0;
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "Centrality (%)";
      nbins[dim] = fNCentHistBins;
      binEdges[dim] = fCentHistBins;
      min[dim] = fCentHistBins[0];
      max[dim] = fCentHistBins[fNCentHistBins];
      dim++;
    }
    
    axisTitle[dim] = "LeadingHadronRequired";
    nbins[dim] = 2;
    min[dim] = -0.5;
    max[dim] = 1.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,min}^{trig}";
    nbins[dim] = fNDijetPtThresholds;
    min[dim] = -0.5;
    max[dim] = fNDijetPtThresholds - 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,min}^{ass}";
    nbins[dim] = fNDijetPtThresholds;
    min[dim] = -0.5;
    max[dim] = fNDijetPtThresholds - 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "isAccepted";
    nbins[dim] = 2;
    min[dim] = -0.5;
    max[dim] = 1.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "LeadingSubleading";
    nbins[dim] = 2;
    min[dim] = -0.5;
    max[dim] = 1.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,jet} (GeV/#it{c})";
    nbins[dim] = fMaxPt/3;
    min[dim] = -fMaxPt/2 + 25;
    max[dim] = fMaxPt/2 + 25;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#phi_{jet}";
    nbins[dim] = fMaxPt/3;
    min[dim] = 0;
    max[dim] = TMath::TwoPi();
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#eta_{jet}";
    nbins[dim] = fMaxPt/3;
    min[dim] = -1;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "N_{tracks}";
    nbins[dim] = fMaxPt/5;
    min[dim] = 0;
    max[dim] = 100;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "A_{jet}";
    nbins[dim] = fMaxPt/3;
    min[dim] = 0;
    max[dim] = 1.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    TString thnname = TString::Format("%s/DijetObservables", jets->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(axisTitle[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
  }
}

/*
 * This function allocates the histograms for dijet imbalance.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateDijetImbalanceHistograms()
{
  // Allocate dijet imbalance THnSparse
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    TString axisTitle[30]= {""};
    Int_t nbins[30]  = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    Double_t *binEdges[20] = {0};
    Int_t dim = 0;
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "Centrality (%)";
      nbins[dim] = fNCentHistBins;
      binEdges[dim] = fCentHistBins;
      min[dim] = fCentHistBins[0];
      max[dim] = fCentHistBins[fNCentHistBins];
      dim++;
    }
    
    axisTitle[dim] = "LeadingHadronRequired";
    nbins[dim] = 2;
    min[dim] = -0.5;
    max[dim] = 1.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,min}^{trig}";
    nbins[dim] = fNDijetPtThresholds;
    min[dim] = -0.5;
    max[dim] = fNDijetPtThresholds - 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,min}^{ass}";
    nbins[dim] = fNDijetPtThresholds;
    min[dim] = -0.5;
    max[dim] = fNDijetPtThresholds - 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;

    axisTitle[dim] = "#Delta#phi";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 4;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#Delta#eta";
    nbins[dim] = 100;
    min[dim] = -2;
    max[dim] = 2;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "A_{J}";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "x_{J}";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
  
    axisTitle[dim] = "k_{Ty} (GeV)";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 100;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;

    TString thnname = TString::Format("%s/DijetImbalanceObservables", jets->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(axisTitle[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
  }
}

/*
 * This function allocates the histograms for the momentum balance study.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateMomentumBalanceHistograms()
{
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    // Allocate THnSparse
    TString axisTitle[30]= {""};
    Int_t nbins[30]  = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    Double_t *binEdges[20] = {0};
    Int_t dim = 0;
    
    axisTitle[dim] = "LeadingHadronRequired";
    nbins[dim] = 2;
    min[dim] = -0.5;
    max[dim] = 1.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,min}^{trig}";
    nbins[dim] = fNDijetPtThresholds;
    min[dim] = -0.5;
    max[dim] = fNDijetPtThresholds - 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,min}^{ass}";
    nbins[dim] = fNDijetPtThresholds;
    min[dim] = -0.5;
    max[dim] = fNDijetPtThresholds - 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "A_{J}";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#Delta#phi";
    nbins[dim] = 100;
    min[dim] = -4;
    max[dim] = 4;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,particle} (GeV/#it{c})";
    nbins[dim] = 9;
    Double_t* pTParticleBins = new Double_t[nbins[dim]+1];
    GenerateFixedBinArray(1, 0.15, 0.3, pTParticleBins);
    GenerateFixedBinArray(1, 0.3, 0.5, pTParticleBins+1);
    GenerateFixedBinArray(1, 0.5, 1, pTParticleBins+2);
    GenerateFixedBinArray(2, 1, 5, pTParticleBins+3);
    GenerateFixedBinArray(3, 5, 20, pTParticleBins+5);
    GenerateFixedBinArray(1, 20, 150, pTParticleBins+8);
    min[dim] = 0;
    max[dim] = pTParticleBins[nbins[dim]];
    binEdges[dim] = pTParticleBins;
    dim++;
    
    axisTitle[dim] = "#it{p}_{T}#parallel (GeV/#it{c})";
    nbins[dim] = 80;
    Double_t* pTParallelBins = new Double_t[nbins[dim]+1];
    GenerateFixedBinArray(20, 0, 2, pTParallelBins);
    GenerateFixedBinArray(16, 2, 10, pTParallelBins+20);
    GenerateFixedBinArray(10, 10, 20, pTParallelBins+36);
    GenerateFixedBinArray(10, 20, 40, pTParallelBins+46);
    GenerateFixedBinArray(24, 40, 150, pTParallelBins+56);
    min[dim] = 0;
    max[dim] = pTParallelBins[nbins[dim]];
    binEdges[dim] = pTParallelBins;
    dim++;
    
    TString thnname = TString::Format("%s/MomentumBalance", jets->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(axisTitle[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
  }
}

/*
 * This function allocates the histograms for the constituent dijet study.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateGeometricalMatchingHistograms()
{
  // Allocate geometrical matching THnSparse
    TString axisTitle[30]= {""};
    Int_t nbins[30]  = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    Double_t *binEdges[20] = {0};
    Int_t dim = 0;
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "Centrality (%)";
      nbins[dim] = fNCentHistBins;
      binEdges[dim] = fCentHistBins;
      min[dim] = fCentHistBins[0];
      max[dim] = fCentHistBins[fNCentHistBins];
      dim++;
    }
    
    axisTitle[dim] = "#it{p}_{T,min}^{trig}";
    nbins[dim] = fNDijetPtThresholds;
    min[dim] = -0.5;
    max[dim] = fNDijetPtThresholds - 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,min}^{ass}";
    nbins[dim] = fNDijetPtThresholds;
    min[dim] = -0.5;
    max[dim] = fNDijetPtThresholds - 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
  
    axisTitle[dim] = "isSwitched";
    nbins[dim] = 2;
    min[dim] = -0.5;
    max[dim] = 1.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
  
    axisTitle[dim] = "#DeltaR_{trig}";
    nbins[dim] = 50;
    min[dim] = 0;
    max[dim] = 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#DeltaR_{ass}";
    nbins[dim] = 50;
    min[dim] = 0;
    max[dim] = 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
  
    axisTitle[dim] = "trig #it{p}_{T,low-thresh} - #it{p}_{T,hard-core}";
    nbins[dim] = 100;
    min[dim] = -50;
    max[dim] = 50;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
  
    axisTitle[dim] = "ass #it{p}_{T,low-thresh} - #it{p}_{T,hard-core}";
    nbins[dim] = 100;
    min[dim] = -50;
    max[dim] = 50;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
  
    axisTitle[dim] = "A_{J} low-threshold";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "A_{J} hard-core";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    TString thnname = "GeometricalMatching";
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(axisTitle[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
  
  // Allocate other histograms
  TString histname;
  TString title;
  histname = "GeometricalMatchingEfficiency";
  title = histname + ";isMatched;counts";
  TH1* hist = fHistManager.CreateTH1(histname.Data(), title.Data(), 2, -0.5, 1.5);
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalDijetImbalance::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
  
  AliInfo(Form("Number of pT thresholds: %d", fNDijetPtThresholds));
  for (Int_t i = 0; i < fNDijetPtThresholds; i++) {
    AliInfo(Form("Trigger jet threshold %d = %f, Associated jet threshold %d = %f", i, fMinTrigJetPt[i], i, fMinAssJetPt[i]));
  }
  AliInfo(Form("Leading hadron threshold (for dijet leading jet): %f GeV", fDijetLeadingHadronPt));
  AliInfo(Form("Momentum balance study: %d", fDoMomentumBalance));
  AliInfo(Form("Geometrical matching study: %d", fDoGeometricalMatching));
  
}

/**
 * This function (overloading the base class) uses AliEventCuts to perform event selection.
 */
Bool_t AliAnalysisTaskEmcalDijetImbalance::IsEventSelected()
{
  if (!fEventCuts.AcceptEvent(InputEvent()))
  {
    PostData(1, fOutput);
    return kFALSE;
  }
  return kTRUE;
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalDijetImbalance::Run()
{
  TString histname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    TString jetContName = jetCont->GetName();
    if (jetContName.Contains("HardCore")) continue;

    // Loop over leading hadron cut or not
    for (Int_t leadingHadronCutType=0; leadingHadronCutType<2; leadingHadronCutType++) {
      
      // Loop through leading jet pT thresholds
      for (Int_t trigJetMinPtType = 0; trigJetMinPtType < fNDijetPtThresholds; trigJetMinPtType++) {
        
        // Loop through subleading jet pT thresholds
        for (Int_t assJetMinPtType=0; assJetMinPtType < fNDijetPtThresholds; assJetMinPtType++) {
          
          // Find the dijet candidate of the event and store its info in struct fDijet
          FindDijet(jetCont, leadingHadronCutType, trigJetMinPtType, assJetMinPtType);

          // If we find a dijet candidate (i.e. acceptable trig jet; ass jet accepted or not), fill the di-jet jet histograms
          if (fDijet.trigJet && fPlotDijetJetHistograms) {
            TString histname = TString::Format("%s/DijetObservables", jetCont->GetArrayName().Data());
            FillDijetJetHistograms(histname, fDijet.isAccepted, 0, fDijet.trigJetPt, fDijet.trigJetPhi, fDijet.trigJetEta,
                                   fDijet.trigJet->GetNumberOfTracks(), fDijet.trigJet->Area());
            if (fDijet.assJet)
              FillDijetJetHistograms(histname, fDijet.isAccepted, 1, fDijet.assJetPt, fDijet.assJetPhi, fDijet.assJetEta,
                                     fDijet.assJet->GetNumberOfTracks(), fDijet.assJet->Area());
          }
          
          // If we find an accepted dijet, fill the dijet imbalance histogram
          if (fDijet.isAccepted && fPlotDijetImbalanceHistograms) {
            TString histname = TString::Format("%s/DijetImbalanceObservables", jetCont->GetArrayName().Data());
            FillDijetImbalanceHistograms(histname);
          }

          // If we find an accepted dijet, perform momentum-balance study (if requested)
          if (fDijet.isAccepted && fDoMomentumBalance) {
            histname = TString::Format("%s/MomentumBalance", jetCont->GetArrayName().Data());
            DoMomentumBalance(histname);
          }
        }
      }
    }
  }

  // Do the constituent threshold and geometrical matching study (if requested)
  if (fDoGeometricalMatching)
    DoGeometricalMatching();

  return kTRUE;
}

/**
 * Find the leading dijet in an event (background subtracted, unless hard-core jet container).
 * Fills dijet to fDijet.
 */
void AliAnalysisTaskEmcalDijetImbalance::FindDijet(AliJetContainer* jetCont, Int_t leadingHadronCutBin, Int_t trigJetMinPtBin, Int_t assJetMinPtBin)
{
  fDijet.clear();
  fDijet.leadingHadronCutType = leadingHadronCutBin;
  fDijet.trigJetMinPtType = trigJetMinPtBin;
  fDijet.assJetMinPtType = assJetMinPtBin;
  
  // Get trigger jet
  AliEmcalJet* trigJet = 0;
  if (jetCont->GetRhoParameter())
    trigJet = jetCont->GetLeadingJet("rho");
  else
    trigJet = jetCont->GetLeadingJet();
  if(!trigJet) return;
  
  // Skip the event if the leading jet doesn't satisfy the pT threshold
  Double_t trigJetPt = GetJetPt(jetCont, trigJet);
  if ( trigJetPt < fMinTrigJetPt[trigJetMinPtBin] ) return;
  
  // Skip the event if the leading jet doesn't satisfy the leading hadron threshold
  if (jetCont->GetLeadingHadronPt(trigJet) < fDijetLeadingHadronPt*leadingHadronCutBin) return;
  
  // Fill the dijet struct
  fDijet.trigJet = trigJet;
  fDijet.trigJetPt = trigJetPt;
  fDijet.trigJetEta = trigJet->Eta();
  fDijet.trigJetPhi = trigJet->Phi();
  
  // Find the subleading jet in the opposite hemisphere
  AliEmcalJet *assJet = 0;
  for(auto assJetCand : jetCont->accepted()) {
    if (!assJetCand) continue;
    Double_t assJetCandPt = GetJetPt(jetCont, assJetCand);
    if ( TMath::Abs(trigJet->Phi() - assJetCand->Phi()) < fDeltaPhiMin ) continue;
    if (assJet) {
      Double_t assJetPt = GetJetPt(jetCont, assJet);
      if ( assJetCandPt < assJetPt ) continue;
    }
    assJet = assJetCand;
  }
  if (!assJet) return;
  
  // Fill the dijet struct
  fDijet.assJet = assJet;
  fDijet.assJetPt = GetJetPt(jetCont, assJet);
  fDijet.assJetPhi  = assJet->Phi();
  fDijet.assJetEta  = assJet->Eta();
  fDijet.isAccepted = fDijet.assJetPt > fMinAssJetPt[assJetMinPtBin];
  
  fDijet.deltaPhi = TMath::Abs(trigJet->Phi() - assJet->Phi());
  fDijet.deltaEta = trigJet->Eta() - assJet->Eta();
  fDijet.AJ = (fDijet.trigJetPt - fDijet.assJetPt)/(fDijet.trigJetPt + fDijet.assJetPt);
  fDijet.xJ = fDijet.assJetPt / fDijet.trigJetPt;
  fDijet.kTy = TMath::Abs( fDijet.trigJetPt * TMath::Sin(fDijet.deltaPhi) );
}

/**
 * Do momentum balance study.
 * Loop through tracks in event, and plot them relative to the dijet axes.
 */
void AliAnalysisTaskEmcalDijetImbalance::DoMomentumBalance(TString histname)
{
  
  AliTrackContainer* trackCont = dynamic_cast<AliTrackContainer*>(GetParticleContainer("tracks"));
  
  AliVTrack* track;
  for (auto trackIterator : trackCont->accepted_momentum() ) {
    
    track = trackIterator.second;
    
    // Compute the delta phi between the track and its nearest jet (of the two jets in the dijet),
    // as well as its pT-parallel projection onto the nearest jet's axis.
    
    Double_t trackPt = track->Pt();
    Double_t trackPhi = track->Phi();
    Double_t trigJetPhi = fDijet.trigJet->Phi();
    Double_t assJetPhi = fDijet.assJet->Phi();
    
    Double_t deltaPhiTrigJet = TMath::Abs(trackPhi - trigJetPhi);
    Double_t deltaPhiAssJet = TMath::Abs(trackPhi - assJetPhi);
    Bool_t isNearside = deltaPhiTrigJet < deltaPhiAssJet;
    
    Double_t deltaPhi;
    Double_t balancePt;
    if (isNearside) {
      deltaPhi = trackPhi - trigJetPhi;
      balancePt = trackPt * TMath::Cos(deltaPhi);
    }
    else {
      deltaPhi = trackPhi - assJetPhi;
      balancePt = -trackPt * TMath::Cos(deltaPhi);
    }
    
    FillMomentumBalanceHistograms(histname, deltaPhi, trackPt, balancePt);
    
  }
}

/**
 * Do the constituent threshold and geometrical matching study.
 */
void AliAnalysisTaskEmcalDijetImbalance::DoGeometricalMatching()
{
  // Get jet container with minimum constituent pT,E thresholds
  TString jetContAllName = Form("Jet_AKTFullR0%d0_tracks_pT0150_caloClusters_E0300_pt_scheme", (int) (fMatchingJetR*10) );
  AliJetContainer* jetContAll = GetJetContainer(jetContAllName.Data());
  
  // Get jet container with X GeV constituent pT,E thresholds
  Int_t trackThreshold = (int) (fTrackConstituentThreshold*1000); // in MeV
  Int_t clusThreshold = (int) (fClusterConstituentThreshold*1000); // in MeV
  TString jetContHardCoreName = Form("JetHardCore_AKTFullR0%d0_tracks_pT%d_caloClusters_E%d_pt_scheme", (int) (fMatchingJetR*10), trackThreshold, clusThreshold);
  AliJetContainer* jetContHardCore = GetJetContainer(jetContHardCoreName.Data());
  
  // Find the di-jet in the hard-core jet sample, then find the matching di-jet and fill histograms
  FindDijet(jetContHardCore, 0, 0, 0);
  if (fDijet.isAccepted) {
    FindMatchingDijet(jetContAll, 0);
    FillGeometricalMatchingHistograms();
  }
}

/**
 * Find the matching leading dijet in an event (background subtracted, unless hard-core jet container).
 * Fills matched dijet to fMatchingDijet.
 */
void AliAnalysisTaskEmcalDijetImbalance::FindMatchingDijet(AliJetContainer* jetCont, Int_t assJetMinPtBin)
{
  fMatchingDijet.clear();
  fMatchingDijet.assJetMinPtType = assJetMinPtBin;
  
  // Loop over jets and find leading jet within R of fDijet.trigJet
  AliEmcalJet *matchingTrigJet = 0;
  for(auto matchingTrigJetCand : jetCont->accepted()) {
    if (!matchingTrigJetCand) continue;
    if ( GetDeltaR(matchingTrigJetCand, fDijet.trigJet) > fMatchingJetR ) continue;
    if (matchingTrigJet) {
      if ( GetJetPt(jetCont, matchingTrigJetCand) < GetJetPt(jetCont, matchingTrigJet) ) continue;
    }
    matchingTrigJet = matchingTrigJetCand;
  }
  if (!matchingTrigJet) return;
  
  // Loop over jets and find leading jet within R of fDijet.assJet
  AliEmcalJet *matchingAssJet = 0;
  for(auto matchingAssJetCand : jetCont->accepted()) {
    if (!matchingAssJetCand) continue;
    if ( GetDeltaR(matchingAssJetCand, fDijet.assJet) > fMatchingJetR ) continue;
    if (matchingAssJet) {
      if ( GetJetPt(jetCont, matchingAssJetCand) < GetJetPt(jetCont, matchingAssJet) ) continue;
    }
    matchingAssJet = matchingAssJetCand;
  }
  
  // Determine which matching jet is the leading jet (i.e. allow them to flip)
  if (matchingAssJet) {
    AliEmcalJet* trigJet = matchingTrigJet;
    AliEmcalJet* assJet = matchingAssJet;
    if ( GetJetPt(jetCont, matchingTrigJet) < GetJetPt(jetCont, matchingAssJet) ) {
      AliEmcalJet* trigJet = matchingAssJet;
      AliEmcalJet* assJet = matchingTrigJet;
    }
    
    // Fill the dijet struct
    fMatchingDijet.trigJet = trigJet;
    fMatchingDijet.trigJetPt = GetJetPt(jetCont, trigJet);
    fMatchingDijet.trigJetEta = trigJet->Eta();
    fMatchingDijet.trigJetPhi = trigJet->Phi();
    
    fMatchingDijet.assJet = assJet;
    fMatchingDijet.assJetPt = GetJetPt(jetCont, assJet);
    fMatchingDijet.assJetPhi  = assJet->Phi();
    fMatchingDijet.assJetEta  = assJet->Eta();
    fMatchingDijet.isAccepted = fMatchingDijet.assJetPt > fMinAssJetPt[assJetMinPtBin];
    
    fMatchingDijet.deltaPhi = TMath::Abs(trigJet->Phi() - assJet->Phi());
    fMatchingDijet.deltaEta = trigJet->Eta() - assJet->Eta();
    fMatchingDijet.AJ = (fMatchingDijet.trigJetPt - fMatchingDijet.assJetPt)/(fMatchingDijet.trigJetPt + fMatchingDijet.assJetPt);
    fMatchingDijet.xJ = fMatchingDijet.assJetPt / fMatchingDijet.trigJetPt;
    fMatchingDijet.kTy = TMath::Abs( fMatchingDijet.trigJetPt * TMath::Sin(fMatchingDijet.deltaPhi) );
  }
}

/**
 * Get pT of jet -- background subtracted, unless hard-core jet
 */
Double_t AliAnalysisTaskEmcalDijetImbalance::GetJetPt(AliJetContainer* jetCont, AliEmcalJet* jet)
{
  Double_t pT = jet->Pt() - jetCont->GetRhoVal() * jet->Area();
  TString jetContName = jetCont->GetName();
  if (jetContName.Contains("HardCore")) pT = jet->Pt();
  return pT;
}

/**
 * Get deltaR of two jets
 */
Double_t AliAnalysisTaskEmcalDijetImbalance::GetDeltaR(AliEmcalJet* jet1, AliEmcalJet* jet2)
{
  Double_t deltaPhi = TMath::Abs(jet1->Phi() - jet2->Phi());
  Double_t deltaEta = TMath::Abs(jet1->Eta() - jet2->Eta());
  Double_t deltaR = TMath::Sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta );
  return deltaR;
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalDijetImbalance::FillHistograms()
{
  if (fPlotJetHistograms) FillJetHistograms();
  
  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalDijetImbalance::FillJetHistograms()
{
  TString histname;
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    TString jetContName = jets->GetName();
    if (jetContName.Contains("HardCore")) continue;
    Double_t rhoVal = 0;
    if (jets->GetRhoParameter()) {
      rhoVal = jets->GetRhoVal();
      histname = TString::Format("%s/fHistRhoVsCent", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), fCent, rhoVal);
    }
    
    for (auto jet : jets->all()) {
      
      UInt_t rejectionReason = 0;
      if (!jets->AcceptJet(jet, rejectionReason)) {
        histname = TString::Format("%s/fHistJetRejectionReason", jets->GetArrayName().Data());
        fHistManager.FillTH2(histname.Data(), jets->GetRejectionReasonBitPosition(rejectionReason), jet->Pt());
        continue;
      }
      
      Float_t ptLeading = jets->GetLeadingHadronPt(jet);
      Float_t corrPt = jet->Pt() - rhoVal * jet->Area();
      
      TLorentzVector leadPart;
      jets->GetLeadingHadronMomentum(leadPart, jet);
      Double_t z = GetParallelFraction(leadPart.Vect(), jet);
      if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
      
      Double_t contents[30]={0};
      histname = TString::Format("%s/JetObservables", jets->GetArrayName().Data());
      THnSparse* histJetObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));
      if (!histJetObservables) return;
      for (Int_t i = 0; i < histJetObservables->GetNdimensions(); i++) {
        TString title(histJetObservables->GetAxis(i)->GetTitle());
        if (title=="Centrality (%)")
          contents[i] = fCent;
        else if (title=="#eta_{jet}")
          contents[i] = jet->Eta();
        else if (title=="#phi_{jet} (rad)")
          contents[i] = jet->Phi_0_2pi();
        else if (title=="#it{p}_{T} (GeV/#it{c})")
          contents[i] = jet->Pt();
        else if (title=="#it{p}_{T}^{corr} (GeV/#it{c})")
          contents[i] = corrPt;
        else if (title=="#it{p}_{T,particle}^{leading} (GeV/#it{c})")
          contents[i] = ptLeading;
        else if (title=="#it{A}_{jet}")
          contents[i] = jet->Area();
        else if (title=="NEF")
          contents[i] = jet->NEF();
        else if (title=="#it{z}_{leading}")
          contents[i] = z;
        else if (title=="No. of constituents")
          contents[i] = jet->GetNumberOfConstituents();
        else
          AliWarning(Form("Unable to fill dimension %s!",title.Data()));
      }
      histJetObservables->Fill(contents);
      
    } //jet loop
  }
}

/**
 * Fill dijet THnSparse.
 */
void AliAnalysisTaskEmcalDijetImbalance::FillDijetJetHistograms(TString histname, Int_t isAccepted, Int_t IsAssJet, Double_t jetPt, Double_t jetPhi, Double_t jetEta, Int_t nTracksJet, Double_t jetArea)
{
  Double_t contents[30]={0};
  THnSparse* histJetObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));
  if (!histJetObservables) return;
  for (Int_t n = 0; n < histJetObservables->GetNdimensions(); n++) {
    TString title(histJetObservables->GetAxis(n)->GetTitle());
    if (title=="Centrality (%)")
      contents[n] = fCent;
    else if (title=="LeadingHadronRequired")
      contents[n] = fDijet.leadingHadronCutType;
    else if (title=="#it{p}_{T,min}^{trig}")
      contents[n] = fDijet.trigJetMinPtType;
    else if (title=="#it{p}_{T,min}^{ass}")
      contents[n] = fDijet.assJetMinPtType;
    else if (title=="isAccepted")
      contents[n] = isAccepted;
    else if (title=="LeadingSubleading")
      contents[n] = IsAssJet;
    else if (title=="#it{p}_{T,jet} (GeV/#it{c})")
      contents[n] = jetPt;
    else if (title=="#phi_{jet}")
      contents[n] = jetPhi;
    else if (title=="#eta_{jet}")
      contents[n] = jetEta;
    else if (title=="N_{tracks}")
      contents[n] = nTracksJet;
    else if (title=="A_{jet}")
      contents[n] = jetArea;
    else
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }
  histJetObservables->Fill(contents);
}

/**
 * Fill dijet imbalance THnSparse.
 */
void AliAnalysisTaskEmcalDijetImbalance::FillDijetImbalanceHistograms(TString histname)
{
  Double_t contents[30]={0};
  THnSparse* histJetObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));
  if (!histJetObservables) return;
  for (Int_t n = 0; n < histJetObservables->GetNdimensions(); n++) {
    TString title(histJetObservables->GetAxis(n)->GetTitle());
    if (title=="Centrality (%)")
      contents[n] = fCent;
    else if (title=="LeadingHadronRequired")
      contents[n] = fDijet.leadingHadronCutType;
    else if (title=="#it{p}_{T,min}^{trig}")
      contents[n] = fDijet.trigJetMinPtType;
    else if (title=="#it{p}_{T,min}^{ass}")
      contents[n] = fDijet.assJetMinPtType;
    else if (title=="#Delta#phi")
      contents[n] = fDijet.deltaPhi;
    else if (title=="#Delta#eta")
      contents[n] = fDijet.deltaEta;
    else if (title=="A_{J}")
      contents[n] = fDijet.AJ;
    else if (title=="x_{J}")
      contents[n] = fDijet.xJ;
    else if (title=="k_{Ty} (GeV)")
      contents[n] = fDijet.kTy;
    else
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }
  histJetObservables->Fill(contents);
}

/**
 * Fill momentum balance histogram.
 */
void AliAnalysisTaskEmcalDijetImbalance::FillMomentumBalanceHistograms(TString histname, Double_t deltaPhi, Double_t trackPt, Double_t balancePt)
{
  Double_t contents[30]={0};
  THnSparse* histJetObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));
  if (!histJetObservables) return;
  for (Int_t n = 0; n < histJetObservables->GetNdimensions(); n++) {
    TString title(histJetObservables->GetAxis(n)->GetTitle());
    if (title=="LeadingHadronRequired")
      contents[n] = fDijet.leadingHadronCutType;
    else if (title=="#it{p}_{T,min}^{trig}")
      contents[n] = fDijet.trigJetMinPtType;
    else if (title=="#it{p}_{T,min}^{ass}")
      contents[n] = fDijet.assJetMinPtType;
    else if (title=="A_{J}")
      contents[n] = fDijet.AJ;
    else if (title=="#Delta#phi")
      contents[n] = deltaPhi;
    else if (title=="#it{p}_{T,particle} (GeV/#it{c})")
      contents[n] = trackPt;
    else if (title=="#it{p}_{T}#parallel (GeV/#it{c})")
      contents[n] = balancePt;
    else
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }
  histJetObservables->Fill(contents);
}

/**
 * Fill momentum balance histogram.
 */
void AliAnalysisTaskEmcalDijetImbalance::FillGeometricalMatchingHistograms()
{
  // Matching efficiency histogram
  TString histname = "GeometricalMatchingEfficiency";
  
  // If we have a matching di-jet, fill the geometrical matching histogram
  if (fMatchingDijet.assJet) {
    fHistManager.FillTH1(histname.Data(), 1);
    
    Double_t trigDeltaR = GetDeltaR(fDijet.trigJet, fMatchingDijet.trigJet);
    Double_t assDeltaR = GetDeltaR(fDijet.assJet, fMatchingDijet.assJet);
    Bool_t isSwitched = trigDeltaR > fMatchingJetR;
    
    TString thnname = "GeometricalMatching";
    Double_t contents[30]={0};
    THnSparse* histJetObservables = static_cast<THnSparse*>(fHistManager.FindObject(thnname.Data()));
    if (!histJetObservables) return;
    for (Int_t n = 0; n < histJetObservables->GetNdimensions(); n++) {
      TString title(histJetObservables->GetAxis(n)->GetTitle());
      if (title=="Centrality (%)")
        contents[n] = fCent;
      else if (title=="#it{p}_{T,min}^{trig}")
        contents[n] = fDijet.trigJetMinPtType;
      else if (title=="#it{p}_{T,min}^{ass}")
        contents[n] = fDijet.assJetMinPtType;
      else if (title=="isSwitched")
        contents[n] = isSwitched;
      else if (title=="#DeltaR_{trig}")
        contents[n] = trigDeltaR;
      else if (title=="#DeltaR_{ass}")
        contents[n] = assDeltaR;
      else if (title=="trig #it{p}_{T,low-thresh} - #it{p}_{T,hard-core}")
        contents[n] = fMatchingDijet.trigJetPt - fDijet.trigJetPt;
      else if (title=="ass #it{p}_{T,low-thresh} - #it{p}_{T,hard-core}")
        contents[n] = fMatchingDijet.assJetPt - fDijet.assJetPt;
      else if (title=="A_{J} low-threshold")
        contents[n] = fMatchingDijet.AJ;
      else if (title=="A_{J} hard-core")
        contents[n] = fDijet.AJ;
      else
        AliWarning(Form("Unable to fill dimension %s!",title.Data()));
    }
    histJetObservables->Fill(contents);
  }
  else {
    fHistManager.FillTH1(histname.Data(), 0.);
  }
}

