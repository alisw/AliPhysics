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
#include <TRandom3.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalDownscaleFactorsOCDB.h"

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
  fMinTrigJetPt(0),
  fMinAssJetPt(0),
  fDijetLeadingHadronPt(0),
  fMaxPt(250),
  fNCentHistBins(0),
  fCentHistBins(0),
  fPlotJetHistograms(kFALSE),
  fPlotDijetCandHistograms(kFALSE),
  fPlotDijetImbalanceHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoMomentumBalance(kFALSE),
  fDoGeometricalMatching(kFALSE),
  fMatchingJetR(0.2),
  fTrackConstituentThreshold(0),
  fClusterConstituentThreshold(0),
  fMBUpscaleFactor(1.)
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
  fMinTrigJetPt(0),
  fMinAssJetPt(0),
  fDijetLeadingHadronPt(0),
  fMaxPt(250),
  fNCentHistBins(0),
  fCentHistBins(0),
  fPlotJetHistograms(kFALSE),
  fPlotDijetCandHistograms(kFALSE),
  fPlotDijetImbalanceHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoMomentumBalance(kFALSE),
  fDoGeometricalMatching(kFALSE),
  fMatchingJetR(0.2),
  fTrackConstituentThreshold(0),
  fClusterConstituentThreshold(0),
  fMBUpscaleFactor(1.)
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
  if (fPlotDijetCandHistograms) AllocateDijetCandHistograms();
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
 * A set of histograms is allocated per each jet container.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateJetHistograms()
{
  TString histname;
  TString title;
  
  Int_t nPtBins = TMath::CeilNint(fMaxPt/2);
  Int_t minPtBin = -fMaxPt/2 + 25;
  Int_t maxPtBin = fMaxPt/2 + 25;
  
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    // Jet rejection reason
    histname = TString::Format("%s/JetHistograms/hJetRejectionReason", jets->GetArrayName().Data());
    title = histname + ";Rejection reason;#it{p}_{T,jet} (GeV/#it{c});counts";
    TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 50, 0, 250);
    SetRejectionReasonLabels(hist->GetXaxis());
    
    // Rho vs. Centrality
    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/JetHistograms/hRhoVsCent", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 101, 0, 101, 100, 0, 500);
    }
    
    // Centrality vs. pT
    histname = TString::Format("%s/JetHistograms/hCentVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});Centrality (%);counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, minPtBin, maxPtBin, 10, 0, 100);
    
    // pT vs. eta vs. phi
    histname = TString::Format("%s/JetHistograms/hPtVsEtaVsPhi", jets->GetArrayName().Data());
    title = histname + ";#eta_{jet} (rad);#phi_{jet} (rad);#it{p}_{T}^{corr} (GeV/#it{c})";
    fHistManager.CreateTH3(histname.Data(), title.Data(), 50, -0.5, 0.5, 101, 0, TMath::Pi() * 2.02, 75, 0, maxPtBin);
    
    // pT upscaled
    histname = TString::Format("%s/JetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c})";
    fHistManager.CreateTH1(histname.Data(), title.Data(), 75, 0, maxPtBin, "s");
    
    // pT-leading vs. pT
    histname = TString::Format("%s/JetHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, minPtBin, maxPtBin, 150, 0, 150);
    
    // A vs. pT
    histname = TString::Format("%s/JetHistograms/hAreaVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, minPtBin, maxPtBin, fMaxPt/3, 0, 1.5);
    
    // Cent vs. NEF vs. pT
    histname = TString::Format("%s/JetHistograms/hCentVsNEFVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});NEF;Centrality (%)";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nPtBins, minPtBin, maxPtBin, fMaxPt/5, 0, 1.0, 50, 0, 100);
    
    // Cent vs. z-leading (charged) vs. pT
    histname = TString::Format("%s/JetHistograms/hCentVsZLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading};Centrality (%)";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nPtBins, minPtBin, maxPtBin, fMaxPt/5, 0, 1.0, 50, 0, 100);
    
    // Cent vs. z (charged) vs. pT
    histname = TString::Format("%s/JetHistograms/hCentVsZVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{z};Centrality (%)";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nPtBins, minPtBin, maxPtBin, fMaxPt/5, 0, 1.0, 50, 0, 100);

    // Cent vs. Nconst vs. pT
    histname = TString::Format("%s/JetHistograms/hCentVsNConstVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents;Centrality (%)";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nPtBins, minPtBin, maxPtBin, fMaxPt/5, 0, fMaxPt, 50, 0, 100);
    
    // Allocate background subtraction histograms, if enabled
    if (fComputeBackground) {
      
      histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCal", jets->GetArrayName().Data());
      title = histname + ";Centrality;Scale factor;counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 5);
      
      histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCal", jets->GetArrayName().Data());
      title = histname + ";#delta#it{p}_{T} (GeV/#it{c});counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), 400, -100, 100);
      
      histname = TString::Format("%s/BackgroundHistograms/hScaleFactorDCalRegion", jets->GetArrayName().Data());
      title = histname + ";Centrality;Scale factor;Eta bin";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 50, 0, 100, 200, 0, 10, 28, -0.5, 27.5);
      
      histname = TString::Format("%s/BackgroundHistograms/hDeltaPtDCalRegion", jets->GetArrayName().Data());
      title = histname + ";#delta#it{p}_{T} (GeV/#it{c});Eta bin;counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 400, -100, 100, 20, -0.5, 19.5);
      
    }
    
  }
  
  // MB downscale factor histogram
  histname = "Trigger/hMBDownscaleFactor";
  title = histname + ";Downscale factor;counts";
  TH1* hist = fHistManager.CreateTH1(histname.Data(), title.Data(), 200, 0, 200);
  
}

/*
 * This function allocates the histograms for dijet candidates, i.e. dijet pairs with the leading jet
 * passing the trigger condition, but no condition on the subleading jet. In particular this histogram is
 * designed to study the kinematic selection of dijets.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateDijetCandHistograms()
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
    
    axisTitle[dim] = "#it{p}_{T,trig jet} (GeV/#it{c})";
    nbins[dim] = fMaxPt/3;
    min[dim] = -fMaxPt/2 + 25;
    max[dim] = fMaxPt/2 + 25;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,ass jet} (GeV/#it{c})";
    nbins[dim] = fMaxPt/3;
    min[dim] = -fMaxPt/2 + 25;
    max[dim] = fMaxPt/2 + 25;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#phi_{trig jet}";
    nbins[dim] = fMaxPt/3;
    min[dim] = 0;
    max[dim] = TMath::TwoPi();
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#phi_{ass jet}";
    nbins[dim] = fMaxPt/3;
    min[dim] = 0;
    max[dim] = TMath::TwoPi();
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#eta_{trig jet}";
    nbins[dim] = fMaxPt/3;
    min[dim] = -1;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#eta_{ass jet}";
    nbins[dim] = fMaxPt/3;
    min[dim] = -1;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    TString thnname = TString::Format("%s/DijetCandObservables", jets->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(axisTitle[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
  }
}

/*
 * This function allocates the histograms for accepted dijets.
 * The purpose is to study in detail the imbalance properties of the accepted dijets.
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
    
    axisTitle[dim] = "N_{tracks, trig jet}";
    nbins[dim] = fMaxPt/5;
    min[dim] = 0;
    max[dim] = 100;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "N_{tracks, ass jet}";
    nbins[dim] = fMaxPt/5;
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
  
  fNeedEmcalGeom = kTRUE;
  
  AliInfo(Form("Trigger jet threshold = %f, Associated jet threshold = %f", fMinTrigJetPt, fMinAssJetPt));
  AliInfo(Form("Leading hadron threshold (for dijet leading jet): %f GeV", fDijetLeadingHadronPt));
  AliInfo(Form("Momentum balance study: %d", fDoMomentumBalance));
  AliInfo(Form("Geometrical matching study: %d", fDoGeometricalMatching));
  
}

/**
 * This function is called automatically when the run number changes.
 */
void AliAnalysisTaskEmcalDijetImbalance::RunChanged(Int_t run){
  
  // Get the downscaling factors for MB triggers (to be used to calculate trigger efficiency)
  
  // Get instance of the downscale factor helper class
  AliEmcalDownscaleFactorsOCDB *downscaleOCDB = AliEmcalDownscaleFactorsOCDB::Instance();
  downscaleOCDB->SetRun(InputEvent()->GetRunNumber());
  
  // There are two possible min bias triggers for LHC15o
  TString triggerNameMB1 = "CINT7-B-NOPF-CENT";
  TString triggerNameMB2 = "CV0L7-B-NOPF-CENT";
  TString triggerNameJE = "CINT7EJ1-B-NOPF-CENTNOPMD";
  
  // Get the downscale factor for whichever MB trigger exists in the given run
  std::vector<TString> runtriggers = downscaleOCDB->GetTriggerClasses();
  Double_t downscalefactor;
  for (auto i : runtriggers) {
    if (i.EqualTo(triggerNameMB1) || i.EqualTo(triggerNameMB2)) {
      downscalefactor = downscaleOCDB->GetDownscaleFactorForTriggerClass(i.Data());
      break;
    }
  }
  
  // Store the inverse of the downscale factor, used later to weight the pT spectrum
  fMBUpscaleFactor = 1/downscalefactor;
  
  TString histname = "Trigger/hMBDownscaleFactor";
  fHistManager.FillTH1(histname.Data(), fMBUpscaleFactor);
  
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

    //-----------------------------------------------------------------------------
    // Find the leading di-jet candidate in each event, and if it satisfies the
    // trig jet pT threshold, then fill di-jet candidate histogram (regardless of ass jet).
    // The idea is to study the kinematic selections in post-processing.
    
    // Loop over leading hadron cut or not
    for (Int_t leadingHadronCutType=0; leadingHadronCutType<2; leadingHadronCutType++) {
      
      // Find the dijet candidate of the event and store its info in struct fDijet
      FindDijet(jetCont, leadingHadronCutType);

      // If we find a dijet candidate (i.e. acceptable trig jet; ass jet accepted or not), fill the di-jet candidate histogram
      if (fDijet.trigJet && fPlotDijetCandHistograms) {
        TString histname = TString::Format("%s/DijetCandObservables", jetCont->GetArrayName().Data());
        FillDijetCandHistograms(histname);
      }
      
    }
    
    //---------------------------------------------------------------------------------------------------
    // Now, study the accepted dijet selection -- specified by the trig/ass jet pT conditions
    
    // Find the dijet candidate of the event and store its info in struct fDijet
    FindDijet(jetCont, 0);
    
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

  //---------------------------------------------------------------------------
  // Do the constituent threshold and geometrical matching study (if requested)
  if (fDoGeometricalMatching)
    DoGeometricalMatching();

  return kTRUE;
}

/**
 * Find the leading dijet in an event (background subtracted, unless hard-core jet container).
 * The trig jet is required to be above pT threshold, or else empty dijet is returned.
 * The ass jet is the leading jet in the opposite hemisphere.
 * Fills dijet to fDijet.
 * The field fDijet.isAccepted is true if the ass jet is above its corresponding pT threshold.
 */
void AliAnalysisTaskEmcalDijetImbalance::FindDijet(AliJetContainer* jetCont, Int_t leadingHadronCutBin)
{
  fDijet.clear();
  fDijet.leadingHadronCutType = leadingHadronCutBin;
  
  // Get trigger jet
  AliEmcalJet* trigJet = 0;
  if (jetCont->GetRhoParameter())
    trigJet = jetCont->GetLeadingJet("rho");
  else
    trigJet = jetCont->GetLeadingJet();
  if(!trigJet) return;
  
  // Skip the event if the leading jet doesn't satisfy the pT threshold
  Double_t trigJetPt = GetJetPt(jetCont, trigJet);
  if ( trigJetPt < fMinTrigJetPt ) return;
  
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
  fDijet.isAccepted = fDijet.assJetPt > fMinAssJetPt;
  
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
  FindDijet(jetContHardCore, 0);
  if (fDijet.isAccepted) {
    FindMatchingDijet(jetContAll);
    FillGeometricalMatchingHistograms();
  }
}

/**
 * Find the matching leading dijet in an event (background subtracted, unless hard-core jet container).
 * Fills matched dijet to fMatchingDijet.
 */
void AliAnalysisTaskEmcalDijetImbalance::FindMatchingDijet(AliJetContainer* jetCont)
{
  fMatchingDijet.clear();
  
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
    fMatchingDijet.isAccepted = fMatchingDijet.assJetPt > fMinAssJetPt;
    
    fMatchingDijet.deltaPhi = TMath::Abs(trigJet->Phi() - assJet->Phi());
    fMatchingDijet.deltaEta = trigJet->Eta() - assJet->Eta();
    fMatchingDijet.AJ = (fMatchingDijet.trigJetPt - fMatchingDijet.assJetPt)/(fMatchingDijet.trigJetPt + fMatchingDijet.assJetPt);
    fMatchingDijet.xJ = fMatchingDijet.assJetPt / fMatchingDijet.trigJetPt;
    fMatchingDijet.kTy = TMath::Abs( fMatchingDijet.trigJetPt * TMath::Sin(fMatchingDijet.deltaPhi) );
  }
}

/**
 * This function performs a study of the heavy-ion background.
 */
void AliAnalysisTaskEmcalDijetImbalance::ComputeBackground(AliJetContainer* jetCont)
{
  // Loop over tracks and clusters in order to:
  //   (1) Compute scale factor for full jets
  //   (2) Compute delta-pT for full jets, with the random cone method
  // For both the scale factor and delta-pT, we compute only one histogram each for EMCal.
  // But for DCal, we bin in eta, in order to study DCal vs. PHOS vs. gap
  
  // Define the acceptance boundaries for the TPC and EMCal/DCal/PHOS
  Double_t etaTPC = 0.9;
  Double_t etaEMCal = 0.7;
  Double_t etaMinDCal = 0.22;
  Double_t etaMaxPHOS = 0.13;
  Double_t phiMinEMCal = fGeom->GetArm1PhiMin() * TMath::DegToRad(); // 80
  Double_t phiMaxEMCal = fGeom->GetEMCALPhiMax() * TMath::DegToRad(); // ~188
  Double_t phiMinDCal = fGeom->GetDCALPhiMin() * TMath::DegToRad(); // 260
  Double_t phiMaxDCal = fGeom->GetDCALPhiMax() * TMath::DegToRad(); // ~327 (1/3 SMs start at 320)
  Double_t phiMinPHOS = 250 * TMath::DegToRad();
  Double_t phiMaxPHOS = 320 * TMath::DegToRad();

  Double_t accTPC = 2 * etaTPC * 2 * TMath::Pi();
  Double_t accEMCal = 2 * etaEMCal * (phiMaxEMCal - phiMinEMCal);
  Double_t accDCalRegion = 2 * etaEMCal * (phiMaxDCal - phiMinDCal);
  
  // Define fiducial acceptances, to be used to generate random cones
  TRandom3* r = new TRandom3(0);
  Double_t jetR = jetCont->GetJetRadius();
  Double_t etaEMCalfid = etaEMCal - jetR;
  Double_t phiMinEMCalfid = phiMinEMCal + jetR;
  Double_t phiMaxEMCalfid = phiMaxEMCal - jetR;
  Double_t phiMinDCalRegionfid = phiMinDCal + jetR;
  Double_t phiMaxDCalRegionfid = phiMaxDCal - jetR;
  
  // Generate EMCal random cone eta-phi
  Double_t etaEMCalRC = r->Uniform(-etaEMCalfid, etaEMCalfid);
  Double_t phiEMCalRC = r->Uniform(phiMinEMCalfid, phiMaxEMCalfid);
  
  // Generate DCalRegion random cone eta-phi inside each eta slice (same phi used for all)
  Double_t etaStep = 0.05;
  const Int_t nEtaBinsSF = 28; // 2 * 0.7 / 0.05
  const Int_t nEtaBinsRC = 20; // 2 * 0.5 / 0.05
  
  Double_t phiDCalRC = r->Uniform(phiMinDCalRegionfid, phiMaxDCalRegionfid);
  Double_t etaDCalRC[nEtaBinsRC];
  Double_t etaMin;
  Double_t etaMax;
  for (Int_t bin=0; bin < nEtaBinsRC; bin++) {
    etaMin = -etaEMCalfid + bin*etaStep;
    etaMax = etaMin + etaStep;
    etaDCalRC[bin] = r->Uniform(etaMin, etaMax);
  }
  
  // Initialize the various sums to 0
  Double_t trackPtSumTPC = 0;
  Double_t trackPtSumEMCal = 0;
  Double_t trackPtSumEMCalRC = 0;
  Double_t clusESumEMCal = 0;
  Double_t clusESumEMCalRC = 0;
  Double_t trackPtSumDCal[nEtaBinsSF] = {0.};
  Double_t trackPtSumDCalRC[nEtaBinsRC] = {0.};
  Double_t clusESumDCal[nEtaBinsSF] = {0.};
  Double_t clusESumDCalRC[nEtaBinsRC] = {0.};
  
  // Loop over tracks. Sum the track pT:
  // (1) in the entire TPC, (2) in the EMCal, (3) in the EMCal random cone,
  // (4) in the DCalRegion at each eta, (5) in the DCalRegion random cone at each eta
  AliTrackContainer* trackCont = dynamic_cast<AliTrackContainer*>(GetParticleContainer("tracks"));
  AliTLorentzVector track;
  Double_t trackEta;
  Double_t trackPhi;
  Double_t trackPt;
  Double_t deltaR;
  for (auto trackIterator : trackCont->accepted_momentum() ) {
    
    track.Clear();
    track = trackIterator.first;
    trackEta = track.Eta();
    trackPhi = track.Phi_0_2pi();
    trackPt = track.Pt();

    // (1)
    if (TMath::Abs(trackEta) < etaTPC) {
      trackPtSumTPC += trackPt;
    }
    
    // (2)
    if (TMath::Abs(trackEta) < etaEMCal && trackPhi > phiMinEMCal && trackPhi < phiMaxEMCal) {
      trackPtSumEMCal += trackPt;
    }
    
    // (3)
    deltaR = GetDeltaR(&track, etaEMCalRC, phiEMCalRC);
    if (deltaR < jetR) {
      trackPtSumEMCalRC += trackPt;
    }
    
    // (4)
    for (Int_t bin=0; bin < nEtaBinsSF; bin++) {
      if (trackPhi > phiMinDCal && trackPhi < phiMaxDCal) {
        etaMin = -etaEMCal+ bin*etaStep;
        etaMax = etaMin + etaStep;
        if (trackEta > etaMin && trackEta < etaMax) {
          trackPtSumDCal[bin] += trackPt;
        }
      }
    }
    
    // (5)
    for (Int_t bin=0; bin < nEtaBinsRC; bin++) {
      deltaR = GetDeltaR(&track, etaDCalRC[bin], phiDCalRC);
      if (deltaR < jetR) {
        trackPtSumDCalRC[bin] += trackPt;
      }
    }
    
  }
  
  // Loop over clusters. Sum the cluster ET:
  // (1) in the EMCal, (2) in the EMCal random cone, (3) in the DCalRegion at each eta,
  // (4) in the DCalRegion random cone at each eta
  AliClusterContainer* clusCont = GetClusterContainer(0);
  AliTLorentzVector clus;
  Double_t clusEta;
  Double_t clusPhi;
  Double_t clusE;
  for (auto clusIterator : clusCont->accepted_momentum() ) {
   
    clus.Clear();
    clus = clusIterator.first;
    clusEta = clus.Eta();
    clusPhi = clus.Phi_0_2pi();
    clusE = clus.E();
    
    // (1)
    if (TMath::Abs(clusEta) < etaEMCal && clusPhi > phiMinEMCal && clusPhi < phiMaxEMCal) {
      clusESumEMCal += clusE;
    }
    
    // (2)
    deltaR = GetDeltaR(&clus, etaEMCalRC, phiEMCalRC);
    if (deltaR < jetR) {
      clusESumEMCalRC += clusE;
    }
    
    // (3)
    for (Int_t bin=0; bin < nEtaBinsSF; bin++) {
      if (clusPhi > phiMinDCal && clusPhi < phiMaxDCal) {
        etaMin = -etaEMCal+ bin*etaStep;
        etaMax = etaMin + etaStep;
        if (clusEta > etaMin && clusEta < etaMax) {
          clusESumDCal[bin] += clusE;
        }
      }
    }
    
    // (4)
    for (Int_t bin=0; bin < nEtaBinsRC; bin++) {
      deltaR = GetDeltaR(&clus, etaDCalRC[bin], phiDCalRC);
      if (deltaR < jetR) {
        clusESumDCalRC[bin] += clusE;
      }
    }
    
  }
  
  // Compute the scale factor for EMCal
  Double_t numerator = (trackPtSumEMCal + clusESumEMCal) / accEMCal;
  Double_t denominator = trackPtSumTPC / accTPC;
  Double_t scaleFactor = numerator / denominator;
  TString histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCal", jetCont->GetArrayName().Data());
  fHistManager.FillTH2(histname, fCent, scaleFactor);
  
  // Compute the scale factor for DCalRegion
  Double_t accDCalRegionBin = accDCalRegion / nEtaBinsSF;
  for (Int_t bin=0; bin < nEtaBinsSF; bin++) {
    numerator = (trackPtSumDCal[bin] + clusESumDCal[bin]) / accDCalRegionBin;
    scaleFactor = numerator / denominator;
    histname = TString::Format("%s/BackgroundHistograms/hScaleFactorDCalRegion", jetCont->GetArrayName().Data());
    fHistManager.FillTH3(histname, fCent, scaleFactor, bin);
  }
  
  // Compute delta pT for EMCal
  Double_t rho = jetCont->GetRhoVal();
  Double_t deltaPt = trackPtSumEMCalRC + clusESumEMCalRC - rho * TMath::Pi() * jetR * jetR;
  histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCal", jetCont->GetArrayName().Data());
  fHistManager.FillTH1(histname, deltaPt);
  
  // Compute delta pT for DCalRegion
  for (Int_t bin=0; bin < nEtaBinsRC; bin++) {
    deltaPt = trackPtSumDCalRC[bin] + clusESumDCalRC[bin] - rho * TMath::Pi() * jetR * jetR;
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtDCalRegion", jetCont->GetArrayName().Data());
    fHistManager.FillTH2(histname, deltaPt, bin);
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
 * Get deltaR of a track/cluster and a reference point.
 */
Double_t AliAnalysisTaskEmcalDijetImbalance::GetDeltaR(AliTLorentzVector* part, Double_t etaRef, Double_t phiRef)
{
  Double_t deltaPhi = TMath::Abs(part->Phi_0_2pi() - phiRef);
  Double_t deltaEta = TMath::Abs(part->Eta() - etaRef);
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
      histname = TString::Format("%s/JetHistograms/hRhoVsCent", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), fCent, rhoVal);
    }
    
    for (auto jet : jets->all()) {
      
      Float_t ptLeading = jets->GetLeadingHadronPt(jet);
      Float_t corrPt = jet->Pt() - rhoVal * jet->Area();
      
      // A vs. pT (fill before area cut)
      histname = TString::Format("%s/JetHistograms/hAreaVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), corrPt, jet->Area());
      
      
      // Rejection reason
      UInt_t rejectionReason = 0;
      if (!jets->AcceptJet(jet, rejectionReason)) {
        histname = TString::Format("%s/JetHistograms/hJetRejectionReason", jets->GetArrayName().Data());
        fHistManager.FillTH2(histname.Data(), jets->GetRejectionReasonBitPosition(rejectionReason), jet->Pt());
        continue;
      }

      // Centrality vs. pT
      histname = TString::Format("%s/JetHistograms/hCentVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), corrPt, fCent);
      
      // pT vs. eta vs. phi
      histname = TString::Format("%s/JetHistograms/hPtVsEtaVsPhi", jets->GetArrayName().Data());
      fHistManager.FillTH3(histname.Data(), jet->Eta(), jet->Phi_0_2pi(), corrPt);
      
      // pT un-downscaled
      histname = TString::Format("%s/JetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
      fHistManager.FillTH1(histname.Data(), corrPt, fMBUpscaleFactor);
      
      // pT-leading vs. pT
      histname = TString::Format("%s/JetHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), corrPt, ptLeading);
      
      // Cent vs. NEF vs. pT
      histname = TString::Format("%s/JetHistograms/hCentVsNEFVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH3(histname.Data(), corrPt, jet->NEF(), fCent);
      
      // Cent vs. z-leading (charged) vs. pT
      TLorentzVector leadPart;
      jets->GetLeadingHadronMomentum(leadPart, jet);
      Double_t z = GetParallelFraction(leadPart.Vect(), jet);
      if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
      histname = TString::Format("%s/JetHistograms/hCentVsZLeadingVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH3(histname.Data(), corrPt, z, fCent);
      
      // Cent vs. z (charged) vs. pT
      histname = TString::Format("%s/JetHistograms/hCentVsZVsPt", jets->GetArrayName().Data());
      AliVTrack* track;
      for (Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
        track = static_cast<AliVTrack*>(jet->Track(i));
        z = track->Pt() / corrPt;
        fHistManager.FillTH3(histname.Data(), corrPt, z, fCent);
      }
      
      // Cent vs. Nconst vs. pT
      histname = TString::Format("%s/JetHistograms/hCentVsNConstVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH3(histname.Data(), corrPt, jet->GetNumberOfConstituents(), fCent);
      
    } //jet loop
    
    //---------------------------------------------------------------------------
    // Do study of background (if requested)
    if (fComputeBackground) ComputeBackground(jets);
  }
}

/**
 * Fill dijet candidate THnSparse.
 */
void AliAnalysisTaskEmcalDijetImbalance::FillDijetCandHistograms(TString histname)
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
    else if (title=="#it{p}_{T,trig jet} (GeV/#it{c})")
      contents[n] = fDijet.trigJetPt;
    else if (title=="#it{p}_{T,ass jet} (GeV/#it{c})")
      contents[n] = fDijet.assJetPt;
    else if (title=="#phi_{trig jet}")
      contents[n] = fDijet.trigJetPhi;
    else if (title=="#phi_{ass jet}")
      contents[n] = fDijet.assJetPhi;
    else if (title=="#eta_{trig jet}")
      contents[n] = fDijet.trigJetEta;
    else if (title=="#eta_{ass jet}")
      contents[n] = fDijet.assJetEta;
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
    else if (title=="N_{tracks, trig jet}")
      contents[n] = fDijet.trigJet->GetNumberOfTracks();
    else if (title=="N_{tracks, ass jet}")
      contents[n] = fDijet.assJet->GetNumberOfTracks();
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
    if (title=="A_{J}")
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

