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

#include <vector>

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3.h>
#include <TList.h>
#include <THnSparse.h>
#include <TRandom3.h>
#include <TGrid.h>
#include <TFile.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliAnalysisManager.h"
#include <AliVEventHandler.h>
#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

#include "AliAnalysisTaskEmcalDijetImbalance.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalDijetImbalance);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalDijetImbalance::AliAnalysisTaskEmcalDijetImbalance() : 
  AliAnalysisTaskEmcalJet(),
  fDeltaPhiMin(0),
  fMinTrigJetPt(0),
  fMinAssJetPt(0),
  fDijetLeadingHadronPt(0),
  fMatchingJetR(0.2),
  fTrackConstituentThreshold(0),
  fClusterConstituentThreshold(0),
  fNEtaBins(40),
  fNPhiBins(200),
  fBackgroundScalingWeights(0),
  fGapJetScalingWeights(0),
  fPlotDijetCandHistograms(kFALSE),
  fPlotDijetImbalanceHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoMomentumBalance(kFALSE),
  fDoGeometricalMatching(kFALSE),
  fLoadBackgroundScalingWeights(kTRUE),
  fComputeMBDownscaling(kFALSE),
  fDoTriggerSimulation(kFALSE),
  fMaxPt(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fMBUpscaleFactor(1.),
  fMedianEMCal(0),
  fMedianDCal(0),
  fkEMCEJE(kFALSE),
  fEmbeddingQA(),
  fHistManager()
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
  fDeltaPhiMin(0),
  fMinTrigJetPt(0),
  fMinAssJetPt(0),
  fDijetLeadingHadronPt(0),
  fMatchingJetR(0.2),
  fTrackConstituentThreshold(0),
  fClusterConstituentThreshold(0),
  fNEtaBins(40),
  fNPhiBins(200),
  fBackgroundScalingWeights(0),
  fGapJetScalingWeights(0),
  fPlotDijetCandHistograms(kFALSE),
  fPlotDijetImbalanceHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoMomentumBalance(kFALSE),
  fDoGeometricalMatching(kFALSE),
  fLoadBackgroundScalingWeights(kTRUE),
  fComputeMBDownscaling(kFALSE),
  fDoTriggerSimulation(kFALSE),
  fMaxPt(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fMBUpscaleFactor(1.),
  fMedianEMCal(0),
  fMedianDCal(0),
  fkEMCEJE(kFALSE),
  fEmbeddingQA(),
  fHistManager(name)
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
  
  fNPtHistBins = 82;
  fPtHistBins = new Double_t[fNPtHistBins+1];
  GenerateFixedBinArray(6, 0, 0.3, fPtHistBins);
  GenerateFixedBinArray(7, 0.3, 1, fPtHistBins+6);
  GenerateFixedBinArray(10, 1, 3, fPtHistBins+13);
  GenerateFixedBinArray(14, 3, 10, fPtHistBins+23);
  GenerateFixedBinArray(10, 10, 20, fPtHistBins+37);
  GenerateFixedBinArray(15, 20, 50, fPtHistBins+47);
  GenerateFixedBinArray(20, 50, 150, fPtHistBins+62);
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalDijetImbalance::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  if (fComputeBackground) {
    AllocateBackgroundHistograms();
  }
  if (fPlotDijetCandHistograms) {
    AllocateDijetCandHistograms();
  }
  if (fPlotDijetImbalanceHistograms) {
    AllocateDijetImbalanceHistograms();
  }
  if (fDoMomentumBalance) {
    AllocateMomentumBalanceHistograms();
  }
  if (fDoGeometricalMatching) {
    AllocateGeometricalMatchingHistograms();
  }
  if (fDoTriggerSimulation) {
    AllocateTriggerSimHistograms();
  }

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  // Intialize AliEventCuts
  if (fUseAliEventCuts) {
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
  }
  
  // Load eta-phi background scale factors from histogram on AliEn
  if (fLoadBackgroundScalingWeights) {
    LoadBackgroundScalingHistogram();
  }
  
  // Initialize embedding QA
  const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if (embeddingHelper) {
    bool res = fEmbeddingQA.Initialize();
    if (res) {
      fEmbeddingQA.AddQAPlotsToList(fOutput);
    }
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates background subtraction histograms, if enabled.
 * A set of histograms is allocated per each jet container.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateBackgroundHistograms()
{
  TString histname;
  TString title;
  
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality;Scale factor;counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 5);
    
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 400, -50, 150);
    
    histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEtaPhi", jets->GetArrayName().Data());
    title = histname + ";#eta;#phi;Centrality;Scale factor;";
    Int_t nbins[4]  = {fNEtaBins, fNPhiBins, 50, 400};
    Double_t min[4] = {-0.5,1., 0, 0};
    Double_t max[4] = {0.5,6., 100, 20};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins, min, max);
    
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEtaPhi", jets->GetArrayName().Data());
    title = histname + ";#eta;#phi;Centrality;#delta#it{p}_{T} (GeV/#it{c})";
    Int_t nbinsDpT[4]  = {fNEtaBins, fNPhiBins, 10, 400};
    Double_t minDpT[4] = {-0.5,1., 0, -50};
    Double_t maxDpT[4] = {0.5,6., 100, 150};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbinsDpT, minDpT, maxDpT);
    
  }
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
    nbins[dim] = TMath::CeilNint(fMaxPt/2);
    min[dim] = 0;
    max[dim] = fMaxPt;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#it{p}_{T,ass jet} (GeV/#it{c})";
    nbins[dim] = TMath::CeilNint(fMaxPt/2);
    min[dim] = 0;
    max[dim] = fMaxPt;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#phi_{trig jet}";
    nbins[dim] = TMath::CeilNint(fMaxPt/2);
    min[dim] = 0;
    max[dim] = TMath::TwoPi();
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#phi_{ass jet}";
    nbins[dim] = TMath::CeilNint(fMaxPt/2);
    min[dim] = 0;
    max[dim] = TMath::TwoPi();
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#eta_{trig jet}";
    nbins[dim] = TMath::CeilNint(fMaxPt/2);
    min[dim] = -1;
    max[dim] = 1;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    axisTitle[dim] = "#eta_{ass jet}";
    nbins[dim] = TMath::CeilNint(fMaxPt/2);
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
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    // Allocate dijet imbalance THnSparse, unless doing trigger simulation
    if (!fDoTriggerSimulation) {
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
    
    // Now, allocate 2d pt spectrum, upscaled according to MB downscaling factors (for trigger efficiency)
    Int_t nPtBins = TMath::CeilNint(fMaxPt/2);
    
    // (Centrality, pT1, pT2) (upscaled)
    TString histname = TString::Format("%s/DijetJetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
    TString title = histname + ";Centrality (%);#it{p}_{T,1}^{corr} (GeV/#it{c});#it{p}_{T,2}^{corr} (GeV/#it{c})";
    fHistManager.CreateTH3(histname.Data(), title.Data(), 20, 0, 100, nPtBins, 0, fMaxPt, nPtBins, 0, fMaxPt, "s");
    
    // Now, plot jet histograms for the leading jet within the dijet (for comparison to single jets, and triggered to MB)
    
    // (Centrality, pT, NEF, calo type)
    histname = TString::Format("%s/DijetJetHistograms/hNEFVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});NEF;type";
    Int_t nbins1[4]  = {20, nPtBins, 50, 3};
    Double_t min1[4] = {0, 0, 0, -0.5};
    Double_t max1[4] = {100, fMaxPt, 1., 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins1, min1, max1);
    
    // (Centrality, pT, z-leading (charged), calo type)
    histname = TString::Format("%s/DijetJetHistograms/hZLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading};type";
    Int_t nbins2[4]  = {20, nPtBins, 50, 3};
    Double_t min2[4] = {0, 0, 0, -0.5};
    Double_t max2[4] = {100, fMaxPt, 1., 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins2, min2, max2);
    
    // (Centrality, pT, z (charged), calo type)
    histname = TString::Format("%s/DijetJetHistograms/hZVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z};type";
    Int_t nbins3[4]  = {20, nPtBins, 50, 3};
    Double_t min3[4] = {0, 0, 0, -0.5};
    Double_t max3[4] = {100, fMaxPt, 1., 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins3, min3, max3);
    
    // (Centrality, pT, Nconst, calo type)
    histname = TString::Format("%s/DijetJetHistograms/hNConstVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents;type";
    Int_t nbins4[4]  = {20, nPtBins, 50, 3};
    Double_t min4[4] = {0, 0, 0, -0.5};
    Double_t max4[4] = {100, fMaxPt, fMaxPt, 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins4, min4, max4);
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
  fHistManager.CreateTH1(histname.Data(), title.Data(), 2, -0.5, 1.5);
}

/*
 * This function allocates the histograms for single jets, when the "simulated" trigger has been fired.
 * A set of histograms is allocated per each jet container.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateTriggerSimHistograms()
{
  TString histname;
  TString title;
  
  //----------------------------------------------
  // Trigger patch histograms
  
  // Median patch energy vs. centrality, for dijets
  histname = "TriggerSimHistograms/hMedPatchDijet";
  title = histname + ";Centrality (%);#it{p}_{T,trig}^{corr} (GeV/#it{c});#it{E}_{patch,med} (GeV);type";
  Int_t nbinsD[4]  = {50, 40, 100, 2};
  Double_t minD[4] = {0, 0, 0, -0.5};
  Double_t maxD[4] = {100, 200, 50, 1.5};
  fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbinsD, minD, maxD);
  
}

/**
 * Load histograms of eta-phi background scale factors from AliEn
 */
void AliAnalysisTaskEmcalDijetImbalance::LoadBackgroundScalingHistogram(const char* path, const char* name1, const char* name2)
{
  
  TString fname(path);
  if (fname.BeginsWith("alien://")) {
    TGrid::Connect("alien://");
  }
  
  TFile* file = TFile::Open(path);
  
  if (!file || file->IsZombie()) {
    ::Error("AliAnalysisTaskEmcalDijetImbalance", "Could not open background scaling histogram");
    return;
  }
  
  // Open background scale factor histogram
  TH2D* h1 = dynamic_cast<TH2D*>(file->Get(name1));
  
  if (h1) {
    ::Info("AliAnalysisTaskEmcalDijetImbalance::LoadBackgroundScalingHistogram", "Background histogram %s loaded from file %s.", name1, path);
  }
  else {
    ::Error("AliAnalysisTaskEmcalDijetImbalance::LoadBackgroundScalingHistogram", "Background histogram  %s not found in file %s.", name1, path);
    return;
  }
  
  fBackgroundScalingWeights = static_cast<TH2D*>(h1->Clone());
  
  // Open jet pT scale factor histogram
  TH2D* h2 = dynamic_cast<TH2D*>(file->Get(name2));
  
  if (h2) {
    ::Info("AliAnalysisTaskEmcalDijetImbalance::LoadBackgroundScalingHistogram", "Jet pT scaling histogram %s loaded from file %s.", name2, path);
  }
  else {
    ::Error("AliAnalysisTaskEmcalDijetImbalance::LoadBackgroundScalingHistogram", "Jet pT scaling histogram  %s not found in file %s.", name2, path);
    return;
  }
  
  fGapJetScalingWeights = static_cast<TH2D*>(h2->Clone());
  
  file->Close();
  delete file;
  
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalDijetImbalance::ExecOnce()
{
  
  // Configure base class to set fTriggerPatchInfo to array of trigger patches, each event
  // (Need to call this before base class ExecOnce)
  if (fDoTriggerSimulation) {
    this->SetCaloTriggerPatchInfoName("EmcalTriggers");
  }
  
  AliAnalysisTaskEmcalJet::ExecOnce();

  fNeedEmcalGeom = kTRUE;
  
  // Check if trigger patches are loaded
  if (fDoTriggerSimulation) {
    if (fTriggerPatchInfo) {
      TString objname(fTriggerPatchInfo->GetClass()->GetName());
      TClass cls(objname);
      if (!cls.InheritsFrom("AliEMCALTriggerPatchInfo")) {
        AliError(Form("%s: Objects of type %s in %s are not inherited from AliEMCALTriggerPatchInfo!",
                      GetName(), cls.GetName(), "EmcalTriggers"));
        fTriggerPatchInfo = 0;
      }
    }
    if (!fTriggerPatchInfo) {
      AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), "EmcalTriggers"));
      return;
    }
  }
  
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
  
  if (fComputeMBDownscaling) {
  
    // Get instance of the downscale factor helper class
    PWG::EMCAL::AliEmcalDownscaleFactorsOCDB *downscaleOCDB = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
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
    
  }
  
}

/**
 * This function (overloading the base class) uses AliEventCuts to perform event selection.
 */
Bool_t AliAnalysisTaskEmcalDijetImbalance::IsEventSelected()
{
  if (fUseAliEventCuts) {
    if (!fEventCuts.AcceptEvent(InputEvent()))
    {
      PostData(1, fOutput);
      return kFALSE;
    }
  }
  else {
    AliAnalysisTaskEmcal::IsEventSelected();
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
        FillDijetCandHistograms(jetCont);
      }
      
    }

    //---------------------------------------------------------------------------------------------------
    // Now, study the accepted dijet selection -- specified by the trig/ass jet pT conditions
    
    // Find the dijet candidate of the event and store its info in struct fDijet
    FindDijet(jetCont, 0);

    // If we find an accepted dijet, fill the dijet imbalance histogram
    if (fDijet.isAccepted && fPlotDijetImbalanceHistograms) {
      FillDijetImbalanceHistograms(jetCont);
    }
    // If we find an accepted dijet, perform momentum-balance study (if requested)
    if (fDijet.isAccepted && fDoMomentumBalance) {
      histname = TString::Format("%s/MomentumBalance", jetCont->GetArrayName().Data());
      DoMomentumBalance(histname);
    }
    
    //---------------------------------------------------------------------------
    // Do a simple trigger simulation (if requested)
    if (fDoTriggerSimulation) {
      DoTriggerSimulation();
    }
    
  }

  //---------------------------------------------------------------------------
  // Do the constituent threshold and geometrical matching study (if requested)
  if (fDoGeometricalMatching) {
    DoGeometricalMatching();
  }
  
  // Only fill the embedding qa plots if:
  //  - We are using the embedding helper
  //  - The class has been initialized
  //  - Both jet collections are available
  if (fEmbeddingQA.IsInitialized()) {
    fEmbeddingQA.RecordEmbeddedEventProperties();
  }

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
  AliEmcalJet* trigJet = GetLeadingJet(jetCont);
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
 * Do a simple trigger simulation, mimicking the median-subtraction method using cell amplitudes, in order to examine 
 * whether there is a bias in the patch median for dijets.
 */
void AliAnalysisTaskEmcalDijetImbalance::DoTriggerSimulation()
{
  TString histname;

  // Check if trigger patches are loaded
  if (fTriggerPatchInfo) {
    TString objname(fTriggerPatchInfo->GetClass()->GetName());
    TClass cls(objname);
    if (!cls.InheritsFrom("AliEMCALTriggerPatchInfo")) {
      AliError(Form("%s: Objects of type %s in %s are not inherited from AliEMCALTriggerPatchInfo!",
                    GetName(), cls.GetName(), "EmcalTriggers"));
      fTriggerPatchInfo = 0;
    }
  }
  if (!fTriggerPatchInfo) {
    AliError(Form("%s: Unable to get trigger patch container with name %s. Aborting", GetName(), "EmcalTriggers"));
    return;
  }

  // Compute patches in EMCal, DCal (I want offline simple trigger patch, i.e. patch calculated using FEE energy)
  std::vector<Double_t> vecEMCal;
  std::vector<Double_t> vecDCal;
  for(auto p : *fTriggerPatchInfo){
    AliEMCALTriggerPatchInfo *recpatch = static_cast<AliEMCALTriggerPatchInfo *>(p);
    if (recpatch) {
      
      if(!recpatch->IsJetHighSimple()) continue;

      if (recpatch->IsEMCal()) {
        vecEMCal.push_back(recpatch->GetPatchE());
      } else {
        vecDCal.push_back(recpatch->GetPatchE());
      }
      
    }
  }
  
  // Compute the median in each calorimeter
  const Int_t nBkgPatchesEMCal = vecEMCal.size(); // 6*8;
  const Int_t nBkgPatchesDCal = vecDCal.size();   // 4*5;
  fMedianEMCal = TMath::Median(nBkgPatchesEMCal, &vecEMCal[0]); // point to array used internally by vector
  fMedianDCal = TMath::Median(nBkgPatchesDCal, &vecDCal[0]);
  
  // Median patch energy vs. pT for dijets
  if (fDijet.isAccepted) {
    histname = "TriggerSimHistograms/hMedPatchDijet";
    Double_t x[4] = {fCent, fDijet.trigJetPt, fMedianEMCal, kEMCal};
    fHistManager.FillTHnSparse(histname, x);
    Double_t y[4] = {fCent, fDijet.trigJetPt, fMedianDCal, kDCal};
    fHistManager.FillTHnSparse(histname, y);
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
      trigJet = matchingAssJet;
      assJet = matchingTrigJet;
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
void AliAnalysisTaskEmcalDijetImbalance::ComputeBackground()
{
  // Loop over tracks and clusters in order to:
  //   (1) Compute scale factor for full jets
  //   (2) Compute delta-pT for full jets, with the random cone method
  // For both the scale factor and delta-pT, we compute only one histogram each for EMCal.
  // And then we bin in eta-phi, in order to compute and perform a corretion to account for the DCal vs. PHOS vs. gap
  
  // Define the acceptance boundaries for the TPC and EMCal/DCal/PHOS
  Double_t etaTPC = 0.9;
  Double_t etaEMCal = 0.7;
  //Double_t etaMinDCal = 0.22;
  //Double_t etaMaxPHOS = 0.13;
  Double_t phiMinEMCal = fGeom->GetArm1PhiMin() * TMath::DegToRad(); // 80
  Double_t phiMaxEMCal = fGeom->GetEMCALPhiMax() * TMath::DegToRad(); // ~188
  //Double_t phiMinDCal = fGeom->GetDCALPhiMin() * TMath::DegToRad(); // 260
  //Double_t phiMaxDCal = fGeom->GetDCALPhiMax() * TMath::DegToRad(); // ~327 (1/3 SMs start at 320)
  //Double_t phiMinPHOS = 250 * TMath::DegToRad();
  //Double_t phiMaxPHOS = 320 * TMath::DegToRad();

  Double_t accTPC = 2 * etaTPC * 2 * TMath::Pi();
  Double_t accEMCal = 2 * etaEMCal * (phiMaxEMCal - phiMinEMCal);
  //Double_t accDCalRegion = 2 * etaEMCal * (phiMaxDCal - phiMinDCal);
  
  // Loop over jet containers
  AliJetContainer* jetCont = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(nextJetColl()))) {
  
    // Define fiducial acceptances, to be used to generate random cones
    TRandom3* r = new TRandom3(0);
    Double_t jetR = jetCont->GetJetRadius();
    Double_t accRC = TMath::Pi() * jetR * jetR;
    Double_t etaEMCalfid = etaEMCal - jetR;
    Double_t phiMinEMCalfid = phiMinEMCal + jetR;
    Double_t phiMaxEMCalfid = phiMaxEMCal - jetR;
    
    // Generate EMCal random cone eta-phi
    Double_t etaEMCalRC = r->Uniform(-etaEMCalfid, etaEMCalfid);
    Double_t phiEMCalRC = r->Uniform(phiMinEMCalfid, phiMaxEMCalfid);
    
    // For eta-phi correction, generate random eta, phi in each eta/phi bin, to be used as center of random cone
    Double_t etaDCalRC[fNEtaBins]; // array storing the RC eta values
    Double_t etaStep = 1./fNEtaBins;
    Double_t etaMin;
    Double_t etaMax;
    for (Int_t etaBin=0; etaBin < fNEtaBins; etaBin++) {
      etaMin = -etaEMCalfid + etaBin*etaStep;
      etaMax = etaMin + etaStep;
      etaDCalRC[etaBin] = r->Uniform(etaMin, etaMax);
    }
    
    Double_t phiDCalRC[fNPhiBins]; // array storing the RC phi values
    Double_t phiStep = 5./fNPhiBins; // phi axis is [1,6] in order to have simple binning
    Double_t phiMin;
    Double_t phiMax;
    for (Int_t phiBin=0; phiBin < fNPhiBins; phiBin++) {
      phiMin = 1 + phiBin*phiStep;
      phiMax = phiMin + phiStep;
      phiDCalRC[phiBin] = r->Uniform(phiMin, phiMax);
    }
    
    // Initialize the various sums to 0
    Double_t trackPtSumTPC = 0;
    Double_t trackPtSumEMCal = 0;
    Double_t trackPtSumEMCalRC = 0;
    Double_t clusESumEMCal = 0;
    Double_t clusESumEMCalRC = 0;
    
    // Define a 2D vector (initialized to 0) to store the sum of track pT, and another for cluster ET
    std::vector<std::vector<Double_t>> trackPtSumDCalRC(fNEtaBins, std::vector<Double_t>(fNPhiBins));
    std::vector<std::vector<Double_t>> clusESumDCalRC(fNEtaBins, std::vector<Double_t>(fNPhiBins));
    
    // Loop over tracks. Sum the track pT:
    // (1) in the entire TPC, (2) in the EMCal, (3) in the EMCal random cone,
    // (4) in a random cone at each eta-phi
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
      for (Int_t etaBin=0; etaBin < fNEtaBins; etaBin++) {
        for (Int_t phiBin=0; phiBin < fNPhiBins; phiBin++) {
          deltaR = GetDeltaR(&track, etaDCalRC[etaBin], phiDCalRC[phiBin]);
          if (deltaR < jetR) {
            trackPtSumDCalRC[etaBin][phiBin] += trackPt;
          }
        }
      }
      
    }
    
    // Loop over clusters. Sum the cluster ET:
    // (1) in the EMCal, (2) in the EMCal random cone, (3) in a random cone at each eta-phi
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
      for (Int_t etaBin=0; etaBin < fNEtaBins; etaBin++) {
        for (Int_t phiBin=0; phiBin < fNPhiBins; phiBin++) {
          deltaR = GetDeltaR(&clus, etaDCalRC[etaBin], phiDCalRC[phiBin]);
          if (deltaR < jetR) {
            clusESumDCalRC[etaBin][phiBin] += clusE;
          }
        }
      }
      
    }
    
    // Compute the scale factor for EMCal, as a function of centrality
    Double_t numerator = (trackPtSumEMCal + clusESumEMCal) / accEMCal;
    Double_t denominator = trackPtSumTPC / accTPC;
    Double_t scaleFactor = numerator / denominator;
    TString histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCal", jetCont->GetArrayName().Data());
    fHistManager.FillTH2(histname, fCent, scaleFactor);
    
    // Compute the scale factor in each eta-phi bin, as a function of centrality
    for (Int_t etaBin=0; etaBin < fNEtaBins; etaBin++) {
      for (Int_t phiBin=0; phiBin < fNPhiBins; phiBin++) {
        numerator = (trackPtSumDCalRC[etaBin][phiBin] + clusESumDCalRC[etaBin][phiBin]) / accRC;
        scaleFactor = numerator / denominator;
        histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEtaPhi", jetCont->GetArrayName().Data());
        Double_t x[4] = {etaDCalRC[etaBin], phiDCalRC[phiBin], fCent, scaleFactor};
        fHistManager.FillTHnSparse(histname, x);
      }
    }
    
    // Compute delta pT for EMCal, as a function of centrality
    Double_t rho = jetCont->GetRhoVal();
    Double_t deltaPt = trackPtSumEMCalRC + clusESumEMCalRC - rho * TMath::Pi() * jetR * jetR;
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCal", jetCont->GetArrayName().Data());
    fHistManager.FillTH2(histname, fCent, deltaPt);
    
    // Compute delta pT in each eta-phi bin, as a function of centrality
    Double_t sf;
    for (Int_t etaBin=0; etaBin < fNEtaBins; etaBin++) {
      for (Int_t phiBin=0; phiBin < fNPhiBins; phiBin++) {
        if (fBackgroundScalingWeights) {
          sf = fBackgroundScalingWeights->GetBinContent(fBackgroundScalingWeights->FindBin(etaDCalRC[etaBin], phiDCalRC[phiBin]));
          rho = sf * jetCont->GetRhoVal();
        }
        deltaPt = trackPtSumDCalRC[etaBin][phiBin] + clusESumDCalRC[etaBin][phiBin] - rho * accRC;
        histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEtaPhi", jetCont->GetArrayName().Data());
        Double_t x[4] = {etaDCalRC[etaBin], phiDCalRC[phiBin], fCent, deltaPt};
        fHistManager.FillTHnSparse(histname, x);
      }
    }
    
    delete r;
    
  }

}

/**
 * Get pT of jet -- background subtracted, unless hard-core jet
 */
Double_t AliAnalysisTaskEmcalDijetImbalance::GetJetPt(AliJetContainer* jetCont, AliEmcalJet* jet)
{
  
  // Get eta-phi dependent jet pT scale factor
  Double_t jetPt = jet->Pt();
  if (fGapJetScalingWeights) {
    Double_t sf = fGapJetScalingWeights->GetBinContent(fGapJetScalingWeights->FindBin(jet->Eta(), jet->Phi_0_2pi()));
    jetPt = jetPt * (1 + sf * jet->NEF());
  }
  
  // Compute pTcorr
  Double_t rho = jetCont->GetRhoVal();
  Double_t pT = jetPt - rho * jet->Area();
  
  // If hard-core jet, don't subtract background
  TString jetContName = jetCont->GetName();
  if (jetContName.Contains("HardCore")) pT = jet->Pt();
  
  return pT;
}

/**
 * Get leading jet
 */
AliEmcalJet* AliAnalysisTaskEmcalDijetImbalance::GetLeadingJet(AliJetContainer* jetCont)
{
  AliEmcalJet* leadingJet = 0;
  
  if (jetCont->GetRhoParameter()) {
    for(auto jetCand : jetCont->accepted()) {
      if (!jetCand) continue;
      Double_t jetCandPt = GetJetPt(jetCont, jetCand);
      if (leadingJet) {
        Double_t leadingJetPt = GetJetPt(jetCont, leadingJet);
        if ( jetCandPt < leadingJetPt ) continue;
      }
      leadingJet = jetCand;
    }
  }
  else {
    leadingJet = jetCont->GetLeadingJet();
  }
  
  return leadingJet;
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
 * Get calo acceptance type of jet
 */
Double_t AliAnalysisTaskEmcalDijetImbalance::GetJetType(AliEmcalJet* jet)
{
  UInt_t jetType = jet->GetJetAcceptanceType();
  Double_t type = -1;
  if (jetType & AliEmcalJet::kEMCAL) {
    type = kEMCal;
  }
  else if (jetType & AliEmcalJet::kDCALonly) {
    type = kDCal;
  }
  else if (jetType & AliEmcalJet::kPHOS) {
    type = kPHOS;
  }

  return type;
}


/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalDijetImbalance::FillHistograms()
{

  if (fComputeBackground) {
    ComputeBackground();
  }
  
  return kTRUE;
}

/**
 * Fill dijet candidate THnSparse.
 */
void AliAnalysisTaskEmcalDijetImbalance::FillDijetCandHistograms(AliJetContainer* jets)
{
  Double_t contents[30]={0};
  TString histname = TString::Format("%s/DijetCandObservables", jets->GetArrayName().Data());
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
void AliAnalysisTaskEmcalDijetImbalance::FillDijetImbalanceHistograms(AliJetContainer* jets)
{
  // Fill the dijet imbalance histogram (unless doing trigger simulation, in which case bypass this)
  if (!fDoTriggerSimulation) {
    Double_t contents[30]={0};
    TString histname = TString::Format("%s/DijetImbalanceObservables", jets->GetArrayName().Data());
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
  
  // (Centrality, pT1, pT2) (upscaled)
  TString histname = TString::Format("%s/DijetJetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
  fHistManager.FillTH3(histname.Data(), fCent, fDijet.trigJetPt, fDijet.assJetPt, fMBUpscaleFactor);
  
  // Get jet acceptance type
  AliEmcalJet* trigJet = fDijet.trigJet;
  Double_t type = GetJetType(trigJet);
  
  // (Centrality, pT, NEF, calo type)
  histname = TString::Format("%s/DijetJetHistograms/hNEFVsPt", jets->GetArrayName().Data());
  Double_t x[4] = {fCent, fDijet.trigJetPt, trigJet->NEF(), type};
  fHistManager.FillTHnSparse(histname, x);
  
  // (Centrality, pT, z-leading (charged), calo type)
  histname = TString::Format("%s/DijetJetHistograms/hZLeadingVsPt", jets->GetArrayName().Data());
  TLorentzVector leadPart;
  jets->GetLeadingHadronMomentum(leadPart, trigJet);
  Double_t z = GetParallelFraction(leadPart.Vect(), trigJet);
  if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
  Double_t y[4] = {fCent, fDijet.trigJetPt, z, type};
  fHistManager.FillTHnSparse(histname, y);
  
  // (Centrality, pT, z (charged), calo type)
  histname = TString::Format("%s/DijetJetHistograms/hZVsPt", jets->GetArrayName().Data());
  AliVTrack* track;
  for (Int_t i=0; i<trigJet->GetNumberOfTracks(); i++) {
    track = static_cast<AliVTrack*>(trigJet->Track(i));
    z = track->Pt() / TMath::Abs(fDijet.trigJetPt);
    Double_t y2[4] = {fCent, fDijet.trigJetPt, z, type};
    fHistManager.FillTHnSparse(histname, y2);
  }
  
  // (Centrality, pT, Nconst, calo type)
  histname = TString::Format("%s/DijetJetHistograms/hNConstVsPt", jets->GetArrayName().Data());
  Double_t a[4] = {fCent, fDijet.trigJetPt, 1.*trigJet->GetNumberOfConstituents(), type};
  fHistManager.FillTHnSparse(histname, a);

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

/**
 * DijetImbalance AddTask.
 */
AliAnalysisTaskEmcalDijetImbalance* AliAnalysisTaskEmcalDijetImbalance::AddTaskEmcalDijetImbalance(const char *ntracks,
                                                                                                   const char *nclusters,
                                                                                                   const Double_t deltaPhiMin,
                                                                                                   const Bool_t doGeomMatching,
                                                                                                   const Double_t minTrPtHardCore,
                                                                                                   const Double_t minClPtHardCore,
                                                                                                   const Double_t jetR,
                                                                                                   const Bool_t includePHOS,
                                                                                                   const Double_t minTrPt,
                                                                                                   const Double_t minClPt,
                                                                                                   const char *suffix)
{
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalDijetImbalance", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalDijetImbalance", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(ntracks);
  TString clusName(nclusters);

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

  TString name("AliAnalysisTaskEmcalDijetImbalance");
  if (!trackName.IsNull()) {
    name += "_";
    name += trackName;
  }
  if (!clusName.IsNull()) {
    name += "_";
    name += clusName;
  }
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  /////////////////////////////////////////////////////////////
  // Configure di-jet task
  AliAnalysisTaskEmcalDijetImbalance* dijetTask = new AliAnalysisTaskEmcalDijetImbalance(name);
  dijetTask->SetDeltaPhiCut(deltaPhiMin);
  if (doGeomMatching) dijetTask->SetDoGeometricalMatching(doGeomMatching, jetR, minTrPtHardCore, minClPtHardCore);

  /////////////////////////////////////////////////////////////
  // Create track and cluster containers with the standard cuts

  AliParticleContainer* partCont = 0;
  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = new AliMCParticleContainer(trackName);
    partCont = mcpartCont;
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = new AliTrackContainer(trackName);
    partCont = trackCont;
  }
  if (partCont) partCont->SetParticlePtCut(minTrPt);
  if (partCont) dijetTask->AdoptParticleContainer(partCont);

  AliClusterContainer* clusCont = 0;
  if (!clusName.IsNull()) {
    clusCont = new AliClusterContainer(clusName);
    clusCont->SetClusECut(0.);
    clusCont->SetClusPtCut(0.);
    clusCont->SetClusNonLinCorrEnergyCut(0.);
    clusCont->SetClusHadCorrEnergyCut(minClPt);
    clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    if (includePHOS) {
      clusCont->SetIncludePHOS(kTRUE);
      clusCont->SetPhosMinNcells(3);
      clusCont->SetPhosMinM02(0.2);
    }
  }
  if (clusCont) dijetTask->AdoptClusterContainer(clusCont);

  /////////////////////////////////////////////////////////////
  // Create track and cluster containers for constituent study with geometrical matching (if enabled)

  if (doGeomMatching) {
    AliParticleContainer* partContThresh = 0;
    if (trackName == "mcparticles") {
      AliMCParticleContainer* mcpartCont = new AliMCParticleContainer(trackName);
      partContThresh = mcpartCont;
    }
    else if (trackName == "tracks" || trackName == "Tracks") {
      AliTrackContainer* trackCont = new AliTrackContainer(trackName);
      partContThresh = trackCont;
    }
    if (partContThresh) {
      partContThresh->SetParticlePtCut(minTrPtHardCore);
      partContThresh->SetName("tracksThresh");
      dijetTask->AdoptParticleContainer(partContThresh);
    }
    
    AliClusterContainer* clusContThresh = 0;
    if (!clusName.IsNull()) {
      clusContThresh = new AliClusterContainer(clusName);
      clusContThresh->SetName("caloClustersThresh");
      clusContThresh->SetClusECut(0.);
      clusContThresh->SetClusPtCut(0.);
      clusContThresh->SetClusNonLinCorrEnergyCut(0.);
      clusContThresh->SetClusHadCorrEnergyCut(minClPtHardCore);
      clusContThresh->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      if (includePHOS) {
        clusContThresh->SetIncludePHOS(kTRUE);
        clusContThresh->SetPhosMinNcells(3);
        clusContThresh->SetPhosMinM02(0.2);
      }
    }
    if (clusContThresh) dijetTask->AdoptClusterContainer(clusContThresh);
      }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(dijetTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                            TList::Class(),AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (dijetTask, 0,  cinput1 );
  mgr->ConnectOutput (dijetTask, 1, coutput1 );

  return dijetTask;
}
