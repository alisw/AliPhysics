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

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliPHOSGeometry.h"
#include "AliOADBContainer.h"
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
  fHistManager(),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fUseAliEventCuts(kTRUE),
  fDeltaPhiMin(0),
  fMinTrigJetPt(0),
  fMinAssJetPt(0),
  fDijetLeadingHadronPt(0),
  fMaxPt(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fPlotJetHistograms(kFALSE),
  fPlotDijetCandHistograms(kFALSE),
  fPlotDijetImbalanceHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoMomentumBalance(kFALSE),
  fDoGeometricalMatching(kFALSE),
  fMatchingJetR(0.2),
  fTrackConstituentThreshold(0),
  fClusterConstituentThreshold(0),
  fMBUpscaleFactor(1.),
  fNEtaBins(40),
  fNPhiBins(200),
  fLoadBackgroundScalingWeights(kTRUE),
  fBackgroundScalingWeights(0),
  fGapJetScalingWeights(0),
  fComputeMBDownscaling(kFALSE),
  fDoCaloStudy(kFALSE),
  fDoTriggerSimulation(kFALSE),
  fMedianEMCal(0),
  fMedianDCal(0),
  fkEMCEJE(kFALSE),
  fPlotNeutralJets(kFALSE),
  fPlotClustersInJets(kFALSE),
  fPlotClusterTHnSparse(kTRUE),
  fPlotClusWithoutNonLinCorr(kFALSE),
  fPlotExotics(kFALSE),
  fEmbeddingQA(),
  fPHOSGeo(nullptr)
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
  fUseAliEventCuts(kTRUE),
  fDeltaPhiMin(0),
  fMinTrigJetPt(0),
  fMinAssJetPt(0),
  fDijetLeadingHadronPt(0),
  fMaxPt(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fPlotJetHistograms(kFALSE),
  fPlotDijetCandHistograms(kFALSE),
  fPlotDijetImbalanceHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoMomentumBalance(kFALSE),
  fDoGeometricalMatching(kFALSE),
  fMatchingJetR(0.2),
  fTrackConstituentThreshold(0),
  fClusterConstituentThreshold(0),
  fMBUpscaleFactor(1.),
  fNEtaBins(40),
  fNPhiBins(200),
  fLoadBackgroundScalingWeights(kTRUE),
  fBackgroundScalingWeights(0),
  fGapJetScalingWeights(0),
  fComputeMBDownscaling(kFALSE),
  fDoCaloStudy(kFALSE),
  fDoTriggerSimulation(kFALSE),
  fMedianEMCal(0),
  fMedianDCal(0),
  fkEMCEJE(kFALSE),
  fPlotNeutralJets(kFALSE),
  fPlotClustersInJets(kFALSE),
  fPlotClusterTHnSparse(kTRUE),
  fPlotClusWithoutNonLinCorr(kFALSE),
  fPlotExotics(kFALSE),
  fEmbeddingQA(),
  fPHOSGeo(nullptr)
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

  if (fPlotJetHistograms) {
    AllocateJetHistograms();
  }
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
  if (fDoCaloStudy) {
    AllocateCaloHistograms();
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
 * This function allocates the histograms for single jets.
 * A set of histograms is allocated per each jet container.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateJetHistograms()
{
  TString histname;
  TString title;
  
  Int_t nPtBins = TMath::CeilNint(fMaxPt/2);
  
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    // Jet rejection reason
    histname = TString::Format("%s/JetHistograms/hJetRejectionReason", jets->GetArrayName().Data());
    title = histname + ";Rejection reason;#it{p}_{T,jet} (GeV/#it{c});counts";
    TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 50, 0, fMaxPt);
    SetRejectionReasonLabels(hist->GetXaxis());
    
    // Rho vs. Centrality
    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/JetHistograms/hRhoVsCent", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 500);
    }
    
    // (Centrality, pT, calo type)
    histname = TString::Format("%s/JetHistograms/hCentVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});type";
    fHistManager.CreateTH3(histname.Data(), title.Data(), 20, 0, 100, nPtBins, 0, fMaxPt, 3, -0.5, 2.5);
    
    // (Centrality, eta, phi, pT, NEF, calo type)
    histname = TString::Format("%s/JetHistograms/hNEFVsPtVsEtaVsPhi", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#eta_{jet} (rad);#phi_{jet} (rad);#it{p}_{T}^{corr} (GeV/#it{c});NEF;type";
    Int_t nbins1[6]  = {20, fNEtaBins, fNPhiBins, nPtBins, 50, 3};
    Double_t min1[6] = {0, -0.5,1., 0, 0, -0.5};
    Double_t max1[6] = {100, 0.5,6., fMaxPt, 1., 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 6, nbins1, min1, max1);
    
    // (Centrality, pT upscaled, calo type)
    histname = TString::Format("%s/JetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});type";
    fHistManager.CreateTH3(histname.Data(), title.Data(), 20, 0, 100, nPtBins, 0, fMaxPt, 3, -0.5, 2.5, "s");
    
    // pT-leading vs. pT
    histname = TString::Format("%s/JetHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nPtBins, 0, fMaxPt);
    
    // A vs. pT
    histname = TString::Format("%s/JetHistograms/hAreaVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, fMaxPt/3, 0, 0.5);
    
    // (Centrality, pT, z-leading (charged), calo type)
    histname = TString::Format("%s/JetHistograms/hZLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading};type";
    Int_t nbins2[4]  = {20, nPtBins, 50, 3};
    Double_t min2[4] = {0, 0, 0, -0.5};
    Double_t max2[4] = {100, fMaxPt, 1., 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins2, min2, max2);
    
    // (Centrality, pT, z (charged), calo type)
    histname = TString::Format("%s/JetHistograms/hZVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z};type";
    Int_t nbins3[4]  = {20, nPtBins, 50, 3};
    Double_t min3[4] = {0, 0, 0, -0.5};
    Double_t max3[4] = {100, fMaxPt, 1., 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins3, min3, max3);

    // (Centrality, pT, Nconst, calo type)
    histname = TString::Format("%s/JetHistograms/hNConstVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents;type";
    Int_t nbins4[4]  = {20, nPtBins, 50, 3};
    Double_t min4[4] = {0, 0, 0, -0.5};
    Double_t max4[4] = {100, fMaxPt, fMaxPt, 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins4, min4, max4);
    
    // (Median patch energy, calo type, jet pT, centrality)
    if (fDoTriggerSimulation) {
      histname = TString::Format("%s/JetHistograms/hMedPatchJet", jets->GetArrayName().Data());
      title = histname + ";#it{E}_{patch,med};type;#it{p}_{T}^{corr} (GeV/#it{c});Centrality (%)";
      Int_t nbins5[4]  = {100, 2, nPtBins, 50};
      Double_t min5[4] = {0,-0.5, 0, 0};
      Double_t max5[4] = {50,1.5, fMaxPt, 100};
      fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins5, min5, max5);
    }
    
  }
  
  // MB downscale factor histogram
  if (fComputeMBDownscaling) {
    histname = "Trigger/hMBDownscaleFactor";
    title = histname + ";Downscale factor;counts";
    TH1* hist = fHistManager.CreateTH1(histname.Data(), title.Data(), 200, 0, 200);
  }
  
}

/*
 * This function allocates background subtraction histograms, if enabled.
 * A set of histograms is allocated per each jet container.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateBackgroundHistograms()
{
  TString histname;
  TString title;
  
  Int_t nPtBins = TMath::CeilNint(fMaxPt/2);
  
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
  TH1* hist = fHistManager.CreateTH1(histname.Data(), title.Data(), 2, -0.5, 1.5);
}

/*
 * This function allocates the histograms for the calorimeter performance study.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateCaloHistograms()
{
  TString histname;
  TString htitle;
  Int_t nPtBins = TMath::CeilNint(fMaxPt/2);
  
  Double_t* clusType = new Double_t[3+1];
  GenerateFixedBinArray(3, -0.5, 2.5, clusType);
  const Int_t nRejBins = 32;
  Double_t* rejReasonBins = new Double_t[nRejBins+1];
  GenerateFixedBinArray(nRejBins, 0, nRejBins, rejReasonBins);
  const Int_t nExBins = 200;
  Double_t* exBins = new Double_t[nExBins+1];
  GenerateFixedBinArray(nExBins, 0, 1, exBins);
  
  AliEmcalContainer* cont = 0;
  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextClusColl()))) {
    
    // rejection reason plot, to make efficiency correction
    histname = TString::Format("%s/hClusterRejectionReasonEMCal", cont->GetArrayName().Data());
    htitle = histname + ";Rejection reason;#it{E}_{clus} (GeV/)";
    TH2* hist = fHistManager.CreateTH2(histname.Data(), htitle.Data(), nRejBins, rejReasonBins, fNPtHistBins, fPtHistBins);
    SetRejectionReasonLabels(hist->GetXaxis());
    
    histname = TString::Format("%s/hClusterRejectionReasonPHOS", cont->GetArrayName().Data());
    htitle = histname + ";Rejection reason;#it{E}_{clus} (GeV/)";
    TH2* histPhos = fHistManager.CreateTH2(histname.Data(), htitle.Data(), nRejBins, rejReasonBins, fNPtHistBins, fPtHistBins);
    SetRejectionReasonLabels(histPhos->GetXaxis());
    
    // plot by SM
    const Int_t nEmcalSM = 20;
    for (Int_t sm = 0; sm < nEmcalSM; sm++) {
      histname = TString::Format("%s/BySM/hEmcalClusEnergy_SM%d", cont->GetArrayName().Data(), sm);
      htitle = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins);
    }
    
    for (Int_t sm = 1; sm < 5; sm++) {
      histname = TString::Format("%s/BySM/hPhosClusEnergy_SM%d", cont->GetArrayName().Data(), sm);
      htitle = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins);
    }
  
    // Plot cluster THnSparse (centrality, cluster type, E, E-hadcorr, has matched track, M02, Ncells)
    if (fPlotClusterTHnSparse) {
      Int_t dim = 0;
      TString title[20];
      Int_t nbins[20] = {0};
      Double_t min[30] = {0.};
      Double_t max[30] = {0.};
      Double_t *binEdges[20] = {0};
      
      if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
        title[dim] = "Centrality %";
        nbins[dim] = fNCentHistBins;
        binEdges[dim] = fCentHistBins;
        min[dim] = fCentHistBins[0];
        max[dim] = fCentHistBins[fNCentHistBins];
        dim++;
      }
      
      title[dim] = "#eta";
      nbins[dim] = 28;
      min[dim] = -0.7;
      max[dim] = 0.7;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "#phi";
      nbins[dim] = 100;
      min[dim] = 1.;
      max[dim] = 6.;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "#it{E}_{clus} (GeV)";
      nbins[dim] = fNPtHistBins;
      binEdges[dim] = fPtHistBins;
      min[dim] = fPtHistBins[0];
      max[dim] = fPtHistBins[fNPtHistBins];
      dim++;
      
      title[dim] = "#it{E}_{clus, hadcorr} or #it{E}_{core} (GeV)";
      nbins[dim] = fNPtHistBins;
      binEdges[dim] = fPtHistBins;
      min[dim] = fPtHistBins[0];
      max[dim] = fPtHistBins[fNPtHistBins];
      dim++;
      
      title[dim] = "Matched track";
      nbins[dim] = 2;
      min[dim] = -0.5;
      max[dim] = 1.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "M02";
      nbins[dim] = 50;
      min[dim] = 0;
      max[dim] = 5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "Ncells";
      nbins[dim] = 30;
      min[dim] = -0.5;
      max[dim] = 29.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "Dispersion cut";
      nbins[dim] = 2;
      min[dim] = -0.5;
      max[dim] = 1.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      TString thnname = TString::Format("%s/clusterObservables", cont->GetArrayName().Data());
      THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
      for (Int_t i = 0; i < dim; i++) {
        hn->GetAxis(i)->SetTitle(title[i]);
        hn->SetBinEdges(i, binEdges[i]);
      }
    }
    
    if (fPlotExotics) {
      histname = TString::Format("%s/hFcrossEMCal", cont->GetArrayName().Data());
      htitle = histname + ";Centrality (%);Fcross;#it{E}_{clus} (GeV/)";
      TH3* hist = fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, nExBins, exBins, fNPtHistBins, fPtHistBins);
    }
  }

  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    // plot neutral jets
    if (fPlotNeutralJets) {
      TString axisTitle[30]= {""};
      Int_t nbinsJet[30]  = {0};
      Double_t minJet[30] = {0.};
      Double_t maxJet[30] = {0.};
      Double_t *binEdgesJet[20] = {0};
      Int_t dimJet = 0;
      
      if (fForceBeamType != kpp) {
        axisTitle[dimJet] = "Centrality (%)";
        nbinsJet[dimJet] = fNCentHistBins;
        binEdgesJet[dimJet] = fCentHistBins;
        minJet[dimJet] = fCentHistBins[0];
        maxJet[dimJet] = fCentHistBins[fNCentHistBins];
        dimJet++;
      }
      
      axisTitle[dimJet] = "#eta_{jet}";
      nbinsJet[dimJet] = 28;
      minJet[dimJet] = -0.7;
      maxJet[dimJet] = 0.7;
      binEdgesJet[dimJet] = GenerateFixedBinArray(nbinsJet[dimJet], minJet[dimJet], maxJet[dimJet]);
      dimJet++;
      
      axisTitle[dimJet] = "#phi_{jet} (rad)";
      nbinsJet[dimJet] = 100;
      minJet[dimJet] = 1.;
      maxJet[dimJet] = 6.;
      binEdgesJet[dimJet] = GenerateFixedBinArray(nbinsJet[dimJet], minJet[dimJet], maxJet[dimJet]);
      dimJet++;
      
      axisTitle[dimJet] = "#it{E}_{T} (GeV)";
      nbinsJet[dimJet] = fNPtHistBins;
      binEdgesJet[dimJet] = fPtHistBins;
      minJet[dimJet] = fPtHistBins[0];
      maxJet[dimJet] = fPtHistBins[fNPtHistBins];
      dimJet++;
      
      axisTitle[dimJet] = "#rho (GeV/#it{c})";
      nbinsJet[dimJet] = 100;
      minJet[dimJet] = 0.;
      maxJet[dimJet] = 1000.;
      binEdgesJet[dimJet] = GenerateFixedBinArray(nbinsJet[dimJet], minJet[dimJet], maxJet[dimJet]);
      dimJet++;
      
      axisTitle[dimJet] = "N_{clusters}";
      nbinsJet[dimJet] = 20;
      minJet[dimJet] = -0.5;
      maxJet[dimJet] = 19.5;
      binEdgesJet[dimJet] = GenerateFixedBinArray(nbinsJet[dimJet], minJet[dimJet], maxJet[dimJet]);
      dimJet++;
      
      TString thnname = TString::Format("%s/hClusterJetObservables", jets->GetArrayName().Data());
      THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dimJet, nbinsJet, minJet, maxJet);
      for (Int_t i = 0; i < dimJet; i++) {
        hn->GetAxis(i)->SetTitle(axisTitle[i]);
        hn->SetBinEdges(i, binEdgesJet[i]);
      }
    }
    
    // Plot cluster spectra within jets
    // (centrality, cluster energy, jet pT, jet eta, jet phi)
    if (fPlotClustersInJets) {
      Int_t dim = 0;
      TString title[20];
      Int_t nbins[20] = {0};
      Double_t min[30] = {0.};
      Double_t max[30] = {0.};
      Double_t *binEdges[20] = {0};
      
      if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
        title[dim] = "Centrality %";
        nbins[dim] = fNCentHistBins;
        binEdges[dim] = fCentHistBins;
        min[dim] = fCentHistBins[0];
        max[dim] = fCentHistBins[fNCentHistBins];
        dim++;
      }
      
      title[dim] = "#it{E}_{clus} (GeV)";
      nbins[dim] = fNPtHistBins;
      binEdges[dim] = fPtHistBins;
      min[dim] = fPtHistBins[0];
      max[dim] = fPtHistBins[fNPtHistBins];
      dim++;
      
      title[dim] = "#it{p}_{T,jet}^{corr}";
      nbins[dim] = nPtBins;
      min[dim] = 0;
      max[dim] = fMaxPt;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "#eta_{jet}";
      nbins[dim] = fNEtaBins;
      min[dim] = -0.5;
      max[dim] = 0.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "#phi_{jet}";
      nbins[dim] = fNPhiBins;
      min[dim] = 1.;
      max[dim] = 6.;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      TString thnname = TString::Format("%s/hClustersInJets", jets->GetArrayName().Data());
      THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
      for (Int_t i = 0; i < dim; i++) {
        hn->GetAxis(i)->SetTitle(title[i]);
        hn->SetBinEdges(i, binEdges[i]);
      }
      
      // (jet type, jet pT, cluster shift)
      histname = TString::Format("%s/hCaloJESshift", jets->GetArrayName().Data());
      htitle = histname + ";type;#it{p}_{T}^{corr} (GeV/#it{c});#Delta_{JES}";
      fHistManager.CreateTH3(histname.Data(), htitle.Data(), 3, -0.5, 2.5, nPtBins, 0, fMaxPt, 100, 0, 20);
    }

  }
  
  // Plot cell histograms
  
  // centrality vs. cell energy vs. cell type (for all cells)
  histname = TString::Format("Cells/hCellEnergyAll");
  htitle = histname + ";#it{E}_{cell} (GeV);Centrality (%); Cluster type";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNCentHistBins, fCentHistBins, 3, clusType);
  
  // centrality vs. cell energy vs. cell type (for cells in accepted clusters)
  histname = TString::Format("Cells/hCellEnergyAccepted");
  htitle = histname + ";#it{E}_{cell} (GeV);Centrality (%); Cluster type";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNCentHistBins, fCentHistBins, 3, clusType);
  
  // centrality vs. cell energy vs. cell type (for leading cells in accepted clusters)
  histname = TString::Format("Cells/hCellEnergyLeading");
  htitle = histname + ";#it{E}_{cell} (GeV);Centrality (%); Cluster type";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNCentHistBins, fCentHistBins, 3, clusType);
  
  // plot cell patches by SM
  const Int_t nEmcalSM = 20;
  for (Int_t sm = 0; sm < nEmcalSM; sm++) {
    histname = TString::Format("Cells/BySM/hEmcalPatchEnergy_SM%d", sm);
    htitle = histname + ";#it{E}_{cell patch} (GeV);Centrality (%)";
    fHistManager.CreateTH2(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNCentHistBins, fCentHistBins);
  }
  
  for (Int_t sm = 1; sm < 5; sm++) {
    histname = TString::Format("Cells/BySM/hPhosPatchEnergy_SM%d", sm);
    htitle = histname + ";#it{E}_{cell patch} (GeV);Centrality (%)";
    fHistManager.CreateTH2(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNCentHistBins, fCentHistBins);
  }
  
}

/*
 * This function allocates the histograms for single jets, when the "simulated" trigger has been fired.
 * A set of histograms is allocated per each jet container.
 */
void AliAnalysisTaskEmcalDijetImbalance::AllocateTriggerSimHistograms()
{
  TString histname;
  TString title;
  Int_t nPtBins = TMath::CeilNint(fMaxPt/2);
  
  //----------------------------------------------
  // Trigger patch histograms
  
  // patch eta vs. phi
  histname = "TriggerSimHistograms/hEtaVsPhi";
  title = histname + ";#eta_{patch} (rad);#phi_{patch} (rad)";
  fHistManager.CreateTH2(histname.Data(), title.Data(), 140, -0.7, 0.7, 500, 1., 6.);
  
  // N patches
  histname = "TriggerSimHistograms/hNPatches";
  title = histname + ";#it{N}_{patches};type";
  fHistManager.CreateTH2(histname.Data(), title.Data(), 200, 0, 200, 2, -0.5, 1.5);
  
  // patch E vs. centrality
  histname = "TriggerSimHistograms/hPatchE";
  title = histname + ";Centrality (%);#it{E}_{patch} (GeV)";
  fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, nPtBins, 0, fMaxPt);
  
  // patch median vs. Centrality
  histname = "TriggerSimHistograms/hPatchMedianE";
  title = histname + ";Centrality (%);#it{E}_{patch,med} (GeV);type";
  fHistManager.CreateTH3(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 50, 2, -0.5, 1.5);
  
  // Median patch energy vs. centrality, for dijets
  histname = "TriggerSimHistograms/hMedPatchDijet";
  title = histname + ";Centrality (%);#it{p}_{T,trig}^{corr} (GeV/#it{c});#it{E}_{patch,med} (GeV);type";
  Int_t nbinsD[4]  = {50, 40, 100, 2};
  Double_t minD[4] = {0, 0, 0, -0.5};
  Double_t maxD[4] = {100, 200, 50, 1.5};
  fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbinsD, minD, maxD);
  
  //----------------------------------------------
  // Jet histograms for "triggered" events
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    // Jet rejection reason
    histname = TString::Format("%s/TriggerSimHistograms/hJetRejectionReason", jets->GetArrayName().Data());
    title = histname + ";Rejection reason;#it{p}_{T,jet} (GeV/#it{c});counts";
    TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 50, 0, fMaxPt);
    SetRejectionReasonLabels(hist->GetXaxis());
    
    // Rho vs. Centrality
    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/TriggerSimHistograms/hRhoVsCent", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 500);
    }
    
    // (Centrality, pT, calo type)
    histname = TString::Format("%s/TriggerSimHistograms/hCentVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});type";
    fHistManager.CreateTH3(histname.Data(), title.Data(), 20, 0, 100, nPtBins, 0, fMaxPt, 3, -0.5, 2.5);
    
    // (Centrality, eta, phi, pT, NEF, calo type)
    histname = TString::Format("%s/TriggerSimHistograms/hNEFVsPtVsEtaVsPhi", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#eta_{jet} (rad);#phi_{jet} (rad);#it{p}_{T}^{corr} (GeV/#it{c});NEF;type";
    Int_t nbins1[6]  = {20, fNEtaBins, fNPhiBins, nPtBins, 50, 3};
    Double_t min1[6] = {0, -0.5,1., 0, 0, -0.5};
    Double_t max1[6] = {100, 0.5,6., fMaxPt, 1., 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 6, nbins1, min1, max1);
    
    // pT-leading vs. pT
    histname = TString::Format("%s/TriggerSimHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nPtBins, 0, fMaxPt);
    
    // A vs. pT
    histname = TString::Format("%s/TriggerSimHistograms/hAreaVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, fMaxPt/3, 0, 0.5);
    
    // (Centrality, pT, z-leading (charged), calo type)
    histname = TString::Format("%s/TriggerSimHistograms/hZLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading};type";
    Int_t nbins2[4]  = {20, nPtBins, 50, 3};
    Double_t min2[4] = {0, 0, 0, -0.5};
    Double_t max2[4] = {100, fMaxPt, 1., 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins2, min2, max2);
    
    // z (charged) vs. pT
    histname = TString::Format("%s/TriggerSimHistograms/hZVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z};type";
    Int_t nbins3[4]  = {20, nPtBins, 50, 3};
    Double_t min3[4] = {0, 0, 0, -0.5};
    Double_t max3[4] = {100, fMaxPt, 1., 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins3, min3, max3);

    // (Centrality, pT, Nconst, calo type)
    histname = TString::Format("%s/TriggerSimHistograms/hNConstVsPt", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents;type";
    Int_t nbins4[4]  = {20, nPtBins, 50, 3};
    Double_t min4[4] = {0, 0, 0, -0.5};
    Double_t max4[4] = {100, fMaxPt, fMaxPt, 2.5};
    fHistManager.CreateTHnSparse(histname.Data(), title.Data(), 4, nbins4, min4, max4);
    
  }
  
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

  // Load the PHOS geometry
  if (fDoCaloStudy) {
    fPHOSGeo = AliPHOSGeometry::GetInstance();
    if (fPHOSGeo) {
      AliInfo("Found instance of PHOS geometry!");
    }
    else {
      AliInfo("Creating PHOS geometry!");
      Int_t runNum = InputEvent()->GetRunNumber();
      if(runNum<209122) //Run1
        fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP");
      else
        fPHOSGeo =  AliPHOSGeometry::GetInstance("Run2");
      
      if (fPHOSGeo) {
        AliOADBContainer geomContainer("phosGeo");
        geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCGeometry.root","PHOSMCRotationMatrixes");
        TObjArray* matrixes = (TObjArray*)geomContainer.GetObject(runNum,"PHOSRotationMatrixes");
        for(Int_t mod=0; mod<6; mod++) {
          if(!matrixes->At(mod)) continue;
          fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod);
          printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo);
          ((TGeoHMatrix*)matrixes->At(mod))->Print();
        }
      }
    }
  }
  
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
  
  if (fPlotJetHistograms && fComputeMBDownscaling) {
  
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
 * Do a simple trigger simulation, mimicking the median-subtraction method using cell amplitudes.
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
      
      histname = "TriggerSimHistograms/hEtaVsPhi";
      fHistManager.FillTH2(histname.Data(), recpatch->GetEtaGeo(), recpatch->GetPhiGeo());
      
      histname = "TriggerSimHistograms/hPatchE";
      fHistManager.FillTH2(histname.Data(), fCent, recpatch->GetPatchE());

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
  
  histname = "TriggerSimHistograms/hPatchMedianE";
  fHistManager.FillTH3(histname.Data(), fCent, fMedianEMCal, kEMCal);
  fHistManager.FillTH3(histname.Data(), fCent, fMedianDCal, kDCal);
  
  histname = "TriggerSimHistograms/hNPatches";
  fHistManager.FillTH2(histname.Data(), nBkgPatchesEMCal, kEMCal);
  fHistManager.FillTH2(histname.Data(), nBkgPatchesDCal, kDCal);
  
  // Median patch energy vs. pT for dijets
  if (fDijet.isAccepted) {
    histname = "TriggerSimHistograms/hMedPatchDijet";
    Double_t x[4] = {fCent, fDijet.trigJetPt, fMedianEMCal, kEMCal};
    fHistManager.FillTHnSparse(histname, x);
    Double_t y[4] = {fCent, fDijet.trigJetPt, fMedianDCal, kDCal};
    fHistManager.FillTHnSparse(histname, y);
  }

  // Then compute background subtracted patches, by subtracting from each patch the median patch E from the opposite hemisphere
  // If a patch is above threshold, the event is "triggered"
  Bool_t fkEMCEJE = kFALSE;
  Double_t threshold = 20;
  for(auto p : *fTriggerPatchInfo){
    AliEMCALTriggerPatchInfo *recpatch = static_cast<AliEMCALTriggerPatchInfo *>(p);
    if (recpatch) {
      
      if(!recpatch->IsJetHighSimple()) continue;
      
      if (recpatch->IsEMCal()) {
        if ((recpatch->GetPatchE() - fMedianDCal) > threshold) {
          fkEMCEJE = kTRUE;
          break;
        }
      } else {
        if ((recpatch->GetPatchE() - fMedianEMCal) > threshold) {
          fkEMCEJE = kTRUE;
          break;
        }
      }
    }
  }
  
  if (fkEMCEJE) {
    FillTriggerSimHistograms();
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
    Double_t phiMinDCalRegionfid = phiMinDCal + jetR;
    Double_t phiMaxDCalRegionfid = phiMaxDCal - jetR;
    
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
  if (fPlotJetHistograms) {
    FillJetHistograms();
  }
  if (fComputeBackground) {
    ComputeBackground();
  }
  if (fDoCaloStudy) {
    FillCaloHistograms();
  }
  
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
      Float_t corrPt = GetJetPt(jets, jet);
      
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
      
      // compute jet acceptance type
      Double_t type = GetJetType(jet);

      // (Centrality, pT, calo type)
      histname = TString::Format("%s/JetHistograms/hCentVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH3(histname.Data(), fCent, corrPt, type);
      
      // (Centrality, eta, phi, pT, NEF, calo type)
      histname = TString::Format("%s/JetHistograms/hNEFVsPtVsEtaVsPhi", jets->GetArrayName().Data());
      Double_t x[6] = {fCent, jet->Eta(), jet->Phi_0_2pi(), corrPt, jet->NEF(), type};
      fHistManager.FillTHnSparse(histname, x);
      
      // (Centrality, pT upscaled, calo type)
      histname = TString::Format("%s/JetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
      fHistManager.FillTH3(histname.Data(), fCent, corrPt, type, fMBUpscaleFactor);
      
      // pT-leading vs. pT
      histname = TString::Format("%s/JetHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), corrPt, ptLeading);
      
      // (Centrality, pT, z-leading (charged), calo type)
      TLorentzVector leadPart;
      jets->GetLeadingHadronMomentum(leadPart, jet);
      Double_t z = GetParallelFraction(leadPart.Vect(), jet);
      if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
      histname = TString::Format("%s/JetHistograms/hZLeadingVsPt", jets->GetArrayName().Data());
      Double_t y[4] = {fCent, corrPt, z, type};
      fHistManager.FillTHnSparse(histname, y);
      
      // (Centrality, pT, z (charged), calo type)
      histname = TString::Format("%s/JetHistograms/hZVsPt", jets->GetArrayName().Data());
      AliVTrack* track;
      for (Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
        track = static_cast<AliVTrack*>(jet->Track(i));
        z = track->Pt() / TMath::Abs(corrPt);
        Double_t y2[4] = {fCent, corrPt, z, type};
        fHistManager.FillTHnSparse(histname, y2);
      }
      
      // (Centrality, pT, Nconst, calo type)
      histname = TString::Format("%s/JetHistograms/hNConstVsPt", jets->GetArrayName().Data());
      Double_t a[4] = {fCent, corrPt, 1.*jet->GetNumberOfConstituents(), type};
      fHistManager.FillTHnSparse(histname, a);
      
      // (Median patch energy, calo type, jet pT, centrality)
      if (fDoTriggerSimulation) {
        histname = TString::Format("%s/JetHistograms/hMedPatchJet", jets->GetArrayName().Data());
        Double_t x[4] = {fMedianEMCal, kEMCal, corrPt, fCent};
        fHistManager.FillTHnSparse(histname, x);
        Double_t y[4] = {fMedianDCal, kDCal, corrPt, fCent};
        fHistManager.FillTHnSparse(histname, y);
      }
      
    } //jet loop
    
  }
}

/**
 * This function performs a loop over the reconstructed jets, when the "simulated" trigger has been fired.
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalDijetImbalance::FillTriggerSimHistograms()
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
      histname = TString::Format("%s/TriggerSimHistograms/hRhoVsCent", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), fCent, rhoVal);
    }
    
    for (auto jet : jets->all()) {
      
      Float_t ptLeading = jets->GetLeadingHadronPt(jet);
      Float_t corrPt = GetJetPt(jets, jet);
      
      // A vs. pT (fill before area cut)
      histname = TString::Format("%s/TriggerSimHistograms/hAreaVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), corrPt, jet->Area());
      
      
      // Rejection reason
      UInt_t rejectionReason = 0;
      if (!jets->AcceptJet(jet, rejectionReason)) {
        histname = TString::Format("%s/TriggerSimHistograms/hJetRejectionReason", jets->GetArrayName().Data());
        fHistManager.FillTH2(histname.Data(), jets->GetRejectionReasonBitPosition(rejectionReason), jet->Pt());
        continue;
      }
      
      // compute jet acceptance type
      Double_t type = GetJetType(jet);
      
      // Centrality vs. pT
      histname = TString::Format("%s/TriggerSimHistograms/hCentVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH3(histname.Data(), fCent, corrPt, type);
      
      // (Centrality, eta, phi, pT, NEF, calo type)
      histname = TString::Format("%s/TriggerSimHistograms/hNEFVsPtVsEtaVsPhi", jets->GetArrayName().Data());
      Double_t x[6] = {fCent, jet->Eta(), jet->Phi_0_2pi(), corrPt, jet->NEF(), type};
      fHistManager.FillTHnSparse(histname, x);
      
      // pT-leading vs. pT
      histname = TString::Format("%s/TriggerSimHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), corrPt, ptLeading);
      
      // (Centrality, pT, z-leading (charged), calo type)
      TLorentzVector leadPart;
      jets->GetLeadingHadronMomentum(leadPart, jet);
      Double_t z = GetParallelFraction(leadPart.Vect(), jet);
      if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
      histname = TString::Format("%s/TriggerSimHistograms/hZLeadingVsPt", jets->GetArrayName().Data());
      Double_t y[4] = {fCent, corrPt, z, type};
      fHistManager.FillTHnSparse(histname, y);
      
      // (Centrality, pT, z (charged), calo type)
      histname = TString::Format("%s/TriggerSimHistograms/hZVsPt", jets->GetArrayName().Data());
      AliVTrack* track;
      for (Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
        track = static_cast<AliVTrack*>(jet->Track(i));
        z = track->Pt() / TMath::Abs(corrPt);
        Double_t y2[4] = {fCent, corrPt, z, type};
        fHistManager.FillTHnSparse(histname, y2);
      }
      
      // (Centrality, pT, Nconst, calo type)
      histname = TString::Format("%s/TriggerSimHistograms/hNConstVsPt", jets->GetArrayName().Data());
      Double_t a[4] = {fCent, corrPt, 1.*jet->GetNumberOfConstituents(), type};
      fHistManager.FillTHnSparse(histname, a);
      
    } //jet loop
    
  }
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

/*
 * This function fills the histograms for the calorimeter performance study.
 */
void AliAnalysisTaskEmcalDijetImbalance::FillCaloHistograms()
{
  // Define some vars
  TString histname;
  Double_t Enonlin;
  Double_t Ehadcorr;
  Int_t absId;
  Double_t ecell;
  Double_t leadEcell;
  
  // Get cells from event
  fCaloCells = InputEvent()->GetEMCALCells();
  AliVCaloCells* phosCaloCells = InputEvent()->GetPHOSCells();
    
  // Loop through clusters and plot cluster THnSparse (centrality, cluster type, E, E-hadcorr, has matched track, M02, Ncells)
  AliClusterContainer* clusters = 0;
  TIter nextClusColl(&fClusterCollArray);
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    AliClusterIterableMomentumContainer itcont = clusters->all_momentum();
    for (AliClusterIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
    
      // Determine cluster type (EMCal/DCal/Phos)
      ClusterType clusType = kNA;
      if (it->second->IsEMCAL()) {
        Double_t phi = it->first.Phi_0_2pi();
        Int_t isDcal = Int_t(phi > fgkEMCalDCalPhiDivide);
        if (isDcal == 0) {
          clusType = kEMCal;
        } else if (isDcal == 1) {
          clusType = kDCal;
        }
      } else if (it->second->GetType() == AliVCluster::kPHOSNeutral){
        clusType = kPHOS;
      }
      
      // rejection reason plots, to make efficiency correction
      if (it->second->IsEMCAL()) {
        histname = TString::Format("%s/hClusterRejectionReasonEMCal", clusters->GetArrayName().Data());
        UInt_t rejectionReason = 0;
        if (!clusters->AcceptCluster(it.current_index(), rejectionReason)) {
          fHistManager.FillTH2(histname, clusters->GetRejectionReasonBitPosition(rejectionReason), it->first.E());
          continue;
        }
      } else if (it->second->GetType() == AliVCluster::kPHOSNeutral){
        histname = TString::Format("%s/hClusterRejectionReasonPHOS", clusters->GetArrayName().Data());
        UInt_t rejectionReason = 0;
        if (!clusters->AcceptCluster(it.current_index(), rejectionReason)) {
          fHistManager.FillTH2(histname, clusters->GetRejectionReasonBitPosition(rejectionReason), it->first.E());
          continue;
        }
      } else {
        continue;
      }
      
      // Fill cluster spectra by SM, and fill cell histograms
      Enonlin = 0;
      Ehadcorr = 0;
      if (it->second->IsEMCAL()) {
        
        Ehadcorr = it->second->GetHadCorrEnergy();
        Enonlin = it->second->GetNonLinCorrEnergy();
        if (fPlotClusWithoutNonLinCorr) {
          Enonlin = it->second->E();
        }
        
        if (fPlotExotics) {
          histname = TString::Format("%s/hFcrossEMCal", clusters->GetArrayName().Data());
          Double_t Fcross = GetFcross(it->second, fCaloCells);
          fHistManager.FillTH3(histname, fCent, Fcross, it->second->E());
        }
        
        Int_t sm = fGeom->GetSuperModuleNumber(it->second->GetCellAbsId(0));
        if (sm >=0 && sm < 20) {
          histname = TString::Format("%s/BySM/hEmcalClusEnergy_SM%d", clusters->GetArrayName().Data(), sm);
          fHistManager.FillTH1(histname, it->second->E());
        }
        else {
          AliError(Form("Supermodule %d does not exist!", sm));
        }
        
        // Get cells from each accepted cluster, and plot centrality vs. cell energy vs. cell type
        histname = TString::Format("Cells/hCellEnergyAccepted");
        leadEcell = 0;
        for (Int_t iCell = 0; iCell < it->second->GetNCells(); iCell++){
          absId = it->second->GetCellAbsId(iCell);
          ecell = fCaloCells->GetCellAmplitude(absId);
          fHistManager.FillTH3(histname, ecell, fCent, kEMCal); // Note: I don't distinguish EMCal from DCal cells
          if (ecell > leadEcell) {
            leadEcell = ecell;
          }
        }
        // Plot also the leading cell
        histname = TString::Format("Cells/hCellEnergyLeading");
        fHistManager.FillTH3(histname, leadEcell, fCent, kEMCal);
        
      } else if (it->second->GetType() == AliVCluster::kPHOSNeutral){
        
        Ehadcorr = it->second->GetCoreEnergy();
        Enonlin = it->second->E();
        
        Int_t relid[4];
        if (fPHOSGeo) {
          fPHOSGeo->AbsToRelNumbering(it->second->GetCellAbsId(0), relid);
          Int_t sm = relid[0];
          if (sm >=1 && sm < 5) {
            histname = TString::Format("%s/BySM/hPhosClusEnergy_SM%d", clusters->GetArrayName().Data(), sm);
            fHistManager.FillTH1(histname, it->second->E());
          }
          else {
            AliError(Form("Supermodule %d does not exist!", sm));
          }
        }
        
        // Get cells from each accepted cluster, and plot centrality vs. cell energy vs. cell type
        histname = TString::Format("Cells/hCellEnergyAccepted");
        leadEcell = 0;
        for (Int_t iCell = 0; iCell < it->second->GetNCells(); iCell++){
          absId = it->second->GetCellAbsId(iCell);
          ecell = phosCaloCells->GetCellAmplitude(absId);
          fHistManager.FillTH3(histname, ecell, fCent, kPHOS);
          if (ecell > leadEcell) {
            leadEcell = ecell;
          }
        }
        // Plot also the leading cell
        histname = TString::Format("Cells/hCellEnergyLeading");
        fHistManager.FillTH3(histname, leadEcell, fCent, kPHOS);
      }
      
      // Check if the cluster has a matched track
      Int_t hasMatchedTrack = -1;
      Int_t nMatchedTracks = it->second->GetNTracksMatched();
      if (nMatchedTracks == 0) {
        hasMatchedTrack = 0;
      } else if (nMatchedTracks > 0) {
        hasMatchedTrack = 1;
      }
      
      // Check if the cluster passes the dispersion cut for photon-like cluster (meaningful only for PHOS)
      Int_t passedDispersionCut = 0;
      if (it->second->Chi2() < 2.5*2.5) {
        passedDispersionCut = 1;
      }
      
      if (fPlotClusterTHnSparse) {
        Double_t contents[30]={0};
        histname = TString::Format("%s/clusterObservables", clusters->GetArrayName().Data());
        THnSparse* histClusterObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));
        if (!histClusterObservables) return;
        for (Int_t i = 0; i < histClusterObservables->GetNdimensions(); i++) {
          TString title(histClusterObservables->GetAxis(i)->GetTitle());
          if (title=="Centrality %")
            contents[i] = fCent;
          else if (title=="#eta")
            contents[i] = it->first.Eta();
          else if (title=="#phi")
            contents[i] = it->first.Phi_0_2pi();
          else if (title=="#it{E}_{clus} (GeV)")
            contents[i] = Enonlin;
          else if (title=="#it{E}_{clus, hadcorr} or #it{E}_{core} (GeV)")
            contents[i] = Ehadcorr;
          else if (title=="Matched track")
            contents[i] = hasMatchedTrack;
          else if (title=="M02")
            contents[i] = it->second->GetM02();
          else if (title=="Ncells")
            contents[i] = it->second->GetNCells();
          else if (title=="Dispersion cut")
            contents[i] = passedDispersionCut;
          else
            AliWarning(Form("Unable to fill dimension %s!",title.Data()));
        }
        histClusterObservables->Fill(contents);
      }
      
    }

  }
  
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    for (auto jet : jets->accepted()) {

      // plot neutral jets THnSparse (centrality, eta, phi, E, Nclusters)
      if (fPlotNeutralJets) {
        Double_t contents[30]={0};
        histname = TString::Format("%s/hClusterJetObservables", jets->GetArrayName().Data());
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
          else if (title=="#it{E}_{T} (GeV)")
            contents[i] = jet->Pt();
          else if (title=="#rho (GeV/#it{c})")
            contents[i] = jet->Pt() / jet->Area();
          else if (title=="N_{clusters}")
            contents[i] = jet->GetNumberOfClusters();
          else
            AliWarning(Form("Unable to fill dimension %s!",title.Data()));
        }
        histJetObservables->Fill(contents);
      }
      
      // Fill cluster spectra of clusters within jets
      //(centrality, cluster energy, jet pT, jet eta, jet phi)
      if (fPlotClustersInJets) {
        histname = TString::Format("%s/hClustersInJets", jets->GetArrayName().Data());
        Int_t nClusters = jet->GetNumberOfClusters();
        AliVCluster* clus;
        for (Int_t iClus = 0; iClus < nClusters; iClus++) {
          clus = jet->Cluster(iClus);
          Double_t x[5] = {fCent, clus->E(), GetJetPt(jets, jet), jet->Eta(), jet->Phi_0_2pi()};
          fHistManager.FillTHnSparse(histname, x);
        }
        
        // Loop through clusters, and plot estimated shift in JES due to cluster bump
        // Only do for 0-10% centrality, and for EMCal/DCal
        Double_t eclus;
        Double_t shift;
        Double_t shiftSum = 0;
        if (fCent < 10.) {
          if (GetJetType(jet) > -0.5 && GetJetType(jet) < 1.5) {
            for (Int_t iClus = 0; iClus < nClusters; iClus++) {
              clus = jet->Cluster(iClus);
              eclus = clus->E();
              if (eclus > 0.5) {
                shift = 0.79 * TMath::Exp(-0.5 * ((eclus - 3.81) / 1.50)*((eclus - 3.81) / 1.50) );
                shiftSum += shift;
              }
            }
            histname = TString::Format("%s/hCaloJESshift", jets->GetArrayName().Data());
            fHistManager.FillTH3(histname, GetJetType(jet), GetJetPt(jets, jet), shiftSum);
          }
        }
      }
      
    }
    
  }
  
  // Loop through all cells and fill histos
  Int_t sm;
  Int_t relid[4];
  Double_t patchSumEMCal[20] = {0.};
  Double_t patchSumPHOS[4] = {0.};
  for (Int_t i=0; i<fCaloCells->GetNumberOfCells(); i++) {
    
    absId = fCaloCells->GetCellNumber(i);
    ecell = fCaloCells->GetCellAmplitude(absId);
    
    // Fill cell histo
    histname = TString::Format("Cells/hCellEnergyAll");
    fHistManager.FillTH3(histname, ecell, fCent, kEMCal); // Note: I don't distinguish EMCal from DCal cells
    
    // Fill cell patch histo, per SM
    sm = fGeom->GetSuperModuleNumber(absId);
    if (sm >=0 && sm < 20) {
      patchSumEMCal[sm] += ecell;
    }
    
  }
  
  for (Int_t i=0; i<phosCaloCells->GetNumberOfCells(); i++) {
    
    absId = phosCaloCells->GetCellNumber(i);
    ecell = phosCaloCells->GetCellAmplitude(absId);
    
    // Fill cell histo
    histname = TString::Format("Cells/hCellEnergyAll");
    fHistManager.FillTH3(histname, ecell, fCent, kPHOS);
    
    // Fill cell patch histo, per SM
    fPHOSGeo->AbsToRelNumbering(absId, relid);
    sm = relid[0];
    if (sm >=1 && sm < 5) {
      patchSumPHOS[sm-1] += ecell;
    }
    
  }
  
  for (Int_t sm = 0; sm < 20; sm++) {
    histname = TString::Format("Cells/BySM/hEmcalPatchEnergy_SM%d", sm);
    fHistManager.FillTH2(histname, patchSumEMCal[sm], fCent);
  }
  
  for (Int_t sm = 1; sm < 5; sm++) {
    histname = TString::Format("Cells/BySM/hPhosPatchEnergy_SM%d", sm);
    fHistManager.FillTH2(histname, patchSumPHOS[sm-1], fCent);
  }

}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalDijetImbalance::GetFcross(AliVCluster *cluster, AliVCaloCells *cells)
{
  Int_t    AbsIdseed  = -1;
  Double_t Eseed      = 0;
  for (Int_t i = 0; i < cluster->GetNCells(); i++) {
    if (cells->GetCellAmplitude(cluster->GetCellAbsId(i)) > Eseed) {
      Eseed     = cells->GetCellAmplitude(cluster->GetCellAbsId(i));
      AbsIdseed = cluster->GetCellAbsId(i);
    }
  }
  
  if (Eseed < 1e-9) {
    return 100;
  }
  
  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
  fGeom->GetCellIndex(AbsIdseed,imod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,iphi,ieta);
  
  //Get close cells index and energy, not in corners
  
  Int_t absID1 = -1;
  Int_t absID2 = -1;
  
  if (iphi < AliEMCALGeoParams::fgkEMCALRows-1) {
    absID1 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  }
  if (iphi > 0) {
    absID2 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);
  }
    
  // In case of cell in eta = 0 border, depending on SM shift the cross cell index
  
  Int_t absID3 = -1;
  Int_t absID4 = -1;
  
  if (ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2)) {
    absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
    absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1);
  }
  else if (ieta == 0 && imod%2) {
    absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
    absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1);
  }
  else {
    if (ieta < AliEMCALGeoParams::fgkEMCALCols-1) {
      absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
    }
    if (ieta > 0) {
      absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1);
    }
  }
  
  Double_t  ecell1 = cells->GetCellAmplitude(absID1);
  Double_t  ecell2 = cells->GetCellAmplitude(absID2);
  Double_t  ecell3 = cells->GetCellAmplitude(absID3);
  Double_t  ecell4 = cells->GetCellAmplitude(absID4);
  
  Double_t Ecross = ecell1 + ecell2 + ecell3 + ecell4;
  
  Double_t Fcross = 1 - Ecross/Eseed;
  
  return Fcross;
}

