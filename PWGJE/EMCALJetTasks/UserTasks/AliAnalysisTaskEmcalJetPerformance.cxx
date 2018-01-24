/**********************************************************************************
* Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                *
* All rights reserved.                                                            *
*                                                                                 *
* Redistribution and use in source and binary forms, with or without              *
* modification, are permitted provided that the following conditions are met:     *
*   * Redistributions of source code must retain the above copyright              *
*     notice, this list of conditions and the following disclaimer.               *
*   * Redistributions in binary form must reproduce the above copyright           *
*     notice, this list of conditions and the following disclaimer in the         *
*     documentation and/or other materials provided with the distribution.        *
*   * Neither the name of the <organization> nor the                              *
*     names of its contributors may be used to endorse or promote products        *
*     derived from this software without specific prior written permission.       *
*                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          *
* DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY             *
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    *
* *********************************************************************************/

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
#include "AliOADBContainer.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliMCEvent.h"

#include "AliAnalysisTaskEmcalJetPerformance.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetPerformance);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalJetPerformance::AliAnalysisTaskEmcalJetPerformance() : 
  AliAnalysisTaskEmcalJet(),
  fPlotJetHistograms(kFALSE),
  fPlotClusterHistograms(kFALSE),
  fPlotParticleCompositionHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoTriggerSimulation(kFALSE),
  fPlotMatchedJetHistograms(kFALSE),
  fComputeMBDownscaling(kFALSE),
  fMaxPt(200),
  fNEtaBins(40),
  fNPhiBins(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNM02HistBins(0),
  fM02HistBins(0),
  fNEoverPBins(0),
  fEoverPBins(0),
  fTrackMatchingDeltaEtaMax(0.015),
  fTrackMatchingDeltaPhiMax(0.030),
  fMBUpscaleFactor(1.),
  fMedianEMCal(0.),
  fMedianDCal(0.),
  fkEMCEJE(kFALSE),
  fEmbeddingQA(),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fGeneratorLevel(0),
  fHistManager()
{
  GenerateHistoBins();
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetPerformance::AliAnalysisTaskEmcalJetPerformance(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fPlotJetHistograms(kFALSE),
  fPlotClusterHistograms(kFALSE),
  fPlotParticleCompositionHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoTriggerSimulation(kFALSE),
  fPlotMatchedJetHistograms(kFALSE),
  fComputeMBDownscaling(kFALSE),
  fMaxPt(200),
  fNEtaBins(40),
  fNPhiBins(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNM02HistBins(0),
  fM02HistBins(0),
  fNEoverPBins(0),
  fEoverPBins(0),
  fTrackMatchingDeltaEtaMax(0.015),
  fTrackMatchingDeltaPhiMax(0.030),
  fMBUpscaleFactor(1.),
  fMedianEMCal(0.),
  fMedianDCal(0.),
  fkEMCEJE(kFALSE),
  fEmbeddingQA(),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fGeneratorLevel(0),
  fHistManager(name)
{
  GenerateHistoBins();
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetPerformance::~AliAnalysisTaskEmcalJetPerformance()
{
}

/**
 * Generate histogram binning arrays
 */
void AliAnalysisTaskEmcalJetPerformance::GenerateHistoBins()
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
  
  fNM02HistBins = 81;
  fM02HistBins = new Double_t[fNM02HistBins+1];
  GenerateFixedBinArray(35, 0, 0.7, fM02HistBins);
  GenerateFixedBinArray(6, 0.7, 1., fM02HistBins+35);
  GenerateFixedBinArray(20, 1., 3., fM02HistBins+41);
  GenerateFixedBinArray(10, 3., 5., fM02HistBins+61);
  GenerateFixedBinArray(10, 5., 10., fM02HistBins+71);
  
  fNEoverPBins = 47;
  fEoverPBins = new Double_t[fNEoverPBins+1];
  GenerateFixedBinArray(30, 0, 1.5, fEoverPBins);
  GenerateFixedBinArray(10, 1.5, 3.5, fEoverPBins+30);
  GenerateFixedBinArray(7, 3.5, 10.5, fEoverPBins+40);
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetPerformance::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  
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
  
  // Get the MC particle branch, in case it exists
  fGeneratorLevel = GetMCParticleContainer("mcparticles");
  
  // Allocate histograms
  if (fPlotJetHistograms) {
    AllocateJetHistograms();
  }
  if (fPlotClusterHistograms) {
    AllocateClusterHistograms();
  }
  if (fPlotParticleCompositionHistograms) {
    AllocateParticleCompositionHistograms();
  }
  if (fComputeBackground) {
    AllocateBackgroundHistograms();
  }
  if (fDoTriggerSimulation) {
    AllocateTriggerSimHistograms();
  }
  if (fPlotMatchedJetHistograms) {
    AllocateMatchedJetHistograms();
  }
  
  // Initialize embedding QA
  const AliAnalysisTaskEmcalEmbeddingHelper * embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance();
  if (embeddingHelper) {
    bool res = fEmbeddingQA.Initialize();
    if (res) {
      fEmbeddingQA.AddQAPlotsToList(fOutput);
    }
  }
  
  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for single jets.
 * A set of histograms is allocated per each jet container.
 */
void AliAnalysisTaskEmcalJetPerformance::AllocateJetHistograms()
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
    
    // (Centrality, pT, NEF)
    Int_t nbinsx = 20; Int_t minx = 0; Int_t maxx = 100;
    Int_t nbinsy = nPtBins; Int_t miny = 0; Int_t maxy = fMaxPt;
    Int_t nbinsz = 50; Int_t minz = 0; Int_t maxz = 1.;
    
    histname = TString::Format("%s/JetHistograms/hNEFVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});NEF";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    histname = TString::Format("%s/JetHistograms/hNEFVsPtDCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});NEF";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    // (Centrality, pT upscaled, calo type)
    if (fComputeMBDownscaling) {
      histname = TString::Format("%s/JetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});type";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 20, 0, 100, nPtBins, 0, fMaxPt, 2, -0.5, 1.5, "s");
    }
    
    // pT-leading vs. pT
    histname = TString::Format("%s/JetHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nPtBins, 0, fMaxPt);
    
    // A vs. pT
    histname = TString::Format("%s/JetHistograms/hAreaVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, fMaxPt/3, 0, 0.5);
    
    // (Centrality, pT, z-leading (charged))
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    
    histname = TString::Format("%s/JetHistograms/hZLeadingVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    histname = TString::Format("%s/JetHistograms/hZLeadingVsPtDCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    // (Centrality, pT, z (charged))
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    
    histname = TString::Format("%s/JetHistograms/hZVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    histname = TString::Format("%s/JetHistograms/hZVsPtDCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    // (Centrality, pT, Nconst, calo type)
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = fMaxPt;
    
    histname = TString::Format("%s/JetHistograms/hNConstVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    histname = TString::Format("%s/JetHistograms/hNConstVsPtDCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    // (Centrality, jet pT, Enonlincorr - Ehadcorr)
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
    nbinsz = nPtBins; minz = 0; maxz = fMaxPt;
    
    histname = TString::Format("%s/JetHistograms/hDeltaEHadCorr", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#sum#it{E}_{nonlincorr} - #it{E}_{hadcorr}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
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
    fHistManager.CreateTH1(histname.Data(), title.Data(), 200, 0, 200);
  }
  
}

/*
 * This function allocates the histograms for the calorimeter performance study.
 */
void AliAnalysisTaskEmcalJetPerformance::AllocateClusterHistograms()
{
  TString histname;
  TString htitle;
  
  const Int_t nRcorrBins = 50;
  Double_t *RcorrBins = GenerateFixedBinArray(nRcorrBins, 0., 1.);
  const Int_t nCellBins = 30;
  Double_t *cellBins = GenerateFixedBinArray(nCellBins, -0.5, 29.5);
  const Int_t nMatchedTrackBins = 5;
  Double_t *matchedTrackBins = GenerateFixedBinArray(nMatchedTrackBins, -0.5, 4.5);
  const Int_t nDeltaEtaBins = 60;
  Double_t *deltaEtaBins = GenerateFixedBinArray(nDeltaEtaBins, -0.015, 0.015);
  
  //////////////////////////////////////////////
  ////// Plot M02 studies
  
  // Plot M02 distribution (centrality, Eclus nonlincorr, M02)
  histname = "ClusterHistograms/hM02";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); M02";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, fNM02HistBins, fM02HistBins);
  
  // Plot Ncell distribution for M02 > 0.4 and 0.1 < M02 < 0.4 (centrality, Eclus nonlincorr, Ncells)
  histname = "ClusterHistograms/hNcellsM02G04";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); Ncells";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, nCellBins, cellBins);
  
  histname = "ClusterHistograms/hNcellsM02L04";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); Ncells";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, nCellBins, cellBins);
  
  //////////////////////////////////////////////
  ////// Plot track matching studies
  
  // Plot matched track pT for all clusters, M02 > 0.4 clusters, and 0.1 < M02 < 0.4 clusters (centrality, Eclus nonlincorr, trackPsum)
  histname = "ClusterHistograms/hMatchedTrackPt";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); #Sigma#it{p}_{track} (GeV/c)";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins);
  
  histname = "ClusterHistograms/hMatchedTrackPtM02G04";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); #Sigma#it{p}_{track} (GeV/c)";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins);
  
  histname = "ClusterHistograms/hMatchedTrackPtM02L04";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); #Sigma#it{p}_{track} (GeV/c)";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins);
  
  // Plot number of matched tracks for all clusters, M02 > 0.4 clusters, and 0.1 < M02 < 0.4 clusters (centrality, Eclus nonlincorr, N matches)
  histname = "ClusterHistograms/hMatchedTrackN";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); N_{tracks}";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, nMatchedTrackBins, matchedTrackBins);
  
  histname = "ClusterHistograms/hMatchedTrackNM02G04";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); N_{tracks}";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, nMatchedTrackBins, matchedTrackBins);
  
  histname = "ClusterHistograms/hMatchedTrackNM02L04";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); N_{tracks}";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, nMatchedTrackBins, matchedTrackBins);
  
  // Plot M02 distribution for clusters with matched tracks (centrality, Eclus nonlincorr, M02)
  histname = "ClusterHistograms/hM02Matched";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); M02";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, fNM02HistBins, fM02HistBins);
  
  // Plot M02 distribution for clusters without matched tracks (centrality, Eclus nonlincorr, M02)
  histname = "ClusterHistograms/hM02Unmatched";
  htitle = histname + ";Centrality (%);#it{E}_{clus} (GeV); M02";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, fNM02HistBins, fM02HistBins);
  
  // Plot clus-track deltaEta of matched tracks (deltaEta, Eclus, M02)
  histname = "ClusterHistograms/hDeltaEtaCentral";
  htitle = histname + ";#eta_{track} - #eta_{clus};#it{E}_{clus} (GeV); M02";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), nDeltaEtaBins, deltaEtaBins, fNPtHistBins, fPtHistBins, fNM02HistBins, fM02HistBins);
  
  histname = "ClusterHistograms/hDeltaEtaPeripheral";
  htitle = histname + ";#eta_{track} - #eta_{clus};#it{E}_{clus} (GeV); M02";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), nDeltaEtaBins, deltaEtaBins, fNPtHistBins, fPtHistBins, fNM02HistBins, fM02HistBins);
  
  //////////////////////////////////////////////
  ////// Plot E/p studies
  
  // Plot E/p vs. M02 for 0-10% and 50-90% (Eclus nonlincorr, Eclus nonlincorr / trackPsum, M02)
  histname = "ClusterHistograms/hEoverPM02Central";
  htitle = histname + ";#it{E}_{clus} (GeV); #it{E}_{clus} / #Sigma#it{p}_{track} (GeV); M02";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNEoverPBins, fEoverPBins, fNM02HistBins, fM02HistBins);
  
  histname = "ClusterHistograms/hEoverPM02Peripheral";
  htitle = histname + ";#it{E}_{clus} (GeV); #it{E}_{clus} / #Sigma#it{p}_{track} (GeV); M02";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNEoverPBins, fEoverPBins, fNM02HistBins, fM02HistBins);
  
  //////////////////////////////////////////////
  ////// Plot hadronic correction studies
  
  // Plot Rcorr distribution (centrality, trackPSum, Rcorr = (Enonlincorr - Ehadcorr) / trackPSum)
  histname = "ClusterHistograms/hRcorrVsCent";
  htitle = histname + ";Centrality (%);#Sigma#it{p}_{track} (GeV); R_{corr} = #frac{#DeltaE_{clus}}{#Sigmap_{track}}";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, nRcorrBins, RcorrBins);
  
  // Plot Rcorr distribution for 0-10% centrality (Eclus nonlincorr, trackPSum, Rcorr)
  histname = "ClusterHistograms/hRcorr0-10";
  htitle = histname + ";#it{E}_{clus} (GeV);#Sigma#it{p}_{track} (GeV); R_{corr} = #frac{#DeltaE_{clus}}{#Sigmap_{track}}";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins, nRcorrBins, RcorrBins);
  
  // Plot Rcorr distribution for 50-90% centrality (Eclus nonlincorr, trackPSum, Rcorr)
  histname = "ClusterHistograms/hRcorr50-90";
  htitle = histname + ";#it{E}_{clus} (GeV);#Sigma#it{p}_{track} (GeV); R_{corr} = #frac{#DeltaE_{clus}}{#Sigmap_{track}}";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins, nRcorrBins, RcorrBins);
  
  // Plot also Rcorr-clus (centrality, trackPSum, Rcorr-clus = (Enonlincorr - Ehadcorr) / Enonlincorr )
  histname = "ClusterHistograms/hRcorrClusVsCent";
  htitle = histname + ";Centrality (%);#Sigma#it{p}_{track} (GeV); #frac{#DeltaE_{clus}}{E_{clus}}";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNPtHistBins, fPtHistBins, nRcorrBins, RcorrBins);
  
  // Rcorr-clus for 0-10% centrality (Eclus nonlincorr, trackPSum, Rcorr-clus = (Enonlincorr - Ehadcorr) / Enonlincorr )
  histname = "ClusterHistograms/hRcorrClus0-10";
  htitle = histname + ";#it{E}_{clus} (GeV);#Sigma#it{p}_{track} (GeV); #frac{#DeltaE_{clus}}{E_{clus}}";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins, nRcorrBins, RcorrBins);
  
  // Rcorr-clus for 50-90% centrality (Eclus nonlincorr, trackPSum, Rcorr-clus = (Enonlincorr - Ehadcorr) / Enonlincorr )
  histname = "ClusterHistograms/hRcorrClus50-90";
  htitle = histname + ";#it{E}_{clus} (GeV);#Sigma#it{p}_{track} (GeV); #frac{#DeltaE_{clus}}{E_{clus}}";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins, nRcorrBins, RcorrBins);
  
  // Plot total track multiplicity
  histname = "ClusterHistograms/hTrackMultiplicity";
  htitle = histname + ";N_{tracks};Centrality (%)";
  fHistManager.CreateTH2(histname.Data(), htitle.Data(), 1000, 0, 10000, 20, 0, 100);
}

/*
 * This function allocates the histograms for the jet composition study.
 */
void AliAnalysisTaskEmcalJetPerformance::AllocateParticleCompositionHistograms()
{
  TString histname;
  TString htitle;
  Int_t nPtBins = TMath::CeilNint(fMaxPt/2);
  
  const Int_t nRejBins = 32;
  Double_t* rejReasonBins = new Double_t[nRejBins+1];
  GenerateFixedBinArray(nRejBins, 0, nRejBins, rejReasonBins);
  const Int_t nContributorTypes = 11;
  Double_t *contributorTypeBins = GenerateFixedBinArray(nContributorTypes, -0.5, 10.5);
  const Int_t nParticleTypes = 17;
  Double_t *particleTypeBins = GenerateFixedBinArray(nParticleTypes, -0.5, 16.5);
  
  AliEmcalContainer* cont = 0;
  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextClusColl()))) {
    
    histname = "ClusterHistogramsMC/hClusterRejectionReasonMC";
    htitle = histname + ";Rejection reason;#it{E}_{clus} (GeV/)";
    TH2* histMC2 = fHistManager.CreateTH2(histname.Data(), htitle.Data(), nRejBins, rejReasonBins, fNPtHistBins, fPtHistBins);
    SetRejectionReasonLabels(histMC2->GetXaxis());
  }
  
  // M02 vs. Energy vs. Particle type
  histname = "ClusterHistogramsMC/hM02VsParticleTypeCentral";
  htitle = histname + ";M02;#it{E}_{clus} (GeV); Particle type";
  TH3* hM02VsParticleTypeCentral = fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNM02HistBins, fM02HistBins, fNPtHistBins, fPtHistBins, nParticleTypes, particleTypeBins);
  SetParticleTypeLabels(hM02VsParticleTypeCentral->GetZaxis());
  
  histname = "ClusterHistogramsMC/hM02VsParticleTypePeripheral";
  htitle = histname + ";M02;#it{E}_{clus} (GeV); Particle type";
  TH3* hM02VsParticleTypePeripheral = fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNM02HistBins, fM02HistBins, fNPtHistBins, fPtHistBins, nParticleTypes, particleTypeBins);
  SetParticleTypeLabels(hM02VsParticleTypePeripheral->GetZaxis());
  
  // Plot photon energy in photon-hadron overlap clusters (Centrality, Photon energy, M02)
  histname = "ClusterHistogramsMC/hPhotonHadronPhotonEnergy";
  htitle = histname + ";Centrality (%);M02;#it{E}_{photon} (GeV)";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNM02HistBins, fM02HistBins, fNPtHistBins, fPtHistBins);
  
  // Plot hadron energy in hadron-photon overlap clusters (Centrality, Photon energy, M02)
  histname = "ClusterHistogramsMC/hHadronPhotonHadronEnergy";
  htitle = histname + ";Centrality (%);M02;#it{E}_{hadron} (GeV)";
  fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, fNM02HistBins, fM02HistBins, fNPtHistBins, fPtHistBins);
  
  if (fPlotJetHistograms) {
  
    // M02 vs. Energy vs. Particle type vs. Jet pT, for particles inside jets
    Int_t dim = 0;
    TString title[20];
    Int_t nbins[20] = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    Double_t *binEdges[20] = {0};
    
    title[dim] = "M02";
    nbins[dim] = fNM02HistBins;
    binEdges[dim] = fM02HistBins;
    min[dim] = fM02HistBins[0];
    max[dim] = fM02HistBins[fNM02HistBins];
    dim++;
    
    title[dim] = "#it{E}_{clus} (GeV)";
    nbins[dim] = fNPtHistBins;
    binEdges[dim] = fPtHistBins;
    min[dim] = fPtHistBins[0];
    max[dim] = fPtHistBins[fNPtHistBins];
    dim++;
    
    title[dim] = "Contributor type";
    nbins[dim] = nContributorTypes;
    min[dim] = -0.5;
    max[dim] = 7.5;
    binEdges[dim] = contributorTypeBins;
    dim++;
    
    title[dim] = "#it{p}_{T,jet}^{corr}";
    nbins[dim] = nPtBins;
    min[dim] = 0;
    max[dim] = fMaxPt;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    TString thnname = "JetPerformanceMC/hM02VsContributorTypeJets";
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(title[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
    
    // Particle composition inside each jet -- jet pT vs. particle type vs. particle number vs. particle pT sum
    // (One entry per jet for each particle type)
    dim = 0;
    
    title[dim] = "#it{p}_{T,jet}^{corr}";
    nbins[dim] = nPtBins;
    min[dim] = 0;
    max[dim] = fMaxPt;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    title[dim] = "Contributor type";
    nbins[dim] = nContributorTypes;
    min[dim] = -0.5;
    max[dim] = 7.5;
    binEdges[dim] = contributorTypeBins;
    dim++;
    
    title[dim] = "N";
    nbins[dim] = 30;
    min[dim] = -0.5;
    max[dim] = 29.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    title[dim] = "#it{p}_{T,sum} (GeV)";
    nbins[dim] = fNPtHistBins;
    binEdges[dim] = fPtHistBins;
    min[dim] = fPtHistBins[0];
    max[dim] = fPtHistBins[fNPtHistBins];
    dim++;
    
    thnname = "JetPerformanceMC/hJetComposition";
    THnSparse* thn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      thn->GetAxis(i)->SetTitle(title[i]);
      thn->SetBinEdges(i, binEdges[i]);
    }
    
    // Hadronic calo energy in each jet
    
    // Jet pT vs. Summed energy of hadronic clusters without a matched track
    histname = "JetPerformance/hHadCaloEnergyUnmatched";
    htitle = histname + ";#it{p}_{T,jet} (GeV);#it{p}_{T,had} (GeV)";
    fHistManager.CreateTH2(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins);
    
    // Jet pT vs. Summed energy of hadronic clusters with a matched track (before hadronic correction)
    histname = "JetPerformance/hHadCaloEnergyMatchedNonlincorr";
    htitle = histname + ";#it{p}_{T,jet} (GeV);#it{p}_{T,had} (GeV)";
    fHistManager.CreateTH2(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins);
    
    // Jet pT vs. Summed energy of hadronic clusters with a matched track (after hadronic correction)
    histname = "JetPerformance/hHadCaloEnergyMatchedHadCorr";
    htitle = histname + ";#it{p}_{T,jet} (GeV);#it{p}_{T,had} (GeV)";
    fHistManager.CreateTH2(histname.Data(), htitle.Data(), fNPtHistBins, fPtHistBins, fNPtHistBins, fPtHistBins);
    
  }
  
}

/*
 * This function sets axis labels for particle type histograms.
 */
void AliAnalysisTaskEmcalJetPerformance::SetParticleTypeLabels(TAxis* axis)
{
  axis->SetBinLabel(1,  "SinglePhoton");
  axis->SetBinLabel(2,  "SingleElectron");
  axis->SetBinLabel(3,  "SingleChargedPion");
  axis->SetBinLabel(4,  "SingleProton");
  axis->SetBinLabel(5,  "SingleAntiProton");
  axis->SetBinLabel(6,  "SingleChargedKaon");
  axis->SetBinLabel(7,  "SingleK0L");
  axis->SetBinLabel(8,  "SingleNeutron");
  axis->SetBinLabel(9,  "SingleAntiNeutron");
  axis->SetBinLabel(10, "SingleOther");
  axis->SetBinLabel(11, "PhotonHadron");
  axis->SetBinLabel(12, "HadronPhoton");
  axis->SetBinLabel(13, "MergedPi0");
  axis->SetBinLabel(14, "PhotonPhotonOther");
  axis->SetBinLabel(15, "HadronHadron");
  axis->SetBinLabel(16, "TwoContributorsOther");
  axis->SetBinLabel(17, "MoreThanTwoContributors");
}
      
/*
 * This function allocates background subtraction histograms, if enabled.
 * A set of histograms is allocated per each jet container.
 */
void AliAnalysisTaskEmcalJetPerformance::AllocateBackgroundHistograms()
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
    
  }
}

/*
 * This function allocates the histograms for single jets, when the "simulated" trigger has been fired.
 * A set of histograms is allocated per each jet container.
 */
void AliAnalysisTaskEmcalJetPerformance::AllocateTriggerSimHistograms()
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
    
    // (Centrality, pT, NEF)
    Int_t nbinsx = 20; Int_t minx = 0; Int_t maxx = 100;
    Int_t nbinsy = nPtBins; Int_t miny = 0; Int_t maxy = fMaxPt;
    Int_t nbinsz = 50; Int_t minz = 0; Int_t maxz = 1.;
    
    histname = TString::Format("%s/TriggerSimHistograms/hNEFVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});NEF";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    histname = TString::Format("%s/TriggerSimHistograms/hNEFVsPtDCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});NEF";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    // pT-leading vs. pT
    histname = TString::Format("%s/TriggerSimHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, nPtBins, 0, fMaxPt);
    
    // A vs. pT
    histname = TString::Format("%s/TriggerSimHistograms/hAreaVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins, 0, fMaxPt, fMaxPt/3, 0, 0.5);
    
    // (Centrality, pT, z-leading (charged))
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    
    histname = TString::Format("%s/TriggerSimHistograms/hZLeadingVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    histname = TString::Format("%s/TriggerSimHistograms/hZLeadingVsPtDCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    // z (charged) vs. pT
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    
    histname = TString::Format("%s/TriggerSimHistograms/hZVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    histname = TString::Format("%s/TriggerSimHistograms/hZVsPtDCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    // (Centrality, pT, Nconst)
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = fMaxPt;
    
    histname = TString::Format("%s/TriggerSimHistograms/hNConstVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    histname = TString::Format("%s/TriggerSimHistograms/hNConstVsPtDCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
  }
}

/*
 * This function allocates histograms for matched truth-det jets in the case of embedding.
 * The jet matching information must be previously filled by another task, such as AliJetResponseMaker.
 */
void AliAnalysisTaskEmcalJetPerformance::AllocateMatchedJetHistograms()
{
  TString histname;
  TString title;
  Int_t nPtBins = TMath::CeilNint(fMaxPt/2);
  
  // Response matrix, (centrality, pT-truth, pT-det)
  Int_t nbinsx = 20; Int_t minx = 0; Int_t maxx = 100;
  Int_t nbinsy = fMaxPt; Int_t miny = 0; Int_t maxy = fMaxPt;
  Int_t nbinsz = fMaxPt; Int_t minz = 0; Int_t maxz = fMaxPt;
  
  histname = "MatchedJetHistograms/hResponseMatrixEMCal";
  title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,corr}^{det} (GeV/#it{c})";
  fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  
  histname = "MatchedJetHistograms/hResponseMatrixDCal";
  title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,corr}^{det} (GeV/#it{c})";
  fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  
  // JES shift, (centrality, pT-truth, (pT-det - pT-truth) / pT-truth)
  nbinsx = 20; minx = 0; maxx = 100;
  nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
  nbinsz = 250; minz = -5.; maxz = 5.;
  
  histname = "MatchedJetHistograms/hJESshiftEMCal";
  title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#frac{#it{p}_{T,corr}^{det} - #it{p}_{T}^{truth}}{#it{p}_{T}^{truth}}";
  fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  
  histname = "MatchedJetHistograms/hJESshiftDCal";
  title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#frac{#it{p}_{T,corr}^{det} - #it{p}_{T}^{truth}}{#it{p}_{T}^{truth}}";
  fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  
  // NEF of det-level matched jets, (centrality, pT-truth, NEF)
  nbinsx = 20; minx = 0; maxx = 100;
  nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
  nbinsz = 50; minz = 0; maxz = 1.;
  
  histname = "MatchedJetHistograms/hNEFVsPt";
  title = histname + ";Centrality (%);#it{p}_{T,corr}^{det} (GeV/#it{c});Calo energy fraction";
  fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  
  // z-leading (charged) of det-level matched jets, (centrality, pT-truth, z-leading)
  nbinsx = 20; minx = 0; maxx = 100;
  nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
  nbinsz = 50; minz = 0; maxz = 1.;
  
  histname = "MatchedJetHistograms/hZLeadingVsPt";
  title = histname + ";Centrality (%);#it{p}_{T,corr}^{det} (GeV/#it{c});#it{z}_{leading}";
  fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  
  // Matching distance, (centrality, pT-truth, R)
  nbinsx = 20; minx = 0; maxx = 100;
  nbinsy = nPtBins; miny = 0; maxy = fMaxPt;
  nbinsz = 50; minz = 0; maxz = 1.;
  
  histname = "MatchedJetHistograms/hMatchingDistance";
  title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});R";
  fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);

}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetPerformance::ExecOnce()
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
}

/**
 * This function is called automatically when the run number changes.
 */
void AliAnalysisTaskEmcalJetPerformance::RunChanged(Int_t run){
  
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
Bool_t AliAnalysisTaskEmcalJetPerformance::IsEventSelected()
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
Bool_t AliAnalysisTaskEmcalJetPerformance::Run()
{
  TString histname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    TString jetContName = jetCont->GetName();
    
    // Do a simple trigger simulation (if requested)
    if (fDoTriggerSimulation) {
      DoTriggerSimulation();
    }
    
  }
  
  // Compute the full jet background scale factor and delta-pt
  if (fComputeBackground) {
    ComputeBackground();
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
 * Do a simple trigger simulation, mimicking the median-subtraction method using cell amplitudes.
 */
void AliAnalysisTaskEmcalJetPerformance::DoTriggerSimulation()
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
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetPerformance::FillHistograms()
{
  
  if (fPlotJetHistograms) {
    FillJetHistograms();
  }
  if (fPlotClusterHistograms) {
    FillClusterHistograms();
  }
  if (fPlotParticleCompositionHistograms) {
    FillParticleCompositionHistograms();
  }
  if (fPlotMatchedJetHistograms) {
    FillMatchedJetHistograms();
  }
  
  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetPerformance::FillJetHistograms()
{
  TString histname;
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    TString jetContName = jets->GetName();
    
    Double_t rhoVal = 0;
    if (jets->GetRhoParameter()) {
      rhoVal = jets->GetRhoVal();
      histname = TString::Format("%s/JetHistograms/hRhoVsCent", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), fCent, rhoVal);
    }
    
    for (auto jet : jets->all()) {
      
      Float_t ptLeading = jets->GetLeadingHadronPt(jet);
      Float_t corrPt = GetJetPt(jet, rhoVal);
      
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
      if ( (type != kEMCal) && (type != kDCal) ) {
        continue;
      }
      
      // (Centrality, pT, NEF)
      if (type == kEMCal) {
        histname = TString::Format("%s/JetHistograms/hNEFVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/JetHistograms/hNEFVsPtDCal", jets->GetArrayName().Data());
      }
      fHistManager.FillTH3(histname, fCent, corrPt, jet->NEF());
      
      // (Centrality, pT upscaled, calo type)
      if (fComputeMBDownscaling) {
        histname = TString::Format("%s/JetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
        fHistManager.FillTH3(histname.Data(), fCent, corrPt, type, fMBUpscaleFactor);
      }
      
      // pT-leading vs. pT
      histname = TString::Format("%s/JetHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), corrPt, ptLeading);
      
      // (Centrality, pT, z-leading (charged))
      if (type == kEMCal) {
        histname = TString::Format("%s/JetHistograms/hZLeadingVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/JetHistograms/hZLeadingVsPtDCal", jets->GetArrayName().Data());
      }
      TLorentzVector leadPart;
      jets->GetLeadingHadronMomentum(leadPart, jet);
      Double_t z = GetParallelFraction(leadPart.Vect(), jet);
      if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin <1
      fHistManager.FillTH3(histname, fCent, corrPt, z);
      
      // (Centrality, pT, z (charged))
      if (type == kEMCal) {
        histname = TString::Format("%s/JetHistograms/hZVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/JetHistograms/hZVsPtDCal", jets->GetArrayName().Data());
      }
      const AliVTrack* track;
      for (Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
        track = static_cast<AliVTrack*>(jet->Track(i));
        z = track->Pt() / TMath::Abs(corrPt);
        fHistManager.FillTH3(histname, fCent, corrPt, z);
      }
 
      // (Centrality, pT, Nconst)
      if (type == kEMCal) {
        histname = TString::Format("%s/JetHistograms/hNConstVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/JetHistograms/hNConstVsPtDCal", jets->GetArrayName().Data());
      }
      fHistManager.FillTH3(histname, fCent, corrPt, 1.*jet->GetNumberOfConstituents());
      
      // (Centrality, jet pT, Enonlincorr - Ehadcorr)
      Double_t deltaEhadcorr = 0;
      const AliVCluster* clus = nullptr;
      Int_t nClusters = jet->GetNumberOfClusters();
      for (Int_t iClus = 0; iClus < nClusters; iClus++) {
        clus = jet->Cluster(iClus);
        deltaEhadcorr += (clus->GetNonLinCorrEnergy() - clus->GetHadCorrEnergy());
      }
      
      histname = TString::Format("%s/JetHistograms/hDeltaEHadCorr", jets->GetArrayName().Data());
      fHistManager.FillTH3(histname, fCent, corrPt, deltaEhadcorr);
      
      
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
      
/*
 * This function fills the histograms for the calorimeter performance study.
 */
void AliAnalysisTaskEmcalJetPerformance::FillClusterHistograms()
{
  TString histname;
  
  // Loop through clusters
  AliClusterContainer* clusters = GetClusterContainer(0);
  const AliVCluster* clus;
  for (auto it : clusters->accepted_momentum()) {
    
    clus = it.second;
    Double_t clusPhi = it.first.Phi_0_2pi();
    Double_t clusEta = it.first.Eta();
    
    // Include only EMCal/DCal clusters
    if (!clus->IsEMCAL()) {
      continue;
    }
    
    // Compute the sum of matched track momentum, and various track matching / hadronic corretion quantities
    Double_t trackPSum = 0;
    Int_t nTracksMatched = 0;
    const AliVTrack* track = nullptr;
    for (Int_t itrack=0; itrack < clus->GetNTracksMatched(); itrack++) {
      track = dynamic_cast<AliVTrack*>(clus->GetTrackMatched(itrack));
      if (!track) {
        continue;
      }
      
      Double_t trackPhi = TVector2::Phi_0_2pi(track->GetTrackPhiOnEMCal());
      Double_t trackEta = track->GetTrackEtaOnEMCal();
      Double_t deta = TMath::Abs(clusEta - trackEta);
      Double_t dphi = TMath::Abs(clusPhi - trackPhi);
      
      if (deta < fTrackMatchingDeltaEtaMax && dphi < fTrackMatchingDeltaPhiMax) {
        trackPSum += track->P();
        nTracksMatched++;
      }
    }
    
    Double_t EoverP = 0;
    if (trackPSum > 1e-3) {
      EoverP = clus->GetNonLinCorrEnergy() / trackPSum;
    }
    
    Double_t deltaE = clus->GetNonLinCorrEnergy() - clus->GetHadCorrEnergy();
    Double_t Rcorr = 0;
    if (trackPSum > 1e-3) {
      Rcorr = deltaE / trackPSum;
    }
    Double_t RcorrClus = deltaE / clus->GetNonLinCorrEnergy();
    
    //////////////////////////////////////////////
    ////// Plot M02 studies

    // Fill M02 distribution (centrality, Eclus nonlincorr, M02)
    histname = "ClusterHistograms/hM02";
    fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), clus->GetM02());
    
    // Plot Ncell distribution for M02 > 0.4 or M02 < 0.4 (centrality, Eclus nonlincorr, Ncells)
    if (clus->GetM02() > 0.4) {
      histname = "ClusterHistograms/hNcellsM02G04";
      fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), clus->GetNCells());
    }
    if (clus->GetM02() > 0.1 && clus->GetM02() < 0.4) {
      histname = "ClusterHistograms/hNcellsM02L04";
      fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), clus->GetNCells());
    }
    
    //////////////////////////////////////////////
    ////// Plot track matching studies
  
    // Plot matched track pT (centrality, Eclus nonlincorr, trackPsum)
    histname = "ClusterHistograms/hMatchedTrackPt";
    fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), trackPSum);
    
    if (clus->GetM02() > 0.4) {
      histname = "ClusterHistograms/hMatchedTrackPtM02G04";
      fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), trackPSum);
    }
    if (clus->GetM02() > 0.1 && clus->GetM02() < 0.4) {
      histname = "ClusterHistograms/hMatchedTrackPtM02L04";
      fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), trackPSum);
    }
    
    // Plot number of matched tracks (centrality, Eclus nonlincorr, N matches)
    histname = "ClusterHistograms/hMatchedTrackN";
    fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), nTracksMatched);
    
    if (clus->GetM02() > 0.4) {
      histname = "ClusterHistograms/hMatchedTrackNM02G04";
      fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), nTracksMatched);
    }
    if (clus->GetM02() > 0.1 && clus->GetM02() < 0.4) {
      histname = "ClusterHistograms/hMatchedTrackNM02L04";
      fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), nTracksMatched);
    }
    
    // Plot M02 distribution for clusters with matched tracks (centrality, Eclus nonlincorr, M02)
    histname = "ClusterHistograms/hM02Matched";
    if (nTracksMatched > 0) {
      fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), clus->GetM02());
    }
    
    // Plot M02 distribution for clusters without matched tracks (centrality, Eclus nonlincorr, M02)
    histname = "ClusterHistograms/hM02Unmatched";
    if (nTracksMatched == 0) {
      fHistManager.FillTH3(histname, fCent, clus->GetNonLinCorrEnergy(), clus->GetM02());
    }
    
    // Plot clus-track deltaEta if there is one matched track (deltaEta, Eclus, M02)
    if (nTracksMatched == 1) {
      
      const AliVTrack* track = dynamic_cast<AliVTrack*>(clus->GetTrackMatched(0));
      if (track) {
        Double_t trackEta = track->GetTrackEtaOnEMCal();
        Double_t deta = trackEta - clusEta;
        
        if (fCent > 0 && fCent < 10) {
          histname = "ClusterHistograms/hDeltaEtaCentral";
          fHistManager.FillTH3(histname.Data(), deta, clus->GetNonLinCorrEnergy(), clus->GetM02());
        }
        
        if (fCent > 50 && fCent < 90) {
          histname = "ClusterHistograms/hDeltaEtaPeripheral";
          fHistManager.FillTH3(histname.Data(), deta, clus->GetNonLinCorrEnergy(), clus->GetM02());
        }
      }
    }
    
    //////////////////////////////////////////////
    ////// Plot E/p studies

    // Plot E/p vs. M02 for 0-10% and 50-90% (Eclus nonlincorr, Eclus nonlincorr / trackPsum, M02)
    if (fCent > 0 && fCent < 10) {
      histname = "ClusterHistograms/hEoverPM02Central";
      fHistManager.FillTH3(histname, clus->GetNonLinCorrEnergy(), EoverP, clus->GetM02());
    }
    if (fCent > 50 && fCent < 90) {
      histname = "ClusterHistograms/hEoverPM02Peripheral";
      fHistManager.FillTH3(histname, clus->GetNonLinCorrEnergy(), EoverP, clus->GetM02());
    }
    
    //////////////////////////////////////////////
    ////// Plot hadronic correction studies
    
    // Fill Rcorr distribution (centrality, trackPSum, Rcorr = (Enonlincorr - Ehadcorr) / trackPSum)
    histname = "ClusterHistograms/hRcorrVsCent";
    fHistManager.FillTH3(histname, fCent, trackPSum, Rcorr);
    
    // Fill Rcorr distribution for 0-10% centrality (Eclus nonlincorr, trackPSum, Rcorr)
    if (fCent > 0 && fCent < 10) {
      histname = "ClusterHistograms/hRcorr0-10";
      fHistManager.FillTH3(histname, clus->GetNonLinCorrEnergy(), trackPSum, Rcorr);
    }
    
    // Fill Rcorr distribution for 50-90% centrality (Eclus nonlincorr, trackPSum, Rcorr)
    if (fCent > 50 && fCent < 90) {
      histname = "ClusterHistograms/hRcorr50-90";
      fHistManager.FillTH3(histname, clus->GetNonLinCorrEnergy(), trackPSum, Rcorr);
    }
    
    // Fill also Rcorr-clus (centrality, trackPSum, Rcorr-clus = (Enonlincorr - Ehadcorr) / Enonlincorr )
    histname = "ClusterHistograms/hRcorrClusVsCent";
    fHistManager.FillTH3(histname, fCent, trackPSum, RcorrClus);
    
    // Fill Rcorr-clus for 0-10% centrality (Eclus nonlincorr, trackPSum, Rcorr-clus = (Enonlincorr - Ehadcorr) / Enonlincorr )
    if (fCent > 0 && fCent < 10) {
      histname = "ClusterHistograms/hRcorrClus0-10";
      fHistManager.FillTH3(histname, clus->GetNonLinCorrEnergy(), trackPSum, RcorrClus);
    }
    
    // Fill Rcorr-clus for 50-90% centrality (Eclus nonlincorr, trackPSum, Rcorr-clus = (Enonlincorr - Ehadcorr) / Enonlincorr )
    if (fCent > 50 && fCent < 90) {
      histname = "ClusterHistograms/hRcorrClus50-90";
      fHistManager.FillTH3(histname, clus->GetNonLinCorrEnergy(), trackPSum, RcorrClus);
    }
    
  }
  
  // Fill total track multiplicity
  histname = "ClusterHistograms/hTrackMultiplicity";
  AliTrackContainer* trackCont = dynamic_cast<AliTrackContainer*>(GetParticleContainer("tracks"));
  Int_t nTracks = trackCont->GetNAcceptedTracks();
  fHistManager.FillTH2(histname.Data(), nTracks, fCent);

}

/*
 * This function fills particle composition histograms for the calorimeter performance study in MC.
 */
void AliAnalysisTaskEmcalJetPerformance::FillParticleCompositionHistograms()
{
  // If MC, get the MC event
  const AliMCEvent* mcevent = nullptr;
  if (fGeneratorLevel) {
    mcevent = MCEvent();
  }
  else {
    return;
  }
  
  // Loop through clusters, and plot M02 for each particle type
  FillParticleCompositionClusterHistograms(mcevent);
  
  // Loop through jets, to fill various histograms
  if (fPlotJetHistograms) {
    FillParticleCompositionJetHistograms(mcevent);
  }
  
}

/*
 * Loop through clusters, and plot M02 for each particle type
 */
void AliAnalysisTaskEmcalJetPerformance::FillParticleCompositionClusterHistograms(const AliMCEvent* mcevent)
{
  TString histname;
  AliClusterContainer* clusters = GetClusterContainer(0);
  const AliVCluster* clus;
  std::vector<ContributorType> vecContributorTypes;
  std::vector<Int_t> vecContributorLabels;
  for (auto it : clusters->accepted_momentum()) {
    
    clus = it.second;
    
    // Include only EMCal/DCal clusters (reject PHOS clusters)
    if (!clus->IsEMCAL()) {
      continue;
    }
    
    // Loop through the cluster's contributors in order to classify its type
    ParticleType particleType = kNotDefined;
    ContributorType contributorType = kUndefined;
    const Int_t nLabels = clus->GetNLabels();
    
    // Create a vector to store the contributor types for PhysicalPrimary particles
    vecContributorTypes.clear();
    vecContributorLabels.clear();
    for (Int_t iLabel=0; iLabel<nLabels; iLabel++) {
      
      Int_t label = clus->GetLabels()[iLabel];
      if (TMath::Abs(label) > 0) { // if the particle has a truth-level match, the label is nonzero
        contributorType = GetContributorType(clus, mcevent, label);
        if (contributorType != kUndefined) {
          vecContributorTypes.push_back(contributorType);
          vecContributorLabels.push_back(label);
        }
      }
    }
    
    Int_t nLabelsPhysPrim = vecContributorTypes.size();
    
    if (nLabelsPhysPrim == 1) {
      
      contributorType = vecContributorTypes[0];
      
      if (contributorType == kPhoton) {
        particleType = kSinglePhoton;
      }
      else if (contributorType == kElectron) {
        particleType = kSingleElectron;
      }
      else if (contributorType == kChargedPion) {
        particleType = kSingleChargedPion;
      }
      else if (contributorType == kProton) {
        particleType = kSingleProton;
      }
      else if (contributorType == kAntiProton) {
        particleType = kSingleAntiProton;
      }
      else if (contributorType == kChargedKaon) {
        particleType = kSingleChargedKaon;
      }
      else if (contributorType == kK0L) {
        particleType = kSingleK0L;
      }
      else if (contributorType == kNeutron) {
        particleType = kSingleNeutron;
      }
      else if (contributorType == kAntiNeutron) {
        particleType = kSingleAntiNeutron;
      }
      else {
        particleType = kSingleOther;
      }
      
    }
    else if (nLabelsPhysPrim == 2) {
      
      // Get the contributor particle types
      ContributorType contributorType1 = vecContributorTypes[0];
      ContributorType contributorType2 = vecContributorTypes[1];
      
      // Get the fraction of cluster energy from each contributor
      //Double_t frac0 = clus->GetClusterMCEdepFraction(0);
      Double_t frac1 = clus->GetClusterMCEdepFraction(1);
      
      // Check whether the leading/subleading contributors are photons/hadrons
      Bool_t isHadron1 = IsHadron(contributorType1);
      Bool_t isHadron2 = IsHadron(contributorType2);
      Bool_t isPhoton1 = contributorType1 == kPhoton;
      Bool_t isPhoton2 = contributorType2 == kPhoton;
      
      if (isHadron1 && isHadron2) {
        particleType = kHadronHadron;
      }
      else if (isPhoton1 && isHadron2) {
        particleType = kPhotonHadron;
        
        // Plot cluster energy when subleading hadron is subtracted
        Double_t photonEnergy = clus->GetNonLinCorrEnergy() * (1 - frac1);
        histname = "ClusterHistogramsMC/hPhotonHadronPhotonEnergy";
        fHistManager.FillTH3(histname.Data(), fCent, clus->GetM02(), photonEnergy);
      }
      else if (isHadron1 && isPhoton2) {
        particleType = kHadronPhoton;
        
        // Plot cluster energy when subleading hadron is subtracted
        Double_t hadronEnergy = clus->GetNonLinCorrEnergy() * (1 - frac1);
        histname = "ClusterHistogramsMC/hHadronPhotonHadronEnergy";
        fHistManager.FillTH3(histname.Data(), fCent, clus->GetM02(), hadronEnergy);
      }
      else if (isPhoton1 && isPhoton2) {
        
        // By default, assume the two photons are not a merged pi0
        particleType = kPhotonPhotonOther;

        // Using the vector of accepted contributor labels, check whether the two photons are from the same pi0
        AliAODMCParticle *part1 = fGeneratorLevel->GetMCParticleWithLabel(vecContributorLabels[0]);
        AliAODMCParticle *part2 = fGeneratorLevel->GetMCParticleWithLabel(vecContributorLabels[1]);
        if (part1 && part2) {
          Int_t iMother1 = part1->GetMother();
          Int_t iMother2 = part2->GetMother();
          AliVParticle *mother1 = mcevent->GetTrack(iMother1);
          AliVParticle *mother2 = mcevent->GetTrack(iMother2);
          
          if (mother1 && mother2) {
            if ( (mother1->PdgCode() == 111) && (mother2->PdgCode() == 111) ) {
              if (iMother1 == iMother2) {
                particleType = kMergedPi0;
              }
            }
          }
        }
      }
      else {
        particleType = kTwoContributorsOther; // this includes partially contained conversion overlaps
      }
      
    }
    else if (nLabelsPhysPrim > 2) {
      particleType = kMoreThanTwoContributors;
    }
    
    // (M02, Eclus, part type)
    if (fCent > 0 && fCent < 10) {
      histname = "ClusterHistogramsMC/hM02VsParticleTypeCentral";
      fHistManager.FillTH3(histname, clus->GetM02(), clus->GetNonLinCorrEnergy(), particleType);
    }
    if (fCent > 50 && fCent < 90) {
      histname = "ClusterHistogramsMC/hM02VsParticleTypePeripheral";
      fHistManager.FillTH3(histname, clus->GetM02(), clus->GetNonLinCorrEnergy(), particleType);
    }
    
  }
}

/*
 * Loop through jets, to fill particle composition histograms.
 */
void AliAnalysisTaskEmcalJetPerformance::FillParticleCompositionJetHistograms(const AliMCEvent* mcevent)
{
  TString histname;
  
  const AliVCluster* clus;
  AliJetContainer* jets = GetJetContainer(0); // there is only a single, det-level jet finder here
  for (const auto jet : jets->accepted()) {
    
    Double_t jetPt = GetJetPt(jet, 0);
    
    // Keep track of hadronic calo energy in each jet
    Double_t hadCaloEnergyUnmatched = 0;
    Double_t hadCaloEnergyMatchedNonlincorr = 0;
    Double_t hadCaloEnergyMatchedHadCorr = 0;
    
    // Loop through clusters in each jet, to plot several histograms
    Int_t nClusters = jet->GetNumberOfClusters();
    for (Int_t iClus = 0; iClus < nClusters; iClus++) {
      
      clus = jet->Cluster(iClus);
      
      // Get the particle type of the cluster
      ContributorType contributorType = kUndefined;
      Int_t label = TMath::Abs(clus->GetLabel());
      if (label > 0) {
        contributorType = GetContributorType(clus, mcevent, label);
      }
      
      // Plot M02 for each particle type
      histname = "JetPerformanceMC/hM02VsContributorTypeJets";
      Double_t x[4] = {clus->GetM02(), clus->GetNonLinCorrEnergy(), contributorType, jetPt};
      fHistManager.FillTHnSparse(histname, x);
      
      // If the cluster is a hadron, sum its energy to compute the jet's hadronic calo energy
      Bool_t isHadron = IsHadron(contributorType);
      if (isHadron) {
        Bool_t hasMatchedTrack = (clus->GetNTracksMatched() > 0);
        //Bool_t hasMatchedTrack = ((clus->GetNonLinCorrEnergy() - clus->GetHadCorrEnergy()) > 1e-3);
        if (hasMatchedTrack) {
          hadCaloEnergyMatchedNonlincorr += clus->GetNonLinCorrEnergy();
          hadCaloEnergyMatchedHadCorr += clus->GetHadCorrEnergy();
        }
        else {
          hadCaloEnergyUnmatched += clus->GetNonLinCorrEnergy();
        }
      }
      
    }
    
    // Fill hadronic calo energy in each jet
    
    // (Jet pT, Summed energy of hadronic clusters without a matched track)
    histname = "JetPerformance/hHadCaloEnergyUnmatched";
    fHistManager.FillTH2(histname, jetPt, hadCaloEnergyUnmatched);
    
    // (Jet pT vs. Summed energy of hadronic clusters with a matched track (before hadronic correction))
    histname = "JetPerformance/hHadCaloEnergyMatchedNonlincorr";
    fHistManager.FillTH2(histname, jetPt, hadCaloEnergyMatchedNonlincorr);
    
    // (Jet pT vs. Summed energy of hadronic clusters with a matched track (after hadronic correction))
    histname = "JetPerformance/hHadCaloEnergyMatchedHadCorr";
    fHistManager.FillTH2(histname, jetPt, hadCaloEnergyMatchedHadCorr);
    
    // Loop through particle types, and plot jet composition for each particle type
    histname = "JetPerformanceMC/hJetComposition";
    for (Int_t type = 0; type < 8; type++) {
      
      ContributorType contributorType = kUndefined;
      Double_t nSum = 0;
      Double_t pTsum = 0;
      
      // Loop through clusters in jet, and add to sum if cluster matches particle type
      for (Int_t iClus = 0; iClus < nClusters; iClus++) {
        
        clus = jet->Cluster(iClus);
        
        Int_t label = TMath::Abs(clus->GetLabel());
        if (label > 0) {
          contributorType = GetContributorType(clus, mcevent, label);
        }
        
        if (type == contributorType) {
          nSum++;
          pTsum += clus->GetNonLinCorrEnergy();
        }
      }
      
      // Fill jet composition histogram
      Double_t x[4] = {jetPt, 1.*type, nSum, pTsum};
      fHistManager.FillTHnSparse(histname, x);
      
    }
  }
}

/**
 * This function performs a study of the heavy-ion background.
 */
void AliAnalysisTaskEmcalJetPerformance::ComputeBackground()
{
  // Loop over tracks and clusters in order to:
  //   (1) Compute scale factor for full jets
  //   (2) Compute delta-pT for full jets, with the random cone method
  
  // Define the acceptance boundaries for the TPC and EMCal/DCal/PHOS
  Double_t etaTPC = 0.9;
  Double_t etaEMCal = 0.7;
  //Double_t etaMinDCal = 0.22;
  Double_t phiMinEMCal = fGeom->GetArm1PhiMin() * TMath::DegToRad(); // 80
  Double_t phiMaxEMCal = fGeom->GetEMCALPhiMax() * TMath::DegToRad(); // ~188
  //Double_t phiMinDCal = fGeom->GetDCALPhiMin() * TMath::DegToRad(); // 260
  //Double_t phiMaxDCal = fGeom->GetDCALPhiMax() * TMath::DegToRad(); // ~327 (1/3 SMs start at 320)
  
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
    Double_t etaEMCalfid = etaEMCal - jetR;
    Double_t phiMinEMCalfid = phiMinEMCal + jetR;
    Double_t phiMaxEMCalfid = phiMaxEMCal - jetR;
    
    // Generate EMCal random cone eta-phi
    Double_t etaEMCalRC = r->Uniform(-etaEMCalfid, etaEMCalfid);
    Double_t phiEMCalRC = r->Uniform(phiMinEMCalfid, phiMaxEMCalfid);
    
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
    // Note: Loops over all det-level track containers. For data there should be only one. For embedding, there should be signal+background tracks.
    AliParticleContainer * partCont = 0;
    AliTLorentzVector track;
    Double_t trackEta;
    Double_t trackPhi;
    Double_t trackPt;
    Double_t deltaR;
    TIter nextPartCont(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
      
      TString partContName = partCont->GetName();
      if (!partContName.CompareTo("tracks")) {
        
        AliTrackContainer* trackCont = dynamic_cast<AliTrackContainer*>(partCont);
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
        }
      }
    }
    
    // Loop over clusters. Sum the cluster ET:
    // (1) in the EMCal, (2) in the EMCal random cone
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
      
    }
    
    // Compute the scale factor for EMCal, as a function of centrality
    Double_t numerator = (trackPtSumEMCal + clusESumEMCal) / accEMCal;
    Double_t denominator = trackPtSumTPC / accTPC;
    Double_t scaleFactor = numerator / denominator;
    TString histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCal", jetCont->GetArrayName().Data());
    fHistManager.FillTH2(histname, fCent, scaleFactor);
    
    // Compute delta pT for EMCal, as a function of centrality
    Double_t rho = jetCont->GetRhoVal();
    Double_t deltaPt = trackPtSumEMCalRC + clusESumEMCalRC - rho * TMath::Pi() * jetR * jetR;
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCal", jetCont->GetArrayName().Data());
    fHistManager.FillTH2(histname, fCent, deltaPt);
    
    delete r;
    
  }
  
}
      
/**
 * This function performs a loop over the reconstructed jets, when the "simulated" trigger has been fired.
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetPerformance::FillTriggerSimHistograms()
{
  TString histname;
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    TString jetContName = jets->GetName();

    Double_t rhoVal = 0;
    if (jets->GetRhoParameter()) {
      rhoVal = jets->GetRhoVal();
      histname = TString::Format("%s/TriggerSimHistograms/hRhoVsCent", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), fCent, rhoVal);
    }
    
    for (auto jet : jets->all()) {
      
      Float_t ptLeading = jets->GetLeadingHadronPt(jet);
      Float_t corrPt = GetJetPt(jet, rhoVal);
      
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
      if ( (type != kEMCal) && (type != kDCal) ) {
        continue;
      }
      
      // (Centrality, pT, NEF, calo type)
      if (type == kEMCal) {
        histname = TString::Format("%s/TriggerSimHistograms/hNEFVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/TriggerSimHistograms/hNEFVsPtDCal", jets->GetArrayName().Data());
      }
      fHistManager.FillTH3(histname, fCent, corrPt, jet->NEF());
      
      // pT-leading vs. pT
      histname = TString::Format("%s/TriggerSimHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), corrPt, ptLeading);
      
      // (Centrality, pT, z-leading (charged), calo type)
      if (type == kEMCal) {
        histname = TString::Format("%s/TriggerSimHistograms/hZLeadingVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/TriggerSimHistograms/hZLeadingVsPtDCal", jets->GetArrayName().Data());
      }
      TLorentzVector leadPart;
      jets->GetLeadingHadronMomentum(leadPart, jet);
      Double_t z = GetParallelFraction(leadPart.Vect(), jet);
      if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin <1
      fHistManager.FillTH3(histname, fCent, corrPt, z);
      
      // (Centrality, pT, z (charged), calo type)
      if (type == kEMCal) {
        histname = TString::Format("%s/TriggerSimHistograms/hZVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/TriggerSimHistograms/hZVsPtDCal", jets->GetArrayName().Data());
      }
      const AliVTrack* track;
      for (Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
        track = static_cast<AliVTrack*>(jet->Track(i));
        z = track->Pt() / TMath::Abs(corrPt);
        fHistManager.FillTH3(histname, fCent, corrPt, z);
      }
      
      // (Centrality, pT, Nconst)
      if (type == kEMCal) {
        histname = TString::Format("%s/TriggerSimHistograms/hNConstVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/TriggerSimHistograms/hNConstVsPtDCal", jets->GetArrayName().Data());
      }
      fHistManager.FillTH3(histname, fCent, corrPt, 1.*jet->GetNumberOfConstituents());
      
    } //jet loop
    
  }
}

/**
 * This function fills histograms for matched truth-det jets in the case of embedding.
 * The jet matching information must be previously filled by another task, such as AliJetResponseMaker.
 */
void AliAnalysisTaskEmcalJetPerformance::FillMatchedJetHistograms()
{
  TString histname;
  AliJetContainer* jets = 0;
  const AliEmcalJet* matchedJet = nullptr;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    TString jetContName = jets->GetName();
    
    // Only loop over jets in the detector-level jet container
    if (jetContName.Contains("mcparticles")) {
      continue;
    }
    
    Double_t rhoVal = 0;
    if (jets->GetRhoParameter()) {
      rhoVal = jets->GetRhoVal();
    }
    
    for (auto jet : jets->accepted()) {
      
      // Get the matched jet, if it exists
      matchedJet = jet->MatchedJet();
      if (!matchedJet) {
        continue;
      }
      
      // compute jet acceptance type
      Double_t type = GetJetType(jet);
      if ( (type != kEMCal) && (type != kDCal) ) {
        continue;
      }
      
      Float_t detPt = GetJetPt(jet, rhoVal);
      Float_t truthPt = matchedJet->Pt();
      
      // Fill response matrix (centrality, pT-truth, pT-det)
      if (type == kEMCal) {
        histname = "MatchedJetHistograms/hResponseMatrixEMCal";
      }
      else if (type == kDCal) {
        histname = "MatchedJetHistograms/hResponseMatrixDCal";
      }
      fHistManager.FillTH3(histname, fCent, truthPt, detPt);
      
      // Fill JES shift (centrality, pT-truth, (pT-det - pT-truth) / pT-truth)
      if (type == kEMCal) {
        histname = "MatchedJetHistograms/hJESshiftEMCal";
      }
      else if (type == kDCal) {
        histname = "MatchedJetHistograms/hJESshiftDCal";
      }
      fHistManager.FillTH3(histname, fCent, truthPt, (detPt-truthPt)/truthPt );
      
      // Fill NEF of det-level matched jets (centrality, pT-truth, NEF)
      histname = "MatchedJetHistograms/hNEFVsPt";
      fHistManager.FillTH3(histname, fCent, truthPt, jet->NEF());

      // Fill z-leading (charged) of det-level matched jets (centrality, pT-truth, z-leading)
      histname = "MatchedJetHistograms/hZLeadingVsPt";
      TLorentzVector leadPart;
      jets->GetLeadingHadronMomentum(leadPart, jet);
      Double_t z = GetParallelFraction(leadPart.Vect(), jet);
      if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin <1
      fHistManager.FillTH3(histname, fCent, truthPt, z);
      
      // Fill matching distance (centrality, pT-truth, R)
      histname = "MatchedJetHistograms/hMatchingDistance";
      fHistManager.FillTH3(histname, fCent, truthPt, jet->ClosestJetDistance());
      
    } //jet loop
  }
}

/*
 * Compute the MC particle type of a given cluster contributor, using the MC particle container
 */
AliAnalysisTaskEmcalJetPerformance::ContributorType AliAnalysisTaskEmcalJetPerformance::GetContributorType(const AliVCluster* clus, const AliMCEvent* mcevent, Int_t label)
{
  ContributorType contributorType = kUndefined;

  AliAODMCParticle *part = fGeneratorLevel->GetMCParticleWithLabel(label);
  if (part) {
  
    TString histname = "ClusterHistogramsMC/hClusterRejectionReasonMC";
    UInt_t rejectionReason = 0;
    if (!fGeneratorLevel->AcceptMCParticle(part, rejectionReason)) {
      fHistManager.FillTH2(histname, fGeneratorLevel->GetRejectionReasonBitPosition(rejectionReason), clus->GetNonLinCorrEnergy());
      return contributorType;
    }

    if (part->GetGeneratorIndex() == 0) { // generator index in cocktail
      
      // select charged pions, protons, kaons, electrons, muons
      Int_t pdg = part->PdgCode();
      
      if (pdg == 22) { // gamma 22
        contributorType = kPhoton;
      }
      else if (TMath::Abs(pdg) == 211) { // pi+ 211 (abs value ensures both particles and antiparticles are included)
        contributorType = kChargedPion;
      }
      else if (pdg == 2212) { // proton 2212
        contributorType = kProton;
      }
      else if (pdg == -2212) {
        contributorType = kAntiProton;
      }
      else if (TMath::Abs(pdg) == 321) {  // K+ 321
        contributorType = kChargedKaon;
      }
      else if (pdg == 130) {  // K0L 130
        contributorType = kK0L;
      }
      else if (pdg == 2112) { // neutron 2112
        contributorType = kNeutron;
      }
      else if (pdg == -2112) {
        contributorType = kAntiNeutron;
      }
      else if (TMath::Abs(pdg) == 11) { // e- 11
        contributorType = kElectron;
      }
      else if (TMath::Abs(pdg) == 13) { // mu- 13
        contributorType = kMuon;
      }
      else {
        contributorType = kOther;
      }
    }
  }
  return contributorType;
}

/**
 * Get deltaR of a track/cluster and a reference point.
 */
Double_t AliAnalysisTaskEmcalJetPerformance::GetDeltaR(const AliTLorentzVector* part, Double_t etaRef, Double_t phiRef)
{
  Double_t deltaPhi = TMath::Abs(part->Phi_0_2pi() - phiRef);
  Double_t deltaEta = TMath::Abs(part->Eta() - etaRef);
  Double_t deltaR = TMath::Sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta );
  return deltaR;
}

/**
 * Get calo acceptance type of jet
 */
Double_t AliAnalysisTaskEmcalJetPerformance::GetJetType(const AliEmcalJet* jet)
{
  UInt_t jetType = jet->GetJetAcceptanceType();
  Double_t type = -1;
  if (jetType & AliEmcalJet::kEMCAL) {
    type = kEMCal;
  }
  else if (jetType & AliEmcalJet::kDCALonly) {
    type = kDCal;
  }
  
  return type;
}

/**
 * Get pT of jet -- background subtracted
 */
Double_t AliAnalysisTaskEmcalJetPerformance::GetJetPt(const AliEmcalJet* jet, Double_t rho)
{
  Double_t pT = jet->Pt() - rho * jet->Area();
  return pT;
}

/**
 * Return whether a contributor is a stable hadron
 */
Bool_t AliAnalysisTaskEmcalJetPerformance::IsHadron(const ContributorType contributor)
{
  return (contributor == kChargedPion) || (contributor == kProton) || (contributor == kAntiProton) || (contributor == kChargedKaon) || (contributor == kK0L) || (contributor == kNeutron) || (contributor == kAntiNeutron);
}

/**
 * JetPerformance AddTask.
 */
AliAnalysisTaskEmcalJetPerformance* AliAnalysisTaskEmcalJetPerformance::AddTaskEmcalJetPerformance(
  const char *ntracks,
  const char *nclusters,
  const char *nGenLev,
  const Double_t minTrPt,
  const Double_t minClPt,
  const char *suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetPerformance", "No analysis manager to connect to.");
    return 0;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetPerformance", "This task requires an input event handler");
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
  
  TString name("AliAnalysisTaskEmcalJetPerformance");
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
  // Configure jet performance task
  AliAnalysisTaskEmcalJetPerformance* task = new AliAnalysisTaskEmcalJetPerformance(name);
  
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
  if (partCont) task->AdoptParticleContainer(partCont);
  
  // Add the generator-level container, if specified
  if (nGenLev && strcmp(nGenLev,"")!=0) {
    AliMCParticleContainer* mcpartCont = task->AddMCParticleContainer(nGenLev);
    mcpartCont->SelectPhysicalPrimaries(kTRUE);
    mcpartCont->SetParticlePtCut(0);
  }
  
  // Add the cluster container
  AliClusterContainer* clusCont = 0;
  if (!clusName.IsNull()) {
    clusCont = new AliClusterContainer(clusName);
    clusCont->SetClusECut(0.);
    clusCont->SetClusPtCut(0.);
  }
  if (clusCont) task->AdoptClusterContainer(clusCont);
  
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                            TList::Class(),AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput1 );
  mgr->ConnectOutput (task, 1, coutput1 );
  
  return task;
}

