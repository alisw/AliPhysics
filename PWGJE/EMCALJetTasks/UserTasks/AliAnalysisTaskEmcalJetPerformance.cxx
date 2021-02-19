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
  fPlotCellNonlinearityHistograms(kFALSE),
  fPlotParticleCompositionHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoTriggerSimulation(kFALSE),
  fPlotMatchedJetHistograms(kFALSE),
  fComputeMBDownscaling(kFALSE),
  fPlotDCal(kFALSE),
  fDoClosureTest(kFALSE),
  fFillChargedFluctuations(kFALSE),
  fMinPt(-100),
  fMaxPt(250),
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
  fDoTriggerResponse(kFALSE),
  fDoJetMatchingGeometrical(kFALSE),
  fDoJetMatchingMCFraction(kFALSE),
  fDoDifferentialRM(kFALSE),
  fEmbeddingQA(),
  fMCJetContainer(nullptr),
  fDetJetContainer(nullptr),
  fDetJetContainerPPIntermediate(nullptr),
  fRequireMatchedJetAccepted(kFALSE),
  fJetMatchingR(0.),
  fMinSharedMomentumFraction(0.5),
  fMCJetMinMatchingPt(0.15),
  fDetJetMinMatchingPt(0.15),
  fPlotJetMatchCandThresh(1.),
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
  fPlotCellNonlinearityHistograms(kFALSE),
  fPlotParticleCompositionHistograms(kFALSE),
  fComputeBackground(kFALSE),
  fDoTriggerSimulation(kFALSE),
  fPlotMatchedJetHistograms(kFALSE),
  fComputeMBDownscaling(kFALSE),
  fPlotDCal(kFALSE),
  fDoClosureTest(kFALSE),
  fFillChargedFluctuations(kFALSE),
  fMinPt(-100),
  fMaxPt(250),
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
  fDoTriggerResponse(kFALSE),
  fDoJetMatchingGeometrical(kFALSE),
  fDoJetMatchingMCFraction(kFALSE),
  fDoDifferentialRM(kFALSE),
  fEmbeddingQA(),
  fMCJetContainer(nullptr),
  fDetJetContainer(nullptr),
  fDetJetContainerPPIntermediate(nullptr),
  fRequireMatchedJetAccepted(kFALSE),
  fJetMatchingR(0.),
  fMinSharedMomentumFraction(0.5),
  fMCJetMinMatchingPt(0.15),
  fPlotJetMatchCandThresh(1.),
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
      fEventCuts.fMC = false;
      fEventCuts.SetupLHC15o();
      fEventCuts.fUseVariablesCorrelationCuts = true;
    }
    fEventCuts.AddQAplotsToList(fEventCutList);
    fOutput->Add(fEventCutList);
  }
  
  // Get the MC particle branch, in case it exists
  fGeneratorLevel = GetMCParticleContainer("mcparticles");
  
  // Get jet containers, in order to check the jet acceptance criteria
  // For geometrical matching, it is expected that there is one det-level and one truth-level container
  // For MC fraction matching, it is expected that there is one Pb-Pb det-level, one pp det-level, and one pp truth-level container
  if (fPlotMatchedJetHistograms) {
    
    if (fDoJetMatchingGeometrical) {
      Printf("Geometrical jet matching enabled.");
      
      for (Int_t i=0; i<2; i++) {
        auto jetCont = GetJetContainer(i);
        TString jetContName = jetCont->GetName();
        if (jetContName.Contains("mcparticles")) {
          fMCJetContainer = jetCont;
        }
        else {
          fDetJetContainer = jetCont;
        }
      }
      
    }
    else if (fDoJetMatchingMCFraction) {
      Printf("MC-fraction jet matching enabled.");
      
      for (Int_t i=0; i<3; i++) {
        auto jetCont = GetJetContainer(i);
        TString jetContName = jetCont->GetName();
        if (jetContName.Contains("mcparticles")) {
          fMCJetContainer = jetCont;
        }
        else if (jetContName.Contains("Combined")) {
          fDetJetContainer = jetCont;
        }
        else {
          fDetJetContainerPPIntermediate = jetCont;
        }
      }
      
    }
    
    if (!fMCJetContainer) {
      Printf("No MC jet container found!");
    }
    Printf("mcJetContainer: %s", fMCJetContainer->GetName());
    
    if (!fDetJetContainer) {
      Printf("No det-level jet container found!");
    }
    Printf("det-level JetContainer: %s", fDetJetContainer->GetName());
    
    if (fDoJetMatchingMCFraction) {
      if (!fDetJetContainerPPIntermediate) {
        Printf("No intermediate pp det-level jet container found, despite MC-fraction matching enabled!");
      }
      Printf("Intermediate pp det-level JetContainer: %s", fDetJetContainerPPIntermediate->GetName());
    }
  
  }
  
  // Allocate histograms
  if (fPlotJetHistograms) {
    AllocateJetHistograms();
  }
  if (fPlotClusterHistograms) {
    AllocateClusterHistograms();
  }
  if (fPlotCellNonlinearityHistograms) {
    AllocateCellNonlinearityHistograms();
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
  
  Int_t nPtBins1 = TMath::CeilNint(fMaxPt-fMinPt);
  Int_t nPtBins2 = TMath::CeilNint((fMaxPt-fMinPt)/2);
  Int_t nPtBins5 = TMath::CeilNint((fMaxPt-fMinPt)/5);
  
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    if (fPlotMatchedJetHistograms && fDoJetMatchingMCFraction && jets != fMCJetContainer) {
      if (!fDoClosureTest) {
        continue; // don't plot det-level histograms if embedding
      }
    }
    
    // Jet rejection reason
    histname = TString::Format("%s/JetHistograms/hJetRejectionReason", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);Rejection reason;#it{p}_{T,corr} (GeV/#it{c});counts";
      TH3* hist = fHistManager.CreateTH3(histname.Data(), title.Data(), 10, 0, 100, 32, 0, 32, nPtBins5, fMinPt, fMaxPt);
      SetRejectionReasonLabels(hist->GetYaxis());
    } else {
      title = histname + ";Rejection reason;#it{p}_{T,corr} (GeV/#it{c});counts";
      TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, nPtBins5, fMinPt, fMaxPt);
      SetRejectionReasonLabels(hist->GetXaxis());
    }
    
    // Rho vs. Centrality
    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/JetHistograms/hRhoVsCent", jets->GetArrayName().Data());
      if (fForceBeamType == kAA) {
    	  title = histname + ";Centrality (%);#rho (GeV/#it{c});counts";
    	  fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 500);
      }
      else{
    	  title = histname + ";Centrality (%);#rho (GeV/#it{c});counts";
    	  fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 250, 0, 50);
      }
    }
    
    // (Centrality, pT, NEF)
    Int_t nbinsx = 20; Int_t minx = 0; Int_t maxx = 100;
    Int_t nbinsy = nPtBins1; Int_t miny = fMinPt; Int_t maxy = fMaxPt;
    Int_t nbinsz = 50; Int_t minz = 0; Int_t maxz = 1.;
    
    histname = TString::Format("%s/JetHistograms/hNEFVsPtEMCal", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});NEF";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    else {
      title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});NEF";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    
    if (fPlotDCal) {
      histname = TString::Format("%s/JetHistograms/hNEFVsPtDCal", jets->GetArrayName().Data());
      if (fForceBeamType == kAA) {
        title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});NEF";
        fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
      }
      else {
        title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});NEF";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
      }
    }
    
    // (Centrality, pT upscaled, calo type)
    if (fComputeMBDownscaling) {
      histname = TString::Format("%s/JetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});type";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 20, 0, 100, nPtBins2, fMinPt, fMaxPt, 2, -0.5, 1.5, "s");
    }
    
    // pT-leading vs. pT
    histname = TString::Format("%s/JetHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c})";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 10, 0, 100, nPtBins1, fMinPt, fMaxPt, fMaxPt, 0, fMaxPt);
    }
    else {
      title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c})";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins2, fMinPt, fMaxPt, nPtBins1, fMinPt, fMaxPt);
    }
    
    // A vs. pT
    histname = TString::Format("%s/JetHistograms/hAreaVsPt", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}/#pi#it{R}^{2}";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 10, 0, 100, nPtBins2, fMinPt, fMaxPt, 75, 0, 3);
    }
    else {
      title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}/#pi#it{R}^{2}";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins2, fMinPt, fMaxPt, 100, 0, 3);
    }
    
    // (Centrality, pT, z-leading (charged))
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins2; miny = fMinPt; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    
    histname = TString::Format("%s/JetHistograms/hZLeadingVsPtEMCal", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    else {
      title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    
    if (fPlotDCal) {
      histname = TString::Format("%s/JetHistograms/hZLeadingVsPtDCal", jets->GetArrayName().Data());
      if (fForceBeamType == kAA) {
        title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
        fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
      }
      else {
        title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
      }
    }
    
    // (Centrality, pT, z (charged))
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins2; miny = fMinPt; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    
    histname = TString::Format("%s/JetHistograms/hZVsPtEMCal", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    else {
      title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    
    if (fPlotDCal) {
      histname = TString::Format("%s/JetHistograms/hZVsPtDCal", jets->GetArrayName().Data());
      if (fForceBeamType == kAA) {
        title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
        fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
      }
      else {
        title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
      }
    }
    
    // (Centrality, pT, Nconst)
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins2; miny = fMinPt; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = fMaxPt;
    
    histname = TString::Format("%s/JetHistograms/hNConstVsPtEMCal", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    else {
      title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, nbinsz);
    }
    
    if (fPlotDCal) {
      histname = TString::Format("%s/JetHistograms/hNConstVsPtDCal", jets->GetArrayName().Data());
      if (fForceBeamType == kAA) {
        title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
        fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
      }
      else {
        title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
        fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, nbinsz);
      }
    }
    
    // (Centrality, pT) for eta<0 and eta>0
    if (fForceBeamType == kAA) {
      nbinsx = 20; minx = 0; maxx = 100;
      nbinsy = nPtBins1; miny = fMinPt; maxy = fMaxPt;
      
      histname = TString::Format("%s/JetHistograms/hEtaPosVsPtEMCal", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c})";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy);
      
      histname = TString::Format("%s/JetHistograms/hEtaNegVsPtEMCal", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c})";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy);
    }
    
    // (Centrality, jet pT, Enonlincorr - Ehadcorr)
    if (fForceBeamType == kAA) {
      nbinsx = 20; minx = 0; maxx = 100;
      nbinsy = nPtBins2; miny = fMinPt; maxy = fMaxPt;
      nbinsz = nPtBins2; minz = fMinPt; maxz = fMaxPt;
      
      histname = TString::Format("%s/JetHistograms/hDeltaEHadCorr", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#sum#it{E}_{nonlincorr} - #it{E}_{hadcorr}";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    
    // (Median patch energy, calo type, jet pT, centrality)
    if (fDoTriggerSimulation) {
      histname = TString::Format("%s/JetHistograms/hMedPatchJet", jets->GetArrayName().Data());
      title = histname + ";#it{E}_{patch,med};type;#it{p}_{T}^{corr} (GeV/#it{c});Centrality (%)";
      Int_t nbins5[4]  = {100, 2, nPtBins2, 50};
      Double_t min5[4] = {0,-0.5, fMinPt, 0};
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
 * This function allocates the histograms for the cell-level non-linearity study for embedding.
 * This should be run over MC.
 */
void AliAnalysisTaskEmcalJetPerformance::AllocateCellNonlinearityHistograms()
{
  TString histname;
  TString htitle;
  
  // Plot cluster non-linearity scale factor
  histname = "CellNonlinearityHistograms/hClusterNonlinearity";
  htitle = histname + ";E_{cluster}; E_{NonLinCorr}/E";
  fHistManager.CreateTH2(histname.Data(), htitle.Data(), 1500, 0, 150, 2000, 0.5, 1.5);
  
  // Plot cell non-linearity scale factor
  histname = "CellNonlinearityHistograms/hCellNonlinearity";
  htitle = histname + ";E_{cell}; E_{NonLinCorr}/E";
  fHistManager.CreateTH2(histname.Data(), htitle.Data(), 1500, 0, 150, 2000, 0.5, 1.5);
  
  // Plot cell non-linearity scale factor for leading cell
  histname = "CellNonlinearityHistograms/hLeadingCellNonlinearity";
  htitle = histname + ";E_{leading cell}; E_{NonLinCorr}/E";
  fHistManager.CreateTH2(histname.Data(), htitle.Data(), 1500, 0, 150, 2000, 0.5, 1.5);
  
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
    
    // EMCal
    histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality;Scale factor;counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 5);
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCal", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
    	title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
    	fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 400, -50, 150);
    }
    else{
    	title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
    	fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 1300, -25, 40);
    }
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCalExcl", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
    	title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
    	fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 400, -50, 150);
    }
    else{
    	title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
    	fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 1300, -25, 40);
    }
    if(fFillChargedFluctuations){
    	histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCalCharged", jets->GetArrayName().Data());
    	if (fForceBeamType == kAA) {
    	title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
    	fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 400, -50, 150);
    }
    else{
    	title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
    	fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 1300, -25, 40);
    }
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCalChargedExcl", jets->GetArrayName().Data());
    if (fForceBeamType == kAA) {
    	title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
    	fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 400, -50, 150);
    }
    else{
    	title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
    	fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 1300, -25, 40);
    }
    }
    histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCalFid", jets->GetArrayName().Data());
    title = histname + ";Centrality;Scale factor;counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 5);
    
    // DCal
    if (fPlotDCal) {
      histname = TString::Format("%s/BackgroundHistograms/hScaleFactorDCal", jets->GetArrayName().Data());
      title = histname + ";Centrality;Scale factor;counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 5);
      
      histname = TString::Format("%s/BackgroundHistograms/hDeltaPtDCal", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#delta#it{p}_{T} (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 10, 0, 100, 400, -50, 150);
      
      histname = TString::Format("%s/BackgroundHistograms/hScaleFactorDCalFid", jets->GetArrayName().Data());
      title = histname + ";Centrality;Scale factor;counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 5);
    }
    
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
  Int_t nPtBins1 = TMath::CeilNint(fMaxPt-fMinPt);
  Int_t nPtBins2 = TMath::CeilNint((fMaxPt-fMinPt)/2);
  Int_t nPtBins5 = TMath::CeilNint((fMaxPt-fMinPt)/5);
  
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
  fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, nPtBins2, 0, fMaxPt);
  
  // patch median vs. Centrality
  histname = "TriggerSimHistograms/hPatchMedianE";
  title = histname + ";Centrality (%);#it{E}_{patch,med} (GeV);type";
  fHistManager.CreateTH3(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 50, 2, -0.5, 1.5);
  
  if (fDoTriggerResponse) {
    histname = "TriggerSimHistograms/hMaxPatchResponseMatrix";
    title = histname + ";#it{E}_{maxpatch,det} (GeV);#it{E}_{maxpatch,truth} (GeV)";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 400, 0, 200, 400, 0, 200);
    
    histname = "TriggerSimHistograms/hEMCalEnergyResponseMatrix";
    title = histname + ";#Sigma#it{E}_{clus,EMCal} (GeV);#it{E}_{EMCal,truth} (GeV)";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 500, 0, 500, 500, 0, 500);
  }
  
  //----------------------------------------------
  // Jet histograms for "triggered" events
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    // Jet rejection reason
    histname = TString::Format("%s/TriggerSimHistograms/hJetRejectionReason", jets->GetArrayName().Data());
    title = histname + ";Rejection reason;#it{p}_{T,jet} (GeV/#it{c});counts";
    TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, nPtBins5, fMinPt, fMaxPt);
    SetRejectionReasonLabels(hist->GetXaxis());
    
    // Rho vs. Centrality
    if (!jets->GetRhoName().IsNull()) {
      histname = TString::Format("%s/TriggerSimHistograms/hRhoVsCent", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#rho (GeV/#it{c});counts";
      fHistManager.CreateTH2(histname.Data(), title.Data(), 50, 0, 100, 100, 0, 500);
    }
    
    // (Centrality, pT, NEF)
    Int_t nbinsx = 20; Int_t minx = 0; Int_t maxx = 100;
    Int_t nbinsy = nPtBins1; Int_t miny = fMinPt; Int_t maxy = fMaxPt;
    Int_t nbinsz = 50; Int_t minz = 0; Int_t maxz = 1.;
    
    histname = TString::Format("%s/TriggerSimHistograms/hNEFVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});NEF";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    if (fPlotDCal) {
      histname = TString::Format("%s/TriggerSimHistograms/hNEFVsPtDCal", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});NEF";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    
    // pT-leading vs. pT
    histname = TString::Format("%s/TriggerSimHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins2, fMinPt, fMaxPt, nPtBins2, fMinPt, fMaxPt);
    
    // A vs. pT
    histname = TString::Format("%s/TriggerSimHistograms/hAreaVsPt", jets->GetArrayName().Data());
    title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});#it{A}_{jet}/#pi#it{R}^{2}";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nPtBins2, fMinPt, fMaxPt, 75, 0, 3);
    
    // (Centrality, pT, z-leading (charged))
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins2; miny = fMinPt; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    
    histname = TString::Format("%s/TriggerSimHistograms/hZLeadingVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    if (fPlotDCal) {
      histname = TString::Format("%s/TriggerSimHistograms/hZLeadingVsPtDCal", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}_{leading}";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    
    // z (charged) vs. pT
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins2; miny = fMinPt; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = 1.;
    
    histname = TString::Format("%s/TriggerSimHistograms/hZVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    if (fPlotDCal) {
      histname = TString::Format("%s/TriggerSimHistograms/hZVsPtDCal", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});#it{z}";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    
    // (Centrality, pT, Nconst)
    nbinsx = 20; minx = 0; maxx = 100;
    nbinsy = nPtBins2; miny = fMinPt; maxy = fMaxPt;
    nbinsz = 50; minz = 0; maxz = fMaxPt;
    
    histname = TString::Format("%s/TriggerSimHistograms/hNConstVsPtEMCal", jets->GetArrayName().Data());
    title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    if (fPlotDCal) {
      histname = TString::Format("%s/TriggerSimHistograms/hNConstVsPtDCal", jets->GetArrayName().Data());
      title = histname + ";Centrality (%);#it{p}_{T}^{corr} (GeV/#it{c});No. of constituents";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    
  }
}

/*
 * This function allocates histograms for matched truth-det jets in the case of embedding.
 */
void AliAnalysisTaskEmcalJetPerformance::AllocateMatchedJetHistograms()
{
  TString histname;
  Double_t jetR=0;
  auto jetCont1 = GetJetContainer(0);
  jetR     = jetCont1->GetJetRadius();

  TString title;
  Int_t nPtBins1 = TMath::CeilNint(fMaxPt-fMinPt);
  Int_t nPtBinsTruth2 = TMath::CeilNint(fMaxPt/2);
  
  // Response matrix, (centrality, pT-truth, pT-det)
  Int_t nbinsx = 20; Int_t minx = 0; Int_t maxx = 100;
  Int_t nbinsy = fMaxPt; Int_t miny = 0; Int_t maxy = fMaxPt;
  Int_t nbinsz = nPtBins1; Int_t minz = fMinPt; Double_t maxz = fMaxPt;
  
  //This is a 5-dim RM with information on the angularity and matching distance
  if (fDoDifferentialRM) {
    //setup the THnSparse
    Int_t nCentBins=20;
    if (fForceBeamType != kAA) {
      nCentBins=1;
    }
    TString titleThn[6]= {"#it{p}_{T}^{truth} (GeV/#it{c})", "#it{p}_{T,corr}^{det} (GeV/#it{c})", "#Delta#it{R}", "shared mom fraction" ,"angularity", "Centrality (%)"};
    Int_t nbinsThn[6]  = {(Int_t)fMaxPt, (Int_t)fMaxPt, 15, 100, 100, nCentBins};
    Double_t minThn[6] = {0., 0., 0., 0.,0., 0.};
    Double_t maxThn[6] = {fMaxPt, fMaxPt, 1.5*jetR, 1, jetR, 100};
    histname = "MatchedJetHistograms/hResponseMatrixEMCalDiff";
    // (1) pt part LvL, (2) pt det LvL, (3) Matching distance (4) shared momentum fraction (5) angularity (6) centrality
    THnSparse* thn = fHistManager.CreateTHnSparse(histname.Data(), histname.Data(), 6, nbinsThn, minThn, maxThn);
    for (Int_t i = 0; i < 6; i++) {
      thn->GetAxis(i)->SetTitle(titleThn[i]);
      //thn->SetBinEdges(i, binEdges[i]);
    }
  }
  //This is a 3D RM for PbPb and a 2D RM for pp
  else {
    histname = "MatchedJetHistograms/hResponseMatrixEMCal";
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,corr}^{det} (GeV/#it{c})";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    else {
      title = histname + ";#it{p}_{T,corr}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c})";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsy, miny, maxy);
    }
  }
  
  if (fPlotDCal) {
    histname = "MatchedJetHistograms/hResponseMatrixDCal";
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,corr}^{det} (GeV/#it{c})";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    else {
      title = histname + ";#it{p}_{T}^{truth} (GeV/#it{c});#it{p}_{T,corr}^{det} (GeV/#it{c})";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
  }
  
  // JES shift, (centrality, pT-truth, (pT-det - pT-truth) / pT-truth)
  nbinsx = 20; minx = 0; maxx = 100;
  nbinsy = nPtBinsTruth2; miny = 0; maxy = fMaxPt;
  nbinsz = 250; minz = -5.; maxz = 5.;
  
  histname = "MatchedJetHistograms/hJESshiftEMCal";
  if (fForceBeamType == kAA) {
    title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#frac{#it{p}_{T,corr}^{det} - #it{p}_{T}^{truth}}{#it{p}_{T}^{truth}}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  }
  else {
    Int_t nbinsxR = 15;
    Int_t minxR = 0;
    Double_t maxxR = 1.5*jetR;
    title = histname + ";#Delta#it{R};#it{p}_{T}^{truth} (GeV/#it{c});#frac{#it{p}_{T,corr}^{det} - #it{p}_{T}^{truth}}{#it{p}_{T}^{truth}}";
    fHistManager.CreateTH3(histname.Data(), title.Data(),nbinsxR, minxR, maxxR, nbinsy, miny, maxy, nbinsz, minz, maxz);
  }
  
  if (fPlotDCal) {
    histname = "MatchedJetHistograms/hJESshiftDCal";
    if (fForceBeamType == kAA) {
      title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#frac{#it{p}_{T,corr}^{det} - #it{p}_{T}^{truth}}{#it{p}_{T}^{truth}}";
      fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
    else {
      title = histname + ";#it{p}_{T}^{truth} (GeV/#it{c});#frac{#it{p}_{T,corr}^{det} - #it{p}_{T}^{truth}}{#it{p}_{T}^{truth}}";
      fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
    }
  }
  
  // NEF of det-level matched jets, (centrality, pT-truth, NEF)
  nbinsx = 20; minx = 0; maxx = 100;
  nbinsy = fMaxPt; miny = 0; maxy = fMaxPt;
  nbinsz = 50; minz = 0; maxz = 1.;
  
  histname = "MatchedJetHistograms/hNEFVsPt";
  if (fForceBeamType == kAA) {
    title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});Calo energy fraction";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  }
  else {
    title = histname + ";#it{p}_{T}^{truth} (GeV/#it{c});Calo energy fraction";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
  }
  
  // z-leading (charged) of det-level matched jets, (centrality, pT-truth, z-leading)
  nbinsx = 20; minx = 0; maxx = 100;
  nbinsy = nPtBinsTruth2; miny = 0; maxy = fMaxPt;
  nbinsz = 50; minz = 0; maxz = 1.;
  
  histname = "MatchedJetHistograms/hZLeadingVsPt";
  if (fForceBeamType == kAA) {
    title = histname + ";Centrality (%);#it{p}_{T}^{truth} (GeV/#it{c});#it{z}_{leading}";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
  }
  else {
    title = histname + ";#it{p}_{T}^{truth} (GeV/#it{c});#it{z}_{leading}";
    fHistManager.CreateTH2(histname.Data(), title.Data(), nbinsy, miny, maxy, nbinsz, minz, maxz);
  }
  
  if (fDoJetMatchingGeometrical) {
   
    // Matching distance, (pT-det, pT-truth, deltaR)
    nbinsx = nPtBins1; minx = fMinPt; maxx = fMaxPt;
    nbinsy = nPtBinsTruth2; miny = 0; maxy = fMaxPt;
    nbinsz = 15; minz = 0; maxz = 1.5*jetR;
    histname = "MatchedJetHistograms/hMatchingDistance";
    title = histname + ";#it{p}_{T}^{det} (GeV/#it{c});#it{p}_{T}^{truth} (GeV/#it{c});R";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, fMaxPt, miny, maxy, nbinsz, minz, maxz);
    
    // Jet matching QA (copied from AliAnalysisTaskEmcalJetHCorrelations.cxx)
    histname = "MatchedJetHistograms/fHistJetMatchingQA";
    title = histname;
    std::vector<std::string> binLabels = {"noMatch", "matchedJet", "uniqueMatch", "jetDistance", "passedAllCuts"};
    auto histMatchedJetCuts = fHistManager.CreateTH1(histname.Data(), title.Data(), binLabels.size(), 0, binLabels.size());
    // Set label names
    for (unsigned int i = 1; i <= binLabels.size(); i++) {
      histMatchedJetCuts->GetXaxis()->SetBinLabel(i, binLabels.at(i-1).c_str());
    }
    histMatchedJetCuts->GetYaxis()->SetTitle("Number of jets");
  }
  
  if (fDoJetMatchingMCFraction) {
    
    // (det pT, shared MC fraction, deltaR) of closest jets
    nbinsx = nPtBins1; minx = fMinPt; maxx = fMaxPt;
    nbinsy = 20; miny = 0; maxy = 1.;
    nbinsz = 15; minz = 0; maxz = 1.5*jetR;
    histname = "MatchedJetHistograms/hMatchingDistanceVsMCFraction";
    title = histname + ";#it{p}_{T}^{det} (GeV/#it{c});MC fraction;#DeltaR";
    fHistManager.CreateTH3(histname.Data(), title.Data(), nbinsx, minx, maxx, nbinsy, miny, maxy, nbinsz, minz, maxz);
    
    // Jet matching QA (copied from AliAnalysisTaskEmcalJetHCorrelations.cxx)
    histname = "MatchedJetHistograms/fHistJetMatchingQA";
    title = histname;
    std::vector<std::string> binLabels = {"noMatch", "matchedJet", "sharedMomentumFraction", "partLevelMatchedJet", "jetDistancePPdet", "jetDistancePPtruth", "passedAllCuts"};
    auto histMatchedJetCuts = fHistManager.CreateTH1(histname.Data(), title.Data(), binLabels.size(), 0, binLabels.size());
    // Set label names
    for (unsigned int i = 1; i <= binLabels.size(); i++) {
      histMatchedJetCuts->GetXaxis()->SetBinLabel(i, binLabels.at(i-1).c_str());
    }
    histMatchedJetCuts->GetYaxis()->SetTitle("Number of jets");
  }

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
    PWG::EMCAL::AliEmcalDownscaleFactorsOCDB *downscaleOCDB = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
    downscaleOCDB->SetRun(InputEvent()->GetRunNumber());
    
    // There are two possible min bias triggers for LHC15o
    TString triggerNameMB1 = "CINT7-B-NOPF-CENT";
    TString triggerNameMB2 = "CV0L7-B-NOPF-CENT";
    TString triggerNameJE = "CINT7EJ1-B-NOPF-CENTNOPMD";
    
    // Get the downscale factor for whichever MB trigger exists in the given run
    std::vector<TString> runtriggers = downscaleOCDB->GetTriggerClasses();
    Double_t downscalefactor=1;
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
    Bool_t answer = AliAnalysisTaskEmcal::IsEventSelected();
    return answer;
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
  Double_t maxPatchEnergy = 0.;
  Double_t patchE;
  AliEMCALTriggerPatchInfo *maxPatch = nullptr;
  for(auto p : *fTriggerPatchInfo){
    AliEMCALTriggerPatchInfo *recpatch = static_cast<AliEMCALTriggerPatchInfo *>(p);
    if (recpatch) {
      
      if(!recpatch->IsJetHighSimple()) continue;
      
      patchE = recpatch->GetPatchE();
      
      histname = "TriggerSimHistograms/hEtaVsPhi";
      fHistManager.FillTH2(histname.Data(), recpatch->GetEtaGeo(), recpatch->GetPhiGeo());
      
      histname = "TriggerSimHistograms/hPatchE";
      fHistManager.FillTH2(histname.Data(), fCent, patchE);
      
      if (recpatch->IsEMCal()) {
        vecEMCal.push_back(patchE);
        if (patchE > maxPatchEnergy) {
          maxPatchEnergy = patchE;
          maxPatch = recpatch;
        }
      } else {
        vecDCal.push_back(patchE);
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
  
  // Fill max patch response and total EMCal energy response, if requested
  if (fGeneratorLevel && fDoTriggerResponse) {
    
    Double_t phiMinEMCal = fGeom->GetArm1PhiMin() * TMath::DegToRad(); // 80
    Double_t phiMaxEMCal = fGeom->GetEMCALPhiMax() * TMath::DegToRad(); // ~188
    Double_t detEnergyEMCal = 0;
    Double_t truthEnergyEMCal = 0;
    Double_t truthEnergyPatch = 0;
    
    AliTLorentzVector part;
    Double_t partEta;
    Double_t partPhi;
    Double_t partE;
    for (auto mcpart:fGeneratorLevel->accepted_momentum()) {
      part.Clear();
      part = mcpart.first;
      partEta = part.Eta();
      partPhi = part.Phi_0_2pi();
      partE = part.E();
      if (maxPatch) {
        if (partPhi < maxPatch->GetPhiMax() && partPhi > maxPatch->GetPhiMin()) {
          if (partEta < maxPatch->GetEtaMax() && partEta > maxPatch->GetEtaMin()) {
            truthEnergyPatch += partE;
          }
        }
      }
      if (partPhi < phiMaxEMCal && partPhi > phiMinEMCal) {
        if (partEta < 0.7 && partEta > -0.7) {
          truthEnergyEMCal += partE;
        }
      }
    }
    
    AliTLorentzVector clusFourVec;
    const AliVCluster* clus;
    Double_t clusEta;
    Double_t clusPhi;
    AliClusterContainer* clusters = GetClusterContainer(0);
    for (auto clusIterator : clusters->accepted_momentum() ) {
      clusFourVec.Clear();
      clusFourVec = clusIterator.first;
      clusEta = clusFourVec.Eta();
      clusPhi = clusFourVec.Phi_0_2pi();
      clus = clusIterator.second;
      if (clus->IsEMCAL()) {
        if (clusPhi < phiMaxEMCal && clusPhi > phiMinEMCal) {
          if (clusEta < 0.7 && clusEta > -0.7) {
            detEnergyEMCal += clusFourVec.E();
          }
        }
      }
    }
    
    if (maxPatch) {
      histname = "TriggerSimHistograms/hMaxPatchResponseMatrix";
      fHistManager.FillTH2(histname.Data(), maxPatchEnergy, truthEnergyPatch);
    }
    
    histname = "TriggerSimHistograms/hEMCalEnergyResponseMatrix";
    fHistManager.FillTH2(histname.Data(), detEnergyEMCal, truthEnergyEMCal);

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
  if (fPlotCellNonlinearityHistograms) {
    FillCellNonlinearityHistograms();
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
    
    if (fPlotMatchedJetHistograms && fDoJetMatchingMCFraction && jets != fMCJetContainer) {
      if (!fDoClosureTest) {
        continue; // don't plot det-level histograms if embedding
      }
    }
    
    Double_t rhoVal = 0;
    if (jets->GetRhoParameter()) {
      rhoVal = jets->GetRhoVal();
      histname = TString::Format("%s/JetHistograms/hRhoVsCent", jets->GetArrayName().Data());
      fHistManager.FillTH2(histname.Data(), fCent, rhoVal);
    }
    
    for (auto jet : jets->all()) {
      
      Float_t ptLeading = jets->GetLeadingHadronPt(jet);
      Float_t corrPt = GetJetPt(jet, rhoVal);
      Double_t jetR = jets->GetJetRadius();

      // compute jet acceptance type
      Double_t type = GetJetType(jet);
      if ( type != kEMCal ) {
        if ( type != kDCal || !fPlotDCal ) {
          continue;
        }
      }

      // (Centrality, Area, pT) (fill before area cut)
      histname = TString::Format("%s/JetHistograms/hAreaVsPt", jets->GetArrayName().Data());
      if (fForceBeamType == kAA) {
        fHistManager.FillTH3(histname.Data(), fCent, corrPt, 1.0*jet->Area()/(1.0*TMath::Pi()*pow(jetR,2)));
      }
      else {
        fHistManager.FillTH2(histname.Data(), corrPt, 1.0*jet->Area()/(1.0*TMath::Pi()*pow(jetR,2)));
      }
      
      // (Centrality, pT-leading, pT) (before leading hadron cuts)
      histname = TString::Format("%s/JetHistograms/hPtLeadingVsPt", jets->GetArrayName().Data());
      if (fForceBeamType == kAA) {
        fHistManager.FillTH3(histname.Data(), fCent, corrPt, ptLeading);
      }
      else {
        fHistManager.FillTH2(histname.Data(), corrPt, ptLeading);
      }
      
      // (Centrality, pT, z-leading (charged)) (before leading hadron cuts)
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
      if (fForceBeamType == kAA) {
        fHistManager.FillTH3(histname, fCent, corrPt, z);
      }
      else {
        fHistManager.FillTH2(histname, corrPt, z);
      }
      
      // Rejection reason
      UInt_t rejectionReason = 0;
      if (!jets->AcceptJet(jet, rejectionReason)) {
        histname = TString::Format("%s/JetHistograms/hJetRejectionReason", jets->GetArrayName().Data());
        if (fForceBeamType == kAA) {
          fHistManager.FillTH3(histname.Data(), fCent, jets->GetRejectionReasonBitPosition(rejectionReason), corrPt);
        }
        else {
          fHistManager.FillTH2(histname.Data(), jets->GetRejectionReasonBitPosition(rejectionReason), corrPt);
        }
        continue;
      }
      
      // (Centrality, pT, NEF)
      if (type == kEMCal) {
        histname = TString::Format("%s/JetHistograms/hNEFVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/JetHistograms/hNEFVsPtDCal", jets->GetArrayName().Data());
      }
      if (fForceBeamType == kAA) {
        fHistManager.FillTH3(histname, fCent, corrPt, jet->NEF());
      }
      else {
        fHistManager.FillTH2(histname, corrPt, jet->NEF());
      }
      
      // (Centrality, pT upscaled, calo type)
      if (fComputeMBDownscaling) {
        histname = TString::Format("%s/JetHistograms/hPtUpscaledMB", jets->GetArrayName().Data());
        fHistManager.FillTH3(histname.Data(), fCent, corrPt, type, fMBUpscaleFactor);
      }
      
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
        if (fForceBeamType == kAA) {
          fHistManager.FillTH3(histname, fCent, corrPt, z);
        }
        else {
          fHistManager.FillTH2(histname, corrPt, z);
        }
      }
 
      // (Centrality, pT, Nconst)
      if (type == kEMCal) {
        histname = TString::Format("%s/JetHistograms/hNConstVsPtEMCal", jets->GetArrayName().Data());
      }
      else if (type == kDCal) {
        histname = TString::Format("%s/JetHistograms/hNConstVsPtDCal", jets->GetArrayName().Data());
      }
      if (fForceBeamType == kAA) {
        fHistManager.FillTH3(histname, fCent, corrPt, 1.*jet->GetNumberOfConstituents());
      }
      else {
        fHistManager.FillTH2(histname, corrPt, 1.*jet->GetNumberOfConstituents());
      }
      
      // (Centrality, pT) for eta<0 and eta>0
      if (fForceBeamType == kAA) {
        if (type == kEMCal) {
          if (jet->Eta() > 0) {
            histname = TString::Format("%s/JetHistograms/hEtaPosVsPtEMCal", jets->GetArrayName().Data());
            fHistManager.FillTH2(histname, fCent, corrPt);
          }
          else if (jet->Eta() < 0) {
            histname = TString::Format("%s/JetHistograms/hEtaNegVsPtEMCal", jets->GetArrayName().Data());
            fHistManager.FillTH2(histname, fCent, corrPt);
          }
        }
      }
      
      // (Centrality, jet pT, Enonlincorr - Ehadcorr)
      if (fForceBeamType == kAA) {
        Double_t deltaEhadcorr = 0;
        const AliVCluster* clus = nullptr;
        Int_t nClusters = jet->GetNumberOfClusters();
        for (Int_t iClus = 0; iClus < nClusters; iClus++) {
          clus = jet->Cluster(iClus);
          deltaEhadcorr += (clus->GetNonLinCorrEnergy() - clus->GetHadCorrEnergy());
        }
        
        histname = TString::Format("%s/JetHistograms/hDeltaEHadCorr", jets->GetArrayName().Data());
        fHistManager.FillTH3(histname, fCent, corrPt, deltaEhadcorr);
      }
      
      
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

void AliAnalysisTaskEmcalJetPerformance::FillCellNonlinearityHistograms()
{
  // Get cells from event
  fCaloCells = InputEvent()->GetEMCALCells();
  
  // Loop through clusters
  AliClusterContainer* clusters = GetClusterContainer(0);
  AliTLorentzVector clusVec;
  const AliVCluster* clus;
  for (auto it : clusters->accepted_momentum()) {
    
    clusVec = it.first;
    clus = it.second;
    
    // Include only EMCal clusters
    if (!clus->IsEMCAL()) {
      continue;
    }
    if (clusVec.Phi_0_2pi() > 4.) {
      continue;
    }
  
    Double_t clusEnergy = clus->E();
    Double_t clusEnergyNonLinCorr = clus->GetNonLinCorrEnergy();
    Double_t correctionFactor = clusEnergyNonLinCorr/clusEnergy;
    
    // Plot cluster non-linearity scale factor
    TString histname = "CellNonlinearityHistograms/hClusterNonlinearity";
    fHistManager.FillTH2(histname.Data(), clusEnergy, correctionFactor);
    
    // Plot cell non-linearity scale factor for all cells
    histname = "CellNonlinearityHistograms/hCellNonlinearity";
    Double_t leadEcell = 0;
    for (Int_t iCell = 0; iCell < clus->GetNCells(); iCell++) {
      
      Int_t absId = clus->GetCellAbsId(iCell);
      Double_t eCell = fCaloCells->GetCellAmplitude(absId);
      if (eCell > leadEcell) {
        leadEcell = eCell;
      }

      fHistManager.FillTH2(histname.Data(), eCell, correctionFactor);
    }
    
    // Plot cell non-linearity scale factor for leading cell
    histname = "CellNonlinearityHistograms/hLeadingCellNonlinearity";
    fHistManager.FillTH2(histname.Data(), leadEcell, correctionFactor);
    
  }
    
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
      Double_t x[4] = {clus->GetM02(), clus->GetNonLinCorrEnergy(), static_cast<Double_t>(contributorType), jetPt};
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
  Double_t etaMinDCal = 0.22;
  Double_t phiMinEMCal = fGeom->GetArm1PhiMin() * TMath::DegToRad(); // 80
  Double_t phiMaxEMCal = fGeom->GetEMCALPhiMax() * TMath::DegToRad(); // ~188
  Double_t phiMinDCal = fGeom->GetDCALPhiMin() * TMath::DegToRad(); // 260
  Double_t phiMaxDCal = fGeom->GetDCALPhiMax() * TMath::DegToRad(); // ~327 (1/3 SMs start at 320)
  
  Double_t accTPC = 2 * etaTPC * 2 * TMath::Pi();
  Double_t accEMCal = 2 * etaEMCal * (phiMaxEMCal - phiMinEMCal);
  Double_t accDCal = 2 * (etaEMCal - etaMinDCal) * (phiMaxDCal - phiMinDCal);
  
  // Loop over jet containers
  AliJetContainer* jetCont = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(nextJetColl()))) {
	  if(jetCont->GetNJets()<1){
	    continue;//Don't enter if there are no jets in the container
	  }
	  Int_t maxJetIds[]   = {-1, -1};
	  Float_t maxJetPts[] = { 0,  0};
	  //loop over the jets in the jet container to find the leading and subleading jet
	  //alternative: AliEmcalJet* leadingJet =jetCont->GetLeadingJet(); but there is no GetSubLeadingJet() function
	  for (Int_t ij = 0; ij < jetCont->GetNJets(); ++ij) {
		//AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(ij));
		AliEmcalJet *jet = static_cast<AliEmcalJet*>(jetCont->GetJet(ij));
		//if (!AcceptJet(jet)) continue;
		if (jet->Pt() > maxJetPts[0]) {
		  maxJetPts[1] = maxJetPts[0];
		  maxJetIds[1] = maxJetIds[0];
		  maxJetPts[0] = jet->Pt();
		  maxJetIds[0] = ij;
		}
		else if (jet->Pt() > maxJetPts[1]) {
		  maxJetPts[1] = jet->Pt();
		  maxJetIds[1] = ij;
		}
	  }

    // Define fiducial acceptances, to be used to generate random cones, and for scale factor studies
    TRandom3* r = new TRandom3(0);
    Double_t jetR = jetCont->GetJetRadius();
    Double_t etaEMCalfid = etaEMCal - jetR;
    Double_t etaMinDCalfid = etaMinDCal + jetR;
    Double_t phiMinEMCalfid = phiMinEMCal + jetR;
    Double_t phiMaxEMCalfid = phiMaxEMCal - jetR;
    Double_t phiMinDCalfid = phiMinDCal + jetR;
    Double_t phiMaxDCalfid = phiMaxDCal - jetR;
    Double_t accEMCalfid = 2 * etaEMCalfid * (phiMaxEMCalfid - phiMinEMCalfid);
    Double_t accDCalfid = 2 * (etaEMCalfid - etaMinDCalfid) * (phiMaxDCalfid - phiMinDCalfid);
    if ( (etaEMCalfid - etaMinDCalfid) < 0) {
      accDCalfid = 0;
    }
    
    // Generate EMCal random cone eta-phi
    Double_t etaEMCalRC = r->Uniform(-etaEMCalfid, etaEMCalfid);
    Double_t phiEMCalRC = r->Uniform(phiMinEMCalfid, phiMaxEMCalfid);
    
    // Generate DCal random cone eta-phi
    Double_t etaDCalRC = r->Uniform(etaMinDCalfid, etaEMCalfid);
    Double_t sign = r->Uniform(-1., 1.);
    if (sign < 0) {
      etaDCalRC = -1*etaDCalRC;
    }
    Double_t phiDCalRC = r->Uniform(phiMinDCalfid, phiMaxDCalfid);

    // Initialize the various sums to 0
    Double_t trackPtSumTPC = 0;
    Double_t trackPtSumEMCal = 0;
    Double_t trackPtSumEMCalfid = 0;
    Double_t trackPtSumDCal = 0;
    Double_t trackPtSumDCalfid = 0;
    Double_t trackPtSumEMCalRC = 0;
    Double_t trackPtSumDCalRC = 0;
    Double_t clusPtSumEMCal = 0;
    Double_t clusPtSumEMCalfid = 0;
    Double_t clusPtSumDCal = 0;
    Double_t clusPtSumDCalfid = 0;
    Double_t clusPtSumEMCalRC = 0;
    Double_t clusPtSumDCalRC = 0;
    Bool_t overlap=0;
    // Loop over tracks. Sum the track pT:
    // (1) in the entire TPC, (2) in the EMCal, (3) in the EMCal fiducial volume, (4) in the DCal, (5) in the DCal fiducial volume, (6) in the EMCal random cone, (7) in the DCal random cone
    // Note: Loops over all det-level track containers. For data there should be only one. For embedding, there should be signal+background tracks.
    AliParticleContainer * partCont = 0;
    AliTLorentzVector trackVec;
    AliVParticle* track;
    Double_t trackEta;
    Double_t trackPhi;
    Double_t trackPt;
    Double_t deltaR;
    Double_t trackID=-1;
    TIter nextPartCont(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
      
      TString partContName = partCont->GetName();
      if (!partContName.CompareTo("tracks") || !partContName.CompareTo("thermalparticles")) {
        for (auto trackIterator : partCont->accepted_momentum() ) {
          
          trackVec.Clear();
          trackVec = trackIterator.first;
          trackEta = trackVec.Eta();
          trackPhi = trackVec.Phi_0_2pi();
          trackPt  = trackVec.Pt();
          //To get track ID and particle pointer
          track = trackIterator.second;
          TClonesArray* fTracksContArray = partCont->GetArray();
          trackID = fTracksContArray->IndexOf(track);

          // (1)
          if (TMath::Abs(trackEta) < etaTPC) {
            trackPtSumTPC += trackPt;
          }
          
          // (2)
          if (TMath::Abs(trackEta) < etaEMCal && trackPhi > phiMinEMCal && trackPhi < phiMaxEMCal) {
            trackPtSumEMCal += trackPt;
          }
          
          // (3)
          if (TMath::Abs(trackEta) < etaEMCalfid && trackPhi > phiMinEMCalfid && trackPhi < phiMaxEMCalfid) {
            trackPtSumEMCalfid += trackPt;
          }
          
          // (4)
          if (fPlotDCal) {
            if (TMath::Abs(trackEta) > etaMinDCal && TMath::Abs(trackEta) < etaEMCal && trackPhi > phiMinDCal && trackPhi < phiMaxDCal) {
              trackPtSumDCal += trackPt;
            }
          }
          
          // (5)
          if (fPlotDCal) {
            if (TMath::Abs(trackEta) > etaMinDCalfid && TMath::Abs(trackEta) < etaEMCalfid && trackPhi > phiMinDCalfid && trackPhi < phiMaxDCalfid) {
              trackPtSumDCalfid += trackPt;
            }
          }
          
          // (6)
          deltaR = GetDeltaR(&trackVec, etaEMCalRC, phiEMCalRC);
          if (deltaR < jetR) {
            trackPtSumEMCalRC += trackPt;
            //Check if there is an overlap with a leading/subleading signal jet
            //if so set the flag to 1 and not take this random cone into accound
            if(IsSignalJetOverlap(1,trackID,jetCont,maxJetIds))
            {
            	  overlap=1;
            }
          }
          
          // (7)
          if (fPlotDCal) {
            deltaR = GetDeltaR(&trackVec, etaDCalRC, phiDCalRC);
            if (deltaR < jetR) {
              trackPtSumDCalRC += trackPt;
            }
          }
        }
      }
    }
    
    // Loop over clusters, if the jet container is for full jets. Sum the cluster ET:
    // (1) in the EMCal, (2) in the EMCal fiducial volume, (3) in the DCal, (4), in the DCal fiducial volume, (5) in the EMCal random cone, (6) in the DCal random cone
    TString jetContName = jetCont->GetName();
    if (!jetContName.Contains("Charged")) {
    	  AliClusterContainer* clusCont = GetClusterContainer(0);
      AliTLorentzVector clusVec;
      AliVCluster* clus;
      Double_t clusEta;
      Double_t clusPhi;
      Double_t clusPt;
      Double_t clusID=-1;
      for (auto clusIterator : clusCont->accepted_momentum() ) {

        clusVec.Clear();
        clusVec = clusIterator.first;
        clusEta = clusVec.Eta();
        clusPhi = clusVec.Phi_0_2pi();
        clusPt  = clusVec.Pt();
        clus    = clusIterator.second;
        TClonesArray* fClusContArray = clusCont->GetArray();
        clusID = fClusContArray->IndexOf(clus);

        // (1)
        if (TMath::Abs(clusEta) < etaEMCal && clusPhi > phiMinEMCal && clusPhi < phiMaxEMCal) {
          clusPtSumEMCal += clusPt;
        }
        
        // (2)
        if (TMath::Abs(clusEta) < etaEMCalfid && clusPhi > phiMinEMCalfid && clusPhi < phiMaxEMCalfid) {
          clusPtSumEMCalfid += clusPt;
        }
        
        // (3)
        if (fPlotDCal) {
          if (TMath::Abs(clusEta) > etaMinDCal && TMath::Abs(clusEta) < etaEMCal && clusPhi > phiMinDCal && clusPhi < phiMaxDCal) {
            clusPtSumDCal += clusPt;
          }
        }
        
        // (4)
        if (fPlotDCal) {
          if (TMath::Abs(clusEta) > etaMinDCalfid && TMath::Abs(clusEta) < etaEMCalfid && clusPhi > phiMinDCalfid && clusPhi < phiMaxDCalfid) {
            clusPtSumDCalfid += clusPt;
          }
        }
        
        // (5)
        deltaR = GetDeltaR(&clusVec, etaEMCalRC, phiEMCalRC);
        if (deltaR < jetR) {
          clusPtSumEMCalRC += clusPt;
          //Check if there is an overlap with a leading/subleading signal jet
          //if so set the flag to 1 and not take this random cone into accound
          if(IsSignalJetOverlap(0,clusID,jetCont,maxJetIds))
          {
        	  overlap=1;
          }
        }
        // (6)
        if (fPlotDCal) {
          deltaR = GetDeltaR(&clusVec, etaDCalRC, phiDCalRC);
          if (deltaR < jetR) {
            clusPtSumDCalRC += clusPt;
          }
        }
      }
    }
    
    // Compute the scale factor, as a function of centrality, for (1) EMCal, (2) EMCalfid, (3) DCal, (4) DCalfid
    // (1)
    Double_t numerator = (trackPtSumEMCal + clusPtSumEMCal) / accEMCal;
    Double_t denominator = trackPtSumTPC / accTPC;
    Double_t scaleFactor = numerator / denominator;
    TString histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCal", jetCont->GetArrayName().Data());
    fHistManager.FillTH2(histname, fCent, scaleFactor);
    
    // (2)
    if (accEMCalfid > 1e-3) {
      numerator = (trackPtSumEMCalfid + clusPtSumEMCalfid) / accEMCalfid;
      scaleFactor = numerator / denominator;
      histname = TString::Format("%s/BackgroundHistograms/hScaleFactorEMCalFid", jetCont->GetArrayName().Data());
      fHistManager.FillTH2(histname, fCent, scaleFactor);
    }
    
    // (3)
    if (fPlotDCal) {
      numerator = (trackPtSumDCal + clusPtSumDCal) / accDCal;
      scaleFactor = numerator / denominator;
      histname = TString::Format("%s/BackgroundHistograms/hScaleFactorDCal", jetCont->GetArrayName().Data());
      fHistManager.FillTH2(histname, fCent, scaleFactor);
    }
    
    // (4)
    if (fPlotDCal) {
      if (accDCalfid > 1e-3) {
        numerator = (trackPtSumDCalfid + clusPtSumDCalfid) / accDCalfid;
        scaleFactor = numerator / denominator;
        histname = TString::Format("%s/BackgroundHistograms/hScaleFactorDCalFid", jetCont->GetArrayName().Data());
        fHistManager.FillTH2(histname, fCent, scaleFactor);
      }
    }
    
    // Compute delta pT, as a function of centrality
    // EMCal acceptance only charged component
    Double_t rho = jetCont->GetRhoVal();
    Double_t deltaPt = trackPtSumEMCalRC - rho * TMath::Pi() * jetR * jetR;
    if(fFillChargedFluctuations){
    	histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCalCharged", jetCont->GetArrayName().Data());
    	fHistManager.FillTH2(histname, fCent, deltaPt);

    	// EMCal acceptance only charged component
    	// excluding random cones that overlap with signal jets
    	histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCalChargedExcl", jetCont->GetArrayName().Data());
    	if(overlap==0)fHistManager.FillTH2(histname, fCent, deltaPt);
    }
    // EMCal acceptance charged and neutral component
    deltaPt += clusPtSumEMCalRC;
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCal", jetCont->GetArrayName().Data());
    fHistManager.FillTH2(histname, fCent, deltaPt);

    // EMCal acceptance charged and neutral component
    // excluding random cones that overlap with signal jets
    histname = TString::Format("%s/BackgroundHistograms/hDeltaPtEMCalExcl", jetCont->GetArrayName().Data());
    if(overlap==0)fHistManager.FillTH2(histname, fCent, deltaPt);
    
    // DCal
    if (fPlotDCal) {
      if (accDCalfid > 1e-3) {
        deltaPt = trackPtSumDCalRC + clusPtSumDCalRC - rho * TMath::Pi() * jetR * jetR;
        histname = TString::Format("%s/BackgroundHistograms/hDeltaPtDCal", jetCont->GetArrayName().Data());
        fHistManager.FillTH2(histname, fCent, deltaPt);
      }
    }
    
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
    Double_t jetR = jets->GetJetRadius();

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
      fHistManager.FillTH2(histname.Data(), corrPt, 1.0*jet->Area()/(1.0*TMath::Pi()*pow(jetR,2)));
      
      // Rejection reason
      UInt_t rejectionReason = 0;
      if (!jets->AcceptJet(jet, rejectionReason)) {
        histname = TString::Format("%s/TriggerSimHistograms/hJetRejectionReason", jets->GetArrayName().Data());
        fHistManager.FillTH2(histname.Data(), jets->GetRejectionReasonBitPosition(rejectionReason), jet->Pt());
        continue;
      }
      
      // compute jet acceptance type
      Double_t type = GetJetType(jet);
      if ( type != kEMCal ) {
        if ( type != kDCal || !fPlotDCal ) {
          continue;
        }
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
 * This function performs matching of det-level jets to truth-level jets, and fills relevant histograms.
 * There are two matching approaches implemented, each of which expect a certain set of jet containers to be attached:
 * (1) fDoJetMatchingGeometrical: Do geometrical matching, with two jet containers fMCJetContainer and fDetJetContainer.
 *     This is appropriate for pp and p-Pb.
 * (2) fDoJetMatchingMCFraction: Do MC fraction based jet matching, with pp-truth, pp-det, and combined-det jet containers.
 *     This is appropriate for Pb-Pb.
 */
void AliAnalysisTaskEmcalJetPerformance::FillMatchedJetHistograms()
{
  // Loop over all jets and fill the ClosestJet(), i.e. the matching candidate.
  // Note: Allow truth jets to be outside of EMCALfid or fail 5 GeV requirement, since these can still contribute accepted det-jets
  //       (but for the jet reconstruction efficiency denominator, the criteria should be enforced).
  if (fDoJetMatchingGeometrical) {
    if (fRequireMatchedJetAccepted) {
      fMCJetContainer->SetJetAcceptanceType(AliEmcalJet::kTPC);
      ComputeJetMatches(fDetJetContainer, fMCJetContainer, kTRUE);
      fMCJetContainer->SetJetAcceptanceType(AliEmcalJet::kEMCALfid);
    }
    else {
      ComputeJetMatches(fDetJetContainer, fMCJetContainer, kFALSE);
    }
  }
  else if (fDoJetMatchingMCFraction) {
    // First match PbPb-det to pp-det
    ComputeJetMatches(fDetJetContainer, fDetJetContainerPPIntermediate, kTRUE);
    
    // Then match pp-det to pp-truth
    if (fRequireMatchedJetAccepted) { // Require pp-truth be accepted (i.e. leading track req), but still allow geometrical acceptance migration
      fMCJetContainer->SetJetAcceptanceType(AliEmcalJet::kTPC);
      ComputeJetMatches(fDetJetContainerPPIntermediate, fMCJetContainer, kTRUE);
      fMCJetContainer->SetJetAcceptanceType(AliEmcalJet::kEMCALfid);
    }
    else{ // Don't require pp-truth jet to be accepted
      ComputeJetMatches(fDetJetContainerPPIntermediate, fMCJetContainer, kFALSE);
    }

  }
  
  // Loop through accepted det-level jets, and retrieve matching candidate.
  // It match passes criteria (i.e. matching distance, uniqueness, MC fraction), fill matching histos.
  Double_t rhoVal = 0;
  if (fDetJetContainer->GetRhoParameter()) {
    rhoVal = fDetJetContainer->GetRhoVal();
  }
  for (auto jet : fDetJetContainer->accepted()) {
    
    Float_t detPt = GetJetPt(jet, rhoVal);
    
    // Get the matched part-level jet
    const AliEmcalJet* matchedPartLevelJet = GetMatchedPartLevelJet(jet, detPt);
    if (!matchedPartLevelJet) {
      continue;
    }
    Float_t truthPt = matchedPartLevelJet->Pt();
    
    // compute jet acceptance type
    Double_t type = GetJetType(jet);
    if ( type != kEMCal ) {
      if ( type != kDCal || !fPlotDCal ) {
        continue;
      }
    }
    
    // Fill response matrix (centrality, pT-truth, pT-det)
    TString histname;
    //This is a 5-dim RM with information on the angularity and matching distance
    if (fDoDifferentialRM) {
      histname = "MatchedJetHistograms/hResponseMatrixEMCalDiff";
      // (1) pt part LvL, (2) pt det LvL, (3) Matching distance (4) angularity (5) centrality
      Double_t angularity     = GetAngularity(matchedPartLevelJet);
      Double_t matchDistance  = matchedPartLevelJet->ClosestJetDistance();
      Double_t sharedFraction = fDetJetContainer->GetFractionSharedPt(jet, nullptr);
      Double_t x[6] = {truthPt, detPt, matchDistance, sharedFraction, angularity, fCent};
      fHistManager.FillTHnSparse(histname, x);
    }
    //This is a 3D RM for PbPb and a 2D RM for pp
    else {
      if (type == kEMCal) {
        histname = "MatchedJetHistograms/hResponseMatrixEMCal";
      }
      else if (type == kDCal) {
        histname = "MatchedJetHistograms/hResponseMatrixDCal";
      }

      if (fForceBeamType == kAA) {
        fHistManager.FillTH3(histname, fCent, truthPt, detPt);
      }
      else {
        fHistManager.FillTH2(histname, detPt, truthPt);
      }
    }
    // Fill JES shift (centrality, pT-truth, (pT-det - pT-truth) / pT-truth)
    if (type == kEMCal) {
      histname = "MatchedJetHistograms/hJESshiftEMCal";
    }
    else if (type == kDCal) {
      histname = "MatchedJetHistograms/hJESshiftDCal";
    }
    if (fForceBeamType == kAA) {
      fHistManager.FillTH3(histname, fCent, truthPt, (detPt-truthPt)/truthPt );
    }
    
    // Fill NEF of det-level matched jets (centrality, pT-truth, NEF)
    histname = "MatchedJetHistograms/hNEFVsPt";
    if (fForceBeamType == kAA) {
      fHistManager.FillTH3(histname, fCent, truthPt, jet->NEF());
    }
    else {
      fHistManager.FillTH2(histname, truthPt, jet->NEF());
    }

    // Fill z-leading (charged) of det-level matched jets (centrality, pT-truth, z-leading)
    histname = "MatchedJetHistograms/hZLeadingVsPt";
    TLorentzVector leadPart;
    fDetJetContainer->GetLeadingHadronMomentum(leadPart, jet);
    Double_t z = GetParallelFraction(leadPart.Vect(), jet);
    if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin <1
    if (fForceBeamType == kAA) {
      fHistManager.FillTH3(histname, fCent, truthPt, z);
    }
    else {
      fHistManager.FillTH2(histname, truthPt, z);
    }
    
  } //jet loop
}

/*
 * Loop over jets of two specified jet collections, and fill the ClosestJet(), i.e. the matching candidate.
 * The first collection always uses the container acceptance criteria.
 * The second collection can be configured to use the container acceptance criteria or not.
 */
void AliAnalysisTaskEmcalJetPerformance::ComputeJetMatches(AliJetContainer* jetCont1, AliJetContainer* jetCont2, Bool_t bUseJetCont2Acceptance) {

  for (auto jet1 : jetCont1->all()) {
    jet1->ResetMatching();
  }
  for (auto jet2 : jetCont2->all()) {
    jet2->ResetMatching();
  }

  for (auto jet1 : jetCont1->accepted()) {//detector level
    if (jet1->Pt() < fDetJetMinMatchingPt) {
        continue;
    }
    if (bUseJetCont2Acceptance) {
      for (auto jet2 : jetCont2->accepted()) {
        SetJetClosestCandidate(jet1, jet2);
      }
    }
    else {
      for (auto jet2 : jetCont2->all()) {//truth level
        if (jet2->Pt() < fMCJetMinMatchingPt) {
          continue;
        }
        SetJetClosestCandidate(jet1, jet2);
      }
    }
  }
}
          
/*
 * Given two jets, set them as closest if they are closer than the current closest jets.
 */
void AliAnalysisTaskEmcalJetPerformance::SetJetClosestCandidate(AliEmcalJet* jet1, AliEmcalJet* jet2) {
          
  Double_t deltaR = jet1->DeltaR(jet2);
  if (deltaR > 0.) {
    
    if (deltaR < jet1->ClosestJetDistance()) {
      jet1->SetClosestJet(jet2, deltaR);
    }
    
    if (deltaR < jet2->ClosestJetDistance()) {
      jet2->SetClosestJet(jet1, deltaR);
    }
  }
}
/*
 * Return a pointer to the matched truth-level jet, if it passes the matching criteria.
 * For fDoJetMatchingGeometrical, this means (1) within R = fJetMatchingR, (2) unique match.
 * For fDoJetMatchingMCFraction, this means also shared MC fraction requirement. That is, if the jet (combined jet)
 * is matched to a pp det-level jet, which is matched to a pp truth-level jet -- it must satisfy
 * (1) The shared momentum fraction being larger than some minimum value fMinSharedMomentumFraction, and
 * (2) Their matched distance being below the max matching distance fMaxMatchedJetDistance
 *
 * @param[in] detJet det-level jet to be checked for a successful truth-level match.
 * @return Pointer to truth-level matched jet, if it exists. False otherwise.
*/
const AliEmcalJet* AliAnalysisTaskEmcalJetPerformance::GetMatchedPartLevelJet(const AliEmcalJet* detJet, Double_t detJetPt) {

  // Track in histogram how many matches pass each distinct matching criteria
  TString histNameQA = "MatchedJetHistograms/fHistJetMatchingQA";
  
  // Geometrical matching case (pp, p-Pb)
  if (fDoJetMatchingGeometrical) {
    
    bool returnValue = false;
    const AliEmcalJet* partLevelJet = detJet->ClosestJet();
    if (partLevelJet) {
      fHistManager.FillTH1(histNameQA, "matchedJet");
      returnValue = true;
      
      // Check if match is unique
      if (partLevelJet->ClosestJet() != detJet) {
        returnValue = false;
      }
      else {
        Double_t truthPt=partLevelJet->Pt();
        fHistManager.FillTH1(histNameQA, "uniqueMatch");
        
        // Fill matching distance between unique matches, without imposing deltaR cut (centrality, pT-truth, R)
        TString histname = "MatchedJetHistograms/hMatchingDistance";
        fHistManager.FillTH3(histname, detJetPt, truthPt, detJet->ClosestJetDistance());
        if (fForceBeamType != kAA)
        {
          histname = "MatchedJetHistograms/hJESshiftEMCal";
          fHistManager.FillTH3(histname, detJet->ClosestJetDistance(), truthPt, (detJetPt-truthPt)/truthPt);
        }
      }
      
      // Check if the matching distance cut is passed
      double matchedJetDistance = detJet->ClosestJetDistance();
      if (matchedJetDistance > fJetMatchingR) {
        returnValue = false;
      }
      else {
        fHistManager.FillTH1(histNameQA, "jetDistance");
      }
      
      // Record all cuts passed
      if (returnValue == true) {
        fHistManager.FillTH1(histNameQA, "passedAllCuts");
      }
    }
    else {
      fHistManager.FillTH1(histNameQA, "noMatch");
      returnValue = false;
    }
    
    if (returnValue) {
      return partLevelJet;
    }
  }
  
  // Shared MC fraction matching case (Pb-Pb)
  else if (fDoJetMatchingMCFraction) { // This function is essentially copied from AliAnalysisTaskEmcalJetHCorrelations::CheckForMatchedJet
    
    bool returnValue = false;
    const AliEmcalJet* partLevelJet = nullptr;
    
    // First, check if combined jet has a pp det-level match assigned
    if (detJet->ClosestJet()) {
      fHistManager.FillTH1(histNameQA, "matchedJet");
      returnValue = true;
      
      // Check shared momentum fraction.
      double sharedFraction = fDetJetContainer->GetFractionSharedPt(detJet, nullptr);
      if (sharedFraction < fMinSharedMomentumFraction) {
        returnValue = false;
      }
      else {
        fHistManager.FillTH1(histNameQA, "sharedMomentumFraction");
      }
    
      // Check that the combined jet has a particle-level match
      AliEmcalJet * detLevelJetPP = detJet->ClosestJet();
      partLevelJet = detLevelJetPP->ClosestJet();
      if (!partLevelJet) {
        returnValue = false;
      }
      else {
        fHistManager.FillTH1(histNameQA, "partLevelMatchedJet");
      }
      
      // Check the matching distance between the combined and pp det-level jets
      double matchedJetDistance = detJet->ClosestJetDistance();
      if (matchedJetDistance > fJetMatchingR) {
        returnValue = false;
      }
      else {
        fHistManager.FillTH1(histNameQA, "jetDistancePPdet");
      }
      
      // Check the matching distance between the combined and pp truth-level jets
      if (partLevelJet) {
        Double_t deltaR = detJet->DeltaR(partLevelJet);
        if (deltaR > fJetMatchingR) {
          returnValue = false;
        }
        else {
          fHistManager.FillTH1(histNameQA, "jetDistancePPtruth");
        }
      }
      
      // Record all cuts passed
      if (returnValue == true) {
        fHistManager.FillTH1(histNameQA, "passedAllCuts");
      }
      
      // Fill (det pT, shared MC fraction, deltaR) of closest jets
      TString histname = "MatchedJetHistograms/hMatchingDistanceVsMCFraction";
      fHistManager.FillTH3(histname, detJetPt, sharedFraction, matchedJetDistance);
      
    }
    else {
      fHistManager.FillTH1(histNameQA, "noMatch");
      returnValue = false;
    }
    
    if (returnValue) {
      return partLevelJet;
    }
  }
  return 0;
}
/*
 * Compute the Angularity of a jet - based on particle tracks (for particle level info)
 */
Double_t AliAnalysisTaskEmcalJetPerformance::GetAngularity(const AliEmcalJet* jet)
{
  //
  Double_t angularity=-1;

  if (jet->GetNumberOfTracks()== 0) return 0;

  Double_t den =0.;
  Double_t num = 0.;
  AliVParticle *vp1 = 0x0;
  //loop over all tracks in the jet
  for(UInt_t i = 0; i < jet->GetNumberOfTracks(); i++) {

    vp1 = static_cast<AliVParticle*>(jet->Track(i));

    if (!vp1){
      Printf("AliVParticle associated to constituent not found");
      continue;
    }

    Double_t dphi = GetRelativePhi(vp1->Phi(),jet->Phi());
    Double_t dr2  = (vp1->Eta()-jet->Eta())*(vp1->Eta()-jet->Eta()) + dphi*dphi;
    Double_t dr   = TMath::Sqrt(dr2);
    num=num+vp1->Pt()*dr;
    den=den+vp1->Pt();
  }
  if (den>0) angularity=num/den;

  return angularity;
}
/*
 * Compute dPhi distance
 */
Double_t AliAnalysisTaskEmcalJetPerformance::GetRelativePhi(Double_t mphi, Double_t vphi){

     if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
     else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
     if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
     else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
     double dphi = mphi-vphi;
     if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
     else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
     return dphi;//dphi in [-Pi, Pi]
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
  if (jetType & AliEmcalJet::kEMCALfid) {
    type = kEMCal;
  }
  else if (jetType & AliEmcalJet::kDCALonlyfid) {
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
 * Return whether a cluster/track is part of a signal jet
 * Check if the particle in the random cone is part of
 * any of the two leading Anti-Kt jets in the event
 */
Bool_t AliAnalysisTaskEmcalJetPerformance::IsSignalJetOverlap(Bool_t isTrack, Int_t particleID, const AliJetContainer* jetCont, Int_t maxJetIds[])
{
	//loop over the two leading anti-kT jets in the event
	for (Int_t i=0; i<2; i++)
	{
		if(maxJetIds[i]==-1)continue;
		AliEmcalJet* jet= jetCont->GetJet(maxJetIds[i]);
		Int_t ClusterConstituentID=-1;
		Int_t TrackConstituentID  =-1;
		//check if the particle in the random cone is a constituent of a signal-Anti-KT jet
		if(isTrack)
		{
			TrackConstituentID   = jet->ContainsTrack(particleID);
		}
		else
		{
			ClusterConstituentID = jet->ContainsCluster(particleID);
		}
		if(ClusterConstituentID>-1 || TrackConstituentID>-1)
		{
			return kTRUE;
		}
	}
	return kFALSE;
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

