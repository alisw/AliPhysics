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
#include "AliPHOSGeometry.h"
#include "AliOADBContainer.h"

#include "AliAnalysisTaskEmcalVsPhos.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalVsPhos);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalVsPhos::AliAnalysisTaskEmcalVsPhos() : 
  AliAnalysisTaskEmcalJet(),
  fPlotClusterHistograms(kTRUE),
  fPlotNeutralJets(kFALSE),
  fPlotClustersInJets(kFALSE),
  fPlotCellHistograms(kTRUE),
  fPlotClusWithoutNonLinCorr(kFALSE),
  fPlotExotics(kFALSE),
  fPlotStandardClusterTHnSparse(kTRUE),
  fPlotNearestNeighborDistribution(kFALSE),
  fPlotClusterCone(kFALSE),
  fPlotCaloCentrality(kFALSE),
  fPlotFineGrainedEtaPhi(kFALSE),
  fPlotEvenOddEta(kFALSE),
  fPlotCellSMDensity(kFALSE),
  fExcludeRejectedCells(kFALSE),
  fPlotFineGrainedCentrality(kFALSE),
  fPlotEventHistograms(kFALSE),
  fMaxPt(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNM02HistBins(0),
  fM02HistBins(0),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fPHOSGeo(nullptr),
  fHistManager()
{
  GenerateHistoBins();
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalVsPhos::AliAnalysisTaskEmcalVsPhos(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fPlotClusterHistograms(kTRUE),
  fPlotNeutralJets(kFALSE),
  fPlotClustersInJets(kFALSE),
  fPlotCellHistograms(kTRUE),
  fPlotClusWithoutNonLinCorr(kFALSE),
  fPlotExotics(kFALSE),
  fPlotStandardClusterTHnSparse(kTRUE),
  fPlotNearestNeighborDistribution(kFALSE),
  fPlotClusterCone(kFALSE),
  fPlotCaloCentrality(kFALSE),
  fPlotFineGrainedEtaPhi(kFALSE),
  fPlotEvenOddEta(kFALSE),
  fPlotCellSMDensity(kFALSE),
  fExcludeRejectedCells(kFALSE),
  fPlotFineGrainedCentrality(kFALSE),
  fPlotEventHistograms(kFALSE),
  fMaxPt(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNM02HistBins(0),
  fM02HistBins(0),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fPHOSGeo(nullptr),
  fHistManager(name)
{
  GenerateHistoBins();
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalVsPhos::~AliAnalysisTaskEmcalVsPhos()
{
}

/**
 * Generate histogram binning arrays
 */
void AliAnalysisTaskEmcalVsPhos::GenerateHistoBins()
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
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalVsPhos::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  
  AllocateCaloHistograms();

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
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for the calorimeter performance study.
 */
void AliAnalysisTaskEmcalVsPhos::AllocateCaloHistograms()
{
  
  if (fPlotClusterHistograms) {
    AllocateClusterHistograms();
  }
  
  if (fPlotCellHistograms) {
    AllocateCellHistograms();
  }
  
  if (fPlotNeutralJets) {
    AllocateNeutralJetHistograms();
  }
  
  if (fPlotClustersInJets) {
    AllocateClustersInJetsHistograms();
  }
  
  if (fPlotEventHistograms) {
    AllocateEventHistograms();
  }
  
}

/*
 * This function allocates the histograms for the calorimeter performance study.
 */
void AliAnalysisTaskEmcalVsPhos::AllocateClusterHistograms()
{
  TString histname;
  TString htitle;
  
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
    
    histname = TString::Format("%s/hClusterRejectionReasonMC", cont->GetArrayName().Data());
    htitle = histname + ";Rejection reason;#it{E}_{clus} (GeV/)";
    TH2* histMC = fHistManager.CreateTH2(histname.Data(), htitle.Data(), nRejBins, rejReasonBins, fNPtHistBins, fPtHistBins);
    SetRejectionReasonLabels(histMC->GetXaxis());
    
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
  
    // Plot cluster THnSparse (centrality, eta, phi, E, E-hadcorr, has matched track, M02, Ncells, passed dispersion cut)
    Int_t dim = 0;
    TString title[20];
    Int_t nbins[20] = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    Double_t *binEdges[20] = {0};
    
    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
      title[dim] = "Centrality %";
      if (fPlotFineGrainedCentrality) {
        nbins[dim] = 18;
        min[dim] = 0;
        max[dim] = 90;
        binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      }
      else {
        nbins[dim] = fNCentHistBins;
        binEdges[dim] = fCentHistBins;
        min[dim] = fCentHistBins[0];
        max[dim] = fCentHistBins[fNCentHistBins];
      }
      dim++;
    }
    
    title[dim] = "#eta";
    nbins[dim] = 28;
    if (fPlotFineGrainedEtaPhi) {
      nbins[dim] = 100;
    }
    min[dim] = -0.7;
    max[dim] = 0.7;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    title[dim] = "#phi";
    nbins[dim] = 100;
    if (fPlotFineGrainedEtaPhi) {
      nbins[dim] = 357;
    }
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
    
    if (fPlotStandardClusterTHnSparse) {
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
      nbins[dim] = fNM02HistBins;
      binEdges[dim] = fM02HistBins;
      min[dim] = fM02HistBins[0];
      max[dim] = fM02HistBins[fNM02HistBins];
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
      
    }
    
    if (fPlotEvenOddEta) {
      title[dim] = "Even/odd eta";
      nbins[dim] = 2;
      min[dim] = -0.5;
      max[dim] = 1.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
    }
    
    if (fPlotNearestNeighborDistribution) {
      title[dim] = "#DeltaR_{NN}";
      nbins[dim] = 100;
      min[dim] = 0.;
      max[dim] = 1.;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
    }
    
    if (fPlotClusterCone) {

      title[dim] = "Cone type";
      nbins[dim] = 2;
      min[dim] = -0.5;
      max[dim] = 1.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "R";
      nbins[dim] = 2;
      min[dim] = 0;
      max[dim] = 0.15;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "#it{E}_{cone} (GeV)";
      nbins[dim] = fNPtHistBins;
      binEdges[dim] = fPtHistBins;
      min[dim] = fPtHistBins[0];
      max[dim] = fPtHistBins[fNPtHistBins];
      dim++;
    
    }
    
    if (fPlotCaloCentrality) {
      
      title[dim] = "#it{#rho}_{cell cone} (GeV)";
      nbins[dim] = 100;
      min[dim] = 0.;
      max[dim] = 1000.;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      title[dim] = "Ncells cone";
      nbins[dim] = 100;
      min[dim] = -0.5;
      max[dim] = 99.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
      
      if (fPlotCellSMDensity) {
        title[dim] = "#it{E}_{cell SM} (GeV)";
        nbins[dim] = fNPtHistBins;
        binEdges[dim] = fPtHistBins;
        min[dim] = fPtHistBins[0];
        max[dim] = fPtHistBins[fNPtHistBins];
        dim++;
        
        title[dim] = "Ncells SM";
        nbins[dim] = 100;
        min[dim] = -0.5;
        max[dim] = 999.5;
        binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
        dim++;
      }
      
    }
    
    TString thnname = TString::Format("%s/clusterObservables", cont->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(title[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
    
    // Plot Fcross distribution
    if (fPlotExotics) {
      histname = TString::Format("%s/hFcrossEMCal", cont->GetArrayName().Data());
      htitle = histname + ";Centrality (%);Fcross;#it{E}_{clus} (GeV/)";
      fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, nExBins, exBins, fNPtHistBins, fPtHistBins);
    }
    
  }
  
}

/*
 * This function allocates the histograms for calo cells.
 */
void AliAnalysisTaskEmcalVsPhos::AllocateCellHistograms()
{
  TString histname;
  TString htitle;
  
  Double_t* clusType = new Double_t[3+1];
  GenerateFixedBinArray(3, -0.5, 2.5, clusType);
  
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
 * This function allocates the histograms for neutral jets.
 */
void AliAnalysisTaskEmcalVsPhos::AllocateNeutralJetHistograms()
{
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
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
    
    axisTitle[dimJet] = "#it{E}_{T}, acc clus within R (GeV)";
    nbinsJet[dimJet] = fNPtHistBins;
    binEdgesJet[dimJet] = fPtHistBins;
    minJet[dimJet] = fPtHistBins[0];
    maxJet[dimJet] = fPtHistBins[fNPtHistBins];
    dimJet++;
    
    axisTitle[dimJet] = "#it{E}_{T}, acc cell within R (GeV)";
    nbinsJet[dimJet] = fNPtHistBins;
    binEdgesJet[dimJet] = fPtHistBins;
    minJet[dimJet] = fPtHistBins[0];
    maxJet[dimJet] = fPtHistBins[fNPtHistBins];
    dimJet++;
    
    TString thnname = TString::Format("%s/hClusterJetObservables", jets->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dimJet, nbinsJet, minJet, maxJet);
    for (Int_t i = 0; i < dimJet; i++) {
      hn->GetAxis(i)->SetTitle(axisTitle[i]);
      hn->SetBinEdges(i, binEdgesJet[i]);
    }
    
  }
}

/*
 * This function allocates the histograms for clusters within jets.
 */
void AliAnalysisTaskEmcalVsPhos::AllocateClustersInJetsHistograms()
{
  TString histname;
  TString htitle;
  Int_t nPtBins = TMath::CeilNint(fMaxPt/2);
  
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {

    // Plot cluster spectra within jets
    // (centrality, cluster energy, jet pT, jet eta, jet phi)
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
    nbins[dim] = 40;
    min[dim] = -0.5;
    max[dim] = 0.5;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    title[dim] = "#phi_{jet}";
    nbins[dim] = 200;
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

/*
 * This function allocates the histograms for the calorimeter performance study.
 */
void AliAnalysisTaskEmcalVsPhos::AllocateEventHistograms()
{
  TString histname;
  TString htitle;
  
  AliEmcalContainer* cont = 0;
  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextClusColl()))) {
    
    // Plot cluster THnSparse (centrality, EMCal cluster energy, PHOS cluster energy, track pT)
    Int_t dim = 0;
    TString title[20];
    Int_t nbins[20] = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    Double_t *binEdges[20] = {0};
    
    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
      title[dim] = "Centrality %";
      if (fPlotFineGrainedCentrality) {
        nbins[dim] = 18;
        min[dim] = 0;
        max[dim] = 90;
        binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      }
      else {
        nbins[dim] = fNCentHistBins;
        binEdges[dim] = fCentHistBins;
        min[dim] = fCentHistBins[0];
        max[dim] = fCentHistBins[fNCentHistBins];
      }
      dim++;
    }
    
    title[dim] = "#Sigma#it{E}_{clus,EMCal} (GeV)";
    nbins[dim] = 200;
    min[dim] = 0;
    max[dim] = 1500;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    title[dim] = "#Sigma#it{E}_{clus,PHOS} (GeV)";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 200;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    title[dim] = "#Sigma#it{p}_{T,tracks} (GeV)";
    nbins[dim] = 200;
    min[dim] = 0;
    max[dim] = 1500;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    TString thnname = "eventObservables";
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(title[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
    
  }
  
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalVsPhos::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();

  fNeedEmcalGeom = kTRUE;

  // Load the PHOS geometry
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

/**
 * This function (overloading the base class) uses AliEventCuts to perform event selection.
 */
Bool_t AliAnalysisTaskEmcalVsPhos::IsEventSelected()
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
Bool_t AliAnalysisTaskEmcalVsPhos::Run()
{
  return kTRUE;
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalVsPhos::FillHistograms()
{
  
  if (fPlotClusterHistograms) {
    FillClusterHistograms();
  }
  
  if (fPlotNeutralJets) {
    FillNeutralJetHistograms();
  }
  
  if (fPlotClustersInJets) {
    FillClustersInJetsHistograms();
  }
  
  if (fPlotCellHistograms) {
    FillCellHistograms();
  }
  
  if (fPlotEventHistograms) {
    FillEventHistograms();
  }
  
  return kTRUE;
}

/*
 * This function fills the histograms for the calorimeter performance study.
 */
void AliAnalysisTaskEmcalVsPhos::FillClusterHistograms()
{
  // Define some vars
  TString histname;
  Double_t Enonlin;
  Double_t Ehadcorr;
  Int_t absId;
  Double_t ecell;
  Double_t leadEcell;
  Int_t leadAbsId;
  
  // Get cells from event
  fCaloCells = InputEvent()->GetEMCALCells();
  AliVCaloCells* phosCaloCells = InputEvent()->GetPHOSCells();
  
  // Loop through clusters and plot cluster THnSparse (centrality, cluster type, E, E-hadcorr, has matched track, M02, Ncells)
  AliClusterContainer* clusters = 0;
  const AliVCluster* clus;
  TString clustersName;
  TIter nextClusColl(&fClusterCollArray);
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    AliClusterIterableMomentumContainer itcont = clusters->all_momentum();
    clustersName = clusters->GetArrayName();
    for (AliClusterIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
      
      clus = it->second;
    
      // Determine cluster type (EMCal/DCal/Phos)
      ClusterType clusType = kNA;
      if (clus->IsEMCAL()) {
        Double_t phi = it->first.Phi_0_2pi();
        Int_t isDcal = Int_t(phi > fgkEMCalDCalPhiDivide);
        if (isDcal == 0) {
          clusType = kEMCal;
        } else if (isDcal == 1) {
          clusType = kDCal;
        }
      } else if (clus->GetType() == AliVCluster::kPHOSNeutral){
        clusType = kPHOS;
      }
      
      // rejection reason plots, to make efficiency correction
      if (clus->IsEMCAL()) {
        histname = TString::Format("%s/hClusterRejectionReasonEMCal", clusters->GetArrayName().Data());
        UInt_t rejectionReason = 0;
        if (!clusters->AcceptCluster(it.current_index(), rejectionReason)) {
          fHistManager.FillTH2(histname, clusters->GetRejectionReasonBitPosition(rejectionReason), it->first.E());
          continue;
        }
      } else if (clus->GetType() == AliVCluster::kPHOSNeutral){
        histname = TString::Format("%s/hClusterRejectionReasonPHOS", clusters->GetArrayName().Data());
        UInt_t rejectionReason = 0;
        if (!clusters->AcceptCluster(it.current_index(), rejectionReason)) {
          fHistManager.FillTH2(histname, clusters->GetRejectionReasonBitPosition(rejectionReason), it->first.E());
          continue;
        }
      } else {
        continue; // avoid CPV clusters
      }
      
      // Fill cluster spectra by SM, and fill cell histograms
      Enonlin = 0;
      Ehadcorr = 0;
      if (clus->IsEMCAL()) {
        
        Ehadcorr = clus->GetHadCorrEnergy();
        Enonlin = clus->GetNonLinCorrEnergy();
        if (fPlotClusWithoutNonLinCorr) {
          Enonlin = clus->E();
        }
        
        if (fPlotExotics) {
          histname = TString::Format("%s/hFcrossEMCal", clusters->GetArrayName().Data());
          Double_t Fcross = GetFcross(clus, fCaloCells);
          fHistManager.FillTH3(histname, fCent, Fcross, Enonlin);
        }
        
        Int_t sm = fGeom->GetSuperModuleNumber(clus->GetCellAbsId(0));
        if (sm >=0 && sm < 20) {
          histname = TString::Format("%s/BySM/hEmcalClusEnergy_SM%d", clusters->GetArrayName().Data(), sm);
          fHistManager.FillTH1(histname, Enonlin);
        }
        else {
          AliError(Form("Supermodule %d does not exist!", sm));
        }
        
        // Get cells from each accepted cluster, and plot centrality vs. cell energy vs. cell type
        histname = TString::Format("Cells/hCellEnergyAccepted");
        leadEcell = 0;
        for (Int_t iCell = 0; iCell < clus->GetNCells(); iCell++){
          absId = clus->GetCellAbsId(iCell);
          ecell = fCaloCells->GetCellAmplitude(absId);
          fHistManager.FillTH3(histname, ecell, fCent, kEMCal); // Note: I don't distinguish EMCal from DCal cells
          if (ecell > leadEcell) {
            leadEcell = ecell;
            leadAbsId = absId;
          }
        }
        // Plot also the leading cell
        histname = TString::Format("Cells/hCellEnergyLeading");
        fHistManager.FillTH3(histname, leadEcell, fCent, kEMCal);
        
      } else if (clus->GetType() == AliVCluster::kPHOSNeutral){
        
        Ehadcorr = clus->GetCoreEnergy();
        Enonlin = clus->E();
        
        Int_t relid[4];
        if (fPHOSGeo) {
          fPHOSGeo->AbsToRelNumbering(clus->GetCellAbsId(0), relid);
          Int_t sm = relid[0];
          if (sm >=1 && sm < 5) {
            histname = TString::Format("%s/BySM/hPhosClusEnergy_SM%d", clusters->GetArrayName().Data(), sm);
            fHistManager.FillTH1(histname, clus->E());
          }
          else {
            AliError(Form("Supermodule %d does not exist!", sm));
          }
        }
        
        // Get cells from each accepted cluster, and plot centrality vs. cell energy vs. cell type
        histname = TString::Format("Cells/hCellEnergyAccepted");
        leadEcell = 0;
        for (Int_t iCell = 0; iCell < clus->GetNCells(); iCell++){
          absId = clus->GetCellAbsId(iCell);
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
      Int_t nMatchedTracks = clus->GetNTracksMatched();
      if (nMatchedTracks == 0) {
        hasMatchedTrack = 0;
      } else if (nMatchedTracks > 0) {
        hasMatchedTrack = 1;
      }
      
      // Check if the cluster passes the dispersion cut for photon-like cluster (meaningful only for PHOS)
      Int_t passedDispersionCut = 0;
      if (clus->Chi2() < 2.5*2.5) {
        passedDispersionCut = 1;
      }
      
      // Fill info about the cluster
      Double_t eta = it->first.Eta();
      Double_t phi = it->first.Phi_0_2pi();
      Double_t M02 = clus->GetM02();
      Int_t nCells = clus->GetNCells();
      Double_t distNN = FindNearestNeighborDistance(it->first);
      
      // If cluster is EMCal, find whether the eta column is even or odd
      Int_t isOddEta = -1;
      Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;
      if (clus->IsEMCAL()) {
        fGeom->GetCellIndex(leadAbsId, nSupMod, nModule, nIphi, nIeta);
        fGeom->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphi, ieta);
        isOddEta = ieta % 2;
      }

      // Standard option: fill once per cluster
      if (!fPlotClusterCone && !fPlotCaloCentrality) {
          FillClusterTHnSparse(clustersName, eta, phi, Enonlin, Ehadcorr, hasMatchedTrack, M02, nCells, passedDispersionCut, distNN, isOddEta);
      }
      
      if (fPlotCaloCentrality) {
        
        Double_t eCellCone = 0.;
        Int_t nCellsCone = 0;
        Double_t eCellSM = 0.;
        Int_t nCellsSM = 0;
        Double_t areaCone = TMath::Pi() * 0.07 * 0.07;
        Double_t areaCell;
        Double_t eDensityCone;
        
        // Get the SM number
        Int_t sm = -1;
        if (clusType == kEMCal) {
          sm = fGeom->GetSuperModuleNumber(clus->GetCellAbsId(0));
          areaCell = 0.014*0.014;
        }
        if (clusType == kPHOS) {
          Int_t relid[4];
          fPHOSGeo->AbsToRelNumbering(clus->GetCellAbsId(0), relid);
          sm = relid[0];
          areaCell = 0.014*0.014*(2.2/6.0)*(2.2/6.0); // approximating same r
        }
        
        // Only fill the THnSparse if the cluster is located in a full SM of EMCal or PHOS
        if ( (clusType == kEMCal && sm < 10 ) || (clusType == kPHOS && sm < 4) ) {
          
          eCellCone = GetConeCellEnergy(eta, phi, 0.07) - Enonlin;
          nCellsCone = (Int_t)GetConeCellEnergy(eta, phi, 0.07, kTRUE) - nCells;

          eDensityCone = eCellCone / (areaCone - nCells*areaCell);
       
          if (fPlotCellSMDensity) {
            eCellSM = GetSMCellEnergy(sm, clusType) - Enonlin;
            nCellsSM = (Int_t)GetSMCellEnergy(sm, clusType, kTRUE) - nCells;
          }
        
          FillClusterTHnSparse(clustersName, eta, phi, Enonlin, eDensityCone, eCellSM, nCellsCone, nCellsSM);
          
        }
        
      }
      
      // If cluster cone option enabled, fill for each R and cone type
      if (fPlotClusterCone) {

        // cluster cone, R=0.05
        FillClusterTHnSparse(clustersName, eta, phi, Enonlin, Ehadcorr, hasMatchedTrack, M02, nCells, passedDispersionCut, distNN, isOddEta, 0, 0.05, GetConeClusterEnergy(eta, phi, 0.05));
        // cluster cone, R=0.1
        FillClusterTHnSparse(clustersName, eta, phi, Enonlin, Ehadcorr, hasMatchedTrack, M02, nCells, passedDispersionCut, distNN, isOddEta, 0, 0.1, GetConeClusterEnergy(eta, phi, 0.1));
        // cell cone, R=0.05
        FillClusterTHnSparse(clustersName, eta, phi, Enonlin, Ehadcorr, hasMatchedTrack, M02, nCells, passedDispersionCut, distNN, isOddEta, 1, 0.05, GetConeCellEnergy(eta, phi, 0.05));
        // cell cone, R=0.1
        FillClusterTHnSparse(clustersName, eta, phi, Enonlin, Ehadcorr, hasMatchedTrack, M02, nCells, passedDispersionCut, distNN, isOddEta, 1, 0.1, GetConeCellEnergy(eta, phi, 0.1));
        
      }

    }
  }
}

/*
 * This function fills the cluster THnSparse.
 */
void AliAnalysisTaskEmcalVsPhos::FillClusterTHnSparse(TString clustersName, Double_t eta, Double_t phi, Double_t Enonlin, Double_t Ehadcorr, Int_t hasMatchedTrack, Double_t M02, Int_t nCells, Int_t passedDispersionCut, Double_t distNN, Int_t isOddEta, Int_t coneType, Double_t R, Double_t Econe)
{
  Double_t contents[30]={0};
  TString histname = TString::Format("%s/clusterObservables", clustersName.Data());
  THnSparse* histClusterObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));
  if (!histClusterObservables) return;
  for (Int_t i = 0; i < histClusterObservables->GetNdimensions(); i++) {
    TString title(histClusterObservables->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = fCent;
    else if (title=="#eta")
      contents[i] = eta;
    else if (title=="#phi")
      contents[i] = phi;
    else if (title=="#it{E}_{clus} (GeV)")
      contents[i] = Enonlin;
    else if (title=="#it{E}_{clus, hadcorr} or #it{E}_{core} (GeV)")
      contents[i] = Ehadcorr;
    else if (title=="Matched track")
      contents[i] = hasMatchedTrack;
    else if (title=="M02")
      contents[i] = M02;
    else if (title=="Ncells")
      contents[i] = nCells;
    else if (title=="Dispersion cut")
      contents[i] = passedDispersionCut;
    else if (title=="#DeltaR_{NN}")
      contents[i] = distNN;
    else if (title=="Even/odd eta")
      contents[i] = isOddEta;
    else if (title=="Cone type")
      contents[i] = coneType;
    else if (title=="R")
      contents[i] = R;
    else if (title=="#it{E}_{cone} (GeV)")
      contents[i] = Econe;
    else
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
                              }
  histClusterObservables->Fill(contents);

}

/*
 * This function fills the cluster THnSparse (alternate signature, used for local density option).
 */
void AliAnalysisTaskEmcalVsPhos::FillClusterTHnSparse(TString clustersName, Double_t eta, Double_t phi, Double_t Enonlin, Double_t eCellCone, Double_t eCellSM, Int_t nCellsCone, Int_t nCellsSM)
{
  Double_t contents[30]={0};
  TString histname = TString::Format("%s/clusterObservables", clustersName.Data());
  THnSparse* histClusterObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));
  if (!histClusterObservables) return;
  for (Int_t i = 0; i < histClusterObservables->GetNdimensions(); i++) {
    TString title(histClusterObservables->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = fCent;
    else if (title=="#eta")
      contents[i] = eta;
    else if (title=="#phi")
      contents[i] = phi;
    else if (title=="#it{E}_{clus} (GeV)")
      contents[i] = Enonlin;
    else if (title=="#it{#rho}_{cell cone} (GeV)")
      contents[i] = eCellCone;
    else if (title=="Ncells cone")
      contents[i] = nCellsCone;
    else if (title=="#it{E}_{cell SM} (GeV)")
      contents[i] = eCellSM;
    else if (title=="Ncells SM")
      contents[i] = nCellsSM;
    else
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }
  histClusterObservables->Fill(contents);
  
}

/*
 * This function fills cell histograms.
 */
void AliAnalysisTaskEmcalVsPhos::FillCellHistograms()
{
  // Define some vars
  TString histname;
  Int_t absId;
  Double_t ecell;
  
  // Get cells from event
  fCaloCells = InputEvent()->GetEMCALCells();
  AliVCaloCells* phosCaloCells = InputEvent()->GetPHOSCells();
  
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
    
    if (absId < 0) {
      continue; // skip CPV cells
    }
    
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

/*
 * This function fills neutral jet histograms.
 */
void AliAnalysisTaskEmcalVsPhos::FillNeutralJetHistograms()
{
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
  
    // plot neutral jets THnSparse (centrality, eta, phi, E, Nclusters)
    Double_t contents[30]={0};
    TString histname = TString::Format("%s/hClusterJetObservables", jets->GetArrayName().Data());
    THnSparse* histJetObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));
    if (!histJetObservables) return;
    
    for (auto jet : jets->accepted()) {
    
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
        else if (title=="#it{E}_{T}, acc clus within R (GeV)")
          contents[i] = GetConeClusterEnergy(jet->Eta(), jet->Phi_0_2pi(), jets->GetJetRadius());
        else if (title=="#it{E}_{T}, acc cell within R (GeV)")
          contents[i] = GetConeCellEnergy(jet->Eta(), jet->Phi_0_2pi(), jets->GetJetRadius());
        else
          AliWarning(Form("Unable to fill dimension %s!",title.Data()));
      }
      histJetObservables->Fill(contents);
      
    }
  }
}

/*
 * This function fills clusters within jets and estimates the JES shift due to the bump.
 */
void AliAnalysisTaskEmcalVsPhos::FillClustersInJetsHistograms()
{
  TString histname;
  Double_t rho;
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    rho = jets->GetRhoVal();
    
    for (const auto jet : jets->accepted()) {

      // Fill cluster spectra of clusters within jets
      //(centrality, cluster energy, jet pT, jet eta, jet phi)
      histname = TString::Format("%s/hClustersInJets", jets->GetArrayName().Data());
      Int_t nClusters = jet->GetNumberOfClusters();
      AliVCluster* clus;
      for (Int_t iClus = 0; iClus < nClusters; iClus++) {
        
        clus = jet->Cluster(iClus);
        
        if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
          Double_t x[5] = {fCent, clus->GetNonLinCorrEnergy(), GetJetPt(jet, rho), jet->Eta(), jet->Phi_0_2pi()};
          fHistManager.FillTHnSparse(histname, x);
        }
        else {
          Double_t x[4] = {clus->GetNonLinCorrEnergy(), GetJetPt(jet, rho), jet->Eta(), jet->Phi_0_2pi()};
          fHistManager.FillTHnSparse(histname, x);
        }

      }
        
      // Loop through clusters, and plot estimated shift in JES due to cluster bump
      // Only do for 0-10% centrality, and for EMCal/DCal
      Double_t eclus;
      Double_t shift;
      Double_t shiftSum = 0;
      if (fCent < 10. || fForceBeamType == AliAnalysisTaskEmcal::kpp) {
        if (GetJetType(jet) > -0.5 && GetJetType(jet) < 1.5) {
          for (Int_t iClus = 0; iClus < nClusters; iClus++) {
            clus = jet->Cluster(iClus);
            eclus = clus->GetNonLinCorrEnergy();
            if (eclus > 0.5) {
              shift = 0.79 * TMath::Exp(-0.5 * ((eclus - 3.81) / 1.50)*((eclus - 3.81) / 1.50) );
              shiftSum += shift;
            }
          }
          histname = TString::Format("%s/hCaloJESshift", jets->GetArrayName().Data());
          fHistManager.FillTH3(histname, GetJetType(jet), GetJetPt(jet, rho), shiftSum);
        }
      }
    
    }
  }
}


/*
 * This function fills the histograms for the event thnsparse (centrality, EMCal clus energy, PHOS clus energy, track pT)
 */
void AliAnalysisTaskEmcalVsPhos::FillEventHistograms()
{
  TString histname;
  Double_t sumEMCal = 0;
  Double_t sumPHOS = 0;
  Double_t sumTracks = 0;
  
  // Loop through clusters
  const AliVCluster* clus;
  AliClusterContainer* clusters = GetClusterContainer(0);
  for (auto clusIterator : clusters->accepted_momentum() ) {
    
    clus = clusIterator.second;
  
    if (clus->IsEMCAL()) {
      sumEMCal += clus->GetNonLinCorrEnergy();
    } else if (clus->GetType() == AliVCluster::kPHOSNeutral){
      sumPHOS += clus->E();
    }
    
  }
  
  // Loop through tracks
  AliTrackContainer* trackCont = dynamic_cast<AliTrackContainer*>(GetParticleContainer("tracks"));
  const AliVTrack* track;
  for (auto trackIterator : trackCont->accepted_momentum() ) {
    
    track = trackIterator.second;
    
    // Sum if track has matched EMCal cluster (within R=0.1, by default)
    if (track->IsEMCAL()) {
      sumTracks += track->Pt();
    }
    
  }
  
  // Fill the THnSparse
  Double_t contents[30]={0};
  THnSparse* histEventObservables = static_cast<THnSparse*>(fHistManager.FindObject("eventObservables"));
  if (!histEventObservables) return;
  for (Int_t i = 0; i < histEventObservables->GetNdimensions(); i++) {
    TString title(histEventObservables->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
    contents[i] = fCent;
    else if (title=="#Sigma#it{E}_{clus,EMCal} (GeV)")
    contents[i] = sumEMCal;
    else if (title=="#Sigma#it{E}_{clus,PHOS} (GeV)")
    contents[i] = sumPHOS;
    else if (title=="#Sigma#it{p}_{T,tracks} (GeV)")
    contents[i] = sumTracks;
    else
    AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }
  histEventObservables->Fill(contents);
  
}

/**
 * Get pT of jet -- background subtracted
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetJetPt(const AliEmcalJet* jet, Double_t rho)
{
  Double_t pT = jet->Pt() - rho * jet->Area();
  return pT;
}

/**
 * Get deltaR of a track/cluster and a reference point.
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetDeltaR(AliTLorentzVector part, Double_t etaRef, Double_t phiRef)
{
  Double_t deltaPhi = TMath::Abs(part.Phi_0_2pi() - phiRef);
  Double_t deltaEta = TMath::Abs(part.Eta() - etaRef);
  Double_t deltaR = TMath::Sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta );
  return deltaR;
}

/**
 * Get deltaR of two points.
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetDeltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  Double_t deltaPhi = TMath::Abs(phi1-phi2);
  Double_t deltaEta = TMath::Abs(eta1-eta2);
  Double_t deltaR = TMath::Sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta );
  return deltaR;
}

/**
 * Get calo acceptance type of jet
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetJetType(const AliEmcalJet* jet)
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
 * Compute Fcross of a cluster
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetFcross(const AliVCluster *cluster, AliVCaloCells *cells)
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

/**
 * Compute the distance to the nearest accepted cluster
 */
Double_t AliAnalysisTaskEmcalVsPhos::FindNearestNeighborDistance(AliTLorentzVector clusterRef)
{
  Double_t distNN = 10.;
  Double_t etaRef = clusterRef.Eta();
  Double_t phiRef = clusterRef.Phi_0_2pi();
  
  AliClusterContainer* clusters = GetClusterContainer(0);
  AliTLorentzVector clusNNcand;
  for (auto clusIterator : clusters->accepted_momentum() ) {
    
    clusNNcand.Clear();
    clusNNcand = clusIterator.first;
    
    Double_t distNNcand = GetDeltaR(clusNNcand, etaRef, phiRef);
    
    if (distNNcand < distNN && distNNcand > 0.001) {
      distNN = distNNcand;
    }
    
  }
  
  return distNN;
  
}

/**
 * Compute the cluster energy within a cone of radius R centered at etaRef,phiRef.
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetConeClusterEnergy(Double_t etaRef, Double_t phiRef, Double_t R)
{
  AliClusterContainer* clusCont = GetClusterContainer(0);
  AliTLorentzVector clus;
  Double_t energy = 0.;
  for (auto clusIterator : clusCont->accepted_momentum() ) {
    
    clus.Clear();
    clus = clusIterator.first;
    
    if (GetDeltaR(clus, etaRef, phiRef) < R) {
      energy += clus.E();
    }
  }
  return energy;
}

/**
 * Compute the cell energy within a cone centered at the jet axis, excluding cells from rejected clusters if requested.
 * Optionally, can instead return the number of cells.
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetConeCellEnergy(Double_t etaRef, Double_t phiRef, Double_t R, Bool_t returnNcells)
{
  Double_t energy = 0.;
  Double_t nCells = 0.;
  
  // Get cells from event
  fCaloCells = InputEvent()->GetEMCALCells();
  AliVCaloCells* phosCaloCells = InputEvent()->GetPHOSCells();
  
  Double_t eta;
  Double_t phi;
  Int_t absId;
  TVector3 pos;
  Int_t relid[4];
  Int_t sm;
  
  for (Int_t i=0; i<fCaloCells->GetNumberOfCells(); i++) {
    
    absId = fCaloCells->GetCellNumber(i);
    
    fGeom->EtaPhiFromIndex(absId, eta, phi);
    phi = TVector2::Phi_0_2pi(phi);

    if (GetDeltaR(eta, phi, etaRef, phiRef) < R) {

      if (fExcludeRejectedCells) {
        if (IsCellRejected(absId, kEMCal)) {
          continue;
        }
      }

      energy += fCaloCells->GetCellAmplitude(absId);
      nCells += 1.;
    }
    
  }
  
  for (Int_t i=0; i<phosCaloCells->GetNumberOfCells(); i++) {
    
    absId = phosCaloCells->GetCellNumber(i);

    fPHOSGeo->AbsToRelNumbering(absId, relid);
    sm = relid[0];
    if (sm < 1 || sm > 4) {
      continue;
    }

    fPHOSGeo->RelPosInAlice(absId, pos); // pos then contains (x,y,z) coordinate of cell
    eta = pos.Eta();
    phi = pos.Phi();
    phi = TVector2::Phi_0_2pi(phi);
    
    if (GetDeltaR(eta, phi, etaRef, phiRef) < R) {

      if (fExcludeRejectedCells) {
        if (IsCellRejected(absId, kPHOS)) {
          continue;
        }
      }

      energy += phosCaloCells->GetCellAmplitude(absId);
      nCells += 1.;
    }
    
  }
  
  if (returnNcells) {
    return nCells;
  }
  
  return energy;
}

/**
 * Compute the cell energy within a SM, excluding cells from rejected clusters if requested.
 * Optionally, can instead return the number of cells.
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetSMCellEnergy(Int_t sm, Int_t clusType, Bool_t returnNcells)
{
  Double_t energy = 0.;
  Double_t nCells = 0.;
  Int_t absId;
  Int_t cellSM;
  Int_t relid[4];
  
  if (clusType == kEMCal) {
    
    fCaloCells = InputEvent()->GetEMCALCells();
  
    for (Int_t i=0; i<fCaloCells->GetNumberOfCells(); i++) {
      
      absId = fCaloCells->GetCellNumber(i);
      cellSM = fGeom->GetSuperModuleNumber(absId);
      
      if (cellSM == sm) {
        
        if (fExcludeRejectedCells) {
          if (IsCellRejected(absId, kEMCal)) {
            continue;
          }
        }
        
        energy += fCaloCells->GetCellAmplitude(absId);
        nCells += 1.;
      
      }
    }
  }
  
  if (clusType == kPHOS) {
    
    AliVCaloCells* phosCaloCells = InputEvent()->GetPHOSCells();
  
    for (Int_t i=0; i<phosCaloCells->GetNumberOfCells(); i++) {
      
      absId = phosCaloCells->GetCellNumber(i);

      fPHOSGeo->AbsToRelNumbering(absId, relid);
      cellSM = relid[0];
    
      if (cellSM == sm) {
        
        if (fExcludeRejectedCells) {
          if (IsCellRejected(absId, kPHOS)) {
            continue;
          }
        }
        
        energy += phosCaloCells->GetCellAmplitude(absId);
        nCells += 1.;
        
      }
    }
  }
  
  if (returnNcells) {
    return nCells;
  }
  
  return energy;
}

/**
 * Determine whether a cell belongs to a rejected cluster
 */
Bool_t AliAnalysisTaskEmcalVsPhos::IsCellRejected(Int_t absId, Int_t cellType)
{
  AliClusterContainer* clusCont = GetClusterContainer(0);
  AliVCluster* cluster;
  
  UInt_t rejectionReason;
  Bool_t skipCell = kFALSE;
  Int_t clusType;
  
  AliClusterIterableMomentumContainer itcont = clusCont->all_momentum();
  for (AliClusterIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
  
    if (!clusCont->AcceptCluster(it.current_index(), rejectionReason)) {
      
      cluster = it->second;
      
      // check that the cell type and cluster type are the same
      if (cluster->IsEMCAL()) {
        clusType = kEMCal;
      }
      if (cluster->GetType() == AliVCluster::kPHOSNeutral) {
        clusType = kPHOS;
      }
      if (clusType != cellType) {
        continue;
      }
      
      // skip the cell if it belongs to a rejected cluster
      for (Int_t i = 0; i < cluster->GetNCells(); i++) {
        
        if (absId == cluster->GetCellAbsId(i)) {
          skipCell = kTRUE;
        }
        
      }
    }
  }
  return skipCell;
}

