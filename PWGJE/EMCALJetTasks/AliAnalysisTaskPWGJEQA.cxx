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

#include <cstring>

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TString.h>
#include <Riostream.h>
#include <TChain.h>
#include <TFile.h>

#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliTrackContainer.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCParticleContainer.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliPHOSGeometry.h"
#include "AliOADBContainer.h"

#include "AliAnalysisTaskPWGJEQA.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskPWGJEQA);
/// \endcond

//________________________________________________________________________
/// Default constructor for ROOT I/O purposes
AliAnalysisTaskPWGJEQA::AliAnalysisTaskPWGJEQA() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskPWGJEQA", kTRUE),
  fCellEnergyCut(0.05),
  fPtBinWidth(0.5),
  fMaxPt(250),
  fSeparateEMCalDCal(kTRUE),
  fNTotTracks(0),
  fLeadingTrack(),
  fGeneratorLevel(0),
  fGeneratorLevelName(),
  fDetectorLevel(0),
  fDetectorLevelName(),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNEtaHistBins(0),
  fEtaHistBins(0),
  fNPhiHistBins(0),
  fPhiHistBins(0),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtRelDiffHistBins(0),
  fPtRelDiffHistBins(0),
  fNPtResHistBins(0),
  fPtResHistBins(0),
  f1OverPtResHistBins(0),
  fN1OverPtResHistBins(0),
  fNIntegerHistBins(0),
  fIntegerHistBins(0),
  fTracks(0),
  fParticlesPhysPrim(0),
  fParticlesMatched(0),
  fHistManager("AliAnalysisTaskPWGJEQA"),
  fDoTrackQA(kTRUE),
  fDoEmcalQA(kTRUE),
  fDoJetQA(kTRUE),
  fDoEventQA(kTRUE),
  fRejectOutlierEvents(kFALSE),
  fIsPtHard(kFALSE),
  fPHOSGeo(nullptr)
{
  // Default constructor.
  memset(fNTotClusters, 0, sizeof(Int_t)*3);
  SetMakeGeneralHistograms(kTRUE);
  GenerateHistoBins();
}

//________________________________________________________________________
/// Standard named constructor
///
/// \param name Name of the task
AliAnalysisTaskPWGJEQA::AliAnalysisTaskPWGJEQA(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fCellEnergyCut(0.05),
  fPtBinWidth(0.5),
  fMaxPt(250),
  fSeparateEMCalDCal(kTRUE),
  fNTotTracks(0),
  fLeadingTrack(),
  fGeneratorLevel(0),
  fGeneratorLevelName(),
  fDetectorLevel(0),
  fDetectorLevelName(),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNEtaHistBins(0),
  fEtaHistBins(0),
  fNPhiHistBins(0),
  fPhiHistBins(0),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtRelDiffHistBins(0),
  fPtRelDiffHistBins(0),
  fNPtResHistBins(0),
  fPtResHistBins(0),
  f1OverPtResHistBins(0),
  fN1OverPtResHistBins(0),
  fNIntegerHistBins(0),
  fIntegerHistBins(0),
  fTracks(0),
  fParticlesPhysPrim(0),
  fParticlesMatched(0),
  fHistManager(name),
  fDoTrackQA(kTRUE),
  fDoEmcalQA(kTRUE),
  fDoJetQA(kTRUE),
  fDoEventQA(kTRUE),
  fRejectOutlierEvents(kFALSE),
  fIsPtHard(kFALSE),
  fPHOSGeo(nullptr)
{
  // Standard
  memset(fNTotClusters, 0, sizeof(Int_t)*3);
  SetMakeGeneralHistograms(kTRUE);
  GenerateHistoBins();
}

//________________________________________________________________________
AliAnalysisTaskPWGJEQA::~AliAnalysisTaskPWGJEQA()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::UserCreateOutputObjects()
{
  // Create histograms
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  
  // Set track container pointers
  fDetectorLevel = GetTrackContainer(fDetectorLevelName);
  fGeneratorLevel = GetMCParticleContainer(fGeneratorLevelName);
  
  // Allocate histograms for tracks, cells, clusters, jets
  if (fDoTrackQA) AllocateTrackHistograms();
  if (fDoEmcalQA) AllocateCellHistograms();
  if (fDoEmcalQA) AllocateClusterHistograms();
  if (fDoJetQA) AllocateJetHistograms();
  if (fDoEventQA) AllocateEventQAHistograms();
  
  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateTrackHistograms() {
  
  AllocateDetectorLevelTHnSparse();
  if (fGeneratorLevel) {
    AllocateGeneratorLevelTHnSparse();
    AllocateMatchedParticlesTHnSparse();
  }
}
  
//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateCellHistograms() {
  
  TString histname;
  TString title;
  
  if (!fCaloCellsName.IsNull()) {

    histname = TString::Format("%s/fHistCellEnergy", fCaloCellsName.Data());
    title = histname + ";#it{E}_{cell} (GeV);counts";
    fHistManager.CreateTH1(histname.Data(), title.Data(), 500,0, 100);
    
    histname = TString::Format("%s/fProfCellAbsIdEnergy", fCaloCellsName.Data());
    title = histname + ";cell absId;<#it{E}_{cell}> (GeV)";
    fHistManager.CreateTProfile(histname.Data(), title.Data(), 18000,0,18000);
     
    histname = TString::Format("%s/fHistCellTime", fCaloCellsName.Data());
    title = histname + ";#it{t}_{cell} (s);counts";
    fHistManager.CreateTH1(histname.Data(), title.Data(), 500,-5e-6, 5e-6);
    
    histname = TString::Format("%s/fProfCellAbsIdTime", fCaloCellsName.Data());
    title = histname + ";cell absId;<#it{t}_{cell}> (s)";
    fHistManager.CreateTProfile(histname.Data(), title.Data(), 18000,0,18000);
    
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateClusterHistograms() {
  
  AliEmcalContainer* cont = 0;
  TString histname;
  TString title;
  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);
  
  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextClusColl()))) {

    histname = TString::Format("%s/fHistClusterRejectionReason", cont->GetArrayName().Data());
    title = histname + ";Rejection reason;#it{E}_{clus} (GeV/);counts";
    TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 50, 0, 250);
    SetRejectionReasonLabels(hist->GetXaxis());
    
    const Int_t nEmcalSM = 20;
    for (Int_t sm = 0; sm < nEmcalSM; sm++) {
      histname = TString::Format("%s/BySM/hEmcalClusEnergy_SM%d", cont->GetArrayName().Data(), sm);
      title = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);
    }
    
    for (Int_t sm = 0; sm < 6; sm++) {
      histname = TString::Format("%s/BySM/hPhosClusEnergy_SM%d", cont->GetArrayName().Data(), sm);
      title = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);
    }
   
    // Cluster THnSparse
    Int_t dim = 0;
    TString title[20];
    Int_t nbins[20] = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    
    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
      title[dim] = "Centrality %";
      nbins[dim] = 5;
      min[dim] = 0;
      max[dim]= 100;
      dim++;
    }
    
    title[dim] = "#it{E}_{clus} (GeV)";
    nbins[dim] = nPtBins/2;
    min[dim] = 0;
    max[dim] = fMaxPt;
    dim++;
    
    title[dim] = "#eta";
    nbins[dim] = fNEtaHistBins;
    min[dim] = -1;
    max[dim] = 1;
    dim++;
    
    title[dim] = "#phi";
    nbins[dim] = fNPhiHistBins;
    min[dim] = 0;
    max[dim] = 2*TMath::Pi();
    dim++;
    
    title[dim] = "cluster type";
    nbins[dim] = 3;
    min[dim] = -0.5;
    max[dim] = 2.5;
    dim++;
    
    TString thnname = TString::Format("%s/clusterObservables", cont->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(title[i]);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateJetHistograms() {
  
  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);
  
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    // Allocate THnSparse
    Double_t jetRadius = jets->GetJetRadius();
    
    TString axisTitle[30]= {""};
    Int_t nbins[30]  = {0};
    Double_t min[30] = {0.};
    Double_t max[30] = {0.};
    Int_t dim = 0;
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "Centrality (%)";
      nbins[dim] = 5;
      min[dim] = 0;
      max[dim] = 100;
      dim++;
    }
    
    axisTitle[dim] = "#eta_{jet}";
    nbins[dim] = nPtBins/10;
    min[dim] = -1;
    max[dim] = 1;
    dim++;
    
    axisTitle[dim] = "#phi_{jet} (rad)";
    nbins[dim] = nPtBins/10*3;
    min[dim] = 0;
    max[dim] = 2*TMath::Pi();
    dim++;
    
    axisTitle[dim] = "#it{p}_{T} (GeV/#it{c})";
    nbins[dim] = nPtBins/2;
    min[dim] = 0;
    max[dim] = fMaxPt;
    dim++;
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "#it{p}_{T}^{corr} (GeV/#it{c})";
      nbins[dim] = nPtBins/2;
      min[dim] = -fMaxPt/2;
      max[dim] = fMaxPt/2;
      dim++;
    }
    
    if (fClusterCollArray.GetEntriesFast() > 0 && fParticleCollArray.GetEntriesFast() > 0) {
      axisTitle[dim] = "NEF";
      nbins[dim] = nPtBins/20;
      min[dim] = 0;
      max[dim] = 1.0;
      dim++;
    }
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "No. of constituents";
      nbins[dim] = 125;
      min[dim] = 0;
      max[dim] = 250;
      dim++;
    }
    else {
      axisTitle[dim] = "No. of constituents";
      nbins[dim] = 50;
      min[dim] = -0.5;
      max[dim] = 49.5;
      dim++;
    }
    
    axisTitle[dim] = "#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    nbins[dim] = nPtBins/10*3;
    min[dim] = 0;
    max[dim] = 150;
    dim++;
    
    TString thnname = TString::Format("%s/fHistJetObservables", jets->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(axisTitle[i]);
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

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateEventQAHistograms() {

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);
  
  Int_t dim = 0;
  TString axistitle[40];
  Int_t nbins[40] = {0};
  Double_t min[40] = {0};
  Double_t max[40] = {0};
  
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    axistitle[dim] = "Centrality %";
    nbins[dim] = 5;
    min[dim] = 0;
    max[dim] = 100;
    dim++;
  }
  
  if (fParticleCollArray.GetEntriesFast()>0) {
    axistitle[dim] = "No. of tracks";
    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
      nbins[dim] = 100;
      min[dim] = -0.5;
      max[dim] = 6000-0.5;
    }
    else {
      nbins[dim] = 100;
      min[dim] = 0;
      max[dim] = 200;
    }
    dim++;
    
    axistitle[dim] = "#it{p}_{T,track}^{leading} (GeV/c)";
    nbins[dim] = nPtBins/2;
    min[dim] = 0;
    max[dim] = fMaxPt;
    dim++;
  }
  
  if (fClusterCollArray.GetEntriesFast()>0) {
    axistitle[dim] = "No. of clusters";
    
    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
      nbins[dim] = 100;
      min[dim] = -0.5;
      max[dim] = 4000-0.5;
    }
    else {
      nbins[dim] = 100;
      min[dim] = 0;
      max[dim] = 200;
    }
    dim++;
    
    if (fSeparateEMCalDCal) {
      axistitle[dim] = "#it{E}_{EMCal cluster}^{leading} (GeV)";
      nbins[dim] = nPtBins/2;
      min[dim] = 0;
      max[dim] = fMaxPt;
      dim++;
      
      axistitle[dim] = "#it{E}_{DCal cluster}^{leading} (GeV)";
      nbins[dim] = nPtBins/2;
      min[dim] = 0;
      max[dim] = fMaxPt;
      dim++;
      
      axistitle[dim] = "#it{E}_{PHOS cluster}^{leading} (GeV)";
      nbins[dim] = nPtBins/2;
      min[dim] = 0;
      max[dim] = fMaxPt;
      dim++;
    }
    else {
      axistitle[dim] = "#it{E}_{cluster}^{leading} (GeV)";
      nbins[dim] = nPtBins/2;
      min[dim] = 0;
      max[dim] = fMaxPt;
      dim++;
    }
  }
  
  THnSparse* hn = fHistManager.CreateTHnSparse("eventQA","eventQA",dim,nbins,min,max);
  for (Int_t i = 0; i < dim; i++)
    hn->GetAxis(i)->SetTitle(axistitle[i]);
  
  if (fIsPtHard) {
    TString histname = "hPtHard";
    TString title = histname + ";#it{p}_{T,hard} (GeV/c);counts";
    fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);
    
    TH1* hTrials = fHistManager.CreateTH1("hNtrials", "hNtrials", 1, 0, 1);
    hTrials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
    
    TH1* hxsection = fHistManager.CreateTH1("hXsec", "hXsec", 1, 0, 1);
    hxsection->GetXaxis()->SetBinLabel(1,"<#sigma>");
  }
  
  fHistEventRejection->GetXaxis()->SetBinLabel(15,"PtHardBinJetOutlier");

}

//________________________________________________________________________
/// Generate histogram binning arrays for track histograms
void AliAnalysisTaskPWGJEQA::GenerateHistoBins()
{
  fNPtHistBins = 82;
  fPtHistBins = new Double_t[fNPtHistBins+1];
  GenerateFixedBinArray(6, 0, 0.3, fPtHistBins);
  GenerateFixedBinArray(7, 0.3, 1, fPtHistBins+6);
  GenerateFixedBinArray(10, 1, 3, fPtHistBins+13);
  GenerateFixedBinArray(14, 3, 10, fPtHistBins+23);
  GenerateFixedBinArray(10, 10, 20, fPtHistBins+37);
  GenerateFixedBinArray(15, 20, 50, fPtHistBins+47);
  GenerateFixedBinArray(20, 50, 150, fPtHistBins+62);
  
  fNEtaHistBins = 100;
  fEtaHistBins = new Double_t[fNEtaHistBins+1];
  GenerateFixedBinArray(fNEtaHistBins, -1, 1, fEtaHistBins);
  
  fNPhiHistBins = 101;
  fPhiHistBins = new Double_t[fNPhiHistBins+1];
  GenerateFixedBinArray(fNPhiHistBins, 0, TMath::Pi() * 2.02, fPhiHistBins);
  
  fNCentHistBins = 4;
  fCentHistBins = new Double_t[fNCentHistBins+1];
  fCentHistBins[0] = 0;
  fCentHistBins[1] = 10;
  fCentHistBins[2] = 30;
  fCentHistBins[3] = 50;
  fCentHistBins[4] = 90;
  
  fNPtResHistBins = 175;
  fPtResHistBins = new Double_t[fNPtResHistBins+1];
  GenerateFixedBinArray(50, 0, 0.05, fPtResHistBins);
  GenerateFixedBinArray(25, 0.05, 0.10, fPtResHistBins+50);
  GenerateFixedBinArray(25, 0.10, 0.20, fPtResHistBins+75);
  GenerateFixedBinArray(30, 0.20, 0.50, fPtResHistBins+100);
  GenerateFixedBinArray(25, 0.50, 1.00, fPtResHistBins+130);
  GenerateFixedBinArray(20, 1.00, 2.00, fPtResHistBins+155);
  
  fNPtRelDiffHistBins = 200;
  fPtRelDiffHistBins = new Double_t[fNPtRelDiffHistBins+1];
  GenerateFixedBinArray(fNPtRelDiffHistBins, -2, 2, fPtRelDiffHistBins);
  
  fN1OverPtResHistBins = 385;
  f1OverPtResHistBins = new Double_t[fN1OverPtResHistBins+1];
  GenerateFixedBinArray(100, 0, 0.02, f1OverPtResHistBins);
  GenerateFixedBinArray(60, 0.02, 0.05, f1OverPtResHistBins+100);
  GenerateFixedBinArray(50, 0.05, 0.1, f1OverPtResHistBins+160);
  GenerateFixedBinArray(50, 0.1, 0.2, f1OverPtResHistBins+210);
  GenerateFixedBinArray(75, 0.2, 0.5, f1OverPtResHistBins+260);
  GenerateFixedBinArray(50, 0.5, 1.5, f1OverPtResHistBins+335);
  
  fNIntegerHistBins = 10;
  fIntegerHistBins = new Double_t[fNIntegerHistBins+1];
  GenerateFixedBinArray(fNIntegerHistBins, -0.5, 9.5, fIntegerHistBins);
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateDetectorLevelTHnSparse()
{
  Int_t dim = 0;
  TString title[20];
  Int_t nbins[20] = {0};
  Double_t *binEdges[20] = {0};
  
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    title[dim] = "Centrality %";
    nbins[dim] = fNCentHistBins;
    binEdges[dim] = fCentHistBins;
    dim++;
  }
  
  title[dim] = "#it{p}_{T} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  dim++;
  
  title[dim] = "#eta";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  dim++;
  
  title[dim] = "#phi";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  dim++;
  
  title[dim] = "track type";
  nbins[dim] = 3;
  binEdges[dim] = fIntegerHistBins;
  dim++;
  
  title[dim] = "#sigma(#it{p}_{T}) / #it{p}_{T}";
  nbins[dim] = fNPtResHistBins;
  binEdges[dim] = fPtResHistBins;
  dim++;
  
  fTracks = new THnSparseF("tracks","tracks",dim,nbins);
  for (Int_t i = 0; i < dim; i++) {
    fTracks->GetAxis(i)->SetTitle(title[i]);
    fTracks->SetBinEdges(i, binEdges[i]);
  }
  
  fOutput->Add(fTracks);
  
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateGeneratorLevelTHnSparse()
{
  Int_t dim = 0;
  TString title[20];
  Int_t nbins[20] = {0};
  Double_t *binEdges[20] = {0};
  
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    title[dim] = "Centrality %";
    nbins[dim] = fNCentHistBins;
    binEdges[dim] = fCentHistBins;
    dim++;
  }
  
  title[dim] = "#it{p}_{T} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  dim++;
  
  title[dim] = "#eta";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  dim++;
  
  title[dim] = "#phi";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  dim++;
  
  title[dim] = "Findable";
  nbins[dim] = 2;
  binEdges[dim] = fIntegerHistBins;
  dim++;
  
  fParticlesPhysPrim = new THnSparseF("tracks_PhysPrim","tracks_PhysPrim",dim,nbins);
  for (Int_t i = 0; i < dim; i++) {
    fParticlesPhysPrim->GetAxis(i)->SetTitle(title[i]);
    fParticlesPhysPrim->SetBinEdges(i, binEdges[i]);
  }
  
  fOutput->Add(fParticlesPhysPrim);
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateMatchedParticlesTHnSparse()
{
  Int_t dim = 0;
  TString title[20];
  Int_t nbins[20] = {0};
  Double_t *binEdges[20] = {0};
  
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    title[dim] = "Centrality %";
    nbins[dim] = fNCentHistBins;
    binEdges[dim] = fCentHistBins;
    dim++;
  }
  
  title[dim] = "#it{p}_{T}^{gen} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  dim++;
  
  title[dim] = "#eta^{gen}";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  dim++;
  
  title[dim] = "#phi^{gen}";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  dim++;
  
  title[dim] = "#it{p}_{T}^{det} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  dim++;
  
  title[dim] = "#eta^{det}";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  dim++;
  
  title[dim] = "#phi^{det}";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  dim++;
  
  title[dim] = "(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{gen}";
  nbins[dim] = fNPtRelDiffHistBins;
  binEdges[dim] = fPtRelDiffHistBins;
  dim++;

  title[dim] = "track type";
  nbins[dim] = 3;
  binEdges[dim] = fIntegerHistBins;
  dim++;
  
  fParticlesMatched = new THnSparseF("tracks_Matched","tracks_Matched",dim,nbins);
  for (Int_t i = 0; i < dim; i++) {
    fParticlesMatched->GetAxis(i)->SetTitle(title[i]);
    fParticlesMatched->SetBinEdges(i, binEdges[i]);
  }
  
  fOutput->Add(fParticlesMatched);
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::ExecOnce()
{
  if (fClusterCollArray.GetEntriesFast() == 0 && fCaloCellsName.IsNull()) {
    fNeedEmcalGeom = kFALSE;
  }
  else {
    fNeedEmcalGeom = kTRUE;
  }
  
  AliAnalysisTaskEmcalJet::ExecOnce();
  
  // Load the PHOS geometry
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

//________________________________________________________________________
Bool_t AliAnalysisTaskPWGJEQA::RetrieveEventObjects()
{
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects()) return kFALSE;
  
  // If Pt-hard production, get the Pt-hard of the event, and have possibility to reject the event for jet outliers
  if (fIsPtHard) {
    AliGenPythiaEventHeader* pygen;
    if (MCEvent()) {
      pygen = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
      if (!pygen) {
        // Check if AOD
        AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));
        if (aodMCH) {
          for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
            pygen = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
            if (pygen) break;
          }
        }
      }
    }
    if (pygen) {
      fPtHard = pygen->GetPtHard();
      //fNTrials = pygen->Trials();
      //fXsection = pygen->GetXsection();
      
      // reject outlier events, where the jet Pt is much larger than Pt-hard
      if (fRejectOutlierEvents) {
        AliTLorentzVector jet;
        Int_t nTriggerJets =  pygen->NTriggerJets();
        Float_t tmpjet[]={0,0,0,0};
        for (Int_t ijet = 0; ijet< nTriggerJets; ijet++) {
          pygen->TriggerJet(ijet, tmpjet);
          jet.SetPxPyPzE(tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3]);

          if (jet.Pt() > 4. * fPtHard) {
            //AliInfo(Form("Reject jet event with: pT Hard %2.2f, pycell jet pT %2.2f, rejection factor %1.1f\n", fPtHard, jet.Pt(), 4.));
            if (fGeneralHistograms) fHistEventRejection->Fill("PtHardBinJetOutlier",1);
            return kFALSE;
          }
        }
      }
    }
  }

  return kTRUE;
}

/**
 * Notifying the user that the input data file has
 * changed. Fill pt-hard histograms if applicable.
 *
 * @return False if the data tree or the data file
 * doesn't exist, true otherwise
 */
Bool_t AliAnalysisTaskPWGJEQA::UserNotify()
{
  if (fIsPtHard) {
    TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
    if (!tree) {
      AliError(Form("%s - UserNotify: No current tree!",GetName()));
      return kFALSE;
    }
    
    Float_t xsection    = 0;
    Float_t trials      = 0;
    Int_t   pthardbin   = 0;
    
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      AliError(Form("%s - UserNotify: No current file!",GetName()));
      return kFALSE;
    }
    
    TChain *chain = dynamic_cast<TChain*>(tree);
    if (chain) tree = chain->GetTree();
    
    PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthardbin);
    
    fHistManager.FillTH1("hNtrials", "#sum{ntrials}", trials);
    fHistManager.FillTH1("hXsec", "<#sigma>", xsection);
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPWGJEQA::FillHistograms()
{
  if (fDoTrackQA) FillTrackHistograms();
  if (fDoEmcalQA) FillCellHistograms();
  if (fDoEmcalQA) FillClusterHistograms();
  if (fDoJetQA) FillJetHistograms();
  if (fDoEventQA) FillEventQAHistograms();
  
  return kTRUE;
}
  
//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillTrackHistograms() {

  fNTotTracks = 0;
  fLeadingTrack.SetPxPyPzE(0,0,0,0);
  
  AliVTrack* track;
  for (auto trackIterator : fDetectorLevel->accepted_momentum() ) {

    fNTotTracks++;
    if (fLeadingTrack.Pt() < trackIterator.first.Pt()) fLeadingTrack = trackIterator.first;
    
    track = trackIterator.second;
    Byte_t type = fDetectorLevel->GetTrackType(track);
    if (type <= 2) {
      Double_t sigma = 0;
      
      if (fIsEsd) {
        
        AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
        if (esdTrack) sigma = TMath::Sqrt(esdTrack->GetSigma1Pt2());
        
      } else { // AOD
        
        AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
        if(!aodtrack) AliFatal("Not a standard AOD");
        
        AliExternalTrackParam exParam;
        
        //get covariance matrix
        Double_t cov[21] = {0,};
        aodtrack->GetCovMatrix(cov);
        Double_t pxpypz[3] = {0,};
        aodtrack->PxPyPz(pxpypz);
        Double_t xyz[3] = {0,};
        aodtrack->GetXYZ(xyz);
        Short_t sign = aodtrack->Charge();
        exParam.Set(xyz,pxpypz,cov,sign);
        
        sigma = TMath::Sqrt(exParam.GetSigma1Pt2());
        
      }
      
      Int_t label = TMath::Abs(track->GetLabel());
      Int_t mcGen = 1;
      // reject particles generated from other generators in the cocktail but keep fake tracks (label == 0)
      if (label==0 || track->GetGeneratorIndex() == 0) mcGen = 0;
      
      FillDetectorLevelTHnSparse(fCent, track->Eta(), track->Phi(), track->Pt(), sigma, type);
      
      if (fGeneratorLevel && label > 0) {
        AliAODMCParticle *part =  fGeneratorLevel->GetAcceptMCParticleWithLabel(label);
        if (part) {
          if (part->GetGeneratorIndex() == 0) {
            Int_t pdg = TMath::Abs(part->PdgCode());
            // select charged pions, protons, kaons, electrons, muons
            if (pdg == 211 || pdg == 2212 || pdg == 321 || pdg == 11 || pdg == 13) {
              FillMatchedParticlesTHnSparse(fCent, part->Eta(), part->Phi(), part->Pt(), track->Eta(), track->Phi(), track->Pt(), type);
            }
          }
        }
      }
    }
    else {
      AliError(Form("Track %d has type %d not recognized!", fDetectorLevel->GetCurrentID(), type));
    }
  }
  
  if (fGeneratorLevel) {
    AliAODMCParticle* part;
    for (auto partIterator : fGeneratorLevel->accepted_momentum() ) {
      part = partIterator.second;
      
      Int_t mcGen = 1;
      Byte_t findable = 0;
      
      if (part->GetGeneratorIndex() == 0) mcGen = 0;
      
      Int_t pdg = TMath::Abs(part->PdgCode());
      // select charged pions, protons, kaons, electrons, muons
      if (pdg == 211 || pdg == 2212 || pdg == 321 || pdg == 11 || pdg == 13) findable = 1;
      
      FillGeneratorLevelTHnSparse(fCent, part->Eta(), part->Phi(), part->Pt(), findable);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillCellHistograms() {
  
  if (!fCaloCells) return;
  
  TString histname_energy = TString::Format("%s/fHistCellEnergy", fCaloCellsName.Data());
  TString histname_time = TString::Format("%s/fHistCellTime", fCaloCellsName.Data());
  TString histname_energyProf = TString::Format("%s/fProfCellAbsIdEnergy", fCaloCellsName.Data());
  TString histname_timeProf = TString::Format("%s/fProfCellAbsIdTime", fCaloCellsName.Data());
  
  const Int_t ncells = fCaloCells->GetNumberOfCells();
  for (Int_t pos = 0; pos < ncells; pos++) {
    Float_t amp   = fCaloCells->GetAmplitude(pos);
    Float_t time   = fCaloCells->GetTime(pos);
    Int_t   absId = fCaloCells->GetCellNumber(pos);
    
    if (amp < fCellEnergyCut) continue;
    
    fHistManager.FillTH1(histname_energy, amp);
    fHistManager.FillTH1(histname_time, time);
    
    fHistManager.FillProfile(histname_energyProf, absId, amp);
    fHistManager.FillProfile(histname_timeProf, absId, time);
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillClusterHistograms() {
  
  memset(fNTotClusters, 0, sizeof(Int_t)*3);
  for (Int_t i = 0; i < 3; i++) fLeadingCluster[i].SetPxPyPzE(0,0,0,0);
  
  TString histname;
  AliClusterContainer* clusters = 0;
  TIter nextClusColl(&fClusterCollArray);
  while ((clusters = static_cast<AliClusterContainer*>(nextClusColl()))) {
    // Cluster loop
    AliClusterIterableMomentumContainer itcont = clusters->all_momentum();
    for (AliClusterIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
      
      UInt_t rejectionReason = 0;
      if (!clusters->AcceptCluster(it.current_index(), rejectionReason)) {
        histname = TString::Format("%s/fHistClusterRejectionReason", clusters->GetArrayName().Data());
        fHistManager.FillTH2(histname, clusters->GetRejectionReasonBitPosition(rejectionReason), it->first.E());
        continue;
      }
      
      // Determine cluster type (EMCal/DCal/Phos) and record relevant eventQA info
      ClusterType clusType = kNA;
      if (it->second->IsEMCAL()) {
        Double_t phi = it->first.Phi_0_2pi();
        Int_t isDcal = Int_t(phi > fgkEMCalDCalPhiDivide);
        if (isDcal == 0) {
          clusType = kEMCal;
        } else if (isDcal == 1) {
          clusType = kDCal;
        }
        
        if (fLeadingCluster[isDcal].E() < it->first.E()) fLeadingCluster[isDcal] = it->first;
        fNTotClusters[isDcal]++;
        
        Int_t sm = fGeom->GetSuperModuleNumber(it->second->GetCellAbsId(0));
        if (sm >=0 && sm < 20) {
          histname = TString::Format("%s/BySM/hEmcalClusEnergy_SM%d", clusters->GetArrayName().Data(), sm);
          fHistManager.FillTH1(histname, it->second->E());
        }
        else {
          AliError(Form("Supermodule %d does not exist!", sm));
        }
        
      } else if (it->second->GetType() == AliVCluster::kPHOSNeutral){
        clusType = kPHOS;
        fNTotClusters[2]++;
        if (fLeadingCluster[2].E() < it->first.E()) fLeadingCluster[2] = it->first;
        
        // Fill Phos spectra by module
        Int_t relid[4];
        if (fPHOSGeo) {
          fPHOSGeo->AbsToRelNumbering(it->second->GetCellAbsId(0), relid);
          Int_t sm = relid[0];
          if (sm >=0 && sm < 6) {
            histname = TString::Format("%s/BySM/hPhosClusEnergy_SM%d", clusters->GetArrayName().Data(), sm);
            fHistManager.FillTH1(histname, it->second->E());
          }
          else {
            AliError(Form("Supermodule %d does not exist!", sm));
          }
        }
      }
      
      Double_t contents[30]={0};
      histname = TString::Format("%s/clusterObservables", clusters->GetArrayName().Data());
      THnSparse* histClusterObservables = static_cast<THnSparse*>(fHistManager.FindObject(histname));
      if (!histClusterObservables) return;
      for (Int_t i = 0; i < histClusterObservables->GetNdimensions(); i++) {
        TString title(histClusterObservables->GetAxis(i)->GetTitle());
        if (title=="Centrality %")
          contents[i] = fCent;
        else if (title=="#it{E}_{clus} (GeV)")
          contents[i] = it->first.E();
        else if (title=="#eta")
          contents[i] = it->first.Eta();
        else if (title=="#phi")
          contents[i] = it->first.Phi_0_2pi();
        else if (title=="cluster type")
          contents[i] = clusType;
        else
          AliWarning(Form("Unable to fill dimension %s!",title.Data()));
      }
      histClusterObservables->Fill(contents);
      
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillJetHistograms() {
  
  TString histname;
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
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
      
      Double_t contents[30]={0};
      histname = TString::Format("%s/fHistJetObservables", jets->GetArrayName().Data());
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
        else if (title=="#it{p}_{T}^{MC} (GeV/#it{c})")
          contents[i] = jet->MCPt();
        else if (title=="#it{p}_{T}^{corr} (GeV/#it{c})")
          contents[i] = corrPt;
        else if (title=="NEF")
          contents[i] = jet->NEF();
        else if (title=="No. of constituents")
          contents[i] = jet->GetNumberOfConstituents();
        else if (title=="#it{p}_{T,particle}^{leading} (GeV/#it{c})")
          contents[i] = ptLeading;
        else
          AliWarning(Form("Unable to fill dimension %s!",title.Data()));
      }
      histJetObservables->Fill(contents);
      
    } //jet loop
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillEventQAHistograms() {
  
  EventQA_t eventQA;
  for (Int_t i = 0; i < 3; i++) {
    eventQA.fMaxCluster[i] = fLeadingCluster[i];
  }
  eventQA.fCent = fCent;
  eventQA.fNTracks = fNTotTracks;
  eventQA.fNClusters[0] = fNTotClusters[0];
  eventQA.fNClusters[1] = fNTotClusters[1];
  eventQA.fNClusters[2] = fNTotClusters[2];
  eventQA.fMaxTrack = fLeadingTrack;
  
  Int_t globalNclusters = eventQA.fNClusters[0] + eventQA.fNClusters[1] + eventQA.fNClusters[2];
  AliTLorentzVector globalMaxCluster;
  for (Int_t i = 0; i < 3; i++) {
    if (globalMaxCluster.E() < eventQA.fMaxCluster[i].E()) globalMaxCluster = eventQA.fMaxCluster[i];
  }
  
  THnSparse* histEventQA = static_cast<THnSparse*>(fHistManager.FindObject("eventQA"));
  Double_t contents[40]={0};
  for (Int_t i = 0; i < histEventQA->GetNdimensions(); i++) {
    TString title(histEventQA->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = eventQA.fCent;
    else if (title=="No. of tracks")
      contents[i] = eventQA.fNTracks;
    else if (title=="No. of clusters")
      contents[i] = globalNclusters;
    else if (title=="#it{p}_{T,track}^{leading} (GeV/c)")
      contents[i] = eventQA.fMaxTrack.Pt();
    else if (title=="#it{E}_{cluster}^{leading} (GeV)")
      contents[i] = globalMaxCluster.E();
    else if (title=="#it{E}_{EMCal cluster}^{leading} (GeV)")
      contents[i] = eventQA.fMaxCluster[0].E();
    else if (title=="#it{E}_{DCal cluster}^{leading} (GeV)")
      contents[i] = eventQA.fMaxCluster[1].E();
    else if (title=="#it{E}_{PHOS cluster}^{leading} (GeV)")
      contents[i] = eventQA.fMaxCluster[2].E();
    else
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }
  
  histEventQA->Fill(contents);
  
  // Fill pythia pt hard histograms, if applicable
  if (fPtHard > 1e-6) {
    fHistManager.FillTH1("hPtHard", fPtHard);
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillDetectorLevelTHnSparse(Double_t cent, Double_t trackEta, Double_t trackPhi, Double_t trackPt,
                                                        Double_t sigma1OverPt, Byte_t trackType)
{
  Double_t contents[20]={0};
  
  for (Int_t i = 0; i < fTracks->GetNdimensions(); i++) {
    TString title(fTracks->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T} (GeV/#it{c})")
      contents[i] = trackPt;
    else if (title=="#eta")
      contents[i] = trackEta;
    else if (title=="#phi")
      contents[i] = trackPhi;
    else if (title=="#sigma(1/#it{p}_{T}) (GeV/#it{c})^{-1}")
      contents[i] = sigma1OverPt;
    else if (title=="#sigma(#it{p}_{T}) / #it{p}_{T}")
      contents[i] = sigma1OverPt*trackPt;
    else if (title=="track type")
      contents[i] = trackType;
    else
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fTracks->GetName()));
  }
  
  fTracks->Fill(contents);
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillGeneratorLevelTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt, Byte_t findable)
{
  Double_t contents[20]={0};
  
  for (Int_t i = 0; i < fParticlesPhysPrim->GetNdimensions(); i++) {
    TString title(fParticlesPhysPrim->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T} (GeV/#it{c})")
      contents[i] = partPt;
    else if (title=="#eta")
      contents[i] = partEta;
    else if (title=="#phi")
      contents[i] = partPhi;
    else if (title=="Findable")
      contents[i] = findable;
    else
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fParticlesPhysPrim->GetName()));
  }
  
  fParticlesPhysPrim->Fill(contents);
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillMatchedParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt,
                                                           Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType)
{
  Double_t contents[20]={0};
  
  for (Int_t i = 0; i < fParticlesMatched->GetNdimensions(); i++) {
    TString title(fParticlesMatched->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title=="#it{p}_{T}^{gen} (GeV/#it{c})")
      contents[i] = partPt;
    else if (title=="#eta^{gen}")
      contents[i] = partEta;
    else if (title=="#phi^{gen}")
      contents[i] = partPhi;
    else if (title=="#it{p}_{T}^{det} (GeV/#it{c})")
      contents[i] = trackPt;
    else if (title=="#eta^{det}")
      contents[i] = trackEta;
    else if (title=="#phi^{det}")
      contents[i] = trackPhi;
    else if (title=="(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{gen}")
      contents[i] = (partPt - trackPt) / partPt;
    else if (title=="(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{det}")
      contents[i] = (partPt - trackPt) / trackPt;
    else if (title=="track type")
      contents[i] = (Double_t)trackType;
    else
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), fParticlesMatched->GetName()));
  }
  
  fParticlesMatched->Fill(contents);
}
