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

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskPWGJEQA);
/// \endcond

//________________________________________________________________________
/// Default constructor for ROOT I/O purposes
AliAnalysisTaskPWGJEQA::AliAnalysisTaskPWGJEQA() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskPWGJEQA", kTRUE),
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fCellEnergyCut(0.05),
  fMaxPt(250),
  fNTotTracks(0),
  fLeadingTrack(),
  fDoTrackQA(kTRUE),
  fDoCaloQA(kTRUE),
  fDoJetQA(kTRUE),
  fDoEventQA(kTRUE),
  fGeneratorLevelName(),
  fDetectorLevelName(),
  fRejectOutlierEvents(kFALSE),
  fIsPtHard(kFALSE),
  fGeneratorLevel(0),
  fDetectorLevel(0),
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
  fNIntegerHistBins(0),
  fIntegerHistBins(0),
  fPHOSGeo(nullptr),
  fHistManager("AliAnalysisTaskPWGJEQA")
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
  fUseAliEventCuts(kTRUE),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fCellEnergyCut(0.05),
  fMaxPt(250),
  fNTotTracks(0),
  fLeadingTrack(),
  fDoTrackQA(kTRUE),
  fDoCaloQA(kTRUE),
  fDoJetQA(kTRUE),
  fDoEventQA(kTRUE),
  fGeneratorLevelName(),
  fDetectorLevelName(),
  fRejectOutlierEvents(kFALSE),
  fIsPtHard(kFALSE),
  fGeneratorLevel(0),
  fDetectorLevel(0),
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
  fNIntegerHistBins(0),
  fIntegerHistBins(0),
  fPHOSGeo(nullptr),
  fHistManager(name)
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

  // Intialize AliEventCuts
  if (fUseAliEventCuts) {
    fEventCutList = new TList();
    fEventCutList ->SetOwner();
    fEventCutList ->SetName("EventCutOutput");

    fEventCuts.OverrideAutomaticTriggerSelection(fOffTrigger);
    if(fUseManualEventCuts==1)
    {
    	  fEventCuts.SetManualMode();
      fEventCuts.fMC = MCEvent(); //before was= false
      if(fForceBeamType != kpp)
      {
    		fEventCuts.SetupLHC15o();
    		fEventCuts.fUseVariablesCorrelationCuts = true;
      }
      else if(fForceBeamType == kpp)
      {
    		fEventCuts.SetupRun2pp();
    		//no other cuts known so far
      }
      else
      {
    		printf("No implementation of manuel event cuts for pPb yet!");
      }
    }
    fEventCuts.AddQAplotsToList(fEventCutList);
    fOutput->Add(fEventCutList);
  }
  
  // Set track container pointers
  fDetectorLevel  = GetTrackContainer(fDetectorLevelName);
  fGeneratorLevel = GetMCParticleContainer(fGeneratorLevelName);
  
  // Allocate histograms for tracks, cells, clusters, jets
  if (fDoTrackQA) AllocateTrackHistograms();
  if (fDoCaloQA)  AllocateCellHistograms();
  if (fDoCaloQA)  AllocateClusterHistograms();
  if (fDoJetQA)   AllocateJetHistograms();
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
    fHistManager.CreateTH1(histname.Data(), title.Data(), 500,-3e-6, 3e-6);
    
    histname = TString::Format("%s/fProfCellAbsIdTime", fCaloCellsName.Data());
    title = histname + ";cell absId;<#it{t}_{cell}> (s)";
    fHistManager.CreateTProfile(histname.Data(), title.Data(), 18000,0,18000);
    
    histname = TString::Format("%s/fHistCellEvsTime", fCaloCellsName.Data());
    title = histname + ";<#it{t}_{cell}> (s);<#it{E}_{cell}> (GeV)";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 1000, -1e-6, 1e-6, 200, 0, 100);
    
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateClusterHistograms() {

  AliEmcalContainer* cont = 0;
  TString histname;
  TString title;
  Int_t nPtBins = 2*fMaxPt;
  
  TIter nextClusColl(&fClusterCollArray);
  while ((cont = static_cast<AliEmcalContainer*>(nextClusColl()))) {

    histname = TString::Format("%s/fHistClusterRejectionReason", cont->GetArrayName().Data());
    title = histname + ";Rejection reason;#it{E}_{clus} (GeV/);counts";
    TH2* hist = fHistManager.CreateTH2(histname.Data(), title.Data(), 32, 0, 32, 50, 0, 250);
    SetRejectionReasonLabels(hist->GetXaxis());
    
    histname = TString::Format("%s/hPhosNCellsVsEnergy", cont->GetArrayName().Data());
    title = histname + ";N cells;#it{E}_{clus} (GeV/);counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 40, 0, 40, nPtBins, 0, fMaxPt);
    
    histname = TString::Format("%s/hPhosM02VsEnergy", cont->GetArrayName().Data());
    title = histname + ";M02;#it{E}_{clus} (GeV/);counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 100, 0, 20, nPtBins, 0, fMaxPt);
    
    histname = TString::Format("%s/hPhosCellIdVsEnergy", cont->GetArrayName().Data());
    title = histname + ";Cell ID;#it{E}_{clus} (GeV/);counts";
    fHistManager.CreateTH2(histname.Data(), title.Data(), 20000, 0, 20000, nPtBins, 0, fMaxPt);
    
    const Int_t nEmcalSM = 20;
    for (Int_t sm = 0; sm < nEmcalSM; sm++) {
      histname = TString::Format("%s/BySM/hEmcalClusEnergy_SM%d", cont->GetArrayName().Data(), sm);
      title = histname + ";#it{E}_{cluster} (GeV);counts";
      fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);
    }
    
    for (Int_t sm = 1; sm < 5; sm++) {
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
    
    title[dim] = "#eta";
    nbins[dim] = fNEtaHistBins;
    binEdges[dim] = fEtaHistBins;
    min[dim] = fEtaHistBins[0];
    max[dim] = fEtaHistBins[fNEtaHistBins];
    dim++;
    
    title[dim] = "#phi";
    nbins[dim] = fNPhiHistBins;
    binEdges[dim] = fPhiHistBins;
    min[dim] = fPhiHistBins[0];
    max[dim] = fPhiHistBins[fNPhiHistBins];
    dim++;
    
    title[dim] = "cluster type";
    nbins[dim] = 3;
    binEdges[dim] = fIntegerHistBins;
    min[dim] = fIntegerHistBins[0];
    max[dim] = fIntegerHistBins[3];
    dim++;
    
    TString thnname = TString::Format("%s/clusterObservables", cont->GetArrayName().Data());
    THnSparse* hn = fHistManager.CreateTHnSparse(thnname.Data(), thnname.Data(), dim, nbins, min, max);
    for (Int_t i = 0; i < dim; i++) {
      hn->GetAxis(i)->SetTitle(title[i]);
      hn->SetBinEdges(i, binEdges[i]);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateJetHistograms() {
  
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
    nbins[dim] = fNEtaHistBins;
    binEdges[dim] = fEtaHistBins;
    min[dim] = fEtaHistBins[0];
    max[dim] = fEtaHistBins[fNEtaHistBins];
    dim++;
    
    axisTitle[dim] = "#phi_{jet} (rad)";
    nbins[dim] = fNPhiHistBins;
    binEdges[dim] = fPhiHistBins;
    min[dim] = fPhiHistBins[0];
    max[dim] = fPhiHistBins[fNPhiHistBins];
    dim++;
    
    axisTitle[dim] = "#it{p}_{T} (GeV/#it{c})";
    nbins[dim] = TMath::CeilNint(fMaxPt/3);
    min[dim] = 0;
    max[dim] = fMaxPt;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    if (fForceBeamType != kpp) {
      axisTitle[dim] = "#it{p}_{T}^{corr} (GeV/#it{c})";
      nbins[dim] = TMath::CeilNint(fMaxPt/3);
      min[dim] = -fMaxPt/2 + 25;
      max[dim] = fMaxPt/2 + 25;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
      dim++;
    }
    
    axisTitle[dim] = "#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    nbins[dim] = fMaxPt/5;
    min[dim] = 0;
    max[dim] = 150;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
    
    TString thnname = TString::Format("%s/fHistJetObservables", jets->GetArrayName().Data());
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
    if (fForceBeamType != kpp) {
    	  histname = TString::Format("%s/hNEFVsPtEMC", jets->GetArrayName().Data());
      title = histname + ";Centrality (%); #it{p}_{T}^{corr} (GeV/#it{c});NEF";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 20, 0, 100, 250, 0, 250, 50, 0, 1);

      histname = TString::Format("%s/hNEFVsPtDCal", jets->GetArrayName().Data());
      title = histname + ";Centrality (%); #it{p}_{T}^{corr} (GeV/#it{c});NEF";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 20, 0, 100, 250, 0, 250, 50, 0, 1);

      histname = TString::Format("%s/hDeltaEHadCorr", jets->GetArrayName().Data());
      title = histname + ";Centrality (%); #it{p}_{T}^{corr} (GeV/#it{c}); #sum #it{E_{nonlincorr}} - #it{E}_{hadcorr}";
      fHistManager.CreateTH3(histname.Data(), title.Data(), 20, 0, 100, 250, 0, 250, 250, 0, 50);
    }
    else{
    	  histname = TString::Format("%s/hNEFVsPtEMC", jets->GetArrayName().Data());
    	  title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});NEF";
    	  fHistManager.CreateTH2(histname.Data(), title.Data(), 250, 0, 250, 50, 0, 1);

    	  histname = TString::Format("%s/hNEFVsPtDCal", jets->GetArrayName().Data());
    	  title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c});NEF";
    	  fHistManager.CreateTH2(histname.Data(), title.Data(), 250, 0, 250, 50, 0, 1);

    	  histname = TString::Format("%s/hDeltaEHadCorr", jets->GetArrayName().Data());
    	  title = histname + ";#it{p}_{T}^{corr} (GeV/#it{c}); #sum#it{E_{nonlincorr}} - #it{E}_{hadcorr}";
    	  fHistManager.CreateTH2(histname.Data(), title.Data(), 250, 0, 250, 250, 0, 50);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateEventQAHistograms() {
  
  Int_t nPtBins = 2*fMaxPt;
  
  Int_t dim = 0;
  TString axistitle[40];
  Int_t nbins[40] = {0};
  Double_t min[40] = {0};
  Double_t max[40] = {0};
  Double_t *binEdges[20] = {0};
  
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
    axistitle[dim] = "Centrality %";
    nbins[dim] = fNCentHistBins;
    binEdges[dim] = fCentHistBins;
    min[dim] = fCentHistBins[0];
    max[dim] = fCentHistBins[fNCentHistBins];
    dim++;
  }
  
  if (fParticleCollArray.GetEntriesFast()>0) {
    axistitle[dim] = "No. of tracks";
    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
      nbins[dim] = 50;
      min[dim] = -0.5;
      max[dim] = 10000-0.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    }
    else {
      nbins[dim] = 50;
      min[dim] = 0;
      max[dim] = 200;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    }
    dim++;
    
    axistitle[dim] = "#it{p}_{T,track}^{leading} (GeV/c)";
    nbins[dim] = fMaxPt/5;
    min[dim] = 0;
    max[dim] = 150;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
  }
  
  if (fClusterCollArray.GetEntriesFast()>0) {
    axistitle[dim] = "No. of clusters";
    
    if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
      nbins[dim] = 50;
      min[dim] = -0.5;
      max[dim] = 4000-0.5;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    }
    else {
      nbins[dim] = 50;
      min[dim] = 0;
      max[dim] = 200;
      binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    }
    dim++;
    
    axistitle[dim] = "#it{E}_{cluster}^{leading} (GeV)";
    nbins[dim] = fMaxPt/5;
    min[dim] = 0;
    max[dim] = 150;
    binEdges[dim] = GenerateFixedBinArray(nbins[dim], min[dim], max[dim]);
    dim++;
  }
  
  THnSparse* hn = fHistManager.CreateTHnSparse("eventQA","eventQA",dim,nbins,min,max);
  for (Int_t i = 0; i < dim; i++) {
    hn->GetAxis(i)->SetTitle(axistitle[i]);
    hn->SetBinEdges(i, binEdges[i]);
  }
    
  if (fIsPtHard) {
    TString histname = "hPtHard";
    TString title = histname + ";#it{p}_{T,hard} (GeV/c);counts";
    fHistManager.CreateTH1(histname.Data(), title.Data(), nPtBins, 0, fMaxPt);
    
    TH1* hTrials = fHistManager.CreateTH1("hNtrials", "hNtrials", 1, 0, 1);
    hTrials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
    
    TH1* hxsection = fHistManager.CreateTH1("hXsec", "hXsec", 1, 0, 1);
    hxsection->GetXaxis()->SetBinLabel(1,"<#sigma>");
  }
  
  if(fUseBuiltinEventSelection)
    fHistEventRejection->GetXaxis()->SetBinLabel(15,"PtHardBinJetOutlier");

}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateDetectorLevelTHnSparse()
{
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
  
  title[dim] = "#it{p}_{T} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  min[dim] = fPtHistBins[0];
  max[dim] = fPtHistBins[fNPtHistBins];
  dim++;
  
  title[dim] = "#eta";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  min[dim] = fEtaHistBins[0];
  max[dim] = fEtaHistBins[fNEtaHistBins];
  dim++;
  
  title[dim] = "#phi";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  min[dim] = fPhiHistBins[0];
  max[dim] = fPhiHistBins[fNPhiHistBins];
  dim++;
  
  title[dim] = "track type";
  nbins[dim] = 3;
  binEdges[dim] = fIntegerHistBins;
  min[dim] = fIntegerHistBins[0];
  max[dim] = fIntegerHistBins[3];
  dim++;
  
  title[dim] = "#sigma(#it{p}_{T}) / #it{p}_{T}";
  nbins[dim] = fNPtResHistBins;
  binEdges[dim] = fPtResHistBins;
  min[dim] = fPtResHistBins[0];
  max[dim] = fPtResHistBins[fNPtResHistBins];
  dim++;

  THnSparse* hn = fHistManager.CreateTHnSparse("tracks","tracks", dim, nbins, min, max);
  for (Int_t i = 0; i < dim; i++) {
    hn->GetAxis(i)->SetTitle(title[i]);
    hn->SetBinEdges(i, binEdges[i]);
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateGeneratorLevelTHnSparse()
{
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
  
  title[dim] = "#it{p}_{T} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  min[dim] = fPtHistBins[0];
  max[dim] = fPtHistBins[fNPtHistBins];
  dim++;
  
  title[dim] = "#eta";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  min[dim] = fEtaHistBins[0];
  max[dim] = fEtaHistBins[fNEtaHistBins];
  dim++;
  
  title[dim] = "#phi";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  min[dim] = fPhiHistBins[0];
  max[dim] = fPhiHistBins[fNPhiHistBins];
  dim++;
  
  title[dim] = "Findable";
  nbins[dim] = 2;
  binEdges[dim] = fIntegerHistBins;
  min[dim] = fIntegerHistBins[0];
  max[dim] = fIntegerHistBins[2];
  dim++;
  
  THnSparse* hn = fHistManager.CreateTHnSparse("tracks_PhysPrim","tracks_PhysPrim", dim, nbins, min, max);
  for (Int_t i = 0; i < dim; i++) {
    hn->GetAxis(i)->SetTitle(title[i]);
    hn->SetBinEdges(i, binEdges[i]);
  }
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::AllocateMatchedParticlesTHnSparse()
{
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
  
  title[dim] = "#it{p}_{T}^{gen} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  min[dim] = fPtHistBins[0];
  max[dim] = fPtHistBins[fNPtHistBins];
  dim++;
  
  title[dim] = "#eta^{gen}";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  min[dim] = fEtaHistBins[0];
  max[dim] = fEtaHistBins[fNEtaHistBins];
  dim++;
  
  title[dim] = "#phi^{gen}";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  min[dim] = fPhiHistBins[0];
  max[dim] = fPhiHistBins[fNPhiHistBins];
  dim++;
  
  title[dim] = "#it{p}_{T}^{det} (GeV/#it{c})";
  nbins[dim] = fNPtHistBins;
  binEdges[dim] = fPtHistBins;
  min[dim] = fPtHistBins[0];
  max[dim] = fPtHistBins[fNPtHistBins];
  dim++;
  
  title[dim] = "#eta^{det}";
  nbins[dim] = fNEtaHistBins;
  binEdges[dim] = fEtaHistBins;
  min[dim] = fEtaHistBins[0];
  max[dim] = fEtaHistBins[fNEtaHistBins];
  dim++;
  
  title[dim] = "#phi^{det}";
  nbins[dim] = fNPhiHistBins;
  binEdges[dim] = fPhiHistBins;
  min[dim] = fPhiHistBins[0];
  max[dim] = fPhiHistBins[fNPhiHistBins];
  dim++;
  
  title[dim] = "(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{gen}";
  nbins[dim] = fNPtRelDiffHistBins;
  binEdges[dim] = fPtRelDiffHistBins;
  min[dim] = fPtRelDiffHistBins[0];
  max[dim] = fPtRelDiffHistBins[fNPtRelDiffHistBins];
  dim++;

  title[dim] = "track type";
  nbins[dim] = 3;
  binEdges[dim] = fIntegerHistBins;
  min[dim] = fIntegerHistBins[0];
  max[dim] = fIntegerHistBins[3];
  dim++;

  THnSparse* hn = fHistManager.CreateTHnSparse("tracks_Matched","tracks_Matched", dim, nbins, min, max);
  for (Int_t i = 0; i < dim; i++) {
    hn->GetAxis(i)->SetTitle(title[i]);
    hn->SetBinEdges(i, binEdges[i]);
  }
}


//________________________________________________________________________
/// Generate histogram binning arrays
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
  
  fNEtaHistBins = 70;
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
  
  fNPtResHistBins = 100;
  fPtResHistBins = new Double_t[fNPtResHistBins+1];
  GenerateFixedBinArray(50, 0, 0.05, fPtResHistBins);
  GenerateFixedBinArray(25, 0.05, 0.10, fPtResHistBins+50);
  GenerateFixedBinArray(15, 0.10, 0.20, fPtResHistBins+75);
  GenerateFixedBinArray(10, 0.20, 0.50, fPtResHistBins+90);
  
  fNPtRelDiffHistBins = 50;
  fPtRelDiffHistBins = new Double_t[fNPtRelDiffHistBins+1];
  GenerateFixedBinArray(fNPtRelDiffHistBins, -0.5, 0.5, fPtRelDiffHistBins);
  
  fNIntegerHistBins = 3;
  fIntegerHistBins = new Double_t[fNIntegerHistBins+1];
  GenerateFixedBinArray(fNIntegerHistBins, -0.5, 2.5, fIntegerHistBins);
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

///
/// This function (overloading the base class) uses AliEventCuts to perform event selection.
///________________________________________________________________________
Bool_t AliAnalysisTaskPWGJEQA::IsEventSelected()
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
//________________________________________________________________________
Bool_t AliAnalysisTaskPWGJEQA::RetrieveEventObjects()
{
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects()) return kFALSE;
  
  // If Pt-hard production, get the Pt-hard of the event, and have possibility to reject the event for jet outliers
  if (fIsPtHard) {
    AliGenPythiaEventHeader* pygen = nullptr;
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
            if (fGeneralHistograms && fHistEventRejection) fHistEventRejection->Fill("PtHardBinJetOutlier",1);
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
//________________________________________________________________________
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
  if (fDoCaloQA) FillCellHistograms();
  if (fDoCaloQA) FillClusterHistograms();
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
      
      Byte_t findable = 0;
      
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
  TString histname_EnergyvsTime = TString::Format("%s/fHistCellEvsTime", fCaloCellsName.Data());
  
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
    
    fHistManager.FillTH2(histname_EnergyvsTime, time, amp);
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
        
        Int_t nCells = it->second->GetNCells();
        Double_t M02 = it->second->GetM02();
        Double_t energy = it->second->E();
        
        histname = TString::Format("%s/hPhosNCellsVsEnergy", clusters->GetArrayName().Data());
        fHistManager.FillTH2(histname, nCells, energy);
        
        histname = TString::Format("%s/hPhosM02VsEnergy", clusters->GetArrayName().Data());
        fHistManager.FillTH2(histname, M02, energy);
        
        clusType = kPHOS;
        fNTotClusters[2]++;
        if (fLeadingCluster[2].E() < it->first.E()) fLeadingCluster[2] = it->first;
        
        // Fill Phos spectra by module
        Int_t relid[4];
        if (fPHOSGeo) {
          fPHOSGeo->AbsToRelNumbering(it->second->GetCellAbsId(0), relid);
          Int_t sm = relid[0];
          if (sm >=1 && sm < 5) {
            histname = TString::Format("%s/BySM/hPhosClusEnergy_SM%d", clusters->GetArrayName().Data(), sm);
            fHistManager.FillTH1(histname, energy);
          }
          else {
            AliError(Form("Supermodule %d does not exist!", sm));
          }
        }
        
        // Loop through cells in the cluster, and histogram their energy
        histname = TString::Format("%s/hPhosCellIdVsEnergy", clusters->GetArrayName().Data());
        for (Int_t i=0; i < nCells; i++) {
          fHistManager.FillTH2(histname, it->second->GetCellAbsId(i), energy);
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
        else if (title=="#it{p}_{T,particle}^{leading} (GeV/#it{c})")
          contents[i] = ptLeading;
        else
          AliWarning(Form("Unable to fill dimension %s!",title.Data()));
      }
      histJetObservables->Fill(contents);
      
      UInt_t jetType = jet->GetJetAcceptanceType();
      if(jetType & AliEmcalJet::kEMCALfid)
      {
    	  	histname = TString::Format("%s/hNEFVsPtEMC", jets->GetArrayName().Data());
    	    if (fForceBeamType != kpp){fHistManager.FillTH3(histname.Data(), fCent, jet->Pt(), jet->NEF());}
    	    else                      {fHistManager.FillTH2(histname.Data(), jet->Pt(), jet->NEF());}
      }
      else if(jetType & AliEmcalJet::kDCALonlyfid)
      {
    	    histname = TString::Format("%s/hNEFVsPtDCal", jets->GetArrayName().Data());
    	    if (fForceBeamType != kpp){
    	    	  fHistManager.FillTH3(histname.Data(), fCent, jet->Pt(), jet->NEF());
    	    }
    	    else{
    	    	  fHistManager.FillTH2(histname.Data(), jet->Pt(), jet->NEF());
    	    }
      }
      Double_t deltaEhadcorr = 0;
      const AliVCluster* clus = nullptr;
      Int_t nClusters = jet->GetNumberOfClusters();
      for (Int_t iClus = 0; iClus < nClusters; iClus++) {
    	    clus = jet->Cluster(iClus);
    	    deltaEhadcorr += (clus->GetNonLinCorrEnergy() - clus->GetHadCorrEnergy());
      }

      histname = TString::Format("%s/hDeltaEHadCorr", jets->GetArrayName().Data());
      if (fForceBeamType != kpp){
    	    fHistManager.FillTH3(histname, fCent, jet->Pt(), deltaEhadcorr);
      }
      else{
    	    fHistManager.FillTH2(histname, jet->Pt(), deltaEhadcorr);
      }
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
  
  THnSparse* thnTracks = static_cast<THnSparse*>(fHistManager.FindObject("tracks"));
  if (!thnTracks) return;
  for (Int_t i = 0; i < thnTracks->GetNdimensions(); i++) {
    TString title(thnTracks->GetAxis(i)->GetTitle());
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
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), thnTracks->GetName()));
  }
  
  thnTracks->Fill(contents);
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillGeneratorLevelTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt, Byte_t findable)
{
  Double_t contents[20]={0};
  
  THnSparse* thnTracks_PhysPrim = static_cast<THnSparse*>(fHistManager.FindObject("tracks_PhysPrim"));
  if (!thnTracks_PhysPrim) return;
  for (Int_t i = 0; i < thnTracks_PhysPrim->GetNdimensions(); i++) {
    TString title(thnTracks_PhysPrim->GetAxis(i)->GetTitle());
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
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), thnTracks_PhysPrim->GetName()));
  }
  
  thnTracks_PhysPrim->Fill(contents);
}

//________________________________________________________________________
void AliAnalysisTaskPWGJEQA::FillMatchedParticlesTHnSparse(Double_t cent, Double_t partEta, Double_t partPhi, Double_t partPt,
                                                           Double_t trackEta, Double_t trackPhi, Double_t trackPt, Byte_t trackType)
{
  Double_t contents[20]={0};
  
  THnSparse* thnTracks_Matched = static_cast<THnSparse*>(fHistManager.FindObject("tracks_Matched"));
  if (!thnTracks_Matched) return;
  for (Int_t i = 0; i < thnTracks_Matched->GetNdimensions(); i++) {
    TString title(thnTracks_Matched->GetAxis(i)->GetTitle());
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
      AliWarning(Form("Unable to fill dimension %s of histogram %s!", title.Data(), thnTracks_Matched->GetName()));
  }
  
  thnTracks_Matched->Fill(contents);
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskPWGJEQA* AliAnalysisTaskPWGJEQA::AddTaskPWGJEQA(
                        const char* ntracks,
                        const char* nclusters,
                        const char* ncells,
                        const char *nGenLev,
                        Bool_t      doTrackQA,
                        Bool_t      doCaloQA,
                        Bool_t      doJetQA,
                        Bool_t      doEventQA,
                        Double_t    trackPtCut,
                        Double_t    clusECut,
                        const char* suffix)
{
	// Get the pointer to the existing analysis manager via the static access method
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		::Error("AddTaskPWGJEQA", "No analysis manager to connect to.");
		return NULL;
	}

	// Check the analysis type using the event handlers connected to the analysis manager
	AliVEventHandler* handler = mgr->GetInputEventHandler();
	if (!handler) {
		::Error("AddTaskPWGJEQA", "This task requires an input event handler");
		return NULL;
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

	// Init the task and do settings
	TString trackName(ntracks);
	TString clusName(nclusters);
	TString cellName(ncells);

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

	if (cellName == "usedefault") {
		if (dataType == kESD) {
			cellName = "EMCALCells";
		}
		else if (dataType == kAOD) {
			cellName = "emcalCells";
		}
		else {
			cellName = "";
		}
	}

	TString name("AliAnalysisTaskPWGJEQA");
	if (!trackName.IsNull()) {
		name += "_";
		name += trackName;
	}
	if (!clusName.IsNull()) {
		name += "_";
		name += clusName;
	}

	if (!cellName.IsNull()) {
		name += "_";
		name += cellName;
	}

	if (strcmp(suffix,"")) {
		name += "_";
		name += suffix;
	}

	if (nGenLev && strcmp(nGenLev,"")!=0) cout<<"MC Part: "<<  nGenLev<<endl;
	else  cout<<"No MC Part: "<<  nGenLev<<endl;

	AliAnalysisTaskPWGJEQA* qaTask = new AliAnalysisTaskPWGJEQA(name);
	qaTask->SetVzRange(-10,10);
	qaTask->SetNeedEmcalGeom(kFALSE);
	qaTask->SetCaloCellsName(cellName);
	qaTask->SetDetectorLevelName(trackName);
	if (nGenLev && strcmp(nGenLev,"")!=0) qaTask->SetGeneratorLevelName(nGenLev);
	qaTask->SetDoTrackQA(doTrackQA);
	qaTask->SetDoCaloQA(doCaloQA);
	qaTask->SetDoJetQA(doJetQA);
	qaTask->SetDoEventQA(doEventQA);

	// Add the detector-level track container
	if (trackName == "tracks" || trackName == "Tracks") {
		AliTrackContainer* trackCont = qaTask->AddTrackContainer(trackName);
		trackCont->SetFilterHybridTracks(kTRUE);
	}
	else if (!trackName.IsNull()) {
		qaTask->AddParticleContainer(trackName);
	}
	AliParticleContainer *partCont = qaTask->GetParticleContainer(trackName);
	if (partCont) {
		partCont->SetParticlePtCut(trackPtCut);
	}

	// Add the generator-level container, if specified
	if (nGenLev && strcmp(nGenLev,"")!=0) {
		AliMCParticleContainer* mcpartCont = qaTask->AddMCParticleContainer(nGenLev);
		mcpartCont->SelectPhysicalPrimaries(kTRUE);
		mcpartCont->SetParticlePtCut(0);
	}

	// Add the cluster container
	AliClusterContainer *clusterCont = qaTask->AddClusterContainer(clusName);
	if (clusterCont) {
		clusterCont->SetClusECut(0.);
		clusterCont->SetClusPtCut(0.);
		clusterCont->SetClusHadCorrEnergyCut(clusECut);
		clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
	}

	// Final settings, pass to manager and set the containers
	mgr->AddTask(qaTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

	TString contName = TString::Format("%s_histos", name.Data());
	TString commonoutput;
	commonoutput = mgr->GetCommonFileName();
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(),
			TList::Class(),AliAnalysisManager::kOutputContainer,
			commonoutput);
	mgr->ConnectInput  (qaTask, 0,  cinput1 );
	mgr->ConnectOutput (qaTask, 1, coutput1 );

	return qaTask;

}
