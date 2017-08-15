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
  fHistManager(),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fUseAliEventCuts(kTRUE),
  fMaxPt(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNEtaBins(40),
  fNPhiBins(200),
  fPlotNeutralJets(kFALSE),
  fPlotClustersInJets(kFALSE),
  fPlotClusterHistograms(kTRUE),
  fPlotCellHistograms(kTRUE),
  fPlotClusWithoutNonLinCorr(kFALSE),
  fPlotExotics(kFALSE),
  fPHOSGeo(nullptr)
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
  fHistManager(name),
  fEventCuts(0),
  fEventCutList(0),
  fUseManualEventCuts(kFALSE),
  fUseAliEventCuts(kTRUE),
  fMaxPt(200),
  fNCentHistBins(0),
  fCentHistBins(0),
  fNPtHistBins(0),
  fPtHistBins(0),
  fNEtaBins(40),
  fNPhiBins(200),
  fPlotNeutralJets(kFALSE),
  fPlotClustersInJets(kFALSE),
  fPlotClusterHistograms(kTRUE),
  fPlotCellHistograms(kTRUE),
  fPlotClusWithoutNonLinCorr(kFALSE),
  fPlotExotics(kFALSE),
  fPHOSGeo(nullptr)
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
  
}

/*
 * This function allocates the histograms for the calorimeter performance study.
 */
void AliAnalysisTaskEmcalVsPhos::AllocateClusterHistograms()
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
    
    if (fPlotExotics) {
      histname = TString::Format("%s/hFcrossEMCal", cont->GetArrayName().Data());
      htitle = histname + ";Centrality (%);Fcross;#it{E}_{clus} (GeV/)";
      TH3* hist = fHistManager.CreateTH3(histname.Data(), htitle.Data(), fNCentHistBins, fCentHistBins, nExBins, exBins, fNPtHistBins, fPtHistBins);
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
  AliJetContainer* jets = 0;
  TIter nextJetColl(&fJetCollArray);
  while ((jets = static_cast<AliJetContainer*>(nextJetColl()))) {
    
    for (auto jet : jets->accepted()) {

      // Fill cluster spectra of clusters within jets
      //(centrality, cluster energy, jet pT, jet eta, jet phi)
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

/**
 * Get pT of jet -- background subtracted, unless hard-core jet
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetJetPt(AliJetContainer* jetCont, AliEmcalJet* jet)
{
  // Get eta-phi dependent jet pT scale factor
  Double_t jetPt = jet->Pt();
  
  // Compute pTcorr
  Double_t rho = jetCont->GetRhoVal();
  Double_t pT = jetPt - rho * jet->Area();
  
  return pT;
}

/**
 * Get deltaR of a track/cluster and a reference point.
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetDeltaR(AliTLorentzVector* part, Double_t etaRef, Double_t phiRef)
{
  Double_t deltaPhi = TMath::Abs(part->Phi_0_2pi() - phiRef);
  Double_t deltaEta = TMath::Abs(part->Eta() - etaRef);
  Double_t deltaR = TMath::Sqrt( deltaPhi*deltaPhi + deltaEta*deltaEta );
  return deltaR;
}

/**
 * Get calo acceptance type of jet
 */
Double_t AliAnalysisTaskEmcalVsPhos::GetJetType(AliEmcalJet* jet)
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
Double_t AliAnalysisTaskEmcalVsPhos::GetFcross(AliVCluster *cluster, AliVCaloCells *cells)
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
