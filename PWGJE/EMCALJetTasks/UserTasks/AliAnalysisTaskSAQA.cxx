// $Id$
//
// General QA task.
//
// Author: S.Aiola

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliAODEvent.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliPicoTrack.h"
#include "AliVVZERO.h"
#include "AliESDUtils.h"

#include "AliAnalysisTaskSAQA.h"

ClassImp(AliAnalysisTaskSAQA)

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA() : 
  AliAnalysisTaskEmcalJetDev("AliAnalysisTaskSAQA", kTRUE),
  fCellEnergyCut(0.1),
  fParticleLevel(kFALSE),
  fIsMC(kFALSE),
  fCentMethod2(""),
  fCentMethod3(""),
  fDoV0QA(0),
  fDoEPQA(0),
  fMaxCellsInCluster(30),
  fCent2(0),
  fCent3(0),
  fVZERO(0),
  fV0ATotMult(0),
  fV0CTotMult(0),
  fHistEventQA(0),
  fHistTrNegativeLabels(0),
  fHistTrZeroLabels(0),
  fHistNCellsEnergy(0),
  fHistFcrossEnergy(0),
  fHistClusTimeEnergy(0),
  fHistCellsAbsIdEnergy(0),
  fHistChVSneCells(0),
  fHistChVSneClus(0),
  fHistChVSneCorrCells(0)
{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < 4; j++) fHistTrPhiEtaPt[i][j] = 0;
    fHistTrPhiEtaZeroLab[i] = 0;
    fHistTrPtZeroLab[i] = 0;
    fHistTrEmcPhiEta[i] = 0;
    fHistTrEmcPt[i] = 0;
    fHistTrPhiEtaNonProp[i] = 0;
    fHistTrPtNonProp[i] = 0;
    fHistDeltaEtaPt[i] = 0;
    fHistDeltaPhiPt[i] = 0;
    fHistDeltaPtvsPtvsMass[i] = 0;
    fHistClusPhiEtaEnergy[i] = 0;
    fHistClusMCEnergyFraction[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA(const char *name) : 
  AliAnalysisTaskEmcalJetDev(name, kTRUE),
  fCellEnergyCut(0.1),
  fParticleLevel(kFALSE),
  fIsMC(kFALSE),
  fCentMethod2(""),
  fCentMethod3(""),
  fDoV0QA(0),
  fDoEPQA(0),
  fMaxCellsInCluster(30),
  fCent2(0),
  fCent3(0),
  fVZERO(0),
  fV0ATotMult(0),
  fV0CTotMult(0),
  fHistEventQA(0),
  fHistTrNegativeLabels(0),
  fHistTrZeroLabels(0),
  fHistNCellsEnergy(0),
  fHistFcrossEnergy(0),
  fHistClusTimeEnergy(0),
  fHistCellsAbsIdEnergy(0),
  fHistChVSneCells(0),
  fHistChVSneClus(0),
  fHistChVSneCorrCells(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < 4; j++) fHistTrPhiEtaPt[i][j] = 0;
    fHistTrPhiEtaZeroLab[i] = 0;
    fHistTrPtZeroLab[i] = 0;
    fHistTrEmcPhiEta[i] = 0;
    fHistTrEmcPt[i] = 0;
    fHistTrPhiEtaNonProp[i] = 0;
    fHistTrPtNonProp[i] = 0;
    fHistDeltaEtaPt[i] = 0;
    fHistDeltaPhiPt[i] = 0;
    fHistDeltaPtvsPtvsMass[i] = 0;
    fHistClusPhiEtaEnergy[i] = 0;
    fHistClusMCEnergyFraction[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSAQA::~AliAnalysisTaskSAQA()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::UserCreateOutputObjects()
{
  // Create histograms

  AliAnalysisTaskEmcalJetDev::UserCreateOutputObjects();

  if (fParticleCollArray.GetEntriesFast()>0) {
    if (!fParticleLevel && fIsMC) {
      fHistTrNegativeLabels = new TH1F("fHistTrNegativeLabels","fHistTrNegativeLabels", 500, 0, 1);
      fHistTrNegativeLabels->GetXaxis()->SetTitle("% of negative labels");
      fHistTrNegativeLabels->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistTrNegativeLabels);

      fHistTrZeroLabels = new TH1F("fHistTrZeroLabels","fHistTrZeroLabels", 500, 0, 1);
      fHistTrZeroLabels->GetXaxis()->SetTitle("% of negative labels");
      fHistTrZeroLabels->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistTrZeroLabels);
    }

    TString histname;

    Int_t nlabels = 4;
    if (fParticleLevel)
      nlabels = 1;

    for (Int_t i = 0; i < fNcentBins; i++) {
      for (Int_t j = 0; j < nlabels; j++) {
	histname = Form("fHistTrPhiEtaPt_%d_%d",i,j);
	fHistTrPhiEtaPt[i][j] = new TH3F(histname,histname, 100, -1, 1, 201, 0, TMath::Pi() * 2.01, fNbins, fMinBinPt, fMaxBinPt);
	fHistTrPhiEtaPt[i][j]->GetXaxis()->SetTitle("#eta");
	fHistTrPhiEtaPt[i][j]->GetYaxis()->SetTitle("#phi");
	fHistTrPhiEtaPt[i][j]->GetZaxis()->SetTitle("p_{T} (GeV/c)");
	fOutput->Add(fHistTrPhiEtaPt[i][j]);
      }

      if (!fParticleLevel) {
	if (fIsMC) {
	  histname = Form("fHistTrPhiEtaZeroLab_%d",i);
	  fHistTrPhiEtaZeroLab[i] = new TH2F(histname,histname, 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
	  fHistTrPhiEtaZeroLab[i]->GetXaxis()->SetTitle("#eta");
	  fHistTrPhiEtaZeroLab[i]->GetYaxis()->SetTitle("#phi");
	  fHistTrPhiEtaZeroLab[i]->GetZaxis()->SetTitle("counts");
	  fOutput->Add(fHistTrPhiEtaZeroLab[i]);

	  histname = Form("fHistTrPtZeroLab_%d",i);
	  fHistTrPtZeroLab[i] = new TH1F(histname,histname, fNbins, fMinBinPt, fMaxBinPt);
	  fHistTrPtZeroLab[i]->GetZaxis()->SetTitle("p_{T} (GeV/c)");
	  fHistTrPtZeroLab[i]->GetYaxis()->SetTitle("counts");
	  fOutput->Add(fHistTrPtZeroLab[i]);
	}
	
	histname = Form("fHistTrEmcPhiEta_%d",i);
	fHistTrEmcPhiEta[i] = new TH2F(histname,histname, 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
	fHistTrEmcPhiEta[i]->GetXaxis()->SetTitle("#eta");
	fHistTrEmcPhiEta[i]->GetYaxis()->SetTitle("#phi");
	fOutput->Add(fHistTrEmcPhiEta[i]);
	
	histname = Form("fHistTrEmcPt_%d",i);
	fHistTrEmcPt[i] = new TH1F(histname,histname, fNbins, fMinBinPt, fMaxBinPt);
	fHistTrEmcPt[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	fHistTrEmcPt[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistTrEmcPt[i]);
	
	histname = Form("fHistTrPhiEtaNonProp_%d",i);
	fHistTrPhiEtaNonProp[i] = new TH2F(histname,histname, 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
	fHistTrPhiEtaNonProp[i]->GetXaxis()->SetTitle("#eta");
	fHistTrPhiEtaNonProp[i]->GetYaxis()->SetTitle("#phi");
	fOutput->Add(fHistTrPhiEtaNonProp[i]);

	histname = Form("fHistTrPtNonProp_%d",i);
	fHistTrPtNonProp[i] = new TH1F(histname,histname, fNbins, fMinBinPt, fMaxBinPt);
	fHistTrPtNonProp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	fHistTrPtNonProp[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistTrPtNonProp[i]);
	
	histname = Form("fHistDeltaEtaPt_%d",i);
	fHistDeltaEtaPt[i] = new TH2F(histname,histname, fNbins, fMinBinPt, fMaxBinPt, 50, -0.5, 0.5);
	fHistDeltaEtaPt[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	fHistDeltaEtaPt[i]->GetYaxis()->SetTitle("#delta#eta");
	fOutput->Add(fHistDeltaEtaPt[i]);
	
	histname = Form("fHistDeltaPhiPt_%d",i);
	fHistDeltaPhiPt[i] = new TH2F(histname,histname, fNbins, fMinBinPt, fMaxBinPt, 200, -2, 2);
	fHistDeltaPhiPt[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	fHistDeltaPhiPt[i]->GetYaxis()->SetTitle("#delta#phi");
	fOutput->Add(fHistDeltaPhiPt[i]);
	
	histname = Form("fHistDeltaPtvsPtvsMass_%d",i);
	fHistDeltaPtvsPtvsMass[i] = new TH3F(histname,histname, fNbins, fMinBinPt, fMaxBinPt, fNbins, -fMaxBinPt/2, fMaxBinPt/2, 30, 0, 3);
	fHistDeltaPtvsPtvsMass[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	fHistDeltaPtvsPtvsMass[i]->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
	fHistDeltaPtvsPtvsMass[i]->GetZaxis()->SetTitle("mass (GeV/c^{2})");
	fOutput->Add(fHistDeltaPtvsPtvsMass[i]);
      }
    }
  }

  if (fClusterCollArray.GetEntriesFast()>0) {
    TString histname;

    for (Int_t i = 0; i < fNcentBins; i++) {
      histname = "fHistClusPhiEtaEnergy_";
      histname += i;
      fHistClusPhiEtaEnergy[i] = new TH3F(histname, histname, 100, -1.2, 1.2, 201, 0, TMath::Pi() * 2.01, fNbins, fMinBinPt, fMaxBinPt);
      fHistClusPhiEtaEnergy[i]->GetXaxis()->SetTitle("#eta");
      fHistClusPhiEtaEnergy[i]->GetYaxis()->SetTitle("#phi");
      fHistClusPhiEtaEnergy[i]->GetZaxis()->SetTitle("E_{cluster} (GeV)");
      fOutput->Add(fHistClusPhiEtaEnergy[i]);

      if (fIsEmbedded) {
	histname = "fHistClusMCEnergyFraction_";
	histname += i;
	fHistClusMCEnergyFraction[i] = new TH1F(histname, histname, fNbins, 0, 1.2);
	fHistClusMCEnergyFraction[i]->GetXaxis()->SetTitle("MC fraction");
	fHistClusMCEnergyFraction[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistClusMCEnergyFraction[i]);
      }
    }

    fHistClusTimeEnergy = new TH2F("fHistClusTimeEnergy","Time vs. energy of clusters", fNbins, fMinBinPt, fMaxBinPt, fNbins,  -1e-6, 1e-6);
    fHistClusTimeEnergy->GetXaxis()->SetTitle("E_{cluster} (GeV)");
    fHistClusTimeEnergy->GetYaxis()->SetTitle("Time");
    fOutput->Add(fHistClusTimeEnergy);

    Int_t nbins = fMaxCellsInCluster;
    while (nbins > fNbins) nbins /= 2;
    fHistNCellsEnergy = new TH2F("fHistNCellsEnergy","Number of cells vs. energy of clusters", fNbins, fMinBinPt, fMaxBinPt, nbins, -0.5, fMaxCellsInCluster-0.5);
    fHistNCellsEnergy->GetXaxis()->SetTitle("E_{cluster} (GeV)");
    fHistNCellsEnergy->GetYaxis()->SetTitle("N_{cells}");
    fOutput->Add(fHistNCellsEnergy);

    fHistFcrossEnergy = new TH2F("fHistFcrossEnergy","fHistFcrossEnergy", fNbins, fMinBinPt, fMaxBinPt, 200, -3.5, 1.5);
    fHistFcrossEnergy->GetXaxis()->SetTitle("E_{cluster} (GeV)");
    fHistFcrossEnergy->GetYaxis()->SetTitle("F_{cross}");
    fOutput->Add(fHistFcrossEnergy); 
     
    fHistCellsAbsIdEnergy = new TH2F("fHistCellsAbsIdEnergy","fHistCellsAbsIdEnergy", 11600,0,11599,(Int_t)(fNbins / 2), fMinBinPt, fMaxBinPt / 2);
    fHistCellsAbsIdEnergy->GetXaxis()->SetTitle("cell abs. Id");
    fHistCellsAbsIdEnergy->GetYaxis()->SetTitle("E_{cluster} (GeV)");
    fHistCellsAbsIdEnergy->GetZaxis()->SetTitle("counts");    
    fOutput->Add(fHistCellsAbsIdEnergy);
    
    fHistChVSneCells = new TH2F("fHistChVSneCells","Charged energy vs. neutral (cells) energy", 
				(Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5);
    fHistChVSneCells->GetXaxis()->SetTitle("Energy (GeV)");
    fHistChVSneCells->GetYaxis()->SetTitle("Momentum (GeV/c)");
    fOutput->Add(fHistChVSneCells);
    
    fHistChVSneClus = new TH2F("fHistChVSneClus","Charged energy vs. neutral (clusters) energy", 
			       (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5);
    fHistChVSneClus->GetXaxis()->SetTitle("Energy (GeV)");
    fHistChVSneClus->GetYaxis()->SetTitle("Momentum (GeV/c)");
    fOutput->Add(fHistChVSneClus);
    
    fHistChVSneCorrCells = new TH2F("fHistChVSneCorrCells","Charged energy vs. neutral (corrected cells) energy", 
				    (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, (Int_t)(fNbins * 2.5), fMinBinPt , fMaxBinPt * 2.5);
    fHistChVSneCorrCells->GetXaxis()->SetTitle("Energy (GeV)");
    fHistChVSneCorrCells->GetYaxis()->SetTitle("Momentum (GeV/c)");
    fOutput->Add(fHistChVSneCorrCells);
  }
       
  if (fJetCollArray.GetEntriesFast()>0) {

    TString histname;

    for (Int_t i = 0; i < fNcentBins; i++) {
      histname = "fHistJetsPhiEta_";
      histname += i;
      fHistJetsPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 100, -1.2, 1.2, 201, 0, TMath::Pi() * 2.01);
      fHistJetsPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistJetsPhiEta[i]->GetYaxis()->SetTitle("#phi");
      fHistJetsPhiEta[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetsPhiEta[i]);

      histname = "fHistJetsPtArea_";
      histname += i;
      fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, 50, 0, 1.5);
      fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistJetsPtArea[i]);
    }
  }
  
  Int_t dim = 0;
  TString title[20];
  Int_t nbins[20] = {0};
  Double_t min[20] = {0};
  Double_t max[20] = {0};
  
  if (fForceBeamType != AliAnalysisTaskEmcalDev::kpp) {
    title[dim] = "Centrality %";
    nbins[dim] = 101;
    min[dim] = 0;
    max[dim] = 101;
    dim++;

    if (!fCentMethod2.IsNull()) {
      title[dim] = Form("Centrality %s %%", fCentMethod2.Data());
      nbins[dim] = 101;
      min[dim] = 0;
      max[dim] = 101;
      dim++;
    }

    if (!fCentMethod3.IsNull()) {
      title[dim] = Form("Centrality %s %%", fCentMethod3.Data());
      nbins[dim] = 101;
      min[dim] = 0;
      max[dim] = 101;
      dim++;
    }
    
    if (fDoV0QA==1) {
      title[dim] = "V0A total multiplicity";
      nbins[dim] = 200;
      min[dim] = 0;
      max[dim] = 20000;
      dim++;

      title[dim] = "V0C total multiplicity";
      nbins[dim] = 200;
      min[dim] = 0;
      max[dim] = 20000;
      dim++;
    }
    else if (fDoV0QA==2) {
      title[dim] = "V0A+V0C total multiplicity";
      nbins[dim] = 300;
      min[dim] = 0;
      max[dim] = 30000;
      dim++;
    }

    if (!fRhoName.IsNull()) {
      title[dim] = "#rho (GeV/c)";
      nbins[dim] = fNbins*4;
      min[dim] = 0;
      max[dim] = 400;
      dim++;
    }

    if (fDoEPQA) {
      title[dim] = "#psi_{RP}";
      nbins[dim] = 200;
      min[dim] = -TMath::Pi();
      max[dim] = TMath::Pi();
      dim++;
    }
  }

  if (fParticleCollArray.GetEntriesFast()>0) {
    title[dim] = "No. of tracks";
    nbins[dim] = 3000;
    min[dim] = -0.5;
    max[dim] = 6000-0.5;
    dim++;

    title[dim] = "p_{T,track}^{leading} (GeV/c)";
    nbins[dim] = fNbins;
    min[dim] = fMinBinPt;
    max[dim] = fMaxBinPt;
    dim++;
  }

  if (fClusterCollArray.GetEntriesFast()>0) {
    title[dim] = "No. of clusters";
    nbins[dim] = 2000;
    min[dim] = 0;
    max[dim] = 4000-0.5;
    dim++;

    title[dim] = "E_{cluster}^{leading} (GeV)";
    nbins[dim] = fNbins;
    min[dim] = fMinBinPt;
    max[dim] = fMaxBinPt;
    dim++;
  }

  if (!fCaloCellsName.IsNull()) {
    title[dim] = "No. of cells";
    nbins[dim] = 3000;
    min[dim] = 0;
    max[dim] = 6000-0.5;
    dim++;
  }

  if (fJetCollArray.GetEntriesFast()>0) {
    title[dim] = "No. of jets";
    nbins[dim] = 200;
    min[dim] = 0;
    max[dim] = 200-0.5;
    dim++;

    title[dim] = "p_{T,jet}^{leading} (GeV/c)";
    nbins[dim] = fNbins;
    min[dim] = fMinBinPt;
    max[dim] = fMaxBinPt;
    dim++;
  }

  fHistEventQA = new THnSparseF("fHistEventQA","fHistEventQA",dim,nbins,min,max);
  for (Int_t i = 0; i < dim; i++)
    fHistEventQA->GetAxis(i)->SetTitle(title[i]);
  fOutput->Add(fHistEventQA);

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::ExecOnce()
{
  AliAnalysisTaskEmcalJetDev::ExecOnce();
  
  if (fDoV0QA) {
    fVZERO = InputEvent()->GetVZEROData();
    if (!fVZERO) {
      AliError("AliVVZERO not available");
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAQA::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJetDev::RetrieveEventObjects())
    return kFALSE;

  if (!fCentMethod2.IsNull() || !fCentMethod3.IsNull()) {
    if (fBeamType == kAA || fBeamType == kpA ) {
      AliCentrality *aliCent = InputEvent()->GetCentrality();
      if (aliCent) {
	if (!fCentMethod2.IsNull()) 
	  fCent2 = aliCent->GetCentralityPercentile(fCentMethod2); 
	if (!fCentMethod3.IsNull()) 
	  fCent3 = aliCent->GetCentralityPercentile(fCentMethod3);
      }
    }
  }

  if (fVZERO) {
    fV0ATotMult = AliESDUtils::GetCorrV0A(fVZERO->GetMTotV0A(),fVertex[2]);
    fV0CTotMult = AliESDUtils::GetCorrV0C(fVZERO->GetMTotV0C(),fVertex[2]);
  }

  return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSAQA::FillHistograms()
{
  // Fill histograms.

  Float_t trackSum = 0;
  Float_t clusSum = 0;
  Float_t cellSum = 0;
  Float_t cellCutSum = 0;

  Int_t ntracks = 0;
  Int_t nclusters = 0;
  Int_t ncells = 0;
  Int_t njets = 0;

  Float_t leadingTrack = 0;
  Float_t leadingClus = 0;
  Float_t leadingJet = 0;
    
  if (fTracks) {
    ntracks = DoTrackLoop(trackSum, leadingTrack);
    AliDebug(2,Form("%d tracks found in the event", ntracks));
  } 

  if (fCaloClusters) {
    nclusters = DoClusterLoop(clusSum, leadingClus);
    AliDebug(2,Form("%d clusters found in the event", nclusters));

    fHistChVSneClus->Fill(clusSum, trackSum);
  }
  
  if (fCaloCells) {
    ncells = DoCellLoop(cellSum, cellCutSum);
    AliDebug(2,Form("%d cells found in the event", ncells));
    
    fHistChVSneCells->Fill(cellSum, trackSum);
    fHistChVSneCorrCells->Fill(cellCutSum, trackSum);
  }

  if (fJets) {
    njets = DoJetLoop(leadingJet);
    AliDebug(2,Form("%d jets found in the event", njets));
  }

  FillEventQAHisto(fCent, fCent2, fCent3, fV0ATotMult, fV0CTotMult, fEPV0, fRhoVal, ntracks, nclusters, ncells, njets, leadingTrack, leadingClus, leadingJet);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::FillEventQAHisto(Float_t cent, Float_t cent2, Float_t cent3, Float_t v0a, Float_t v0c, 
					   Float_t ep, Float_t rho, Int_t ntracks, Int_t nclusters, Int_t ncells, Int_t njets, 
					   Float_t maxtrack, Float_t maxcluster, Float_t maxjet)
{
  Double_t contents[20]={0};

  for (Int_t i = 0; i < fHistEventQA->GetNdimensions(); i++) {
    TString title(fHistEventQA->GetAxis(i)->GetTitle());
    if (title=="Centrality %")
      contents[i] = cent;
    else if (title==Form("Centrality %s %%", fCentMethod2.Data()))
      contents[i] = cent2;
    else if (title==Form("Centrality %s %%", fCentMethod3.Data()))
      contents[i] = cent3;
    else if (title=="V0A total multiplicity")
      contents[i] = v0a;
    else if (title=="V0C total multiplicity")
      contents[i] = v0c;
    else if (title=="V0A+V0C total multiplicity")
      contents[i] = v0a+v0c;
    else if (title=="#psi_{RP}")
      contents[i] = ep;
    else if (title=="#rho (GeV/c)")
      contents[i] = rho;
    else if (title=="No. of tracks")
	contents[i] = ntracks;
    else if (title=="No. of clusters")
      contents[i] = nclusters;
    else if (title=="No. of cells")
      contents[i] = ncells;
    else if (title=="No. of jets")
      contents[i] = njets;
    else if (title=="p_{T,track}^{leading} (GeV/c)")
      contents[i] = maxtrack;
    else if (title=="E_{cluster}^{leading} (GeV)")
      contents[i] = maxcluster;
    else if (title=="p_{T,jet}^{leading} (GeV/c)")
      contents[i] = maxjet;
    else 
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }

  fHistEventQA->Fill(contents);
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAQA::DoCellLoop(Float_t &sum, Float_t &sum_cut)
{
  // Do cell loop.

  AliVCaloCells *cells = InputEvent()->GetEMCALCells();

  if (!cells)
    return 0;

  const Int_t ncells = cells->GetNumberOfCells();

  for (Int_t pos = 0; pos < ncells; pos++) {
    Float_t amp   = cells->GetAmplitude(pos);
    Int_t   absId = cells->GetCellNumber(pos);
    fHistCellsAbsIdEnergy->Fill(absId,amp);
    sum += amp;
    if (amp < fCellEnergyCut)
      continue;
    sum_cut += amp;
  } 

  return ncells;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSAQA::GetFcross(AliVCluster *cluster, AliVCaloCells *cells)
{
  Int_t    AbsIdseed  = -1;
  Double_t Eseed      = 0;
  for (Int_t i = 0; i < cluster->GetNCells(); i++) {
    if (cells->GetCellAmplitude(cluster->GetCellAbsId(i)) > AbsIdseed) {
      Eseed     = cells->GetCellAmplitude(cluster->GetCellAbsId(i));
      AbsIdseed = cluster->GetCellAbsId(i);
    }
  }

  if (Eseed < 1e-9)
    return 100;

  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1; 
  fGeom->GetCellIndex(AbsIdseed,imod,iTower,iIphi,iIeta); 
  fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,iphi,ieta);  
  
  //Get close cells index and energy, not in corners
  
  Int_t absID1 = -1;
  Int_t absID2 = -1;
  
  if (iphi < AliEMCALGeoParams::fgkEMCALRows-1) absID1 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  if (iphi > 0)                                 absID2 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);
  
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
    if (ieta < AliEMCALGeoParams::fgkEMCALCols-1) 
      absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
    if (ieta > 0)                                 
      absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1); 
  }
  
  Double_t  ecell1 = cells->GetCellAmplitude(absID1);
  Double_t  ecell2 = cells->GetCellAmplitude(absID2);
  Double_t  ecell3 = cells->GetCellAmplitude(absID3);
  Double_t  ecell4 = cells->GetCellAmplitude(absID4);

  Double_t Ecross = ecell1 + ecell2 + ecell3 + ecell4;
  
  Double_t Fcross = 1 - Ecross/Eseed;

  return Fcross;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAQA::DoClusterLoop(Float_t &sum, Float_t &leading)
{
  // Do cluster loop.

  if (!fCaloClusters)
    return 0;

  Int_t nAccClusters = 0;

  AliVCaloCells *cells = InputEvent()->GetEMCALCells();

  sum = 0;
  leading = 0;

  // Cluster loop
  const Int_t nclusters = fCaloClusters->GetEntriesFast();

  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = static_cast<AliVCluster*>(fCaloClusters->At(iClusters));
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  

    if (!AcceptCluster(cluster))
      continue;

    sum += cluster->E();

    if (leading < cluster->E()) leading = cluster->E();

    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);

    fHistClusPhiEtaEnergy[fCentBin]->Fill(nPart.Eta(), nPart.Phi(), cluster->E());
    fHistNCellsEnergy->Fill(cluster->E(), cluster->GetNCells());

    fHistClusTimeEnergy->Fill(cluster->E(), cluster->GetTOF());

    if (cells)
      fHistFcrossEnergy->Fill(cluster->E(), GetFcross(cluster, cells));

    if (fHistClusMCEnergyFraction[fCentBin])
      fHistClusMCEnergyFraction[fCentBin]->Fill(cluster->GetMCEnergyFraction());

    nAccClusters++;
  }

  return nAccClusters;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAQA::DoTrackLoop(Float_t &sum, Float_t &leading)
{
  // Do track loop.

  if (!fTracks)
    return 0;

  Int_t nAccTracks = 0;

  sum = 0;
  leading = 0;

  const Int_t ntracks = fTracks->GetEntriesFast();
  Int_t neg = 0;
  Int_t zero = 0;

  for (Int_t i = 0; i < ntracks; i++) {

    AliVParticle* track = static_cast<AliVParticle*>(fTracks->At(i)); // pointer to reconstructed to track  

    if (!track) {
      AliError(Form("Could not retrieve track %d",i)); 
      continue; 
    }

    if (!AcceptTrack(track)) 
      continue;

    nAccTracks++;

    sum += track->P();

    if (leading < track->Pt()) leading = track->P();

    if (fParticleLevel) {
      fHistTrPhiEtaPt[fCentBin][0]->Fill(track->Eta(), track->Phi(), track->Pt());
    }
    else {
      fHistTrPhiEtaPt[fCentBin][3]->Fill(track->Eta(), track->Phi(), track->Pt());
      if (track->GetLabel() == 0) {
	zero++;
	if (fHistTrPhiEtaZeroLab[fCentBin]) {
	  fHistTrPhiEtaZeroLab[fCentBin]->Fill(track->Eta(), track->Phi());
	  fHistTrPtZeroLab[fCentBin]->Fill(track->Pt());
	}
      }

      if (track->GetLabel() < 0)
	neg++;

      Int_t type = 0;

      AliPicoTrack* ptrack = dynamic_cast<AliPicoTrack*>(track);
      if (ptrack)
	type = ptrack->GetTrackType();

      if (type >= 0 && type < 3)
	fHistTrPhiEtaPt[fCentBin][type]->Fill(track->Eta(), track->Phi(), track->Pt());
      else
	AliDebug(2,Form("%s: track type %d not recognized!", GetName(), type));
    }

    AliVTrack* vtrack = dynamic_cast<AliVTrack*>(track); 

    if (!vtrack)
      continue;

    if ((vtrack->GetTrackEtaOnEMCal() == -999 || vtrack->GetTrackPhiOnEMCal() == -999) && fHistTrPhiEtaNonProp[fCentBin]) {
      fHistTrPhiEtaNonProp[fCentBin]->Fill(vtrack->Eta(), vtrack->Phi());
      fHistTrPtNonProp[fCentBin]->Fill(vtrack->Pt());
    }

    if (fHistTrEmcPhiEta[fCentBin])
      fHistTrEmcPhiEta[fCentBin]->Fill(vtrack->GetTrackEtaOnEMCal(), vtrack->GetTrackPhiOnEMCal());   
    if (fHistTrEmcPt[fCentBin])
      fHistTrEmcPt[fCentBin]->Fill(vtrack->GetTrackPtOnEMCal());   
    if (fHistDeltaEtaPt[fCentBin])
      fHistDeltaEtaPt[fCentBin]->Fill(vtrack->Pt(), vtrack->Eta() - vtrack->GetTrackEtaOnEMCal());
    if (fHistDeltaPhiPt[fCentBin])
      fHistDeltaPhiPt[fCentBin]->Fill(vtrack->Pt(), vtrack->Phi() - vtrack->GetTrackPhiOnEMCal());
    if (fHistDeltaPtvsPtvsMass[fCentBin])
      fHistDeltaPtvsPtvsMass[fCentBin]->Fill(vtrack->Pt(), vtrack->Pt() - vtrack->GetTrackPtOnEMCal(), vtrack->M());
  }

  if (fHistTrNegativeLabels)
    fHistTrNegativeLabels->Fill(1. * neg / ntracks);

  if (fHistTrZeroLabels)
    fHistTrZeroLabels->Fill(1. * zero / ntracks);

  return nAccTracks;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAQA::DoJetLoop(Float_t &leading)
{
  // Do jet loop.

  if (!fJets)
    return 0;

  Int_t nAccJets = 0;

  leading = 0;

  Int_t njets = fJets->GetEntriesFast();

  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    if (leading < jet->Pt()) leading = jet->Pt();

    nAccJets++;

    fHistJetsPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());
    fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
  }

  return nAccJets;
}
