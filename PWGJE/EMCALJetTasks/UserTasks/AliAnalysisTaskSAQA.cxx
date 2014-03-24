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
  AliAnalysisTaskEmcalJet("AliAnalysisTaskSAQA", kTRUE),
  fCellEnergyCut(0.05),
  fParticleLevel(kFALSE),
  fIsMC(kFALSE),
  fCentMethod2(""),
  fCentMethod3(""),
  fDoV0QA(0),
  fDoEPQA(0),
  fDoLeadingObjectPosition(0),
  fMaxCellsInCluster(30),
  fCent2(0),
  fCent3(0),
  fVZERO(0),
  fV0ATotMult(0),
  fV0CTotMult(0),
  fHistEventQA(0)
{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < 4; j++) fHistTrPhiEtaPt[i][j] = 0;

    fHistTrNegativeLabels[i] = 0;
    fHistTrZeroLabels[i] = 0;
    fHistTrPhiEtaZeroLab[i] = 0;
    fHistTrPtZeroLab[i] = 0;
    fHistTrEmcPhiEta[i] = 0;
    fHistTrEmcPt[i] = 0;
    fHistTrPhiEtaNonProp[i] = 0;
    fHistTrPtNonProp[i] = 0;
    fHistDeltaEtaPt[i] = 0;
    fHistDeltaPhiPt[i] = 0;
    fHistDeltaPtvsPt[i] = 0;

    fHistClusPhiEtaEnergy[i] = 0;
    fHistClusDeltaPhiEPEnergy[i] = 0;
    fHistNCellsEnergy[i] = 0;
    fHistFcrossEnergy[i] = 0;
    fHistClusTimeEnergy[i] = 0;
    fHistClusMCEnergyFraction[i] = 0;

    fHistCellsAbsIdEnergy[i] = 0;

    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fCellEnergyCut(0.05),
  fParticleLevel(kFALSE),
  fIsMC(kFALSE),
  fCentMethod2(""),
  fCentMethod3(""),
  fDoV0QA(0),
  fDoEPQA(0),
  fDoLeadingObjectPosition(0),
  fMaxCellsInCluster(30),
  fCent2(0),
  fCent3(0),
  fVZERO(0),
  fV0ATotMult(0),
  fV0CTotMult(0),
  fHistEventQA(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < 4; j++) fHistTrPhiEtaPt[i][j] = 0;

    fHistTrNegativeLabels[i] = 0;
    fHistTrZeroLabels[i] = 0;
    fHistTrPhiEtaZeroLab[i] = 0;
    fHistTrPtZeroLab[i] = 0;
    fHistTrEmcPhiEta[i] = 0;
    fHistTrEmcPt[i] = 0;
    fHistTrPhiEtaNonProp[i] = 0;
    fHistTrPtNonProp[i] = 0;
    fHistDeltaEtaPt[i] = 0;
    fHistDeltaPhiPt[i] = 0;
    fHistDeltaPtvsPt[i] = 0;

    fHistClusPhiEtaEnergy[i] = 0;
    fHistClusDeltaPhiEPEnergy[i] = 0;
    fHistNCellsEnergy[i] = 0;
    fHistFcrossEnergy[i] = 0;
    fHistClusTimeEnergy[i] = 0;
    fHistClusMCEnergyFraction[i] = 0;

    fHistCellsAbsIdEnergy[i] = 0;

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

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  TString histname;

  if (fParticleCollArray.GetEntriesFast()>0) {
    if (!fParticleLevel && fIsMC) {
      for (Int_t i = 0; i < fNcentBins; i++) {
	histname = Form("fHistTrNegativeLabels_%d",i);
	fHistTrNegativeLabels[i] = new TH1F(histname,histname, 500, 0, 1);
	fHistTrNegativeLabels[i]->GetXaxis()->SetTitle("% of negative labels");
	fHistTrNegativeLabels[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistTrNegativeLabels[i]);

	histname = Form("fHistTrZeroLabels_%d",i);	
	fHistTrZeroLabels[i] = new TH1F(histname,histname, 500, 0, 1);
	fHistTrZeroLabels[i]->GetXaxis()->SetTitle("% of negative labels");
	fHistTrZeroLabels[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistTrZeroLabels[i]);
      }
    }

    Int_t nlabels = 4;
    if (fParticleLevel)
      nlabels = 1;

    for (Int_t i = 0; i < fNcentBins; i++) {
      for (Int_t j = 0; j < nlabels; j++) {
	histname = Form("fHistTrPhiEtaPt_%d_%d",i,j);
	fHistTrPhiEtaPt[i][j] = new TH3F(histname,histname, 100, -1, 1, 101, 0, TMath::Pi() * 2.02, fNbins, fMinBinPt, fMaxBinPt);
	fHistTrPhiEtaPt[i][j]->GetXaxis()->SetTitle("#eta");
	fHistTrPhiEtaPt[i][j]->GetYaxis()->SetTitle("#phi");
	fHistTrPhiEtaPt[i][j]->GetZaxis()->SetTitle("p_{T} (GeV/c)");
	fOutput->Add(fHistTrPhiEtaPt[i][j]);
      }

      if (!fParticleLevel) {
	if (fIsMC) {
	  histname = Form("fHistTrPhiEtaZeroLab_%d",i);
	  fHistTrPhiEtaZeroLab[i] = new TH2F(histname,histname, 100, -1, 1, 101, 0, TMath::Pi() * 2.02);
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
	fHistTrEmcPhiEta[i] = new TH2F(histname,histname, 100, -1, 1, 101, 0, TMath::Pi() * 2.02);
	fHistTrEmcPhiEta[i]->GetXaxis()->SetTitle("#eta");
	fHistTrEmcPhiEta[i]->GetYaxis()->SetTitle("#phi");
	fOutput->Add(fHistTrEmcPhiEta[i]);
	
	histname = Form("fHistTrEmcPt_%d",i);
	fHistTrEmcPt[i] = new TH1F(histname,histname, fNbins, fMinBinPt, fMaxBinPt);
	fHistTrEmcPt[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	fHistTrEmcPt[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistTrEmcPt[i]);
	
	histname = Form("fHistTrPhiEtaNonProp_%d",i);
	fHistTrPhiEtaNonProp[i] = new TH2F(histname,histname, 100, -1, 1, 101, 0, TMath::Pi() * 2.02);
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
	
	histname = Form("fHistDeltaPtvsPt_%d",i);
	fHistDeltaPtvsPt[i] = new TH2F(histname,histname, fNbins, fMinBinPt, fMaxBinPt, fNbins, -fMaxBinPt/2, fMaxBinPt/2);
	fHistDeltaPtvsPt[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	fHistDeltaPtvsPt[i]->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
	fHistDeltaPtvsPt[i]->GetZaxis()->SetTitle("counts");
	fOutput->Add(fHistDeltaPtvsPt[i]);
      }
    }
  }

  if (fClusterCollArray.GetEntriesFast()>0) {
    for (Int_t i = 0; i < fNcentBins; i++) {
      histname = "fHistClusPhiEtaEnergy_";
      histname += i;
      fHistClusPhiEtaEnergy[i] = new TH3F(histname, histname, 100, -1, 1, 101, 0, TMath::Pi() * 2.02, fNbins, fMinBinPt, fMaxBinPt);
      fHistClusPhiEtaEnergy[i]->GetXaxis()->SetTitle("#eta");
      fHistClusPhiEtaEnergy[i]->GetYaxis()->SetTitle("#phi");
      fHistClusPhiEtaEnergy[i]->GetZaxis()->SetTitle("E_{cluster} (GeV)");
      fOutput->Add(fHistClusPhiEtaEnergy[i]);

      histname = "fHistClusDeltaPhiEPEnergy_";
      histname += i;
      fHistClusDeltaPhiEPEnergy[i] = new TH2F(histname, histname, fNbins, fMinBinPt, fMaxBinPt, 100, 0, TMath::Pi());
      fHistClusDeltaPhiEPEnergy[i]->GetXaxis()->SetTitle("E_{cluster} (GeV)");
      fHistClusDeltaPhiEPEnergy[i]->GetYaxis()->SetTitle("#phi_{cluster} - #psi_{RP}");
      fOutput->Add(fHistClusDeltaPhiEPEnergy[i]);

      if (fIsEmbedded) {
	histname = "fHistClusMCEnergyFraction_";
	histname += i;
	fHistClusMCEnergyFraction[i] = new TH1F(histname, histname, fNbins, 0, 1.2);
	fHistClusMCEnergyFraction[i]->GetXaxis()->SetTitle("MC fraction");
	fHistClusMCEnergyFraction[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistClusMCEnergyFraction[i]);
      }

      histname = "fHistClusTimeEnergy_";
      histname += i;
      fHistClusTimeEnergy[i] = new TH2F(histname,histname, fNbins, fMinBinPt, fMaxBinPt, fNbins,  -1e-6, 1e-6);
      fHistClusTimeEnergy[i]->GetXaxis()->SetTitle("E_{cluster} (GeV)");
      fHistClusTimeEnergy[i]->GetYaxis()->SetTitle("Time");
      fOutput->Add(fHistClusTimeEnergy[i]);

      Int_t nbins = fMaxCellsInCluster;
      while (nbins > fNbins) nbins /= 2;
      histname = "fHistNCellsEnergy_";
      histname += i;
      fHistNCellsEnergy[i] = new TH2F(histname,histname, fNbins, fMinBinPt, fMaxBinPt, nbins, -0.5, fMaxCellsInCluster-0.5);
      fHistNCellsEnergy[i]->GetXaxis()->SetTitle("E_{cluster} (GeV)");
      fHistNCellsEnergy[i]->GetYaxis()->SetTitle("N_{cells}");
      fOutput->Add(fHistNCellsEnergy[i]);

      histname = "fHistFcrossEnergy_";
      histname += i;
      fHistFcrossEnergy[i] = new TH2F(histname,histname, fNbins, fMinBinPt, fMaxBinPt, 200, -3.5, 1.5);
      fHistFcrossEnergy[i]->GetXaxis()->SetTitle("E_{cluster} (GeV)");
      fHistFcrossEnergy[i]->GetYaxis()->SetTitle("F_{cross}");
      fOutput->Add(fHistFcrossEnergy[i]); 

      histname = "fHistCellsAbsIdEnergy_";
      histname += i;
      fHistCellsAbsIdEnergy[i] = new TH2F(histname,histname, 11600,0,11599,(Int_t)(fNbins / 2), fMinBinPt, fMaxBinPt / 2);
      fHistCellsAbsIdEnergy[i]->GetXaxis()->SetTitle("cell abs. Id");
      fHistCellsAbsIdEnergy[i]->GetYaxis()->SetTitle("E_{cluster} (GeV)");
      fHistCellsAbsIdEnergy[i]->GetZaxis()->SetTitle("counts");    
      fOutput->Add(fHistCellsAbsIdEnergy[i]);
    }
  }
       
  if (fJetCollArray.GetEntriesFast()>0) {
    for (Int_t i = 0; i < fNcentBins; i++) {
      histname = "fHistJetsPhiEta_";
      histname += i;
      fHistJetsPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 100, -1, 1, 101, 0, TMath::Pi() * 2.02);
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
  
  if (fForceBeamType != AliAnalysisTaskEmcal::kpp) {
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

    if (fDoLeadingObjectPosition) {
      title[dim] = "#eta_{track}^{leading}";
      nbins[dim] = 100;
      min[dim] = -1;
      max[dim] = 1;
      dim++;

      title[dim] = "#phi_{track}^{leading}";
      nbins[dim] = 101;
      min[dim] = 0;
      max[dim] = TMath::Pi() * 2.02;
      dim++;
    }
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

    if (fDoLeadingObjectPosition) {
      title[dim] = "#eta_{cluster}^{leading}";
      nbins[dim] = 100;
      min[dim] = -1;
      max[dim] = 1;
      dim++;

      title[dim] = "#phi_{cluster}^{leading}";
      nbins[dim] = 101;
      min[dim] = 0;
      max[dim] = TMath::Pi() * 2.02;
      dim++;
    }
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

    if (fDoLeadingObjectPosition) {
      title[dim] = "#eta_{jet}^{leading}";
      nbins[dim] = 100;
      min[dim] = -1;
      max[dim] = 1;
      dim++;

      title[dim] = "#phi_{jet}^{leading}";
      nbins[dim] = 101;
      min[dim] = 0;
      max[dim] = TMath::Pi() * 2.02;
      dim++;
    }
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
  AliAnalysisTaskEmcalJet::ExecOnce();
  
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

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
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

  Int_t ntracks = 0;
  Int_t nclusters = 0;
  Int_t ncells = 0;
  Int_t njets = 0;

  Float_t leadingClusE = 0;
  Float_t leadingClusEta = 0;
  Float_t leadingClusPhi = 0;

  Float_t leadingTrackPt = 0;
  Float_t leadingTrackEta = 0;
  Float_t leadingTrackPhi = 0;

  Float_t leadingJetPt = 0;
  Float_t leadingJetEta = 0;
  Float_t leadingJetPhi = 0;
    
  if (fTracks) {
    AliVParticle *leadingTrack = 0;

    ntracks = DoTrackLoop(trackSum, leadingTrack);
    AliDebug(2,Form("%d tracks found in the event", ntracks));

    if (leadingTrack) {
      leadingTrackPt = leadingTrack->Pt();
      leadingTrackEta = leadingTrack->Eta();
      leadingTrackPhi = leadingTrack->Phi();
    }
  } 

  if (fCaloClusters) {
    AliVCluster  *leadingClus = 0;

    nclusters = DoClusterLoop(clusSum, leadingClus);
    AliDebug(2,Form("%d clusters found in the event", nclusters));

    if (leadingClus) {
      TLorentzVector leadingClusVect;
      leadingClus->GetMomentum(leadingClusVect, fVertex);
      leadingClusE = leadingClus->E();
      leadingClusEta = leadingClusVect.Eta();
      leadingClusPhi = leadingClusVect.Phi();
    }
  }
  
  if (fCaloCells) {
    ncells = DoCellLoop(cellSum);
    AliDebug(2,Form("%d cells found in the event", ncells));
  }

  if (fJets) {
    AliEmcalJet  *leadingJet = 0;

    njets = DoJetLoop(leadingJet);
    AliDebug(2,Form("%d jets found in the event", njets));

    if (leadingJet) {
      leadingJetPt = leadingJet->Pt();
      leadingJetEta = leadingJet->Eta();
      leadingJetPhi = leadingJet->Phi();
    }
  }

  FillEventQAHisto(fCent, fCent2, fCent3, fV0ATotMult, fV0CTotMult, fEPV0, fRhoVal, 
		   ntracks, nclusters, ncells, njets, 
		   leadingTrackPt, leadingTrackEta, leadingTrackPhi, 
		   leadingClusE, leadingClusEta, leadingClusPhi, 
		   leadingJetPt, leadingJetEta, leadingJetPhi);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::FillEventQAHisto(Float_t cent, Float_t cent2, Float_t cent3, Float_t v0a, Float_t v0c, 
					   Float_t ep, Float_t rho, Int_t ntracks, Int_t nclusters, Int_t ncells, Int_t njets, 
					   Float_t maxTrackPt, Float_t maxTrackEta, Float_t maxTrackPhi,
					   Float_t maxClusterE, Float_t maxClusterEta, Float_t maxClusterPhi,
					   Float_t maxJetPt, Float_t maxJetEta, Float_t maxJetPhi)
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
      contents[i] = maxTrackPt;
    else if (title=="#eta_{track}^{leading}")
      contents[i] = maxTrackEta;
    else if (title=="#phi_{track}^{leading}")
      contents[i] = maxTrackPhi;
    else if (title=="E_{cluster}^{leading} (GeV)")
      contents[i] = maxClusterE;
    else if (title=="#eta_{cluster}^{leading}")
      contents[i] = maxClusterEta;
    else if (title=="#phi_{cluster}^{leading}")
      contents[i] = maxClusterPhi;
    else if (title=="p_{T,jet}^{leading} (GeV/c)")
      contents[i] = maxJetPt;
    else if (title=="#eta_{jet}^{leading}")
      contents[i] = maxJetEta;
    else if (title=="#phi_{jet}^{leading}")
      contents[i] = maxJetPhi;
    else 
      AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }

  fHistEventQA->Fill(contents);
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAQA::DoCellLoop(Float_t &sum)
{
  // Do cell loop.

  AliVCaloCells *cells = InputEvent()->GetEMCALCells();

  if (!cells)
    return 0;

  const Int_t ncells = cells->GetNumberOfCells();
  Int_t nAccCells = 0;

  for (Int_t pos = 0; pos < ncells; pos++) {
    Float_t amp   = cells->GetAmplitude(pos);
    Int_t   absId = cells->GetCellNumber(pos);

    if (amp < fCellEnergyCut)
      continue;

    fHistCellsAbsIdEnergy[fCentBin]->Fill(absId,amp);
    nAccCells++;
    sum += amp;
  } 

  return nAccCells;
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
Int_t AliAnalysisTaskSAQA::DoClusterLoop(Float_t &sum, AliVCluster* &leading)
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

    if (!leading || leading->E() < cluster->E()) leading = cluster;

    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);

    fHistClusPhiEtaEnergy[fCentBin]->Fill(nPart.Eta(), nPart.Phi(), cluster->E());

    Double_t ep = nPart.Phi() - fEPV0;
    while (ep < 0) ep += TMath::Pi();
    while (ep >= TMath::Pi()) ep -= TMath::Pi();
    fHistClusDeltaPhiEPEnergy[fCentBin]->Fill(cluster->E(), ep);

    fHistNCellsEnergy[fCentBin]->Fill(cluster->E(), cluster->GetNCells());

    fHistClusTimeEnergy[fCentBin]->Fill(cluster->E(), cluster->GetTOF());

    if (cells) fHistFcrossEnergy[fCentBin]->Fill(cluster->E(), GetFcross(cluster, cells));

    if (fHistClusMCEnergyFraction[fCentBin])
      fHistClusMCEnergyFraction[fCentBin]->Fill(cluster->GetMCEnergyFraction());

    nAccClusters++;
  }

  return nAccClusters;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAQA::DoTrackLoop(Float_t &sum, AliVParticle* &leading)
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

    if (!leading || leading->Pt() < track->Pt()) leading = track;

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
    if (fHistDeltaPtvsPt[fCentBin])
      fHistDeltaPtvsPt[fCentBin]->Fill(vtrack->Pt(), vtrack->Pt() - vtrack->GetTrackPtOnEMCal());
  }

  if (fHistTrNegativeLabels[fCentBin])
    fHistTrNegativeLabels[fCentBin]->Fill(1. * neg / ntracks);

  if (fHistTrZeroLabels[fCentBin])
    fHistTrZeroLabels[fCentBin]->Fill(1. * zero / ntracks);

  return nAccTracks;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAQA::DoJetLoop(AliEmcalJet* &leading)
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

    if (!leading || leading->Pt() < jet->Pt()) leading = jet;

    nAccJets++;

    fHistJetsPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());
    fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
  }

  return nAccJets;
}
