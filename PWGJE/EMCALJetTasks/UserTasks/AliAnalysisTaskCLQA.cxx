// $Id: AliAnalysisTaskCLQA.cxx 60694 2013-02-04 15:35:56Z morsch $
//
// Constantin's Task
//
// Author: C.Loizides

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeoParams.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalJet.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliTrackerBase.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliAnalysisTaskCLQA.h"

ClassImp(AliAnalysisTaskCLQA)

//________________________________________________________________________
AliAnalysisTaskCLQA::AliAnalysisTaskCLQA() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskCLQA", kTRUE)
{
  // Default constructor.
}

//________________________________________________________________________
AliAnalysisTaskCLQA::AliAnalysisTaskCLQA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE)
{
  // Standard constructor.
}

//________________________________________________________________________
AliAnalysisTaskCLQA::~AliAnalysisTaskCLQA()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskCLQA::UserCreateOutputObjects()
{
  // Create histograms

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

#if 0
  if (!fTracksName.IsNull()) {
    fHistTracksCent = new TH2F("fHistTracksCent","Tracks vs. centrality", 100, 0, 100, fNbins, 0, 4000);
    fHistTracksCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistTracksCent->GetYaxis()->SetTitle("No. of tracks");
    fOutput->Add(fHistTracksCent);

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
    }

    fHistTrEmcPhiEta = new TH2F("fHistTrEmcPhiEta","Phi-Eta emcal distribution of tracks", 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
    fHistTrEmcPhiEta->GetXaxis()->SetTitle("#eta");
    fHistTrEmcPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistTrEmcPhiEta);

    fHistTrPhiEtaNonProp = new TH2F("fHistTrPhiEtaNonProp","fHistTrPhiEtaNonProp", 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
    fHistTrPhiEtaNonProp->GetXaxis()->SetTitle("#eta");
    fHistTrPhiEtaNonProp->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistTrPhiEtaNonProp);
    
    fHistDeltaEtaPt = new TH2F("fHistDeltaEtaPt","fHistDeltaEtaPt", fNbins, fMinBinPt, fMaxBinPt, 80, -0.5, 0.5);
    fHistDeltaEtaPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistDeltaEtaPt->GetYaxis()->SetTitle("#delta#eta");
    fOutput->Add(fHistDeltaEtaPt);
    
    fHistDeltaPhiPt = new TH2F("fHistDeltaPhiPt","fHistDeltaPhiPt", fNbins, fMinBinPt, fMaxBinPt, 256, -1.6, 4.8);
    fHistDeltaPhiPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistDeltaPhiPt->GetYaxis()->SetTitle("#delta#phi");
    fOutput->Add(fHistDeltaPhiPt);
  }

  if (!fCaloName.IsNull()) {
    fHistClusCent = new TH2F("fHistClusCent","Clusters vs. centrality", 100, 0, 100, fNbins, 0, 2000);
    fHistClusCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistClusCent->GetYaxis()->SetTitle("No. of clusters");
    fOutput->Add(fHistClusCent);

    fHistClusTracks = new TH2F("fHistClusTracks","Clusters vs. tracks", fNbins, 0, 4000, fNbins, 0, 2000);
    fHistClusTracks->GetXaxis()->SetTitle("No. of tracks");
    fHistClusTracks->GetYaxis()->SetTitle("No. of clusters");
    fOutput->Add(fHistClusTracks);

    fHistCellsCent = new TH2F("fHistCellsCent","Cells vs. centrality", 100, 0, 100, fNbins, 0, 6000);
    fHistCellsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistCellsCent->GetYaxis()->SetTitle("No. of EMCal cells");
    fOutput->Add(fHistCellsCent);

    fHistCellsTracks = new TH2F("fHistCellsTracks","Cells vs. tracks", fNbins, 0, 4000, fNbins, 0, 6000);
    fHistCellsTracks->GetXaxis()->SetTitle("No. of tracks");
    fHistCellsTracks->GetYaxis()->SetTitle("No. of EMCal cells");
    fOutput->Add(fHistCellsTracks);

    TString histname;

    for (Int_t i = 0; i < fNcentBins; i++) {
      histname = "fHistClusPhiEtaEnergy_";
      histname += i;
      fHistClusPhiEtaEnergy[i] = new TH3F(histname, histname, 100, -1.2, 1.2, 201, 0, TMath::Pi() * 2.01, fNbins, fMinBinPt, fMaxBinPt);
      fHistClusPhiEtaEnergy[i]->GetXaxis()->SetTitle("#eta");
      fHistClusPhiEtaEnergy[i]->GetYaxis()->SetTitle("#phi");
      fHistClusPhiEtaEnergy[i]->GetZaxis()->SetTitle("E_{cluster} (GeV)");
      fOutput->Add(fHistClusPhiEtaEnergy[i]);
    }

    fHistClusTimeEnergy = new TH2F("fHistClusTimeEnergy","Time vs. energy of clusters", fNbins, fMinBinPt, fMaxBinPt, fNbins,  -1e-6, 1e-6);
    fHistClusTimeEnergy->GetXaxis()->SetTitle("E_{cluster} (GeV)");
    fHistClusTimeEnergy->GetYaxis()->SetTitle("Time");
    fOutput->Add(fHistClusTimeEnergy);

    fHistNCellsEnergy = new TH2F("fHistNCellsEnergy","Number of cells vs. energy of clusters", fNbins, fMinBinPt, fMaxBinPt, 30, 0, 30);
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
       
  if (!fJetsName.IsNull()) {
    fHistJetsCent = new TH2F("fHistJetsCent","Jets vs. centrality", 100, 0, 100, 150, -0.5, 149.5);
    fHistJetsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetsCent->GetYaxis()->SetTitle("No. of jets");
    fOutput->Add(fHistJetsCent);
    
    fHistJetsParts = new TH2F("fHistJetsParts","Jets vs. centrality", fNbins, 0, 6000, 150, -0.5, 149.5);
    fHistJetsParts->GetXaxis()->SetTitle("No. of particles");
    fHistJetsParts->GetYaxis()->SetTitle("No. of jets");
    fOutput->Add(fHistJetsParts);

    TString histname;

    for (Int_t i = 0; i < fNcentBins; i++) {
      histname = "fHistJetsPhiEtaPt_";
      histname += i;
      fHistJetsPhiEtaPt[i] = new TH3F(histname.Data(), histname.Data(), 100, -1.2, 1.2, 201, 0, TMath::Pi() * 2.01,  (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5);
      fHistJetsPhiEtaPt[i]->GetXaxis()->SetTitle("#eta");
      fHistJetsPhiEtaPt[i]->GetYaxis()->SetTitle("#phi");
      fHistJetsPhiEtaPt[i]->GetZaxis()->SetTitle("p_{T} (GeV/c)");
      fOutput->Add(fHistJetsPhiEtaPt[i]);

      histname = "fHistJetsPtArea_";
      histname += i;
      fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, 30, 0, fJetRadius * fJetRadius * TMath::Pi() * 3);
      fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistJetsPtArea[i]);
    }
  }
#endif

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskCLQA::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskCLQA::FillHistograms()
{
  // Fill histograms.

#if 0
  Float_t clusSum = 0;
  Float_t trackSum = 0;

  if (fTracks) {
    trackSum = DoTrackLoop();

    fHistTracksCent->Fill(fCent, fNtracks);
  } 

  if (fCaloClusters) {
    clusSum = DoClusterLoop();

    fHistClusCent->Fill(fCent, fNclusters);
    fHistClusTracks->Fill(fNtracks, fNclusters);

    Float_t cellSum = 0, cellCutSum = 0;
    
    Int_t ncells = DoCellLoop(cellSum, cellCutSum);

    if (fTracks)
      fHistCellsTracks->Fill(fNtracks, ncells);

    fHistCellsCent->Fill(fCent, ncells);
    
    fHistChVSneCells->Fill(cellSum, trackSum);
    fHistChVSneClus->Fill(clusSum, trackSum);
    fHistChVSneCorrCells->Fill(cellCutSum, trackSum);
  }

  if (fJets) {
    DoJetLoop();
    
    fHistJetsCent->Fill(fCent, fNjets);
    fHistJetsParts->Fill(fNtracks + fNclusters, fNjets);
  }
#endif

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskCLQA::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
