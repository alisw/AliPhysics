// $Id$
//
// General QA task.
//
// Author: S.Aiola

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
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

#include "AliAnalysisTaskSAQA.h"

ClassImp(AliAnalysisTaskSAQA)

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskSAQA", kTRUE),
  fCellEnergyCut(0.1),
  fDoTrigger(kFALSE),
  fRepropagateTracks(kFALSE),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTrgClusters(0),
  fNclusters(0),
  fNtracks(0),
  fNjets(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistTracksCent(0),
  fHistClusCent(0),
  fHistJetsCent(0),
  fHistClusTracks(0),
  fHistJetsParts(0),
  fHistCellsCent(0),
  fHistCellsTracks(0),
  fHistMaxL1FastORCent(0),
  fHistMaxL1ClusCent(0),
  fHistMaxL1ThrCent(0),
  fHistTracksPt(0),
  fHistTrPhiEta(0),
  fHistTrEmcPhiEta(0),
  fHistTrPhiEtaNonProp(0),
  fHistDeltaEtaPt(0),
  fHistDeltaPhiPt(0),
  fHistDeltaEtaNewProp(0),
  fHistDeltaPhiNewProp(0),
  fHistClusPhiEtaEnergy(0),
  fHistNCellsEnergy(0),
  fHistClusTimeEnergy(0),
  fHistCellsEnergy(0),
  fHistChVSneCells(0),
  fHistChVSneClus(0),
  fHistChVSneCorrCells(0)
{
  // Default constructor.

  for (Int_t i = 0; i < 5; i++) {
    fHistTrackPhi[i] = 0;
    fHistTrackEta[i] = 0;
  }

  for (Int_t i = 0; i < 6; i++) {
    fHistTrackPhiPt[i] = 0;
    fHistTrackEtaPt[i] = 0;
  }

  for (Int_t i = 0; i < 4; i++) {
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtNonBias[i] = 0;
    fHistJetsPtTrack[i] = 0;
    fHistJetsPtClus[i] = 0;
    fHistJetsPt[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtAreaNonBias[i] = 0;
  }
}

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fCellEnergyCut(0.1),
  fDoTrigger(kFALSE),
  fRepropagateTracks(kFALSE),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTrgClusters(0),
  fNclusters(0),
  fNtracks(0),
  fNjets(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistTracksCent(0),
  fHistClusCent(0),
  fHistJetsCent(0),
  fHistClusTracks(0),
  fHistJetsParts(0),
  fHistCellsCent(0),
  fHistCellsTracks(0),
  fHistMaxL1FastORCent(0),
  fHistMaxL1ClusCent(0),
  fHistMaxL1ThrCent(0),
  fHistTracksPt(0),
  fHistTrPhiEta(0),
  fHistTrEmcPhiEta(0),
  fHistTrPhiEtaNonProp(0),
  fHistDeltaEtaPt(0),
  fHistDeltaPhiPt(0),
  fHistDeltaEtaNewProp(0),
  fHistDeltaPhiNewProp(0),
  fHistClusPhiEtaEnergy(0),
  fHistNCellsEnergy(0),
  fHistClusTimeEnergy(0),
  fHistCellsEnergy(0),
  fHistChVSneCells(0),
  fHistChVSneClus(0),
  fHistChVSneCorrCells(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 5; i++) {
    fHistTrackPhi[i] = 0;
    fHistTrackEta[i] = 0;
  }

  for (Int_t i = 0; i < 6; i++) {
    fHistTrackPhiPt[i] = 0;
    fHistTrackEtaPt[i] = 0;
  }

  for (Int_t i = 0; i < 4; i++) {
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtNonBias[i] = 0;
    fHistJetsPtTrack[i] = 0;
    fHistJetsPtClus[i] = 0;
    fHistJetsPt[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtAreaNonBias[i] = 0;
  }
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

  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();

  fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", 100, 0, 100);
  fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentrality->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCentrality);

  fHistZVertex = new TH1F("fHistZVertex","Z vertex position", 60, -30, 30);
  fHistZVertex->GetXaxis()->SetTitle("z");
  fHistZVertex->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistZVertex);

  fHistTracksCent = new TH2F("fHistTracksCent","Tracks vs. centrality", 100, 0, 100, fNbins, 0, 4000);
  fHistTracksCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistTracksCent->GetYaxis()->SetTitle("No. of tracks");
  fOutput->Add(fHistTracksCent);

  fHistJetsCent = new TH2F("fHistJetsCent","Jets vs. centrality", 100, 0, 100, 60, 0, 60);
  fHistJetsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistJetsCent->GetYaxis()->SetTitle("No. of jets");
  fOutput->Add(fHistJetsCent);

  fHistJetsParts = new TH2F("fHistJetsParts","Jets vs. centrality", fNbins, 0, 6000, 60, 0, 60);
  fHistJetsParts->GetXaxis()->SetTitle("No. of particles");
  fHistJetsParts->GetYaxis()->SetTitle("No. of jets");
  fOutput->Add(fHistJetsParts);

  if (fAnaType == kEMCAL || fAnaType == kEMCALOnly) {
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

    fHistClusTimeEnergy = new TH2F("fHistClusTimeEnergy","Time vs. energy of clusters", fNbins, fMinBinPt, fMaxBinPt, fNbins,  -1e-6, 1e-6);
    fHistClusTimeEnergy->GetXaxis()->SetTitle("Energy (GeV)");
    fHistClusTimeEnergy->GetYaxis()->SetTitle("Time");
    fOutput->Add(fHistClusTimeEnergy);

    if (fDoTrigger) {
      fHistMaxL1FastORCent = new TH2F("fHistMaxL1FastORCent","fHistMaxL1ClusCent", 100, 0, 100, 250, 0, 250);
      fHistMaxL1FastORCent->GetXaxis()->SetTitle("Centrality [%]");
      fHistMaxL1FastORCent->GetYaxis()->SetTitle("Maximum L1 FastOR");
      fOutput->Add(fHistMaxL1FastORCent);
      
      fHistMaxL1ClusCent = new TH2F("fHistMaxL1ClusCent","fHistMaxL1ClusCent", 100, 0, 100, 250, 0, 250);
      fHistMaxL1ClusCent->GetXaxis()->SetTitle("Centrality [%]");
      fHistMaxL1ClusCent->GetYaxis()->SetTitle("Maximum L1 trigger cluster");
      fOutput->Add(fHistMaxL1ClusCent);
      
      fHistMaxL1ThrCent = new TH2F("fHistMaxL1ThrCent","fHistMaxL1ThrCent", 100, 0, 100, 250, 0, 250);
      fHistMaxL1ThrCent->GetXaxis()->SetTitle("Centrality [%]");
      fHistMaxL1ThrCent->GetYaxis()->SetTitle("Maximum L1 threshold");
      fOutput->Add(fHistMaxL1ThrCent);
    }
  }

  fHistTracksPt = new TH1F("fHistTracksPt","p_{T} spectrum of reconstructed tracks", fNbins, fMinBinPt, fMaxBinPt);
  fHistTracksPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistTracksPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistTracksPt);
       
  fHistTrPhiEta = new TH2F("fHistTrPhiEta","Phi-Eta distribution of tracks", 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
  fHistTrPhiEta->GetXaxis()->SetTitle("#eta");
  fHistTrPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistTrPhiEta);

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

  if (fRepropagateTracks) {
    fHistDeltaEtaNewProp = new TH1F("fHistDeltaEtaNewProp","fHistDeltaEtaNewProp", 800, -0.4, 0.4);
    fHistDeltaEtaNewProp->GetXaxis()->SetTitle("#delta#eta");
    fHistDeltaEtaNewProp->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaEtaNewProp);
    
    fHistDeltaPhiNewProp = new TH1F("fHistDeltaPhiNewProp","fHistDeltaPhiNewProp", 2560, -1.6, 3.2);
    fHistDeltaPhiNewProp->GetXaxis()->SetTitle("#delta#phi");
    fHistDeltaPhiNewProp->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPhiNewProp);
  }

  if (fAnaType == kEMCAL || fAnaType == kEMCALOnly) {
    fHistClusPhiEtaEnergy = new TH3F("fHistClusPhiEtaEnergy","Phi-Eta-Energy distribution of clusters", 
				     fNbins, fMinBinPt, fMaxBinPt, 100, -1.2, 1.2, 201, 0, TMath::Pi() * 2.01);
    fHistClusPhiEtaEnergy->GetXaxis()->SetTitle("E [GeV]");
    fHistClusPhiEtaEnergy->GetYaxis()->SetTitle("#eta");
    fHistClusPhiEtaEnergy->GetZaxis()->SetTitle("#phi");
    fOutput->Add(fHistClusPhiEtaEnergy);

    fHistNCellsEnergy = new TH2F("fHistNCellsEnergy","Number of cells vs. energy of clusters", fNbins, fMinBinPt, fMaxBinPt, 30, 0, 30);
    fHistNCellsEnergy->GetXaxis()->SetTitle("E [GeV]");
    fHistNCellsEnergy->GetYaxis()->SetTitle("N_{cells}");
    fOutput->Add(fHistNCellsEnergy);
  }

  if (fAnaType == kEMCAL || fAnaType == kEMCALOnly) {
   
    fHistCellsEnergy = new TH1F("fHistCellsEnergy","Energy spectrum of cells", fNbins, fMinBinPt, fMaxBinPt);
    fHistCellsEnergy->GetXaxis()->SetTitle("E [GeV]");
    fHistCellsEnergy->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCellsEnergy);
    
    fHistChVSneCells = new TH2F("fHistChVSneCells","Charged energy vs. neutral (cells) energy", 
				(Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5);
    fHistChVSneCells->GetXaxis()->SetTitle("E [GeV]");
    fHistChVSneCells->GetYaxis()->SetTitle("P [GeV/c]");
    fOutput->Add(fHistChVSneCells);
    
    fHistChVSneClus = new TH2F("fHistChVSneClus","Charged energy vs. neutral (clusters) energy", 
			       (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5);
    fHistChVSneClus->GetXaxis()->SetTitle("E [GeV]");
    fHistChVSneClus->GetYaxis()->SetTitle("P [GeV/c]");
    fOutput->Add(fHistChVSneClus);
    
    fHistChVSneCorrCells = new TH2F("fHistChVSneCorrCells","Charged energy vs. neutral (corrected cells) energy", 
				    (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, (Int_t)(fNbins * 2.5), fMinBinPt , fMaxBinPt * 2.5);
    fHistChVSneCorrCells->GetXaxis()->SetTitle("E [GeV]");
    fHistChVSneCorrCells->GetYaxis()->SetTitle("P [GeV/c]");
    fOutput->Add(fHistChVSneCorrCells);
  }
 
  for (Int_t i = 0; i < 5; i++) {
    TString histnamephi("fHistTrackPhi_");
    histnamephi += i;
    fHistTrackPhi[i] = new TH1F(histnamephi.Data(),histnamephi.Data(), 201, 0, TMath::Pi() * 2.01);
    fHistTrackPhi[i]->GetXaxis()->SetTitle("Phi");
    fOutput->Add(fHistTrackPhi[i]);

    TString histnameeta("fHistTrackEta_");
    histnameeta += i;
    fHistTrackEta[i] = new TH1F(histnameeta.Data(),histnameeta.Data(), 100, -1, 1);
    fHistTrackEta[i]->GetXaxis()->SetTitle("Eta");
    fOutput->Add(fHistTrackEta[i]);
  }

  for (Int_t i = 0; i < 6; i++) {
    TString histnamephipt("fHistTrackPhiPt_");
    histnamephipt += i;
    fHistTrackPhiPt[i] = new TH1F(histnamephipt.Data(),histnamephipt.Data(), 201, 0, TMath::Pi() * 2.01);
    fHistTrackPhiPt[i]->GetXaxis()->SetTitle("Phi");
    fOutput->Add(fHistTrackPhiPt[i]);

    TString histnameetapt("fHistTrackEtaPt_");
    histnameetapt += i;
    fHistTrackEtaPt[i] = new TH1F(histnameetapt.Data(),histnameetapt.Data(), 100, -1, 1);
    fHistTrackEtaPt[i]->GetXaxis()->SetTitle("Eta");
    fOutput->Add(fHistTrackEtaPt[i]);
  }
  
  fHistTrackPhi[0]->SetLineColor(kRed);
  fHistTrackEta[0]->SetLineColor(kRed);
  fHistTrackPhi[1]->SetLineColor(kBlue);
  fHistTrackEta[1]->SetLineColor(kBlue);
  fHistTrackPhi[2]->SetLineColor(kGreen);
  fHistTrackEta[2]->SetLineColor(kGreen);
  fHistTrackPhi[3]->SetLineColor(kOrange);
  fHistTrackEta[3]->SetLineColor(kOrange);
  fHistTrackPhi[4]->SetLineColor(kBlack);
  fHistTrackEta[4]->SetLineColor(kBlack);

  if (!fJetsName.IsNull()) {

    TString histname;

    for (Int_t i = 0; i < 4; i++) {
      histname = "fHistJetsPhiEta_";
      histname += i;
      fHistJetsPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 100, -1.2, 1.2, 201, 0, TMath::Pi() * 2.01);
      fHistJetsPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistJetsPhiEta[i]->GetYaxis()->SetTitle("#phi");
      fHistJetsPhiEta[i]->GetZaxis()->SetTitle("p_{T} [GeV/c]");
      fOutput->Add(fHistJetsPhiEta[i]);

      histname = "fHistJetsPtNonBias_";
      histname += i;
      fHistJetsPtNonBias[i] = new TH1F(histname.Data(), histname.Data(), (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5);
      fHistJetsPtNonBias[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      fHistJetsPtNonBias[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistJetsPtNonBias[i]);

      histname = "fHistJetsPtTrack_";
      histname += i;
      fHistJetsPtTrack[i] = new TH1F(histname.Data(), histname.Data(), (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5);
      fHistJetsPtTrack[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      fHistJetsPtTrack[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistJetsPtTrack[i]);

      if (fAnaType == kEMCAL || fAnaType == kEMCALOnly) {
	histname = "fHistJetsPtClus_";
	histname += i;
	fHistJetsPtClus[i] = new TH1F(histname.Data(), histname.Data(), (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5);
	fHistJetsPtClus[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	fHistJetsPtClus[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistJetsPtClus[i]);
      }

      histname = "fHistJetsPt_";
      histname += i;
      fHistJetsPt[i] = new TH1F(histname.Data(), histname.Data(), (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5);
      fHistJetsPt[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      fHistJetsPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistJetsPt[i]);

      histname = "fHistJetsPtAreaNonBias_";
      histname += i;
      fHistJetsPtAreaNonBias[i] = new TH2F(histname.Data(), histname.Data(), (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, 20, 0, fJetRadius * fJetRadius * TMath::Pi() * 1.5);
      fHistJetsPtAreaNonBias[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      fHistJetsPtAreaNonBias[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistJetsPtAreaNonBias[i]);

      histname = "fHistJetsPtArea_";
      histname += i;
      fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), (Int_t)(fNbins * 2.5), fMinBinPt, fMaxBinPt * 2.5, 20, 0, fJetRadius * fJetRadius * TMath::Pi() * 1.5);
      fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistJetsPtArea[i]);
    }
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAQA::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  fNclusters = 0;
  fNtracks = 0;
  fNjets = 0;

  if (!fTrgClusName.IsNull() && fDoTrigger && !fTrgClusters) {
    fTrgClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrgClusName));
    if (!fTrgClusters) {
      AliError(Form("%s: Could not retrieve trigger clusters %s!", GetName(), fTrgClusName.Data())); 
      return kFALSE;
    }
    else {
      TClass *cl = fTrgClusters->GetClass();
      if (!cl->GetBaseClass("AliVCluster") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fTrgClusName.Data())); 
	fTrgClusters = 0;
	return kFALSE;
      }
    }
  }

  return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSAQA::FillHistograms()
{
  // Fill histograms.

  fHistCentrality->Fill(fCent);
  fHistZVertex->Fill(fVertex[2]);

  Float_t trackSum = DoTrackLoop();

  DoJetLoop();

  fHistTracksCent->Fill(fCent, fNtracks);

  if (fAnaType == kEMCAL || fAnaType == kEMCALOnly) {

    Float_t clusSum = DoClusterLoop();

    fHistClusCent->Fill(fCent, fNclusters);
    fHistClusTracks->Fill(fNtracks, fNclusters);

    Float_t cellSum = 0, cellCutSum = 0;
    
    Int_t ncells = DoCellLoop(cellSum, cellCutSum);

    if (fTracks)
      fHistCellsTracks->Fill(fTracks->GetEntriesFast(), ncells);

    fHistCellsCent->Fill(fCent, ncells);
    
    fHistChVSneCells->Fill(cellSum, trackSum);
    fHistChVSneClus->Fill(clusSum, trackSum);
    fHistChVSneCorrCells->Fill(cellCutSum, trackSum);

    if (fDoTrigger) {
      Float_t maxTrgClus = DoTriggerClusLoop();
      fHistMaxL1ClusCent->Fill(fCent, maxTrgClus);
    
      Int_t maxL1amp = -1;
      Int_t maxL1thr = -1;
    
      DoTriggerPrimitives(maxL1amp, maxL1thr);
      
      if (maxL1amp > -1) 
	fHistMaxL1FastORCent->Fill(fCent, maxL1amp);
      
      if (maxL1thr > -1) 
	fHistMaxL1ThrCent->Fill(fCent, maxL1thr);
    }
  }

  fHistJetsCent->Fill(fCent, fNjets);
  fHistJetsParts->Fill(fNtracks + fNclusters, fNjets);

  return kTRUE;
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
    Float_t amp = cells->GetAmplitude(pos);
    fHistCellsEnergy->Fill(amp);
    sum += amp;
    if (amp < fCellEnergyCut)
      continue;
    sum_cut += amp;
  } 

  return ncells;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAQA::DoClusterLoop()
{
  // Do cluster loop.

  if (!fCaloClusters)
    return 0;

  Float_t sum = 0;

  // Cluster loop
  Int_t nclusters =  fCaloClusters->GetEntriesFast();

  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = static_cast<AliVCluster*>(fCaloClusters->At(iClusters));
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  

    if (!AcceptCluster(cluster, kTRUE))
      continue;

    sum += cluster->E();

    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);

    fHistClusPhiEtaEnergy->Fill(cluster->E(), nPart.Eta(), nPart.Phi());
    fHistNCellsEnergy->Fill(cluster->E(), cluster->GetNCells());

    fHistClusTimeEnergy->Fill(cluster->E(), cluster->GetTOF());

    fNclusters++;
  }

  return sum;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAQA::DoTrackLoop()
{
  // Do track loop.

  if (!fTracks)
    return 0;

  Float_t sum = 0;

  Int_t ntracks = fTracks->GetEntriesFast();
  Int_t nclusters = 0;
  if (fCaloClusters)
    nclusters = fCaloClusters->GetEntriesFast();

  for (Int_t i = 0; i < ntracks; i++) {

    AliVParticle* track = static_cast<AliVParticle*>(fTracks->At(i)); // pointer to reconstructed to track  

    if (!track) {
      AliError(Form("Could not retrieve track %d",i)); 
      continue; 
    }

    AliVTrack* vtrack = dynamic_cast<AliVTrack*>(track); 
    
    if (vtrack && !AcceptTrack(vtrack, kTRUE)) 
      continue;

    fNtracks++;
    
    fHistTracksPt->Fill(track->Pt());

    sum += track->P();

    Int_t label = track->GetLabel();
      
    fHistTrPhiEta->Fill(track->Eta(), track->Phi());
    
    fHistTrackEta[4]->Fill(track->Eta());
    fHistTrackPhi[4]->Fill(track->Phi());

    if (label >= 0 && label < 4) {
      fHistTrackEta[label]->Fill(track->Eta());
      fHistTrackPhi[label]->Fill(track->Phi());
    }

    if (track->Pt() < 0.5) {
      fHistTrackPhiPt[0]->Fill(track->Phi());
      fHistTrackEtaPt[0]->Fill(track->Eta());
    }
    else if (track->Pt() < 1) {
      fHistTrackPhiPt[1]->Fill(track->Phi());
      fHistTrackEtaPt[1]->Fill(track->Eta());
    }
    else if (track->Pt() < 2) {
      fHistTrackPhiPt[2]->Fill(track->Phi());
      fHistTrackEtaPt[2]->Fill(track->Eta());
    }
    else if (track->Pt() < 3) {
      fHistTrackPhiPt[3]->Fill(track->Phi());
      fHistTrackEtaPt[3]->Fill(track->Eta());
    }
    else if (track->Pt() < 5) {
      fHistTrackPhiPt[4]->Fill(track->Phi());
      fHistTrackEtaPt[4]->Fill(track->Eta());
    }
    else {
      fHistTrackPhiPt[5]->Fill(track->Phi());
      fHistTrackEtaPt[5]->Fill(track->Eta());
    }

    if (!vtrack)
      continue;

    if (vtrack->GetTrackEtaOnEMCal() == -999 || vtrack->GetTrackPhiOnEMCal() == -999)
      fHistTrPhiEtaNonProp->Fill(vtrack->Eta(), vtrack->Phi());

    fHistTrEmcPhiEta->Fill(vtrack->GetTrackEtaOnEMCal(), vtrack->GetTrackPhiOnEMCal());
    fHistDeltaEtaPt->Fill(vtrack->Pt(), vtrack->Eta() - vtrack->GetTrackEtaOnEMCal());
    fHistDeltaPhiPt->Fill(vtrack->Pt(), vtrack->Phi() - vtrack->GetTrackPhiOnEMCal());

    if (fRepropagateTracks && vtrack->GetTrackEtaOnEMCal() > -2) {    
      Float_t propeta = -999, propphi = -999;
      PropagateTrack(vtrack, propeta, propphi);
      fHistDeltaEtaNewProp->Fill(propeta - vtrack->GetTrackEtaOnEMCal());
      fHistDeltaPhiNewProp->Fill(propphi - vtrack->GetTrackPhiOnEMCal());
    }
  }
  
  return sum;
}

//____________________________________________________________________________
void AliAnalysisTaskSAQA::PropagateTrack(AliVTrack *track, Float_t &eta, Float_t &phi)
{
  eta = -999;
  phi = -999;

  if (!track)
    return;

  if (track->Pt() == 0) 
    return;

  // init the magnetic field if not already on
  if(!TGeoGlobalMagField::Instance()->GetField()) {
    AliInfo("Init the magnetic field\n");
    AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(InputEvent());
    if (aodevent) {
      Double_t curSol = 30000*aodevent->GetMagneticField()/5.00668;
      Double_t curDip = 6000 *aodevent->GetMuonMagFieldScale();
      AliMagF *field  = AliMagF::CreateFieldMap(curSol,curDip);
      TGeoGlobalMagField::Instance()->SetField(field);
    }
  }
    
  Double_t cv[21];
  for (Int_t i = 0; i < 21; i++) cv[i] = 0;

  Double_t pos[3], mom[3];
  track->GetXYZ(pos);
  track->GetPxPyPz(mom);
  AliExternalTrackParam *trackParam = new AliExternalTrackParam(pos, mom, cv, track->Charge());
 
  if(!AliTrackerBase::PropagateTrackToBxByBz(trackParam, 430., 0.139, 20, kTRUE, 0.8, -1)) return;
  Double_t trkPos[3] = {0., 0., 0.};
  if(!trackParam->GetXYZ(trkPos)) return;
  TVector3 trkPosVec(trkPos[0], trkPos[1], trkPos[2]);
  eta = trkPosVec.Eta();
  phi = trkPosVec.Phi();
  if(phi < 0)
    phi += 2 * TMath::Pi();
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::DoJetLoop()
{
  // Do jet loop.

  if (!fJets)
    return;

  Int_t njets = fJets->GetEntriesFast();

  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  

    if (!AcceptJet(jet, kFALSE))
      continue;

    fHistJetsPtNonBias[fCentBin]->Fill(jet->Pt());
    fHistJetsPtAreaNonBias[fCentBin]->Fill(jet->Pt(), jet->Area());

    fNjets++;

    if (jet->MaxTrackPt() > fPtBiasJetTrack)
      fHistJetsPtTrack[fCentBin]->Fill(jet->Pt());
    
    if (fAnaType == kEMCAL && jet->MaxClusterPt() > fPtBiasJetClus)
      fHistJetsPtClus[fCentBin]->Fill(jet->Pt());
    
    if (!AcceptBiasJet(jet))
      continue;

    fHistJetsPt[fCentBin]->Fill(jet->Pt());
    fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());

    fHistJetsPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());
  }
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAQA::DoTriggerClusLoop()
{
  // Do trigger cluster loop.

  if (!fTrgClusters)
    return 0;

  Int_t ntrgclusters = fTrgClusters->GetEntriesFast();
  Float_t maxTrgClus = 0;

  for (Int_t iClusters = 0; iClusters < ntrgclusters; iClusters++) {
    AliVCluster* cluster = static_cast<AliVCluster*>(fTrgClusters->At(iClusters));
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) continue;

    if (cluster->E() > maxTrgClus)
      maxTrgClus = cluster->E();

  }
  return maxTrgClus;
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::DoTriggerPrimitives(Int_t &maxL1amp, Int_t &maxL1thr)
{
  // Do trigger primitives loop.

  AliVCaloTrigger *triggers = InputEvent()->GetCaloTrigger("EMCAL");

  if (!triggers || triggers->GetEntries() == 0)
    return;
    
  triggers->Reset();
  Int_t L1amp = 0;
  Int_t L1thr = 0;
  maxL1amp = -1;
  maxL1thr = -1;

  while (triggers->Next()) {
    triggers->GetL1TimeSum(L1amp);
    if (maxL1amp < L1amp) 
      maxL1amp = L1amp;

    triggers->GetL1Threshold(L1thr);
    if (maxL1thr < L1thr) 
      maxL1thr = L1thr;
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
