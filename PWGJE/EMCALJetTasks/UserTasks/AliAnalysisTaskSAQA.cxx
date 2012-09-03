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
  fNclusters(0),
  fNtracks(0),
  fNjets(0),
  fHistTracksCent(0),
  fHistClusCent(0),
  fHistJetsCent(0),
  fHistClusTracks(0),
  fHistJetsParts(0),
  fHistCellsCent(0),
  fHistCellsTracks(0),
  fHistTrEmcPhiEta(0),
  fHistTrPhiEtaNonProp(0),
  fHistDeltaEtaPt(0),
  fHistDeltaPhiPt(0),
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

  for (Int_t i = 0; i < 4; i++) {
    fHistTrPhiEtaPt[i] = 0;
    fHistClusPhiEtaEnergy[i] = 0;
    fHistJetsPhiEtaPt[i] = 0;
    fHistJetsPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fCellEnergyCut(0.1),
  fNclusters(0),
  fNtracks(0),
  fNjets(0),
  fHistTracksCent(0),
  fHistClusCent(0),
  fHistJetsCent(0),
  fHistClusTracks(0),
  fHistJetsParts(0),
  fHistCellsCent(0),
  fHistCellsTracks(0),
  fHistTrEmcPhiEta(0),
  fHistTrPhiEtaNonProp(0),
  fHistDeltaEtaPt(0),
  fHistDeltaPhiPt(0),
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

  for (Int_t i = 0; i < 4; i++) {
    fHistTrPhiEtaPt[i] = 0;
    fHistClusPhiEtaEnergy[i] = 0;
    fHistJetsPhiEtaPt[i] = 0;
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

  if (!fTracksName.IsNull()) {
    fHistTracksCent = new TH2F("fHistTracksCent","Tracks vs. centrality", 100, 0, 100, fNbins, 0, 4000);
    fHistTracksCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistTracksCent->GetYaxis()->SetTitle("No. of tracks");
    fOutput->Add(fHistTracksCent);

    TString histname;

    for (Int_t i = 0; i < 4; i++) {
      histname = "fHistTrPhiEtaPt_";
      histname += i;
      fHistTrPhiEtaPt[i] = new TH3F(histname,histname, 100, -1, 1, 201, 0, TMath::Pi() * 2.01, fNbins, fMinBinPt, fMaxBinPt);
      fHistTrPhiEtaPt[i] ->GetXaxis()->SetTitle("#eta");
      fHistTrPhiEtaPt[i] ->GetYaxis()->SetTitle("#phi");
      fHistTrPhiEtaPt[i] ->GetZaxis()->SetTitle("p_{T} (GeV/c)");
      fOutput->Add(fHistTrPhiEtaPt[i]);
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

    for (Int_t i = 0; i < 4; i++) {
      histname = "fHistClusPhiEtaEnergy_";
      histname += i;
      fHistClusPhiEtaEnergy[i] = new TH3F(histname, histname, 100, -1.2, 1.2, 201, 0, TMath::Pi() * 2.01, fNbins, fMinBinPt, fMaxBinPt);
      fHistClusPhiEtaEnergy[i]->GetXaxis()->SetTitle("#eta");
      fHistClusPhiEtaEnergy[i]->GetYaxis()->SetTitle("#phi");
      fHistClusPhiEtaEnergy[i]->GetZaxis()->SetTitle("Energy (GeV)");
      fOutput->Add(fHistClusPhiEtaEnergy[i]);
    }

    fHistClusTimeEnergy = new TH2F("fHistClusTimeEnergy","Time vs. energy of clusters", fNbins, fMinBinPt, fMaxBinPt, fNbins,  -1e-6, 1e-6);
    fHistClusTimeEnergy->GetXaxis()->SetTitle("Energy (GeV)");
    fHistClusTimeEnergy->GetYaxis()->SetTitle("Time");
    fOutput->Add(fHistClusTimeEnergy);

    fHistNCellsEnergy = new TH2F("fHistNCellsEnergy","Number of cells vs. energy of clusters", fNbins, fMinBinPt, fMaxBinPt, 30, 0, 30);
    fHistNCellsEnergy->GetXaxis()->SetTitle("Energy (GeV)");
    fHistNCellsEnergy->GetYaxis()->SetTitle("N_{cells}");
    fOutput->Add(fHistNCellsEnergy); 
     
    fHistCellsEnergy = new TH1F("fHistCellsEnergy","Energy spectrum of cells", fNbins, fMinBinPt, fMaxBinPt);
    fHistCellsEnergy->GetXaxis()->SetTitle("Energy (GeV)");
    fHistCellsEnergy->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCellsEnergy);
    
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

    for (Int_t i = 0; i < 4; i++) {
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

  return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSAQA::FillHistograms()
{
  // Fill histograms.

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

    fHistClusPhiEtaEnergy[fCentBin]->Fill(nPart.Eta(), nPart.Phi(), cluster->E());
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

    sum += track->P();
      
    fHistTrPhiEtaPt[fCentBin]->Fill(track->Eta(), track->Phi(), track->Pt());

    Int_t label = track->GetLabel();

    if (label >= 0 && label < 4) {
      fHistTrackEta[label]->Fill(track->Eta());
      fHistTrackPhi[label]->Fill(track->Phi());
    }

    fHistTrackEta[4]->Fill(track->Eta());
    fHistTrackPhi[4]->Fill(track->Phi());

    if (!vtrack)
      continue;

    if (vtrack->GetTrackEtaOnEMCal() == -999 || vtrack->GetTrackPhiOnEMCal() == -999)
      fHistTrPhiEtaNonProp->Fill(vtrack->Eta(), vtrack->Phi());

    fHistTrEmcPhiEta->Fill(vtrack->GetTrackEtaOnEMCal(), vtrack->GetTrackPhiOnEMCal());
    fHistDeltaEtaPt->Fill(vtrack->Pt(), vtrack->Eta() - vtrack->GetTrackEtaOnEMCal());
    fHistDeltaPhiPt->Fill(vtrack->Pt(), vtrack->Phi() - vtrack->GetTrackPhiOnEMCal());
  }
  
  return sum;
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

    if (!AcceptJet(jet))
      continue;

    fNjets++;

    fHistJetsPhiEtaPt[fCentBin]->Fill(jet->Eta(), jet->Phi(), jet->Pt());
    fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
  }
}


//________________________________________________________________________
void AliAnalysisTaskSAQA::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
