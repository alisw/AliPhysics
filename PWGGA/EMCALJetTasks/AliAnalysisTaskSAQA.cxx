// $Id$
//
// General QA task (S.Aiola).
//
//

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliLog.h"

#include "AliAnalysisTaskSAQA.h"

ClassImp(AliAnalysisTaskSAQA)

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskSAQA"),
  fCellEnergyCut(0.1),
  fDoEoverP(kFALSE),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTrgClusters(0),
  fHistCentrality(0),
  fHistTracksCent(0),
  fHistClusCent(0),
  fHistMaxL1FastORCent(0),
  fHistMaxL1ClusCent(0),
  fHistMaxL1ThrCent(0),
  fHistTracksPt(0),
  fHistTrPhiEta(0),
  fHistClustersEnergy(0),
  fHistClusPhiEta(0),
  fHistJetsPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistEoverP(0),
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
}

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA(const char *name) : 
  AliAnalysisTaskEmcal(name),
  fCellEnergyCut(0.1),
  fDoEoverP(kFALSE),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTrgClusters(0),
  fHistCentrality(0),
  fHistTracksCent(0),
  fHistClusCent(0),
  fHistMaxL1FastORCent(0),
  fHistMaxL1ClusCent(0),
  fHistMaxL1ThrCent(0),
  fHistTracksPt(0),
  fHistTrPhiEta(0),
  fHistClustersEnergy(0),
  fHistClusPhiEta(0),
  fHistJetsPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistEoverP(0),
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
  fOutput->SetOwner();  // IMPORTANT!

  fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", fNbins, 0, 100);
  fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentrality->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCentrality);

  fHistTracksCent = new TH2F("fHistTracksCent","Tracks vs. centrality", fNbins, 0, 100, fNbins, 0, 4000);
  fHistTracksCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistTracksCent->GetYaxis()->SetTitle("No. of tracks");
  fOutput->Add(fHistTracksCent);

  fHistClusCent = new TH2F("fHistClusCent","Clusters vs. centrality", fNbins, 0, 100, fNbins, 0, 2000);
  fHistClusCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistClusCent->GetYaxis()->SetTitle("No. of clusters");
  fOutput->Add(fHistClusCent);

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
    
  fHistTracksPt = new TH1F("fHistTracksPt","P_{T} spectrum of reconstructed tracks", fNbins, fMinPt, fMaxPt);
  fHistTracksPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  fHistTracksPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistTracksPt);
       
  fHistTrPhiEta = new TH2F("fHistTrPhiEta","Phi-Eta distribution of tracks", 20, -2, 2, 32, 0, 6.4);
  fHistTrPhiEta->GetXaxis()->SetTitle("#eta");
  fHistTrPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistTrPhiEta);

  fHistClustersEnergy = new TH1F("fHistClustersEnergy","Energy spectrum of clusters", fNbins, fMinPt, fMaxPt);
  fHistClustersEnergy->GetXaxis()->SetTitle("E [GeV]");
  fHistClustersEnergy->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistClustersEnergy);

  fHistClusPhiEta = new TH2F("fHistClusPhiEta","Phi-Eta distribution of clusters", 20, -2, 2, 32, 0, 6.4);
  fHistClusPhiEta->GetXaxis()->SetTitle("#eta");
  fHistClusPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistClusPhiEta);

  fHistJetsPt = new TH1F("fHistJetsPt","P_{T} spectrum of reconstructed jets", fNbins, fMinPt, fMaxPt);
  fHistJetsPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  fHistJetsPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsPt);
       
  fHistJetsPhiEta = new TH2F("fHistJetsPhiEta","Phi-Eta distribution of jets", 20, -2, 2, 32, 0, 6.4);
  fHistJetsPhiEta->GetXaxis()->SetTitle("#eta");
  fHistJetsPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJetsPhiEta);

  fHistJetsPtArea = new TH2F("fHistJetsPtArea","P_{T} vs. area of reconstructed jets", fNbins, fMinPt, fMaxPt, 20, 0, fJetRadius * fJetRadius * TMath::Pi() * 1.5);
  fHistJetsPtArea->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  fHistJetsPtArea->GetYaxis()->SetTitle("area");
  fOutput->Add(fHistJetsPtArea);

  if (fDoEoverP) {
    fHistEoverP = new TH2F("fHistEoverP","E/P vs. E", fNbins, fMinPt, fMaxPt, fNbins, 0, 10);
    fHistEoverP->GetXaxis()->SetTitle("E [GeV]");
    fHistEoverP->GetYaxis()->SetTitle("E/P [c]");
    fOutput->Add(fHistEoverP);
  }

  fHistCellsEnergy = new TH1F("fHistCellsEnergy","Energy spectrum of cells", fNbins, fMinPt, fMaxPt);
  fHistCellsEnergy->GetXaxis()->SetTitle("E [GeV]");
  fHistCellsEnergy->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCellsEnergy);

  fHistChVSneCells = new TH2F("fHistChVSneCells","Charged energy vs. neutral (cells) energy", fNbins, fMinPt * 4, fMaxPt * 4, fNbins, fMinPt * 4, fMaxPt * 4);
  fHistChVSneCells->GetXaxis()->SetTitle("E [GeV]");
  fHistChVSneCells->GetYaxis()->SetTitle("P [GeV/c]");
  fOutput->Add(fHistChVSneCells);

  fHistChVSneClus = new TH2F("fHistChVSneClus","Charged energy vs. neutral (clusters) energy", fNbins, fMinPt * 4, fMaxPt * 4, fNbins, fMinPt * 4, fMaxPt * 4);
  fHistChVSneClus->GetXaxis()->SetTitle("E [GeV]");
  fHistChVSneClus->GetYaxis()->SetTitle("P [GeV/c]");
  fOutput->Add(fHistChVSneClus);

  fHistChVSneCorrCells = new TH2F("fHistChVSneCorrCells","Charged energy vs. neutral (corrected cells) energy", fNbins, fMinPt * 4, fMaxPt * 4, fNbins, fMinPt * 4, fMaxPt * 4);
  fHistChVSneCorrCells->GetXaxis()->SetTitle("E [GeV]");
  fHistChVSneCorrCells->GetYaxis()->SetTitle("P [GeV/c]");
  fOutput->Add(fHistChVSneCorrCells);
 
  for (Int_t i = 0; i < 5; i++) {
    TString histnamephi("fHistTrackPhi_");
    histnamephi += i;
    fHistTrackPhi[i] = new TH1F(histnamephi.Data(),histnamephi.Data(), 128, 0, 6.4);
    fHistTrackPhi[i]->GetXaxis()->SetTitle("Phi");
    fOutput->Add(fHistTrackPhi[i]);

    TString histnameeta("fHistTrackEta_");
    histnameeta += i;
    fHistTrackEta[i] = new TH1F(histnameeta.Data(),histnameeta.Data(), 100, -2, 2);
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

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::RetrieveEventObjects()
{
  AliAnalysisTaskEmcal::RetrieveEventObjects();

  if (strcmp(fTrgClusName,"")) {
    fTrgClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrgClusName));
    if (!fTrgClusters) {
      AliWarning(Form("Could not retrieve trigger clusters!")); 
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::FillHistograms()
{
  fHistCentrality->Fill(fCent);
  if (fTracks)
    fHistTracksCent->Fill(fCent, fTracks->GetEntriesFast());
  if (fCaloClusters)
    fHistClusCent->Fill(fCent, fCaloClusters->GetEntriesFast());

  Float_t clusSum = DoClusterLoop();
  
  Float_t trackSum = DoTrackLoop();

  Float_t cellSum = 0, cellCutSum = 0;
  DoCellLoop(cellSum, cellCutSum);

  fHistChVSneCells->Fill(cellSum, trackSum);
  fHistChVSneClus->Fill(clusSum, trackSum);
  fHistChVSneCorrCells->Fill(cellCutSum, trackSum);

  DoJetLoop();

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

//________________________________________________________________________
void AliAnalysisTaskSAQA::DoCellLoop(Float_t &sum, Float_t &sum_cut)
{
  AliVCaloCells *cells = InputEvent()->GetEMCALCells();

  if (!cells)
    return;

  Int_t ncells = cells->GetNumberOfCells();

  for (Int_t pos = 0; pos < ncells; pos++) {

    Float_t amp = cells->GetAmplitude(pos);

    fHistCellsEnergy->Fill(amp);

    sum += amp;

    if (amp < fCellEnergyCut)
      continue;

    sum_cut += amp;

  } 
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAQA::DoClusterLoop()
{
  if (!fCaloClusters)
    return 0;

  Float_t sum = 0;

  // Cluster loop
  Int_t nclusters =  fCaloClusters->GetEntriesFast();

  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = dynamic_cast<AliVCluster*>(fCaloClusters->At(iClusters));
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) continue;

    if (!AcceptCluster(cluster, kTRUE)) continue;

    fHistClustersEnergy->Fill(cluster->E());

    sum += cluster->E();

    Float_t pos[3];
    cluster->GetPosition(pos);
    TVector3 clusVec(pos);
    fHistClusPhiEta->Fill(clusVec.Eta(), clusVec.Phi());

  } //cluster loop 

  return sum;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAQA::DoTrackLoop()
{
  if (!fTracks)
    return 0;

  Float_t sum = 0;

  // Track loop 
  Int_t ntracks = fTracks->GetEntriesFast();
  Int_t nclusters = 0;
  if (fCaloClusters)
    nclusters = fCaloClusters->GetEntriesFast();

  for(Int_t i = 0; i < ntracks; i++) {

    AliVTrack* track = dynamic_cast<AliVTrack*>(fTracks->At(i)); // pointer to reconstructed to track          
    if(!track) {
      AliError(Form("Could not retrieve track %d",i)); 
      continue; 
    }
    
    if (!AcceptTrack(track, kTRUE)) continue;

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

    if (!fDoEoverP)
      continue;

    Int_t clId = track->GetEMCALcluster();
    if (clId > -1 && clId < nclusters && fCaloClusters) {
      AliVCluster* cluster = dynamic_cast<AliVCluster*>(fCaloClusters->At(i));
      if (cluster) {
	fHistEoverP->Fill(cluster->E(), cluster->E() / track->P());
      }
    } 
  }
  
  return sum;
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::DoJetLoop()
{
  if (!fJets)
    return;

  Int_t njets = fJets->GetEntriesFast();

  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fJets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;
    
    fHistJetsPt->Fill(jet->Pt());
    fHistJetsPtArea->Fill(jet->Pt(), jet->Area());
    fHistJetsPhiEta->Fill(jet->Eta(), jet->Phi());
  }
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAQA::DoTriggerClusLoop()
{
  if (!fTrgClusters)
    return 0;

  Int_t ntrgclusters = fTrgClusters->GetEntriesFast();
  Float_t maxTrgClus = 0;

  for (Int_t iClusters = 0; iClusters < ntrgclusters; iClusters++) {
    AliVCluster* cluster = dynamic_cast<AliVCluster*>(fTrgClusters->At(iClusters));
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
