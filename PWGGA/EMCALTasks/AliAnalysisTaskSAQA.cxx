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
#include <TParticle.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliESDtrack.h"
#include "AliEmcalJet.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVEventHandler.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskSAQA.h"

ClassImp(AliAnalysisTaskSAQA)

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA() : 
  AliAnalysisTaskSE("AliAnalysisTaskSAQA"),
  fOutput(0),
  fCellEnergyCut(0.1),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fTrgClusters(0),
  fCent(0),
  fHistCentrality(0),
  fHistTracksCent(0),
  fHistClusCent(0),
  fHistMaxL1FastORCent(0),
  fHistMaxL1ClusCent(0),
  fHistMaxL1ThrCent(0),
  fHistTracksPt(0),
  fHistCellsEnergy(0),
  fHistClustersEnergy(0),
  fHistEoverP(0),
  fHistTrPhiEta(0),
  fHistClusPhiEta(0),
  fHistChVSneCells(0),
  fHistChVSneClus(0),
  fHistChVSneCorrCells(0),
  fNbins(250),
  fMinPt(0),
  fMaxPt(25)
{
  // Default constructor.

  for (Int_t i = 0; i < 5; i++) {
    fHistTrackPhi[i] = 0;
    fHistTrackEta[i] = 0;
  }

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
}

//________________________________________________________________________
AliAnalysisTaskSAQA::AliAnalysisTaskSAQA(const char *name) : 
  AliAnalysisTaskSE(name),
  fOutput(0),
  fCellEnergyCut(0.1),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fTrgClusters(0),
  fCent(0),
  fHistCentrality(0),
  fHistTracksCent(0),
  fHistClusCent(0),
  fHistMaxL1FastORCent(0),
  fHistMaxL1ClusCent(0),
  fHistMaxL1ThrCent(0),
  fHistTracksPt(0),
  fHistCellsEnergy(0),
  fHistClustersEnergy(0),
  fHistEoverP(0),
  fHistTrPhiEta(0),
  fHistClusPhiEta(0),
  fHistChVSneCells(0),
  fHistChVSneClus(0),
  fHistChVSneCorrCells(0),
  fNbins(250),
  fMinPt(0),
  fMaxPt(25)
{
  // Standard constructor.

  for (Int_t i = 0; i < 5; i++) {
    fHistTrackPhi[i] = 0;
    fHistTrackEta[i] = 0;
  }

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
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

  fHistCellsEnergy = new TH1F("fHistCellsEnergy","Energy spectrum of cells", fNbins, fMinPt, fMaxPt);
  fHistCellsEnergy->GetXaxis()->SetTitle("E [GeV]");
  fHistCellsEnergy->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCellsEnergy);
	
  fHistClustersEnergy = new TH1F("fHistClustersEnergy","Energy spectrum of clusters", fNbins, fMinPt, fMaxPt);
  fHistClustersEnergy->GetXaxis()->SetTitle("E [GeV]");
  fHistClustersEnergy->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistClustersEnergy);

  fHistEoverP = new TH2F("fHistEoverP","E/P vs. E", fNbins, fMinPt, fMaxPt, fNbins, 0, 10);
  fHistEoverP->GetXaxis()->SetTitle("E [GeV]");
  fHistEoverP->GetYaxis()->SetTitle("E/P [c]");
  fOutput->Add(fHistEoverP);

  fHistTrPhiEta = new TH2F("fHistTrPhiEta","Phi-Eta distribution of tracks", 20, -2, 2, 32, 0, 6.4);
  fHistTrPhiEta->GetXaxis()->SetTitle("Eta");
  fHistTrPhiEta->GetYaxis()->SetTitle("Phi");
  fOutput->Add(fHistTrPhiEta);

  fHistClusPhiEta = new TH2F("fHistClusPhiEta","Phi-Eta distribution of clusters", 20, -2, 2, 32, 0, 6.4);
  fHistClusPhiEta->GetXaxis()->SetTitle("Eta");
  fHistClusPhiEta->GetYaxis()->SetTitle("Phi");
  fOutput->Add(fHistClusPhiEta);

  fHistChVSneCells = new TH2F("fHistChVSneCells","Charged energy vs. neutral (cells) energy", fNbins, fMinPt * 10, fMaxPt * 10, fNbins, fMinPt * 10, fMaxPt * 10);
  fHistChVSneCells->GetXaxis()->SetTitle("E [GeV]");
  fHistChVSneCells->GetYaxis()->SetTitle("P [GeV/c]");
  fOutput->Add(fHistChVSneCells);

  fHistChVSneClus = new TH2F("fHistChVSneClus","Charged energy vs. neutral (clusters) energy", fNbins, fMinPt * 10, fMaxPt * 10, fNbins, fMinPt * 10, fMaxPt * 10);
  fHistChVSneClus->GetXaxis()->SetTitle("E [GeV]");
  fHistChVSneClus->GetYaxis()->SetTitle("P [GeV/c]");
  fOutput->Add(fHistChVSneClus);

  fHistChVSneCorrCells = new TH2F("fHistChVSneCorrCells","Charged energy vs. neutral (corrected cells) energy", fNbins, fMinPt * 10, fMaxPt * 10, fNbins, fMinPt * 10, fMaxPt * 10);
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

void AliAnalysisTaskSAQA::RetrieveEventObjects()
{
  fCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
  if (!fCaloClusters) {
    AliWarning(Form("Could not retrieve clusters!")); 
  }

  fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
  if (!fTracks) {
    AliWarning(Form("Could not retrieve tracks!")); 
  }

  if (strcmp(fTrgClusName,"")) {
    fTrgClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrgClusName));
    if (!fTrgClusters) {
      AliWarning(Form("Could not retrieve trigger clusters!")); 
    }
  }

  fCent = InputEvent()->GetCentrality();
}

AliVTrack* AliAnalysisTaskSAQA::GetTrack(const Int_t i) const
{
  if (fTracks)
    return dynamic_cast<AliVTrack*>(fTracks->At(i));
  else
    return 0;
}

Int_t AliAnalysisTaskSAQA::GetNumberOfTracks() const
{
  if (fTracks)
    return fTracks->GetEntriesFast();
  else
    return 0;
}

AliVCluster* AliAnalysisTaskSAQA::GetCaloCluster(const Int_t i) const
{ 
  if (fCaloClusters)
    return dynamic_cast<AliVCluster*>(fCaloClusters->At(i));
  else
    return 0;
}

Int_t AliAnalysisTaskSAQA::GetNumberOfCaloClusters() const
{ 
  if (fCaloClusters)
    return fCaloClusters->GetEntriesFast();
  else
    return 0;
}

AliVCluster* AliAnalysisTaskSAQA::GetTrgCluster(const Int_t i) const
{
  if (fTrgClusters)
    return dynamic_cast<AliVCluster*>(fTrgClusters->At(i));
  else
    return 0;
}

Int_t AliAnalysisTaskSAQA::GetNumberOfTrgClusters() const
{
  if (fTrgClusters)
    return fTrgClusters->GetEntriesFast();
  else
    return 0;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAQA::GetAcceptanceNormFactor() const
{
  const Float_t EmcalEtaMin = -0.7;
  const Float_t EmcalEtaMax = 0.7;
  const Float_t EmcalPhiMin = 80 * TMath::DegToRad();
  const Float_t EmcalPhiMax = 180 * TMath::DegToRad();

  const Float_t TpcEtaMin = -0.9;
  const Float_t TpcEtaMax = 0.9;
  const Float_t TpcPhiMin = 0;
  const Float_t TpcPhiMax = 2 * TMath::Pi();

  Float_t emcalArea = (EmcalEtaMax - EmcalEtaMin) * (EmcalPhiMax - EmcalPhiMin);
  Float_t tpcArea   = (TpcEtaMax - TpcEtaMin) * (TpcPhiMax - TpcPhiMin);

  return emcalArea / tpcArea;
}

void AliAnalysisTaskSAQA::FillHistograms()
{
  Float_t cent = 100;
  
  if (fCent)
    cent = fCent->GetCentralityPercentile("V0M");
  else
    AliWarning("Centrality not available!");

  fHistCentrality->Fill(cent);
  fHistTracksCent->Fill(cent, GetNumberOfTracks());
  fHistClusCent->Fill(cent, GetNumberOfCaloClusters());

  Float_t clusSum = DoClusterLoop();
  
  Float_t trackSum = DoTrackLoop();

  //Normalization to EMCal acceptance
  trackSum *= GetAcceptanceNormFactor();

  Float_t cellSum = 0, cellCutSum = 0;
  DoCellLoop(cellSum, cellCutSum);

  fHistChVSneCells->Fill(cellSum, trackSum);
  fHistChVSneClus->Fill(clusSum, trackSum);
  fHistChVSneCorrCells->Fill(cellCutSum, trackSum);

  Float_t maxTrgClus = DoTriggerClusLoop();
  fHistMaxL1ClusCent->Fill(cent, maxTrgClus);
  
  Int_t maxL1amp = -1;
  Int_t maxL1thr = -1;

  DoTriggerPrimitives(maxL1amp, maxL1thr);

  if (maxL1amp > -1) 
    fHistMaxL1FastORCent->Fill(cent, maxL1amp);

  if (maxL1thr > -1) 
    fHistMaxL1ThrCent->Fill(cent, maxL1thr);
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
  Float_t sum = 0;

  // get primary vertex
  //Double_t vertex[3] = {0, 0, 0};
  //InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Cluster loop
  Int_t nclusters =  GetNumberOfCaloClusters();

  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = GetCaloCluster(iClusters);
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) continue;

    fHistClustersEnergy->Fill(cluster->E());

    //TLorentzVector nPart;
    //cluster->GetMomentum(nPart, vertex);
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
  Float_t sum = 0;

  // Track loop 
  Int_t nclusters =  GetNumberOfCaloClusters();
  Int_t ntracks = GetNumberOfTracks();

  for(Int_t i = 0; i < ntracks; i++) {

    AliVTrack* track = GetTrack(i); // pointer to reconstructed to track          
    if(!track) {
      AliError(Form("Could not retrieve esdtrack %d",i)); 
      continue; 
    }
    
    if (!AcceptTrack(track)) continue;

    fHistTracksPt->Fill(track->Pt());

    sum += track->P();

    Int_t clId = track->GetEMCALcluster();
    if (clId > -1 && clId < nclusters) {
      AliVCluster* cluster = GetCaloCluster(clId);
      if (cluster) {
	fHistEoverP->Fill(cluster->E(), cluster->E() / track->P());
      }
    } 

    Int_t label = track->GetLabel();

    fHistTrPhiEta->Fill(track->Eta(), track->Phi());
    
    fHistTrackEta[4]->Fill(track->Eta());
    fHistTrackPhi[4]->Fill(track->Phi());

    if (label >= 0 && label < 4) {
      fHistTrackEta[label]->Fill(track->Eta());
      fHistTrackPhi[label]->Fill(track->Phi());
    }
    else {
      AliWarning(Form("Track label %d not recognized!",label));
    }
  }
  
  return sum;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAQA::DoTriggerClusLoop()
{
  Int_t ntrgclusters =  GetNumberOfTrgClusters();
  Float_t maxTrgClus = 0;

  for (Int_t iClusters = 0; iClusters < ntrgclusters; iClusters++) {
    AliVCluster* cluster = GetTrgCluster(iClusters);
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
  //Float_t L0FastORamp = 0;
  Int_t L1amp = 0;
  Int_t L1thr = 0;
  maxL1amp = -1;
  maxL1thr = -1;

  while (triggers->Next()) {
    /*
    triggers->GetAmplitude(L0FastORamp);
    
    if (L0FastORamp < 0)
      continue;
    
    Int_t ntimes = 0;
      
    triggers->GetNL0Times(ntimes);
    
    if (ntimes > 0) {
      Int_t trgtimes[25];
      triggers->GetL0Times(trgtimes);
      
      Int_t mintime = trgtimes[0];
      Int_t maxtime = trgtimes[0];
        
      for (Int_t i = 0; i < ntimes; ++i) {
	if (trgtimes[i] < mintime)
	  mintime = trgtimes[i];
	if (maxtime < trgtimes[i])
	  maxtime = trgtimes[i];
      }
    }
    
    Int_t gCol = 0, gRow = 0;
    triggers->GetPosition(gCol, gRow);
    
    Int_t find = -1;
    fGeom->GetAbsFastORIndexFromPositionInEMCAL(gCol, gRow, find);
      
    if (find < 0)
      continue;
    
    Int_t cidx[4] = {-1};
    Bool_t ret = fGeom->GetCellIndexFromFastORIndex(find, cidx);
    
    if (!ret)
      continue;
    */

    triggers->GetL1TimeSum(L1amp);
    if (maxL1amp < L1amp) 
      maxL1amp = L1amp;

    triggers->GetL1Threshold(L1thr);
    if (maxL1thr < L1thr) 
      maxL1thr = L1thr;
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAQA::AcceptTrack(AliVTrack* /*track*/) const
{
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  // Add jets to event if not yet there

  RetrieveEventObjects();

  FillHistograms();
    
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskSAQA::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
