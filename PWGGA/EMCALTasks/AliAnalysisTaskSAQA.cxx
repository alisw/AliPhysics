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
  fHistTracksPt(0),
  fHistClustersEnergy(0),
  fHistEPcorrelation(0),
  fHistTrPhiEta(0),
  fHistClusPhiEta(0),
  fHistMaxTrgCluster(0),
  Ptbins(100),
  Ptlow(0),
  Ptup(50),
  Ebins(100),
  Elow(0),
  Eup(50)
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
  fHistTracksPt(0),
  fHistClustersEnergy(0),
  fHistEPcorrelation(0),
  fHistTrPhiEta(0),
  fHistClusPhiEta(0),
  fHistMaxTrgCluster(0),
  Ptbins(100),
  Ptlow(0),
  Ptup(50),
  Ebins(100),
  Elow(0),
  Eup(50)
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

  fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", Ebins, 0, 100);
  fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentrality->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCentrality);

  fHistTracksCent = new TH2F("fHistTracksCent","Tracks vs. centrality", Ebins, 0, 100, Ebins, 0, 4000);
  fHistTracksCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistTracksCent->GetYaxis()->SetTitle("No. of tracks");
  fOutput->Add(fHistTracksCent);

  fHistClusCent = new TH2F("fHistClusCent","Clusters vs. centrality", Ebins, 0, 100, Ebins, 0, 2000);
  fHistClusCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistClusCent->GetYaxis()->SetTitle("No. of clusters");
  fOutput->Add(fHistClusCent);
    
  fHistTracksPt = new TH1F("fHistTracksPt","P_{T} spectrum of reconstructed tracks", Ptbins, Ptlow, Ptup);
  fHistTracksPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  fHistTracksPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistTracksPt);
	
  fHistClustersEnergy = new TH1F("fHistClustersEnergy","Energy spectrum of clusters", Ebins, Elow, Eup);
  fHistClustersEnergy->GetXaxis()->SetTitle("E [GeV]");
  fHistClustersEnergy->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistClustersEnergy);

  fHistEPcorrelation = new TH2F("fHistEPcorrelation","Energy-momentum correlation", Ptbins, Ptlow, Ptup, Ebins, Elow, Eup);
  fHistEPcorrelation->GetXaxis()->SetTitle("P [GeV/c]");
  fHistEPcorrelation->GetYaxis()->SetTitle("E [GeV]");
  fOutput->Add(fHistEPcorrelation);

  fHistTrPhiEta = new TH2F("fHistTrPhiEta","Phi-Eta distribution of tracks", 20, -2, 2, 32, 0, 6.4);
  fHistTrPhiEta->GetXaxis()->SetTitle("Eta");
  fHistTrPhiEta->GetYaxis()->SetTitle("Phi");
  fOutput->Add(fHistTrPhiEta);

  fHistClusPhiEta = new TH2F("fHistClusPhiEta","Phi-Eta distribution of clusters", 20, -2, 2, 32, 0, 6.4);
  fHistClusPhiEta->GetXaxis()->SetTitle("Eta");
  fHistClusPhiEta->GetYaxis()->SetTitle("Phi");
  fOutput->Add(fHistClusPhiEta);

  fHistMaxTrgCluster = new TH1F("fHistMaxTrgCluster","Energy distribution of max trigger clusters", Ebins, Elow, Eup);
  fHistMaxTrgCluster->GetXaxis()->SetTitle("E [GeV]");
  fHistMaxTrgCluster->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistMaxTrgCluster);
 
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

  // Cluster loop
  Int_t nclusters =  GetNumberOfCaloClusters();
  //cout << nclusters << " clusters" << endl;
  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = GetCaloCluster(iClusters);
    if (!cluster) {
      printf("ERROR: Could not receive cluster %d\n", iClusters);
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) continue;

    fHistClustersEnergy->Fill(cluster->E());

    Float_t pos[3];
    cluster->GetPosition(pos);
    TVector3 clusVec(pos);
    fHistClusPhiEta->Fill(clusVec.Eta(), clusVec.Phi());
  } //cluster loop 
  

  // Track loop 
  Int_t ntracks = GetNumberOfTracks();
  //cout << ntracks << " tracks" << endl;
  for(Int_t i = 0; i < ntracks; i++) {
    //cout << "track n. " << i << endl;
    AliVTrack* track = GetTrack(i); // pointer to reconstructed to track          
    if(!track) {
      AliError(Form("ERROR: Could not retrieve esdtrack %d",i)); 
      continue; 
    }
    
    if (!AcceptTrack(track)) continue;

    fHistTracksPt->Fill(track->Pt());

    Int_t clId = track->GetEMCALcluster();
    if (clId > -1) {
      AliVCluster* cluster = GetCaloCluster(clId);
      if (cluster)
	fHistEPcorrelation->Fill(track->P(),cluster->E());
    } 

    Float_t eta,phi;
    Int_t label;

    if(track->InheritsFrom("AliESDtrack")) {
      AliESDtrack *esdtrack = dynamic_cast<AliESDtrack*>(track);
      eta = esdtrack->Eta();
      phi = esdtrack->Phi();
      label = esdtrack->GetLabel();
    }
    else if (track->InheritsFrom("AliAODTrack")) {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
      eta = aodtrack->Eta();
      phi = aodtrack->Phi();
      label = aodtrack->GetLabel();
    }
    else if (track->InheritsFrom("AliPicoTrack")) {
      AliPicoTrack *picotrack = dynamic_cast<AliPicoTrack*>(track);
      eta = picotrack->Eta();
      phi = picotrack->Phi();
      label = picotrack->GetLabel();
    }
    else {
      AliWarning("Track type not recognized! Will not plot phi\eta distributions!");
      continue;
    }

    fHistTrPhiEta->Fill(eta, phi);
    
    fHistTrackEta[4]->Fill(eta);
    fHistTrackPhi[4]->Fill(phi);

    if (label >= 0 && label < 4) {
      fHistTrackEta[label]->Fill(eta);
      fHistTrackPhi[label]->Fill(phi);
    }
    else {
      AliWarning(Form("Track label %d not recognized!",label));
    }
   
  }

  Int_t ntrgclusters =  GetNumberOfTrgClusters();
  Float_t maxe = 0;
  //cout << ntrgclusters << " clusters" << endl;
  for (Int_t iClusters = 0; iClusters < ntrgclusters; iClusters++) {
    AliVCluster* cluster = GetTrgCluster(iClusters);
    if (!cluster) {
      printf("ERROR: Could not receive cluster %d\n", iClusters);
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) continue;

    if (cluster->E() > maxe)
      maxe = cluster->E();

  } //cluster loop 
  fHistMaxTrgCluster->Fill(maxe);
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAQA::AcceptTrack(AliVTrack* /*track*/)
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
