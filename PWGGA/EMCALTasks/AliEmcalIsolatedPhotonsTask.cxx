// $Id: AliEmcalIsolatedPhotonsTask.cxx  $

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDCaloCluster.h"
#include "AliESDtrack.h"
#include "AliEmcalJet.h"
#include "AliFJWrapper.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliPicoTrack.h"

#include "AliEmcalIsolatedPhotonsTask.h"

ClassImp(AliEmcalIsolatedPhotonsTask)

//________________________________________________________________________
AliEmcalIsolatedPhotonsTask::AliEmcalIsolatedPhotonsTask() : 
  AliAnalysisTaskSE("AliEmcalIsolatedPhotonsTask"),
  fOutput(0),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fESDTrackCuts(0),
  fFilterBit(16),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fHistTracksPt(0),
  fHistClustersEnergy(0),
  fHistEPcorrelation(0),
  fHistJetsEnergy(0),
  fHistJetsNE(0),
  fHistJetsNEF(0),
  fHistJetsZ(0),
  fHistTrPhiEta(0),
  fHistClusPhiEta(0),
  fHistJetPhiEta(0),
  Ptbins(100),
  Ptlow(0),
  Ptup(50),
  Ebins(100),
  Elow(0),
  Eup(50)
{
  // Default constructor.
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
}

//________________________________________________________________________
AliEmcalIsolatedPhotonsTask::AliEmcalIsolatedPhotonsTask(const char *name) : 
  AliAnalysisTaskSE("AliEmcalIsolatedPhotonsTask"),
  fOutput(0),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fESDTrackCuts(0),
  fFilterBit(16),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fHistTracksPt(0),
  fHistClustersEnergy(0),
  fHistEPcorrelation(0),
  fHistJetsEnergy(0),
  fHistJetsNE(0),
  fHistJetsNEF(0),
  fHistJetsZ(0),
  fHistTrPhiEta(0),
  fHistClusPhiEta(0),
  fHistJetPhiEta(0),
  Ptbins(100),
  Ptlow(0),
  Ptup(50),
  Ebins(100),
  Elow(0),
  Eup(50)
{
  // Standard constructor.

  if (!name)
    return;

  SetName(name);
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//________________________________________________________________________
AliEmcalIsolatedPhotonsTask::~AliEmcalIsolatedPhotonsTask()
{
  // Destructor
}

//________________________________________________________________________
void AliEmcalIsolatedPhotonsTask::UserCreateOutputObjects()
{
   // Create histograms
  // Called once (on the worker node)
  
  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
  
  AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  
  if( handler && handler->InheritsFrom("AliESDInputHandler") ) {

    if (!fESDTrackCuts) fESDTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    
    if (!fESDTrackCuts) {
      AliFatal("Invalid pointer to ESDtrackCuts");
    }
    
    if (fESDTrackCuts) {
      fESDTrackCuts->DefineHistograms(kRed);
      fOutput->Add(fESDTrackCuts);
    }
  }

  /*
  fGeom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
  if (!fGeom) {
    AliFatal("Unable to get the EMCAL_COMPLETEV1 geom");
  }
  */
    
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
  
  fHistJetsEnergy = new TH1F("fHistJetsEnergy","Energy spectrum of jets", Ebins, Elow, Eup);
  fHistJetsEnergy->GetXaxis()->SetTitle("E [GeV]");
  fHistJetsEnergy->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsEnergy);
  
  fHistJetsNE = new TH1F("fHistJetsNE","Neutral energy spectrum of jets", Ebins, Elow, Eup);
  fHistJetsNE->GetXaxis()->SetTitle("E [GeV]");
  fHistJetsNE->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsNE);
  
  fHistJetsNEF = new TH1F("fHistJetsNEF","Jets neutral energy fraction", Ebins, 0, 1.2);
  fHistJetsNEF->GetXaxis()->SetTitle("E [GeV]");
  fHistJetsNEF->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsNEF);

  fHistJetsZ = new TH1F("fHistJetsZ","Z of jet constituents", Ebins, 0, 1.2);
  fHistJetsZ->GetXaxis()->SetTitle("Z");
  fHistJetsZ->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsZ);

  fHistTrPhiEta = new TH2F("fHistTrPhiEta","Phi-Eta distribution of tracks", 20, -2, 2, 32, 0, 6.4);
  fHistTrPhiEta->GetXaxis()->SetTitle("Eta");
  fHistTrPhiEta->GetYaxis()->SetTitle("Phi");
  fOutput->Add(fHistTrPhiEta);

  fHistClusPhiEta = new TH2F("fHistClusPhiEta","Phi-Eta distribution of clusters", 20, -2, 2, 32, 0, 6.4);
  fHistClusPhiEta->GetXaxis()->SetTitle("Eta");
  fHistClusPhiEta->GetYaxis()->SetTitle("Phi");
  fOutput->Add(fHistClusPhiEta);

  fHistJetPhiEta = new TH2F("fHistJetPhiEta","Phi-Eta distribution of jets", 20, -2, 2, 32, 0, 6.4);
  fHistJetPhiEta->GetXaxis()->SetTitle("Eta");
  fHistJetPhiEta->GetYaxis()->SetTitle("Phi");
  fOutput->Add(fHistJetPhiEta);
	
  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

void AliEmcalIsolatedPhotonsTask::RetrieveEventObjects()
{
  fCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
  fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
  fJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));

  if (!fTracks) {
    AliError(Form("ERROR: Could not retrieve tracks!")); 
  }
  if (!fCaloClusters) {
    AliError(Form("ERROR: Could not retrieve clusters!")); 
  }
  if (!fJets) {
    AliError(Form("ERROR: Could not retrieve jets!")); 
  }
}

AliVTrack* AliEmcalIsolatedPhotonsTask::GetTrack(const Int_t i) const
{
  return dynamic_cast<AliVTrack*>(fTracks->At(i));
}

Int_t AliEmcalIsolatedPhotonsTask::GetNumberOfTracks() const
{
  return fTracks->GetEntriesFast();
}

AliVCluster* AliEmcalIsolatedPhotonsTask::GetCaloCluster(const Int_t i) const
{
  return dynamic_cast<AliVCluster*>(fCaloClusters->At(i));
}

Int_t AliEmcalIsolatedPhotonsTask::GetNumberOfCaloClusters() const
{
  return fCaloClusters->GetEntriesFast();
}

AliEmcalJet* AliEmcalIsolatedPhotonsTask::GetJet(const Int_t i) const
{
  return dynamic_cast<AliEmcalJet*>(fJets->At(i));
}

Int_t AliEmcalIsolatedPhotonsTask::GetNumberOfJets() const
{
  return fJets->GetEntriesFast();
}

void AliEmcalIsolatedPhotonsTask::FillHistograms()
{
  
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
    
    if(track->InheritsFrom("AliExternalTrackParam")) {
      AliExternalTrackParam *trackparam = dynamic_cast<AliExternalTrackParam*>(track);
      fHistTrPhiEta->Fill(trackparam->Eta(), trackparam->Phi());
    }
    else if (track->InheritsFrom("AliAODTrack")) {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
      fHistTrPhiEta->Fill(aodtrack->Eta(), aodtrack->Phi());
    }
    else if (track->InheritsFrom("AliPicoTrack")) {
      AliPicoTrack *picotrack = dynamic_cast<AliPicoTrack*>(track);
      fHistTrPhiEta->Fill(picotrack->Eta(), picotrack->Phi());
    }

    Int_t clId = track->GetEMCALcluster();
    if (clId > -1) {
      AliVCluster* cluster = GetCaloCluster(clId);
      if (cluster)
	fHistEPcorrelation->Fill(track->P(),cluster->E());
    } 
    
  }

  // Jet loop
  Int_t njets =  GetNumberOfJets();
  //cout << njets << " jets" << endl;
  for (Int_t ij = 0; ij < njets; ij++) {
    AliEmcalJet* jet = GetJet(ij);
    if (!jet) {
      printf("ERROR: Could not receive jet %d\n", ij);
      continue;
    }  
    
    fHistJetPhiEta->Fill(jet->Eta(), jet->Phi());
    fHistJetsEnergy->Fill(jet->E());
    fHistJetsNEF->Fill(jet->NEF());
    fHistJetsNE->Fill(jet->E() * jet->NEF());
    //if (jet->E() <= 0)
    //continue;
    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      Int_t trackid = jet->TrackAt(it);
      AliVTrack *track = GetTrack(trackid);
      if (track)
	fHistJetsZ->Fill(track->Pt() / jet->E());
    }
    for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
      Int_t clusterid = jet->ClusterAt(ic);
      AliVCluster *cluster = GetCaloCluster(clusterid);
      if (cluster)
	fHistJetsZ->Fill(cluster->E() / jet->E());
    }
  } //jet loop 
}

//________________________________________________________________________
Bool_t AliEmcalIsolatedPhotonsTask::AcceptTrack(AliVTrack *track)
{
  if (!strcmp(track->ClassName(), "AliESDTrack") && fESDTrackCuts) {
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack*>(track);
    if(esdtrack) {
      return fESDTrackCuts->AcceptTrack(esdtrack);
    }
    //else {
    //cout << "no esdtrack!" << endl;
    //}
  }
  else if (!strcmp(track->ClassName(), "AliAODTrack")) {
    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(track);
    if (aodtrack) {
      //cout << "filter bit = " << fFilterBit << ", filter map = " << aodtrack->GetFilterMap() << endl;
      return aodtrack->TestFilterBit(fFilterBit);
      
    }
  }
  return 1;
}

//________________________________________________________________________
void AliEmcalIsolatedPhotonsTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  // Add jets to event if not yet there

  RetrieveEventObjects();

  FillHistograms();
    
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliEmcalIsolatedPhotonsTask::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
