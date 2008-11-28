#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TIterator.h"

#include "AliStack.h"
#include "AliMagFMaps.h"
#include "AliTracker.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDMuonTrack.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliAnalysisTaskRecoCheck.h"

#include "AliMUONTrackLight.h"
#include "AliMUONPairLight.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTrack.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONTrackExtrap.h"

// analysis task for decoding reconstructed tracks and kinematics (AliMUONRecoCheck)
// Authors: Bogdan Vulpescu

ClassImp(AliAnalysisTaskRecoCheck)

//________________________________________________________________________
AliAnalysisTaskRecoCheck::AliAnalysisTaskRecoCheck(const char *name) 
  : AliAnalysisTask(name, ""), fESDEvent(0), fTree(0), fArray1Mu(0), fArray2Mu(0), fL3Current(30000.0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TTree container
  DefineOutput(0, TTree::Class());

}

//________________________________________________________________________
void AliAnalysisTaskRecoCheck::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("MuonTracks.*", kTRUE);
    tree->SetBranchStatus("AliESDHeader.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESDEvent = esdH->GetEvent();
  }
  
  // calculate the filed map in the L3 magnet using the current value
  AliMagFMaps* field = 0x0;
  if (fL3Current == 30000.0) {
    field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  }

  // set the tracker field map
  AliTracker::SetFieldMap(field, kFALSE);
  AliMUONTrackExtrap::SetField(AliTracker::GetFieldMap());

}

//________________________________________________________________________
void AliAnalysisTaskRecoCheck::CreateOutputObjects() 
{
  // Create output tree with branches of arrays
  // Called once

  fArray1Mu = new TClonesArray("AliMUONTrackLight",0); 
  fArray2Mu = new TClonesArray("AliMUONPairLight",0); 
  fTree = new TTree("tree","tree",0); 
  fTree->Branch("muons",  &fArray1Mu); 
  fTree->Branch("dimuons",&fArray2Mu); 
  fTree->SetDirectory(0);
  
}

//________________________________________________________________________
void AliAnalysisTaskRecoCheck::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
  // This handler can return the current MC event
  
  Int_t nTracks = fESDEvent->GetNumberOfMuonTracks();
  if (nTracks == 0) return;
  
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    Printf("ERROR: Could not retrieve MC event handler");
    return;
  }

  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if (!mcEvent) {
     Printf("ERROR: Could not retrieve MC event");
     return;
  }

  Int_t eventNumber = fESDEvent->GetEventNumberInFile();

  fArray1Mu->Clear();
  fArray2Mu->Clear();

  TLorentzVector v; 

  AliMUONRecoCheck *rc = new AliMUONRecoCheck(fESDEvent, eventHandler);

  AliMUONVTrackStore* trackRefs = rc->TrackRefs(eventNumber);
  AliMUONVTrackStore* recoTracks = rc->ReconstructedTracks(eventNumber);
  AliStack* pstack = mcEvent->Stack();
  
  TIter next(recoTracks->CreateIterator());
  AliMUONTrack* trackReco = NULL;
  
  Int_t nTrackReco = recoTracks->GetSize();
  Int_t nTracksESD = fESDEvent->GetNumberOfMuonTracks();
  //printf("Tracks in recoTracks (%d) and in ESD (%d).\n", nTrackReco, nTracksESD); 

  if (nTrackReco != nTracksESD) printf ("Tracks in recoTracks (%d) and in ESD (%d) do not match!\n", nTrackReco, nTracksESD);
  
  Int_t nreftracks = 0;
  Int_t itrRec = 0;
  while ( (trackReco = static_cast<AliMUONTrack*>(next())) != NULL ) {
    // assign parameters concerning the reconstructed tracks
    AliMUONTrackLight muLight;
    
    muLight.FillFromESD(fESDEvent->GetMuonTrack(itrRec));
    
    // find the reference track and store further information	
    TParticle *part = muLight.FindRefTrack(trackReco, trackRefs, pstack);
    if (part) { 
      v.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
      muLight.SetPGen(v); 
      muLight.FillMuonHistory(pstack, part);
      //muLight.PrintInfo("A");
      //store the referenced track in the muonArray:
      TClonesArray &muons = *fArray1Mu;
      new (muons[nreftracks++]) AliMUONTrackLight(muLight);
    } 
    
    itrRec++;
  }  // end reco track

  // now loop over muon pairs to build dimuons
  Int_t nmuons = fArray1Mu->GetEntriesFast(); 
  Int_t ndimuons = 0; 
  for(Int_t itrRec1 = 0; itrRec1 < nmuons-1; itrRec1++) { 
    AliMUONTrackLight* mu1 = (AliMUONTrackLight*) fArray1Mu->At(itrRec1); 
    for(Int_t itrRec2 = itrRec1+1; itrRec2 < nmuons; itrRec2++){
      AliMUONTrackLight* mu2 = (AliMUONTrackLight*) fArray1Mu->At(itrRec2); 
      AliMUONPairLight dimuLight;
      dimuLight.SetMuons(*mu1, *mu2);
      //dimuLight.PrintInfo("A");
      TClonesArray &dimuons = *fArray2Mu;
      new (dimuons[ndimuons++]) AliMUONPairLight(dimuLight);
    }
  }
    
  delete rc;

  // Post output data.
  if (nreftracks != 0) {
    fTree->Fill();
  }
  
  PostData(0, fTree);
  
  
}      

//________________________________________________________________________
void AliAnalysisTaskRecoCheck::Terminate(Option_t *) 
{
  // the end

}
