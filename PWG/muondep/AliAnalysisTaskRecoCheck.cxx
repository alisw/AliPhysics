
/* $Id$ */

#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TIterator.h"

#include "AliStack.h"

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
#include "AliMUONESDInterface.h"
#include "AliMUONRecoCheck.h"

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

  fArray1Mu->Clear();
  fArray2Mu->Clear();

  TLorentzVector v; 

  AliMUONRecoCheck *rc = new AliMUONRecoCheck(fESDEvent, eventHandler);

  AliMUONVTrackStore* trackRefs = rc->TrackRefs(-1);
  AliStack* pstack = mcEvent->Stack();
  
  // loop over ESD tracks
  Int_t nreftracks = 0;
  for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {
    
    AliESDMuonTrack* esdTrack = fESDEvent->GetMuonTrack(iTrack);
    
    // skip ghosts
    if (!esdTrack->ContainTrackerData()) continue;
    
    // convert ESD track to MUON track (without recomputing track parameters at each clusters)
    AliMUONTrack muonTrack;
    AliMUONESDInterface::ESDToMUON(*esdTrack, muonTrack, kFALSE);
    
    // try to match the reconstructed track with a simulated one
    Int_t nMatchClusters = 0;
    AliMUONTrack* matchedTrackRef = rc->FindCompatibleTrack(muonTrack, *trackRefs, nMatchClusters, kFALSE, 10.);
    
    if (matchedTrackRef) {
      
      //store new referenced track in the muonArray with parameters from the reconstructed tracks
      AliMUONTrackLight* muLight = new ((*fArray1Mu)[nreftracks++]) AliMUONTrackLight(esdTrack);
      
      // store further information related to the simulated track
      muLight->SetTrackPythiaLine(matchedTrackRef->GetUniqueID());
      TParticle *part = pstack->Particle(matchedTrackRef->GetUniqueID());
      muLight->SetTrackPDGCode(part->GetPdgCode());
      v.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
      muLight->SetPGen(v); 
      muLight->FillMuonHistory(pstack, part);
      
    }
    
  } // end esd track
  
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
