#include "AliLog.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "AliMuonForwardTrack.h"

TFile *inputFileWithoutBranson=0, *inputFileWithBranson=0, *outputFile=0;
TTree *fInputTreeWithoutBranson=0, *fInputTreeWithBranson=0;
TTree *fOutputTreeWithoutBranson=0, *fOutputTreeWithBranson=0;
TClonesArray *fInputMuonForwardTracksWithoutBranson=0, *fInputMuonForwardTracksWithBranson=0;
TClonesArray *fOutputMuonForwardTracksWithoutBranson=0, *fOutputMuonForwardTracksWithBranson=0;

AliMuonForwardTrack* GetTrackWithoutBranson(Int_t event, Int_t trackMCId);

//====================================================================================================================================================

void FilterMuonGlobalTracks() {

  // ---- reading file of tracks reconstructed without Branson correction

  inputFileWithoutBranson = new TFile("MuonGlobalTracks.withoutBransonCorrection.root");
  if (!inputFileWithoutBranson || !inputFileWithoutBranson->IsOpen()) {
    printf("Error opening file MuonGlobalTracks.withoutBransonCorrection.root");
    return;
  }
  fInputTreeWithoutBranson = (TTree*) inputFileWithoutBranson->Get("AliMuonForwardTracks");
  if (!fInputTreeWithoutBranson) {
    printf("Error reading input tree from MuonGlobalTracks.withoutBransonCorrection.root");
    return;
  }
  fInputTreeWithoutBranson->SetName("AliInputMuonForwardTracksWithoutBranson");
  fInputMuonForwardTracksWithoutBranson = new TClonesArray("AliMuonForwardTrack");
  fInputTreeWithoutBranson->SetBranchAddress("tracks", &fInputMuonForwardTracksWithoutBranson);  

  // ---- reading file of tracks reconstructed with Branson correction

  inputFileWithBranson = new TFile("MuonGlobalTracks.withBransonCorrection.root");
  if (!inputFileWithBranson || !inputFileWithBranson->IsOpen()) {
    printf("Error opening file MuonGlobalTracks.withBransonCorrection.root");
    return;
  }
  fInputTreeWithBranson = (TTree*) inputFileWithBranson->Get("AliMuonForwardTracks");
  if (!fInputTreeWithBranson) {
    printf("Error reading input tree from MuonGlobalTracks.withBransonCorrection.root");
    return;
  }
  fInputTreeWithBranson->SetName("AliInputMuonForwardTracksWithBranson");
  fInputMuonForwardTracksWithBranson = new TClonesArray("AliMuonForwardTrack");
  fInputTreeWithBranson->SetBranchAddress("tracks", &fInputMuonForwardTracksWithBranson);  

  // ---- preparing output file and trees
  
  outputFile = new TFile("MuonGlobalTracks.root","recreate");
  fOutputTreeWithoutBranson = new TTree("AliMuonForwardTracksWithoutBranson", "Tree of AliMuonForwardTracks");
  fOutputTreeWithBranson    = new TTree("AliMuonForwardTracksWithBranson",    "Tree of AliMuonForwardTracks");
  fOutputMuonForwardTracksWithoutBranson = new TClonesArray("AliMuonForwardTrack");
  fOutputMuonForwardTracksWithBranson    = new TClonesArray("AliMuonForwardTrack");
  fOutputTreeWithoutBranson -> Branch("tracks", &fOutputMuonForwardTracksWithoutBranson);
  fOutputTreeWithBranson    -> Branch("tracks", &fOutputMuonForwardTracksWithBranson);

  // ---- 

  Int_t nEvents = fInputTreeWithBranson->GetEntries();

  for (Int_t iEv=0; iEv<nEvents; iEv++) {

    fInputTreeWithBranson -> GetEvent(iEv);

    printf("reading event %d\n", iEv);
    
    Int_t nTracksOfEvent = fInputMuonForwardTracksWithBranson->GetEntries();
    for (Int_t iTr=0; iTr<nTracksOfEvent; iTr++) {
      AliMuonForwardTrack *trackWithBranson = (AliMuonForwardTrack*) fInputMuonForwardTracksWithBranson->At(iTr);
      AliMuonForwardTrack *trackWithoutBranson = GetTrackWithoutBranson(iEv, trackWithBranson->GetTrackMCId());
      if (!trackWithoutBranson) {
	printf("Event %d, track %d : track without Branson not available!\n", iEv, iTr);
	continue;
      }
      new ((*fOutputMuonForwardTracksWithBranson)[fOutputMuonForwardTracksWithBranson->GetEntries()]) AliMuonForwardTrack(*trackWithBranson);      
      new ((*fOutputMuonForwardTracksWithoutBranson)[fOutputMuonForwardTracksWithoutBranson->GetEntries()]) AliMuonForwardTrack(*trackWithoutBranson);    
    }

    fOutputTreeWithBranson    -> Fill();
    fOutputTreeWithoutBranson -> Fill();

    fOutputMuonForwardTracksWithBranson    -> Delete();
    fOutputMuonForwardTracksWithoutBranson -> Delete();

    fInputMuonForwardTracksWithBranson    -> Delete();
    fInputMuonForwardTracksWithoutBranson -> Delete();

  }

  outputFile -> cd();
  fOutputTreeWithBranson    -> Write();
  fOutputTreeWithoutBranson -> Write();
  outputFile -> Close();

}

//====================================================================================================================================================

AliMuonForwardTrack* GetTrackWithoutBranson(Int_t event, Int_t trackMCId) {

  fInputTreeWithoutBranson -> GetEvent(event);
  
  Int_t nTracksOfEvent = fInputMuonForwardTracksWithoutBranson->GetEntries();
  for (Int_t iTr=0; iTr<nTracksOfEvent; iTr++) {
    AliMuonForwardTrack *trackWithoutBranson = (AliMuonForwardTrack*) fInputMuonForwardTracksWithoutBranson->At(iTr);
    if (trackWithoutBranson->GetTrackMCId() == trackMCId) return trackWithoutBranson;
  }
  
  return NULL;
  
}

//====================================================================================================================================================
