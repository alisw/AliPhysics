#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

//#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliVEventHandler.h"
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
#include "AliVfriendEvent.h"
#include "AliVfriendTrack.h"
#include "AliESDInputHandler.h"

#include "AliAnalysisTaskPt.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

ClassImp(AliAnalysisTaskPt)

//________________________________________________________________________
AliAnalysisTaskPt::AliAnalysisTaskPt(const char *name) 
: AliAnalysisTask(name, ""), fESD(0), fESDfriend(0), fHistPt(0), fCuts(0), fEv(0), fHistQ(0), fListOut(0), fUseFriends(kFALSE), fHistNTPCCl(0), fHistNESDtracks(0),   fHistNESDfriendtracks(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskPt::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  printf("----> AliAnalysisTaskPt::ConnectInputData\n");
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
    /*
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("fTracks.*", kTRUE);
    */

    AliVEventHandler *esdH = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    TString classInputHandler = esdH->ClassName();

    Printf("----> AliAnalysisTaskPt: ClassName of handler = %s", classInputHandler.Data());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } 
    else {
      Printf("----> AliAnalysisTaskPt::ConnectInputData Getting the event from handler %p", esdH);
      fESD = esdH->GetEvent();
      if (fUseFriends){
	Printf("...We have to use the friends...");
	if (classInputHandler.Contains("HLT")) { // we are running in HLT
	  fESDfriend = esdH->GetVfriendEvent();
	}
	else { /// we are running offline
	  if (esdH && esdH->GetTree()) {
	    Printf("...We got the tree...");
	    if (esdH->GetTree()->GetBranch("ESDfriend.")){
	      Printf("Yu-huuuu!!! friend branch found");
	      fESDfriend = ((AliESDInputHandler*)esdH)->GetESDfriend();
	    }
	    else {
	      Printf("No friend branch found");
	    }
	  }
	}	
	Printf("and the result is: fESDfriend = %p", fESDfriend);
      }
      else {
	Printf("The friends are not requested");
      }
    }
    if (!fESD) {
      Printf("ERROR, no ESD event");
    }
    if (fUseFriends && !fESDfriend) {
      Printf("ERROR, no ESD friend");
    }
  }

  Printf("fESD = %p, fESDfriend = %p", fESD, fESDfriend);
  printf("<---- AliAnalysisTaskPt::ConnectInputData\n");
}

//________________________________________________________________________
void AliAnalysisTaskPt::CreateOutputObjects()
{
  // Create histograms
  // Called once

  fListOut = new TList();
  fListOut->SetOwner();
  fListOut->SetName("listHistos");

  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 3.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);

  fHistQ = new TH1F("fHistQ", "TPC clusters: Q distribution", 1000, 0, 10000);
  fHistQ->GetXaxis()->SetTitle("Q");
  fHistQ->GetYaxis()->SetTitle("dN/dQ");
  fHistQ->SetMarkerStyle(kFullCircle);

  fHistNTPCCl = new TH1F("fHistNTPCCl", "Number of TPC clusters", 160, -0.5, 159.5);
  fHistNTPCCl->GetXaxis()->SetTitle("n. TPC Cl.");
  fHistNTPCCl->GetYaxis()->SetTitle("dN/d(n. TPC Cl)");
  fHistNTPCCl->SetMarkerStyle(kFullCircle);

  fHistNESDtracks = new TH1F("fHistNESDtracks", "Number of ESD tracks", 1000, -0.5, 999.5);
  fHistNESDtracks->GetXaxis()->SetTitle("n. ESD tracks");
  fHistNESDtracks->GetYaxis()->SetTitle("dN/d(n. ESD tracks)");
  fHistNESDtracks->SetMarkerStyle(kFullCircle);

  fHistNESDfriendtracks = new TH1F("fHistNESDfriendtracks", "Number of ESD friend tracks", 1000, -0.5, 999.5);
  fHistNESDfriendtracks->GetXaxis()->SetTitle("n. ESD friend tracks");
  fHistNESDfriendtracks->GetYaxis()->SetTitle("dN/d(n. ESD friend tracks)");
  fHistNESDfriendtracks->SetMarkerStyle(kFullCircle);

  fListOut->Add(fHistPt);
  fListOut->Add(fHistQ);
  fListOut->Add(fHistNTPCCl);
  fListOut->Add(fHistNESDtracks);
  fListOut->Add(fHistNESDfriendtracks);

  PostData(0, fListOut);

  fCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(1);
}

//________________________________________________________________________
void AliAnalysisTaskPt::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  /*
  if (fUseFriends){
    Printf("In Exec: ...We have to use the friends...");
    fESDfriend = fESD->FindFriend();
    Printf("...and we got friends = %p", fESDfriend);
    if (!fESDfriend) {
      Printf("ERROR: fESDfriend not available");
	return;
    }
  }
  */

  if (fUseFriends){
    Printf("In Exec: ...We have to use the friends...");
    Printf("...and we got friends = %p", fESDfriend);
    if (!fESDfriend) {
      Printf("ERROR: fESDfriend not available");
	return;
    }
  }

  Int_t nESDtracks = fESD->GetNumberOfTracks();
  Int_t nESDfriendtracks = 0;
  if (fUseFriends) nESDfriendtracks = fESDfriend->GetNumberOfTracks();
  Printf("There are %d tracks in this event", nESDtracks);
  Printf("... and there are %d friends in this event", nESDfriendtracks);

  fHistNESDtracks->Fill(nESDtracks);
  fHistNESDfriendtracks->Fill(nESDfriendtracks);

  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < nESDtracks; iTracks++) {
    Printf("Checking track %d", iTracks);
    const AliVTrack* track = dynamic_cast<AliVTrack*>(fESD->GetTrack(iTracks));
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    Printf("track %d has pt = %f", iTracks, track->Pt());
    fHistPt->Fill(track->Pt());
    fHistNTPCCl->Fill(track->GetTPCNcls());
  } //track loop 


  if (fUseFriends){
    Printf("In the loop over the friends");
    // Friend Track loop
    for (Int_t iFriend = 0; iFriend < nESDfriendtracks; iFriend++) {
      Printf("Getting friend %d", iFriend);
      const AliVfriendTrack* friendTrack = fESDfriend->GetTrack(iFriend);
      if (!friendTrack) {
	Printf("ERROR: Could not receive track %d", iFriend);
	continue;
      }
      else {
	Printf("friend track = %p", friendTrack);
      }
      
      AliTPCseed seed;
      Int_t err = friendTrack->GetTPCseed( seed );
      Printf("err = %d", err);
      if( err==0 ){
	Printf("Found TPC seed" );
	for (Int_t irow = 0; irow < kMaxRow; irow++){
	  AliTPCclusterMI* cluMI = seed.GetClusterPointer(irow);
	  if (cluMI){
	    Printf("Found cluster at row %d", irow);
	    Float_t q = cluMI->GetQ();
	    Printf("Q = %f", q);
	    fHistQ->Fill(q);
	  }
	  else {
	    Printf("Row %d does not contain clusters", irow);
	  }      	 
	}
      }	    
      else {
	//Printf("Schade... seed is %p", seed);
      }
    }
  }

  // Post output data.
  PostData(0, fListOut);
  fEv++;
}      

//________________________________________________________________________
void AliAnalysisTaskPt::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  Printf("Terminate called: fESD = %p", fESD);

  fListOut = dynamic_cast<TList*> (GetOutputData(0)); 

  if (fListOut) {
    fHistPt = dynamic_cast<TH1F*>(fListOut->FindObject("fHistPt")); 
    if (!fHistPt) {
      Printf("ERROR: fHistPt not available");
      return;
    }
   
    TCanvas *c1 = new TCanvas("AliAnalysisTaskPt","Pt",10,10,510,510);
    c1->cd(1)->SetLogy();
    fHistPt->DrawCopy("E");
  }
  else {
    Printf("In Terminate: no TList found");
  }

}
