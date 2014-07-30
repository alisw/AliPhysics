#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliVVevent.h"
#include "AliVVtrack.h"
#include "AliESDtrackCuts.h"
#include "AliVEventHandler.h"
#include "../TPC/Rec/AliTPCseed.h"
#include "../TPC/Rec/AliTPCclusterMI.h"

#include "AliAnalysisTaskPt.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

ClassImp(AliAnalysisTaskPt)

//________________________________________________________________________
AliAnalysisTaskPt::AliAnalysisTaskPt(const char *name) 
: AliAnalysisTask(name, ""), fESD(0), fESDfriend(0), fHistPt(0), fCuts(0), fEv(0), fHistQ(0), fListOut(0)
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
  } else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
    /*
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("fTracks.*", kTRUE);
    */

    AliVEventHandler *esdH = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else {
      Printf("----> AliAnalysisTaskPt::ConnectInputData Getting the event from handler %p", esdH);
      //fESD = dynamic_cast<AliESDEvent*>(esdH->GetEvent());
      fESD = esdH->GetEvent();
      fESDfriend = esdH->GetFriendEvent();
    }
    if (!fESD) {
      Printf("ERROR, no ESD event");
    }
    if (!fESDfriend) {
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

  fListOut->Add(fHistPt);
  fListOut->Add(fHistQ);

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
  if (!fESDfriend) {
    Printf("ERROR: fESDfriend not available");
    return;
  }

  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
  Printf("... and there are %d friends in this event", fESDfriend->GetNumberOfTracks());

  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliVVtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }

    fHistPt->Fill(track->Pt());
  } //track loop 


  // Friend Track loop
  for (Int_t iFriend = 0; iFriend < fESDfriend->GetNumberOfTracks(); iFriend++) {
    AliESDfriendTrack* friendTrack = fESDfriend->GetTrack(iFriend);
    if (!friendTrack) {
      Printf("ERROR: Could not receive track %d", iFriend);
      continue;
    } 
    TObject* calibObject;
    AliTPCseed* seed = NULL;
    for (Int_t idx = 0; (calibObject = friendTrack->GetCalibObject(idx)); ++idx) {
      Printf(" |Cal %d = %p", idx, calibObject); 
      if ((seed = dynamic_cast<AliTPCseed*>(calibObject))) {
	Printf("Found TPC seed");
	for (Int_t irow = 0; irow < 160; irow++){
	  AliTPCclusterMI* cluMI = seed->GetClusterPointer(irow);
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
