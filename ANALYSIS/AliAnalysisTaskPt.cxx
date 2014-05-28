#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"

#include "AliAnalysisTaskPt.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

ClassImp(AliAnalysisTaskPt)

//________________________________________________________________________
AliAnalysisTaskPt::AliAnalysisTaskPt(const char *name) 
: AliAnalysisTask(name, ""), fESD(0), fHistPt(0), fCuts(0), fEv(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TH1F::Class());
}

//________________________________________________________________________
void AliAnalysisTaskPt::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  printf("AliAnalysisTaskPt::ConnectInputData\n");
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

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskPt::CreateOutputObjects()
{
  // Create histograms
  // Called once


  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 3.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);

  fCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(1);
  PostData(0, fHistPt);
}

//________________________________________________________________________
void AliAnalysisTaskPt::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

   AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

   if (!esdH) {
     Printf("ERROR: Could not get ESDInputHandler");
   } else
     fESD = esdH->GetEvent();
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());

  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }

    fHistPt->Fill(track->Pt());
  } //track loop 

  TObjArray* tmp = fCuts->GetAcceptedTracks(fESD);
  Printf("Event %d has %d accepted tracks", fEv, tmp->GetEntries());
  // Post output data.
  PostData(0, fHistPt);
  fEv++;
}      

//________________________________________________________________________
void AliAnalysisTaskPt::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  Printf("Terminate called");

  fHistPt = dynamic_cast<TH1F*> (GetOutputData(0));
  if (!fHistPt) {
    Printf("ERROR: fHistPt not available");
    return;
  }
   
  TCanvas *c1 = new TCanvas("AliAnalysisTaskPt","Pt",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistPt->DrawCopy("E");
}
