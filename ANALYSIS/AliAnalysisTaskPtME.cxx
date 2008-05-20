
#include "TH1F.h"
#include "TCanvas.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskPtME.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

ClassImp(AliAnalysisTaskPtME)

//________________________________________________________________________
AliAnalysisTaskPtME::AliAnalysisTaskPtME(const char *name) 
  : AliAnalysisTaskME(name), fHistPt(0)
{
  // Constructor
  // Define input and output slots here
  DefineOutput(1, TH1F::Class());
}


//________________________________________________________________________
void AliAnalysisTaskPtME::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 3.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
}

//________________________________________________________________________
void AliAnalysisTaskPtME::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

    for (Int_t jev = 0; jev < 3; jev++) {
	AliAODEvent* aod = GetEvent(jev);
	Int_t ntracks =  aod->GetNumberOfTracks();
	printf("Number of tracks %5d %5d\n", jev, ntracks);
	
	for (Int_t iTracks = 0; iTracks < aod->GetNumberOfTracks(); iTracks++) 
	{
	 	    AliAODTrack* track = aod->GetTrack(iTracks);
	    if (!track) {
		Printf("ERROR: Could not receive track %d", iTracks);
		continue;
	    }
	    
	    fHistPt->Fill(track->Pt());
	} //track loop
    } // event loop
  // Post output data.
  PostData(0, fHistPt);
}      

//________________________________________________________________________
void AliAnalysisTaskPtME::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  TCanvas *c1 = new TCanvas("AliAnalysisTaskPt","Pt",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistPt->DrawCopy("E");
}
