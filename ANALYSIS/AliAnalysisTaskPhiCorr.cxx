
#include "TH1F.h"
#include "TList.h"
#include "TCanvas.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskPhiCorr.h"
#include "AliMultiEventInputHandler.h"

ClassImp(AliAnalysisTaskPhiCorr)

//________________________________________________________________________
AliAnalysisTaskPhiCorr::AliAnalysisTaskPhiCorr(const char *name) 
    : AliAnalysisTaskME(name), fHists(0), fHistDphiCO(0),  fHistDphiUC(0)
{
  // Constructor
  // Define input and output slots here
  DefineOutput(1, TList::Class());
}


//________________________________________________________________________
void AliAnalysisTaskPhiCorr::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fHistDphiCO = new TH1F("fHistDPhiCO", "#Delta #Phi distribution", 64, 0., 3.2);
  fHistDphiCO->GetXaxis()->SetTitle("#Delta#Phi [rad]");
  fHistDphiCO->GetYaxis()->SetTitle("dN/d#Phi");
  fHistDphiCO->SetMarkerStyle(kFullCircle);

  fHistDphiUC = new TH1F("fHistDPhiUC", "#Delta #Phi distribution", 64, 0., 3.2);
  fHistDphiUC->GetXaxis()->SetTitle("#Delta#Phi [rad]");
  fHistDphiUC->GetYaxis()->SetTitle("dN/d#Phi");
  fHistDphiUC->SetMarkerStyle(kOpenCircle);
  
  fHists = new TList();
  fHists->Add(fHistDphiCO);
  fHists->Add(fHistDphiUC);  
}

//________________________________________________________________________
void AliAnalysisTaskPhiCorr::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
    // Uncorrelated tracks
    Int_t nev = fInputHandler->GetBufferSize();
    Float_t wgt = 1./(nev*(nev-1));
    
    for (Int_t iev = 0; iev < nev; iev++) {
	for (Int_t jev = (iev + 1); jev < nev; jev++) {
	    AliAODEvent* aod1 = (AliAODEvent*)GetEvent(iev);
	    AliAODEvent* aod2 = (AliAODEvent*)GetEvent(jev);
	    
	    Int_t ntracks1 =  aod1->GetNumberOfTracks();
	    Int_t ntracks2 =  aod2->GetNumberOfTracks();
	    
	    printf("Number of tracks %5d:%5d %5d:%5d\n", iev, ntracks1, jev, ntracks2);
	    
	    for (Int_t iTracks = 0; iTracks < ntracks1; iTracks++) {
		for (Int_t jTracks = 0; jTracks < ntracks2; jTracks++) {
		    
		    
		    AliAODTrack* track1 = aod1->GetTrack(iTracks);
		    AliAODTrack* track2 = aod2->GetTrack(jTracks);
		    
		    Float_t phi1 = track1->Phi();
		    Float_t phi2 = track2->Phi();
		    Float_t dphi = TMath::Abs(phi1 - phi2);
		    if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
		    fHistDphiUC->Fill(dphi, wgt);
		} // tracks
	    } // tracks 
	} // event loop
    } // event loop
    
//    Correlated 
    AliAODEvent* aod = (AliAODEvent*)fInputHandler->GetLatestEvent();

    Int_t ntracks = aod->GetNumberOfTracks();
    printf("Number of tracks %5d: \n", ntracks);
    
    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
	for (Int_t jTracks = (iTracks+1); jTracks < ntracks; jTracks++) {
	      AliAODTrack* track1 = aod->GetTrack(iTracks);
	      AliAODTrack* track2 = aod->GetTrack(jTracks);
	      Float_t phi1 = track1->Phi();
	      Float_t phi2 = track2->Phi();
	      Float_t dphi = TMath::Abs(phi1 - phi2);
	      if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
	      fHistDphiCO->Fill(dphi);
	} // tracks
    }

  // Post output data.
  PostData(1, fHists);
}      

//________________________________________________________________________
void AliAnalysisTaskPhiCorr::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  TCanvas *c1 = new TCanvas("AliAnalysisTaskPt","corr1",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistDphiUC->DrawCopy("E");
  fHistDphiCO->DrawCopy("Esame");
}
