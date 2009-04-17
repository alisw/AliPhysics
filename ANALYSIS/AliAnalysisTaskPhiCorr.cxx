
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
    fMixedEvent.Reset();
    
    for (Int_t iev = 0; iev < nev; iev++) {
	AliAODEvent* aod = (AliAODEvent*) GetEvent(iev);
	fMixedEvent.AddEvent(aod);
    }
    fMixedEvent.Init();
    Int_t ntrack = fMixedEvent.GetNumberOfTracks();
    if (ntrack > 1) {
      for (Int_t itr = 0; itr < ntrack -1; itr++) {
	for (Int_t jtr = itr+1; jtr < ntrack; jtr++) {
	  AliVParticle* track1 = fMixedEvent.GetTrack(itr);
	  AliVParticle* track2 = fMixedEvent.GetTrack(jtr);
	  Int_t iev1 =  fMixedEvent.EventIndex(itr);
	  Int_t iev2 =  fMixedEvent.EventIndex(jtr);

	  Float_t phi1 = track1->Phi();
	  Float_t phi2 = track2->Phi();
	  Float_t dphi = TMath::Abs(phi1 - phi2);
	  if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
	  if (iev1 != iev2) {
	    fHistDphiUC->Fill(dphi, wgt);
	  } else {
	    fHistDphiCO->Fill(dphi, 0.5);	    
	  }
	} // tarcks
      } // tracks
    } // more than 1 
  PostData(1, fHists);
}      

//________________________________________________________________________
void AliAnalysisTaskPhiCorr::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  TCanvas *c1 = new TCanvas("AliAnalysisTaskPt","Pt",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistDphiCO->DrawCopy("E");
  fHistDphiUC->DrawCopy("Esame");
}
