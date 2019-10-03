//-------------------------------------------------------------------------
//
// This is the PROOF-enabled version of TOF/AliTOFComparison.C macro.
// Origin:  Andrei.Zalite@cern.ch
//
//-------------------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TText.h"
#include "TFile.h"

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliTrackReference.h"
#include "AliMCComparisonTrack.h"

#include "AliTOFComparisonTask.h"

ClassImp(AliTOFComparisonTask)

extern TStyle *gStyle;

AliTOFComparisonTask::AliTOFComparisonTask()
  : AliAnalysisTaskSE("AliTOFComparisonTask"),
    fListOfHistos(0),
    fGood(0),
    fFound(0),
    fFake(0),
    fGoodPhi(0),
    fFoundPhi(0),
    fGoodl(0),
    fFakel(0),
    fFoundl(0)
{
  // Default constructor
  AliInfo("Defaulut Constructor AliTOFComparisonTask");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}

AliTOFComparisonTask::AliTOFComparisonTask(const char* name)
  : AliAnalysisTaskSE(name),
    fListOfHistos(0),
    fGood(0),
    fFound(0),
    fFake(0),
    fGoodPhi(0),
    fFoundPhi(0),
    fGoodl(0),
    fFakel(0),
    fFoundl(0)
{
  // Constructor
  AliInfo("Constructor AliTOFComparisonTask");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}



void AliTOFComparisonTask::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  AliInfo("AliTOFComparisonTask::UserCreateOutputObjects");
  // Create output container
  fListOfHistos = new TList();
  
  fGood = new TH1F("fGood", "Pt for good tracks", 34, 0.2, 7.0);
  fFound = new TH1F("fFound", "Pt for found tracks", 34, 0.2, 7.0);
  fFake = new TH1F("fFake",  "Pt for fake tracks", 34, 0.2, 7.0);
  fGoodPhi = new TH1F("fGoodPhi", "Phi for Good tracks", 90, 0., 2.*TMath::Pi());
  fFoundPhi = new TH1F("fFoundPhi", "Phi for Found tracks", 90, 0., 2.*TMath::Pi());
  fGoodl = new TH1F("fGoodl", "Good tracks", 30, -1., 1.);
  fFakel = new TH1F("fFakel", "Mismatched tracks", 30, -1., 1.);
  fFoundl = new TH1F("fFoundl", "Matched tracks", 30, -1., 1.);
  
  fListOfHistos->Add(fGood);
  fListOfHistos->Add(fFound);
  fListOfHistos->Add(fFake);
  fListOfHistos->Add(fGoodPhi);
  fListOfHistos->Add(fFoundPhi);
  fListOfHistos->Add(fGoodl);
  fListOfHistos->Add(fFakel);
  fListOfHistos->Add(fFoundl);
}


void AliTOFComparisonTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // MC information
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
 
  Int_t nt = 0;
    
  TClonesArray dummy("AliMCComparisonTrack",1000), *refs=&dummy;
    
  // Loop over all MC tracks and select "good" tracks
  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
    AliMCParticle* track = (AliMCParticle*)mcEvent->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d (mc loop)", iTracks);
      continue;
    }
    
    // Track selection
    if (track->Pt() < 0.2) continue; 
    if (track->Pt() > 7.0) continue;   
    if (TMath::Abs(track->Pz()/track->Pt()) > 0.999) continue;
    
    Double_t vx = track->Xv(), vy = track->Yv(), vz = track->Zv();
    if (TMath::Sqrt(vx*vx+vy*vy) > 3.5) continue;
    if (TMath::Abs(vz) > 50.) continue; 
    
    if (TMath::Abs(vx) > 0.1) continue;
    if (TMath::Abs(vy) > 0.1) continue;
	
	
    // Loop over Track References
    Bool_t labelTPC = kFALSE;
    Bool_t labelTOF = kFALSE;
    AliTrackReference* trackRef = 0;
    UInt_t iTSLayerMap = 0;
    for (Int_t iTrackRef = 0; iTrackRef  < track->GetNumberOfTrackReferences(); iTrackRef++) {
      trackRef = track->GetTrackReference(iTrackRef);
      if(trackRef) {
        Int_t detectorId = trackRef->DetectorId();
	if (detectorId == AliTrackReference::kITS) {
          Float_t Radius = trackRef->R();
	  if (Radius > 2.5 && Radius < 5.5)
            iTSLayerMap |= 1;
          else if (Radius > 5.5 && Radius < 8.5)
            iTSLayerMap |= 1 << 1;
          else if (Radius > 13.5 && Radius < 16.5)
            iTSLayerMap |= 1 << 2;
          else if (Radius > 22. && Radius < 26.)
            iTSLayerMap |= 1 << 3;
          else if (Radius > 36. && Radius < 41.)
            iTSLayerMap |= 1 << 4;
          else if (Radius > 42. && Radius < 46.)
            iTSLayerMap |= 1 << 5;
          else {
            Printf("Wrong radius %f ", Radius);
            return;
          }

        }
	if (detectorId == AliTrackReference::kTPC) {	    
	  labelTPC = kTRUE;
	}
	if (detectorId == AliTrackReference::kTOF) {
	  labelTOF = kTRUE;
	}
      }      
    }    // track references loop   
    
    // Skip tracks that passed not all ITS layers 
    if (iTSLayerMap != 0x3F) continue;
    
    // "Good" tracks
    if (labelTPC && labelTOF) {
      AliMCComparisonTrack* ref = new((*refs)[nt]) AliMCComparisonTrack();
      ref->SetLabel(iTracks);
      TParticle* particle = track->Particle();
      Int_t pdg = particle->GetPdgCode();
      ref->SetPDG(pdg); 
      ref->SetPz(track->Pz());
      Float_t Pt = track->Pt();
      ref->SetPt(Pt);
      fGood->Fill(Pt);
      Float_t phig = track->Phi();
      ref->SetPhi(phig);
      fGoodPhi->Fill(phig);
      Double_t tgl = track->Pz()/Pt;
      fGoodl->Fill(tgl);
      nt++;  
    }      
  }    // track loop 
  
  // ESD information  
  AliVEvent* event = InputEvent();
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }
    
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  
  Bool_t iFound;
  Int_t nfound = 0;
  Int_t nfake = 0;
  Int_t nlost = 0;
  
  Int_t mcGoods = refs->GetEntriesFast();
    
  // Loop over all "good" MC tracks
  for (Int_t k = 0; k < mcGoods; k++) {
    AliMCComparisonTrack* ref = (AliMCComparisonTrack*)refs->UncheckedAt(k); 
    if (!ref) continue;
    Int_t MCLabel = ref->GetLabel();
    Float_t ptg = ref->GetPt(); 
    Float_t Phig = ref->GetPhi();
    iFound = kFALSE;
    Int_t iTrack;
    AliESDtrack* track = 0;
    for (iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
      track = esd->GetTrack(iTrack);  
      if (!track) {
        Printf("ERROR: Could not receive track %d", iTrack);
        continue;
      }
      
      Int_t Label =  track->GetLabel();
      if (MCLabel != TMath::Abs(Label)) continue;
  
      if (track->GetTOFsignal() <= 0.) continue;
      
      Double_t tgl = ref->GetPz()/ref->GetPt();      
      Int_t TOFLabel[3];  
      track->GetTOFLabel(TOFLabel);            
      if (TOFLabel[0] != MCLabel)
        if (TOFLabel[1] != MCLabel)
	  if (TOFLabel[2] != MCLabel) {
	    fFake->Fill(ptg); 
	    fFakel->Fill(tgl);
	    nfake++;
	    break;
	  }
      fFound->Fill(ptg); 
      fFoundl->Fill(tgl);
      fFoundPhi->Fill(Phig);
      nfound++;	  
      break;      
    } 
    if(iTrack == esd->GetNumberOfTracks()) {
      Printf("Not matched: MCLabel %d ", MCLabel);
      nlost++;
      if (track) {
        Printf("  kITSout %d kTPCout %d kTRDout %d kTIME %d",
	       (track->GetStatus()&AliESDtrack::kITSout),
	       (track->GetStatus()&AliESDtrack::kTPCout),
	       (track->GetStatus()&AliESDtrack::kTRDout),
	       (track->GetStatus()&AliESDtrack::kTIME));
      }
      else
        Printf(" No ESD track !");
    }  
  }    
     
  Printf(" Results: " );
  Printf(" Good %d Found %d Fake %d Lost %d", nt, nfound, nfake, nlost);
    
  refs->Clear();

  // Post output data.
  PostData(1, fListOfHistos);

}      


void AliTOFComparisonTask::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query 
  fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
  if (!fListOfHistos) {
    Printf("ERROR: fListOfHistos not available");
    return;
  }

  fGood = dynamic_cast<TH1F*>(fListOfHistos->At(0));
  fFound = dynamic_cast<TH1F*>(fListOfHistos->At(1));
  fFake = dynamic_cast<TH1F*>(fListOfHistos->At(2));  
  fGoodPhi = dynamic_cast<TH1F*>(fListOfHistos->At(3));
  fFoundPhi = dynamic_cast<TH1F*>(fListOfHistos->At(4));
  fGoodl = dynamic_cast<TH1F*>(fListOfHistos->At(5));
  fFakel = dynamic_cast<TH1F*>(fListOfHistos->At(6));
  fFoundl = dynamic_cast<TH1F*>(fListOfHistos->At(7));
  
  gStyle->SetOptStat(111110);
  gStyle->SetOptFit(1); 
  
  TCanvas* c1 = new TCanvas("c1","",0,0,600,900);
  
  TPad* p1 = new TPad("p1", "", 0., 0.5, 1., 1.); p1->Draw();  
  p1->cd(); p1->SetFillColor(42); p1->SetFrameFillColor(10);
  
  fFound->Sumw2(); fGood->Sumw2(); fFake->Sumw2();
  TH1F* hgp = new TH1F("hgp", " ", 34, 0.2, 7.0);
  TH1F* hfp = new TH1F("hfp", "Probability of mismatching", 34, 0.2, 7.0);
  hgp->Divide(fFound, fGood, 1., 1., "b");
  hfp->Divide(fFake, fGood, 1., 1., "b");
  hgp->SetLineColor(4); hgp->SetLineWidth(2);
  hgp->SetMaximum(1.4);
  hgp->SetYTitle("Tracking efficiency");
  hgp->SetXTitle("Pt (GeV/c)");
  hfp->SetFillColor(1); hfp->SetFillStyle(3013); hfp->SetLineWidth(2);
  hfp->SetLineColor(2);
  
  hgp->Draw();
  hfp->Draw("histsame");
  TLine *line1 = new TLine(0.2,1.0,7.0,1.0); line1->SetLineStyle(4);
  line1->Draw("same");
  TLine *line2 = new TLine(0.2,0.9,7.0,0.9); line2->SetLineStyle(4);
  line2->Draw("same");
  c1->cd();
  
  TPad* p2 = new TPad("p2", "", 0., 0., 1., 0.5); p2->Draw();
  p2->cd(); p2->SetFillColor(42); p2->SetFrameFillColor(10);  
  
  fFoundl->Sumw2(); fGoodl->Sumw2(); fFakel->Sumw2();
  TH1F* hgl = new TH1F("hgl", "", 30, -1., 1.);
  TH1F* hfl = new TH1F("hfl", "Probability of mismatching", 30, -1., 1.);
  hgl->Divide(fFoundl, fGoodl, 1., 1., "b");
  hfl->Divide(fFakel, fGoodl, 1., 1., "b");
  hgl->SetLineColor(4); hgl->SetLineWidth(2);
  hgl->SetMaximum(1.4);
  hgl->SetYTitle("Tracking efficiency");
  hgl->SetXTitle("Tan(lambda)");
  hfl->SetFillColor(1); hfl->SetFillStyle(3013); hfl->SetLineWidth(2); 
  hfl->SetLineColor(2);
  
  hgl->Draw();
  hfl->Draw("histsame");
  TLine *line3 = new TLine(-1,1.0,1,1.0); line3->SetLineStyle(4);
  line3->Draw("same");
  TLine *line4 = new TLine(-1,0.9,1,0.9); line4->SetLineStyle(4);
  line4->Draw("same");
  
  c1->Update();
  
  TCanvas* c2 = new TCanvas("c2", "", 10, 10, 510, 510);
  c2->SetFillColor(42); c2->SetFrameFillColor(10);
  fFoundPhi->Sumw2(); fGoodPhi->Sumw2();
  TH1F* hgphi = new TH1F("hgphi", "Efficiency for good tracks (Phi)",
			 90, 0., 2.*TMath::Pi());
  hgphi->Divide(fFoundPhi, fGoodPhi, 1., 1., "B");
  hgphi->SetYTitle("Tracking efficiency");
  hgphi->SetXTitle("Phi (rad)");
  hgphi->Draw();
  
  TFile fc("AliTOFComparison.root","RECREATE");
  c1->Write();
  c2->Write();
  fc.Close();  
}
