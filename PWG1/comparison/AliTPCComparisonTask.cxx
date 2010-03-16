//-------------------------------------------------------------------------
//
// This is the PROOF-enabled version of TPC/AliTPCComparison.C macro.
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
#include "AliMCComparisonTrack.cxx"

#include "AliTPCComparisonTask.h"

ClassImp(AliTPCComparisonTask)


AliTPCComparisonTask::AliTPCComparisonTask()
  : AliAnalysisTaskSE("AliTPCComparisonTask"),
    fListOfHistos(0),
    fGood(0),
    fFound(0),
    fFake(0),
    fP(0),
    fL(0),
    fPt(0),
    fHmpt(0),
    fE(0),
    fEp(0),
    fGoodPhi(0),
    fFoundPhi(0)
{
  // Default constructor
  AliInfo("Default constructor AliTPCComparisonTask");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}


AliTPCComparisonTask::AliTPCComparisonTask(const char* name)
  : AliAnalysisTaskSE(name),
    fListOfHistos(0),
    fGood(0),
    fFound(0),
    fFake(0),
    fP(0),
    fL(0),
    fPt(0),
    fHmpt(0),
    fE(0),
    fEp(0),
    fGoodPhi(0),
    fFoundPhi(0)
{
  // Constructor
  AliInfo("Constructor AliTPCComparisonTask");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}


void AliTPCComparisonTask::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  AliInfo("AliTPCComparisonTask::UserCreateOutputObjects");
  // Create output container
  fListOfHistos = new TList();
  
  fGood = new TH1F("fGood", "Pt for good tracks", 34, 0.2, 7.0);
  fFound = new TH1F("fFound", "Pt for found tracks", 34, 0.2, 7.0);
  fFake = new TH1F("fFake",  "Pt for fake tracks", 34, 0.2, 7.0);
  fP = new TH1F("fP", "PHI resolution", 50, -20., 20.);
  fL = new TH1F("fL", "LAMBDA resolution", 50, -20., 20.);
  fPt = new TH1F("fPt", "Relative Pt resolution", 30, -10., 10.);
  fHmpt = new TH1F("fHmpt", "Relative Pt resolution (pt>4GeV/c)", 30, -60., 60.);
  fE = new TH1F("fE", "dE/dX for pions with 0.4<p<0.5 GeV/c", 50, 0., 100.);
  fEp = new TH2F("fEp", "dE/dX vs momentum", 50, 0., 2., 50, 0., 400.);
  fGoodPhi = new TH1F("fGoodPhi", "Phi for good tracks", 90, 0., 2.*TMath::Pi());
  fFoundPhi = new TH1F("fFoundPhi", "Phi for found tracks", 90, 0., 2.*TMath::Pi());
  
  fListOfHistos->Add(fGood);
  fListOfHistos->Add(fFound);
  fListOfHistos->Add(fFake);
  fListOfHistos->Add(fP);
  fListOfHistos->Add(fL);
  fListOfHistos->Add(fPt);
  fListOfHistos->Add(fHmpt);
  fListOfHistos->Add(fE);
  fListOfHistos->Add(fEp);
  fListOfHistos->Add(fGoodPhi);
  fListOfHistos->Add(fFoundPhi);
}


void AliTPCComparisonTask::UserExec(Option_t *) 
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
    
    Double_t vx = track->Xv(),vy = track->Yv(),vz = track->Zv();
    if (TMath::Sqrt(vx*vx+vy*vy) > 3.5) continue;
    if (TMath::Abs(vz) > 50.) continue; 
    
    // Loop over Track References
    Bool_t LabelTPC = kFALSE;
    AliTrackReference* trackRef = 0;
    for (Int_t iTrackRef = 0; iTrackRef  < track->GetNumberOfTrackReferences(); iTrackRef++) {
      trackRef = track->GetTrackReference(iTrackRef);
      if(trackRef) {
	Int_t detectorId = trackRef->DetectorId(); 
	if (detectorId == AliTrackReference::kTPC) {	    
	  LabelTPC = kTRUE;
	  break;
	}
      }      
    }
    // "Good" tracks
    if (LabelTPC) {
      AliMCComparisonTrack* ref = new((*refs)[nt]) AliMCComparisonTrack();
      ref->SetLabel(iTracks);
      TParticle* particle = track->Particle();
      Int_t pdg = particle->GetPdgCode();
      ref->SetPDG(pdg); 
      ref->SetPz(trackRef->Pz());
      Float_t pt = trackRef->Pt();
      ref->SetPt(pt);
      fGood->Fill(pt);
      Float_t phig = trackRef->Phi();
      ref->SetPhi(phig);
      fGoodPhi->Fill(phig);
      nt++;  
    }    
  } //track loop 
  
  
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
  
  for (Int_t k = 0; k < mcGoods; k++) {
    AliMCComparisonTrack* ref = (AliMCComparisonTrack*)refs->UncheckedAt(k);
    if (!ref) continue;
    Int_t mcLabel = ref->GetLabel();
    Float_t ptg = ref->GetPt();
    Float_t Phig = ref->GetPhi();

    iFound = kFALSE;
    for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
      AliESDtrack* track = esd->GetTrack(iTrack);  
      if (!track) {
        Printf("ERROR: Could not receive track %d", iTrack);
        continue;
      }
      
      if (! track->IsOn(AliESDtrack::kTPCrefit)) continue;     
      
      Int_t tpcLabel =  track->GetTPCLabel();      
      if (mcLabel == TMath::Abs(tpcLabel)) {	  
        if (mcLabel == tpcLabel) {
	  nfound++;
	  fFound->Fill(ptg);
	  fFoundPhi->Fill(Phig);
	} 
	else {
	  nfake++;
	  fFake->Fill(ptg);
	}
	iFound = kTRUE;
	
      
        Double_t pxpypz[3]; track->GetInnerPxPyPz(pxpypz);
        Float_t phi = TMath::ATan2(pxpypz[1],pxpypz[0]);
	if (phi < 0) phi += 2*TMath::Pi();
        Double_t pt = TMath::Sqrt(pxpypz[0]*pxpypz[0]+pxpypz[1]*pxpypz[1]);
        Float_t lam = TMath::ATan2(pxpypz[2],pt); 
        Float_t pt1 = 1/pt;

        Int_t pdg = ref->GetPDG();  
        if (TMath::Abs(pdg) == 11 && ptg>4.) {
	  //high momentum electrons
          fHmpt->Fill((pt1 - 1/ptg)/(1/ptg)*100.);
        } 
	else {
          fP->Fill((phi - Phig)*1000.);          
	  Float_t lamg = TMath::ATan2(ref->GetPz(),ptg);
          fL->Fill((lam - lamg)*1000.);
          fPt->Fill((pt1 - 1/ptg)/(1/ptg)*100.);
        }
	
	Float_t mom = pt/TMath::Cos(lam);
        Float_t dedx = track->GetTPCsignal();
        fEp->Fill(mom, dedx, 1.);
        if (TMath::Abs(pdg) == 211) { // pions
	  if (mom > 0.4 && mom < 0.5) {
            fE->Fill(dedx,1.);
          }
	}   
	   
        break; 
      }
    }  
    if (!iFound) {
      nlost++;
    } 
  }
     
  Printf(" Results: " );
  Printf(" Found %d Fake %d Lost %d ", nfound, nfake, nlost);
  
  refs->Clear();
  
  // Post output data.
  PostData(1, fListOfHistos);
}      


void AliTPCComparisonTask::Terminate(Option_t *) 
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
  fP = dynamic_cast<TH1F*>(fListOfHistos->At(3));
  fL = dynamic_cast<TH1F*>(fListOfHistos->At(4));
  fPt = dynamic_cast<TH1F*>(fListOfHistos->At(5));
  fHmpt = dynamic_cast<TH1F*>(fListOfHistos->At(6));
  fE = dynamic_cast<TH1F*>(fListOfHistos->At(7));
  fEp = dynamic_cast<TH2F*>(fListOfHistos->At(8));
  fGoodPhi = dynamic_cast<TH1F*>(fListOfHistos->At(9));
  fFoundPhi = dynamic_cast<TH1F*>(fListOfHistos->At(10));

  gStyle->SetOptStat(111110);
  gStyle->SetOptFit(1);
   
  TCanvas* c1 = new TCanvas("c1", "", 0, 0, 700, 850);

  Int_t minc = 33;
  
  TPad* p1 = new TPad("p1", "", 0., 0.3, 0.5, 0.6); p1->Draw();
  p1->cd(); p1->SetFillColor(42); p1->SetFrameFillColor(10);
  fP->SetFillColor(4); fP->SetXTitle("(mrad)");
  if (fP->GetEntries() < minc) fP->Draw(); else fP->Fit("gaus"); c1->cd();

  TPad* p2 = new TPad("p2", "", 0.5, 0.3, 1, 0.6); p2->Draw();
  p2->cd(); p2->SetFillColor(42); p2->SetFrameFillColor(10);
  fL->SetFillColor(4); fL->SetXTitle("(mrad)");
  if (fL->GetEntries() < minc) fL->Draw(); else fL->Fit("gaus"); c1->cd();
  
  TPad* p3 = new TPad("p3", "", 0., 0., 0.5, 0.3); p3->Draw();
  p3->cd(); p3->SetFillColor(42); p3->SetFrameFillColor(10);
  fPt->SetFillColor(2); fPt->SetXTitle("(%)");
  if (fPt->GetEntries() < minc) fPt->Draw(); else fPt->Fit("gaus"); c1->cd();
  
  TPad* p4 = new TPad("p4", "", 0.5, 0., 1., 0.3); p4->Draw();
  p4->cd(); p4->SetFillColor(42); p4->SetFrameFillColor(10);
  fHmpt->SetFillColor(6); fHmpt->SetXTitle("(%)");
  if (fHmpt->GetEntries() < minc) fHmpt->Draw(); else fHmpt->Fit("gaus"); c1->cd(); 
  
  TPad* p5 = new TPad("p5", "", 0., 0.6, 1., 1.); p5->Draw(); p5->cd();
  p5->SetFillColor(41); p5->SetFrameFillColor(10);
  fFound->Sumw2(); fGood->Sumw2(); fFake->Sumw2();
  TH1F* hg = new TH1F("hg", "Efficiency for good tracks", 34, 0.2, 7.0);
  TH1F* hf = new TH1F("hf", "Efficiency for fake tracks", 34, 0.2, 7.0);
  hg->Divide(fFound, fGood, 1., 1., "B");
  hf->Divide(fFake, fGood, 1., 1., "B");
  hg->SetLineColor(4); hg->SetLineWidth(2);
  hg->SetMaximum(1.4);
  hg->SetYTitle("Tracking efficiency");
  hg->SetXTitle("Pt (GeV/c)");
  hg->Draw();
  
  TLine* line1 = new TLine(0.2, 1.0, 7.0, 1.0); line1->SetLineStyle(4);
  line1->Draw("same");
  TLine* line2 = new TLine(0.2, 0.9, 7.0, 0.9); line2->SetLineStyle(4);
  line2->Draw("same");
  
  hf->SetFillColor(1);
  hf->SetFillStyle(3013);
  hf->SetLineColor(2);
  hf->SetLineWidth(2);
  hf->Draw("histsame");
  TText* text = new TText(0.461176, 0.248448, "Fake tracks");
  text->SetTextSize(0.05);
  text->Draw();
  text = new TText(0.453919, 1.11408, "Good tracks");
  text->SetTextSize(0.05);
  text->Draw();
  
  TCanvas* c2 = new TCanvas("c2", "", 320, 32, 530, 590);
  
  TPad* p6 = new TPad("p6", "", 0., 0., 1., 0.5); p6->Draw();
  p6->cd(); p6->SetFillColor(42); p6->SetFrameFillColor(10);
  fE->SetFillColor(2); fE->SetFillStyle(3005);
  fE->SetXTitle("Arbitrary Units");
  if (fE->GetEntries() < minc) fE->Draw(); else fE->Fit("gaus"); c2->cd();
  
  TPad* p7 = new TPad("p7", "", 0., 0.5, 1., 1.); p7->Draw();
  p7->cd(); p7->SetFillColor(42); p7->SetFrameFillColor(10);
  fEp->SetFillColor(2); fEp->SetFillStyle(3005);
  fEp->SetMarkerStyle(8);
  fEp->SetMarkerSize(0.4);
  fEp->SetXTitle("p (GeV/c)"); fEp->SetYTitle("dE/dX (Arb. Units)");
  fEp->Draw(); c1->cd();
  
  TCanvas* c3 = new TCanvas("c3", "", 10, 10, 510, 510);
  c3->SetFillColor(42); c3->SetFrameFillColor(10);
  fFoundPhi->Sumw2(); fGoodPhi->Sumw2();
  TH1F* hgphi = new TH1F("hgphi", "Efficiency for good tracks (Phi)", 
			 90, 0., 2.*TMath::Pi());
  hgphi->Divide(fFoundPhi, fGoodPhi, 1., 1., "B"); 
  hgphi->SetYTitle("Tracking efficiency");
  hgphi->SetXTitle("Phi (rad)");
  hgphi->Draw();
  
  TFile fc("AliTPCComparison.root", "RECREATE");
  c1->Write();
  c2->Write();
  c3->Write();
  fc.Close();

}
