//-------------------------------------------------------------------------
//
// This is the PROOF-enabled version of TRD/Macros/AliTRDComparisonV2.C macro.
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

#include "AliTRDComparisonTask.h"

ClassImp(AliTRDComparisonTask)

extern TStyle *gStyle;

AliTRDComparisonTask::AliTRDComparisonTask()
  : AliAnalysisTaskSE("AliTRDComparisonTask"),
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
    fFoundPhi(0),
    fZ(0),
    fC(0)
{
  // Default constructor
  AliInfo("Default constructor AliTRDComparisonTask");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}

AliTRDComparisonTask::AliTRDComparisonTask(const char* name)
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
    fFoundPhi(0),
    fZ(0),
    fC(0)
{
  // Constructor
  AliInfo("Constructor AliTRDComparisonTask");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}



void AliTRDComparisonTask::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  AliInfo("AliTRDComparisonTask::UserCreateOutputObjects");
  // Create output container
  fListOfHistos = new TList();
  
  fGood = new TH1F("fGood", "Pt for good tracks", 34, 0.2, 7.0);
  fFound = new TH1F("fFound", "Pt for found tracks", 34, 0.2, 7.0);
  fFake = new TH1F("fFake", "Pt for fake tracks", 34, 0.2, 7.0);
  fP = new TH1F("fP", "PHI resolution", 50, -20., 20.);
  fL = new TH1F("fL", "LAMBDA resolution", 50, -20., 20.);
  fPt = new TH1F("fPt", "Relative Pt resolution", 30, -10., 10.);
  fHmpt = new TH1F("fHmpt", "Y and Z resolution", 30, -30., 30.);
  fE = new TH1F("fE", "dE/dX for pions with 0.4<p<0.5 GeV/c", 50, 0., 1000.);
  fEp = new TH2F("fEp", "dE/dX vs momentum", 50, 0., 2., 50, 0., 2000.);
  fGoodPhi = new TH1F("fGoodPhi", "Phi for Good tracks", 90, 0., 2.*TMath::Pi());
  fFoundPhi = new TH1F("fFoundPhi", "Phi for Found tracks", 90, 0., 2.*TMath::Pi());
  fZ = new TH1F("fZ", "Z resolution", 30, -30., 30.);
  fC = new TH1F("fC", "Number of the assigned clusters", 25, 110., 135.);
  
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
  fListOfHistos->Add(fZ);
  fListOfHistos->Add(fC);  
}


void AliTRDComparisonTask::UserExec(Option_t *) 
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
    
    // Check TRD information       
    Int_t nRefs = track->GetNumberOfTrackReferences();
    if (nRefs <= 1) continue;

    AliTrackReference* trackRef0 = 0;
    for (Int_t iTrackRef = 0; iTrackRef < nRefs; iTrackRef++) {
      trackRef0 = track->GetTrackReference(iTrackRef);
      if(trackRef0) {
        Int_t detectorId = trackRef0->DetectorId();
        if (detectorId == AliTrackReference::kTRD) {
	  break;
        }
	trackRef0 = 0;
      }
    }
    if (!trackRef0) continue;

    if (trackRef0->LocalX() > 300.) continue;
    
    AliTrackReference* trackRef = 0;
    for (Int_t iTrackRef = nRefs - 1; iTrackRef >= 0; --iTrackRef) {
      trackRef = track->GetTrackReference(iTrackRef);
      if (!trackRef) continue;
      if (trackRef->LocalX() > 363. &&
	  trackRef->DetectorId() == AliTrackReference::kTRD) break;       
      trackRef = 0;	       
    }
    if (!trackRef) continue;
     
    if (TMath::Abs(trackRef0->Alpha() - trackRef->Alpha()) > 1e-5) continue;
    
    // Check TPC information
    Bool_t labelTPC = kFALSE;
    for (Int_t iTrackRef = 0; iTrackRef  < nRefs; iTrackRef++) {
      trackRef = track->GetTrackReference(iTrackRef);
      if(trackRef) {
	Int_t detectorId = trackRef->DetectorId();
	if (detectorId == AliTrackReference::kTPC) {	    
	  labelTPC = kTRUE;
	  break;
	}
      }      
    }    // track references loop   
    
    // "Good" tracks
    if (labelTPC) {
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
      ref->SetLocalX(trackRef->LocalX());
      ref->SetLocalY(trackRef->LocalY());
      ref->SetZ(trackRef->Z());
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
    Int_t mcLabel = ref->GetLabel();
    Float_t ptg = ref->GetPt(); 
    Float_t phiG = ref->GetPhi();
    iFound = kFALSE;
    for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
      AliESDtrack* track = esd->GetTrack(iTrack);  
      if (!track) {
        Printf("ERROR: Could not receive track %d", iTrack);
        continue;
      }
      
      if (! track->IsOn(AliESDtrack::kTRDout)) continue;
      const AliExternalTrackParam * outorig = track->GetOuterParam();
      if (!outorig) continue;
      
      Int_t lable =  track->GetTRDLabel();      

      if (mcLabel == TMath::Abs(lable)) {	  
        if (mcLabel == lable) {
	  nfound++;
	  fFound->Fill(ptg);
	  fFoundPhi->Fill(phiG);
	} 
	else {
	  nfake++;
	  fFake->Fill(ptg);
	}
	iFound = kTRUE;
	
	AliExternalTrackParam out(*outorig);

        Double_t bz = esd->GetMagneticField();
        Double_t xg = ref->GetLocalX();
        out.PropagateTo(xg,bz);
	
        Float_t phi = TMath::ASin(out.GetSnp()) + out.GetAlpha();
        if (phi < 0) phi += 2*TMath::Pi();
        if (phi >= 2*TMath::Pi()) phi -= 2*TMath::Pi();
        Float_t phig = ref->GetPhi();
        fP->Fill((phi - phig)*1000.);
	
        Float_t lam = TMath::ATan(out.GetTgl()); 
        Float_t lamg = TMath::ATan2(ref->GetPz(), ptg);
        fL->Fill((lam - lamg)*1000.);
	
        fC->Fill(track->GetTRDclusters(0));
	
        Float_t pt1 = out.OneOverPt();
        fPt->Fill((pt1 - 1/ptg)/(1/ptg)*100.);

        Float_t y = out.GetY();
        Float_t yg = ref->GetLocalY();
        fHmpt->Fill(10*(y - yg));

        Float_t z = out.GetZ();
        Float_t zg = ref->GetZ();
        fZ->Fill(10.*(z - zg));

        Float_t mom = 1./(pt1*TMath::Cos(lam));
        Float_t dedx = track->GetTRDsignal();
	Printf (" dedx = %f ", dedx);
        fEp->Fill(mom, dedx, 1.);
	
        Int_t pdg = ref->GetPDG();
        if (TMath::Abs(pdg)==211) //pions
	  if (mom>0.4 && mom<0.5) fE->Fill(dedx,1.);
	   
        break; 
      }
    }  
    if (!iFound) {
      nlost++;
    } 
  }
     
  Printf(" Results: " );
  Printf(" Good %d Found %d Fake %d Lost %d ", nt, nfound, nfake, nlost);
    
  refs->Clear();

  // Post output data.
  PostData(1, fListOfHistos);  
}      


void AliTRDComparisonTask::Terminate(Option_t *) {
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
  fZ = dynamic_cast<TH1F*>(fListOfHistos->At(11));
  fC = dynamic_cast<TH1F*>(fListOfHistos->At(12));

  gStyle->SetOptStat(111110);
  gStyle->SetOptFit(1);
  
  TCanvas* c1 = new TCanvas("c1", "", 0, 0, 700, 850);

  Int_t minc = 33;
  
  TPad* p1 = new TPad("p1", "", 0., 0.3, 0.5, 0.6); p1->Draw();
  p1->cd(); p1->SetFillColor(42); p1->SetFrameFillColor(10);
  fP->SetFillColor(4); fP->SetXTitle("(mrad)");
  if (fP->GetEntries() < minc) fP->Draw(); else fP->Fit("gaus"); c1->cd();
  
  TPad* p2 = new TPad("p2", "", 0.5, 0.3, 1., 0.6); p2->Draw();
  p2->cd(); p2->SetFillColor(42); p2->SetFrameFillColor(10);
  fL->SetFillColor(4); fL->SetXTitle("(mrad)");
  if (fL->GetEntries() < minc) fL->Draw(); else fL->Fit("gaus"); c1->cd();
  
  TPad* p3 = new TPad("p3", "", 0., 0., 0.5, 0.3); p3->Draw();
  p3->cd(); p3->SetFillColor(42); p3->SetFrameFillColor(10);
  fPt->SetFillColor(2); fPt->SetXTitle("(%)");
  if (fPt->GetEntries() < minc) fPt->Draw(); else fPt->Fit("gaus"); c1->cd();
  
  TPad* p4 = new TPad("p4", "", 0.5, 0., 1., 0.3); p4->Draw();
  p4->cd(); p4->SetFillColor(42); p4->SetFrameFillColor(10);
  fHmpt->SetFillColor(6); fHmpt->SetXTitle("(mm)");
  if (fHmpt->GetEntries() < minc) fHmpt->Draw(); else fHmpt->Fit("gaus");
  fZ->Draw("same"); c1->cd();
  
  TPad* p5 = new TPad("p5", "", 0., 0.6, 1., 1.); p5->Draw(); p5->cd();
  p5->SetFillColor(41); p5->SetFrameFillColor(10);
  fFound->Sumw2(); fGood->Sumw2(); fFake->Sumw2();
  TH1F* hg = new TH1F("hg", "Efficiency for good tracks", 34, 0.2, 7.0);
  TH1F* hf = new TH1F("hf", "Efficiency for fake tracks", 34, 0.2, 7.0);
  hg->Divide(fFound,fGood,1.,1.,"B");
  hf->Divide(fFake,fGood,1.,1.,"B");
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
  
  TPad* p6 = new TPad("p6", "", 0., 0., 1., .5); p6->Draw();
  p6->cd(); p6->SetFillColor(42); p6->SetFrameFillColor(10);
  fE->SetFillColor(2); fE->SetFillStyle(3005);
  fE->SetXTitle("Arbitrary Units");
  if (fE->GetEntries() < minc) fE->Draw(); else fE->Fit("gaus"); c2->cd();

  TPad* p7 = new TPad("p7", "", 0., 0.5, 1., 1.); p7->Draw();
  p7->cd(); p7->SetFillColor(42); p7->SetFrameFillColor(10);
  fEp->SetMarkerStyle(8); fEp->SetMarkerSize(0.4);
  fEp->SetFillColor(2); fEp->SetFillStyle(3005);
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

  TFile fc("AliTRDComparison.root","RECREATE");
  c1->Write();
  c2->Write();
  c3->Write();
  fc.Close();  
}
