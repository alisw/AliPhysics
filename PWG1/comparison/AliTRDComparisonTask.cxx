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
#include "TBenchmark.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliTrackReference.h"
#include "AliMCComparisonTrack.h"

#include "AliTRDComparisonTask.h"

ClassImp(AliTRDComparisonTask)

AliTRDComparisonTask::AliTRDComparisonTask()
  : AliAnalysisTaskSE("AliTRDComparisonTask"),
    fListOfHistos(0),
    hgood(0),
    hfound(0),
    hfake(0),
    hp(0),
    hl(0),
    hpt(0),
    hmpt(0),
    he(0),
    hep(0),
    hgoodPhi(0),
    hfoundPhi(0),
    hz(0),
    hc(0)
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
    hgood(0),
    hfound(0),
    hfake(0),
    hp(0),
    hl(0),
    hpt(0),
    hmpt(0),
    he(0),
    hep(0),
    hgoodPhi(0),
    hfoundPhi(0),
    hz(0),
    hc(0)
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
  
  hgood = new TH1F("hgood", "Pt for good tracks", 34, 0.2, 7.0);
  hfound = new TH1F("hfound", "Pt for found tracks", 34, 0.2, 7.0);
  hfake = new TH1F("hfake", "Pt for fake tracks", 34, 0.2, 7.0);
  hp = new TH1F("hp", "PHI resolution", 50, -20., 20.);
  hl = new TH1F("hl", "LAMBDA resolution", 50, -20., 20.);
  hpt = new TH1F("hpt", "Relative Pt resolution", 30, -10., 10.);
  hmpt = new TH1F("hmpt", "Y and Z resolution", 30, -30., 30.);
  he = new TH1F("he", "dE/dX for pions with 0.4<p<0.5 GeV/c", 50, 0., 1000.);
  hep = new TH2F("hep", "dE/dX vs momentum", 50, 0., 2., 50, 0., 2000.);
  hgoodPhi = new TH1F("hgoodPhi", "Phi for Good tracks", 90, 0., 2.*TMath::Pi());
  hfoundPhi = new TH1F("hfoundPhi", "Phi for Found tracks", 90, 0., 2.*TMath::Pi());
  hz = new TH1F("hz", "Z resolution", 30, -30., 30.);
  hc = new TH1F("hc", "Number of the assigned clusters", 25, 110., 135.);
  
  fListOfHistos->Add(hgood);
  fListOfHistos->Add(hfound);
  fListOfHistos->Add(hfake);
  fListOfHistos->Add(hp);
  fListOfHistos->Add(hl);
  fListOfHistos->Add(hpt);
  fListOfHistos->Add(hmpt);
  fListOfHistos->Add(he);
  fListOfHistos->Add(hep);
  fListOfHistos->Add(hgoodPhi);
  fListOfHistos->Add(hfoundPhi);
  fListOfHistos->Add(hz);
  fListOfHistos->Add(hc);  
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
    AliMCParticle* track = mcEvent->GetTrack(iTracks);
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
    Bool_t LabelTPC = kFALSE;
    for (Int_t iTrackRef = 0; iTrackRef  < nRefs; iTrackRef++) {
      trackRef = track->GetTrackReference(iTrackRef);
      if(trackRef) {
	Int_t detectorId = trackRef->DetectorId();
	if (detectorId == AliTrackReference::kTPC) {	    
	  LabelTPC = kTRUE;
	  break;
	}
      }      
    }    // track references loop   
    
    // "Good" tracks
    if (LabelTPC) {
      AliMCComparisonTrack* ref = new((*refs)[nt]) AliMCComparisonTrack();
      ref->SetLabel(iTracks);
      TParticle* particle = track->Particle();
      Int_t pdg = particle->GetPdgCode();
      ref->SetPDG(pdg); 
      ref->SetPz(trackRef->Pz());
      Float_t Pt = trackRef->Pt();
      ref->SetPt(Pt);
      hgood->Fill(Pt);
      Float_t phig = trackRef->Phi();
      ref->SetPhi(phig);
      hgoodPhi->Fill(phig);
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
  
  Int_t MCgoods = refs->GetEntriesFast();
    
  // Loop over all "good" MC tracks
  for (Int_t k = 0; k < MCgoods; k++) {
    AliMCComparisonTrack* ref = (AliMCComparisonTrack*)refs->UncheckedAt(k); 
    if (!ref) continue;
    Int_t MCLabel = ref->GetLabel();
    Float_t ptg = ref->GetPt(); 
    Float_t PhiG = ref->GetPhi();
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
      
      Int_t Label =  track->GetTRDLabel();      

      if (MCLabel == TMath::Abs(Label)) {	  
        if (MCLabel == Label) {
	  nfound++;
	  hfound->Fill(ptg);
	  hfoundPhi->Fill(PhiG);
	} 
	else {
	  nfake++;
	  hfake->Fill(ptg);
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
        hp->Fill((phi - phig)*1000.);
	
        Float_t lam = TMath::ATan(out.GetTgl()); 
        Float_t lamg = TMath::ATan2(ref->GetPz(), ptg);
        hl->Fill((lam - lamg)*1000.);
	
        hc->Fill(track->GetTRDclusters(0));
	
        Float_t pt_1 = out.OneOverPt();
        hpt->Fill((pt_1 - 1/ptg)/(1/ptg)*100.);

        Float_t y = out.GetY();
        Float_t yg = ref->GetLocalY();
        hmpt->Fill(10*(y - yg));

        Float_t z = out.GetZ();
        Float_t zg = ref->GetZ();
        hz->Fill(10.*(z - zg));

        Float_t mom = 1./(pt_1*TMath::Cos(lam));
        Float_t dedx = track->GetTRDsignal();
	Printf (" dedx = %f ", dedx);
        hep->Fill(mom, dedx, 1.);
	
        Int_t pdg = ref->GetPDG();
        if (TMath::Abs(pdg)==211) //pions
	  if (mom>0.4 && mom<0.5) he->Fill(dedx,1.);
	   
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

  hgood = dynamic_cast<TH1F*>(fListOfHistos->At(0));
  hfound = dynamic_cast<TH1F*>(fListOfHistos->At(1));
  hfake = dynamic_cast<TH1F*>(fListOfHistos->At(2));  
  hp = dynamic_cast<TH1F*>(fListOfHistos->At(3));
  hl = dynamic_cast<TH1F*>(fListOfHistos->At(4));
  hpt = dynamic_cast<TH1F*>(fListOfHistos->At(5));
  hmpt = dynamic_cast<TH1F*>(fListOfHistos->At(6));
  he = dynamic_cast<TH1F*>(fListOfHistos->At(7));
  hep = dynamic_cast<TH2F*>(fListOfHistos->At(8));
  hgoodPhi = dynamic_cast<TH1F*>(fListOfHistos->At(9));
  hfoundPhi = dynamic_cast<TH1F*>(fListOfHistos->At(10));
  hz = dynamic_cast<TH1F*>(fListOfHistos->At(11));
  hc = dynamic_cast<TH1F*>(fListOfHistos->At(12));

  gStyle->SetOptStat(111110);
  gStyle->SetOptFit(1);
  
  TCanvas* c1 = new TCanvas("c1", "", 0, 0, 700, 850);

  Int_t minc = 33;
  
  TPad* p1 = new TPad("p1", "", 0., 0.3, 0.5, 0.6); p1->Draw();
  p1->cd(); p1->SetFillColor(42); p1->SetFrameFillColor(10);
  hp->SetFillColor(4); hp->SetXTitle("(mrad)");
  if (hp->GetEntries() < minc) hp->Draw(); else hp->Fit("gaus"); c1->cd();
  
  TPad* p2 = new TPad("p2", "", 0.5, 0.3, 1., 0.6); p2->Draw();
  p2->cd(); p2->SetFillColor(42); p2->SetFrameFillColor(10);
  hl->SetFillColor(4); hl->SetXTitle("(mrad)");
  if (hl->GetEntries() < minc) hl->Draw(); else hl->Fit("gaus"); c1->cd();
  
  TPad* p3 = new TPad("p3", "", 0., 0., 0.5, 0.3); p3->Draw();
  p3->cd(); p3->SetFillColor(42); p3->SetFrameFillColor(10);
  hpt->SetFillColor(2); hpt->SetXTitle("(%)");
  if (hpt->GetEntries() < minc) hpt->Draw(); else hpt->Fit("gaus"); c1->cd();
  
  TPad* p4 = new TPad("p4", "", 0.5, 0., 1., 0.3); p4->Draw();
  p4->cd(); p4->SetFillColor(42); p4->SetFrameFillColor(10);
  hmpt->SetFillColor(6); hmpt->SetXTitle("(mm)");
  if (hmpt->GetEntries() < minc) hmpt->Draw(); else hmpt->Fit("gaus");
  hz->Draw("same"); c1->cd();
  
  TPad* p5 = new TPad("p5", "", 0., 0.6, 1., 1.); p5->Draw(); p5->cd();
  p5->SetFillColor(41); p5->SetFrameFillColor(10);
  hfound->Sumw2(); hgood->Sumw2(); hfake->Sumw2();
  TH1F* hg = new TH1F("hg", "Efficiency for good tracks", 34, 0.2, 7.0);
  TH1F* hf = new TH1F("hf", "Efficiency for fake tracks", 34, 0.2, 7.0);
  hg->Divide(hfound,hgood,1.,1.,"B");
  hf->Divide(hfake,hgood,1.,1.,"B");
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
  he->SetFillColor(2); he->SetFillStyle(3005);
  he->SetXTitle("Arbitrary Units");
  if (he->GetEntries() < minc) he->Draw(); else he->Fit("gaus"); c2->cd();

  TPad* p7 = new TPad("p7", "", 0., 0.5, 1., 1.); p7->Draw();
  p7->cd(); p7->SetFillColor(42); p7->SetFrameFillColor(10);
  hep->SetMarkerStyle(8); hep->SetMarkerSize(0.4);
  hep->SetFillColor(2); hep->SetFillStyle(3005);
  hep->SetXTitle("p (GeV/c)"); hep->SetYTitle("dE/dX (Arb. Units)");
  hep->Draw(); c1->cd();
   
  TCanvas* c3 = new TCanvas("c3", "", 10, 10, 510, 510);
  c3->SetFillColor(42); c3->SetFrameFillColor(10);
  hfoundPhi->Sumw2(); hgoodPhi->Sumw2();
  TH1F* hgphi = new TH1F("hgphi", "Efficiency for good tracks (Phi)", 
			  90, 0., 2.*TMath::Pi());
  hgphi->Divide(hfoundPhi, hgoodPhi, 1., 1., "B"); 
  hgphi->SetYTitle("Tracking efficiency");
  hgphi->SetXTitle("Phi (rad)");
  hgphi->Draw();

  TFile fc("AliTRDComparison.root","RECREATE");
  c1->Write();
  c2->Write();
  c3->Write();
  fc.Close();  
}
