#ifndef __CINT__
#include "TFile.h"
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TBenchmark.h"       
#include "AliTOFTrackV2.h"
#endif


Int_t AliTOFanalyzeMatchingV2Phenix() 
{
  // Read TPC and TRD reconstructed tracks and matched 
  // with TOF digits
  // origin: F. Pierella | pierella@bo.infn.it
  // Analysis done a la TOF in PHENIX

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(33);
  gStyle->SetFrameFillColor(18);
  

  gBenchmark->Start("AliTOFanalyzeMatchingV2");


  if(gAlice){
    delete gAlice;
    gAlice=0;
  }
  else{
    // Dynamically link some shared libs
    if(gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("$ALICE_ROOT/macros/loadlibs.C");
      loadlibs();
    } // end if
  }

  // ==================>histos
  TH2F *hpvsm     = new TH2F("hpvsm"," ",170,-0.5,1.2,400,-4.,4.);
  hpvsm->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  hpvsm->GetYaxis()->SetTitle("charge*momentum [GeV/c]");

  TH2F *hdEdXvsp  = new TH2F("hdEdXvsp"," ",400,-4.,4.,1000,0.,1000.);
  hdEdXvsp->GetYaxis()->SetTitle("dE/dX [au]");
  hdEdXvsp->GetXaxis()->SetTitle("charge*momentum [GeV/c]");

  TH2F *hdEdXvsabsp  = new TH2F("hdEdXvsabsp"," ",400,0.,4.,100,0.,1000.);
  hdEdXvsabsp->GetYaxis()->SetTitle("dE/dX [au]");
  hdEdXvsabsp->GetXaxis()->SetTitle("momentum [GeV/c]");

  TH2F *hinvpvstof= new TH2F("hinvpvstof"," ",300,10.,40.,400,-4.,4.);
  hinvpvstof->GetYaxis()->SetTitle("charge/momentum [(GeV/c)^{-1}]");
  hinvpvstof->GetXaxis()->SetTitle("Time of flight [ns]");

  TH2F *htofvsp   = new TH2F("htofvsp"," ",400,0.,4.,300,10.,40.);
  htofvsp->GetYaxis()->SetTitle("Time of flight [ns]");
  htofvsp->GetXaxis()->SetTitle("momentum [GeV/c]");

  TH2F *hsqmvsp   = new TH2F("hsqmvsp"," ",400,-4.,4.,130,-0.1,1.2);
  hsqmvsp->GetYaxis()->SetTitle("mass^{2} (GeV/c^{2})^{2}");
  hsqmvsp->GetXaxis()->SetTitle("charge*momentum [GeV/c]");


  TH2F *hinvbgvsp   = new TH2F("hinvbgvsp"," ",400,-4.,4.,50,0.,5);
  hinvbgvsp->GetXaxis()->SetTitle("charge*momentum [GeV/c]");
  hinvbgvsp->GetYaxis()->SetTitle("1/(#beta #gamma)");

  TH2F *hpvstofde     = new TH2F("hpvstofde"," ",350,-5,30.,400,-4.,4.);
  hpvstofde->GetXaxis()->SetTitle("Time-of-flight difference from electron [ns]");
  hpvstofde->GetYaxis()->SetTitle("charge*momentum [GeV/c]");


  TH2F *hinvpvstofde     = new TH2F("hinvpvstofde"," ",350,-5,30.,400,-4.,4.);
  hinvpvstofde->GetXaxis()->SetTitle("Time-of-flight difference from electron [ns]");
  hinvpvstofde->GetYaxis()->SetTitle("charge/momentum [(GeV/c)^{-1}]");
  // =======================================> done



  // chain in the case of more events
  TChain ch("T");
  ch.Add("tofTracks*.root"); // use wildcards

  TClonesArray *arr = new TClonesArray("AliTOFTrackV2"); 
  T->GetBranch("tracks")->SetAutoDelete(kFALSE);               
  T->SetBranchAddress("tracks",&arr);                  
  Int_t nentries = (Int_t)(T->GetEntries());           
  cout << "number of events " << nentries << endl;
  for (Int_t ev=0;ev<nentries;ev++) {
    arr->Clear();
    T->GetEntry(ev);
    Int_t ntracks = arr->GetEntriesFast();   
    cout << ntracks << endl;
    for (Int_t i=0;i<ntracks;i++) {
      AliTOFTrackV2 *toftrack = (AliTOFTrackV2*)arr->At(i);
      if(toftrack->GetTrackLabel()<0) continue; // reject fake

      Float_t momTPC=toftrack->GetPTPC(); // signed momentum
      if(toftrack->GetPdgCode()<0) momTPC=-momTPC;
      Int_t pdgCode=toftrack->GetPdgCode();
      //if(TMath::Abs(pdgCode)!=211) continue;
      hdEdXvsp->Fill(momTPC,toftrack->GetdEdX());
      hdEdXvsabsp->Fill(TMath::Abs(momTPC),toftrack->GetdEdX());

      Int_t matchStatus=toftrack->GetMatchingStatus();
      Float_t tof=toftrack->GetTof();
      if(matchStatus==3 || matchStatus==4){ // track has a tof
	htofvsp->Fill(TMath::Abs(momTPC),tof);
	hinvpvstof->Fill(tof,1./momTPC);
      } // if(matchStatus==3 || matchStatus==4)

      Float_t length=toftrack->GetLength();
      // starting mass calculation/ only possible with track length
      if(length>0. && (matchStatus==3 || matchStatus==4)) // to be skipped when track length will be in the Kalman 
	{
	  Float_t squareMass=momTPC*momTPC*((29.9792*tof/length)*(29.9792*tof/length)-1);
	  hsqmvsp->Fill(momTPC,TMath::Abs(squareMass));
	  Float_t dummy=squareMass;
	  if(dummy<0) dummy=-dummy;
	  Float_t mass=TMath::Sqrt(dummy);
	  if(squareMass<0) mass=-mass;
	  hpvsm->Fill(mass,momTPC);
	  Float_t invbg=TMath::Sqrt(dummy)/TMath::Abs(momTPC);
	  hinvbgvsp->Fill(momTPC,invbg);
	  Float_t betaEl=TMath::Abs(momTPC)/TMath::Sqrt(momTPC*momTPC+0.000510999*0.000510999);
	  Float_t tofEl=length/(betaEl*29.9792); // [ns]
	  hpvstofde->Fill((tof-tofEl),momTPC);
	  hinvpvstofde->Fill((tof-tofEl),1./momTPC);
	} // if(toftrack->GetLength()>0.)

    } // end loop on tracks
  } // end loop on 'events'
  
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);

  //TPaveLabel pl;

  TCanvas* c1 = new TCanvas("c1", "TOF a la PHENIX (I)",210, 210, 740, 690);
  c1->SetFillColor(10);
  c1->cd();
  gPad->SetGrid();
  c1->SetHighLightColor(2);
  c1->SetFillColor(10);
  c1->SetBorderSize(2);

  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hpvsm->Draw("contz");
  
  //Float_t x1=0.67, y1=0.875, x2=0.85, y2=0.95;
  /*hpvsm->Draw("contz");*/ 
  //pl.DrawPaveLabel(x1,y1,x2,y2,"CONTZ","brNDC");
  
  

  c1->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdEdXvsp->Draw();
  c1->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  hinvpvstof->Draw();
  c1->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  htofvsp->Draw();

  TCanvas* c2 = new TCanvas("c2", "TOF a la PHENIX (II)", 210, 210, 740, 690);
  c2->SetFillColor(10);
  c2->cd();

  //gStyle->SetOptStat(0);
  c2->SetHighLightColor(2);
  c2->SetFillColor(10);
  c2->SetBorderSize(2);
  c2->SetGridy();
  c2->Divide(2,2);


  c2->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hsqmvsp->Draw("contz");

  c2->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hinvbgvsp->Draw();


  c2->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  hdEdXvsabsp->Draw();


  c2->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  hpvstofde->Draw();


  TCanvas* c2b = new TCanvas("c2b", " ", 210, 210, 740, 690);
  c2b->SetFillColor(10);
  c2b->cd();
  //gStyle->SetOptStat(0);
  c2b->SetHighLightColor(2);
  c2b->SetFillColor(10);
  c2b->SetBorderSize(2);
  c2b->SetGridy();
  c2b->Divide(2,2);
  c2b->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hinvpvstofde->SetFillColor(33);
  hinvpvstofde->SetFillStyle(3001);
  hinvpvstofde->Draw();
  c2b->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hpvstofde->Draw();



  TCanvas* c3 = new TCanvas("c3", "TOF a la PHENIX (II)", 210, 210, 740, 690);
  c3->SetFillColor(10);
  c3->cd();

  //gStyle->SetOptStat(0);
  c3->SetHighLightColor(2);
  c3->SetFillColor(10);
  c3->SetBorderSize(2);
  c3->SetGridy();
  c3->Divide(2,2);
  c3->cd(1); 
  gPad->SetGridx();
  gPad->SetGridy();
  hpvsm->Draw("contz");

  c3->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hsqmvsp->Draw("contz");

  c3->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  hdEdXvsabsp->Draw();

}
