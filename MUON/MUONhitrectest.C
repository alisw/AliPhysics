#include "iostream.h"

void MUONhitrectest (Int_t evNumber1=0,Int_t evNumber2=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
//  Create some histograms

   TH2F *h21 = new TH2F("h21","Hits",100,-100,100,100,-100,100);
   TH2F *h22 = new TH2F("h22","CoG ",100,-100,100,100,-100,100);
   TH1F *h1 = new TH1F("h1","Multiplicity",30,-0.5,29.5);
   TH1F *hmult = new TH1F("hmult","Multiplicity",30,-0.5,29.5);
   TH1F *hresx = new TH1F("hresx","Residuals",100,-1,1);
   TH1F *hresy = new TH1F("hresy","Residuals",100,-10.,10);
   TH1F *hresym = new TH1F("hresym","Residuals",100,-500,500);
   TH1F *hpos = new TH1F("hnpos","Possibilities",10,-0.5,9.5);
   TH2F *hchi1 = new TH2F("hchi1","Chi2 vs Residuals",100,0,0.2,100,-500,500);
   TH2F *hchi2 = new TH2F("hchi2","Chi2 vs Residuals",100,0,20,100,-500,500);
//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;

   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       cout << "nev         " << nev <<endl;
       cout << "nparticles  " << nparticles <<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       
       TTree *TH = gAlice->TreeH();
       Int_t ntracks = TH->GetEntries();
       cout<<"ntracks "<<ntracks<<endl;

       Int_t nbytes = 0;



// Get pointers to Alice detectors and Digits containers
       AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
       TClonesArray *Particles = gAlice->Particles();
       TTree *TR = gAlice->TreeR();
       Int_t nent=TR->GetEntries();
       printf("Found %d entries in the tree (must be one per cathode per event! + 1empty)\n",nent);
       if (MUON) {
	   Float_t xc[2000], yc[2000], itrc[2000];
	   TClonesArray *MUONrawclust  = MUON->RawClustAddress(0);
	   nbytes += TR->GetEvent(1);
	   Int_t nrawcl = MUONrawclust->GetEntries();
	   AliMUONRawCluster  *mRaw;
	   for (Int_t iraw=0; iraw < nrawcl; iraw++) {
	       mRaw = (AliMUONRawCluster*)MUONrawclust->UncheckedAt(iraw);
	       xc[iraw]=mRaw->fX[1];
	       yc[iraw]=mRaw->fY[0];       
	       itrc[iraw]=mRaw->fTracks[1];
	   } // cluster

	   gAlice->ResetHits();
	   for (Int_t itrack=0; itrack<ntracks ;itrack++) 
	   {
	       Int_t nbytes += TH->GetEvent(itrack);
	       Int_t nhit=MUON->Hits()->GetEntriesFast();
	       for(AliMUONHit* mHit=(AliMUONHit*)MUON->FirstHit(itrack); 
		   mHit;
		   mHit=(AliMUONHit*)MUON->NextHit()) 
	       {
		   Int_t   nch   = mHit->fChamber;  // chamber number
		   Float_t x     = mHit->fX;        // x-pos of hit
		   Float_t y     = mHit->fY;        // y-pos
		   Int_t   ip    = mHit->fParticle;
		   if (ip != kMuonPlus && ip != kMuonMinus) continue;
		   if (nch != 1) continue;
		   Float_t dmin=1.e8;
		   Int_t imin=-1;
//
// Loop over clusters
		   Int_t npos=0;
		   
		   for (Int_t i=0; i<nrawcl; i++) {
		       Float_t dx = xc[i]-x;
		       Float_t dy = yc[i]-y;
		       Float_t dr  = TMath::Sqrt(dx*dx+dy*dy);
		       if (dr<dmin) {
			   imin=i;
			   dmin=dr;
		       }
		       if (TMath::Abs(dy)< 0.05 && TMath::Abs(dx)< 0.2 ) npos++;
		   }
		   mRaw = (AliMUONRawCluster*)MUONrawclust->UncheckedAt(imin);
		   Int_t   track=mRaw->fTracks[1];
		   Float_t xrec=mRaw->fX[1];
		   Float_t yrec=mRaw->fY[0];
		   Float_t x_res=xrec-x;
		   Float_t y_res=yrec-y;
		   hresx->Fill(x_res,1.);
		   hresy->Fill(y_res,1.);
		   hpos->Fill(Float_t(npos), 1.);
		   hresym->Fill(y_res*1.e4,1.);	       
		   if (npos == 0) {
		       printf("No cluster %d %f %f %d\n", 
			      itrack, x, y, nhit);
		       h21->Fill(x,y,1.);
		   }
		   
		   
		   
	       }  // hits
	   } // tracks
       }   // end if MUON
   }   // event loop 
   TCanvas *c1 = new TCanvas("c1","Charge and Residuals",400,10,600,700);
   c1->Divide(2,2);
   c1->cd(1);
   hresx->SetFillColor(42);
   hresx->SetXTitle("xrec-x");
   hresx->Draw();

   c1->cd(2);
   hresy->SetFillColor(42);
   hresy->SetXTitle("yrec-y");
   hresy->Draw();

   c1->cd(3);
   hpos->SetFillColor(42);
   hpos->SetXTitle("Possibilities");
   hpos->Draw();

   c1->cd(4);
   hresym->SetFillColor(42);
   hresym->SetXTitle("yrec-y");
   hresym->Fit("gaus");
   hresym->Draw();

   TCanvas *c2 = new TCanvas("c2","Charge and Residuals",400,10,600,700);
   c2->Divide(2,2);

   c2->cd(1);
   h21->SetFillColor(42);
   h21->SetXTitle("x");
   h21->SetYTitle("y");
   h21->Draw();

   c2->cd(2);
   hmult->SetFillColor(42);
   hmult->SetXTitle("multiplicity");
   hmult->Draw();

   c2->cd(3);
   hchi1->SetFillColor(42);
   hchi1->SetXTitle("chi2");
   hchi1->Draw();

   c2->cd(4);
   hchi2->SetFillColor(42);
   hchi2->SetXTitle("chi2");
   hchi2->Draw();

}



