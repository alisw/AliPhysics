#include "iostream.h"

void MUONmultest (Int_t evNumber1=0,Int_t evNumber2=0) 
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
   TH1F *hresx = new TH1F("hresx","Residuals",100,-4,4);
   TH1F *hresy = new TH1F("hresy","Residuals",100,-.1,.1);
   TH1F *hresym = new TH1F("hresym","Residuals",100,-500,500);
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

       AliMUONRawCluster  *mRaw;

// Get pointers to Alice detectors and Digits containers
       AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
       TClonesArray *Particles = gAlice->Particles();
       TTree *TR = gAlice->TreeR();
       Int_t nent=TR->GetEntries();
       printf("Found %d entries in the tree (must be one per cathode per event! + 1empty)\n",nent);
       if (MUON) {
	   for (Int_t ich=0;ich<1;ich++) {
	       TClonesArray *MUONrawclust  = MUON->RawClustAddress(ich);
	       TClonesArray *MUONdigits  = MUON->DigitsAddress(ich);
	       for (Int_t icat=1; icat<2; icat++) {
		   MUON->ResetRawClusters();
		   nbytes += TR->GetEvent(icat);
		   Int_t nrawcl = MUONrawclust->GetEntries();
		   printf(
		       "Found %d raw-clusters for cathode %d in chamber %d \n"
		      ,nrawcl,icat,ich+1);


		   gAlice->ResetDigits();

		   Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
		   gAlice->TreeD()->GetEvent(nent-2+icat-1);
		   Int_t ndigits = MUONdigits->GetEntriesFast();
		   printf(
		   "Found %d digits for cathode %d in chamber %d \n"
		       ,ndigits,icat,ich+1);


		   Int_t TotalMult =0;

		   for (Int_t iraw=0; iraw < nrawcl; iraw++) {
		       mRaw = (AliMUONRawCluster*)MUONrawclust->UncheckedAt(iraw);
		       Int_t mult=mRaw->fMultiplicity;
		       h1->Fill(mult,float(1));
		       TotalMult+=mult;
		       Int_t itrack=mRaw->fTracks[0];
		       Float_t xrec=mRaw->fX;
		       Float_t yrec=mRaw->fY;
		       Float_t R=TMath::Sqrt(xrec*xrec+yrec*yrec);
		       if (R > 55.2) continue;
		       if (itrack ==1) continue;
		       Float_t res[2];
		       Int_t nres=0;
		       nbytes=0;
		       gAlice->ResetHits();
		       Int_t nbytes += TH->GetEvent(itrack);
		       
		       for(AliMUONHit* mHit=(AliMUONHit*)MUON->FirstHit(-1); 
			   mHit;
			   mHit=(AliMUONHit*)MUON->NextHit()) 
		       {
			   Int_t   nch   = mHit->fChamber;  // chamber number
			   Float_t x     = mHit->fX;        // x-pos of hit
			   Float_t y     = mHit->fY;        // y-pos
			   if (nch==(ich+1)){
			       hresx->Fill(xrec-x,float(1));
			       hresy->Fill(yrec-y,float(1));
			       if ((yrec-y)*1e4 <500 )
				   hresym->Fill((yrec-y)*1e4,float(1));
			       if (TMath::Abs(yrec-y)>.02) {
				   h22->Fill(mRaw->fX,mRaw->fY,float(1));
				   hmult->Fill(mult,float(1));
			       }
			   } // chamber
		       } //hit
		   } //iraw
		   printf("Total Cluster Multiplicity %d \n" , TotalMult);
	       }  // icat
	   } // ich
       }   // end if MUON
   }   // event loop 
   TCanvas *c1 = new TCanvas("c1","Charge and Residuals",400,10,600,700);
   pad11 = new TPad("pad11"," ",0.01,0.51,0.49,0.99);
   pad12 = new TPad("pad12"," ",0.51,0.51,0.99,0.99);
   pad13 = new TPad("pad13"," ",0.01,0.01,0.49,0.49);
   pad14 = new TPad("pad14"," ",0.51,0.01,0.99,0.49);
   pad11->SetFillColor(11);
   pad12->SetFillColor(11);
   pad13->SetFillColor(11);
   pad14->SetFillColor(11);
   pad11->Draw();
   pad12->Draw();
   pad13->Draw();
   pad14->Draw();
   
   pad11->cd();
   hresx->SetFillColor(42);
   hresx->SetXTitle("xrec-x");
   hresx->Draw();

   pad12->cd();
   hresy->SetFillColor(42);
   hresy->SetXTitle("yrec-y");
   hresy->Draw();

   pad13->cd();
   h1->SetFillColor(42);
   h1->SetXTitle("multiplicity");
   h1->Draw();

   pad14->cd();
   hresym->SetFillColor(42);
   hresym->SetXTitle("yrec-y");
   hresym->Fit("gaus");
   hresym->Draw();
   TCanvas *c2 = new TCanvas("c2","Charge and Residuals",400,10,600,700);
   pad21 = new TPad("pad21"," ",0.01,0.51,0.49,0.99);
   pad22 = new TPad("pad22"," ",0.51,0.51,0.99,0.99);
   pad23 = new TPad("pad23"," ",0.01,0.01,0.49,0.49);
   pad24 = new TPad("pad24"," ",0.51,0.01,0.99,0.49);
   pad21->SetFillColor(11);
   pad22->SetFillColor(11);
   pad23->SetFillColor(11);
   pad24->SetFillColor(11);
   pad21->Draw();
   pad22->Draw();
   pad23->Draw();
   pad24->Draw();
   
   pad21->cd();
   h22->SetFillColor(42);
   h22->SetXTitle("x");
   h22->SetYTitle("y");
   h22->Draw();

   pad22->cd();
   hmult->SetFillColor(42);
   hmult->SetXTitle("multiplicity");
   hmult->Draw();

}



