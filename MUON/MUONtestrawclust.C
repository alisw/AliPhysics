#include "iostream.h"

void MUONtestrawclust (Int_t evNumber1=0,Int_t evNumber2=0, Int_t ich1=0, Int_t ich2=0) 
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
   TH2F *hchi1 = new TH2F("hchi1","Chi2 vs Residuals",100,0,0.2,100,-500,500);
   TH2F *hchi2 = new TH2F("hchi2","Chi2 vs Residuals",100,0,20,100,-500,500);
//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (int nev=evNumber1; nev<= evNumber2; nev++) {
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
	   for (Int_t ich=ich1;ich<=ich2;ich++) {
	       TClonesArray *MUONrawclust  = MUON->RawClustAddress(ich);
	       //          printf ("MUONrawclust %p \n",MUONrawclust);

	   for (Int_t icat=1; icat<2; icat++) {
	       MUON->ResetRawClusters();
	       nbytes += TR->GetEvent(icat);
	       Int_t nrawcl = MUONrawclust->GetEntries();
	       printf("Found %d raw clusters for cathode %d in chamber %d \n"
		      ,nrawcl,icat,ich+1);
	       for (Int_t iraw=0; iraw < nrawcl; iraw++) {
		   mRaw = (AliMUONRawCluster*)MUONrawclust->UncheckedAt(iraw);
		   Int_t mult=mRaw->fMultiplicity[1];
		   Int_t itrack=mRaw->fTracks[1];
		   printf("\n mult1 mult2 %d %d chi2 %f itrack %d"
			  ,mRaw->fMultiplicity[0], mRaw->fMultiplicity[1], mRaw->fChi2[0], itrack);
		   h1->Fill(mult,float(1));

		   Float_t xrec=mRaw->fX[1];
		   Float_t yrec=mRaw->fY[0];
		   Float_t R=TMath::Sqrt(xrec*xrec+yrec*yrec);
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
			       hchi1->Fill(mRaw->fChi2[0],(yrec-y)*1e4,float(1));
			       hchi2->Fill(mRaw->fChi2[0],(yrec-y)*1e4,float(1));
			       
			       if ((yrec-y)*1e4 <500 )
				   hresym->Fill((yrec-y)*1e4,float(1));
			       if (mRaw->fChi2[0]>.3) {
				   h22->Fill(mRaw->fX[1],mRaw->fY[0],float(1));
				   hmult->Fill(mult,float(1));
			       }
			   } // chamber
		       } //hit
	       } //iraw
	   }  // icat
       } // ich
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
   h1->SetFillColor(42);
   h1->SetXTitle("multiplicity");
   h1->Draw();

   c1->cd(4);
   hresym->SetFillColor(42);
   hresym->SetXTitle("yrec-y");
   hresym->Fit("gaus");
   hresym->Draw();

   TCanvas *c2 = new TCanvas("c2","Charge and Residuals",400,10,600,700);
   c2->Divide(2,2);

   c2->cd(1);
   h22->SetFillColor(42);
   h22->SetXTitle("x");
   h22->SetYTitle("y");
   h22->Draw();

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



