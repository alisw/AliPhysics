#include "iostream.h"

void MUONdoubles (Int_t evNumber1=0,Int_t evNumber2=0) 
{
    {
    hprof  = new TProfile("hprof","Profile dmin vs r",24,0,120,0,50);
    }
    hdist1  = new TH1F("hdist1","distance",100,0,10);
    hdist2  = new TH1F("hdist2","distance",100,0,10);
    Float_t a1=3.7/TMath::Sqrt(4);
    Float_t a2=26/TMath::Sqrt(4);    
    for (Int_t i=1; i<10000; i++) {
	Float_t x1,x2,y1,y2;
	x1 = (2*gRandom->Rndm(i)-1.)*a1;
	x2 = (2*gRandom->Rndm(i)-1.)*a1;
	y1 = (2*gRandom->Rndm(i)-1.)*a1;
	y2 = (2*gRandom->Rndm(i)-1.)*a1;	
	Float_t d=TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	hdist1->Fill(d,1);
    }
    Float_t xh[100], yh[100];
    for (Int_t k=0; k<1000; k++) {
	for (Int_t i=0; i<100; i++) {
	    xh[i] = (2*gRandom->Rndm(i)-1.)*a2;
	    yh[i] = (2*gRandom->Rndm(i)-1.)*a2;
	}
	Float_t x1,y1;
	
	for (Int_t i=0; i<10; i++) {
	    Float_t dmin=1000;
	    x1 = (2*gRandom->Rndm(i)-1.)*a2;
	    y1 = (2*gRandom->Rndm(i)-1.)*a2;
	    for (Int_t j=0; j<100; j++) {
		Float_t d=TMath::Sqrt((x1-xh[j])*(x1-xh[j])+(y1-yh[j])*(y1-yh[j]));
		if (d<dmin) dmin=d;
	    }
	    hdist2->Fill(dmin,1);
	}
    }
    
    

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
    if (!file) file = new TFile("galice.root","UPDATE");

// Get AliRun object from file or create it if not on file

    if (!gAlice) {
	gAlice = (AliRun*)file->Get("gAlice");
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
//
// Set reconstruction models
//
// Get pointers to Alice detectors and Digits containers
    AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
//
// Book some histograms
    TH1F *Hcentre[10], *Hall[10], *Hcfrac[10], *Htail[10], *Htfrac[10]; 
    TH1F *Hdcut[10], *Hdall[10], *Hdfrac[10];

    for (Int_t i=0; i<10; i++) {
	Hcentre[i]     =  new TH1F("Hcentre","Hit Distribution",26,0,260);
	Hcentre[i]->SetFillColor(0);
	Hcentre[i]->SetXTitle("R (cm)");

	Hall[i]     =  new TH1F("Hall","Hit Distribution",26,0,260);
	Hall[i]->SetFillColor(0);
	Hall[i]->SetXTitle("R (cm)");

	Hcfrac[i]     =  new TH1F("Hcfrag","Hit Distribution",26,0,260);
	Hcfrac[i]->SetFillColor(0);
	Hcfrac[i]->SetXTitle("R (cm)");

	Htail[i]     =  new TH1F("Htail","Hit Distribution",26,0,260);
	Htail[i]->SetFillColor(0);
	Htail[i]->SetXTitle("R (cm)");

	Htfrac[i]     =  new TH1F("Htfrag","Hit Distribution",26,0,260);
	Htfrac[i]->SetFillColor(0);
	Htfrac[i]->SetXTitle("R (cm)");

	Hdall[i]     =  new TH1F("Hdall","Hit Distribution",26,0,260);
	Hdall[i]->SetFillColor(0);
	Hdall[i]->SetXTitle("R (cm)");

	Hdcut[i]     =  new TH1F("Hdcut","Hit Distribution",26,0,260);
	Hdcut[i]->SetFillColor(0);
	Hdcut[i]->SetXTitle("R (cm)");

	Hdfrac[i]     =  new TH1F("Hdfrac","Hit Distribution",26,0,260);
	Hdfrac[i]->SetFillColor(0);
	Hdfrac[i]->SetXTitle("R (cm)");

    }

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


	TClonesArray *Particles = gAlice->Particles();
	TTree *TD = gAlice->TreeD();
	Int_t nent=TD->GetEntries();
	printf("Found %d entries in the tree (must be one per cathode per event!)\n",nent);
	Float_t x[2000];
	Float_t y[2000];
	Int_t nhit=0;
	
	if (MUON) {
//
// Loop on chambers and on cathode planes
//
	    TTree *TH = gAlice->TreeH();
	    Int_t ntracks = TH->GetEntries();

//
//              Loop over Hits
//
	    for (i=0; i<1000; i++) {
		Float_t r  =100.*gRandom->Rndm(i);
		Float_t phi=2.*TMath::Pi()*gRandom->Rndm(i);
		Float_t xr=r*TMath::Sin(phi);
		Float_t yr=r*TMath::Cos(phi);		
		Hdall[0]->Fill(r,1.);
		Float_t dmin=100000.;
		
		if (i==0) {
		    for (Int_t track=0; track<ntracks; track++) {
			gAlice->ResetHits();
			Int_t nbytes += TH->GetEvent(track);
			for(AliMUONHit* mHit=(AliMUONHit*)MUON->FirstHit(-1); 
			    mHit;
			    mHit=(AliMUONHit*)MUON->NextHit()) 
			{
			    Int_t   nch   = mHit->fChamber;  // chamber number
			    if (nch!=1) continue;			
			    
			    x[nhit]     = mHit->fX;        // x-pos of hit
			    y[nhit]     = mHit->fY;        // y-pos
			    nhit++;
			} // hit loop 
		    } // track loop
		} else {
		    for (Int_t nh=0; nh<nhit; nh++) {
			Float_t d=TMath::Sqrt((x[nh]-xr)*(x[nh]-xr)+(y[nh]-yr)*(y[nh]-yr));
//			printf ("\n r,d %f %f",r,d);
			Float_t dx=TMath::Abs(x[nh]-xr);
			Float_t dy=TMath::Abs(y[nh]-yr);			
			if (d < dmin) dmin=d;
			
			if (dx<0.5 && dy < 0.75) Hdcut[0]->Fill(r,1.);			    
		    } // hit loop
		} // i==0
		if (r>20) hprof->Fill(r,dmin,1);
	    } // random loop
	    Int_t icat=0;
	    gAlice->ResetDigits();
	    gAlice->TreeD()->GetEvent(icat+1); // spurious +1 ...

	    for (Int_t ich=0;ich<1;ich++) {
		AliMUONChamber iChamber= MUON->Chamber(ich);
		TClonesArray *MUONdigits  = MUON->DigitsAddress(ich);
		if (MUONdigits == 0) continue;
		//
		// Get ready the current chamber stuff
		//
		AliMUONResponse* response = iChamber.GetResponseModel();
		AliMUONSegmentation*  seg = iChamber.GetSegmentationModel(icat+1);
		HitMap = new AliMUONHitMapA1(seg, MUONdigits);
		HitMap->FillHits();
		Int_t nxmax=seg->Npx();
		Int_t nymax=seg->Npy();
//
// generate random positions on the chamber
//		
		for (Int_t i=0; i<4000; i++) {
		    Int_t ix = 2*Int_t((gRandom->Rndm(i)-0.5)*nxmax);
		    Int_t iy = 2*Int_t((gRandom->Rndm(i)-0.5)*nymax);
		    
// test for centre overlap
//
		    Float_t xp,yp;
		    seg->GetPadCxy(ix,iy,xp,yp);
		    Float_t r=TMath::Sqrt(xp*xp+yp*yp);
		    if (TMath::Abs(xp)>3 && TMath::Abs(yp)>3) 
			Hall[ich]->Fill(r,1.);
		    if (HitMap->TestHit(ix,iy) != 0) {
			Hcentre[ich]->Fill(r,1.);
		    }
// test for neighbour overlap
//		   
		    Int_t nn;
		    Int_t X[20], Y[20];
		    seg->Neighbours(ix,iy,&nn,X,Y);
		    Bool_t hit=false;
		    
		    for (Int_t j=0; j<nn; j++) {
			if (HitMap->TestHit(X[j],Y[j]) != 0) hit=true;
		    }
		    if (hit &&  HitMap->TestHit(ix,iy) == 0) Htail[ich]->Fill(r,1.);
		} //random loop
		delete HitMap;
	    } // chamber loop
	}   // event loop
    } // if MUON
    
    for (Int_t ich=0; ich<10; ich++) {
	Hcfrac[ich]->Divide(Hcentre[ich],Hall[ich]);
	Htfrac[ich]->Divide(Htail[ich],Hall[ich]);
	Hdfrac[ich]->Divide(Hdcut[ich],Hdall[ich]);
    }

    TCanvas *c1 = new TCanvas("c1","Hit Densities",400,10,600,700);
    pad11 = new TPad("pad11"," ",0.01,0.51,0.49,0.99);
    pad12 = new TPad("pad12"," ",0.51,0.51,0.99,0.99);
    pad13 = new TPad("pad13"," ",0.01,0.01,0.49,0.49);
    pad14 = new TPad("pad14"," ",0.51,0.01,0.99,0.49);
    pad11->SetFillColor(0);
    pad12->SetFillColor(0);
    pad13->SetFillColor(0);
    pad14->SetFillColor(0);
    pad11->Draw();
    pad12->Draw();
    pad13->Draw();
    pad14->Draw();
   
    pad11->cd();
    Hcentre[0]->Draw();
    
    pad12->cd();
    Hall[0]->Draw();
    
   pad13->cd();
   Hcfrac[0]->Draw();

   pad14->cd();
   Htfrac[0]->Draw();
   
   TCanvas *c2 = new TCanvas("c2","Hit Densities",400,10,600,700);
   pad21 = new TPad("pad21"," ",0.01,0.51,0.49,0.99);
   pad22 = new TPad("pad22"," ",0.51,0.51,0.99,0.99);
   pad23 = new TPad("pad23"," ",0.01,0.01,0.49,0.49);
   pad24 = new TPad("pad24"," ",0.51,0.01,0.99,0.49);
   pad21->SetFillColor(0);
   pad22->SetFillColor(0);
   pad23->SetFillColor(0);
   pad24->SetFillColor(0);
   pad21->Draw();
   pad22->Draw();
   pad23->Draw();
   pad24->Draw();
   
   pad21->cd();
   hdist1->Draw();

   pad22->cd();
   hdist2->Draw();

   pad23->cd();
   Hdfrac[0]->Draw();

   pad24->cd();
   hprof->Draw();
   file->Close();
}

