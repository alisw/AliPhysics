void MUONocc (Int_t evNumber1=0,Int_t evNumber2=0, Int_t ic=1) 
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

    
    if (!file) {
	printf("\n Creating galice.root \n");
	file = new TFile("galice.root");
    } else {
	printf("\n galice.root found in file list");
    }
    file->ls();
    
// Get AliRun object from file or create it if not on file
    if (!gAlice) {
	gAlice = (AliRun*)(file->Get("gAlice"));
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) {
	    printf("\n create new gAlice object");
	    gAlice = new AliRun("gAlice","Alice test program");
	}
    }
    
//  Create some histograms

    TH1F *Hits[10], *HitDensity[10], *Occ0[10],*Occ1[10], *Mult[10];
    for (Int_t i=0; i<10; i++) {
	Hits[i]     =  new TH1F("hits1","Hit Distribution",26,0,260);
	Hits[i]->SetFillColor(0);
	Hits[i]->SetXTitle("R (cm)");

	HitDensity[i]  =  new TH1F("dhits1","Hit Density Distribution",26,0,260);
	HitDensity[i]->SetFillColor(0);
	HitDensity[i]->SetXTitle("R (cm)");

	Occ0[i]  =  new TH1F("occ0","Occupancy Density",26,0,260);
	Occ0[i] -> SetFillColor(0);
	Occ0[i] -> SetXTitle("R (cm)");
	Occ1[i]  =  new TH1F("occ1","Occupancy",26,0,260);
	Occ1[i] -> SetFillColor(0);
	Occ1[i] -> SetXTitle("R (cm)");
	Occ1[i] -> SetYTitle("Occupancy");
	Mult[i]  =  new TH1F("mult","Mult distribution",26,0,260);
	Mult[i] -> SetFillColor(0);
	Mult[i] -> SetXTitle("R (cm)");
    }


   TH1F *theta   =  new TH1F("theta","Theta distribution",180,0,180);
 
   TH1F *emult   =  new TH1F("emult","Event Multiplicity",100,0,1000);   
   AliMUONChamber*  iChamber;
   AliMUONSegmentation* seg;
   
//
//   Loop over events 
//
   Float_t dpx[2], dpy[2];
   Int_t mdig =0;
   Int_t nhit=0;
   
   for (Int_t nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;

       AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
       printf("\n track %d %d \n ", nev, MUON);

       TTree *TH = gAlice->TreeH();
       Int_t ntracks = TH->GetEntries();
//
//   Loop over events
//
       Int_t EvMult=0;
       
       for (Int_t track=0; track<ntracks;track++) {
	   gAlice->ResetHits();
	   Int_t nbytes += TH->GetEvent(track);
	   if (MUON)  {
//
//  Loop over hits
//
	       Int_t NperCh=0;
	       
	       for(AliMUONHit* mHit=(AliMUONHit*)MUON->FirstHit(-1); 
		   mHit;
		   mHit=(AliMUONHit*)MUON->NextHit()) 
	       {
		   Int_t   nch   = mHit->fChamber;  // chamber number
		   Float_t x     = mHit->fX;        // x-pos of hit
		   Float_t y     = mHit->fY;        // y-pos
		   Float_t Eloss = mHit->fEloss;
		   Float_t Theta = mHit->fTheta;
		   theta->Fill(Theta,(float) 1);
		   if (nch >10) continue;
		   if (nch ==1) EvMult++;
//
//
		   iChamber = & MUON->Chamber(nch-1);
		   response= iChamber.GetResponseModel();
		   seg     = iChamber.GetSegmentationModel(1);
		   NperCh++;
		   Int_t i,j;
		   seg->GetPadIxy(x,y,i,j);
		   Int_t isec = seg->Sector(i,j);
		   if (isec ==1 && nch==ic) nhit++;
		   
		   Float_t a=seg->Dpx(isec)*seg->Dpy(isec);
		   Float_t r=TMath::Sqrt(x*x+y*y);
		   Float_t wgt=1/(2*10*TMath::Pi()*r)/(evNumber2+1);
		   wgt=wgt*(1.+24./(2.*TMath::Pi()*r));
		   Hits[nch-1]->Fill(r,(float) 1);
		   HitDensity[nch-1]->Fill(r,wgt);
		   Occ0[nch-1]->Fill(r,wgt*a);
	       } // hit loop
	   } // if MUON
       } // track loop
       emult->Fill(Float_t(EvMult), (Float_t) 1);
       
       Int_t iseg=1;

       for (Int_t ich=0; ich<10; ich++) {
	   iChamber = & MUON->Chamber(ich);
	   seg=iChamber.GetSegmentationModel(iseg);
	   gAlice->ResetDigits();
	   gAlice->TreeD()->GetEvent(iseg);
	   TClonesArray *MUONDigits  = MUON->DigitsAddress(ich);
	   Int_t Ndigits=MUONDigits->GetEntriesFast();
	   AliMUONDigit* dig;
	   printf("\n Reading %d digits\n", Ndigits);
	   if (MUONDigits)  {
	       for (Int_t ndig=0; ndig<Ndigits; ndig++) 
	       {
		   dig = (AliMUONDigit*)MUONDigits->UncheckedAt(ndig); 
		   Int_t i=dig->fPadX;
		   Int_t j=dig->fPadY;
		   Float_t x,y;
		   seg->GetPadCxy(i,j,x,y);
		   Int_t isec = seg->Sector(i,j);
		   Float_t a=seg->Dpx(isec)*seg->Dpy(isec);
		   Float_t r=TMath::Sqrt(x*x+y*y);
//	       if (r<25) 
//	       printf("\n Sector,a %d %f", isec,a);
		   if (isec==1) mdig++;
		   
		   Float_t wgt;
		   wgt=1/(2.*10.*TMath::Pi()*r)/(evNumber2+1)*a;
// Take into account inefficiency due to frames
		   wgt=wgt*(1.+24./(2.*TMath::Pi()*r));
		   Occ1[ich]->Fill(r,wgt);
	       } // digit loop
	       Mult[ich]->Divide(Occ1[ich],Occ0[ich]);
	   } // chamber loop
       } // if MUONDigits
   } // event loop 
//
   printf("\n hits, digits %d %d\n ", nhit, mdig);
   

   
//Create a canvas, set the view range, show histograms
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
   HitDensity[0]->Draw();
   pad12->cd();
   HitDensity[1]->Draw();
   pad13->cd();
   HitDensity[2]->Draw();
   pad14->cd();
   HitDensity[3]->Draw();

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
   HitDensity[4]->Draw();
   pad22->cd();
   HitDensity[5]->Draw();
   pad23->cd();
   HitDensity[6]->Draw();
   pad24->cd();
   HitDensity[7]->Draw();

   TCanvas *c3 = new TCanvas("c3","Hit Densities",400,10,600,700);
   pad31 = new TPad("pad31"," ",0.01,0.51,0.49,0.99);
   pad32 = new TPad("pad32"," ",0.51,0.51,0.99,0.99);
   pad33 = new TPad("pad33"," ",0.01,0.01,0.49,0.49);
   pad34 = new TPad("pad34"," ",0.51,0.01,0.99,0.49);
   pad31->SetFillColor(0);
   pad32->SetFillColor(0);
   pad33->SetFillColor(0);
   pad34->SetFillColor(0);
   pad31->Draw();
   pad32->Draw();
   pad33->Draw();
   pad34->Draw();

   pad31->cd();
   HitDensity[8]->Draw();
   pad32->cd();
   HitDensity[9]->Draw();

   TCanvas *c4 = new TCanvas("c4","Occupancies",400,10,600,700);
   pad41 = new TPad("pad41"," ",0.01,0.51,0.49,0.99);
   pad42 = new TPad("pad42"," ",0.51,0.51,0.99,0.99);
   pad43 = new TPad("pad43"," ",0.01,0.01,0.49,0.49);
   pad44 = new TPad("pad44"," ",0.51,0.01,0.99,0.49);
   pad41->SetFillColor(0);
   pad42->SetFillColor(0);
   pad43->SetFillColor(0);
   pad44->SetFillColor(0);
   pad41->Draw();
   pad42->Draw();
   pad43->Draw();
   pad44->Draw();


   pad41->cd();
   Occ1[0]->Draw();
   pad42->cd();
   Occ1[2]->Draw();
   pad43->cd();
   Occ1[4]->Draw();
   pad44->cd();
   Occ1[6]->Scale(1.25);
   Occ1[6]->Draw();

   TCanvas *c5 = new TCanvas("c5","Occupancies",400,10,600,700);
   pad51 = new TPad("pad41"," ",0.01,0.51,0.49,0.99);
   pad52 = new TPad("pad42"," ",0.51,0.51,0.99,0.99);
   pad53 = new TPad("pad43"," ",0.01,0.01,0.49,0.49);
   pad54 = new TPad("pad44"," ",0.51,0.01,0.99,0.49);
   pad51->SetFillColor(0);
   pad52->SetFillColor(0);
   pad53->SetFillColor(0);
   pad54->SetFillColor(0);
   pad51->Draw();
   pad52->Draw();
   pad53->Draw();
   pad54->Draw();


   pad51->cd();
   Occ1[8]->Scale(1.25);
   
   Occ1[8]->Draw();
}










