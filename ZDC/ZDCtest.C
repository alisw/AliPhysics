void ZDCtest (Int_t detector=0, Int_t evTot = 0) 
{
   delete gAlice;
   gAlice=0;
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }
      
// Connect the Root Galice file containing Geometry, Kine, Hits and Digits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("Hijing_b2.root");
    if (!file) {
	printf("\n Creating galice.root \n");
	file = new TFile("Hijing_b2.root");
    } else {
	printf("\n galice.root found in file list");
    }
    file->ls();

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) {
	    printf("\n create new gAlice object");
	    gAlice = new AliRun("gAlice","Alice test program");
	}
   }
   file->ls();
   
// Create some histograms

   TH2F *hspotzn  = new TH2F("hspotzn","Y vs X on ZN front face",100,-4.,4.,100,-4.,4.);
   hspotzn -> SetXTitle("X on ZN");
   hspotzn -> SetYTitle("Y on ZN");
   TH2F *hspotzp  = new TH2F("hspotzp","Y vs X on ZP front face",100,-12.,12.,100,-12.,12.);
   hspotzp -> SetXTitle("X on ZP");
   hspotzp -> SetYTitle("Y on ZP");
   TH2F *hspotzem = new TH2F("hspotzem","Y vs X on ZEM front face",100,-4.,4.,100,-4.,4.);
   hspotzem -> SetXTitle("X on ZEM");
   hspotzem -> SetYTitle("Y on ZEM");

   TH1F *hEzn    = new TH1F("hEzn","Energy deposited in ZN",100,0,4.e3);
   hEzn -> SetXTitle("E (GeV)");
   TH1F *hLzn    = new TH1F("hLzn", "Total light in ZN    ",100,0,4.e3);
   hLzn -> SetXTitle("phe");
   TH1F *hPMCzn  = new TH1F("hPMCzn", "Light in common PM    ",100,0,1.e3);
   hPMCzn -> SetXTitle("phe");
   TH1F *hPMQ1zn = new TH1F("hPMQ1zn","Light in quadrant 1 PM",100,0,1.e3);
   hPMQ1zn -> SetXTitle("phe");
   TH1F *hPMQ2zn = new TH1F("hPMQ2zn","Light in quadrant 2 PM",100,0,1.e3);
   hPMQ2zn -> SetXTitle("phe");
   TH1F *hPMQ3zn = new TH1F("hPMQ3zn","Light in quadrant 3 PM",100,0,1.e3);
   hPMQ3zn -> SetXTitle("phe");
   TH1F *hPMQ4zn = new TH1F("hPMQ4zn","Light in quadrant 4 PM",100,0,1.e3);
   hPMQ4zn -> SetXTitle("phe");

   TH1F *hEzp    = new TH1F("hEzp","Energy deposited in ZP",100,0,4.e3);
   hEzp -> SetXTitle("E (GeV)");
   TH1F *hLzp    = new TH1F("hLzp", "Total light in ZP    ",100,0,4.e3);
   hLzp -> SetXTitle("phe");
   TH1F *hPMCzp  = new TH1F("hPMCzp","Light in common PM     ",100,0,1.e3);
   hPMCzp -> SetXTitle("phe");
   TH1F *hPMQ1zp = new TH1F("hPMQ1zp","Light in quadrant 1 PM",100,0,1.e3);
   hPMQ1zp -> SetXTitle("phe");
   TH1F *hPMQ2zp = new TH1F("hPMQ2zp","Light in quadrant 2 PM",100,0,1.e3);
   hPMQ2zp -> SetXTitle("phe");
   TH1F *hPMQ3zp = new TH1F("hPMQ3zp","Light in quadrant 3 PM",100,0,1.e3);
   hPMQ3zp -> SetXTitle("phe");
   TH1F *hPMQ4zp = new TH1F("hPMQ4zp","Light in quadrant 4 PM",100,0,1.e3);
   hPMQ4zp -> SetXTitle("phe");


   TH1F *hEzem   = new TH1F("hEzem","Energy deposited in ZEM",100,0,1.e3);
   hEzem -> SetXTitle("E (GeV)");
   TH1F *hPMzem  = new TH1F("hPMzem","Light produced in ZEM PM",100,0,4.e3);
   hPMzem -> SetXTitle("phe");
   
   //
   
   TH1F *dPMCzn  = new TH1F("dPMCzn","Common PM     ",100,0,1.e3);
   dPMCzn -> SetXTitle("ADC channels");
   TH1F *dPMQ1zn = new TH1F("dPMQ1zn","Quadrant 1 PM",100,0,1.e3);
   dPMQ1zn -> SetXTitle("ADC channels");
   TH1F *dPMQ2zn = new TH1F("dPMQ2zn","Quadrant 2 PM",100,0,1.e3);
   dPMQ2zn -> SetXTitle("ADC channels");
   TH1F *dPMQ3zn = new TH1F("dPMQ3zn","Quadrant 3 PM",100,0,1.e3);
   dPMQ3zn -> SetXTitle("ADC channels");
   TH1F *dPMQ4zn = new TH1F("dPMQ4zn","Quadrant 4 PM",100,0,1.e3);
   dPMQ4zn -> SetXTitle("ADC channels");
   TH1F *dZN     = new TH1F("dZN","Total light in ZN",100,0,1.e3);
   dZN -> SetXTitle("ADC channels");
   
   TH1F *dPMCzp  = new TH1F("dPMCzp","Common PM     ",100,0,1.e3);
   dPMCzp -> SetXTitle("ADC channels");
   TH1F *dPMQ1zp = new TH1F("dPMQ1zp","Quadrant 1 PM",100,0,1.e3);
   dPMQ1zp -> SetXTitle("ADC channels");
   TH1F *dPMQ2zp = new TH1F("dPMQ2zp","Quadrant 2 PM",100,0,1.e3);
   dPMQ2zp -> SetXTitle("ADC channels");
   TH1F *dPMQ3zp = new TH1F("dPMQ3zp","Quadrant 3 PM",100,0,1.e3);
   dPMQ3zp -> SetXTitle("ADC channels");
   TH1F *dPMQ4zp = new TH1F("dPMQ4zp","Quadrant 4 PM",100,0,1.e3);
   dPMQ4zp -> SetXTitle("ADC channels");
   TH1F *dZP     = new TH1F("dZP","Total light in ZP",100,0,1.e3);
   dZP -> SetXTitle("ADC channels");

   TH1F *dZEM     = new TH1F("dZEM","Total light in ZEM",100,0,1.e3);
   dZEM -> SetXTitle("ADC channels");

//
//   Loop over events 
//

   for (Int_t evNumber=0; evNumber<evTot; evNumber++){

// Import the Kine and Hits Trees for the event evNumber in the file
   Int_t nparticles = gAlice->GetEvent(evNumber);
   if (nparticles <= 0) return;
   printf("\n	--- nparticles = %d\n",nparticles);
   
   Float_t energy, EtotZN=0, EtotZP=0, LightCzn=0, LightCzp=0, LtotZN=0, LtotZP=0;
   Int_t nbytes=0, nbytesd=0, ipart, nhits, ndigits, pdgcode, ADCzn, ADCzp;
   TParticle   *particle;
   AliZDCHit   *ZDChit;
   AliZDCDigit *ZDCdigit;

// Get pointers to Alice detectors and Hits containers
   AliDetector *ZDC  = gAlice->GetDetector("ZDC");
   TClonesArray *Particles = gAlice->Particles();
   if (ZDC) {
     TClonesArray *ZDChits    = ZDC->Hits();
     TClonesArray *ZDCdigits  = ZDC->Digits();
   }

// # of entries in Hits tree
   TTree *TH = gAlice->TreeH();
   Int_t ntracks = TH->GetEntries();

// # of entries in Digits tree
   TTree *TD = gAlice->TreeD();
   Int_t ndigen = TD->GetEntries();

   gAlice->ResetDigits();
   nbytesd += TD->GetEvent(ndigen-1);

   if (ZDC) {
     ndigits = ZDCdigits->GetEntries();
     printf("\n	Digits Tree --- # of entries: %d; # of digits: %d\n",ndigen, ndigits);
   }
   for(Int_t digit=0; digit<ndigits; digit++) {
     ZDCdigit = (AliZDCDigit*)ZDCdigits->UncheckedAt(digit);
     printf("\n	Digit# %d, fDetector = %d, fVolume = %d, fADCValue = %f\n",
            digit,ZDCdigit->fDetector,ZDCdigit->fQuadrant,ZDCdigit->fADCValue);
   }
   
// Start loop on tracks in the hits containers
      for (Int_t track=0; track<ntracks; track++) {
         gAlice->ResetHits();
         nbytes += TH->GetEvent(track);

         if (ZDC) {
//          nhits = ZDChits->GetEntries();
//          nhits = ZDChits->GetLast()+1;
           nhits   = ZDChits->GetEntriesFast();
//           printf("\n	Hits Tree --- Event %d track %d nhits %d\n",evNumber,track,nhits);

	   particle = (TParticle*)Particles->UncheckedAt(track);
//	   pdgcode = particle->GetPdgCode();
//	   printf("\nParticle %d\n",pdgcode);

           for(Int_t hit=0; hit<nhits; hit++) {
	     ZDChit   = (AliZDCHit*)ZDChits->UncheckedAt(hit);
	    
	     // Print of the hits
//	     printf("\nHit # %d, fVolume = %d %d\n",hit,ZDChit->fVolume[0],ZDChit->fVolume[1]);
//	     printf("Primary energy = %f, Secondary Flag = %d\n",ZDChit->fPrimKinEn,ZDChit->fSFlag);
//	     printf("Impact point -> %f %f\n",ZDChit->fXImpact,ZDChit->fYImpact);
//	     printf("Energy = %f, Light in quadrant = %f, Light in common PM = %f\n\n",
//	             ZDChit->fEnergy,ZDChit->fLightPMQ,ZDChit->fLightPMC);
	    
	     // Filling histos
	     if(ZDChit->fVolume[0]==1) { //ZN
	       if(ZDChit->fVolume[1]==1)hPMQ1zn->Fill(ZDChit->fLightPMQ);
	       if(ZDChit->fVolume[1]==2)hPMQ2zn->Fill(ZDChit->fLightPMQ);
	       if(ZDChit->fVolume[1]==3)hPMQ3zn->Fill(ZDChit->fLightPMQ);
	       if(ZDChit->fVolume[1]==4)hPMQ4zn->Fill(ZDChit->fLightPMQ);
	       EtotZN   += ZDChit->fEnergy;
	       LtotZN   += (ZDChit->fLightPMQ) + (ZDChit->fLightPMC);
	       LightCzn += ZDChit->fLightPMC;
	       hspotzn->Fill(ZDChit->fXImpact,ZDChit->fYImpact);
	     }
	     if(ZDChit->fVolume[0]==2) { //ZP
	       if(ZDChit->fVolume[1]==1)hPMQ1zp->Fill(ZDChit->fLightPMQ);
	       if(ZDChit->fVolume[1]==2)hPMQ2zp->Fill(ZDChit->fLightPMQ);
	       if(ZDChit->fVolume[1]==3)hPMQ3zp->Fill(ZDChit->fLightPMQ);
	       if(ZDChit->fVolume[1]==4)hPMQ4zp->Fill(ZDChit->fLightPMQ);
	       EtotZP   += ZDChit->fEnergy;
	       LtotZP   += (ZDChit->fLightPMQ) + (ZDChit->fLightPMC);
	       LightCzp += ZDChit->fLightPMC;
	       hspotzp->Fill(ZDChit->fXImpact,ZDChit->fYImpact);
	     }
	     if(ZDChit->fVolume[0]==3) { //ZEM
	       hEzem->Fill(ZDChit->fEnergy);
	       hPMzem->Fill(ZDChit->fLightPMQ);
	       hspotzem->Fill(ZDChit->fXImpact,ZDChit->fYImpact);
	     }
	  	 
           }
           if(nhits!=0){
             hEzn->Fill(EtotZN);
             hLzn->Fill(LtotZN);
             hPMCzn->Fill(LightCzn);
             hEzp->Fill(EtotZP);
             hLzp->Fill(LtotZP);
             hPMCzp->Fill(LightCzp);
             printf("\n	Histos var -> Ezn = %f, Lzn = %f, Ezp = %f, Lzp = %f \n\n",
                     EtotZN, LtotZN, EtotZP, LtotZP);
           }
         }//ZDC        
      }//Track loop
   }//Hit loop 
   
// Control prints


if(detector == 1){ // ZN histos      
   TCanvas *c1 = new TCanvas("c1","ZN hits",0,10,580,700);
   c1->cd();
   TPad *pad1= new TPad("pad1"," ",0.01,0.51,0.99,0.99);
   TPad *pad2= new TPad("pad2"," ",0.01,0.01,0.49,0.49);
   TPad *pad3= new TPad("pad3"," ",0.51,0.01,0.99,0.49);
   pad1->SetFillColor(18);
   pad2->SetFillColor(18);
   pad3->SetFillColor(18);
   pad1->Draw();
   pad2->Draw();
   pad3->Draw();
   pad1->cd();
   hspotzn->Draw();
   pad2->cd();
   hEzn->Draw();
   pad3->cd();
   hPMCzn->Draw();
   
   TCanvas *c2 = new TCanvas("c2","ZN hits",600,10,600,700);
   c2->cd();
   TPad *pad4= new TPad("pad4"," ",0.01,0.51,0.49,0.99);
   TPad *pad5= new TPad("pad5"," ",0.51,0.51,0.99,0.99);
   TPad *pad6= new TPad("pad6"," ",0.01,0.01,0.49,0.49);
   TPad *pad7= new TPad("pad7"," ",0.51,0.01,0.99,0.49);
   pad4->SetFillColor(18);
   pad5->SetFillColor(18);
   pad6->SetFillColor(18);
   pad7->SetFillColor(18);
   pad4->Draw();
   pad5->Draw();
   pad6->Draw();
   pad7->Draw();
   pad4->cd();
   hPMQ1zn->Draw();
   pad5->cd();
   hPMQ2zn->Draw();
   pad6->cd();
   hPMQ3zn->Draw();
   pad7->cd();
   hPMQ4zn->Draw();

   TCanvas *c3 = new TCanvas("c3","ZN hits",300,10,600,700);
   c3->cd();
   hLzn->Draw();
}

if(detector == 2){ // ZP histos      
   TCanvas *c1 = new TCanvas("c1","ZP hits",0,10,580,700);
   c1->cd();
   TPad *pad1= new TPad("pad1"," ",0.01,0.51,0.99,0.99);
   TPad *pad2= new TPad("pad2"," ",0.01,0.01,0.49,0.49);
   TPad *pad3= new TPad("pad3"," ",0.51,0.01,0.99,0.49);
   pad1->SetFillColor(18);
   pad2->SetFillColor(18);
   pad3->SetFillColor(18);
   pad1->Draw();
   pad2->Draw();
   pad3->Draw();
   pad1->cd();
   hspotzp->Draw();
   pad2->cd();
   hEzp->Draw();
   pad3->cd();
   hPMCzp->Draw();
   
   TCanvas *c2 = new TCanvas("c2","ZP hits",600,10,600,700);
   c2->cd();
   TPad *pad4= new TPad("pad4"," ",0.01,0.51,0.49,0.99);
   TPad *pad5= new TPad("pad5"," ",0.51,0.51,0.99,0.99);
   TPad *pad6= new TPad("pad6"," ",0.01,0.01,0.49,0.49);
   TPad *pad7= new TPad("pad7"," ",0.51,0.01,0.99,0.49);
   pad4->SetFillColor(18);
   pad5->SetFillColor(18);
   pad6->SetFillColor(18);
   pad7->SetFillColor(18);
   pad4->Draw();
   pad5->Draw();
   pad6->Draw();
   pad7->Draw();
   pad4->cd();
   hPMQ1zp->Draw();
   pad5->cd();
   hPMQ2zp->Draw();
   pad6->cd();
   hPMQ3zp->Draw();
   pad7->cd();
   hPMQ4zp->Draw();

   TCanvas *c3 = new TCanvas("c3","ZP hits",300,10,600,700);
   c3->cd();
   hLzp->Draw();
}

if(detector == 3){ // ZEM histos      
   TCanvas *c1 = new TCanvas("c1","ZEM hits",0,10,580,700);
   c1->cd();
   TPad *pad1= new TPad("pad1"," ",0.01,0.51,0.99,0.99);
   TPad *pad2= new TPad("pad2"," ",0.01,0.01,0.99,0.49);
   pad1->SetFillColor(18);
   pad2->SetFillColor(18);
   pad1->Draw();
   pad2->Draw();
   pad1->cd();
   hspotzem->Draw();
   pad2->cd();
   hEzem->Draw();
   
   TCanvas *c2 = new TCanvas("c2","ZEM hits",600,10,600,700);
   c2->cd();
   hPMzem->Draw();
}
//   file->Close();   
}
