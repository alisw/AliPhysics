void RICHchamberView (Int_t evNumber=0, Int_t ChamberView = 3) 
{

/////////////////////////////////////////////////////////////////////////

    gClassTable->GetID("AliRun");


// Dynamically link some shared libs
 
    if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
    }else {
      delete gAlice;
      gAlice = 0;
   }

    gAlice=0;
    
// Connect the Root Galice file containing Geometry, Kine and Hits
    
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (!file) file = new TFile("galice.root","UPDATE");
    
// Get AliRun object from file or create it if not on file
    
    if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file %d \n",gAlice);
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
    else {
      delete gAlice;
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }

   TH2F *ChView = new TH2F("ChView","RICH DISPLAY",160,0,160,144,0,144);

   gStyle->SetPalette(1);

   TCanvas *view = new TCanvas("Display","ALICE RICH Display",0,0,1200,750);

   ChView->SetStats(0);

   ChView->SetMaximum(100);

   Int_t Nevents = gAlice->GetEventsPerRun();

   printf(" Events in Run: %d\n",Nevents);

   Int_t ich = ChamberView - 1;

   AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
   
   AliRICHSegmentationV0*  segmentation;
   AliRICHChamber*       chamber;
   
   chamber = &(RICH->Chamber(ich));
   segmentation=(AliRICHSegmentationV0*) chamber->GetSegmentationModel(ich);
   
   Int_t NpadX = segmentation->Npx();                 // number of pads on X
   Int_t NpadY = segmentation->Npy();                 // number of pads on Y
   
   printf(" NpadX: %d NpadY: %d \n",NpadX,NpadY);
   
   //   Start loop over events 

   for (int nev=evNumber; nev<= Nevents; nev++) {

       Int_t nparticles = gAlice->GetEvent(nev);

       if (nparticles == -1) break;

       printf("Particles:%d\n",nparticles);

       gAlice->ResetDigits();

       gAlice->TreeD()->GetEvent(0);

       printf("gAlice D: %x\n",gAlice);

       TClonesArray *Digits = RICH->DigitsAddress(ich);    //  Raw clusters branch

       Int_t ndigits = Digits->GetEntriesFast();

       printf("Digits:%d %d\n",Digits,ndigits);

       for (Int_t hit=0;hit<ndigits;hit++) {

	 AliRICHDigit *dHit = (AliRICHDigit*) Digits->UncheckedAt(hit);

	 //	 printf(" dHit %d\n",dHit);

	 Int_t qtot = dHit->Signal();                       // charge
	 Int_t ipx  = dHit->PadX() + NpadX/2;               // pad number on X
	 Int_t ipy  = dHit->PadY() + NpadY/2;               // pad number on Y

	 ChView -> Fill((Float_t)ipx,(Float_t)ipy,(Float_t)qtot); 

	 //	 printf(" X: %d Y: %d Q: %d\n",ipx,ipy,qtot);

       }

       
              
       gAlice->TreeR()->GetEvent(0);
       Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
       printf("ich: %d nent: %d\n",ich,nent);
       gAlice->TreeR()->GetEvent(nent-1);

       TClonesArray *RecRings = RICH->RecHitsAddress1D(ich);

       printf(" Ring Pointer: %x\n",RecRings);

       Int_t nRecRings = RecRings->GetEntriesFast();
 
       printf(" nRings: %d\n",nRecRings);
       

       for (Int_t ring=0;ring<nRecRings;ring++) {

	 AliRICHRecHit1D *recHit1D = (AliRICHRecHit1D*) RecRings->UncheckedAt(ring);

	 printf(" Pointer to PatRec: %d \n",recHit1D);

	 Float_t r_omega = recHit1D->fOmega;                  // Cerenkov angle
	 Float_t *cer_pho = recHit1D->fCerPerPhoton;        // Cerenkov angle per photon
	 Int_t *padsx = recHit1D->fPadsUsedX;           // Pads Used fo reconstruction (x)
	 Int_t *padsy = recHit1D->fPadsUsedY;           // Pads Used fo reconstruction (y)
	 Int_t goodPhotons = recHit1D->fGoodPhotons;    // Number of pads used for reconstruct

         printf(" Theta Cerenkov %d : %f NgoodPhotons: %d\n",ring,r_omega,goodPhotons);
       }

       
       ChView->Draw("colz");
       
       view->Modified();
       view->Update();
       
       gSystem->Sleep(1000);
       
       cout << endl;
       
       ChView->Reset();
       
   }
   
   printf("\nEnd of macro\n");
   printf("**********************************\n");
   
}


