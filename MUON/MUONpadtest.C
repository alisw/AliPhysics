void MUONpadtest (Int_t evNumber1=0,Int_t evNumber2=1) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

  const Int_t NpadX = 252;                 // number of pads on X
  const Int_t NpadY = 374;                 // number of pads on Y

// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }

//  Create some histograms

   Int_t xmin=-NpadX/2;  
   Int_t xmax= NpadX/2;
   Int_t ymin=-NpadY/2;
   Int_t ymax= NpadY/2;

   TH2F *hc = new TH2F("hc","Chamber1 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *h = new TH2F("h","Chamber1 hit distribution",100,-100,100,100,-100,100);
   TH1F *charge = new TH1F("charge","Charge distribution",100,0.,1000.);

//   Start loop over events 

   Int_t Nh=0;
   Int_t Nh1=0;
   AliMUON *pMUON  = (AliDetector*) gAlice->GetDetector("MUON");
   for (int nev=0; nev<= evNumber2; nev++) {
      Int_t nparticles = gAlice->GetEvent(nev);
      cout<<"nev  "<<nev<<endl;
      cout<<"nparticles  "<<nparticles<<endl;
      if (nev < evNumber1) continue;
      if (nparticles <= 0) return;
     
// Get pointers to MUON detector and Hits containers
      TTree *TH = gAlice->TreeH();
      Int_t ntracks = TH->GetEntries();

// Start loop on tracks in the hits containers

      //      Int_t Nh=0;
      Int_t Nc=0;
      for (Int_t track=0; track<ntracks;track++) {
//	  printf("ntracks %d\n",ntracks);
          gAlice->ResetHits();
          Int_t nbytes += TH->GetEvent(track);
          if (pMUON)  {
              TClonesArray *Clusters = pMUON->PadHits();  // Cluster branch
              TClonesArray *Hits = pMUON->Hits();         // Hits branch
//	      printf("%d %d \n",Clusters,Hits);
          }
          //see hits distribution
          Int_t nhits = Hits->GetEntriesFast();
          if (nhits) Nh+=nhits;
//	  printf("nhits %d\n",nhits);
          for (Int_t hit=0;hit<nhits;hit++) {
//	      printf("hit# %d\n",hit);
              mHit = (AliMUONHit*) Hits->UncheckedAt(hit);
              Int_t nch  = mHit->Chamber();              // chamber number
              Float_t x  = mHit->X();                    // x-pos of hit
              Float_t y  = mHit->Y();                    // y-pos
              // Fill the histograms
	      if( nch==1) {
		  Float_t rhit=TMath::Sqrt(x*x+y*y);
		  if( rhit<= 55 ) Nh1+=nhits;
                  h->Fill(x,y,(float) 1);
              }


          //    see signal distribution
	      for (AliMUONPadHit* mPad =
		       (AliMUONPadHit*)pMUON->FirstPad(mHit,Clusters);
		   mPad;
		   mPad = (AliMUONPadHit*)pMUON->NextPad(Clusters))
	      {
		  Int_t nhit   = mPad->HitNumber();          // hit number
		  Int_t qtot   = mPad->Q();                  // charge
		  Int_t ipx    = mPad->PadX();               // pad number on X
		  Int_t ipy    = mPad->PadY();               // pad number on Y
		  Int_t iqpad  = mPad->QPad();               // charge per pad
//		  printf("%d %d %d\n",ipx,ipy,iqpad);
		  if (qtot > 0 && nch == 1){
		      charge->Fill(float(iqpad),(float) 1);
		  }
		  if(nch == 1) {
		      hc->Fill(ipx,ipy,(float) iqpad);
		  } // if rpad	 	
	      } // padhits
          } // hitloops
      }
      //      cout<<"Nh  "<<Nh<<endl;
   }
//Create a canvas, set the view range, show histograms
   TCanvas *c1 = new TCanvas("c1","Alice MUON pad hits",400,10,600,700);
   hc->SetXTitle("ix (npads)");
   hc->Draw();
   TCanvas *c2 = new TCanvas("c2","Alice MUON hits",400,10,600,700);
   h->SetFillColor(42);
   h->SetXTitle("x (cm)");
   h->Draw();
   TCanvas *c3 = new TCanvas("c3","Charge distribution",400,10,600,700);
   charge->SetFillColor(42);
   charge->SetXTitle("ADC units");
   charge->Draw();

}
