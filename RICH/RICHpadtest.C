void RICHpadtest (Int_t evNumber1=0,Int_t evNumber2=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////


    Int_t NpadX = 252;                 // number of pads on X
    Int_t NpadY = 374;                 // number of pads on Y
    
    Int_t Pad[252][374];
    for (Int_t i=0;i<NpadX;i++) {
	for (Int_t j=0;j<NpadY;j++) {
	    Pad[i][j]=0;
	}
    }
    


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
   for (int nev=0; nev<= evNumber2; nev++) {
      Int_t nparticles = gAlice->GetEvent(nev);
      //cout<<"nev  "<<nev<<endl;
      printf ("Number of Events   : %d\n",nev);
      //cout<<"nparticles  "<<nparticles<<endl;
      printf ("Number of particles: %d\n",nparticles);
      if (nev < evNumber1) continue;
      if (nparticles <= 0) return;
     
// Get pointers to RICH detector and Hits containers

      AliRICH *RICH  = gAlice->GetDetector("RICH");
      TTree *TH = gAlice->TreeH();
      Int_t ntracks = TH->GetEntries();

// Start loop on tracks in the hits containers

      //      Int_t Nh=0;
      Int_t Nc=0;
      for (Int_t track=0; track<ntracks;track++) {
	//          printf("ntracks %d\n",ntracks);
          gAlice->ResetHits();
          Int_t nbytes += TH->GetEvent(track);
          if (RICH)  {
              TClonesArray *Clusters = RICH->Clusters(); // Cluster branch
              TClonesArray *Hits = RICH->Hits();         // Hits branch
	      //printf("%d %d \n",Clusters,Hits);
          }
          //see hits distribution
          Int_t nhits = Hits->GetEntriesFast();
          if (nhits) Nh+=nhits;
	  //          printf("nhits %d\n",nhits);
          for (Int_t hit=0;hit<nhits;hit++) {
              mHit = (AliRICHhit*) Hits->UncheckedAt(hit);
              Int_t nch  = mHit->fChamber;              // chamber number
              Float_t x  = mHit->fX;                    // x-pos of hit
              Float_t y  = mHit->fY;                    // y-pos
              // Fill the histograms
	      if( nch==1) {
		 Float_t rhit=TMath::Sqrt(x*x+y*y);
                 if( rhit<= 55 ) Nh1+=nhits;
                  h->Fill(x,y,(float) 1);
              }
          }
          //    see signal distribution
          Int_t nclust = Clusters->GetEntriesFast();
	  //	   printf("nclust %d\n",nclust);
          if (nclust) {
            Nc+=nclust;
            for (Int_t hit=0;hit<nclust;hit++) {
	      padHit = (AliRICHcluster*) Clusters->UncheckedAt(hit);
	      Int_t nchamber = padHit->fChamber;     // chamber number
	      Int_t nhit = padHit->fHitNumber;          // hit number
	      Int_t qtot = padHit->fQ;                // charge
	      Int_t ipx  = padHit->fPadX;               // pad number on X
	      Int_t ipy  = padHit->fPadY;               // pad number on Y
	      Int_t iqpad  = padHit->fQpad;           // charge per pad
	      Int_t rpad  = padHit->fRpad;            // R-position of pad
              if (qtot > 0 && nchamber==1){
                 charge->Fill(qtot,(float) 1);
              }
      	      if(rpad <= 55 && nchamber==1) {
	      //	      if(nchamber==1) {
                  Pad[ipx+126][ipy+187]+=(iqpad);
                  hc->Fill(ipx,ipy,(float) iqpad);
	      }	 	
            }
          }
       }
      //      cout<<"Nh  "<<Nh<<endl;
      //cout<<"Nc  "<<Nc<<endl;
      printf ("Cluster Number: %d\n",Nc);
   }
//Create a canvas, set the view range, show histograms
   TCanvas *c1 = new TCanvas("c1","Alice RICH pad hits",400,10,600,700);
   hc->SetXTitle("ix (npads)");
   hc->Draw();
   TCanvas *c2 = new TCanvas("c2","Alice RICH hits",400,10,600,700);
   h->SetFillColor(42);
   h->SetXTitle("x (cm)");
   h->Draw();
   TCanvas *c3 = new TCanvas("c3","Charge distribution",400,10,600,700);
   charge->SetFillColor(42);
   charge->SetXTitle("ADC units");
   charge->Draw();

   
   // calculate the number of pads which give a signal


    Int_t Np=0;
    for (Int_t i=0;i< NpadX;i++) {
      for (Int_t j=0;j< NpadY;j++) {
         if (Pad[i][j]>=6){
             Np+=1;
	     //             cout<<"i, j  "<<i<<"  "<<j<<endl;
         }
      }
    }

    //cout<<"The total number of pads which give a signal "<<Np<<endl; 
    //cout<<"Nh  "<<Nh<<endl;
    //cout<<"Nh1  "<<Nh1<<endl;
    printf("The total number of pads which give a signal: %d %d\n",Nh,Nh1);

}
