void RICHpadtest (Int_t diaglevel,Int_t evNumber1=0,Int_t evNumber2=0) 
{

// Int_t diaglevel=3;         // 1->Hits, 2->Spectra, 3->Statistics 

/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////


    Int_t NpadX = 162;                 // number of pads on X
    Int_t NpadY = 162;                 // number of pads on Y
    
    Int_t Pad[162][162];
    for (Int_t i=0;i<NpadX;i++) {
	for (Int_t j=0;j<NpadY;j++) {
	    Pad[i][j]=0;
	}
    }
    gClassTable->GetID("AliRun");


// Dynamically link some shared libs
 
    if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
    }
    else {
      //delete gAlice;
      gAlice = 0;
    }

    gAlice=0;
    
// Connect the Root Galice file containing Geometry, Kine and Hits
    
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (!file) file = new TFile("galice.root");
    
// Get AliRun object from file or create it if not on file
    
    if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
    else {
      delete gAlice;
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }

//  Create some histograms

   Int_t xmin= -NpadX/2;  
   Int_t xmax=  NpadX/2;
   Int_t ymin= -NpadY/2;
   Int_t ymax=  NpadY/2;

   TH2F *hc0 = new TH2F("hc0","Zoom on center of Chamber 1",100,-50,50,100,-50,50);
   TH2F *hc1 = new TH2F("hc1","Chamber 1 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc2 = new TH2F("hc2","Chamber 2 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc3 = new TH2F("hc3","Chamber 3 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc4 = new TH2F("hc4","Chamber 4 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc5 = new TH2F("hc5","Chamber 5 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc6 = new TH2F("hc6","Chamber 6 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc7 = new TH2F("hc7","Chamber 7 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *h = new TH2F("h","Detector hit distribution",150,-300,300,150,-300,300);
   TH1F *Clcharge = new TH1F("Clcharge","Cluster Charge Distribution",500,0.,500.);
   TH2F *cerenkov = new TH2F("cerenkov","Cerenkov hit distribution",150,-300,300,150,-300,300);
   TH1F *ckovangle = new TH1F("ckovangle","Cerenkov angle per photon",200,.5,1);
   TH1F *hckphi = new TH1F("hckphi","Cerenkov phi angle per photon",620,-3.1,3.1);
   TH2F *feedback = new TH2F("feedback","Feedback hit distribution",150,-300,300,150,-300,300);
   TH2F *mip = new TH2F("mip","Mip hit distribution",150,-300,300,150,-300,300);
   TH1F *mother = new TH1F("mother","Cerenkovs per Mip",75,0.,75.);
   TH1F *radius = new TH1F("radius","Mean distance to Mip",100,0.,20.);
   TH1F *phspectra1 = new TH1F("phspectra","Photon Spectra",200,5.,10.);
   TH1F *phspectra2 = new TH1F("phspectra","Photon Spectra",200,5.,10.);
   TH1F *totalphotons = new TH1F("totalphotons","Produced Photons per Mip",100,200,700.);
   TH1F *feedbacks = new TH1F("feedbacks","Produced Feedbacks per Mip",50,0.5,50.);
   TH1F *padnumber = new TH1F("padnumber","Number of pads per cluster",50,-0.5,50.);
   TH1F *padsev = new TH1F("padsev","Number of pads hit per event",50,0.5,100.);
   TH1F *clusev = new TH1F("clusev","Number of clusters per event",50,0.5,50.);
   TH1F *photev = new TH1F("photev","Number of photons per event",50,0.5,50.);
   TH1F *feedev = new TH1F("feedev","Number of feedbacks per event",50,0.5,50.);
   TH1F *padsmip = new TH1F("padsmip","Number of pads per event inside MIP region",50,0.5,50.);
   TH1F *padscl = new TH1F("padscl","Number of pads per event from cluster count",50,0.5,100.);
   TH1F *pionspectra = new TH1F("pionspectra","Pion Spectra",200,.5,10.);
   TH1F *protonspectra = new TH1F("protonspectra","Proton Spectra",200,.5,10.);
   TH1F *kaonspectra = new TH1F("kaonspectra","Kaon Spectra",100,.5,10.);
   TH1F *kaonspectra = new TH1F("kaonspectra","Kaon Spectra",100,.5,10.);
   TH1F *chargedspectra = new TH1F("chargedspectra","Charged particles above 1 GeV Spectra",100,.5,10.);
   TH1F *hitsX = new TH1F("digitsX","Distribution of hits along x-axis",200,-300,300);
   TH1F *hitsY = new TH1F("digitsY","Distribution of hits along z-axis",200,-300,300);
   TH1F *hitsPhi = new TH1F("hitsPhi","Distribution of phi angle of incidence",100,-180,180);
   TH1F *hitsTheta = new TH1F("hitsTheta","Distribution of Theta angle of incidence",100,0,15);
   TH1F *Omega = new TH1F("omega","Reconstructed Cerenkov angle per track",200,.5,1);
   TH1F *Theta = new TH1F("theta","Reconstructed theta incidence angle per track",200,0,15);
   TH1F *Phi = new TH1F("phi","Reconstructed phi incidence per track",200,-180,180);
   
   
   

//   Start loop over events 

   Int_t Nh=0;
   Int_t pads=0;
   Int_t Nh1=0;
   Int_t mothers[80000];
   Int_t mothers2[80000];
   Float_t mom[3];
   Int_t nraw=0;
   Int_t phot=0;
   Int_t feed=0;
   Int_t padmip=0;
   for (Int_t i=0;i<100;i++) mothers[i]=0;
   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       

       //cout<<"nev  "<<nev<<endl;
       printf ("\nProcessing event: %d\n",nev);
       //cout<<"nparticles  "<<nparticles<<endl;
       printf ("Particles       : %d\n",nparticles);
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       
// Get pointers to RICH detector and Hits containers
       
       AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
       Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
       gAlice->TreeR()->GetEvent(nent-1);
       TClonesArray *Rawclusters = RICH->RawClustAddress(2);    //  Raw clusters branch
       //printf ("Rawclusters:%p",Rawclusters);
       Int_t nrawclusters = Rawclusters->GetEntriesFast();
       //printf (" nrawclusters:%d\n",nrawclusters);
       gAlice->TreeR()->GetEvent(nent-1);
       TClonesArray *RecHits = RICH->RecHitsAddress(2);
       Int_t nrechits = RecHits->GetEntriesFast();
       //printf (" nrechits:%d\n",nrechits);
       TTree *TH = gAlice->TreeH(); 
       Int_t ntracks = TH->GetEntries();


       
       Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
       gAlice->TreeD()->GetEvent(nent-1);
       

       TClonesArray *Digits = RICH->DigitsAddress(2);    //  Raw clusters branch
       Int_t ndigits = Digits->GetEntriesFast();
       //printf("Digits:%d\n",ndigits);
       padsev->Fill(ndigits,(float) 1);

       for (Int_t ich=0;ich<7;ich++)
	 {
	   TClonesArray *Digits = RICH->DigitsAddress(ich);    //  Raw clusters branch
	   Int_t ndigits = Digits->GetEntriesFast();
	   //printf("Digits:%d\n",ndigits);
	   padsev->Fill(ndigits,(float) 1); 
	   if (ndigits) {
	     for (Int_t hit=0;hit<ndigits;hit++) {
	       dHit = (AliRICHDigit*) Digits->UncheckedAt(hit);
	       //Int_t nchamber = padHit->fChamber;     // chamber number
	       //Int_t nhit = dHit->fHitNumber;          // hit number
	       Int_t qtot = dHit->fSignal;                // charge
	       Int_t ipx  = dHit->fPadX;               // pad number on X
	       Int_t ipy  = dHit->fPadY;               // pad number on Y
	       //Int_t iqpad  = dHit->fQpad;           // charge per pad
	       //Int_t rpad  = dHit->fRSec;            // R-position of pad
	       //printf ("Pad hit, PadX:%d, PadY:%d\n",ipx,ipy);
	       if( ipx<=100 && ipy <=100 && ich==0) hc0->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==0) hc1->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==1) hc2->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==2) hc3->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==3) hc4->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==4) hc5->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==5) hc6->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==6) hc7->Fill(ipx,ipy,(float) qtot);
	     }
	   }
	 }
       
// Start loop on tracks in the hits containers
       Int_t Nc=0;
       for (Int_t track=0; track<ntracks;track++) {
	   printf ("Processing Track: %d\n",track);
	   gAlice->ResetHits();
	   Int_t nbytes += TH->GetEvent(track);
	   if (RICH)  {
	       //RICH->ResetRawClusters();
	       TClonesArray *PadHits = RICH->PadHits();      // Cluster branch
	       TClonesArray *Hits = RICH->Hits();            // Hits branch
	       TClonesArray *Cerenkovs = RICH->Cerenkovs();  // Cerenkovs branch
	   }
	   //see hits distribution
	   

	   Int_t nhits = Hits->GetEntriesFast();
	   if (nhits) Nh+=nhits;
	   //printf("nhits %d\n",nhits);
	   for (Int_t hit=0;hit<nhits;hit++) {
              mHit = (AliRICHHit*) Hits->UncheckedAt(hit);
              Int_t nch  = mHit->fChamber;              // chamber number
              Float_t x  = mHit->fX;                    // x-pos of hit
              Float_t y  = mHit->fZ;                    // y-pos
	      Float_t phi = mHit->fPhi;                 //Phi angle of incidence
	      Float_t theta = mHit->fTheta;             //Theta angle of incidence
	      Int_t index = mHit->fTrack;
	      Int_t particle = mHit->fParticle;        
	      Int_t freon = mHit->fLoss;    

	     hitsX->Fill(x,(float) 1);
	     hitsY->Fill(y,(float) 1);

	      //printf("Particle:%d\n",particle);
	      
	      TParticle *current = (TParticle*)(*gAlice->Particles())[index];
	      //printf("Particle type: %d\n",current->GetPdgCode());
	      if (TMath::Abs(particle) < 50000000)
		{
		  mip->Fill(x,y,(float) 1);
		  if (current->Energy() - current->GetCalcMass()>1 && freon==1)
		    {
		      hitsPhi->Fill(phi,(float) 1);
		      hitsTheta->Fill(theta,(float) 1);
		      //printf("Theta:%f, Phi:%f\n",theta,phi);
		    }
		}
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  pionspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      if (TMath::Abs(particle)==2212)
		{
		  protonspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  kaonspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      if(TMath::Abs(particle)==211 || TMath::Abs(particle)==2212 || TMath::Abs(particle)==321)
		{
		  if (current->Energy() - current->GetCalcMass()>1)
		    chargedspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      //printf("Hits:%d\n",hit);
	      //printf ("Chamber number:%d x:%f y:%f\n",nch,x,y);
              // Fill the histograms
	      Nh1+=nhits;
	      h->Fill(x,y,(float) 1);
		  //}
              //}
          }
	   
	   Int_t ncerenkovs = Cerenkovs->GetEntriesFast();

	   if (ncerenkovs) {
	       for (Int_t hit=0;hit<ncerenkovs;hit++) {
		   cHit = (AliRICHCerenkov*) Cerenkovs->UncheckedAt(hit);
		   Int_t nchamber = cHit->fChamber;     // chamber number
		   Int_t index =    cHit->fTrack;
		   Int_t pindex =   cHit->fIndex;
		   Int_t cx  =      cHit->fX;                // x-position
		   Int_t cy  =      cHit->fZ;                // y-position
		   Int_t cmother =  cHit->fCMother;      // Index of mother particle
		   Int_t closs =    cHit->fLoss;           // How did the particle get lost? 
		  //printf ("Cerenkov hit, X:%d, Y:%d\n",cx,cy); 

		  
		 
		   TParticle *current = (TParticle*)(*gAlice->Particles())[index];
		   
		   if (current->GetPdgCode() == 50000051)
		   {
		       if (closs==4)
		       {
			   feedback->Fill(cx,cy,(float) 1);
			   feed++;
		       }
		   }
		   if (current->GetPdgCode() == 50000050)
		   {
		     if (closs==4)
		       {
			 cerenkov->Fill(cx,cy,(float) 1);
			 
			 
			 
			 TParticle *MIP = (TParticle*)(*gAlice->Particles())[cmother];
			 mipHit = (AliRICHHit*) Hits->UncheckedAt(0);
			 mom[0] = current->Px();
			 mom[1] = current->Py();
			 mom[2] = current->Pz();
			 Float_t energyckov = current->Energy();
			 /*mom[0] = cHit->fMomX;
			   mom[1] = cHit->fMomZ;
			   mom[2] = cHit->fMomY;*/
			 Float_t energymip = MIP->Energy();
			 Float_t Mip_px = mipHit->fMomX;
			 Float_t Mip_py = mipHit->fMomY;
			 Float_t Mip_pz = mipHit->fMomZ;
			 /*Float_t Mip_px = MIP->Px();
			   Float_t Mip_py = MIP->Py();
			   Float_t Mip_pz = MIP->Pz();*/
			 
			 
			 
			 Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
			 Float_t rt = TMath::Sqrt(r);
			 Float_t Mip_r = Mip_px*Mip_px + Mip_py*Mip_py + Mip_pz*Mip_pz;	
			 Float_t Mip_rt = TMath::Sqrt(Mip_r);
			 Float_t coscerenkov = (mom[0]*Mip_px + mom[1]*Mip_py + mom[2]*Mip_pz)/(rt*Mip_rt);
			 Float_t cherenkov = TMath::ACos(coscerenkov);
			 ckovangle->Fill(cherenkov,(float) 1);                           //Cerenkov angle calculus
			 //printf("Cherenkov: %f\n",cherenkov);
			 Float_t ckphi=TMath::ATan2(mom[0], mom[2]);
			 hckphi->Fill(ckphi,(float) 1);
			 
			 
			 Float_t mix = MIP->Vx();
			 Float_t miy = MIP->Vy();
			 Float_t mx = mipHit->fX;
			 Float_t my = mipHit->fZ;
			 //printf("FX %e, FY %e, VX %e, VY %e\n",cx,cy,mx,my);
			 Float_t dx = cx - mx;
			 Float_t dy = cy - my;
			 //printf("Dx:%f, Dy:%f\n",dx,dy);
			 Float_t final_radius = TMath::Sqrt(dx*dx+dy*dy);
			 //printf("Final radius:%f\n",final_radius);
			 radius->Fill(final_radius,(float) 1);
			 
			 phspectra1->Fill(energyckov*1e9,(float) 1);
			 phot++;
		       }
		     else
		       phspectra2->Fill(energyckov*1e9,(float) 1);
		     for (Int_t nmothers=0;nmothers<=ntracks;nmothers++){
		       if (cmother == nmothers){
			 if (closs == 4)
			   mothers2[cmother]++;
			 mothers[cmother]++;
		       }
		     } 
		   }
	       }
	   }
	   
	   if (nrawclusters) {
	       for (Int_t hit=0;hit<nrawclusters;hit++) {
		   rcHit = (AliRICHRawCluster*) Rawclusters->UncheckedAt(hit);
		   //Int_t nchamber = rcHit->fChamber;     // chamber number
		   //Int_t nhit = cHit->fHitNumber;        // hit number
		   Int_t qtot = rcHit->fQ;                 // charge
		   Int_t fx  =  rcHit->fX;                 // x-position
		   Int_t fy  =  rcHit->fY;                 // y-position
		   Int_t type = rcHit->fCtype;             // cluster type ?   
		   Int_t mult = rcHit->fMultiplicity;      // How many pads form the cluster
		   pads += mult;
		   if (qtot > 0) {
		       if (fx>-4 && fx<4 && fy>-4 && fy<4) {
			   padmip+=mult;
		       } else {
			   padnumber->Fill(mult,(float) 1);
			   nraw++;
			   if (mult<4) Clcharge->Fill(qtot,(float) 1);
		       }
		   }
	       }
	   }

	   if(nrechits)
	     {
	       for (Int_t hit=0;hit<nrechits;hit++) {
		 recHit = (AliRICHRecHit*) RecHits->UncheckedAt(hit);
		 Float_t r_omega = recHit->fOmega;                  // Cerenkov angle
		 Float_t r_theta = recHit->fTheta;                  // Theta angle of incidence
		 Float_t r_phi   = recHit->fPhi;                    // Phi angle if incidence
		 
		 Omega->Fill(r_omega,(float) 1);
		 Theta->Fill(r_theta*180/TMath::Pi(),(float) 1);
		 Phi->Fill(r_phi*180/TMath::Pi(),(float) 1);
		 
		 //printf("Omega: %f, Theta: %f, Phi: %f\n",r_omega,r_theta,r_phi);
	       }
	     }
       }
       
       for (Int_t nmothers=0;nmothers<ntracks;nmothers++){
	   totalphotons->Fill(mothers[nmothers],(float) 1);
	   mother->Fill(mothers2[nmothers],(float) 1);
	   //printf ("Entries in %d : %d\n",nmothers, mothers[nmothers]);
       }
       
       clusev->Fill(nraw,(float) 1);
       photev->Fill(phot,(float) 1);
       feedev->Fill(feed,(float) 1);
       padsmip->Fill(padmip,(float) 1);
       padscl->Fill(pads,(float) 1);
       //printf("Photons:%d\n",phot);
       phot = 0;
       feed = 0;
       pads = 0;
       nraw=0;
       padmip=0;

   }
   
   //Create canvases, set the view range, show histograms

   switch(diaglevel)
     {
     case 1:

       if (ndigits)
	 {
	   TCanvas *c1 = new TCanvas("c1","Alice RICH pad hits",50,10,1200,700);
	   c1->Divide(4,2);
	   c1->cd(1);
	   hc1->SetXTitle("ix (npads)");
	   hc1->Draw();
	   c1->cd(2);
	   hc2->SetXTitle("ix (npads)");
	   hc2->Draw();
	   c1->cd(3);
	   hc3->SetXTitle("ix (npads)");
	   hc3->Draw();
	   c1->cd(4);
	   hc4->SetXTitle("ix (npads)");
	   hc4->Draw();
	   c1->cd(5);
	   hc5->SetXTitle("ix (npads)");
	   hc5->Draw();
	   c1->cd(6);
	   hc6->SetXTitle("ix (npads)");
	   hc6->Draw();
	   c1->cd(7);
	   hc7->SetXTitle("ix (npads)");
	   hc7->Draw();
	   c1->cd(8);
	   hc0->SetXTitle("ix (npads)");
	   hc0->Draw();
	 }
//
       TCanvas *c4 = new TCanvas("c4","Hits per type",400,10,600,700);
       c4->Divide(2,2);
       
       c4->cd(1);
       feedback->SetFillColor(42);
       feedback->SetXTitle("x (pads)");
       feedback->SetYTitle("y (pads)");
       feedback->Draw();
       
       c4->cd(2);
       mip->SetFillColor(42);
       mip->SetXTitle("x (pads)");
       mip->SetYTitle("y (pads)");
       mip->Draw();
       
       c4->cd(3);
       cerenkov->SetFillColor(42);
       cerenkov->SetXTitle("x (pads)");
       cerenkov->SetYTitle("y (pads)"); 
       cerenkov->Draw();
       
       c4->cd(4);
       h->SetFillColor(42);
       h->SetXTitle("x (cm)");
       h->SetYTitle("y (cm)");
       h->Draw();

       TCanvas *c10 = new TCanvas("c10","Hits distribution",400,10,600,350);
       c10->Divide(2,1);
       
       c10->cd(1);
       hitsX->SetFillColor(42);
       hitsX->SetXTitle("(cm)");
       hitsX->Draw();
       
       c10->cd(2);
       hitsY->SetFillColor(42);
       hitsY->SetXTitle("(cm)");
       hitsY->Draw();
       
      
       break;
//
     case 2:
       
       TCanvas *c6 = new TCanvas("c6","Photon Spectra",50,10,600,350);
       c6->Divide(2,1);
       
       c6->cd(1);
       phspectra2->SetFillColor(42);
       phspectra2->SetXTitle("energy (eV)");
       phspectra2->Draw();
       c6->cd(2);
       phspectra1->SetFillColor(42);
       phspectra1->SetXTitle("energy (eV)");
       phspectra1->Draw();
       
       TCanvas *c9 = new TCanvas("c9","Particles Spectra",400,10,600,700);
       c9->Divide(2,2);
       
       c9->cd(1);
       pionspectra->SetFillColor(42);
       pionspectra->SetXTitle("(GeV)");
       pionspectra->Draw();
       
       c9->cd(2);
       protonspectra->SetFillColor(42);
       protonspectra->SetXTitle("(GeV)");
       protonspectra->Draw();
       
       c9->cd(3);
       kaonspectra->SetFillColor(42);
       kaonspectra->SetXTitle("(GeV)");
       kaonspectra->Draw();
       
       c9->cd(4);
       chargedspectra->SetFillColor(42);
       chargedspectra->SetXTitle("(GeV)");
       chargedspectra->Draw();

       break;
       
     case 3:
       
       if (nrawclusters) {
	 TCanvas *c3=new TCanvas("c3","Clusters Statistics",400,10,600,700);
	 c3->Divide(2,2);
	 
	 c3->cd(1);
	 Clcharge->SetFillColor(42);
	 Clcharge->SetXTitle("ADC units");
	 Clcharge->Draw();
	 
	 c3->cd(2);
	 padnumber->SetFillColor(42);
	 padnumber->SetXTitle("(counts)");
	 padnumber->Draw();
	 
	 c3->cd(3);
	 clusev->SetFillColor(42);
	 clusev->SetXTitle("(counts)");
	 clusev->Draw();

	 c3->cd(4);
	 padsmip->SetFillColor(42);
	 padsmip->SetXTitle("(counts)");
	 padsmip->Draw(); 
       }

       if (nev<1)
	 {
	   TCanvas *c11 = new TCanvas("c11","Cherenkov per Mip",400,10,600,700);
	   mother->SetFillColor(42);
	   mother->SetXTitle("counts");
	   mother->Draw();
	 }
       
       TCanvas *c2 = new TCanvas("c2","Angles of incidence",50,10,600,700);
       c2->Divide(2,2);
       
       c2->cd(1);
       hitsPhi->SetFillColor(42);
       hitsPhi->Draw();
       c2->cd(2);
       hitsTheta->SetFillColor(42);
       hitsTheta->Draw();
       c2->cd(3);
       Phi->SetFillColor(42);
       Phi->Draw();
       c2->cd(4);
       Theta->SetFillColor(42);
       Theta->Draw();
       
       
       TCanvas *c5 = new TCanvas("c5","Ring Statistics",50,10,600,700);
       c5->Divide(2,2);
       
       c5->cd(1);
       ckovangle->SetFillColor(42);
       ckovangle->SetXTitle("angle (radians)");
       ckovangle->Draw();
       
       c5->cd(2);
       radius->SetFillColor(42);
       radius->SetXTitle("radius (cm)");
       radius->Draw();
       
       c5->cd(3);
       Omega->SetFillColor(42);
       Omega->SetXTitle("angle (radians)");
       Omega->Draw();

       TCanvas *c7 = new TCanvas("c7","Production Statistics",400,10,600,700);
       c7->Divide(2,2);
       
       c7->cd(1);
       totalphotons->SetFillColor(42);
       totalphotons->SetXTitle("Photons (counts)");
       totalphotons->Draw();
       
       c7->cd(2);
       photev->SetFillColor(42);
       photev->SetXTitle("(counts)");
       photev->Draw();
       
       c7->cd(3);
       feedev->SetFillColor(42);
       feedev->SetXTitle("(counts)");
       feedev->Draw();

       c7->cd(4);
       padsev->SetFillColor(42);
       padsev->SetXTitle("(counts)");
       padsev->Draw();

       break;
       
     }
       

   // calculate the number of pads which give a signal


   Int_t Np=0;
   for (Int_t i=0;i< NpadX;i++) {
       for (Int_t j=0;j< NpadY;j++) {
	   if (Pad[i][j]>=6){
	       Np+=1;
	   }
       }
   }
   //printf("The total number of pads which give a signal: %d %d\n",Nh,Nh1);
   printf("End of macro\n");
}


