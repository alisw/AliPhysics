Int_t diaglevel=2;         // 1->Hits, 2->Spectra, 3->Statistics 


void RICHpadtestC (Int_t evNumber1=0,Int_t evNumber2=0) 
{
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

   Int_t xmin= -NpadX/2;  
   Int_t xmax=  NpadX/2;
   Int_t ymin= -NpadY/2;
   Int_t ymax=  NpadY/2;

   /*TH2F *hc1 = new TH2F("hc1","Chamber 1 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc2 = new TH2F("hc2","Chamber 2 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc3 = new TH2F("hc3","Chamber 3 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc4 = new TH2F("hc4","Chamber 4 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc5 = new TH2F("hc5","Chamber 5 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc6 = new TH2F("hc6","Chamber 6 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc7 = new TH2F("hc7","Chamber 7 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *h = new TH2F("h","Detector hit distribution",150,-300,300,150,-300,300);
   TH1F *Clcharge = new TH1F("Clcharge","Cluster Charge Distribution",500,0.,500.);
   TH2F *cerenkov = new TH2F("cerenkov","Cerenkov hit distribution",150,-300,300,150,-300,300);
   TH1F *ckovangle = new TH1F("ckovangle","Cerenkov angle per photon",200,.6,.85);
   TH1F *hckphi = new TH1F("hckphi","Cerenkov phi angle per photon",620,-3.1,3.1);
   TH2F *feedback = new TH2F("feedback","Feedback hit distribution",150,-300,300,150,-300,300);
   TH2F *mip = new TH2F("mip","Mip hit distribution",150,-300,300,150,-300,300);
   TH1F *mother = new TH1F("mother","Cerenkovs per Mip",75,0.,75.);
   TH1F *radius = new TH1F("radius","Mean distance to Mip",200,0.,20.);
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
   TH1F *padscl = new TH1F("padscl","Number of pads per event from cluster count",50,0.5,100.);*/
   TH1F *pionspectra1 = new TH1F("pionspectra1","Pion Spectra",200,-4,2);
   TH1F *pionspectra2 = new TH1F("pionspectra2","Pion Spectra",200,-4,2);
   TH1F *pionspectra3 = new TH1F("pionspectra3","Pion Spectra",200,-4,2);
   TH1F *protonspectra1 = new TH1F("protonspectra1","Proton Spectra",200,-4,2);
   TH1F *protonspectra2 = new TH1F("protonspectra2","Proton Spectra",200,-4,2);
   TH1F *protonspectra3 = new TH1F("protonspectra3","Proton Spectra",200,-4,2);
   TH1F *kaonspectra1 = new TH1F("kaonspectra1","Kaon Spectra",100,-4,2);
   TH1F *kaonspectra2 = new TH1F("kaonspectra2","Kaon Spectra",100,-4,2);
   TH1F *kaonspectra3 = new TH1F("kaonspectra3","Kaon Spectra",100,-4,2);
   TH1F *electronspectra1 = new TH1F("electronspectra1","Electron Spectra",100,-4,2);
   TH1F *electronspectra2 = new TH1F("electronspectra2","Electron Spectra",100,-4,2);
   TH1F *electronspectra3 = new TH1F("electronspectra3","Electron Spectra",100,-4,2);
   TH1F *muonspectra1 = new TH1F("muonspectra1","Muon Spectra",100,-4,2);
   TH1F *muonspectra2 = new TH1F("muonspectra2","Muon Spectra",100,-4,2);
   TH1F *muonspectra3 = new TH1F("muonspectra3","Muon Spectra",100,-4,2);
   TH1F *neutronspectra1 = new TH1F("neutronspectra1","Neutron Spectra",100,-4,2);
   TH1F *neutronspectra2 = new TH1F("neutronspectra2","Neutron Spectra",100,-4,2);
   TH1F *neutronspectra3 = new TH1F("neutronspectra2","Neutron Spectra",100,-4,2);
   TH1F *chargedspectra1 = new TH1F("chargedspectra1","Charged particles above 1 GeV Spectra",100,-1,3);
   TH1F *chargedspectra2 = new TH1F("chargedspectra2","Charged particles above 1 GeV Spectra",100,-1,3);
   TH1F *chargedspectra3 = new TH1F("chargedspectra2","Charged particles above 1 GeV Spectra",100,-1,3);
/*   TH1F *hitsX = new TH1F("digitsX","Distribution of hits along x-axis",200,-300,300);
   TH1F *hitsY = new TH1F("digitsY","Distribution of hits along z-axis",200,-300,300);*/
   TH2F *production = new TH2F("production","Mother production vertices",100,-300,300,100,0,600);
   
   
   

//   Start loop over events 

   Int_t Nh=0;
   Int_t pads=0;
   Int_t Nh1=0;
   //Int_t mothers[100000];
   //Int_t mothers2[100000];
   Float_t mom[3];
   //Float_t random;
   Int_t nraw=0;
   Int_t phot=0;
   Int_t feed=0;
   Int_t padmip=0;
   Int_t pion=0, kaon=0, proton=0, electron=0, positron=0, neutron=0, highneutrons=0, muon=0;
   Int_t chargedpions=0,primarypions=0,highprimarypions=0,chargedkaons=0,primarykaons=0,highprimarykaons=0;
   Int_t chargedmuons=0, photons=0, primaryphotons=0, highprimaryphotons=0;

   TRandom random;

   //for (Int_t i=0;i<100;i++) mothers[i]=0;
   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       

       //cout<<"nev  "<<nev<<endl;
       printf ("Event number       : %d\n",nev);
       //cout<<"nparticles  "<<nparticles<<endl;
       printf ("Number of particles: %d\n",nparticles);
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
       TTree *TH = gAlice->TreeH(); 
       Int_t ntracks = TH->GetEntries();


       
     /*  Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
       gAlice->TreeD()->GetEvent(nent-1);
       

       TClonesArray *Digits = RICH->DigitsAddress(2);    //  Raw clusters branch
       Int_t ndigits = Digits->GetEntriesFast();
       //printf("Digits:%d\n",ndigits);
       padsev->Fill(ndigits,(float) 1);*/

      /* for (Int_t ich=0;ich<7;ich++)
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
	       if( ipx<=162 && ipy <=162 && ich==0) hc1->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==1) hc2->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==2) hc3->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==3) hc4->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==4) hc5->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==5) hc6->Fill(ipx,ipy,(float) qtot);
	       if( ipx<=162 && ipy <=162 && ich==6) hc7->Fill(ipx,ipy,(float) qtot);
	     }
	   }
	 }*/
       
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
	      Float_t z  = mHit->fY;
	      Int_t index = mHit->fTrack;
	      Int_t particle = mHit->fParticle;    
	      Float_t R;

	      //hitsX->Fill(x,(float) 1);
	      //hitsY->Fill(y,(float) 1);

	      //printf("Particle:%d\n",particle);
	      
	      TParticle *current = (TParticle*)(*gAlice->Particles())[index];
	      
	      R=TMath::Sqrt(current->Vx()*current->Vx() + current->Vy()*current->Vy());

	      //printf("Particle type: %d\n",current->GetPdgCode());
	      if (TMath::Abs(particle) < 50000051)
		{
		  //if (TMath::Abs(particle) == 50000050 || TMath::Abs(particle) == 2112)
		  if (TMath::Abs(particle) == 2112 || TMath::Abs(particle) == 50000050)
		    {
		      //gMC->Rndm(&random, 1);
		      if (random->Rndm() < .1)
			production->Fill(current->Vz(),R,(float) 1);
		      if (TMath::Abs(particle) == 50000050)
			{
			  photons +=1;
			  if (R<.005)
			    {
			      primaryphotons +=1;
			      if (current->Energy()>0.001)
				highprimaryphotons +=1;
			    }
			}	
		      if (TMath::Abs(particle) == 2112)
			{
			  neutron +=1;
			  if (current->Energy()>0.0001)
			    highneutrons +=1;
			}
		    }
		  else 
		    {
		      production->Fill(current->Vz(),R,(float) 1);
		      printf("Adding %d at %f\n",particle,R);
		    }
		  //mip->Fill(x,y,(float) 1);
		}
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  pionspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>.005 && current->Vy()>.005 && current->Vz()>.005)
		    pionspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>2.5 && R<4.5)
		    {
		    pionspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\R:%f\n\n\n\n\n\n\n\n\n",R);
		    }
		  printf("Pion mass: %e\n",current->GetCalcMass());
		  pion +=1;
		  if (TMath::Abs(particle)==211)
		    {
		      chargedpions +=1;
		      if (R<.005)
			{
			  primarypions +=1;
			  if (current->Energy()>1)
			    highprimarypions +=1;
			}
		    }	
		}
	      if (TMath::Abs(particle)==2212)
		{
		  protonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>.005 && current->Vy()>.005 && current->Vz()>.005)
		    protonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>3 && R<4.3)
		    protonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //printf("\n\n\n\n\n\n\nProton mass: %e\n\n\n\n\n\n\n\n\n",current->GetCalcMass());
		  proton +=1;
		}
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  kaonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>.005 && current->Vy()>.005 && current->Vz()>.005)
		    kaonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>2.5 && R<4.5)
		    kaonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  printf("Kaon mass: %e\n",current->GetCalcMass());
		  kaon +=1;
		  if (TMath::Abs(particle)==321)
		    {
		      chargedkaons +=1;
		      if (R<.005)
			{
			  primarykaons +=1;
			  if (current->Energy()>1)
			    highprimarykaons +=1;
			}
		    }
		}
	      if (TMath::Abs(particle)==11)
		{
		  electronspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>.005 && current->Vy()>.005 && current->Vz()>.005)
		    electronspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>2.5 && R<4.5)
		    electronspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  printf("Electron mass: %e\n",current->GetCalcMass());
		  if (particle == -11)
		    electron +=1;
		  if (particle == 11)
		    positron +=1;
		}
	      if (TMath::Abs(particle)==13)
		{
		  muonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>.005 && current->Vy()>.005 && current->Vz()>.005)
		    muonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>2.5 && R<4.5)
		    muonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  printf("Muon mass: %e\n",current->GetCalcMass());
		  muon +=1;
		}
	      if (TMath::Abs(particle)==2112)
		{
		  neutronspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>.005 && current->Vy()>.005 && current->Vz()>.005)
		    neutronspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>2.5 && R<4.5)
		    {
		      neutronspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		      printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\R:%f\n\n\n\n\n\n\n\n\n",R);
		    }
		  printf("Neutron mass: %e\n",current->GetCalcMass());
		  neutron +=1;
		}
	      if(TMath::Abs(particle)==211 || TMath::Abs(particle)==2212 || TMath::Abs(particle)==321)
		{
		  if (current->Energy()-current->GetCalcMass()>1)
		    {
		      chargedspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		      if (current->Vx()>.005 && current->Vy()>.005 && current->Vz()>.005)
			chargedspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		      if (R>2.5 && R<4.5)
			chargedspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		}
	      //printf("Hits:%d\n",hit);
	      //printf ("Chamber number:%d x:%f y:%f\n",nch,x,y);
              // Fill the histograms
	      Nh1+=nhits;
	      //h->Fill(x,y,(float) 1);
		  //}
              //}
	   }          
	   
/*	   Int_t ncerenkovs = Cerenkovs->GetEntriesFast();

	   if (ncerenkovs) {
	       for (Int_t hit=0;hit<ncerenkovs;hit++) {
		   cHit = (AliRICHCerenkov*) Cerenkovs->UncheckedAt(hit);
		   Int_t nchamber = cHit->fChamber;     // chamber number
		   Int_t index =    cHit->fTrack;
		   Int_t pindex =   cHit->fIndex;
		   Int_t cx  =      cHit->fX;                // x-position
		   Int_t cy  =      cHit->fZ;                // y-position
		   Int_t cmother =  cHit->fCMother;      // Index of mother particle
		   Int_t closs =    cHit->fLoss;           // How did the paryicle get lost? 
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
		       cerenkov->Fill(cx,cy,(float) 1);
		     
		     TParticle *MIP = (TParticle*)(*gAlice->Particles())[current->GetFirstMother()];
		     //TParticle *MIP = (TParticle*)(*gAlice->Particles())[MIP1->GetFirstDaughter()];
		     //printf("Second Mother:%d",MIP1->GetFirstDaughter());
		     mom[0] = current->Px();
		     mom[1] = current->Py();
		     mom[2] = current->Pz();
		     Float_t energy = current->Energy();
		     Float_t Mip_px = MIP->Px();
		     Float_t Mip_py = MIP->Py();
		     Float_t Mip_pz = MIP->Pz();
		     
		     Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
		     Float_t rt = TMath::Sqrt(r);
		     Float_t Mip_r = Mip_px*Mip_px + Mip_py*Mip_py + Mip_pz*Mip_pz;	
		     Float_t Mip_rt = TMath::Sqrt(Mip_r);
		     Float_t coscerenkov = (mom[0]*Mip_px + mom[1]*Mip_py + mom[2]*Mip_pz)/(rt*Mip_rt);
		     Float_t cherenkov = TMath::ACos(coscerenkov);
		     ckovangle->Fill(cherenkov,(float) 1);                           //Cerenkov angle calculus
		     Float_t ckphi=TMath::ATan2(mom[0], mom[2]);
		     hckphi->Fill(ckphi,(float) 1);
		     
		     //mipHit = (AliRICHHit*) Hits->UncheckedAt(0);
		     
		     Float_t mx = MIP->Vx();
		     Float_t my = MIP->Vz();
		     Float_t mz = MIP->Vy();
		     
		     //Float_t mx = mipHit->fX;
		     //Float_t my = mipHit->fZ;
		     Float_t dx = cx - mx;
		     Float_t dy = cy - my;
		     Float_t final_radius = TMath::Sqrt(dx*dx+dy*dy);
		     radius->Fill(final_radius,(float) 1);
		     
		     if (closs == 4)
		       {
			 phspectra1->Fill(energy*1e9,(float) 1);
			 phot++;
		       }
		     else
		       phspectra2->Fill(energy*1e9,(float) 1);
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
	   }*/
       }

      /* for (Int_t nmothers=0;nmothers<ntracks;nmothers++){
	   totalphotons->Fill(mothers[nmothers],(float) 1);
	   mother->Fill(mothers2[nmothers],(float) 1);
	   //printf ("Entries in %d : %d\n",nmothers, mothers[nmothers]);
       }*/
       
      /* clusev->Fill(nraw,(float) 1);
       photev->Fill(phot,(float) 1);
       feedev->Fill(feed,(float) 1);
       padsmip->Fill(padmip,(float) 1);
       padscl->Fill(pads,(float) 1);
       printf("Photons:%d\n",phot);
       phot = 0;
       feed = 0;
       pads = 0;
       nraw=0;
       padmip=0;*/

   }
   
   //Create canvases, set the view range, show histograms

   switch(diaglevel)
     {
     case 1:

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
       hitsX->SetXTitle("(GeV)");
       hitsX->Draw();
       
       c10->cd(2);
       hitsY->SetFillColor(42);
       hitsY->SetXTitle("(GeV)");
       hitsY->Draw();
       
      
       break;
//
     case 2:
       
       /*TCanvas *c6 = new TCanvas("c6","Photon Spectra",50,10,600,350);
       c6->Divide(2,1);
       
       c6->cd(1);
       phspectra2->SetFillColor(42);
       phspectra2->SetXTitle("energy (eV)");
       phspectra2->Draw();
       c6->cd(2);
       phspectra1->SetFillColor(42);
       phspectra1->SetXTitle("energy (eV)");
       phspectra1->Draw();*/
       
       //TCanvas *c9 = new TCanvas("c9","Particles Spectra",400,10,600,700);
       TCanvas *c9 = new TCanvas("c9","Pion Spectra",400,10,600,700);
       //c9->Divide(2,2);
       
       //c9->cd(1);
       pionspectra1->SetFillColor(42);
       pionspectra1->SetXTitle("log(GeV)");
       pionspectra2->SetFillColor(46);
       pionspectra2->SetXTitle("log(GeV)");
       pionspectra3->SetFillColor(10);
       pionspectra3->SetXTitle("log(GeV)");
       //c9->SetLogx();
       pionspectra1->Draw();
       pionspectra2->Draw("same");
       pionspectra3->Draw("same");
       
       //c9->cd(2);
       
       TCanvas *c10 = new TCanvas("c10","Proton Spectra",400,10,600,700);
       protonspectra1->SetFillColor(42);
       protonspectra1->SetXTitle("log(GeV)");
       protonspectra2->SetFillColor(46);
       protonspectra2->SetXTitle("log(GeV)");
       protonspectra3->SetFillColor(10);
       protonspectra3->SetXTitle("log(GeV)");
       //c10->SetLogx();
       protonspectra1->Draw();
       protonspectra2->Draw("same");
       protonspectra3->Draw("same");

       //c9->cd(3);
      TCanvas *c11 = new TCanvas("c11","Kaon Spectra",400,10,600,700); 
       kaonspectra1->SetFillColor(42);
       kaonspectra1->SetXTitle("log(GeV)");
       kaonspectra2->SetFillColor(46);
       kaonspectra2->SetXTitle("log(GeV)");
       kaonspectra3->SetFillColor(10);
       kaonspectra3->SetXTitle("log(GeV)");
       //c11->SetLogx();
       kaonspectra1->Draw();
       kaonspectra2->Draw("same");
       kaonspectra3->Draw("same");
       
       //c9->cd(4);
       TCanvas *c12 = new TCanvas("c12","Charged Particles Spectra",400,10,600,700);
       chargedspectra1->SetFillColor(42);
       chargedspectra1->SetXTitle("log(GeV)");
       chargedspectra2->SetFillColor(46);
       chargedspectra2->SetXTitle("log(GeV)");
       chargedspectra3->SetFillColor(10);
       chargedspectra3->SetXTitle("log(GeV)");
       //c12->SetLogx();
       chargedspectra1->Draw();
       chargedspectra2->Draw("same");
       chargedspectra3->Draw("same");

       //TCanvas *c16 = new TCanvas("c16","Particles Spectra II",400,10,600,700);
       //c16->Divide(2,2);
       
       //c16->cd(1);
       TCanvas *c13 = new TCanvas("c13","Electron Spectra",400,10,600,700);
       electronspectra1->SetFillColor(42);
       electronspectra1->SetXTitle("log(GeV)");
       electronspectra2->SetFillColor(46);
       electronspectra2->SetXTitle("log(GeV)");
       electronspectra3->SetFillColor(10);
       electronspectra3->SetXTitle("log(GeV)");
       //c13->SetLogx();
       electronspectra1->Draw();
       electronspectra2->Draw("same");
       electronspectra3->Draw("same");
       
       //c16->cd(2);
       TCanvas *c14 = new TCanvas("c14","Muon Spectra",400,10,600,700);
       muonspectra1->SetFillColor(42);
       muonspectra1->SetXTitle("log(GeV)");
       muonspectra2->SetFillColor(46);
       muonspectra2->SetXTitle("log(GeV)");
       muonspectra3->SetFillColor(10);
       muonspectra3->SetXTitle("log(GeV)");
       //c14->SetLogx();
       muonspectra1->Draw();
       muonspectra2->Draw("same");
       muonspectra3->Draw("same");

       //c16->cd(4);
       TCanvas *c16 = new TCanvas("c16","Neutron Spectra",400,10,600,700);
       neutronspectra1->SetFillColor(42);
       neutronspectra1->SetXTitle("log(GeV)");
       neutronspectra2->SetFillColor(46);
       neutronspectra2->SetXTitle("log(GeV)");
       neutronspectra3->SetFillColor(10);
       neutronspectra3->SetXTitle("log(GeV)");
       //c16->SetLogx();
       neutronspectra1->Draw();
       neutronspectra2->Draw("same");
       neutronspectra3->Draw("same");

       TCanvas *c15 = new TCanvas("c15","Mothers Production Vertices",500,100,800,800);
       production->SetFillColor(42);
       production->SetXTitle("z (m)");
       production->SetYTitle("R (m)");
       production->Draw();

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
       }

       if (nev<1)
	 {
	   TCanvas *c11 = new TCanvas("c11","Cherenkov per Mip",400,10,600,700);
	   mother->SetFillColor(42);
	   mother->SetXTitle("counts");
	   mother->Draw();
	 }
       
       
       TCanvas *c5 = new TCanvas("c5","Ring Statistics",50,10,600,350);
       c5->Divide(2,1);
       
       c5->cd(1);
       ckovangle->SetFillColor(42);
       ckovangle->SetXTitle("angle (radians)");
       ckovangle->Draw();
       
       c5->cd(2);
       radius->SetFillColor(42);
       radius->SetXTitle("radius (cm)");
       radius->Draw();

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
      
   /*
	 TCanvas *c8 = new TCanvas("c25","Number of pads per event inside MIP region",400,10,600,700);
	 padsmip->SetFillColor(42);
	 padsmip->SetXTitle("(counts)");
	 padsmip->Draw(); 
       */

       
       /*TCanvas *c8 = new TCanvas("c8","Number of pads per event inside MIP region",400,10,600,700);
	 hckphi->SetFillColor(42);
	 hckphi->SetXTitle("phi");
	 hckphi->Draw();*/ 
       

   // calculate the number of pads which give a signal


   Int_t Np=0;
   for (Int_t i=0; i< NpadX;i++) {
       for (Int_t j=0;j< NpadY;j++) {
	   if (Pad[i][j]>=6){
	       Np+=1;
	   }
       }
   }
   //printf("The total number of pads which give a signal: %d %d\n",Nh,Nh1);
   
   printf("****************************************\n");
   printf("* Particle                  * Flux(m2) *\n");
   printf("****************************************\n");

   printf("* Pions:                    *   %3.1f   *\n",pion/11.757312);
   printf("* Charged Pions:            *   %3.1f   *\n",chargedpions/11.757312);
   printf("* Primary Pions:            *   %3.1f   *\n",primarypions/11.757312);
   printf("* Primary Pions (p>1GeV):   *   %3.1f   *\n",highprimarypions/11.757312);
   printf("* Kaons:                    *   %3.1f   *\n",kaon/11.757312);
   printf("* Charged Kaons:            *   %3.1f   *\n",chargedkaons/11.757312);
   printf("* Primary Kaons:            *   %3.1f   *\n",primarykaons/11.757312);
   printf("* Primary Kaons (p>1GeV):   *   %3.1f   *\n",highprimarykaons/11.757312);
   printf("* Muons:                    *   %3.1f   *\n",muon/11.757312);
   printf("* Electrons:                *   %3.1f   *\n",electron/11.757312);
   printf("* Positrons:                *   %3.1f   *\n",positron/11.757312);
   printf("* Protons:                  *   %3.1f   *\n",proton/11.757312);
   printf("* All Charged:              *   %3.1f   *\n",(chargedpions+chargedkaons+muon+electron+positron+proton)/11.757312);
   printf("****************************************\n");
   printf("* Photons:                  *   %3.1f   *\n",photons/11.757312); 
   printf("* Primary Photons:          *   %3.1f   *\n",primaryphotons/11.757312);
   printf("* Primary Photons (p>1MeV): *   %3.1f   *\n",highprimaryphotons/11.757312);
   printf("****************************************\n");
   printf("* Neutrons:                 *   %3.1f   *\n",neutron);
   printf("* Neutrons (p>100keV:       *   %3.1f   *\n",highneutrons);
   printf("****************************************\n");

   printf("End of macro\n");
}





