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
    else {
      //delete gAlice;
      gAlice = 0;
    }

    gAlice=0;
    
// Connect the Root Galice file containing Geometry, Kine and Hits
    
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (!file) file = new TFile("galice.root","UPDATE");
    
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
   TH2F *production = new TH2F("production","Mother production vertices",100,-300,300,100,0,600);
   
   
   

//   Start loop over events 

   Int_t Nh=0;
   Int_t pads=0;
   Int_t Nh1=0;
   Float_t mom[3];
   Int_t nraw=0;
   Int_t phot=0;
   Int_t feed=0;
   Int_t padmip=0;
   Int_t pion=0, kaon=0, proton=0, electron=0, positron=0, neutron=0, highneutrons=0, muon=0;
   Int_t chargedpions=0,primarypions=0,highprimarypions=0,chargedkaons=0,primarykaons=0,highprimarykaons=0;
   Int_t chargedmuons=0, photons=0, primaryphotons=0, highprimaryphotons=0;

   TRandom random;

   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       

       printf ("Event number       : %d\n",nev);
       printf ("Number of particles: %d\n",nparticles);
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       
// Get pointers to RICH detector and Hits containers
       
       AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
       Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
       gAlice->TreeR()->GetEvent(nent-1);
       TClonesArray *Rawclusters = RICH->RawClustAddress(2);    //  Raw clusters branch
       Int_t nrawclusters = Rawclusters->GetEntriesFast();
       TTree *TH = gAlice->TreeH(); 
       Int_t ntracks = TH->GetEntries();


       
// Start loop on tracks in the hits containers
       Int_t Nc=0;
       for (Int_t track=0; track<ntracks;track++) {
	   printf ("Processing Track: %d\n",track);
	   gAlice->ResetHits();
	   Int_t nbytes += TH->GetEvent(track);
	   if (RICH)  {
	       TClonesArray *PadHits = RICH->PadHits();      // Cluster branch
	       TClonesArray *Hits = RICH->Hits();            // Hits branch
	       TClonesArray *Cerenkovs = RICH->Cerenkovs();  // Cerenkovs branch
	   }
	   //see hits distribution
	   Int_t nhits = Hits->GetEntriesFast();
	   if (nhits) Nh+=nhits;
	   for (Int_t hit=0;hit<nhits;hit++) {
              mHit = (AliRICHHit*) Hits->UncheckedAt(hit);
              Int_t nch  = mHit->fChamber;              // chamber number
              Float_t x  = mHit->fX;                    // x-pos of hit
              Float_t y  = mHit->fZ;                    // y-pos
	      Float_t z  = mHit->fY;
	      Int_t index = mHit->fTrack;
	      Int_t particle = mHit->fParticle;    
	      Float_t R;

	      TParticle *current = (TParticle*)(*gAlice->Particles())[index];
	      
	      Float_t energy=current->Energy();

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
		      //printf("Adding %d at %f\n",particle,R);
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
		    //printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\R:%f\n\n\n\n\n\n\n\n\n",R);
		    }
		  //printf("Pion mass: %e\n",current->GetCalcMass());
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
		  //printf("Kaon mass: %e\n",current->GetCalcMass());
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
		  //printf("Electron mass: %e\n",current->GetCalcMass());
		  if (particle == 11)
		    electron +=1;
		  if (particle == -11)
		    positron +=1;
		}
	      if (TMath::Abs(particle)==13)
		{
		  muonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>.005 && current->Vy()>.005 && current->Vz()>.005)
		    muonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>2.5 && R<4.5)
		    muonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //printf("Muon mass: %e\n",current->GetCalcMass());
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
		      //printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\R:%f\n\n\n\n\n\n\n\n\n",R);
		    }
		  //printf("Neutron mass: %e\n",current->GetCalcMass());
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
	   
       }

   }
   
   //Create canvases, set the view range, show histograms

   TCanvas *c15 = new TCanvas("c15","Mothers Production Vertices",50,50,600,600);
   production->SetFillColor(42);
   production->SetXTitle("z (m)");
   production->SetYTitle("R (m)");
   production->Draw();

   TCanvas *c16 = new TCanvas("c16","Particles Spectra II",100,100,600,350);
   c16->Divide(2,1);
   
   c16->cd(1);
   //TCanvas *c13 = new TCanvas("c13","Electron Spectra",400,10,600,700);
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
   
   c16->cd(2);
   //TCanvas *c14 = new TCanvas("c14","Muon Spectra",400,10,600,700);
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
   
   //c16->cd(3);
   //TCanvas *c16 = new TCanvas("c16","Neutron Spectra",400,10,600,700);
   //neutronspectra1->SetFillColor(42);
   //neutronspectra1->SetXTitle("log(GeV)");
   //neutronspectra2->SetFillColor(46);
   //neutronspectra2->SetXTitle("log(GeV)");
   //neutronspectra3->SetFillColor(10);
   //neutronspectra3->SetXTitle("log(GeV)");
   //c16->SetLogx();
   //neutronspectra1->Draw();
   //neutronspectra2->Draw("same");
   //neutronspectra3->Draw("same");

   TCanvas *c9 = new TCanvas("c9","Particles Spectra",150,150,600,700);
   //TCanvas *c9 = new TCanvas("c9","Pion Spectra",400,10,600,700);
   c9->Divide(2,2);
   
   c9->cd(1);
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
   
   c9->cd(2);
   //TCanvas *c10 = new TCanvas("c10","Proton Spectra",400,10,600,700);
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
   
   c9->cd(3);
   //TCanvas *c11 = new TCanvas("c11","Kaon Spectra",400,10,600,700); 
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
   
   c9->cd(4);
   //TCanvas *c12 = new TCanvas("c12","Charged Particles Spectra",400,10,600,700);
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
   
   printf("*****************************************\n");
   printf("* Particle                   * Flux(m2) *\n");
   printf("*****************************************\n");

   printf("* Pions:                     *   %3.1f   *\n",pion/11.757312);
   printf("* Charged Pions:             *   %3.1f   *\n",chargedpions/11.757312);
   printf("* Primary Pions:             *   %3.1f   *\n",primarypions/11.757312);
   printf("* Primary Pions (p>1GeV/c):  *   %3.1f   *\n",highprimarypions/11.757312);
   printf("* Kaons:                     *   %3.1f   *\n",kaon/11.757312);
   printf("* Charged Kaons:             *   %3.1f   *\n",chargedkaons/11.757312);
   printf("* Primary Kaons:             *   %3.1f   *\n",primarykaons/11.757312);
   printf("* Primary Kaons (p>1GeV/c):  *   %3.1f   *\n",highprimarykaons/11.757312);
   printf("* Muons:                     *   %3.1f   *\n",muon/11.757312);
   printf("* Electrons:                 *   %3.1f   *\n",electron/11.757312);
   printf("* Positrons:                 *   %3.1f   *\n",positron/11.757312);
   printf("* Protons:                   *   %3.1f   *\n",proton/11.757312);
   printf("* All Charged:               *   %3.1f   *\n",(chargedpions+chargedkaons+muon+electron+positron+proton)/11.757312);
   printf("*****************************************\n");
   //printf("* Photons:                   *   %3.1f   *\n",photons/11.757312); 
   //printf("* Primary Photons:           *   %3.1f   *\n",primaryphotons/11.757312);
   //printf("* Primary Photons (p>1MeV/c):*   %3.1f   *\n",highprimaryphotons/11.757312);
   //printf("*****************************************\n");
   //printf("* Neutrons:                  *   %3.1f   *\n",neutron);
   //printf("* Neutrons (p>100keV/c):     *   %3.1f   *\n",highneutrons);
   //printf("*****************************************\n");

   printf("\nEnd of macro\n");
}





