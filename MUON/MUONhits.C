void MUONhits (Int_t evNumber1=0,Int_t evNumber2=0, Int_t ic=1) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs                    
    Float_t rmin[14] = {17.5, 17.5, 23.5, 23.5, 33.5, 33.5, 43., 43., 50., 50,
			56.1, 56.1, 59.6, 59.6};
    Float_t rmax[14] = {81.6, 81.6, 109.3, 109.3, 154.4, 154.4, 197.8, 
			197.8, 229.5, 229.5, 254., 254, 270., 270.};
    

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
//
    TDatabasePDG*  DataBase = new TDatabasePDG();
    
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

    TH1F *Hits[14],   *HitDensity[14];
    TH1F *Hits_g[14], *HitDensity_g[14];
    TH1F *Hits_n[14], *HitDensity_n[14];
    TH1F *Mom[14], *Mom_g[14], *Mom_n[14];
    
    for (Int_t i=0; i<14; i++) {
	Hits[i]     =  new TH1F("hits1","Hit Distribution",30,0,300);
	Hits[i]->SetFillColor(0);
	Hits[i]->SetXTitle("R (cm)");
	
	HitDensity[i]  =  new TH1F("dhits1","Hit Density",30,0,300);
	HitDensity[i]->SetFillColor(0);
	HitDensity[i]->SetXTitle("R (cm)");

	Hits_g[i]     =  new TH1F("hits1","Hit Distribution",30,0,300);
	Hits_g[i]->SetFillColor(0);
	Hits_g[i]->SetXTitle("R (cm)");

	HitDensity_g[i]  =  new TH1F("dhits1","Hit Density",30,0,300);
	HitDensity_g[i]->SetFillColor(0);
	HitDensity_g[i]->SetXTitle("R (cm)");

	Mom[i]  =  new TH1F("mom","Energy",70,-6,1);
	Mom[i]->SetFillColor(0);
	Mom[i]->SetXTitle("E (GeV)");

	Mom_g[i]  =  new TH1F("mom","Energy",70,-6,1);
	Mom_g[i]->SetFillColor(0);
	Mom_g[i]->SetXTitle("E (GeV)");

	Mom_n[i]  =  new TH1F("mom","Energy",70,-6,1);
	Mom_n[i]->SetFillColor(0);
	Mom_n[i]->SetXTitle("E (GeV)");

	Hits_n[i]     =  new TH1F("hits1","Hit Distribution",30,0,300);
	Hits_n[i]->SetFillColor(0);
	Hits_n[i]->SetXTitle("R (cm)");

	HitDensity_n[i]  =  new TH1F("dhits1","Hit Density",30,0,300);
	HitDensity_n[i]->SetFillColor(0);
	HitDensity_n[i]->SetXTitle("R (cm)");
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
		    Int_t   nch   =    mHit->fChamber;  // chamber number
		    Float_t x     =    mHit->fX;        // x-pos of hit
		    Float_t y     =    mHit->fY;        // y-pos
		    Float_t Eloss =    mHit->fEloss;
		    Float_t Theta =    mHit->fTheta;
		    Float_t Particle = mHit->fParticle;
		    Float_t P    =     
			TMath::Sqrt(mHit->fCxHit*mHit->fCxHit+
				    mHit->fCyHit*mHit->fCyHit+
				    mHit->fCzHit*mHit->fCzHit);
		    TParticlePDG* Part = DataBase->GetParticle(Particle);
		    Double_t mass = Part->Mass();
		    
		    if (nch >13) continue;
		    if (nch ==1) EvMult++;
		    if (mHit->fAge > 5.e-6) continue;
		   
		    Float_t r=TMath::Sqrt(x*x+y*y);
		    if (nch ==0) continue;
		    
		    if (r < rmin[nch-1]) continue;
		    if (r > rmax[nch-1]) continue;
		   
		    Float_t wgt=1/(2*10*TMath::Pi()*r)/(evNumber2+1);

		    Float_t Ekin=TMath::Sqrt(P*P+mass*mass)-mass;
		    if (Particle == 2112) {
			Hits_n[nch-1]->Fill(r,(float) 1);
			HitDensity_n[nch-1]->Fill(r,wgt);
			Mom_n[nch-1]->Fill(TMath::Log10(Ekin),1);
		    } else if (Particle == 22) {
			Hits_g[nch-1]->Fill(r,(float) 1);
			HitDensity_g[nch-1]->Fill(r,wgt);
			Mom_g[nch-1]->Fill(TMath::Log10(Ekin),1);
		    } else {
			Hits[nch-1]->Fill(r,(float) 1);
			HitDensity[nch-1]->Fill(r,wgt);
			Mom[nch-1]->Fill(TMath::Log10(Ekin),1);
		    }
		} // hit loop
	    } // if MUON
	} // track loop
    }
    
//Create a canvas, set the view range, show histograms
    Int_t k;
    TCanvas *c1 = new TCanvas("c1","Hit Densities",400,10,600,700);
    c1->Divide(2,4);
    for (k=0; k<7; k++) {
	c1->cd(k+1);
	HitDensity[2*k]->Draw();
    }

    TCanvas *c2 = new TCanvas("c2","Hit Densities (gamma)",400,10,600,700);
    c2->Divide(2,4);
    for (k=0; k<7; k++) {
	c2->cd(k+1);
	HitDensity_g[2*k]->Draw();
    }
    
    TCanvas *c3 = new TCanvas("c3","Hit Densities (neutron)",400,10,600,700);
    c3->Divide(2,4);
    for (k=0; k<7; k++) {
	c3->cd(k+1);
	HitDensity_n[2*k]->Draw();
    }

    TCanvas *c4 = new TCanvas("c4","Energy (charged)",400,10,600,700);
    c4->Divide(2,4);
    for (k=0; k<7; k++) {
	c4->cd(k+1);
	Mom[2*k]->Draw();
    }

    TCanvas *c5 = new TCanvas("c5","Energy (gamma)",400,10,600,700);
    c5->Divide(2,4);
    for (k=0; k<7; k++) {
	c5->cd(k+1);
	Mom_g[2*k]->Draw();
    }

    TCanvas *c6 = new TCanvas("c6","Energy (gamma)",400,10,600,700);
    c6->Divide(2,4);
    for (k=0; k<7; k++) {
	c6->cd(k+1);
	Mom_n[2*k]->Draw();
    }
}







