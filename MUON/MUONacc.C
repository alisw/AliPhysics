Bool_t SelectThis(AliMUONHit* mHit);

Bool_t SelectThis(Float_t x, Float_t y, Float_t d)
{
//    Float_t x     =    mHit->X();        // x-pos 
//    Float_t y     =    mHit->Y();        // y-pos
    if (TMath::Abs(x) > d/2. && TMath::Abs(y) > d/2.) {
	return 1;
    } else {
	return 0;
    }
}

Bool_t SelectThat(Float_t x, Float_t y, Float_t d)
{
    const Float_t phi0=60.*TMath::Pi()/180.;
    
    if (y<0.) return 1;
    if (x<0.) x=-x;
    
    Float_t phi = TMath::ATan2(y,x);
    Float_t dist = TMath::Sin(phi0-phi)*TMath::Sqrt(x*x+y*y);
    if (TMath::Abs(dist) > d) {
	return 1;
    } else {
	return 0 ;
    }
}

void MUONacc (Float_t d=3., Int_t evNumber1=0, Int_t evNumber2=0) 
{
// Dynamically link some shared libs                    
    Float_t rmin[14] = {17.5, 17.5, 23.5, 23.5, 33.5, 33.5, 43., 43., 50., 50,
			56.1, 56.1, 59.6, 59.6};
    Float_t rmax[14] = {91.5, 91.5, 122.5, 122.5, 158.3, 158.3, 260.0, 
			260.0, 260.0, 260.0, 293., 293, 293., 293.};
    

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

    TH1F *theta   =  new TH1F("theta","Theta distribution",180,0,180);
    TH1F *emult   =  new TH1F("emult","Event Multiplicity",100,0,1000);   
    TH2F *supp    =  new TH2F("supp", "Diselected",100, -200., 200, 100, -200., 200.);   

    AliMUONChamber*  iChamber;
    AliSegmentation* seg;
    
//
//   Loop over events 
//
    
    for (Int_t nev=0; nev<= evNumber2; nev++) {
	Int_t nparticles = gAlice->GetEvent(nev);
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;
	
	AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
	
	
	TTree *TH = gAlice->TreeH();
	Int_t ntracks = TH->GetEntries();
	TClonesArray *fPartArray = gAlice->Particles();       
	Int_t npart = fPartArray->GetEntriesFast();
	
	Int_t EvMult=0;
	Int_t Nsel=0;
//
// Loop over primary particles (jpsi. upsilon, ...)
//       
	Int_t prim=0;
	
	for (Int_t part=0; part<30000; part+=3) {
	    TParticle *MPart = gAlice->Particle(part);
	    Int_t mpart  = MPart->GetPdgCode();
	    Int_t child1 = MPart->GetFirstDaughter();
	    Int_t child2 = MPart->GetLastDaughter();	
	    Int_t mother = MPart->GetFirstMother();	    
//
// Loop over muons
//
	    Int_t muon1=2*prim;
	    Int_t muon2=muon1+1;

//	    printf("\n PDG Code %d %d %d %d %d %d %d %d \n", 
//		   part, mpart, child1, child2, mother,prim, muon1,muon2);


	    Bool_t selected[2]={1,1};
	    Int_t muons=0;
	    for (Int_t track=muon1; track<=muon2;track++) {
		gAlice->ResetHits();
		Int_t nbytes += TH->GetEvent(track);
		
		if (MUON)  {
//
//  Loop over hits
//
		    for(AliMUONHit* mHit=(AliMUONHit*)MUON->FirstHit(-1); 
			mHit;
			mHit=(AliMUONHit*)MUON->NextHit()) 
		    {
			Int_t   nch   =    mHit->Chamber();  // chamber number
			Float_t x     =    mHit->X();        // x-pos of hit
			Float_t y     =    mHit->Y();        // y-pos
			Float_t z     =    mHit->Z();        // z-pos
			Float_t Eloss =    mHit->Eloss();    // energy loss
			Float_t Theta =    mHit->Theta();    // theta
			Float_t Particle = mHit->Particle(); // Particle type
			Int_t itrack   = Int_t(mHit->Track());

		        TParticle *thePart = gAlice->Particle(itrack);
			Float_t pTheta=thePart->Theta();
			
			if (Particle != kMuonPlus && Particle != kMuonMinus) continue;
			
			Float_t P   = mhit->Momentum();
			Float_t R   = TMath::Sqrt(x*x+y*y);
			TParticlePDG* Part = DataBase->GetParticle(Particle);
			Double_t mass = Part->Mass();
			theta->Fill(pTheta*180./TMath::Pi(),1.);
			
			if (nch>14) continue;
//			printf("\n %f %f %f %d\n ", R, rmin[nch-1], rmax[nch-1], nch);
// Geometrical acceptance of frames			
//			selected[muons] = selected[muons]&&SelectThis(x,y,d);
// Cut on default rmin and rmax
//			selected[muons] = selected[muons] &&( R > rmin[nch-1] && R < rmax[nch-1]);
// Cut on rmin
//			selected[muons] = selected[muons] &&( R > rmin[nch-1]);
// Cut on rmax
//			selected[muons] = selected[muons] &&( R < rmax[nch-1]);
//			selected[muons] = selected[muons] &&( R > z*TMath::Tan(2.0*TMath::Pi()/180.));
//			Bool_t ok = (z<970 && ( R > z*TMath::Tan(2.0*TMath::Pi()/180.))) ||
//			    (z>970 && ( R > 970 * TMath::Tan(2.0*TMath::Pi()/180.)
//					+(z-970)* TMath::Tan(1.6*TMath::Pi()/180.)));
//
//			Bool_t ok = (z<970 && ( R < z*TMath::Tan(9.0*TMath::Pi()/180.))) ||
//			    (z>970 && ( R < 970 * TMath::Tan(9.0*TMath::Pi()/180.)
//					+(z-970)* TMath::Tan(12.0*TMath::Pi()/180.)));
			if (nch == 4 ||  nch == 7) {
			    Bool_t ok = SelectThat(x,y,d);
			} else {
			    Bool_t ok = kTRUE;
			}
			
			
			    selected[muons] = selected[muons] && ok;
			if (!ok) {
			    supp->Fill(x,y,1.);
			    
			}
			

		    } // hit loop
		} // if MUON
		muons++;
	    } // track loop
	    if (selected[0] && selected[1]) Nsel++;
	    prim++;
	} // primary loop
	printf("\n Selected %d %d \n", Nsel, prim);
    }
    
//Create a canvas, set the view range, show histograms
    Int_t k;
    TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
    c1->Divide(2,2);
    c1->cd(1);
    theta->Draw();
    c1->cd(2);
    supp->Draw();
    
}









