void MUONtestabso (Int_t evNumber1=0,Int_t evNumber2=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs                    
    static Float_t xmuon, ymuon;
    
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
	    printf("\n Create new gAlice object");
	    gAlice = new AliRun("gAlice","Alice test program");
	}
    }

    
//  Create some histograms


   TH1F *zv = new TH1F("zv","z vertex" ,140,470.,505.);
   TH1F *xv = new TH1F("xv","x vertex" ,140,-70.,70.);
   TH1F *yv = new TH1F("yv","y vertex", 140,-70.,70.);

   TH1F *ip  = new TH1F("ip","geant part",50,0.,10.);
   TH2F *rzv  = new TH2F("rzv","R-z vert",100,502,504.,100,0.,100.);

   TH1F *ptUps  = new TH1F("ptUps","pT Upsilon",50,0.,25.);
   TH1F *hde  = new TH1F("hde","dE",200,0.,50);
   TH1F *hde2  = new TH1F("hde2","dE",100,1.5,5.5);
   TH2F *hdevsn  = new TH2F("hdevsn","dE vs N electron",100,0.,15., 20, 0.,20.);

   TH1F *ekine  = new TH1F("ekine","E_kin electrons",70,-5,2);
   TH1F *etheta = new TH1F("etheta","Theta electrons",90,0,90);
   TH1F *edr = new TH1F("edr","Distance to muon",100,0,10);

   TH1F *de  = new TH1F("de","correction",100,-1,1);

   TH1F *dtheta = new TH1F("dtheta","Delta Theta" ,200,-5.,5.);
   AliMUONChamber*  iChamber;
//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
       Int_t Nel1=0;
       Int_t Nel2=0;
       Int_t Nel3=0;
       Int_t Nel4=0;
       
   for (Int_t nev=0; nev<= evNumber2; nev++) {
       //cout << "nev         " << nev <<endl;
       Int_t nparticles = gAlice->GetEvent(nev);
       //cout << "nparticles  " << nparticles <<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       
       AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");


       TTree *TH = gAlice->TreeH();
       Int_t ntracks = TH->GetEntries();
//
//   Loop over tracks
//

       Float_t dE;
       
       for (Int_t track=0; track<ntracks;track++) {
	   Int_t Nelt=0;	   
	   printf("\n nev %d %d %d %d %d\n ", nev, Nel1, Nel2, Nel3, Nel4);
	   gAlice->ResetHits();
	   Int_t nbytes += TH->GetEvent(track);
	   if (MUON)  {
	       for(AliMUONHit* mHit=(AliMUONHit*)MUON->FirstHit(-1); 
		   mHit;
		   mHit=(AliMUONHit*)MUON->NextHit()) 
	       {
		   Int_t   nch   = mHit->fChamber;  // chamber number
		   Float_t x     = mHit->fX;        // x-pos of hit
		   Float_t y     = mHit->fY;        // y-pos
		   Float_t z     = mHit->fZ;        // y-pos
		   
		   if (nch != 1) continue;
		   
		   Int_t ipart = mHit->fParticle;
		   TClonesArray *fPartArray = gAlice->Particles();
		   TParticle *Part;
		   Int_t ftrack = mHit->fTrack;
		   Part = (TParticle*) fPartArray->UncheckedAt(ftrack);
		   Int_t ipart = Part->GetPdgCode();
		   ip->Fill((float)ipart);
		   TParticle *Mother;
		   
		   Float_t px0=Part->Px();
		   Float_t py0=Part->Py();		   
		   Float_t pz0=Part->Pz();
		   Float_t thetax0=TMath::ATan2(px0,pz0);
		   Float_t thetay0=TMath::ATan2(py0,pz0);		   
		   
		   if (ipart == kMuonPlus || ipart == kMuonMinus) {
//		       Int_t imo = Part->GetFirstMother();
//		       Mother = (TParticle*) fPartArray->UncheckedAt(imo);
//		       Float_t pt = Mother->Pt();
//		       ptUps->Fill(pt, (float) 1);
		       Float_t E=Part->Energy();
		       Float_t Eloc=mHit->fPTot;
		       Float_t corr=Eloc+CorrectP(Eloc,mHit->fTheta);
		       printf("\n %f %f %f", E, Eloc, corr);
		       de->Fill(E-corr,1.);
		       dE = E-Eloc;
		       if (dE<50)  hde->Fill(dE, (float) 1);
		       if (dE<5.5) hde2->Fill(dE, (float) 1);
		       xmuon=mHit->fX;
		       ymuon=mHit->fY;
		       Float_t thetax=TMath::ATan2(mHit->Px(), mHit->Momentum());
		       Float_t thetay=TMath::ATan2(mHit->Py(), mHit->Momentum());
		       dtheta->Fill((thetax-thetax0)*1000., 1.);
		       dtheta->Fill((thetay-thetay0)*1000., 1.);		       
		   }
		   
		   if (ipart == kElectron || ipart == kPositron) {

		   
		       Float_t xvert = Part->Vx();       //  vertex  
		       Float_t yvert = Part->Vy();       //  vertex
		       Float_t zvert = Part->Vz();       // z vertex 
		       if (zvert < 503 && mHit->fTheta<90) {
			   Nelt++;
			   Float_t px = Part->Px();
			   Float_t py = Part->Py();
			   Float_t pz = Part->Pz();
			   Float_t Ek = Part->Energy()-Part->GetMass();
			   
			   Int_t imo = Part->GetFirstMother();
			   Mother = (TParticle*) fPartArray->UncheckedAt(imo);
			   Int_t imot = Mother->GetPdgCode();
			   xv->Fill(xvert);
			   yv->Fill(yvert);
			   zv->Fill(zvert);
			   Float_t rvert=TMath::Sqrt(xvert*xvert+yvert*yvert);
			   rzv->Fill(zvert,rvert);
			   ekine->Fill(TMath::Log10(Ek),1.);
			   etheta->Fill(mHit->fTheta,1.);
			   Float_t ex=mHit->fX;
			   Float_t ey=mHit->fY;
			   dr=TMath::Sqrt((ex-xmuon)*(ex-xmuon)+(ey-ymuon)*(ey-ymuon));
			   edr->Fill(dr,1.);
			   
		       }
		   }
		   
	       } // hits 
	   } // if MUON 
	   if (Nelt == 1) Nel1++;
	   if (Nelt == 2) Nel2++;
	   if (Nelt == 3) Nel3++;
	   if (Nelt  > 3) Nel4++;
	   hdevsn->Fill(dE, (float) Nelt, (float) 1);

       } // tracks

    } // event

   
//Create a canvas, set the view range, show histograms
   TCanvas *c1 = new TCanvas("c1","Vetices from electrons and positrons",400,10,600,700);
   c1->Divide(2,2);
   
   c1->cd(1);
   ip->SetFillColor(42);
   ip->SetXTitle("ipart");
   ip->Draw();

   c1->cd(2);
   xv->SetFillColor(42);
   xv->SetXTitle("xvert");
   xv->Draw();

   c1->cd(3);
   yv->SetFillColor(42);
   yv->SetXTitle("yvert");
   yv->Draw();

   c1->cd(4);
   zv->SetFillColor(42);
   zv->SetXTitle("zvert");
   zv->Draw();

   TCanvas *c2 = new TCanvas("c2"," ",400,10,600,700);
   c2->Divide(2,2);
   c2->cd(1);
   rzv->SetXTitle("zvert");
   rzv->SetYTitle("rvert");
   rzv->Draw();

   c2->cd(2);
   ptUps->SetXTitle("pt");
   ptUps->Draw();

   c2->cd(3);
   hde->SetXTitle("dE");
   hde->Draw();


   c2->cd(4);
   hde2->SetXTitle("dE");
   hde2->Draw();
   TCanvas *c3 = new TCanvas("c3"," ",400,10,600,700);
   c3->Divide(2,2);
   c3->cd(1);
   ekine->SetXTitle("E_kin");
   ekine->Draw();

   c3->cd(2);
   etheta->SetXTitle("Theta");
   etheta->Draw();

   c3->cd(3);
   edr->SetXTitle("Distance to muon");
   edr->Draw();

   c3->cd(4);
   dtheta->SetXTitle(" ");
   dtheta->Draw();

}

Float_t CorrectP(Float_t p, Float_t theta)
{
    printf("\n %f%", theta);
    
    if (theta<3.) {
//W
	if (p<15) {
	    return 2.737+0.0494*p-0.001123*p*p;
	} else {
	    return 3.0643+0.01346*p;
	}
    } else {
//Pb
	if (p<15) {
	    return 2.1380+0.0351*p-0.000853*p*p;
	} else {
	    return 2.407+0.00702*p;
	}
    }
}

