   TH1F *mass1 = new TH1F("mass1","Invariant Mass",120,0,12);
   TH1F *mass2 = new TH1F("mass2","Invariant Mass",100,9.0,10.5);
   TH1F *mass3 = new TH1F("mass3","Invariant Mass",100,2.5,3.5);
void MUONstraggling (Int_t evNumber1=0,Int_t evNumber2=0) 
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


   AliMUONChamber*  iChamber;
//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (Int_t nev=0; nev<= evNumber2; nev++) {
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

       Float_t pups[4]={0,0,0,0};
       for (Int_t track=0; track<ntracks;track++) {
	   Int_t Nelt=0;	   
	   gAlice->ResetHits();
	   Int_t nbytes += TH->GetEvent(track);

	   
	   if (MUON)  {
	       for(AliMUONHit* mHit=(AliMUONHit*)MUON->FirstHit(-1); 
		   mHit;
		   mHit=(AliMUONHit*)MUON->NextHit()) 
	       {
		   Int_t   nch   = mHit->Chamber();  // chamber number
		   Float_t x     = mHit->X();        // x-pos of hit
		   Float_t y     = mHit->Y();        // y-pos
		   Float_t z     = mHit->Z();        // y-pos
		   Float_t p  = mHit->Momentum();
		   Float_t px = mHit->Px();
		   Float_t py = mHit->Py();
		   Float_t pz = mHit->Pz();
		   
		   if (nch != 1) continue;
		   
		   Int_t ipart = mHit->Particle();
		   TParticle *Part;
		   Int_t ftrack = mHit->Track();
		   Part = gAlice->Particle(ftrack);
		   Int_t ipart = Part->GetPdgCode();
		   TParticle *Mother;
		   Float_t px0 = Part->Px();
		   Float_t py0 = Part->Py();		   
		   Float_t pz0 = Part->Pz();
		   Float_t p0  = Part->P();
		   
		   Float_t e0=Part->Energy();
		   
		   if (ipart == kMuonPlus || ipart == kMuonMinus) {
//
// Branson Correction
//
		       Float_t zch = z;
		       Float_t r   = TMath::Sqrt(x*x+y*y);
		       Float_t zb;
		       
		       if (r<26.3611) {
			   zb = 466.;
		       } else {
			   zb = 441.;
		       }

		       Float_t xb  = x-(zch-zb)*px/pz;
		       Float_t yb  = y-(zch-zb)*py/pz;		       
		       Float_t pzB = p * zb/TMath::Sqrt(zb*zb+yb*yb+xb*xb);
		       Float_t pxB = pzB * xb/zb;
		       Float_t pyB = pzB * yb/zb;
		       Float_t theta 
			   = TMath::ATan2(TMath::Sqrt(pxB*pxB+pyB*pyB), pzB) *
			   180./TMath::Pi();
		       
		       
//
// Energy Correction
//
//
		       
		       Float_t corr=(p+CorrectP(p, theta))/p;
		       pups[0]+=p*corr;
		       pups[1]+=pxB*corr;
		       pups[2]+=pyB*corr;
		       pups[3]+=pzB*corr;
 		   }
	       } // hits 
	   } // if MUON 
       } // tracks
       Float_t mass =
	   TMath::Sqrt(pups[0]*pups[0]
		       -pups[1]*pups[1]-pups[2]*pups[2]-pups[3]*pups[3]);
       mass1->Fill(mass, 1.);
       mass2->Fill(mass, 1.);
       mass3->Fill(mass, 1.);	       
    } // event
//Create a canvas, set the view range, show histograms
   TCanvas *c1 = new TCanvas("c1","Mass Spectra",400,10,600,700);

   mass2->SetFillColor(42);
   mass2->SetXTitle("Mass (GeV)");
   mass2->Draw();
}



Float_t CorrectP(Float_t p, Float_t theta)
{
    Float_t c;
    
    if (theta<3.) {
//W
	if (p<15) {
	    c = 2.916+0.0156*p-0.00000866*p*p;
	} else {
	    c = 2.998+0.0137*p;
	}
    } else {
//Cu+Fe+Cc+C
	if (p<15) {
	    c = 2.464+0.00877*p-0.0000106*p*p;
	} else {
	    c = 2.496+0.00817*p;
	}
    }
//
//  
    c*=0.95;
//
    return c;
}






