void readPrimaries(Int_t evNumber1=0, Int_t evNumber2=3500) 
{
// Dynamically link some shared libs
   /*                    
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
        //gSystem->Load("libpythia6_4_28");
    }
   */

        gSystem->Load("libpythia6_4_28");

// Connect the Root Galice file containing Geometry, Kine and Hits

    AliRunLoader* rl = AliRunLoader::Open("galice.root");
//
    TDatabasePDG*  DataBase = new TDatabasePDG();
    
//  Create some histograms

    TH1F *pte      =  new TH1F("pte","Pt distribution",50,0,50);
//
//   Loop over events 
//
    rl->LoadKinematics();
    rl->LoadHeader();    
    for (Int_t nev=0; nev< evNumber2; nev++) {
	rl->GetEvent(nev);
	AliStack* stack = rl->Stack();
	Int_t npart = stack->GetNprimary();
	if (nev < evNumber1) continue;
//
// Loop over primary particles (jpsi. upsilon, ...)
//       
	
	for (Int_t part=0; part<npart; part++) {
	    TParticle *MPart = stack->Particle(part);
	    Int_t mpart  = MPart->GetPdgCode();
            if(TMath::Abs(mpart)!=11)continue;

	    Float_t Pt = MPart->Pt();
	    Float_t E  = MPart->Energy();
	    Float_t Pz = MPart->Pz();
	    Float_t Py = MPart->Py();
	    Float_t Px = MPart->Px();
	    Float_t pT = TMath::Sqrt(Px*Px+Py*Py);
	    Float_t theta = MPart->Theta();
	    Float_t phi   = MPart->Phi()-TMath::Pi();
	    Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));
	    Float_t y     = 0.5*TMath::Log((E+Pz+1.e-13)/(E-Pz+1.e-13));
 
            Int_t iMom = MPart->GetFirstMother();
	    TParticle *MPart_mom = stack->Particle(iMom);
            //cout << MPart_mom->GetPdgCode() << endl;
            if(TMath::Abs(MPart_mom->GetPdgCode())!=24)continue; 

	    if(TMath::Abs(y)<0.6)pte->Fill(Pt);     

	} // primary loop
    }
    
//Create a canvas, set the view range, show histograms
    TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
    pte->Draw();

 TFile *fout = new TFile("We_pp.root","recreate");
 pte->Write("ze");

}


