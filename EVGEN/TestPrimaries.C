void TestPrimaries(Int_t evNumber1=0, Int_t evNumber2=0) 
{
// Dynamically link some shared libs                    
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

    TH1F *thetaH   =  new TH1F("thetaH","Theta distribution",180,0,180);
    TH1F *phiH     =  new TH1F("phiH","Phi distribution"  ,180,-180,180);
    TH1F *etaH     =  new TH1F("etaH","Pseudorapidity",120,-12,12);
    TH1F *yH       =  new TH1F("yH","Rapidity distribution",120,-12,12);
    TH1F *eH       =  new TH1F("eH","Energy distribution",100,0,100);
    TH1F *eetaH    =  new TH1F("eetaH","Pseudorapidity",120,-12,12);
    TH1F *ptH      =  new TH1F("ptH","Pt distribution",150,0,15);
//
//   Loop over events 
//
    
    for (Int_t nev=0; nev<= evNumber2; nev++) {
	Int_t nparticles = gAlice->GetEvent(nev);
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;
	TClonesArray *fPartArray = gAlice->Particles();       
	Int_t npart = fPartArray->GetEntriesFast();
//
// Loop over primary particles (jpsi. upsilon, ...)
//       
	for (Int_t part=0; part<npart; part++) {
	    TParticle *MPart = (TParticle*) fPartArray->UncheckedAt(part);
	    Int_t mpart  = MPart->GetPdgCode();
	    Int_t child1 = MPart->GetFirstDaughter();
	    Int_t child2 = MPart->GetLastDaughter();	
	    Int_t mother = MPart->GetFirstMother();	   
	    
	    Float_t Pt = MPart->Pt();
	    Float_t E  = MPart->Energy();
	    Float_t Pz = MPart->Pz();
	    Float_t Py = MPart->Py();
	    Float_t Px = MPart->Px();
	    Float_t pT = TMath::Sqrt(Px*Px+Py*Py);
	    Float_t theta = MPart->Theta();
	    Float_t phi   = MPart->Phi()-TMath::Pi();
	    Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));
	    Float_t y     = TMath::Log((E+Pz)/(E-Pz+1.e-13));

	    thetaH->Fill(theta*180./TMath::Pi(),1.);
	    phiH->Fill(phi*180./TMath::Pi(),1.);
	    etaH->Fill(eta,1.);    
	    eetaH->Fill(eta,E);    
	    yH->Fill(y,1.);      
	    eH->Fill(E,1.);      
	    ptH->Fill(pT,1.);     
	} // primary loop
    }
    
//Create a canvas, set the view range, show histograms
    TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
    c1->Divide(2,2);
    c1->cd(1); ptH->Draw();
    c1->cd(2); etaH->Draw();
    c1->cd(3); yH->Draw();
    c1->cd(4



); eH->Draw();

    TCanvas *c2 = new TCanvas("c2","Canvas 1",400,10,600,700);
    c2->Divide(2,2);
    c2->cd(1); phiH->Draw();
    c2->cd(2); thetaH->Draw();
    c2->cd(3); eetaH->Draw();
}







