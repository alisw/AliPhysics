void testJetFinder(Int_t evNumber1=0, Int_t evNumber2=0) 
{
//*-- Author: Andreas Morsch (CERN)
// Dynamically link some shared libs                    
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("../macros/loadlibs.C");
	loadlibs();
    }
// Connect the Root Galice file containing Geometry, Kine and Hits

    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");

    
    if (!file) {
	printf("\n Creating galice.root \n");
	file = new TFile("galice.root", "update");
    } else {
	printf("\n galice.root found in file list");
    }
    
// Get AliRun object from file or create it if not on file
    if (!gAlice) {
	gAlice = (AliRun*)(file->Get("gAlice"));
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) {
	    printf("\n create new gAlice object");
	    gAlice = new AliRun("gAlice","Alice test program");
	}
    }
//
//  Create and configure jet finder
    AliEMCALJetFinder* jetFinder = 
	new AliEMCALJetFinder("UA1 Jet Finder", "Test");
    jetFinder->SetDebug(1);
    jetFinder->SetConeRadius(0.5);
    jetFinder->SetEtSeed(4.);
    jetFinder->SetMinJetEt(10.);
    jetFinder->SetMinCellEt(1.);
    jetFinder->SetPtCut(2.);
//    jetFinder->SetMomentumSmearing(1);
//    jetFinder->SetEfficiencySim(1);
    jetFinder->SetHadronCorrection(0);
//    jetFinder->SetHadronCorrector(AliEMCALHadronCorrectionv0::Instance());
    jetFinder->SetSamplingFraction(12.9);
//
//  I/O
    jetFinder->SetOutputFileName("jets.root");
//
//  Initialization    
    jetFinder->Init();
//
//   Loop over events 
//

    Int_t nhit=0;
    for (Int_t nev = evNumber1; nev<= evNumber2; nev++) {
	file->cd();
	Int_t nparticles = gAlice->GetEvent(nev);
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;
// ECAL information	
	jetFinder->FillFromHits();
//	jetFinder->FillFromDigits();
// TPC  information
	jetFinder->FillFromTracks(1, 0);
//                                ^ preserves info from hit

// TPC  information from Hits associated to the EMCAL
//      jetFinder->FillFromHitFlaggedTracks(1);
//  
	jetFinder->Find();

//
// Look at results
	printf("\n Test Jets: %d\n", jetFinder->Njets());
	Int_t njet = jetFinder->Njets();
	for (Int_t nj=0; nj<njet; nj++)
	{
	    printf("\n Jet Energy %5d %8.2f \n", 
		   nj, jetFinder->JetEnergy(nj));
	    printf("\n Jet Phi    %5d %8.2f %8.2f \n", 
		   nj, jetFinder->JetPhiL(nj), jetFinder->JetPhiW(nj));
	    printf("\n Jet Eta    %5d %8.2f %8.2f \n", 
		   nj, jetFinder->JetEtaL(nj), jetFinder->JetEtaW(nj));
	}
//      TCanvas *c1 = new TCanvas("c1","Canvas 1",400,10,600,700);
//      (jetFinder->GetLego())->Draw();
    } // event loop
}
















