void testJetFinderBg(Int_t evNumber1=0, Int_t evNumber2=0) 
{
//*-- Author: Andreas Morsch (CERN)
// Dynamically link some shared libs                    
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("../macros/loadlibs.C");
	loadlibs();
    }
//
//  Create and configure JetFinder
//
    AliEMCALJetFinder* jetFinder = 
	new AliEMCALJetFinder("UA1 Jet Finder", "Test");
//
//  Debug
    jetFinder->SetDebug(1);
//
//  Input and fast simulation
    jetFinder->SetPtCut(2.);
    jetFinder->SetIncludeK0andN(0);
//    jetFinder->SetMomentumSmearing(1);
//    jetFinder->SetEfficiencySim(1);
    jetFinder->SetHadronCorrection(0);
//    jetFinder->SetHadronCorrector(AliEMCALHadronCorrectionv0::Instance());
    jetFinder->SetSamplingFraction(12.9);
//
//  Parameters for jet finding
    jetFinder->SetConeRadius(0.5);
    jetFinder->SetEtSeed(4.);
    jetFinder->SetMinJetEt(10.);
    jetFinder->SetMinCellEt(1.);

//.............
//            |
//            |
//            V
//  This part will go into class
//
//  Open background file
//
// 
    printf("\n Opening Background File !!\n");
    
    TFile* fileB =  new TFile("bg.root");
    gAlice = (AliRun*)(fileB->Get("gAlice"));
    Int_t nparticles = gAlice->GetEvent(0);
    // Read and save background event
    jetFinder->FillFromHits();
    jetFinder->FillFromTracks();
    jetFinder->SaveBackgroundEvent();
    delete gAlice;    


//
//   Open signal file
//
    printf("\n Opening Signal file !!\n");
    TFile* fileS =  new TFile("galice.root", "update");
    gAlice = (AliRun*)(fileS->Get("gAlice"));
    
//
//   Loop over signal events 
//

    Int_t nhit=0;
    for (Int_t nev = evNumber1; nev<= evNumber2; nev++) {
	Int_t nparticles = gAlice->GetEvent(nev);
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;
// Load background
	jetFinder->InitFromBackground();
// ECAL information	
	jetFinder->FillFromHits(1);
//                              ^ preserves info from background
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

//            ^
//            |
//            |
//            |
//.............
    fileB->Close();
    fileS->Close();
}
















