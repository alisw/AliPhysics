#include "iostream.h"
// example of macro which retrieves the Trigger Output from galice.root 

void MUONTestTrigger (Int_t evNumber1=0,Int_t evNumber2=0) 
{
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (file) file->Close(); 
   file = new TFile("galice.root","READ");

// Get AliRun object from file or create it if not on file
  printf ("I'm after Map \n");
  if (!gAlice) { 
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }
  printf ("I'm after gAlice \n");
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  for (int nev=evNumber1; nev<= evNumber2; nev++) {
    Int_t nparticles = gAlice->GetEvent(nev);
    cout << "nev npart =" << nev << " , " << nparticles << "\n";
    if (nev < evNumber1) continue;
    if (nparticles <= 0) return;
                     
    Int_t nbytes = 0;
    AliMUONGlobalTrigger *gloTrg;
    AliMUONLocalTrigger *locTrg;

// Get pointers to Alice detectors and Triggers containers
    AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
    TClonesArray *globalTrigger = MUON->GlobalTrigger();
    TClonesArray *localTrigger = MUON->LocalTrigger();

    TTree *TR = gAlice->TreeR();    
    TBranch *gbranch = TR->GetBranch("MUONGlobalTrigger");
    TBranch *lbranch = TR->GetBranch("MUONLocalTrigger");
    gbranch->SetAddress(&globalTrigger);
    lbranch->SetAddress(&localTrigger);

    Int_t nent = TR->GetEntries();
    cout << ">>> " << nent << " entries found in TreeR of event " 
	 << nev << "\n";
    Int_t nb = 0;

    for (Int_t n=1; n<nent; n++) {
      MUON->ResetTrigger();
      nbytes += TR->GetEvent(n);

      Int_t nglobals = globalTrigger->GetEntries(); // should be = 1
      Int_t nlocals  = localTrigger->GetEntries();  // who knows ?
      
      for (Int_t i=0; i<nglobals; i++) { // inspect Global Trigger
	cout << " >>> Output for Global Trigger " << "\n";
	
	gloTrg = (AliMUONGlobalTrigger*)globalTrigger->UncheckedAt(i);
	
	cout << "fSinglePlusLpt = " << gloTrg->SinglePlusLpt() << "\n"; 
	cout << "fSinglePlusHpt = " << gloTrg->SinglePlusHpt() << "\n";
	cout << "fSinglePlusApt = " << gloTrg->SinglePlusApt() << "\n";

	cout << "fSingleMinusLpt = " << gloTrg->SingleMinusLpt() << "\n";
	cout << "fSingleMinusHpt = " << gloTrg->SingleMinusHpt() << "\n";
	cout << "fSingleMinusApt = " << gloTrg->SingleMinusApt() << "\n";

	cout << "fSingleUndefLpt = " << gloTrg->SingleUndefLpt() << "\n";
	cout << "fSingleUndefHpt = " << gloTrg->SingleUndefHpt() << "\n";
	cout << "fSingleUndefApt = " << gloTrg->SingleUndefApt() << "\n";

	cout << "fPairLikeLpt = " << gloTrg->PairLikeLpt() << "\n";
	cout << "fPairLikeHpt = " << gloTrg->PairLikeHpt() << "\n";
	cout << "fPairLikeApt = " << gloTrg->PairLikeApt() << "\n";

	cout << "fPairUnlikeLpt = " << gloTrg->PairUnlikeLpt() << "\n";
	cout << "fPairUnlikeHpt = " << gloTrg->PairUnlikeHpt() << "\n";
	cout << "fPairUnlikeApt = " << gloTrg->PairUnlikeApt() << "\n";
      } // end of loop on Global Trigger

      for (Int_t i=0; i<nlocals; i++) { // inspect Local Trigger
	cout << " >>> Output for Local Trigger # " << i << "\n";

	locTrg = (AliMUONLocalTrigger*)localTrigger->UncheckedAt(i);	

	cout << "fLoCircuit = " << locTrg->LoCircuit() << "\n";
	cout << "fLoStripX = "  << locTrg->LoStripX() << "\n";
	cout << "fLoDev = "     << locTrg->LoDev() << "\n";
	cout << "fLoStripY = "  << locTrg->LoStripY() << "\n";
	cout << "fLoLpt = "     << locTrg->LoLpt() << "\n";
	cout << "fLoHpt = "     << locTrg->LoHpt() << "\n";
	cout << "fLoApt = "     << locTrg->LoApt() << "\n";

      } // end of loop on Local Trigger
    } // end of loop on entries of TreeR
  } // loop on event
  
} 




