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
	
	cout << "fSinglePlusLpt = " << gloTrg->fSinglePlusLpt << "\n"; 
	cout << "fSinglePlusHpt = " << gloTrg->fSinglePlusHpt << "\n";
	cout << "fSinglePlusApt = " << gloTrg->fSinglePlusApt << "\n";

	cout << "fSingleMinusLpt = " << gloTrg->fSingleMinusLpt << "\n";
	cout << "fSingleMinusHpt = " << gloTrg->fSingleMinusHpt << "\n";
	cout << "fSingleMinusApt = " << gloTrg->fSingleMinusApt << "\n";

	cout << "fSingleUndefLpt = " << gloTrg->fSingleUndefLpt << "\n";
	cout << "fSingleUndefHpt = " << gloTrg->fSingleUndefHpt << "\n";
	cout << "fSingleUndefApt = " << gloTrg->fSingleUndefApt << "\n";

	cout << "fPairLikeLpt = " << gloTrg->fPairLikeLpt << "\n";
	cout << "fPairLikeHpt = " << gloTrg->fPairLikeHpt << "\n";
	cout << "fPairLikeApt = " << gloTrg->fPairLikeApt << "\n";

	cout << "fPairUnlikeLpt = " << gloTrg->fPairUnlikeLpt << "\n";
	cout << "fPairUnlikeHpt = " << gloTrg->fPairUnlikeHpt << "\n";
	cout << "fPairUnlikeApt = " << gloTrg->fPairUnlikeApt << "\n";
      } // end of loop on Global Trigger

      for (Int_t i=0; i<nlocals; i++) { // inspect Local Trigger
	cout << " >>> Output for Local Trigger # " << i << "\n";

	locTrg = (AliMUONLocalTrigger*)localTrigger->UncheckedAt(i);	

	cout << "fLoCircuit = " << locTrg->fLoCircuit << "\n";
	cout << "fLoStripX = "  << locTrg->fLoStripX << "\n";
	cout << "fLoDev = "     << locTrg->fLoDev << "\n";
	cout << "fLoStripY = "  << locTrg->fLoStripY << "\n";
	cout << "fLoLpt = "     << locTrg->fLoLpt << "\n";
	cout << "fLoHpt = "     << locTrg->fLoHpt << "\n";
	cout << "fLoApt = "     << locTrg->fLoApt << "\n";

      } // end of loop on Local Trigger
    } // end of loop on entries of TreeR
  } // loop on event
  
// store histos in ouput file
    hfile->Write();
} 




