void MUONdigit (Int_t evNumber1=0, Int_t evNumber2=0, Int_t ibg=0, Int_t bgr=10) 
{
  //////////////////////////////////////
  //                                  //
  // ROOT macro for ALICE Dimuon Arm: //
  // Digitization                     //
  //                                  //
  //////////////////////////////////////
  //
  // Adds the tree TD for digits (charges deposited on pads)
  // to the ROOT file "galice.root"
  // containing the signal hit coordinates from simulation (tree TH).
  // Eventually (argument "ibg"), background hits are also taken into account,
  // from the ROOT file "bg.root".
  //
  // Arguments:
  //   evNumber1 = first event number to digitize in file "galice.root"
  //   evNumber2 = last event number to digitize in file "galice.root"
  //   ibg       = 0 if no background hits to be taken into account;
  //             = 1 if  background hits to be taken into account
  //                     in file "bg.root"
  //   bgr       : used only if "ibg" = 1
  //             = number of events in the background file "bg.root";
  //               the signal events are divided into "bgr" successive parts,
  //               all events of each part are associated
  //               with the same background event,
  //               starting with event number 0,
  //               incrementing it by 1 for the next part.
  //               Strictly speaking, "bgr" can be smaller than
  //               the number of events in the background file,
  //               in which case one will only use
  //               the first "bgr" events of the background file.
  //               But it SHOULD NOT BE LARGER THAN
  //               THE NUMBER OF EVENTS IN THE BACKGROUND FILE.
  //
  // Input file(s):
  //   "galice.root" for signal
  //   "bg.root" for background (used only if "ibg" = 1)
  //
  // Output file:
  //   "galice.root"
  //
  //__________________________________________________________________________

// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (file) file->Close(); 
   file = new TFile("galice.root","UPDATE");

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
       gAlice = (AliRun*)file->Get("gAlice");
       if (gAlice) printf("AliRun object found on file\n");
       if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   AliMUON *pMUON  = (AliMUON*) gAlice->GetModule("MUON");
   if (pMUON) {
// creation
       AliMUONMerger* merger = new AliMUONMerger();
// configuration
       if (ibg) {
	 merger->SetMode(ibg);
	 merger->SetBackgroundFileName("bg.root");
       }
       // pass
       pMUON->SetMerger(merger);
   }
// Action !
//
//   Loop over events              
//
    for (int nev=evNumber1; nev<= evNumber2; nev++) {
	Int_t nparticles = gAlice->GetEvent(nev);
	cout << "nev         " << nev <<endl;
	cout << "nparticles  " << nparticles <<endl;
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;
	Int_t nbgr_ev = Int_t(nev*bgr/(evNumber2+1));
	
	if (ibg) {
	    merger->SetBackgroundEventNumber(nbgr_ev);
	}

	gAlice->SDigits2Digits();

	char hname[30];
	sprintf(hname,"TreeD%d",nev);
	gAlice->TreeD()->Write(hname);
	// reset tree
	gAlice->TreeD()->Reset();

    }   // event loop 
}














