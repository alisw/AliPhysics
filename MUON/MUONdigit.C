void MUONdigit (Int_t evNumber1=0, Int_t evNumber2=0, Int_t ibg=1, Int_t bgr=10) 
{
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
       merger->SetMode(0);
       merger->SetSignalEventNumber(0);
       merger->SetBackgroundEventNumber(0);
       merger->SetBackgroundFileName("bg.root");
       
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
	gAlice->SDigits2Digits();

	char hname[30];
	sprintf(hname,"TreeD%d",nev);
	gAlice->TreeD()->Write(hname);
	// reset tree
	gAlice->TreeD()->Reset();

    }   // event loop 
}














