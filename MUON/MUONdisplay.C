MUONdisplay (Int_t nevent=0) {
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
      /*
   } else {
      delete gAlice;
      gAlice = 0;
      */
   }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   
// Create Event Display object
   AliMUONDisplay *muondisplay = new AliMUONDisplay(750);

// Display first event
   gAlice->GetEvent(nevent);
   muondisplay->ShowNextEvent(0);
}
