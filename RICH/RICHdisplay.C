RICHdisplay (Int_t nevent=0) {
// Dynamically link some shared libs
 
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    }
    else {
      delete gAlice;
      gAlice = 0;
    }
    
    galice=0;
    
// Connect the Root Galice file containing Geometry, Kine and Hits
    
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (file) file->Close(); 
    file = new TFile("galice.root","UPDATE");
       
    //printf ("I'm after Map \n");
    
// Get AliRun object from file or create it if not on file
    
    if (!gAlice) {
	gAlice = (AliRun*)file->Get("gAlice");
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
    else {
      delete gAlice;
      gAlice = (AliRun*)file->Get("gAlice");
      	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
    
    //printf ("I'm after gAlice \n");
    
    // Create Event Display object
    AliRICHDisplay *richdisplay = new AliRICHDisplay(750);
    
// Display first event
    gAlice->GetEvent(nevent);
    richdisplay->ShowNextEvent(0);
    
    
    //file->Close();
    //delete file;
    //delete richdisplay;
}
