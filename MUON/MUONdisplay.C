MUONdisplay (Int_t nevent=0) {
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
       gSystem->Load("$ALITOP/cern.so/lib/libpdfDUMMY.so");
       gSystem->Load("$ALITOP/cern.so/lib/libPythia.so");
       gSystem->Load("$ROOTSYS/lib/libEG.so");       
       gSystem->Load("$ROOTSYS/lib/libEGPythia.so");    
       gSystem->Load("libGeant3Dummy.so");        //a dummy version of Geant3
       gSystem->Load("PHOS/libPHOSdummy.so");     //the standard Alice classes 
       gSystem->Load("libgalice.so");             //the standard Alice classes 
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
   AliMUONdisplay *muondisplay = new AliMUONdisplay(750);

// Display first event
   gAlice->GetEvent(nevent);
   muondisplay->ShowNextEvent(0);
}
