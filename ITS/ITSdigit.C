void ITSdigit (Int_t evNumber=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   for the ITS digits
//   
//     Root > .L ITSdigit.C   //this loads the macro in memory
//     Root > ITSdigit();     //by default process first event   
//     Root > ITSdigit(2);    //process third event
//End_Html
/////////////////////////////////////////////////////////////////////////


// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gSystem->Load("libGeant3Dummy.so");   // a dummy version of Geant3
      gSystem->Load("PHOS/libPHOSdummy.so");        // the standard Alice classes 
      gSystem->Load("libgalice.so");        // the standard Alice classes 
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
      
// Import the Kine and Hits Trees for the event evNumber in the file
   Int_t nparticles = gAlice->GetEvent(evNumber);
   if (nparticles <= 0) return;
   Float_t x,y,z,mass,e;
   Int_t nbytes = 0;
   Int_t j,hit,ipart;
   Int_t nhits;
   Int_t sector,plane;
   GParticle *particle;
   AliITSdigit  *ITSdigit;

// Get pointers to Alice detectors and Hits containers
   AliDetector *ITS  = gAlice->GetDetector("ITS");
   TClonesArray *Particles = gAlice->Particles();
   if (ITS) TClonesArray *ITSdigits  = ITS->Digits();

   TTree *TD = gAlice->TreeD();
   Int_t nent    = TD->GetEntries();

   printf("Found %d entries in the tree (must be one per event!)\n",nent);
   
   for (Int_t dig=0; dig<nent; dig++) {
     gAlice->ResetDigits();
     nbytes += TD->GetEvent(dig);
     if (ITS) {
       Int_t ndigits = ITSdigits->GetEntriesFast();
       printf("Found %d digits\n",ndigits);
       for (Int_t digit=0;digit<ndigits;digit++) {
	 ITSdigit   = (AliITSdigit*)ITSdigits->UncheckedAt(digit);
	 printf("%d %d %d %d \n",ITSdigit->fEvent,ITSdigit->fLayer,ITSdigit->fDet,ITSdigit->fNoverl);
       }
     }        
   }
}
