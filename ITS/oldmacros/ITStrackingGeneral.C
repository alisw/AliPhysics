#include "iostream.h"

void ITStrackingGeneral(Int_t evNumber1=0,Int_t evNumber2=0,Int_t min_t=-1,Int_t max_t=0,Bool_t flagvert=1) {

  const char *filename="galice.root";
  
  ///////////////// Dynamically link some shared libs ////////////////////////////////
  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice=0;
  }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile("galice.root","UPDATE");
   //if (!file) file = new TFile(filename);

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }


  AliITS* ITS =(AliITS *)gAlice->GetDetector("ITS");
  if (!ITS) return;

//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (int nev=0; nev<= evNumber2; nev++) {
     Int_t nparticles = gAlice->GetEvent(nev);
     cout << "nev         " << nev <<endl;
     cout << "nparticles  " << nparticles <<endl;
     if (nev < evNumber1) continue;
     if (nparticles <= 0) return;

     TTree *TR=gAlice->TreeR();
     Int_t nent=TR->GetEntries();
     //printf("Found %d entries in the TreeR (must be one per module per event!)\n",nent);

     ITS->DoTracking(nev,min_t,max_t,file,flagvert);
   }   // event loop 
   file->Close();   
}

