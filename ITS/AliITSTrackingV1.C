#ifndef __CINT__
#include "iostream.h"
#endif

void AliITSTrackingV1(Int_t evNumber1=0,Int_t evNumber2=0, Int_t min_t=-1, Int_t max_t=0,Bool_t flagvert=1, Bool_t realmass=0, const char *filename="galice.root") {

  
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

  AliITS* IITTSS =(AliITS *)gAlice->GetDetector("ITS");        
  if (!IITTSS) return;                                           

//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (Int_t nev=0; nev<= evNumber2; nev++) {
	AliITSTrackerV1 *ITStracker = new AliITSTrackerV1(IITTSS,nev,flagvert);
     Int_t nparticles = gAlice->GetEvent(nev);
     cout << "nev         " << nev <<endl;
     cout << "nparticles  " << nparticles <<endl;
     if (nev < evNumber1) continue;
     if (nparticles <= 0) return;

     TTree *TR=gAlice->TreeR();
     Int_t nent=TR->GetEntries();
     //printf("Found %d entries in the TreeR (must be one per module per event!)\n",nent);


     TStopwatch timer;
	  
	  timer.Start();
     ITStracker->DoTracking(nev,min_t,max_t,file,realmass);    
     timer.Stop(); timer.Print();
	 AliITSgeom *g1 = IITTSS->GetITSgeom();
    Int_t NumOfModules = g1->GetIndexMax();	   
	  ITStracker->DelMatrix(NumOfModules);
	  delete ITStracker;  
   }   // event loop  
   file->Close();   
}

