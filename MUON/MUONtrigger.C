#include "iostream.h"
//Add the Trigger output in the tree TR of the ROOT file galice.root

void MUONtrigger (Int_t evNumber1=0, Int_t evNumber2=0)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
   printf ("I'm after gAlice \n");
   
   AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   for (int nev=evNumber1; nev<= evNumber2; nev++) { // event loop
     Int_t nparticles = gAlice->GetEvent(nev);
     cout << "nev         " <<nev<<endl;
     cout << "nparticles  " <<nparticles<<endl; 
     if (nparticles <= 0) return;
     
     if (MUON) MUON->Trigger(nev);
     
   } // event loop 
   file->Close();
}














