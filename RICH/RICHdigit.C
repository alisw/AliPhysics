// 0 = all
// 1 = not pion
// 2 = not kaon
// 3 = not proton
// 4 = not muon
// 5 = not electron
// 6 = not neutron


Int_t particle_type=0;

#include "iostream.h"

void RICHdigit (Int_t evNumber1=0,Int_t evNumber2=0, Int_t merging) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs

   
  if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
  }else {
    delete gAlice;
    gAlice = 0;
  }
  
  galice=0;

// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (file) file->Close(); 
   file = new TFile("galice.root","UPDATE");
//   file->ls();
// Get AliRun object from file or create it if not on file

   

   if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    }


   if (!gAlice) {
       gAlice = (AliRun*)file->Get("gAlice");
       if (gAlice) printf("AliRun object found on file\n");
       if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   } else {
      delete gAlice;
      gAlice = (AliRun*)file->Get("gAlice");
      	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
      
   AliRICH *RICH  = (AliRICH*) gAlice->GetDetector("RICH");

   if (merging)
     printf("Merging is ON\n");
   else
     printf("Merging is OFF\n");

// Creation of merger object
   AliRICHMerger* merger = new AliRICHMerger();
   
// Configuration
   merger->SetMode(merging);
   merger->SetSignalEventNumber(0);
   merger->SetBackgroundEventNumber(0);
   merger->SetBackgroundFileName("bg.root");
       
// Pass
   RICH->SetMerger(merger);


//
// Event Loop
//
   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       cout <<endl<< "Processing event:" <<nev<<endl;
       cout << "Particles       :" <<nparticles<<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       if (RICH) 
	 {
	   //gAlice->MakeTree("D");
	   //RICH->MakeBranch("D");
	   //RICH->Digitise(nev, particle_type);
	   //gAlice->SDigits2Digits("RICH");
	   //gAlice->Tree2Tree("D");
	   RICH->MakeBranch("D");
	   RICH->SDigits2Digits(nev, particle_type);
	 }
   } // event loop 
   file->Close();

   //delete gAlice;
   printf("\nEnd of Macro  *************************************\n");
}



