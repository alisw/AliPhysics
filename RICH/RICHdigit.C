// 0 = all
// 1 = not pion
// 2 = not kaon
// 3 = not proton
// 4 = not muon
// 5 = not electron
// 6 = not neutron


Int_t particle_type=0;

#include "iostream.h"

void RICHdigit (Int_t evNumber1=0,Int_t evNumber2=0) 
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

   AliRICHChamber*       iChamber;
   
   printf("Generating tresholds...\n");

   for(Int_t i=0;i<7;i++)
     {
       iChamber = &(RICH->Chamber(i));
       iChamber->GenerateTresholds();
     }
   
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
	   RICH->Digitise(nev, particle_type);
	 }
       //if (RICH) gAlice->SDigits2Digits("RICH");
       //char hname[30];
       //sprintf(hname,"TreeD%d",nev);
       //gAlice->TreeD()->Write(hname);
       //gAlice->TreeD()->Reset();
   } // event loop 
   file->Close();

   //delete gAlice;
   printf("\nEnd of Macro  *************************************\n");
}














