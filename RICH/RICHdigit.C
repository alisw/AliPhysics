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
   }


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
    else {
      //delete gAlice;
      gAlice = 0;
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
//
// Event Loop
//
   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       cout <<endl<< "Processing event:" <<nev<<endl;
       cout << "Particles       :" <<nparticles<<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       if (RICH) RICH->Digitise(nev);
       char hname[30];
       sprintf(hname,"TreeD%d",nev);
       gAlice->TreeD()->Write(hname);
       gAlice->TreeD()->Reset();
   } // event loop 
   file->Close();

   //delete gAlice;
   printf("\nEnd of Macro  *************************************\n");
}














