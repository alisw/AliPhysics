#include "iostream.h"

void ITSHitsToFastPoints (Int_t evNumber1=0,Int_t evNumber2=0,Int_t nsignal=25, Int_t size=-1) 
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
   file->ls();

   printf ("I'm after Map \n");

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
       gAlice = (AliRun*)file->Get("gAlice");
       if (gAlice) printf("AliRun object found on file\n");
       if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   printf ("I'm after gAlice \n");
   
   AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
   if (!ITS) return;

   // Set the simulation model
   AliITSsimulationFastPoints *sim = new AliITSsimulationFastPoints();

   for (Int_t i=0;i<3;i++) {
       ITS->SetSimulationModel(i,sim);
   }

//
// Event Loop
//

   Int_t nbgr_ev=0;

   for (int ev=evNumber1; ev<= evNumber2; ev++) {
       Int_t nparticles = gAlice->GetEvent(ev);
       cout << "event         " <<ev<<endl;
       cout << "nparticles  " <<nparticles<<endl;
       if (ev < evNumber1) continue;
       if (nparticles <= 0) return;

       Int_t bgr_ev=Int_t(ev/nsignal);
       //printf("bgr_ev %d\n",bgr_ev);
       ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","All"," ");
   } // event loop 

   file->Close();
}














