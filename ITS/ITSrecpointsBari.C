#include "iostream.h"

void ITSrecpointsBari (Int_t evNumber1=0,Int_t evNumber2=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
//   version to execute the Bari model
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }


// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   printf("file %p\n",file);
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

//
// Event Loop
//


   AliITSgeom *geom = ITS->GetITSgeom();

   // SPD

   AliITSDetType *iDetType=ITS->DetType(0);
   AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
   TClonesArray *dig0  = ITS->DigitsAddress(0);
   TClonesArray *recp0  = ITS->ClustersAddress(0);
   AliITSClusterFinderSPDbari *rec0=new AliITSClusterFinderSPDbari(seg0,dig0,recp0);
   ITS->SetReconstructionModel(0,rec0);

   // test
   printf("SPD dimensions %f %f \n",seg0->Dx(),seg0->Dz());
   printf("SPD npixels %d %d \n",seg0->Npz(),seg0->Npx());



   for (int nev=evNumber1; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       cout << "nev         " <<nev<<endl;
       cout << "nparticles  " <<nparticles<<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;

       TTree *TD = gAlice->TreeD();
       Int_t nent=TD->GetEntries();
       printf("Found %d entries in the tree (must be one per module per event!)\n",nent);
       Int_t nmodules=geom->GetLastSSD();
       Int_t last_entry=1;
//       ITS->DigitsToRecPoints(nev,last_entry,"All");
       ITS->DigitsToRecPoints(nev,last_entry,"SPD");
   } // event loop 

   file->Close();
}

