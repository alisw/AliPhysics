#include "iostream.h"

void ITSmixedpoints (Int_t evNumber1=0,Int_t evNumber2=0,Int_t nsignal  =25, Int_t size=-1) 
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

   // set the simulation models
   AliITSgeom *geom = ITS->GetITSgeom();

   // SPD - simulation slow points

   AliITSDetType *iDetType=ITS->DetType(0);
   AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
   AliITSresponseSPD *res0 = (AliITSresponseSPD*)iDetType->GetResponseModel();
   AliITSsimulationSPD *sim0=new AliITSsimulationSPD(seg0,res0);
   ITS->SetSimulationModel(0,sim0);

   // SPD - cluster finder
   AliITSDetType *iDetType=ITS->DetType(0);
   AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
   TClonesArray *dig0  = ITS->DigitsAddress(0);
   TClonesArray *recp0  = ITS->ClustersAddress(0);
   AliITSClusterFinderSPD *rec0=new AliITSClusterFinderSPD(seg0,dig0,recp0);
   ITS->SetReconstructionModel(0,rec0);

   // SDD+SSD - fast poinst

   AliITSsimulationFastPoints *sim = new AliITSsimulationFastPoints();

   for (Int_t i=1;i<3;i++) {
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
       Int_t last_entry=1;
       ITS->DigitsToRecPoints(ev,last_entry,"SPD");
       ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","SDD"," ");
       ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","SSD"," ");
   } // event loop 

   file->Close();
}














