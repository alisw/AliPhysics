#include "iostream.h"

void ITSDigitsToClusters (Int_t evNumber1=0,Int_t evNumber2=0) 
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
   } else {
      delete gAlice;
      gAlice=0;
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

   AliITSgeom *geom = ITS->GetITSgeom();


   // Set the models for cluster finding

   // SPD

   AliITSDetType *iDetType=ITS->DetType(0);
   AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
   TClonesArray *dig0  = ITS->DigitsAddress(0);
   TClonesArray *recp0  = ITS->ClustersAddress(0);
   AliITSClusterFinderSPD *rec0=new AliITSClusterFinderSPD(seg0,dig0,recp0);
   ITS->SetReconstructionModel(0,rec0);
   // test
   printf("SPD dimensions %f %f \n",seg0->Dx(),seg0->Dz());
   printf("SPD npixels %d %d \n",seg0->Npz(),seg0->Npx());


   // SDD

   AliITSDetType *iDetType=ITS->DetType(1);
   AliITSgeom *geom = ITS->GetITSgeom();

   AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->GetSegmentationModel();
   if (!seg1) seg1 = new AliITSsegmentationSDD(geom);
   AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
   if (!res1) res1=new AliITSresponseSDD();
   res1->SetNoiseParam(0.,0.);
   TClonesArray *dig1  = ITS->DigitsAddress(1);
   TClonesArray *recp1  = ITS->ClustersAddress(1);
   AliITSClusterFinderSDD *rec1=new AliITSClusterFinderSDD(seg1,res1,dig1,recp1);
   ITS->SetReconstructionModel(1,rec1);



   // SSD

   AliITSDetType *iDetType=ITS->DetType(2);
   AliITSsegmentationSSD *seg2=(AliITSsegmentationSSD*)iDetType->GetSegmentationModel();
   TClonesArray *dig2  = ITS->DigitsAddress(2);
   TClonesArray *recp2  = ITS->ClustersAddress(2);
   AliITSClusterFinderSSD *rec2=new AliITSClusterFinderSSD(seg2,dig2,recp2);
   ITS->SetReconstructionModel(2,rec2);
   // test
   printf("SSD dimensions %f %f \n",seg2->Dx(),seg2->Dz());
   printf("SSD nstrips %d %d \n",seg2->Npz(),seg2->Npx());



//
// Event Loop
//

   for (int nev=evNumber1; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       cout << "nev         " <<nev<<endl;
       cout << "nparticles  " <<nparticles<<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;

       TTree *TD = gAlice->TreeD();
       Int_t nent=TD->GetEntries();
       printf("Found %d entries in the tree (must be one per module per event!)\n",nent);
       //Int_t nmodules=geom->GetLastSSD();
       //Int_t last_entry=nent-(nmodules+1);
       Int_t last_entry=1;
       ITS->DigitsToRecPoints(nev,last_entry,"All");
   } // event loop 

   file->Close();
}














