#include "iostream.h"

void ITSdigitBari (Int_t evNumber1=0,Int_t evNumber2=0,Int_t nsignal  =25, Int_t size=-1) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//  
//   Macro to run the Bari model on the SPDs
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

   AliITSgeom *geom = ITS->GetITSgeom();

//
// Event Loop
//


   // SPD
   printf ("Beginning of SPD processing \n\n");

   AliITSDetType *iDetType=ITS->DetType(0);
   AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
   //   AliITSresponseSPDbari *res0=(AliITSresponseSPDbari*)iDetType->GetResponseModel();
 AliITSresponseSPDbari *res0= new AliITSresponseSPDbari();

// to change or monitor the parameters

//--occhio
   res0->SetThresholds(7.2e-6, 1.e-6);
   Float_t thresh, sigma;
   res0->Thresholds(thresh, sigma);
   printf("SPDbari: threshold %e sigma %e\n",thresh, sigma);

//--occhio
    res0->SetNoiseParam(0., 0.);
//   res0->SetNoiseParam(0.04, 0.08);
   Float_t col, row;
   res0->GetNoiseParam(col, row);
   printf("SPDbari: Couplcol %e Couplrow %e\n",col, row);

   ITS->SetResponseModel(0,res0);

   
   AliITSsimulationSPDbari *sim0=new AliITSsimulationSPDbari(seg0,res0);
   ITS->SetSimulationModel(0,sim0);
   
   // test
   printf("SPD dimensions %f %f \n",seg0->Dx(),seg0->Dz());
   printf("SPD npixels %d %d \n",seg0->Npz(),seg0->Npx());
   printf("SPD pitches %d %d \n",seg0->Dpz(0),seg0->Dpx(0));
   // end test

   printf("End SPD processing \n\n");





   Int_t nbgr_ev=0;

   for (int nev=evNumber1; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       cout << "nev         " <<nev<<endl;
       cout << "nparticles  " <<nparticles<<endl;
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;

       Int_t nbgr_ev=Int_t(nev/nsignal);
       //printf("nbgr_ev %d\n",nbgr_ev);
//       ITS->HitsToDigits(nev,nbgr_ev,evNumber2,size," ","All"," ");
       ITS->HitsToDigits(nev,nbgr_ev,size," ","SPD"," ");
   } // event loop 

   file->Close();
}

