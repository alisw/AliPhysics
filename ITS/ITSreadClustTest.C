#include "iostream.h"

void ITSreadClustTest (Int_t evNumber1=0,Int_t evNumber2=0) 
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
   if (!file) file = new TFile("galice.root");
   file->ls();

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
 
//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (int nev=0; nev<= evNumber2; nev++) {
     Int_t nparticles = gAlice->GetEvent(nev);
     cout << "nev         " << nev <<endl;
     cout << "nparticles  " << nparticles <<endl;
     if (nev < evNumber1) continue;
     if (nparticles <= 0) return;

     TTree *TH = gAlice->TreeH();
     Int_t ntracks = TH->GetEntries();
     cout<<"ntracks "<<ntracks<<endl;

   Int_t nbytes = 0;

   //AliITSRawClusterSDD  *ITSclust;
   AliITSRawClusterSPD  *ITSclust;

// Get pointers to Alice detectors and Digits containers
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if(!ITS) return;
   TClonesArray *Particles = gAlice->Particles();

     ITS->GetTreeC(nev);
     TTree *TC=ITS->TreeC();
     Int_t nent=TC->GetEntries();
     printf("Found %d entries in the tree (must be one per module per event!)\n",nent);
   
     for (Int_t idet=0;idet<3;idet++) {

       TClonesArray *ITSclu  = ITS->ClustersAddress(idet);

          if (idet != 0) continue;
	  for (Int_t mod=0; mod<nent; mod++) {
              ITS->ResetClusters();
              TC->GetEvent(mod);
	      Int_t ncl = ITSclu->GetEntries();
	      if (ncl) printf("Found %d clusters for module %d in det type %d \n",nrecp,mod,idet);

	      if (!ncl) continue;

	      for (Int_t icl=0;icl<ncl;icl++) {
		ITSclust   = (AliITSRawClusterSPD*)ITSclu->UncheckedAt(icl);
		printf("%d %d %f %f %f\n",ITSclust->NclZ(),ITSclust->NclX(),ITSclust->Q(),ITSclust->X(),ITSclust->Z());
		/*
	      for (Int_t icl=0;icl<ncl;icl++) {
		ITSclust   = (AliITSRawClusterSDD*)ITSclu->UncheckedAt(icl);
		printf("%d %d %f %f %f\n",ITSclust->Anodes(),ITSclust->Samples(),ITSclust->Q(),ITSclust->X(),ITSclust->Z());
		*/

	      }
	  }        
     }

   }   // event loop 

     cout<<"END  test for clusters "<<endl;

     file->Close();   
}










