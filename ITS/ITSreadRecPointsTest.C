#include "iostream.h"

void ITSreadRecPointsTest (Int_t evNumber1=0,Int_t evNumber2=0) 
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

     AliITSRecPoint  *recp;

     // Get pointers to Alice detectors and Digits containers
     AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
     if (!ITS) return;
     TClonesArray *Particles = gAlice->Particles();
     TTree *TR = gAlice->TreeR();
     Int_t nent=TR->GetEntries();
     printf("Found %d entries in the tree (must be one per module per event!)\n",nent);

     TClonesArray *ITSrec  = ITS->RecPoints();
     for (Int_t mod=0; mod<nent; mod++) {
       ITS->ResetRecPoints();
       TR->GetEvent(mod);
       Int_t nrecp = ITSrec->GetEntries();
       if (nrecp) printf("Found %d rec points for module %d \n",nrecp,mod);
       if (!nrecp) continue;

       for (Int_t irec=0;irec<nrecp;irec++) {
		recp   = (AliITSRecPoint*)ITSrec->UncheckedAt(irec);
		printf("%d %f %f %d %d %d\n",irec,recp->GetX(),recp->GetZ(),recp->GetLabel(0),recp->GetLabel(1),recp->GetLabel(2));

       }
     }        
   }   // event loop 
   
     cout<<"END  test for rec points "<<endl;

     file->Close();   
}



