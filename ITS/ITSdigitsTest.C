#include "iostream.h"

void ITSdigitsTest (Int_t evNumber1=0,Int_t evNumber2=0) 
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

   AliITSdigitSSD  *ITSdigit;

// Get pointers to Alice detectors and Digits containers
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   TClonesArray *Particles = gAlice->Particles();
   TTree *TD = gAlice->TreeD();
   TD->Print();
   Int_t nent=TD->GetEntries();
   printf("Found %d entries in the tree (must be one per module per event!)\n",nent);
   if (ITS) {
   TObjArray *fBranches=TD->GetListOfBranches();
   
     for (Int_t ich=0;ich<3;ich++) {
       TBranch *branch = (TBranch*)fBranches->UncheckedAt(ich);
       printf ("branch %p \n",branch);
       printf ("branch %p entries %d \n",branch,branch->GetEntries());
       TClonesArray *ITSdigits  = ITS->DigitsAddress(ich);
          printf ("ITSdigits %p \n",ITSdigits);

          if (ich != 2) continue;
	  for (Int_t mod=1; mod<nent; mod++) {
	    //Int_t nmodules=2269;
	    //for (Int_t mod=nent-nmodules; mod<nent; mod++) {
              ITS->ResetDigits();
              nbytes += TD->GetEvent(mod);
              //nbytes += branch->GetEvent(mod); this works as well
	      Int_t ndigits = ITSdigits->GetEntries();
	      if (ndigits) printf("Found %d digits for module %d in det type %d \n",ndigits,mod,ich+1);

	      if (!ndigits) continue;
	      /*
	      for (Int_t digit=0;digit<ndigits;digit++) {
		ITSdigit   = (AliITSdigitSPD*)ITSdigits->UncheckedAt(digit);
		printf("%d %d %d %d \n",ITSdigit->fCoord1,ITSdigit->fCoord2,ITSdigit->fTracks[0],ITSdigit->fTracks[1]);
	      */
	      /*
	      for (Int_t digit=0;digit<ndigits;digit++) {
		ITSdigit   = (AliITSdigitSDD*)ITSdigits->UncheckedAt(digit);
		printf("%d %d %d %d %f %f\n",ITSdigit->fCoord1,ITSdigit->fCoord2,ITSdigit->fSignal,ITSdigit->fTracks[0],ITSdigit->fTcharges[0],ITSdigit->fPhysics);
	      */
	      for (Int_t digit=0;digit<ndigits;digit++) {
		ITSdigit   = (AliITSdigitSSD*)ITSdigits->UncheckedAt(digit);
		printf("%d %d %d %d \n",ITSdigit->fCoord1,ITSdigit->fCoord2,ITSdigit->fSignal,ITSdigit->fTracks[0]);

	      }
	  }        
     }
   }   // end if ITS

   }   // event loop 

     cout<<"END  test for digits "<<endl;

     file->Close();   
}



