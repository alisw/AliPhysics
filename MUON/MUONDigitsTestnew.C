#include "iostream.h"

void MUONdigitsTestnew (Int_t evNumber1=0,Int_t evNumber2=0) 
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

   AliMUONdigit  *MUONdigit;

// Get pointers to Alice detectors and Digits containers
   AliMUON *MUON  = gAlice->GetDetector("MUON");
   TClonesArray *Particles = gAlice->Particles();
   TTree *TD = gAlice->TreeD();
   Int_t nent=TD->GetEntries();
   printf("Found %d entries in the tree (must be one per cathode per event!)\n",nent);
   if (MUON) {
     for (Int_t ich=0;ich<10;ich++) {
        TClonesArray *MUONdigits  = MUON->DigitsAddress(ich);
        //   printf ("MUONdigits %f \n",MUONdigits);

        for (Int_t dig=1; dig<nent; dig++) {
          gAlice->ResetDigits();
          nbytes += TD->GetEvent(dig);
          Int_t ndigits = MUONdigits->GetEntries();
          printf("Found %d digits for cathode %d in chamber %d \n",ndigits,dig,ich+1);
	  for (Int_t digit=0;digit<ndigits;digit++) {
  	      MUONdigit   = (AliMUONdigit*)MUONdigits->UncheckedAt(digit);
	      //	       	 printf("%d %d %d %d %d\n",MUONdigit->fPadX,MUONdigit->fPadY,MUONdigit->fSignal,MUONdigit->fTracks[0],MUONdigit->fTcharges[0]);
          }
        }        
     }
   }   // end if MUON

   }   // event loop 

     cout<<"END  test for digits "<<endl;

     file->Close();   
}



