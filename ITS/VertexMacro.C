#include "iostream.h"

void VertexMacro(Int_t evNumber1=0,Int_t evNumber2=0) {

  const char *filename="galice.root";
  
  ///////////////// Dynamically link some shared libs ////////////////////////////////
  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice=0;
  }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile(filename,"UPDATE");


// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }

//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (int nev=0; nev<= evNumber2; nev++) {
     gAlice->SetEvent(nev);
	 Int_t nparticles = gAlice->GetEvent(nev);
     cout << "nev         " << nev <<endl;
     cout << "nparticles  " << nparticles <<endl;
     if (nev < evNumber1) continue;
     if (nparticles <= 0) return;
  
  
     TStopwatch timer;
     timer.Start();

     AliITSVertex *V = new AliITSVertex();;

     timer.Stop();
     timer.Print();

     cout << endl << "Zv = " << V->GetZv() << " cm" << endl;
     cout << "Z Resolution = " << V->GetZRes()*10000 << " microns" << endl;
     cout << "Signal/Noise for Z = " << V->GetZSNR() <<endl;
     cout << endl << "Yv (MC value) = " << V->GetYv() << " cm"  << endl;
//     cout << "Y resolution = " << V->GetYRes()*10000 << " microns"  << endl;
//     cout << "Signal/Noise for Y = " << V->GetYSNR() << endl;
     cout << endl << "Xv (MC value) = " << V->GetXv() << " cm" << endl;
//     cout << "X resolution = " << V->GetXRes()*10000 << " microns"  << endl;
//     cout << "Signal/Noise for X = " << V->GetXSNR() << endl;
  
     delete V;
   
  } 
  
  file->Close();  
}

