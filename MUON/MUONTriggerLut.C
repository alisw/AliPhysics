#include "iostream.h"
//---------------------------------------------------------
// this macro is used to generate the Trigger Look up Table
// and store the TH3S histos in the MUONTriggerLut.root file
//---------------------------------------------------------
void MUONTriggerLut() 
{
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (file) file->Close(); 
   file = new TFile("galice.root","READ");

// Get AliRun object from file or create it if not on file
  printf ("I'm after Map \n");
  if (!gAlice) { 
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }
  printf ("I'm after gAlice \n");
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
    AliMUONTriggerLut* Lut= new AliMUONTriggerLut;
    Lut->LoadLut(); 
}


