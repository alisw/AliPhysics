#include "iostream.h"

void MUONraw (Int_t evNumber1=0,Int_t evNumber2=0) 
{
  //////////////////////////////////////
  //                                  //
  // ROOT macro for ALICE Dimuon Arm: //
  // Clusterization of digits         //
  //                                  //
  //////////////////////////////////////
  //
  // Adds the tree TR for raw clusters
  // to the ROOT file "galice.root"
  // containing the digits (tree TD).
  //
  // Arguments:
  //   evNumber1 = first event number to act on in file "galice.root"
  //   evNumber2 = last event number to act on in file "galice.root"
  //
  // Input/output file:
  //   "galice.root"
  //
  //__________________________________________________________________________

// Dynamically link some shared libs

    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    }

// Connect the Root Galice file containing Geometry, Kine and Hits

    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (!file) file = new TFile("galice.root","UPDATE");

// Get AliRun object from file or create it if not on file

    if (!gAlice) {
	gAlice = (AliRun*)file->Get("gAlice");
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }

// Set reconstruction models
//
// Get pointers to Alice detectors and Digits containers
    AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
    RecModel = new AliMUONClusterFinderAZ(0,1); //AZ
    //RecModel = new AliMUONClusterFinderAZ(1,0); //AZ
    for (Int_t i=0; i<10; i++) {
      MUON->SetReconstructionModel(i,(AliMUONClusterFinderVS*)RecModel);
    }
//
//   Loop over events              
//
    Int_t Nh=0;
    Int_t Nh1=0;
    gAlice->RunReco("MUON", evNumber1, evNumber2);
}

