#include "iostream.h"

void MUONrawclusters (Int_t evNumber1=0,Int_t evNumber2=0) 
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
    for (Int_t i=0; i<10; i++) {
 	AliMUONChamber* iChamber= &(MUON->Chamber(i));
	AliMUONResponse* response =  iChamber->ResponseModel();
	AliMUONSegmentation*  seg1 = iChamber->SegmentationModel(1);
	AliMUONSegmentation*  seg2 = iChamber->SegmentationModel(2);
//
	RecModel = new AliMUONClusterFinderVS();
	RecModel->SetNperMax(90);
	RecModel->SetClusterSize(100);
	RecModel->SetDeclusterFlag(0);
	RecModel->SetSegmentation(seg1,seg2);
	RecModel->SetResponse(response); 
//	RecModel->SetTracks(16,17);    
//	RecModel->SetTracks(266,267);    
	MUON->SetReconstructionModel(i,RecModel);
    }
//
//   Loop over events              
//
    Int_t Nh=0;
    Int_t Nh1=0;
    for (int nev=evNumber1; nev<= evNumber2; nev++) {
	Int_t nparticles = gAlice->GetEvent(nev);
	cout << "nev         " << nev <<endl;
	cout << "nparticles  " << nparticles <<endl;
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;
	Int_t nbytes = 0;
	TClonesArray *Particles = gAlice->Particles();
	TTree *TD = gAlice->TreeD();
	Int_t nent=TD->GetEntries();
	if (MUON) {
	    MUON->FindClusters(nev,nent-2);
	}   // end if MUON
    }   // event loop 

    file->Close();
}

