#include "iostream.h"

void MUONrawclusters (Int_t evNumber1=0,Int_t evNumber2=0) 
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
    if (!file) file = new TFile("galice.root","UPDATE");

// Get AliRun object from file or create it if not on file

    if (!gAlice) {
	gAlice = (AliRun*)file->Get("gAlice");
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
//
// Set reconstruction models
//
// Get pointers to Alice detectors and Digits containers
    AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");

    RecModel1 = new AliMUONClusterFinder();
    RecModel1->SetNperMax(90);
    //RecModel1->SetClusterSize(12);
    RecModel1->SetClusterSize(100);
    RecModel1->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(0,RecModel1);
    
    RecModel2 = new AliMUONClusterFinder();
    RecModel2->SetNperMax(90);
    //RecModel2->SetClusterSize(12);
    RecModel2->SetClusterSize(100);
    RecModel2->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(1,RecModel2);
 
    RecModel3 = new AliMUONClusterFinder();
    RecModel3->SetNperMax(90);
    //RecModel3->SetClusterSize(12);
    RecModel3->SetClusterSize(100);
    RecModel3->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(2,RecModel3);
    
    RecModel4 = new AliMUONClusterFinder();
    RecModel4->SetNperMax(90);
    //RecModel4->SetClusterSize(12);
    RecModel4->SetClusterSize(100);
    RecModel4->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(3,RecModel4);
    
    //RecModel5 = new AliMUONClusterFinderv0();
    RecModel5 = new AliMUONClusterFinder();
    RecModel5->SetNperMax(90);
    //RecModel5->SetClusterSize(15);
    RecModel5->SetClusterSize(100);
    RecModel5->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(4,RecModel5);
    
    //RecModel6 = new AliMUONClusterFinderv0();
    RecModel6 = new AliMUONClusterFinder();
    RecModel6->SetNperMax(90);
    //RecModel6->SetClusterSize(15);
    RecModel6->SetClusterSize(100);
    RecModel6->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(5,RecModel6);
    
    RecModel7 = new AliMUONClusterFinder();
    RecModel7->SetNperMax(90);
    //RecModel7->SetClusterSize(9);
    RecModel7->SetClusterSize(100);
    RecModel7->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(6,RecModel7);
    
    RecModel8 = new AliMUONClusterFinder();
    RecModel8->SetNperMax(90);
    //RecModel8->SetClusterSize(9);
    RecModel8->SetClusterSize(100);
    RecModel8->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(7,RecModel8);

    RecModel9 = new AliMUONClusterFinder();
    RecModel9->SetNperMax(90);
    //RecModel9->SetClusterSize(9);
    RecModel9->SetClusterSize(100);
    RecModel9->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(8,RecModel9);

    RecModel10 = new AliMUONClusterFinder();
    RecModel10->SetNperMax(90);
    //    RecModel1->SetClusterSize(9);
    RecModel1->SetClusterSize(100);
    RecModel10->SetDeclusterFlag(0);
    MUON->SetReconstructionModel(9,RecModel10);

 
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


	TClonesArray *Particles = gAlice->Particles();
	TTree *TD = gAlice->TreeD();
	Int_t nent=TD->GetEntries();
	//printf("Found %d entries in the tree (must be one per cathode per event!)\n",nent);
	if (MUON) {
	    MUON->FindClusters(nev,nent-2);
	}   // end if MUON
    }   // event loop 
    //file->ls();
    file->Close();
}

