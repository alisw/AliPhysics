/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "iostream.h"

#include <TClassTable.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TTree.h>

#include "AliHeader.h"
#include "AliRun.h"

#include "AliMUON.h"

#include "AliMUONClusterFinderVS.h"
#endif

void MUONrawclusters (char* filename="galice.root", Int_t evNumber1=0,Int_t evNumber2=9999) 
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

//  // Dynamically link some shared libs

  //if (gClassTable->GetID("AliRun") < 0) {
  //	gROOT->LoadMacro("$ALICE_ROOT/macros/loadlibs.C");
  //	loadlibs();
  //    }

 // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","UPDATE");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }

  // Loading AliRun master
  RunLoader->UnloadgAlice();
  RunLoader->LoadgAlice();
  gAlice = RunLoader->GetAliRun();

  // Loading MUON subsystem
  AliMUON * MUON = (AliMUON *) gAlice->GetDetector("MUON");
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  AliMUONData * muondata = MUON->GetMUONData();

  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();

  for (Int_t i=0; i<10; i++) {
    //RecModel = new AliMUONClusterFinderVS();
    AliMUONClusterFinderVS *RecModel = new AliMUONClusterFinderVS();
    //	RecModel->SetTracks(16,17);    
    //	RecModel->SetTracks(266,267);    
    RecModel->SetGhostChi2Cut(10);
    MUON->SetReconstructionModel(i,RecModel);
  } 

  MUONLoader->LoadDigits("READ");
  MUONLoader->LoadRecPoints("UPDATE");

  // Testing if RawClusterisation has already been done
  RunLoader->GetEvent(0);
  if (MUONLoader->TreeR()) {
    if (muondata->IsRawClusterBranchesInTree()) {
      MUONLoader->UnloadRecPoints();
      MUONLoader->LoadRecPoints("RECREATE");
      printf("Recreating RecPoints files\n");
    }
  }
//
//   Loop over events              
  //
    Int_t Nh=0;
    Int_t Nh1=0;
    //    gAlice->RunReco("MUON", evNumber1, evNumber2);
    if (evNumber2>nevents) evNumber2=nevents;
    for(Int_t ievent=evNumber1; ievent<evNumber2; ievent++) {
      printf("event %d\n",ievent);
      RunLoader->GetEvent(ievent);

      // Test if rawcluster has already been done before
      if (MUONLoader->TreeR() == 0x0) 
	MUONLoader->MakeRecPointsContainer();
      else {
	if (muondata->IsRawClusterBranchesInTree()){ // Test if rawcluster has already been done before
	  if (ievent==evNumber1) MUONLoader->UnloadRecPoints();
	  MUONLoader->MakeRecPointsContainer();  // Redoing clusterisation
	  Info("RecPointsContainer","Recreating RecPointsContainer and deleting previous ones");
	}
      }
      muondata->MakeBranch("RC");
      muondata->SetTreeAddress("D,RC");
      MUON->Digits2Reco(); 
    }
    MUONLoader->UnloadDigits();
    MUONLoader->UnloadRecPoints();
}


