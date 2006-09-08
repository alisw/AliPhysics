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

// A. De Falco, H. Woehri, INFN Cagliari, July 2006
// This macro reads the generation/reconstruction files in the 
// input directory (dirname), 
// builds a tree that contains, for each event, an array of muons and dimuons
// (TClonesArrays of AliMUONTrackLight and AliMUONPairLight objects) 
// and writes the tree in an output file (outFileName) 
// Note that if the path for the output file is not explicitly specified, 
// it will be written in the directory containing the generation/reconstruction

#include <iostream>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TSystem.h>
#include <AliMUONRecoCheck.h>
#include <AliMUONTrack.h>
#include <AliMUONTrackParam.h>
#include <AliRunLoader.h>
#include "AliMUONTrackLight.h"
#include "AliMUONPairLight.h"

void DecodeRecoCocktail(char* dirname=".", char* outFileName = "MuonLight.root"){ 
  const char *startingDir = gSystem->pwd(); 
  gSystem->cd(dirname); 
  TFile *fout = new TFile(outFileName,"recreate"); 

  AliMUONRecoCheck *rc = new AliMUONRecoCheck("galice.root");
  AliRunLoader *runLoader = rc->GetRunLoader();
  
  TClonesArray *muonArray   = new TClonesArray("AliMUONTrackLight",10); 
  TClonesArray *dimuonArray = new TClonesArray("AliMUONPairLight",10); 
  TTree *treeOut = new TTree("tree","tree"); 

  fout->cd(); 
  treeOut->Branch("muons",&muonArray); 
  treeOut->Branch("dimuons",&dimuonArray); 
  
  runLoader->LoadKinematics("READ");;
  Int_t nev = runLoader->GetNumberOfEvents(); 
  
  TLorentzVector v; 
  for(Int_t ievent = 0; ievent < nev; ievent++){ // loop over events 
    printf ("Event %d of %d\n",ievent,nev);
    muonArray->Clear();     // clean muon and dimuon arrays 
    dimuonArray->Clear(); 
    runLoader->GetEvent(ievent);
    rc->ResetTracks();
    rc->MakeTrackRef(); // make reconstructable tracks
    
    TClonesArray * trackRecoArray = rc->GetTrackReco();
    TClonesArray * trackRefArray = rc->GetMuonTrackRef();
    Int_t nTrackReco = trackRecoArray->GetEntriesFast();
    Int_t nreftracks = 0; 
    for (Int_t itrRec = 0; itrRec<nTrackReco; itrRec++) { 
      // assign parameters concerning the reconstructed tracks
      AliMUONTrack *trackReco = (AliMUONTrack *)trackRecoArray->At(itrRec);
      AliMUONTrackLight *muLight = new AliMUONTrackLight(); 
      AliMUONTrackParam *trPar = trackReco->GetTrackParamAtVertex(); 
      muLight->SetCharge(Int_t(TMath::Sign(1.,trPar->GetInverseBendingMomentum())));
      muLight->SetPxPyPz(trPar->Px(),trPar->Py(), trPar->Pz()); 
      muLight->SetTriggered(trackReco->GetMatchTrigger()); 
      Double_t xyz[3] = { trPar->GetNonBendingCoor(), 
			  trPar->GetBendingCoor(), 
			  trPar->GetZ()};
      muLight->SetVertex(xyz); 
      // find the reference track and store further information
      TParticle *part = muLight->FindRefTrack(trackReco,trackRefArray,runLoader); 
      if (part) { 
	v.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
	muLight->SetPGen(v); 
	muLight->FillMuonHistory(runLoader, part);
	muLight->PrintInfo("A");
	//store the referenced track in the muonArray:
	TClonesArray &muons = *muonArray;
	new (muons[nreftracks++]) AliMUONTrackLight(*muLight);
      }
    }
    
    // now loop over muon pairs to build dimuons
    Int_t nmuons = muonArray->GetEntriesFast(); 
    Int_t ndimuons = 0; 
    for(Int_t itrRec1 = 0; itrRec1 < nmuons-1; itrRec1++) { 
      AliMUONTrackLight* mu1 = (AliMUONTrackLight*) muonArray->At(itrRec1); 
      for(Int_t itrRec2 = itrRec1+1; itrRec2 < nmuons; itrRec2++){
	AliMUONTrackLight* mu2 = (AliMUONTrackLight*) muonArray->At(itrRec2); 
	AliMUONPairLight *dimuLight = new AliMUONPairLight(); 
	dimuLight->SetMuons(*mu1, *mu2);
	//     	if(dimuLight->GetCreationProcess() == 2){
	// 	  printf ("#dimuon = %d (%d, %d) \n", ndimuons, itrRec1, itrRec2);
	// dimuLight->PrintInfo("A");
	//   	}
	//store the referenced track in the dimuonArray:
	TClonesArray &dimuons = *dimuonArray;
	new (dimuons[ndimuons++]) AliMUONPairLight(*dimuLight);
      }
    }
    Int_t ndimu2 = dimuonArray->GetEntriesFast(); 
    //     printf ("dimuonArray has %d entries\n",ndimu2); 
    //    dimuonArray->Dump(); 
    //     printf ("filling tree\n"); 
    //     // fill the tree
    treeOut->Fill(); 
    //     printf ("done\n"); 
    
  }
  fout->cd(); 
  treeOut->Write(); 
  gSystem->cd(startingDir); 
}
