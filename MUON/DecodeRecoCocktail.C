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
// 27-Nov-2006: modified by in order to loop on files

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TSystem.h>
#include <TGeoManager.h>
#include "AliMUONRecoCheck.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliRunLoader.h"
#include "AliMUONTrackLight.h"
#include "AliMUONPairLight.h"
#include "AliMUONTrackExtrap.h"
#include "AliESD.h"
#include "AliESDVertex.h" 
#include "AliITSVertexerPPZ.h"
#include "AliITSLoader.h"
#include "AliTracker.h"
#include "AliMagFMaps.h"
#endif


void DecodeRecoCocktail(char* dirname=".", char* outFileName = "MuonLight.root", char* geoFilename = "geometry.root"){

  
  char startingDir[200]; 
  sprintf (startingDir,"%s",gSystem->pwd()); 
  gSystem->cd(dirname); 


  TClonesArray *muonArray   = new TClonesArray("AliMUONTrackLight",10); 
  TClonesArray *dimuonArray = new TClonesArray("AliMUONPairLight",10); 
  TTree *treeOut = new TTree("tree","tree"); 

  TFile *fout = new TFile(outFileName,"recreate"); 
  fout->cd(); 

  treeOut->Branch("muons",&muonArray); 
  treeOut->Branch("dimuons",&dimuonArray); 
  
  TFile* esdFile = TFile::Open("AliESDs.root");
  if (!esdFile || !esdFile->IsOpen()) {
    Error("DecodeRecoCocktailNew", "opening ESD file AliESDs.root failed");
    return;
  }
    
  AliESD* esd = new AliESD();
  TTree* treeESD = (TTree*) esdFile->Get("esdTree");
  if (!treeESD) {
    Error("CheckESD", "no ESD tree found");
    return;
  }
  treeESD->SetBranchAddress("ESD", &esd);
   
  // Import TGeo geometry (needed by AliMUONTrackExtrap::ExtrapToVertex)
  if (!gGeoManager) {
    TGeoManager::Import(geoFilename);
    if (!gGeoManager) {
      Error("MUONmass_ESD", "getting geometry from file %s failed", geoFilename);
      return;
    }
  }
  
  // set  mag field 
  // waiting for mag field in CDB 
  printf("Loading field map...\n");
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field, kFALSE);
  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField(AliTracker::GetFieldMap());

  AliMUONRecoCheck *rc = new AliMUONRecoCheck("galice.root", "galice_sim.root");
  AliRunLoader *runLoader = rc->GetRunLoader();
  AliRunLoader *runLoaderSim = rc->GetRunLoaderSim();
  
  runLoaderSim->LoadKinematics("READ");
  Int_t nev = runLoaderSim->GetNumberOfEvents(); 
//   nevent = nevent +nev;  
//     printf(" number of files and events = %d - %d \n",irun,nevent);
    
  AliITSLoader* ITSloader =  (AliITSLoader*) runLoaderSim->GetLoader("ITSLoader");
  AliITSVertexerPPZ *dovert = 0; 
  if (ITSloader) { 
    dovert = new AliITSVertexerPPZ("default",0,0);
    // dovert->SetDebug(0);
    dovert->SetDiffPhiMax(0.05);
    dovert->SetWindow(3.);
  }
  AliESDVertex *vert = 0;
  
  TLorentzVector v; 
  
  for(Int_t ievent = 0; ievent < nev; ievent++){ // loop over events 
    runLoaderSim->GetHeader();
    if (ITSloader) { 
      vert = dovert->FindVertexForCurrentEvent(ievent);
    }
    //printf ("Event %d of %d\n",ievent,nev);
    muonArray->Clear();     // clean muon and dimuon arrays 
    dimuonArray->Clear(); 
    runLoader->GetEvent(ievent);
    runLoaderSim->GetEvent(ievent);
    treeESD->GetEvent(ievent);
    rc->ResetTracks();
    rc->MakeTrackRef(); // make reconstructable tracks
    
    TClonesArray * trackRecoArray = rc->GetTrackReco();
    TClonesArray * trackRefArray = rc->GetMuonTrackRef();
    Int_t nTrackReco = trackRecoArray->GetEntriesFast();
    Int_t nTracksESD = (Int_t)esd->GetNumberOfMuonTracks() ; 
    if (nTrackReco != nTracksESD) printf ("Tracks in AliMUONTrack and in ESD do not match!\n");
    Int_t nreftracks = 0; 
    for (Int_t itrRec = 0; itrRec<nTrackReco; itrRec++) { 
      // assign parameters concerning the reconstructed tracks
      AliMUONTrackLight muLight;
      AliMUONTrack *trackReco = (AliMUONTrack *)trackRecoArray->At(itrRec);
      
      muLight.FillFromESD(esd->GetMuonTrack(itrRec));
      // 	muLight.FillFromAliMUONTrack(trackReco);
      
      // find the reference track and store further information	
      TParticle *part = muLight.FindRefTrack(trackReco,trackRefArray,runLoaderSim); 
      if (part) { 
	v.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
	muLight.SetPGen(v); 
	muLight.FillMuonHistory(runLoaderSim, part);
	// 	  muLight.PrintInfo("A");
	//store the referenced track in the muonArray:
	TClonesArray &muons = *muonArray;
	new (muons[nreftracks++]) AliMUONTrackLight(muLight);
      } 
    }  // end reco track
    
    // now loop over muon pairs to build dimuons
    Int_t nmuons = muonArray->GetEntriesFast(); 
    Int_t ndimuons = 0; 
    for(Int_t itrRec1 = 0; itrRec1 < nmuons-1; itrRec1++) { 
      AliMUONTrackLight* mu1 = (AliMUONTrackLight*) muonArray->At(itrRec1); 
      for(Int_t itrRec2 = itrRec1+1; itrRec2 < nmuons; itrRec2++){
	AliMUONTrackLight* mu2 = (AliMUONTrackLight*) muonArray->At(itrRec2); 
	AliMUONPairLight dimuLight;
	dimuLight.SetMuons(*mu1, *mu2);
	TClonesArray &dimuons = *dimuonArray;
	new (dimuons[ndimuons++]) AliMUONPairLight(dimuLight);
      }
    }
    treeOut->Fill(); 
  } 
  
  fout->cd(); 
  treeOut->Write(); 
  gSystem->cd(startingDir); 
  fout->Close();
  delete muonArray;
  delete dimuonArray;
  delete treeOut; 
  delete fout;
  delete rc;
  delete esd; 
}

