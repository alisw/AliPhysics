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
// input directories (recodir and simdir), 
// builds a tree that contains, for each event, an array of muons and dimuons
// (TClonesArrays of AliMUONTrackLight and AliMUONPairLight objects) 
// and writes the tree in an output file (outFileName) 
// Note that if the path for the output file is not explicitly specified, 
// it will be written in the directory containing the generation/reconstruction
// 27-Nov-2006: modified by in order to loop on files

// 13 Nov 2007:
// Updated this macro to work with new version of AliMUONRecoCheck. Also, we are
// now fetching reconstructed track data from ESD and not the track tree, because
// the AliMUONTrack objects for reconstructed tracks are no longer stored on disk.
//  - Artur Szostak <artursz@iafrica.com>

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
#include "AliMUONTrackLight.h"
#include "AliMUONPairLight.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTrackExtrap.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
/*TODO: need to update this with changes made to ITS
#include "AliITSVertexerPPZ.h"
#include "AliITSLoader.h"
*/
#include "AliTracker.h"
#include "AliMagFMaps.h"
#endif


void DecodeRecoCocktail(
    char* recodir=".",          // The directory containing galice.root for reconstructed data.
    char* simdir="generated/",  // The directory containing galice.root for simulated data.
    char* outFileName = "MuonLight.root", // The output filename containing AliMUONTrackLight and AliMUONPairLight objects.
    char* geoFilename = "generated/geometry.root"  // The filename containing the geometry.
  )
{
  char startingDir[200]; 
  sprintf (startingDir,"%s",gSystem->pwd()); 
  gSystem->cd(recodir); 

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
    
  AliESDEvent* esd = new AliESDEvent();
  TTree* treeESD = (TTree*) esdFile->Get("esdTree");
  if (!treeESD) {
    Error("CheckESD", "no ESD tree found");
    return;
  }
  esd->ReadFromTree(treeESD);
  
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

  AliMUONRecoCheck *rc = new AliMUONRecoCheck("galice.root", simdir);
  Int_t nev = rc->NumberOfEvents();
  
  /*TODO: need to update this with changes made to ITS
  AliITSLoader* ITSloader =  (AliITSLoader*) runLoaderSim->GetLoader("ITSLoader");
  AliITSVertexerPPZ *dovert = 0; 
  if (ITSloader) { 
    dovert = new AliITSVertexerPPZ("default",0,0);
    // dovert->SetDebug(0);
    dovert->SetDiffPhiMax(0.05);
    dovert->SetWindow(3.);
  }
  AliESDVertex *vert = 0;
  */
  
  TLorentzVector v; 
 
  for(Int_t ievent = 0; ievent < nev; ievent++){ // loop over events 
    treeESD->GetEvent(ievent);
    rc->GetMCEventHandler()->GetEvent(ievent);
    AliStack* pstack = rc->GetMCEventHandler()->MCEvent()->Stack();
    
    /*TODO: need to update this with changes made to ITS
    runLoaderSim->GetHeader();
    if (ITSloader) { 
      vert = dovert->FindVertexForCurrentEvent(ievent);
    }
    */
    
    //printf ("Event %d of %d\n",ievent,nev);
    muonArray->Clear();     // clean muon and dimuon arrays 
    dimuonArray->Clear(); 
    
    //AliMUONVTrackStore* recoTracks = rc->ReconstructedTracks(ievent);  // Use tracks from actual reconstruction, but this does not work anymore.
    AliMUONVTrackStore* recoTracks = rc->ReconstructedTracks(esd);
    //AliMUONVTrackStore* trackRefs = rc->ReconstructibleTracks(ievent);  // Only use reconstructible reference tracks.
    AliMUONVTrackStore* trackRefs = rc->TrackRefs(ievent);
    
    TIter next(recoTracks->CreateIterator());
    AliMUONTrack* trackReco = NULL;
    
    Int_t nTrackReco = recoTracks->GetSize();
    Int_t nTracksESD = (Int_t)esd->GetNumberOfMuonTracks();
    if (nTrackReco != nTracksESD) printf ("Tracks in recoTracks (%d) and in ESD (%d) do not match!\n", nTrackReco, nTracksESD);
    Int_t nreftracks = 0;
    Int_t itrRec = 0;
    while ( (trackReco = static_cast<AliMUONTrack*>(next())) != NULL )
    {
      // assign parameters concerning the reconstructed tracks
      AliMUONTrackLight muLight;
      
      muLight.FillFromESD(esd->GetMuonTrack(itrRec));
      // muLight.FillFromAliMUONTrack(trackReco);
      
      // find the reference track and store further information	
      TParticle *part = muLight.FindRefTrack(trackReco, trackRefs, pstack, kTRUE);
      if (part) { 
	v.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
	muLight.SetPGen(v); 
	muLight.FillMuonHistory(pstack, part);
	// 	  muLight.PrintInfo("A");
	//store the referenced track in the muonArray:
	TClonesArray &muons = *muonArray;
	new (muons[nreftracks++]) AliMUONTrackLight(muLight);
      } 
      
      itrRec++;
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
    
    delete recoTracks;
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

