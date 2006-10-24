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

// ROOT includes
#include "TClonesArray.h"
#include "TH1.h"
#include "TParticle.h"
#include "TFile.h"

// STEER includes
#include "AliRun.h"
#include "AliHeader.h"
#include "AliMC.h"
#include "AliStack.h"
#include "AliRunLoader.h"
#include "AliMagFMaps.h"

// MUON includes
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONTrackParam.h"
#include "AliTracker.h"

Int_t TrackCheck( Bool_t *compTrack);

void MUONRecoCheck (Int_t nEvent = 1, char * filename="galice.root"){

  // Utility macro to check the muon reconstruction. Reconstructed tracks are compared
  // to reference tracks. The reference tracks are built from AliTrackReference for the
  // hit in chamber (0..9) and from kinematics (TreeK) for the vertex parameters.     
  
  Int_t nTrackReco, nTrackRef;
  AliMUONTrack *trackReco, *trackRef;
  Bool_t *compTrack;
  Bool_t compTrackOK[10];
  Int_t nHitOK = 0;
  Int_t testTrack = 0;	
  Int_t iTrack = 0;
  Int_t indexOK = 0;
  Int_t trackID = 0;
  Double_t sigma2Cut = 16;  // 4 sigmas cut, sigma2Cut = 4*4
  AliMUONTrackParam *trackParam;
  TClonesArray *trackParamAtHit;
  Double_t x1,y1,z1,pX1,pY1,pZ1,p1;
  Double_t x2,y2,z2,pX2,pY2,pZ2,p2;
  TParticle* particle = new TParticle();

  TClonesArray *trackRecoArray = NULL;
  TClonesArray *trackRefArray = NULL;


  // File for histograms and histogram booking
  TFile *histoFile = new TFile("MUONRecoCheck.root", "RECREATE");

  TH1F *hReconstructible = new TH1F("hReconstructible"," Nb of reconstructible tracks ",15,-0.5,14.5);
  TH1F *hReco = new TH1F("hReco"," Nb of reconstructed tracks / evt",15,-0.5,14.5);
  TH1F *hNHitComp = new TH1F("hNHitComp"," Nb of compatible hits / track ",15,-0.5,14.5);
  TH1F *hTestTrack = new TH1F("hTestTrack"," Reconstruction requirement / track",15,-0.5,14.5);
  TH1F *hTrackRefID = new TH1F("hTrackRefID"," track reference ID ",100,-0.5,99.5);
  
  TH1F *hResMomVertex = new TH1F("hMomVertex"," delta P vertex (GeV/c)",100,-10.,10);
  TH1F *hResMomFirstHit = new TH1F("hMomFirstHit"," delta P first hit (GeV/c)",100,-10.,10);

  // set  mag field 
  // waiting for mag field in CDB 
  printf("Loading field map...\n");
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field, kFALSE);

  AliRunLoader* runLoader = AliRunLoader::Open(filename,"read");
  AliLoader * MUONLoader = runLoader->GetLoader("MUONLoader");
  AliMUONData * MUONData = new AliMUONData(MUONLoader,"MUON","MUON");

  runLoader->LoadKinematics("READ");
  runLoader->LoadTrackRefs("READ");
  MUONLoader->LoadTracks("READ");
  
  AliMUONRecoCheck rc(runLoader,MUONData);
    
  Int_t nevents = runLoader->GetNumberOfEvents();
  
  if (nevents < nEvent) nEvent = nevents;
  
  Int_t ievent;
  Int_t nReconstructibleTracks = 0;
  Int_t nReconstructedTracks = 0;
  Int_t nReconstructibleTracksCheck = 0;

  for (ievent=0; ievent<nEvent; ievent++) {
    if (!(ievent%10)) printf(" **** event # %d  \n",ievent);
    runLoader->GetEvent(ievent);
    rc.ResetTracks();
    rc.MakeTrackRef(); // make reconstructible tracks
//     rc.PrintEvent();

    
    trackRecoArray = rc.GetTrackReco();
    trackRefArray = rc.GetMuonTrackRef();
    
    nTrackRef = trackRefArray->GetEntriesFast();
    nTrackReco = trackRecoArray->GetEntriesFast();

    hReconstructible->Fill(rc.GetNumberOfReconstuctibleTracks());
    hReco->Fill(rc.GetNumberOfRecoTracks());

    nReconstructibleTracks += rc.GetNumberOfReconstuctibleTracks();
    nReconstructedTracks += rc.GetNumberOfRecoTracks();

    for (Int_t index1 = 0; index1 < nTrackRef; index1++) {  
      trackRef = (AliMUONTrack *)trackRefArray->At(index1);

      testTrack = 0;
      indexOK = 0;
      for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) 
	compTrackOK[ch] = kFALSE;
      for (Int_t index2 = 0; index2 < nTrackReco; index2++) {
	trackReco = (AliMUONTrack *)trackRecoArray->At(index2);

	// check if trackRef is compatible with trackReco
	compTrack = trackRef->CompatibleTrack(trackReco,sigma2Cut);

	iTrack = TrackCheck(compTrack);
	
	if (iTrack > testTrack) {
	  nHitOK = 0;
	  for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
	    if (compTrack[ch]) nHitOK++;
	    compTrackOK[ch] = compTrack[ch];
	  }
	  testTrack = iTrack;
	  indexOK = index2;
	}
      }

      hTestTrack->Fill(testTrack);
      trackID = trackRef->GetTrackID();
      hTrackRefID->Fill(trackID);

      if (testTrack == 4) {     // tracking requirements verified, track is found
	nReconstructibleTracksCheck++;
	hNHitComp->Fill(nHitOK);
	particle = runLoader->GetHeader()->Stack()->Particle(trackID);
// 	printf(" trackID: %d , PDG code: %d \n",trackID,particle->GetPdgCode());
	trackParam = trackRef->GetTrackParamAtVertex();
	x1 = trackParam->GetNonBendingCoor();
	y1 = trackParam->GetBendingCoor();
	z1 = trackParam->GetZ();
	pX1 = trackParam->Px();
	pY1 = trackParam->Py();
	pZ1 = trackParam->Pz();
	p1  = trackParam->P();
	
// 	printf(" Ref. track at vertex: x,y,z: %f %f %f px,py,pz,p: %f %f %f %f \n",x1,y1,z1,pX1,pY1,pZ1,p1);
	
	trackParam = ((AliMUONTrack *)trackRecoArray->At(indexOK))->GetTrackParamAtVertex();
	x2 = trackParam->GetNonBendingCoor();
	y2 = trackParam->GetBendingCoor();
	z2 = trackParam->GetZ();
	pX2 = trackParam->Px();
	pY2 = trackParam->Py();
	pZ2 = trackParam->Pz();
	p2  = trackParam->P();
// 	printf(" Reconst. track at vertex: x,y,z: %f %f %f px,py,pz: %f %f %f %f \n",x2,y2,z2,pX2,pY2,pZ2,p2);
	
	hResMomVertex->Fill(p2-p1);

 	trackParamAtHit =  trackRef->GetTrackParamAtHit();
	trackParam = (AliMUONTrackParam*) trackParamAtHit->First();
	x1 = trackParam->GetNonBendingCoor();
	y1 = trackParam->GetBendingCoor();
	z1 = trackParam->GetZ();
	pX1 = trackParam->Px();
	pY1 = trackParam->Py();
	pZ1 = trackParam->Pz();
	p1  = trackParam->P();
// 	printf(" Ref. track at 1st hit: x,y,z: %f %f %f px,py,pz: %f %f %f \n",x1,y1,z1,pX1,pY1,pZ1);
	trackParamAtHit =  ((AliMUONTrack *) trackRecoArray->At(indexOK))->GetTrackParamAtHit();
	trackParam = (AliMUONTrackParam*) trackParamAtHit->First();
	x2 = trackParam->GetNonBendingCoor();
	y2 = trackParam->GetBendingCoor();
	z2 = trackParam->GetZ();
	pX2 = trackParam->Px();
	pY2 = trackParam->Py();
	pZ2 = trackParam->Pz();
	p2  = trackParam->P();
// 	printf(" Reconst. track at 1st hit: x,y,z: %f %f %f px,py,pz: %f %f %f \n",x2,y2,z2,pX2,pY2,pZ2);

	hResMomFirstHit->Fill(p2-p1);
	       
      }
    } // end loop track ref.

  } // end loop on event  

  MUONLoader->UnloadTracks();
  runLoader->UnloadKinematics();
  runLoader->UnloadTrackRefs();
  runLoader->Delete();
  field->Delete();

  printf(" nb of reconstructible tracks: %d \n", nReconstructibleTracks);
  printf(" nb of reconstructed tracks: %d \n", nReconstructedTracks);
  printf(" nb of reconstructible tracks which are reconstructed: %d \n", nReconstructibleTracksCheck);
  
  histoFile->Write();
  histoFile->Close();
}


Int_t TrackCheck( Bool_t *compTrack)
{
  // Apply reconstruction requirements
  // Return number of validated conditions 
  // If all the tests are verified then TrackCheck = 4 (good track)
  Int_t iTrack = 0;
  Int_t hitsInLastStations = 0;
  
  // apply reconstruction requirements
  if (compTrack[0] || compTrack[1]) iTrack++; // at least one hit in st. 0
  if (compTrack[2] || compTrack[3]) iTrack++; // at least one hit in st. 1
  if (compTrack[4] || compTrack[5]) iTrack++; // at least one hit in st. 2
  for (Int_t ch = 6; ch < AliMUONConstants::NTrackingCh(); ch++) {
    if (compTrack[ch]) hitsInLastStations++; 
  }
  if (hitsInLastStations > 2) iTrack++; // at least 3 hits in st. 3 & 4
  
  return iTrack;

}







