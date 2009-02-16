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

//-----------------------------------------------------------------------------
/// \class AliMUONRecoCheck
/// Utility class to check reconstruction
/// Reconstructed tracks are compared to reference tracks. 
/// The reference tracks are built from AliTrackReference for the
/// hit in chamber (0..9) and from kinematics for the vertex parameters.     
//-----------------------------------------------------------------------------

// 13 Nov 2007:
// Added a method to create a list of reconstructed AliMUONTrack objects from
// ESD data. This is necessary since the track objects that are actually created
// during offline reconstruction are no longer stored to disk.
//  - Artur Szostak <artursz@iafrica.com>
// 25 Jan 2008:
// Use the new ESDInterface to create MUON objects from ESD data
// - Philippe Pillot
// 

#include "AliMUONRecoCheck.h"
#include "AliMUONTrack.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONConstants.h"
#include "AliMUONESDInterface.h"
#include "AliMUONTrackParam.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONRecoCheck)
/// \endcond

//_____________________________________________________________________________
AliMUONRecoCheck::AliMUONRecoCheck(const Char_t *esdFileName, const Char_t *pathSim)
: TObject(),
fMCEventHandler(new AliMCEventHandler()),
fESDEvent(new AliESDEvent()),
fESDTree (0x0),
fESDFile (0x0),
fCurrentEvent(0),
fTrackRefStore(0x0),
fRecoTrackRefStore(0x0),
fRecoTrackStore(0x0),
fESDEventOwner(kTRUE)
{
  /// Normal ctor
  
  // TrackRefs and Particules
  fMCEventHandler->SetInputPath(pathSim);
  fMCEventHandler->InitIO("");
  
  // ESD MUON Tracks
  fESDFile = TFile::Open(esdFileName); // open the file
  if (!fESDFile || !fESDFile->IsOpen()) {
    AliError(Form("opening ESD file %s failed", esdFileName));
    fESDFile = 0x0;
    return;
  }
  fESDTree = (TTree*) fESDFile->Get("esdTree"); // get the tree
  if (!fESDTree) {
    AliError("no ESD tree found");
    fESDFile->Close();
    fESDFile = 0x0;
    return;
  }
  fESDEvent->ReadFromTree(fESDTree); // link fESDEvent to the tree
}

//_____________________________________________________________________________
AliMUONRecoCheck::AliMUONRecoCheck(AliESDEvent *esdEvent, AliMCEventHandler *mcEventHandler)
: TObject(),
fMCEventHandler(0),
fESDEvent(0),
fESDTree (0x0),
fESDFile (0x0),
fCurrentEvent(0),
fTrackRefStore(0x0),
fRecoTrackRefStore(0x0),
fRecoTrackStore(0x0),
fESDEventOwner(kFALSE)
{
  /// Normal ctor
  
  // TrackRefs and Particules
  fMCEventHandler = mcEventHandler;
  
  // ESD MUON Tracks
  fESDEvent = esdEvent;
  
}

//_____________________________________________________________________________
AliMUONRecoCheck::~AliMUONRecoCheck()
{
  /// Destructor
  if (fESDEventOwner) {
    delete fMCEventHandler;
    delete fESDEvent;
    if (fESDFile) fESDFile->Close();
  }
  ResetStores();
}

//_____________________________________________________________________________
void AliMUONRecoCheck::ResetStores()
{
  /// Deletes all the store objects that have been created and resets the pointers to 0x0
  delete fTrackRefStore;      fTrackRefStore = 0x0;
  delete fRecoTrackRefStore;  fRecoTrackRefStore = 0x0;
  delete fRecoTrackStore;     fRecoTrackStore = 0x0;
}

//_____________________________________________________________________________
Int_t AliMUONRecoCheck::NumberOfEvents() const
{
  /// Return the number of events
  if (fESDEventOwner && fESDTree) return fESDTree->GetEntries();
  return 0;
}

//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONRecoCheck::ReconstructedTracks(Int_t event)
{
  /// Return a track store containing the reconstructed tracks (converted into 
  /// MUONTrack objects) for a given event
  
  if (!fESDEventOwner) {
    MakeReconstructedTracks();
    return fRecoTrackStore;
  }

  if (event != fCurrentEvent) {
    ResetStores();
    fCurrentEvent = event;
  }
  
  if (fRecoTrackStore != 0x0) return fRecoTrackStore;
  else {
    if (!fESDTree) return 0x0;
    if (fESDTree->GetEvent(event) <= 0) {
      AliError(Form("fails to read ESD object for event %d", event));
      return 0x0;
    }
    MakeReconstructedTracks();
    return fRecoTrackStore;
  }
}

//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONRecoCheck::TrackRefs(Int_t event)
{
  /// Return a track store containing the track references (converted into 
  /// MUONTrack objects) for a given event
  
  if (!fESDEventOwner) {
    MakeTrackRefs();
    return fTrackRefStore;
  }

  if (event != fCurrentEvent) {
    ResetStores();
    fCurrentEvent = event;
  }
  
  if (fTrackRefStore != 0x0) return fTrackRefStore;
  else {
    if (!fMCEventHandler->GetEvent(event)) {
      AliError(Form("fails to read MC objects for event %d", event));
      return 0x0;
    }
    MakeTrackRefs();
    return fTrackRefStore;
  }
}

//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONRecoCheck::ReconstructibleTracks(Int_t event)
{
  /// Return a track store containing the reconstructible tracks for a given event

  if (!fESDEventOwner) {
    if (TrackRefs(event) == 0x0) return 0x0;
    MakeReconstructibleTracks();
    return fRecoTrackRefStore;
  }

  if (event != fCurrentEvent) {
    ResetStores();
    fCurrentEvent = event;
  }
  
  if (fRecoTrackRefStore != 0x0) return fRecoTrackRefStore;
  else {
    if (TrackRefs(event) == 0x0) return 0x0;
    MakeReconstructibleTracks();
    return fRecoTrackRefStore;
  }
}

//_____________________________________________________________________________
void AliMUONRecoCheck::MakeReconstructedTracks()
{
  /// Make reconstructed tracks
  if (!(fRecoTrackStore = AliMUONESDInterface::NewTrackStore())) return;
  
  // loop over all reconstructed tracks and add them to the store (skip ghosts)
  Int_t nTracks = (Int_t) fESDEvent->GetNumberOfMuonTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliESDMuonTrack* esdTrack = fESDEvent->GetMuonTrack(iTrack);
    if (esdTrack->ContainTrackerData()) AliMUONESDInterface::Add(*esdTrack, *fRecoTrackStore);
  }
  
}

//_____________________________________________________________________________
void AliMUONRecoCheck::MakeTrackRefs()
{
  /// Make reconstructible tracks
  AliMUONVTrackStore *tmpTrackRefStore = AliMUONESDInterface::NewTrackStore();
  if (!tmpTrackRefStore) return;
  
  Double_t x, y, z, pX, pY, pZ, bendingSlope, nonBendingSlope, inverseBendingMomentum;
  TParticle* particle;
  TClonesArray* trackRefs;
  Int_t nTrackRef = fMCEventHandler->MCEvent()->GetNumberOfTracks();
  AliMUONVClusterStore* cStore = AliMUONESDInterface::NewClusterStore();
  if (!cStore) return;
  AliMUONVCluster* hit = cStore->CreateCluster(0,0,0);
  
  // loop over simulated tracks
  for (Int_t iTrackRef  = 0; iTrackRef < nTrackRef; ++iTrackRef) {
    Int_t nHits = fMCEventHandler->GetParticleAndTR(iTrackRef, particle, trackRefs);
    
    // skip empty trackRefs
    if (nHits < 1) continue;
    
    // get the particle charge for further calculation
    TParticlePDG* ppdg = particle->GetPDG();
    Int_t charge = ppdg != NULL ? (Int_t)(ppdg->Charge()/3.0) : 0;
    
    AliMUONTrack track;
    
    // loop over simulated track hits
    for (Int_t iHit = 0; iHit < nHits; ++iHit) {        
      AliTrackReference* trackReference = static_cast<AliTrackReference*>(trackRefs->UncheckedAt(iHit));
      
      // skip trackRefs not in MUON
      if (trackReference->DetectorId() != AliTrackReference::kMUON) continue;
    
      // Get track parameters of current hit
      x = trackReference->X();
      y = trackReference->Y();
      z = trackReference->Z();
      pX = trackReference->Px();
      pY = trackReference->Py();
      pZ = trackReference->Pz();
      
      // check chamberId of current trackReference
      Int_t detElemId = trackReference->UserId();
      Int_t chamberId = detElemId / 100 - 1;
      if (chamberId < 0 || chamberId >= AliMUONConstants::NTrackingCh()) continue;
      
      // set hit parameters
      hit->SetUniqueID(AliMUONVCluster::BuildUniqueID(chamberId, detElemId, iHit));
      hit->SetXYZ(x,y,z);
      hit->SetErrXY(0.,0.);
      
      // compute track parameters at hit
      bendingSlope = 0;
      nonBendingSlope = 0;
      inverseBendingMomentum = 0;
      if (TMath::Abs(pZ) > 0) {
	bendingSlope = pY/pZ;
	nonBendingSlope = pX/pZ;
      }
      Double_t pYZ = TMath::Sqrt(pY*pY+pZ*pZ);
      if (pYZ >0) inverseBendingMomentum = 1/pYZ; 
      inverseBendingMomentum *= charge;
      
      // set track parameters at hit
      AliMUONTrackParam trackParam;
      trackParam.SetNonBendingCoor(x);
      trackParam.SetBendingCoor(y);
      trackParam.SetZ(z);
      trackParam.SetBendingSlope(bendingSlope);
      trackParam.SetNonBendingSlope(nonBendingSlope);
      trackParam.SetInverseBendingMomentum(inverseBendingMomentum);
      
      // add track parameters at current hit to the track
      track.AddTrackParamAtCluster(trackParam, *hit, kTRUE);
    }
    
    // if none of the track hits was in MUON, goto the next track
    if (track.GetNClusters() < 1) continue;
    
    // get track parameters at particle's vertex
    x = particle->Vx();
    y = particle->Vy();
    z = particle->Vz();
    pX = particle->Px();
    pY = particle->Py();
    pZ = particle->Pz();
    
    // compute rest of track parameters at particle's vertex
    bendingSlope = 0;
    nonBendingSlope = 0;
    inverseBendingMomentum = 0;
    if (TMath::Abs(pZ) > 0) {
      bendingSlope = pY/pZ;
      nonBendingSlope = pX/pZ;
    }
    Double_t pYZ = TMath::Sqrt(pY*pY+pZ*pZ);
    if (pYZ >0) inverseBendingMomentum = 1/pYZ;
    inverseBendingMomentum *= charge;
    
    // set track parameters at particle's vertex
    AliMUONTrackParam trackParamAtVertex;
    trackParamAtVertex.SetNonBendingCoor(x);
    trackParamAtVertex.SetBendingCoor(y);
    trackParamAtVertex.SetZ(z);
    trackParamAtVertex.SetBendingSlope(bendingSlope);
    trackParamAtVertex.SetNonBendingSlope(nonBendingSlope);
    trackParamAtVertex.SetInverseBendingMomentum(inverseBendingMomentum);
    
    // add track parameters at vertex
    track.SetTrackParamAtVertex(&trackParamAtVertex);
    
    // store the track
    track.SetTrackID(iTrackRef);
    tmpTrackRefStore->Add(track);
  }
  
  CleanMuonTrackRef(tmpTrackRefStore);
  
  delete hit;
  delete cStore;
  delete tmpTrackRefStore;
}

//_____________________________________________________________________________
void AliMUONRecoCheck::CleanMuonTrackRef(const AliMUONVTrackStore *tmpTrackRefStore)
{
  /// Re-calculate hits parameters because two AliTrackReferences are recorded for
  /// each chamber (one when particle is entering + one when particle is leaving 
  /// the sensitive volume) 
  if (!(fTrackRefStore = AliMUONESDInterface::NewTrackStore())) return;
  
  Double_t maxGasGap = 1.; // cm 
  Double_t x, y, z, pX, pY, pZ, x1, y1, z1, pX1, pY1, pZ1, z2;
  Double_t bendingSlope,nonBendingSlope,inverseBendingMomentum;
  AliMUONVClusterStore* cStore = AliMUONESDInterface::NewClusterStore();
  if (!cStore) return;
  AliMUONVCluster* hit = cStore->CreateCluster(0,0,0);
  
  // create iterator
  TIter next(tmpTrackRefStore->CreateIterator());
  
  // loop over tmpTrackRef
  AliMUONTrack* track;
  while ( ( track = static_cast<AliMUONTrack*>(next()) ) ) {
    
    AliMUONTrack newTrack;
    
    // loop over tmpTrackRef's hits
    Int_t iHit1 = 0;
    Int_t nTrackHits = track->GetNClusters();
    while (iHit1 < nTrackHits) {
      AliMUONTrackParam *trackParam1 = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->UncheckedAt(iHit1);
      
      // get track parameters at hit1
      x1  = trackParam1->GetNonBendingCoor();
      y1  = trackParam1->GetBendingCoor();
      z1  = trackParam1->GetZ();
      pX1 = trackParam1->Px();
      pY1 = trackParam1->Py();
      pZ1 = trackParam1->Pz();
      
      // prepare new track parameters
      x  = x1;
      y  = y1;
      z  = z1;
      pX = pX1;
      pY = pY1;
      pZ = pZ1;
      
      // loop over next tmpTrackRef's hits
      Int_t nCombinedHits = 1;
      for (Int_t iHit2 = iHit1+1; iHit2 < nTrackHits; iHit2++) {
	AliMUONTrackParam *trackParam2 = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->UncheckedAt(iHit2);
        
	// get z position of hit2
	z2 = trackParam2->GetZ();
	
	// complete new track parameters if hit2 is on the same detection element
        if ( TMath::Abs(z2-z1) < maxGasGap ) {
	  x  += trackParam2->GetNonBendingCoor();
	  y  += trackParam2->GetBendingCoor();
	  z  += z2;
	  pX += trackParam2->Px();
	  pY += trackParam2->Py();
	  pZ += trackParam2->Pz();
	  nCombinedHits++;
          iHit1 = iHit2;
        }
        
      }
      
      // finalize new track parameters
      x  /= (Double_t)nCombinedHits;
      y  /= (Double_t)nCombinedHits;
      z  /= (Double_t)nCombinedHits;
      pX /= (Double_t)nCombinedHits;
      pY /= (Double_t)nCombinedHits;
      pZ /= (Double_t)nCombinedHits;
      bendingSlope = 0;
      nonBendingSlope = 0;
      inverseBendingMomentum = 0;
      if (TMath::Abs(pZ) > 0) {
	bendingSlope = pY/pZ;
	nonBendingSlope = pX/pZ;
      }
      Double_t pYZ = TMath::Sqrt(pY*pY+pZ*pZ);
      if (pYZ >0) inverseBendingMomentum = 1/pYZ; 
      inverseBendingMomentum *= trackParam1->GetCharge();
      
      // set hit parameters
      hit->SetUniqueID(trackParam1->GetClusterPtr()->GetUniqueID());
      hit->SetXYZ(x,y,z);
      hit->SetErrXY(0.,0.);
      
      // set new track parameters at new hit
      AliMUONTrackParam trackParam;
      trackParam.SetNonBendingCoor(x);
      trackParam.SetBendingCoor(y);
      trackParam.SetZ(z);
      trackParam.SetBendingSlope(bendingSlope);
      trackParam.SetNonBendingSlope(nonBendingSlope);
      trackParam.SetInverseBendingMomentum(inverseBendingMomentum);
      
      // add track parameters at current hit to the track
      newTrack.AddTrackParamAtCluster(trackParam, *hit, kTRUE);
      
      iHit1++;
    }
    
    newTrack.SetTrackID(track->GetTrackID());
    newTrack.SetTrackParamAtVertex(track->GetTrackParamAtVertex());
    fTrackRefStore->Add(newTrack);
    
  }
  
  delete hit;
  delete cStore;
}

//_____________________________________________________________________________
void AliMUONRecoCheck::MakeReconstructibleTracks()
{
  /// Isolate the reconstructible tracks
  if (!(fRecoTrackRefStore = AliMUONESDInterface::NewTrackStore())) return;
  
  // create iterator on trackRef
  TIter next(fTrackRefStore->CreateIterator());
  
  // loop over trackRef
  AliMUONTrack* track;
  while ( ( track = static_cast<AliMUONTrack*>(next()) ) ) {
    
    Bool_t* chamberInTrack = new Bool_t(AliMUONConstants::NTrackingCh());
    for (Int_t iCh = 0; iCh < AliMUONConstants::NTrackingCh(); iCh++) chamberInTrack[iCh] = kFALSE;
    
    // loop over trackRef's hits to get hit chambers
    Int_t nTrackHits = track->GetNClusters();
    for (Int_t iHit = 0; iHit < nTrackHits; iHit++) {
      AliMUONVCluster* hit = ((AliMUONTrackParam*) track->GetTrackParamAtCluster()->UncheckedAt(iHit))->GetClusterPtr(); 
      chamberInTrack[hit->GetChamberId()] = kTRUE;
    } 
    
    // track is reconstructible if the particle is depositing a hit
    // in the following chamber combinations:
    Bool_t trackOK = kTRUE;
    if (!chamberInTrack[0] && !chamberInTrack[1]) trackOK = kFALSE;
    if (!chamberInTrack[2] && !chamberInTrack[3]) trackOK = kFALSE;
    if (!chamberInTrack[4] && !chamberInTrack[5]) trackOK = kFALSE;
    Int_t nHitsInLastStations = 0;
    for (Int_t iCh = 6; iCh < AliMUONConstants::NTrackingCh(); iCh++)
      if (chamberInTrack[iCh]) nHitsInLastStations++; 
    if(nHitsInLastStations < 3) trackOK = kFALSE;
    
    // Add reconstructible tracks to fRecoTrackRefStore
    if (trackOK) fRecoTrackRefStore->Add(*track);
    
    delete [] chamberInTrack;
  }

}

