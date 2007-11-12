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

#include "AliMUONRecoCheck.h"
#include "AliMUONRawClusterV2.h"
#include "AliMUONTrack.h"
#include "AliMUONConstants.h"
#include "AliMUONDataInterface.h"
#include "AliMUONTrackStoreV1.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliLog.h" 

#include <TParticle.h>
#include <TParticlePDG.h>

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONRecoCheck)
/// \endcond

//_____________________________________________________________________________
AliMUONRecoCheck::AliMUONRecoCheck(Char_t *chLoader, Char_t *pathSim)
: TObject(),
fMCEventHandler(new AliMCEventHandler()),
fDataInterface(new AliMUONDataInterface(chLoader)),
fCurrentEvent(0),
fTrackRefStore(0x0),
fRecoTrackRefStore(0x0),
fRecoTrackStore(0x0)
{
  /// Normal ctor
  fMCEventHandler->SetInputPath(pathSim);
  fMCEventHandler->InitIO("");
}

//_____________________________________________________________________________
AliMUONRecoCheck::~AliMUONRecoCheck()
{
  /// Destructor
  delete fMCEventHandler;
  delete fDataInterface;
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
  if (fDataInterface->IsValid()) return fDataInterface->NumberOfEvents();
  return 0;
}

//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONRecoCheck::ReconstructedTracks(Int_t event)
{
  /// Return the reconstructed track store for a given event
  if (!fDataInterface->IsValid()) return new AliMUONTrackStoreV1();
  
  if (event != fCurrentEvent) {
    ResetStores();
    fCurrentEvent = event;
  }
  
  if (fRecoTrackStore != 0x0) return fRecoTrackStore;
  
  fRecoTrackStore = new AliMUONTrackStoreV1();
  
  return fRecoTrackStore;
}

//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONRecoCheck::TrackRefs(Int_t event)
{
  /// Return a track store containing the track references (converted into 
  /// MUONTrack objects) for a given event
  if (event != fCurrentEvent) {
    ResetStores();
    fCurrentEvent = event;
  }
  
  if (fTrackRefStore != 0x0) return fTrackRefStore;
  else {
    if (!fMCEventHandler->GetEvent(event)) return 0x0;
    MakeTrackRefs();
    return fTrackRefStore;
  }
}

//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONRecoCheck::ReconstructibleTracks(Int_t event)
{
  /// Return a track store containing the reconstructible tracks for a given event
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
void AliMUONRecoCheck::MakeTrackRefs()
{
  /// Make reconstructible tracks
  AliMUONVTrackStore *tmpTrackRefStore = new AliMUONTrackStoreV1();
  
  Double_t x, y, z, pX, pY, pZ, bendingSlope, nonBendingSlope, inverseBendingMomentum;
  TParticle* particle;
  TClonesArray* trackRefs;
  Int_t nTrackRef = fMCEventHandler->MCEvent()->GetNumberOfTracks();
  
  // loop over simulated tracks
  for (Int_t iTrackRef  = 0; iTrackRef < nTrackRef; ++iTrackRef) {
    Int_t nHits = fMCEventHandler->GetParticleAndTR(iTrackRef, particle, trackRefs);
    
    // skip empty trackRefs
    if (nHits < 1) continue;
    
    // skip trackRefs not in MUON
    AliTrackReference* trackReference0 = static_cast<AliTrackReference*>(trackRefs->UncheckedAt(0));
    if (trackReference0->DetectorId() != AliTrackReference::kMUON) continue;
    
    AliMUONTrack track;
    
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
    TParticlePDG* ppdg = particle->GetPDG(1);
    Int_t charge = (Int_t)(ppdg->Charge()/3.0);
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
    
    // loop over simulated track hits
    for (Int_t iHit = 0; iHit < nHits; ++iHit) {        
      AliTrackReference* trackReference = static_cast<AliTrackReference*>(trackRefs->UncheckedAt(iHit));
      
      // Get track parameters of current hit
      x = trackReference->X();
      y = trackReference->Y();
      z = trackReference->Z();
      pX = trackReference->Px();
      pY = trackReference->Py();
      pZ = trackReference->Pz();
      
      // check chamberId of current trackReference
      Int_t chamberId = AliMUONConstants::ChamberNumber(z);
      if (chamberId < 0 || chamberId >= AliMUONConstants::NTrackingCh()) continue;
      
      // set hit parameters
      AliMUONRawClusterV2 hit(chamberId, 0, 0);
      hit.SetXYZ(x,y,z);
      hit.SetErrXY(0.,0.);
      
      // compute track parameters at hit
      Double_t bendingSlope = 0;
      Double_t nonBendingSlope = 0;
      Double_t inverseBendingMomentum = 0;
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
      track.AddTrackParamAtCluster(trackParam,hit,kTRUE);
    }
    
    track.GetTrackParamAtCluster()->Sort();
    track.SetTrackID(iTrackRef);
    tmpTrackRefStore->Add(track);
  }
  
  CleanMuonTrackRef(tmpTrackRefStore);
  
  delete tmpTrackRefStore;
}

//_____________________________________________________________________________
void AliMUONRecoCheck::CleanMuonTrackRef(const AliMUONVTrackStore *tmpTrackRefStore)
{
  /// Re-calculate hits parameters because two AliTrackReferences are recorded for
  /// each chamber (one when particle is entering + one when particle is leaving 
  /// the sensitive volume) 
  fTrackRefStore = new AliMUONTrackStoreV1();
  
  Double_t maxGasGap = 1.; // cm 
  Double_t x, y, z, pX, pY, pZ, x1, y1, z1, pX1, pY1, pZ1, z2;
  Double_t bendingSlope,nonBendingSlope,inverseBendingMomentum;
  
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
      AliMUONRawClusterV2 hit(trackParam1->GetClusterPtr()->GetChamberId(), 0, 0);
      hit.SetXYZ(x,y,z);
      hit.SetErrXY(0.,0.);
      
      // set new track parameters at new hit
      AliMUONTrackParam trackParam;
      trackParam.SetNonBendingCoor(x);
      trackParam.SetBendingCoor(y);
      trackParam.SetZ(z);
      trackParam.SetBendingSlope(bendingSlope);
      trackParam.SetNonBendingSlope(nonBendingSlope);
      trackParam.SetInverseBendingMomentum(inverseBendingMomentum);
      
      // add track parameters at current hit to the track
      newTrack.AddTrackParamAtCluster(trackParam,hit,kTRUE);
      
      iHit1++;
    }
    
    newTrack.GetTrackParamAtCluster()->Sort();
    newTrack.SetTrackID(track->GetTrackID());
    newTrack.SetTrackParamAtVertex(track->GetTrackParamAtVertex());
    fTrackRefStore->Add(newTrack);
    
  }
  
}

//_____________________________________________________________________________
void AliMUONRecoCheck::MakeReconstructibleTracks()
{
  /// Calculate the number of reconstructible tracks
  fRecoTrackRefStore = new AliMUONTrackStoreV1();
  
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

