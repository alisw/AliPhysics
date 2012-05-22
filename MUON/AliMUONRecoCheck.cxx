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
#include "AliMUONTriggerTrack.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONTriggerStoreV1.h"
#include "AliMUONLocalTrigger.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliMUONDigitStoreV2S.h"
#include "AliMUONVDigit.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpPad.h"

#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMUONCDB.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONTriggerElectronics.h"


#include "TGeoManager.h"

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <Riostream.h>

#include "AliMUONRecoCheck.h"

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
fRecoTriggerRefStore(0x0),
fRecoTrackStore(0x0),
fRecoTriggerTrackStore(0x0),
fGeometryTransformer(0x0),
fTriggerCircuit(0x0),
fCalibrationData(0x0),
fTriggerElectronics(0x0),
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
fRecoTriggerRefStore(0x0),
fRecoTrackStore(0x0),
fRecoTriggerTrackStore(0x0),
fGeometryTransformer(0x0),
fTriggerCircuit(0x0),
fCalibrationData(0x0),
fTriggerElectronics(0x0),
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
  delete fGeometryTransformer;
  delete fTriggerCircuit;
  delete fTriggerElectronics;
  delete fCalibrationData;
}

//_____________________________________________________________________________
void AliMUONRecoCheck::ResetStores()
{
  /// Deletes all the store objects that have been created and resets the pointers to 0x0
  delete fTrackRefStore;      fTrackRefStore = 0x0;
  delete fRecoTrackRefStore;  fRecoTrackRefStore = 0x0;
  delete fRecoTriggerRefStore;  fRecoTriggerRefStore = 0x0;
  delete fRecoTrackStore;     fRecoTrackStore = 0x0;
  delete fRecoTriggerTrackStore; fRecoTriggerTrackStore = 0x0;
}

//_____________________________________________________________________________
Bool_t AliMUONRecoCheck::InitCircuit()
{
  
  if ( fTriggerCircuit ) return kTRUE;
  
  if ( ! InitGeometryTransformer() ) return kFALSE;
  
  // reset tracker for local trigger to trigger track conversion
  if ( ! AliMUONESDInterface::GetTracker() ) AliMUONESDInterface::ResetTracker();
  
  fTriggerCircuit = new AliMUONTriggerCircuit(fGeometryTransformer);
  
  return kTRUE;
}


//_____________________________________________________________________________
Int_t AliMUONRecoCheck::GetRunNumber()
{
  /// Return the run number of the current ESD event
  
  if (fESDEventOwner && fRecoTrackStore == 0x0 && fRecoTriggerTrackStore == 0x0) {
    if (!fESDTree || fESDTree->GetEvent(fCurrentEvent) <= 0) {
      AliError(Form("fails to read ESD object for event %d: cannot get the run number",fCurrentEvent));
      return -1;
    }
  }
  
  return fESDEvent->GetRunNumber();
}

//_____________________________________________________________________________
Int_t AliMUONRecoCheck::NumberOfEvents() const
{
  /// Return the number of events
  if (fESDEventOwner && fESDTree) return fESDTree->GetEntries();
  return 0;
}

//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONRecoCheck::ReconstructedTracks(Int_t event, Bool_t refit)
{
  /// Return a track store containing the reconstructed tracks (converted into 
  /// MUONTrack objects) for a given event.
  /// Track parameters at each clusters are computed or not depending on the flag "refit".
  /// If not, only the track parameters at first cluster are valid.
  
  if (!fESDEventOwner) {
    if (fRecoTrackStore == 0x0) MakeReconstructedTracks(refit);
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
    MakeReconstructedTracks(refit);
    return fRecoTrackStore;
  }
}


//_____________________________________________________________________________
AliMUONVTriggerTrackStore* AliMUONRecoCheck::TriggeredTracks(Int_t event)
{
  /// Return a track store containing the reconstructed trigger tracks (converted into 
  /// MUONTriggerTrack objects) for a given event.
	
  if (!fESDEventOwner) {
    if (fRecoTriggerTrackStore == 0x0) MakeTriggeredTracks();
    return fRecoTriggerTrackStore;
  }
	
  if (event != fCurrentEvent) {
    ResetStores();
    fCurrentEvent = event;
  }
	
  if (fRecoTriggerTrackStore != 0x0) return fRecoTriggerTrackStore;
  else {
    if (!fESDTree) return 0x0;
    if (fESDTree->GetEvent(event) <= 0) {
      AliError(Form("fails to read ESD object for event %d", event));
      return 0x0;
    }
    MakeTriggeredTracks();
    return fRecoTriggerTrackStore;
  }
}


//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONRecoCheck::TrackRefs(Int_t event)
{
  /// Return a track store containing the track references (converted into 
  /// MUONTrack objects) for a given event
  
  if (!fESDEventOwner) {
    if (fTrackRefStore == 0x0) MakeTrackRefs();
    return fTrackRefStore;
  }

  if (event != fCurrentEvent) {
    ResetStores();
    fCurrentEvent = event;
  }
  
  if (fTrackRefStore != 0x0) return fTrackRefStore;
  else {
    if (!fMCEventHandler->LoadEvent(event)) {
      AliError(Form("fails to read MC objects for event %d", event));
      return 0x0;
    }
    MakeTrackRefs();
    return fTrackRefStore;
  }
}

//_____________________________________________________________________________
AliMUONVTriggerTrackStore* AliMUONRecoCheck::TriggerableTracks(Int_t event)
{
  /// Return a trigger track store containing the triggerable track references (converted into 
  /// AliMUONTriggerTrack objects) for a given event
	
  if (!fESDEventOwner) {
    if (fRecoTriggerRefStore == 0x0) MakeTriggerableTracks();
    return fRecoTriggerRefStore;
  }
	
  if (event != fCurrentEvent) {
    ResetStores();
    fCurrentEvent = event;
  }
	
  if (fRecoTriggerRefStore != 0x0) return fRecoTriggerRefStore;
  else {
    if (!fMCEventHandler->LoadEvent(event)) {
      AliError(Form("fails to read MC objects for event %d", event));
      return 0x0;
    }
    MakeTriggerableTracks();
    return fRecoTriggerRefStore;
  }
}


//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONRecoCheck::ReconstructibleTracks(Int_t event, UInt_t requestedStationMask,
							    Bool_t request2ChInSameSt45, Bool_t hitInFrontOfPad)
{
  /// Return a track store containing the reconstructible tracks for a given event,
  /// according to the mask of requested stations and the minimum number of chambers hit in stations 4 & 5.

  if (!fESDEventOwner) {
    if (fRecoTrackRefStore == 0x0) {
      if (TrackRefs(event) == 0x0) return 0x0;
      MakeReconstructibleTracks(requestedStationMask, request2ChInSameSt45, hitInFrontOfPad);
    }
    return fRecoTrackRefStore;
  }

  if (event != fCurrentEvent) {
    ResetStores();
    fCurrentEvent = event;
  }
  
  if (fRecoTrackRefStore != 0x0) return fRecoTrackRefStore;
  else {
    if (TrackRefs(event) == 0x0) return 0x0;
    MakeReconstructibleTracks(requestedStationMask, request2ChInSameSt45, hitInFrontOfPad);
    return fRecoTrackRefStore;
  }
}


//_____________________________________________________________________________
void AliMUONRecoCheck::MakeReconstructedTracks(Bool_t refit)
{
  /// Make reconstructed tracks
  if (!(fRecoTrackStore = AliMUONESDInterface::NewTrackStore())) return;
  
  // loop over all reconstructed tracks and add them to the store (skip ghosts)
  Int_t nTracks = (Int_t) fESDEvent->GetNumberOfMuonTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliESDMuonTrack* esdTrack = fESDEvent->GetMuonTrack(iTrack);
    if (esdTrack->ContainTrackerData()) AliMUONESDInterface::Add(*esdTrack, *fRecoTrackStore, refit);
  }
  
}


//_____________________________________________________________________________
void AliMUONRecoCheck::MakeTriggeredTracks()
{
  /// Make reconstructed trigger tracks
  if (!(fRecoTriggerTrackStore = AliMUONESDInterface::NewTriggerTrackStore())) return;
	
  AliMUONVTriggerStore* tmpTriggerStore = AliMUONESDInterface::NewTriggerStore();
  if ( ! tmpTriggerStore ) return;
  
  // loop over all reconstructed tracks and add them to the store (include ghosts)
  Int_t nTracks = (Int_t) fESDEvent->GetNumberOfMuonTracks();
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliESDMuonTrack* esdTrack = fESDEvent->GetMuonTrack(iTrack);
    if (esdTrack->ContainTriggerData()) AliMUONESDInterface::Add(*esdTrack, *tmpTriggerStore);
  }
	
  if ( InitCircuit() ) {
  
    AliMUONVTrackReconstructor* tracker = AliMUONESDInterface::GetTracker();
    tracker->EventReconstructTrigger(*fTriggerCircuit, *tmpTriggerStore, *fRecoTriggerTrackStore);
  }
  
  delete tmpTriggerStore;
}

//_____________________________________________________________________________
void AliMUONRecoCheck::TriggerToTrack(const AliMUONLocalTrigger& locTrg, AliMUONTriggerTrack& triggerTrack)
{
  /// Make trigger track from local trigger info
  if ( ! InitCircuit() ) return;
  AliMUONVTrackReconstructor* tracker = AliMUONESDInterface::GetTracker();
  tracker->TriggerToTrack(*fTriggerCircuit, locTrg, triggerTrack);
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
  AliMUONVCluster* hit = AliMUONESDInterface::NewCluster();
  
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
    track.SetUniqueID(iTrackRef);
    tmpTrackRefStore->Add(track);
  }
  
  CleanMuonTrackRef(tmpTrackRefStore);
  
  delete hit;
  delete tmpTrackRefStore;
}

//_____________________________________________________________________________
void AliMUONRecoCheck::MakeTriggerableTracks()
{
  /// Make triggerable tracks
 if (!(fRecoTriggerRefStore = AliMUONESDInterface::NewTriggerTrackStore()))
   return;

  Double_t x, y, z, slopeX, slopeY, pZ, xLoc, yLoc, zLoc;
  TParticle* particle;
  TClonesArray* trackRefs;
  Int_t nTrackRef = fMCEventHandler->MCEvent()->GetNumberOfTracks();
	
  // loop over simulated tracks
  for (Int_t iTrackRef  = 0; iTrackRef < nTrackRef; ++iTrackRef) {
    Int_t nHits = fMCEventHandler->GetParticleAndTR(iTrackRef, particle, trackRefs);
		
    // skip empty trackRefs
    if (nHits < 1) continue;
				
    AliMUONTriggerTrack track;
    AliMUONDigitStoreV2S digitStore;
    Int_t hitsOnTrigger = 0;
    Int_t currCh = -1;
		
    // loop over simulated track hits
    for (Int_t iHit = 0; iHit < nHits; ++iHit) {        
      AliTrackReference* trackReference = static_cast<AliTrackReference*>(trackRefs->UncheckedAt(iHit));
			
      // skip trackRefs not in MUON
      if (trackReference->DetectorId() != AliTrackReference::kMUON) continue;

      // check chamberId of current trackReference
      Int_t detElemId = trackReference->UserId();
      Int_t chamberId = detElemId / 100 - 1;
      if (chamberId < AliMUONConstants::NTrackingCh() || chamberId >= AliMUONConstants::NCh() ) continue;

      // Get track parameters of current hit
      x = trackReference->X();
      y = trackReference->Y();
      z = trackReference->Z();
      
      if ( InitTriggerResponse() ) {
        fGeometryTransformer->Global2Local(detElemId, x, y, z, xLoc, yLoc, zLoc);
      
        Int_t nboard = 0;
        for ( Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; ++cath )
        {
          const AliMpVSegmentation* seg 
          = AliMpSegmentation::Instance()
          ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));
        
          AliMpPad pad = seg->PadByPosition(xLoc,yLoc,kFALSE);
          Int_t ix = pad.GetIx();
          Int_t iy = pad.GetIy();
        
          if ( !pad.IsValid() ) continue;
        
          if ( cath == AliMp::kCath0 ) nboard = pad.GetLocalBoardId(0);
        
          AliMUONVDigit* digit = digitStore.Add(detElemId, nboard, pad.GetLocalBoardChannel(0),
						cath, AliMUONVDigitStore::kDeny);
          if (digit) {
	    digit->SetPadXY(ix,iy);
	    digit->SetCharge(1.);
	  }
        }
      }
    
      if ( hitsOnTrigger == 0 ) {
        pZ = trackReference->Pz();
        slopeX = ( pZ == 0. ) ? 99999. : trackReference->Px() / pZ;
        slopeY = ( pZ == 0. ) ? 99999. : trackReference->Py() / pZ;

        track.SetX11(x);
        track.SetY11(y);
        track.SetZ11(z);
        track.SetSlopeX(slopeX);
        track.SetSlopeY(slopeY);
      }

      if ( currCh != chamberId ) {
        hitsOnTrigger++;
        currCh = chamberId;
      }

    } // loop on hits
    
    if ( hitsOnTrigger < 3 ) continue;

    // Check if the track passes the trigger algorithm
    if ( InitTriggerResponse() ) {
      AliMUONTriggerStoreV1 triggerStore;
      fTriggerElectronics->Digits2Trigger(digitStore,triggerStore);
    
      TIter next(triggerStore.CreateIterator());
      AliMUONLocalTrigger* locTrg(0x0);
    
      Int_t ptCutLevel = 0;
    
      while ( ( locTrg = static_cast<AliMUONLocalTrigger*>(next()) ) )
      {
        if ( locTrg->IsTrigX() && locTrg->IsTrigY() ) 
        {
          ptCutLevel = TMath::Max(ptCutLevel, 1);
          if ( locTrg->LoHpt() ) ptCutLevel = TMath::Max(ptCutLevel, 3);
          else if ( locTrg->LoLpt() ) ptCutLevel = TMath::Max(ptCutLevel, 2);
        } // board is fired 
      } // end of loop on Local Trigger
      track.SetPtCutLevel(ptCutLevel);
    }
    track.SetUniqueID(iTrackRef);
    fRecoTriggerRefStore->Add(track);
  }
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
  AliMUONVCluster* hit = AliMUONESDInterface::NewCluster();
  
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
    
    newTrack.SetUniqueID(track->GetUniqueID());
    newTrack.SetTrackParamAtVertex(track->GetTrackParamAtVertex());
    fTrackRefStore->Add(newTrack);
    
  }
  
  delete hit;
}

//_____________________________________________________________________________
void AliMUONRecoCheck::MakeReconstructibleTracks(UInt_t requestedStationMask, Bool_t request2ChInSameSt45,
						 Bool_t hitInFrontOfPad)
{
  /// Isolate the reconstructible tracks
  if (!(fRecoTrackRefStore = AliMUONESDInterface::NewTrackStore())) return;
  
  // need geometry transformer and segmentation to check if trackrefs are located in front of pad(s)
  if (hitInFrontOfPad) {
    if (!InitGeometryTransformer()) return;
    if (!AliMpSegmentation::Instance(kFALSE) && !AliMUONCDB::LoadMapping(kTRUE)) return;
  }
  
  // create iterator on trackRef
  TIter next(fTrackRefStore->CreateIterator());
  
  // loop over trackRef and add reconstructible tracks to fRecoTrackRefStore
  AliMUONTrack* track;
  while ( ( track = static_cast<AliMUONTrack*>(next()) ) ) {
    
    // check if the track is valid as it is
    if (track->IsValid(requestedStationMask, request2ChInSameSt45)) {
      
      if (hitInFrontOfPad) {
	
	AliMUONTrack trackTmp(*track);
	
	// loop over clusters and remove those which are not in front of pad(s)
	AliMUONTrackParam *param = static_cast<AliMUONTrackParam*>(trackTmp.GetTrackParamAtCluster()->First());
	while (param) {
	  AliMUONTrackParam *nextParam = static_cast<AliMUONTrackParam*>(trackTmp.GetTrackParamAtCluster()->After(param));
	  if (!IsHitInFrontOfPad(param)) trackTmp.RemoveTrackParamAtCluster(param);
	  param = nextParam;
	}
	
	// add the track if it is still valid
	if (trackTmp.IsValid(requestedStationMask, request2ChInSameSt45)) fRecoTrackRefStore->Add(trackTmp);
	
      } else {
	
	// add the track
	fRecoTrackRefStore->Add(*track);
	
      }
      
    }
      
  }

}

//_____________________________________________________________________________
AliMUONTrack* AliMUONRecoCheck::FindCompatibleTrack(AliMUONTrack &track, AliMUONVTrackStore &trackStore,
						    Int_t &nMatchClusters, Bool_t useLabel, Double_t sigmaCut)
{
  /// Return the track from the store matched with the given track (or 0x0) and the number of matched clusters.
  /// Matching is done by using the MC label of by comparing cluster/TrackRef positions according to the flag "useLabel".
  /// WARNING: Who match who matters since the matching algorithm uses the *fraction* of matched clusters of the given track
  
  AliMUONTrack *matchedTrack = 0x0;
  nMatchClusters = 0;
  
  if (useLabel) { // by using the MC label
    
    // get the corresponding simulated track if any
    Int_t label = track.GetMCLabel();
    matchedTrack = (AliMUONTrack*) trackStore.FindObject(label);
    
    // get the fraction of matched clusters
    if (matchedTrack) {
      Int_t nClusters = track.GetNClusters();
      for (Int_t iCl = 0; iCl < nClusters; iCl++)
	if (((AliMUONTrackParam*) track.GetTrackParamAtCluster()->UncheckedAt(iCl))->GetClusterPtr()->GetMCLabel() == label)
	  nMatchClusters++;
    }
    
  } else { // by comparing cluster/TrackRef positions
    
    // look for the corresponding simulated track if any
    TIter next(trackStore.CreateIterator());
    AliMUONTrack* track2;
    while ( ( track2 = static_cast<AliMUONTrack*>(next()) ) ) {
      
      // check compatibility
      Int_t n = 0;
      if (track.Match(*track2, sigmaCut, n)) {
	matchedTrack = track2;
	nMatchClusters = n;
	break;
      }
      
    }
    
  }
  
  return matchedTrack;
  
}


//_____________________________________________________________________________
AliMUONTriggerTrack* AliMUONRecoCheck::FindCompatibleTrack(AliMUONTriggerTrack &track, const AliMUONVTriggerTrackStore &triggerTrackStore,
                                                           Double_t sigmaCut)
{
  /// Return the trigger track from the store matched with the given track (or 0x0).
  /// Matching is done by comparing cluster/TrackRef positions.
  
  AliMUONTriggerTrack *matchedTrack = 0x0;
  
  // look for the corresponding simulated track if any
  TIter next(triggerTrackStore.CreateIterator());
  AliMUONTriggerTrack* track2;
  while ( ( track2 = static_cast<AliMUONTriggerTrack*>(next()) ) ) {
    // check compatibility
    if (track.Match(*track2, sigmaCut)) {
      matchedTrack = track2;
      break;
    }
  }
  
  return matchedTrack;
  
}


//____________________________________________________________________________ 
Bool_t AliMUONRecoCheck::InitTriggerResponse()
{
  /// Initialize trigger electronics
  /// for building of triggerable tracks from MC
  
  if ( fTriggerElectronics ) return kTRUE;
  
  if ( ! InitGeometryTransformer() ) return kFALSE;
  
  if ( ! AliMpDDLStore::Instance(false) && ! AliMpCDB::LoadDDLStore()) return kFALSE;
  
  if ( ! InitCalibrationData() ) return kFALSE;
  
  fTriggerElectronics = new AliMUONTriggerElectronics(fCalibrationData);
  
  return kTRUE;
}


//____________________________________________________________________________ 
Bool_t AliMUONRecoCheck::InitCalibrationData()
{
  /// Initialize calibration data
  if ( ! fCalibrationData ) {
    if ( !AliMUONCDB::CheckOCDB() ) return kFALSE;
    fCalibrationData = new AliMUONCalibrationData(AliCDBManager::Instance()->GetRun());
  }
  return kTRUE;
}


//____________________________________________________________________________ 
Bool_t AliMUONRecoCheck::InitGeometryTransformer()
{
  /// Return geometry transformer
  /// (create it if necessary)
  if ( ! fGeometryTransformer ) {
    
    if ( !AliMUONCDB::CheckOCDB() ) return kFALSE;
    
    if (!AliGeomManager::GetGeometry()) {
      AliGeomManager::LoadGeometry();
      if (!AliGeomManager::GetGeometry()) return kFALSE;  
      if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return kFALSE;
    }
    
    fGeometryTransformer = new AliMUONGeometryTransformer();
    fGeometryTransformer->LoadGeometryData();
  }
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliMUONRecoCheck::IsHitInFrontOfPad(AliMUONTrackParam *param) const
{
  /// Return kTRUE if the hit is located in front of at least 1 pad
  
  Int_t deId = param->GetClusterPtr()->GetDetElemId();
  
  Double_t xL, yL, zL;
  fGeometryTransformer->Global2Local(deId, param->GetNonBendingCoor(), param->GetBendingCoor(), param->GetZ(), xL, yL, zL);
  
  const AliMpVSegmentation* seg1 = AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath0);
  const AliMpVSegmentation* seg2 = AliMpSegmentation::Instance()->GetMpSegmentation(deId, AliMp::kCath1);
  if (!seg1 || !seg2) return kFALSE;
  
  AliMpPad pad1 = seg1->PadByPosition(xL, yL, kFALSE);
  AliMpPad pad2 = seg2->PadByPosition(xL, yL, kFALSE);
  
  return (pad1.IsValid() || pad2.IsValid());
}

