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
// Class AliMUONTrack
//-------------------
// Reconstructed track in ALICE dimuon spectrometer
//-----------------------------------------------------------------------------

#include "AliMUONTrack.h"

#include "AliMUONTrackParam.h" 
#include "AliMUONHitForRec.h" 
#include "AliMUONObjectPair.h" 
#include "AliMUONConstants.h"
#include "AliMUONTrackExtrap.h" 

#include "AliLog.h"

#include <TMath.h>
#include <TMatrixD.h>

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrack) // Class implementation in ROOT context
/// \endcond

//__________________________________________________________________________
AliMUONTrack::AliMUONTrack()
  : TObject(),
    fTrackParamAtVertex(),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(0),
    fFitWithVertex(kFALSE),
    fVertex(0x0),
    fFitWithMCS(kFALSE),
    fHitWeightsNonBending(0x0),
    fHitWeightsBending(0x0),
    fGlobalChi2(-1.),
    fImproved(kFALSE),
    fMatchTrigger(-1),
    floTrgNum(-1),
    fChi2MatchTrigger(0.),
    fTrackID(0),
    fHitsPatternInTrigCh(0),
    fLocalTrigger(0)
{
  /// Default constructor
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack(AliMUONObjectPair *segment)
  : TObject(),
    fTrackParamAtVertex(),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(0),
    fFitWithVertex(kFALSE),
    fVertex(0x0),
    fFitWithMCS(kFALSE),
    fHitWeightsNonBending(0x0),
    fHitWeightsBending(0x0),
    fGlobalChi2(0.),
    fImproved(kFALSE),
    fMatchTrigger(-1),
    floTrgNum(-1),    
    fChi2MatchTrigger(0.),
    fTrackID(0),
    fHitsPatternInTrigCh(0),
    fLocalTrigger(0)
{
  /// Constructor from thw hitForRec's

  fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",10);
  fTrackParamAtHit->SetOwner(kTRUE);
  fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",10);
  fHitForRecAtHit->SetOwner(kTRUE);
  
  if (!segment) return; //AZ
  
  // Pointers to hits from the segment
  AliMUONHitForRec* hit1 = (AliMUONHitForRec*) segment->First();
  AliMUONHitForRec* hit2 = (AliMUONHitForRec*) segment->Second();
  
  // check sorting in -Z (spectro z<0)
  if (hit1->GetZ() < hit2->GetZ()) {
    hit1 = hit2;
    hit2 = (AliMUONHitForRec*) segment->First();
  }
  
  // order the hits into the track according to the station the segment belong to
  //(the hit first attached is the one from which we will start the tracking procedure)
  if (hit1->GetChamberNumber() == 8) {
    AddTrackParamAtHit(0,hit1);
    AddTrackParamAtHit(0,hit2);
  } else {
    AddTrackParamAtHit(0,hit2);
    AddTrackParamAtHit(0,hit1);
  }
  
  AliMUONTrackParam* trackParamAtFirstHit = (AliMUONTrackParam*) fTrackParamAtHit->First();
  AliMUONHitForRec* firstHit = trackParamAtFirstHit->GetHitForRecPtr();
  AliMUONTrackParam* trackParamAtLastHit = (AliMUONTrackParam*) fTrackParamAtHit->Last();
  AliMUONHitForRec* lastHit = trackParamAtLastHit->GetHitForRecPtr();
  
  
  // Compute track parameters
  Double_t dZ = firstHit->GetZ() - lastHit->GetZ();
  // Non bending plane
  Double_t nonBendingCoor1 = firstHit->GetNonBendingCoor();
  Double_t nonBendingCoor2 = lastHit->GetNonBendingCoor();
  Double_t nonBendingSlope = (nonBendingCoor1 - nonBendingCoor2) / dZ;
  // Bending plane
  Double_t bendingCoor1 = firstHit->GetBendingCoor();
  Double_t bendingCoor2 = lastHit->GetBendingCoor();
  Double_t bendingSlope = (bendingCoor1 - bendingCoor2) / dZ;
  // Inverse bending momentum
  Double_t bendingImpact = bendingCoor1 - firstHit->GetZ() * bendingSlope;
  Double_t inverseBendingMomentum = 1. / AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(bendingImpact);
  
  
  // Set track parameters at first hit
  trackParamAtFirstHit->SetNonBendingCoor(nonBendingCoor1);
  trackParamAtFirstHit->SetNonBendingSlope(nonBendingSlope);
  trackParamAtFirstHit->SetBendingCoor(bendingCoor1);
  trackParamAtFirstHit->SetBendingSlope(bendingSlope);
  trackParamAtFirstHit->SetInverseBendingMomentum(inverseBendingMomentum);
  
  
  // Set track parameters at last hit
  trackParamAtLastHit->SetNonBendingCoor(nonBendingCoor2);
  trackParamAtLastHit->SetNonBendingSlope(nonBendingSlope);
  trackParamAtLastHit->SetBendingCoor(bendingCoor2);
  trackParamAtLastHit->SetBendingSlope(bendingSlope);
  trackParamAtLastHit->SetInverseBendingMomentum(inverseBendingMomentum);
  
  
  // Compute and set track parameters covariances at first hit
  TMatrixD paramCov1(5,5);
  paramCov1.Zero();
  // Non bending plane
  paramCov1(0,0) = firstHit->GetNonBendingReso2();
  paramCov1(0,1) = firstHit->GetNonBendingReso2() / dZ;
  paramCov1(1,0) = paramCov1(0,1);
  paramCov1(1,1) = ( firstHit->GetNonBendingReso2() + lastHit->GetNonBendingReso2() ) / dZ / dZ;
  // Bending plane
  paramCov1(2,2) = firstHit->GetBendingReso2();
  paramCov1(2,3) = firstHit->GetBendingReso2() / dZ;
  paramCov1(3,2) = paramCov1(2,3);
  paramCov1(3,3) = ( firstHit->GetBendingReso2() + lastHit->GetBendingReso2() ) / dZ / dZ;
  // Inverse bending momentum (50% error)
  paramCov1(4,4) = 0.5*inverseBendingMomentum * 0.5*inverseBendingMomentum;
  // Set covariances
  trackParamAtFirstHit->SetCovariances(paramCov1);
  
  
  // Compute and set track parameters covariances at last hit (as if the first hit did not exist)
  TMatrixD paramCov2(5,5);
  paramCov2.Zero();
  // Non bending plane
  paramCov2(0,0) = paramCov1(0,0);
  paramCov2(1,1) = 100.*paramCov1(1,1);
  // Bending plane
  paramCov2(2,2) = paramCov1(2,2);
  paramCov2(3,3) = 100.*paramCov1(3,3);
  // Inverse bending momentum
  paramCov2(4,4) = paramCov1(4,4);
  // Set covariances
  trackParamAtLastHit->SetCovariances(paramCov2);
  
  
  // Flag first hit as being removable
  trackParamAtFirstHit->SetRemovable(kTRUE);
  
  // Flag last hit as being removable
  trackParamAtLastHit->SetRemovable(kTRUE);
  
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack (const AliMUONTrack& track)
  : TObject(track),
    fTrackParamAtVertex(track.fTrackParamAtVertex),
    fTrackParamAtHit(0x0),
    fHitForRecAtHit(0x0),
    fNTrackHits(track.fNTrackHits),
    fFitWithVertex(track.fFitWithVertex),
    fVertex(0x0),
    fFitWithMCS(track.fFitWithMCS),
    fHitWeightsNonBending(0x0),
    fHitWeightsBending(0x0),
    fGlobalChi2(track.fGlobalChi2),
    fImproved(track.fImproved),
    fMatchTrigger(track.fMatchTrigger),
    floTrgNum(track.floTrgNum),    
    fChi2MatchTrigger(track.fChi2MatchTrigger),
    fTrackID(track.fTrackID),
    fHitsPatternInTrigCh(track.fHitsPatternInTrigCh),
    fLocalTrigger(track.fLocalTrigger)
{
  ///copy constructor
  Int_t maxIndex = 0;
  
  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (track.fTrackParamAtHit) {
    maxIndex = (track.fTrackParamAtHit)->GetEntriesFast();
    fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",maxIndex);
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fTrackParamAtHit)[index]) AliMUONTrackParam(*(AliMUONTrackParam*)track.fTrackParamAtHit->At(index));
    }
  }
  
  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (track.fHitForRecAtHit) {
    maxIndex = (track.fHitForRecAtHit)->GetEntriesFast();
    fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",maxIndex);
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fHitForRecAtHit)[index]) AliMUONHitForRec(*(AliMUONHitForRec*)track.fHitForRecAtHit->At(index));
    }
  }
  
  // copy vertex used during the tracking procedure if any
  if (track.fVertex) fVertex = new AliMUONHitForRec(*(track.fVertex));
  
  // copy hit weights matrices if any
  if (track.fHitWeightsNonBending) fHitWeightsNonBending = new TMatrixD(*(track.fHitWeightsNonBending));
  if (track.fHitWeightsBending) fHitWeightsBending = new TMatrixD(*(track.fHitWeightsBending));
  
}

  //__________________________________________________________________________
AliMUONTrack & AliMUONTrack::operator=(const AliMUONTrack& track)
{
  /// Asignment operator
  // check assignement to self
  if (this == &track)
    return *this;

  // base class assignement
  TObject::operator=(track);

  fTrackParamAtVertex = track.fTrackParamAtVertex;

  Int_t maxIndex = 0;
  
  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (track.fTrackParamAtHit) {
    if (fTrackParamAtHit) fTrackParamAtHit->Clear();
    else fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",10);
    maxIndex = (track.fTrackParamAtHit)->GetEntriesFast();
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fTrackParamAtHit)[fTrackParamAtHit->GetEntriesFast()])
      	AliMUONTrackParam(*(AliMUONTrackParam*)(track.fTrackParamAtHit)->At(index));
    }
  } else if (fTrackParamAtHit) {
    delete fTrackParamAtHit;
    fTrackParamAtHit = 0x0;
  }

  // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (track.fHitForRecAtHit) {
    if (fHitForRecAtHit) fHitForRecAtHit->Clear();
    else fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",10);
    maxIndex = (track.fHitForRecAtHit)->GetEntriesFast();
    for (Int_t index = 0; index < maxIndex; index++) {
      new ((*fHitForRecAtHit)[fHitForRecAtHit->GetEntriesFast()])
      	AliMUONHitForRec(*(AliMUONHitForRec*)(track.fHitForRecAtHit)->At(index));
    }
  } else if (fHitForRecAtHit) {
    delete fHitForRecAtHit;
    fHitForRecAtHit = 0x0;
  }
  
  // copy vertex used during the tracking procedure if any.
  if (track.fVertex) {
    if (fVertex) *fVertex = *(track.fVertex);
    else fVertex = new AliMUONHitForRec(*(track.fVertex));
  } else if (fVertex) {
    delete fVertex;
    fVertex = 0x0;
  }
  
  // copy hit weights matrix if any
  if (track.fHitWeightsNonBending) {
    if (fHitWeightsNonBending) {
      fHitWeightsNonBending->ResizeTo(*(track.fHitWeightsNonBending));
      *fHitWeightsNonBending = *(track.fHitWeightsNonBending);
    } else fHitWeightsNonBending = new TMatrixD(*(track.fHitWeightsNonBending));
  } else if (fHitWeightsNonBending) {
    delete fHitWeightsNonBending;
    fHitWeightsNonBending = 0x0;
  }
  
  // copy hit weights matrix if any
  if (track.fHitWeightsBending) {
    if (fHitWeightsBending) {
      fHitWeightsBending->ResizeTo(*(track.fHitWeightsBending));
      *fHitWeightsBending = *(track.fHitWeightsBending);
    } else fHitWeightsBending = new TMatrixD(*(track.fHitWeightsBending));
  } else if (fHitWeightsBending) {
    delete fHitWeightsBending;
    fHitWeightsBending = 0x0;
  }
  
  fNTrackHits         =  track.fNTrackHits;
  fFitWithVertex      =  track.fFitWithVertex;
  fFitWithMCS         =  track.fFitWithMCS;
  fGlobalChi2         =  track.fGlobalChi2;
  fImproved           =  track.fImproved;
  fMatchTrigger       =  track.fMatchTrigger;
  floTrgNum           =  track.floTrgNum;
  fChi2MatchTrigger   =  track.fChi2MatchTrigger;
  fTrackID            =  track.fTrackID; 
  fHitsPatternInTrigCh = track.fHitsPatternInTrigCh;
  fLocalTrigger        = track.fLocalTrigger;

  return *this;
}

  //__________________________________________________________________________
AliMUONTrack::~AliMUONTrack()
{
  /// Destructor
  delete fTrackParamAtHit;
  delete fHitForRecAtHit;
  delete fVertex;
  delete fHitWeightsNonBending;
  delete fHitWeightsBending;
}

  //__________________________________________________________________________
void AliMUONTrack::Clear(Option_t* opt)
{
  /// Clear arrays
  if ( fTrackParamAtHit ) fTrackParamAtHit->Clear(opt);
  if ( fHitForRecAtHit ) fHitForRecAtHit->Clear(opt);
  delete fVertex; fVertex = 0x0;
  delete fHitWeightsNonBending; fHitWeightsNonBending = 0x0;
  delete fHitWeightsBending; fHitWeightsBending = 0x0;
}

  //__________________________________________________________________________
void AliMUONTrack::AddTrackParamAtHit(const AliMUONTrackParam *trackParam, AliMUONHitForRec *hitForRec)
{
  /// Add TrackParamAtHit if "trackParam" != NULL
  /// else create empty TrackParamAtHit and set the z position to the one of "hitForRec" if any
  /// Update link to HitForRec if "hitForRec" != NULL
  if (!fTrackParamAtHit) {
    fTrackParamAtHit = new TClonesArray("AliMUONTrackParam",10);  
    fNTrackHits = 0;
  }
  AliMUONTrackParam* trackParamAtHit;
  if (trackParam) {
    trackParamAtHit = new ((*fTrackParamAtHit)[fNTrackHits]) AliMUONTrackParam(*trackParam);
    if (hitForRec) {
      if (hitForRec->GetZ() != trackParam->GetZ())
        AliWarning("Added track parameters at a different z position than the one of the attached hit");
    }
  } else {
    trackParamAtHit = new ((*fTrackParamAtHit)[fNTrackHits]) AliMUONTrackParam();
    if (hitForRec) trackParamAtHit->SetZ(hitForRec->GetZ());
  }
  if (hitForRec) trackParamAtHit->SetHitForRecPtr(hitForRec);
  fNTrackHits++;
}

  //__________________________________________________________________________
void AliMUONTrack::RemoveTrackParamAtHit(AliMUONTrackParam *trackParam)
{
  /// Remove trackParam from the array of TrackParamAtHit
  if (!fTrackParamAtHit) {
    AliWarning("array fTrackParamAtHit does not exist");
    return;
  }
  
  if (!fTrackParamAtHit->Remove(trackParam)) {
    AliWarning("object to remove does not exist in array fTrackParamAtHit");
    return;
  }
  
  fTrackParamAtHit->Compress();
  fNTrackHits--;
}

  //__________________________________________________________________________
void AliMUONTrack::AddHitForRecAtHit(const AliMUONHitForRec *hitForRec) 
{
  /// Add hitForRec to the array of hitForRec at hit
  if (!fHitForRecAtHit)
    fHitForRecAtHit = new TClonesArray("AliMUONHitForRec",10); 
  
  if (!hitForRec)
    AliFatal("AliMUONTrack::AddHitForRecAtHit: hitForRec == NULL");
  
  new ((*fHitForRecAtHit)[fHitForRecAtHit->GetEntriesFast()]) AliMUONHitForRec(*hitForRec);
}

  //__________________________________________________________________________
void AliMUONTrack::UpdateTrackParamAtHit()
{
  /// Update track parameters at each attached hit
  
  if (fNTrackHits == 0) {
    AliWarning("no hit attached to the track");
    return;
  }
  
  Double_t z;
  AliMUONTrackParam* startingTrackParam = (AliMUONTrackParam*) fTrackParamAtHit->First();
  AliMUONTrackParam* trackParamAtHit = (AliMUONTrackParam*) fTrackParamAtHit->After(startingTrackParam);
  while (trackParamAtHit) {
    
    // save current z
    z = trackParamAtHit->GetZ();
    
    // reset track parameters and their covariances
    trackParamAtHit->SetParameters(startingTrackParam->GetParameters());
    trackParamAtHit->SetZ(startingTrackParam->GetZ());
    
    // extrapolation to the given z
    AliMUONTrackExtrap::ExtrapToZ(trackParamAtHit, z);
    
    // prepare next step
    startingTrackParam = trackParamAtHit;
    trackParamAtHit = (AliMUONTrackParam*) (fTrackParamAtHit->After(trackParamAtHit));
  }

}

  //__________________________________________________________________________
void AliMUONTrack::UpdateCovTrackParamAtHit()
{
  /// Update track parameters and their covariances at each attached hit
  
  if (fNTrackHits == 0) {
    AliWarning("no hit attached to the track");
    return;
  }
  
  Double_t z;
  AliMUONTrackParam* startingTrackParam = (AliMUONTrackParam*) fTrackParamAtHit->First();
  AliMUONTrackParam* trackParamAtHit = (AliMUONTrackParam*) fTrackParamAtHit->After(startingTrackParam);
  while (trackParamAtHit) {
    
    // save current z
    z = trackParamAtHit->GetZ();
    
    // reset track parameters and their covariances
    trackParamAtHit->SetParameters(startingTrackParam->GetParameters());
    trackParamAtHit->SetZ(startingTrackParam->GetZ());
    trackParamAtHit->SetCovariances(startingTrackParam->GetCovariances());
    
    // extrapolation to the given z
    AliMUONTrackExtrap::ExtrapToZCov(trackParamAtHit, z);
    
    // prepare next step
    startingTrackParam = trackParamAtHit;
    trackParamAtHit = (AliMUONTrackParam*) (fTrackParamAtHit->After(trackParamAtHit));
  }

}

  //__________________________________________________________________________
void AliMUONTrack::SetVertex(const AliMUONHitForRec* vertex)
{
  /// Set the vertex used during the tracking procedure
  if (!fVertex) fVertex = new AliMUONHitForRec(*vertex);
  else *fVertex = *vertex;
}


  //__________________________________________________________________________
Bool_t AliMUONTrack::ComputeLocalChi2(Bool_t accountForMCS)
{
  /// Compute the removable hit contribution to the chi2 of the track
  /// accounting for multiple scattering or not according to the flag
  /// - Also recompute the weight matrices of the attached hits if accountForMCS=kTRUE
  /// - Assume that track parameters at each hit are corrects
  /// - Return kFALSE if computation failed
  
  // Check hits (if the first one exist, assume that the other ones exit too!)
  AliMUONTrackParam* trackParamAtHit = (AliMUONTrackParam*) fTrackParamAtHit->First();
  if (!trackParamAtHit->GetHitForRecPtr()) {
    AliWarning("hit is missing");
    return kFALSE;
  }
  
  if (accountForMCS) { // Compute local chi2 taking into account multiple scattering effects
      
    // Compute MCS covariance matrix only once
    TMatrixD mcsCovariances(fNTrackHits,fNTrackHits);
    ComputeMCSCovariances(mcsCovariances);
    
    // Make sure hit weights are consistent with following calculations
    if (!ComputeHitWeights(&mcsCovariances)) {
      AliWarning("cannot take into account the multiple scattering effects");
      return ComputeLocalChi2(kFALSE);
    }
    
    // Compute chi2 of the track
    Double_t globalChi2 = ComputeGlobalChi2(kTRUE);
    if (globalChi2 < 0.) return kFALSE;
    
    // Loop over removable hits and compute their local chi2
    AliMUONTrackParam* trackParamAtHit1;
    AliMUONHitForRec *hitForRec, *discardedHit;
    Int_t hitNumber1, hitNumber2, currentHitNumber1, currentHitNumber2;
    TMatrixD hitWeightsNB(fNTrackHits-1,fNTrackHits-1);
    TMatrixD hitWeightsB(fNTrackHits-1,fNTrackHits-1);
    Double_t *dX = new Double_t[fNTrackHits-1];
    Double_t *dY = new Double_t[fNTrackHits-1];
    Double_t globalChi2b;
    while (trackParamAtHit) {
      
      discardedHit = trackParamAtHit->GetHitForRecPtr();
      
      // Recompute hit weights without the current hit
      if (!ComputeHitWeights(hitWeightsNB, hitWeightsB, &mcsCovariances, discardedHit)) {
  	AliWarning("cannot take into account the multiple scattering effects");
  	ComputeLocalChi2(kFALSE);
      }
      
      // Compute track chi2 without the current hit
      globalChi2b = 0.;
      currentHitNumber1 = 0;
      for (hitNumber1 = 0; hitNumber1 < fNTrackHits ; hitNumber1++) { 
    	trackParamAtHit1 = (AliMUONTrackParam*) fTrackParamAtHit->UncheckedAt(hitNumber1);
    	hitForRec = trackParamAtHit1->GetHitForRecPtr();
        
        if (hitForRec == discardedHit) continue;
        
        // Compute and save residuals
    	dX[currentHitNumber1] = hitForRec->GetNonBendingCoor() - trackParamAtHit1->GetNonBendingCoor();
    	dY[currentHitNumber1] = hitForRec->GetBendingCoor() - trackParamAtHit1->GetBendingCoor();
        
        currentHitNumber2 = 0;
    	for (hitNumber2 = 0; hitNumber2 < hitNumber1; hitNumber2++) {
    	  hitForRec = ((AliMUONTrackParam*) fTrackParamAtHit->UncheckedAt(hitNumber2))->GetHitForRecPtr();
          
          if (hitForRec == discardedHit) continue;
          
          // Add contribution from covariances
          globalChi2b += (hitWeightsNB(currentHitNumber1, currentHitNumber2) +
        		  hitWeightsNB(currentHitNumber2, currentHitNumber1)) * dX[currentHitNumber1] * dX[currentHitNumber2] +
        		 (hitWeightsB(currentHitNumber1, currentHitNumber2) +
        		  hitWeightsB(currentHitNumber2, currentHitNumber1)) * dY[currentHitNumber1] * dY[currentHitNumber2];
          
          currentHitNumber2++;
    	}
        
        // Add contribution from variances
    	globalChi2b += hitWeightsNB(currentHitNumber1, currentHitNumber1) * dX[currentHitNumber1] * dX[currentHitNumber1] +
        	       hitWeightsB(currentHitNumber1, currentHitNumber1) * dY[currentHitNumber1] * dY[currentHitNumber1];
    	
        currentHitNumber1++;
      }

      // Set local chi2
      trackParamAtHit->SetLocalChi2(globalChi2 - globalChi2b);
      
      trackParamAtHit = (AliMUONTrackParam*) fTrackParamAtHit->After(trackParamAtHit);
    }
    
    delete [] dX;
    delete [] dY;
    
  } else { // without multiple scattering effects
    
    AliMUONHitForRec *discardedHit;
    Double_t dX, dY;
    while (trackParamAtHit) {
      
      discardedHit = trackParamAtHit->GetHitForRecPtr();
      
      // Compute residuals
      dX = discardedHit->GetNonBendingCoor() - trackParamAtHit->GetNonBendingCoor();
      dY = discardedHit->GetBendingCoor() - trackParamAtHit->GetBendingCoor();
      
      // Set local chi2
      trackParamAtHit->SetLocalChi2(dX * dX / discardedHit->GetNonBendingReso2() + dY * dY / discardedHit->GetBendingReso2());
    
    trackParamAtHit = (AliMUONTrackParam*) fTrackParamAtHit->After(trackParamAtHit);
    }
  
  }
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Double_t AliMUONTrack::ComputeGlobalChi2(Bool_t accountForMCS)
{
  /// Compute the chi2 of the track accounting for multiple scattering or not according to the flag
  /// - Assume that track parameters at each hit are corrects
  /// - Assume the hits weights matrices are corrects
  /// - Return negative value if chi2 computation failed
  
  // Check hits (if the first one exist, assume that the other ones exit too!)
  AliMUONTrackParam* trackParamAtHit = (AliMUONTrackParam*) fTrackParamAtHit->First();
  if (!trackParamAtHit->GetHitForRecPtr()) {
    AliWarning("hit is missing");
    return -1.;
  }
  
  Double_t chi2 = 0.;
  
  if (accountForMCS) {
    
    // Check the weight matrices. If weight matrices are not available compute chi2 without MCS
    if (!fHitWeightsNonBending || !fHitWeightsBending) {
      AliWarning("hit weights including multiple scattering effects are not available\n\t\t --> compute chi2 WITHOUT multiple scattering");
      return ComputeGlobalChi2(kFALSE);
    }
    if (fHitWeightsNonBending->GetNrows() != fNTrackHits || fHitWeightsBending->GetNcols() != fNTrackHits) {
      AliWarning("hit weights including multiple scattering effects are not available\n\t\t --> compute chi2 WITHOUT multiple scattering");
      return ComputeGlobalChi2(kFALSE);
    }
    
    // Compute chi2
    AliMUONHitForRec *hitForRec;
    Double_t *dX = new Double_t[fNTrackHits];
    Double_t *dY = new Double_t[fNTrackHits];
    Int_t hitNumber1, hitNumber2;
    for (hitNumber1 = 0; hitNumber1 < fNTrackHits ; hitNumber1++) { 
      trackParamAtHit = (AliMUONTrackParam*) fTrackParamAtHit->UncheckedAt(hitNumber1);
      hitForRec = trackParamAtHit->GetHitForRecPtr();
      dX[hitNumber1] = hitForRec->GetNonBendingCoor() - trackParamAtHit->GetNonBendingCoor();
      dY[hitNumber1] = hitForRec->GetBendingCoor() - trackParamAtHit->GetBendingCoor();
      for (hitNumber2 = 0; hitNumber2 < hitNumber1; hitNumber2++) {
        chi2 += ((*fHitWeightsNonBending)(hitNumber1, hitNumber2) + (*fHitWeightsNonBending)(hitNumber2, hitNumber1)) * dX[hitNumber1] * dX[hitNumber2] +
		((*fHitWeightsBending)(hitNumber1, hitNumber2) + (*fHitWeightsBending)(hitNumber2, hitNumber1)) * dY[hitNumber1] * dY[hitNumber2];
      }
      chi2 += ((*fHitWeightsNonBending)(hitNumber1, hitNumber1) * dX[hitNumber1] * dX[hitNumber1]) +
	      ((*fHitWeightsBending)(hitNumber1, hitNumber1) * dY[hitNumber1] * dY[hitNumber1]);
    }
    delete [] dX;
    delete [] dY;
    
  } else {
    
    AliMUONHitForRec *hitForRec;
    Double_t dX, dY;
    for (Int_t hitNumber = 0; hitNumber < fNTrackHits ; hitNumber++) { 
      trackParamAtHit = (AliMUONTrackParam*) fTrackParamAtHit->UncheckedAt(hitNumber);
      hitForRec = trackParamAtHit->GetHitForRecPtr();
      dX = hitForRec->GetNonBendingCoor() - trackParamAtHit->GetNonBendingCoor();
      dY = hitForRec->GetBendingCoor() - trackParamAtHit->GetBendingCoor();
      chi2 += dX * dX / hitForRec->GetNonBendingReso2() + dY * dY / hitForRec->GetBendingReso2();
    }
    
  }
  
  return chi2;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrack::ComputeHitWeights(TMatrixD* mcsCovariances)
{
  /// Compute the weight matrices of the attached hits, in non bending and bending direction,
  /// accounting for multiple scattering correlations and hits resolution
  /// - Use the provided MCS covariance matrix if any (otherwise build it temporarily)
  /// - Assume that track parameters at each hit are corrects
  /// - Return kFALSE if computation failed
  
  // Alocate memory
  if (!fHitWeightsNonBending) fHitWeightsNonBending = new TMatrixD(fNTrackHits,fNTrackHits);
  if (!fHitWeightsBending) fHitWeightsBending = new TMatrixD(fNTrackHits,fNTrackHits);
  
  // Check hits (if the first one exist, assume that the other ones exit too!)
  if (!((AliMUONTrackParam*) fTrackParamAtHit->First())->GetHitForRecPtr()) {
    AliWarning("hit is missing");
    fHitWeightsNonBending->ResizeTo(0,0);
    fHitWeightsBending->ResizeTo(0,0);
    return kFALSE;
  }
  
  // Compute weights matrices
  if (!ComputeHitWeights(*fHitWeightsNonBending, *fHitWeightsBending, mcsCovariances)) return kFALSE;
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrack::ComputeHitWeights(TMatrixD& hitWeightsNB, TMatrixD& hitWeightsB, TMatrixD* mcsCovariances, AliMUONHitForRec* discardedHit) const
{
  /// Compute the weight matrices, in non bending and bending direction,
  /// of the other attached hits assuming the discarded one does not exist
  /// accounting for multiple scattering correlations and hits resolution
  /// - Use the provided MCS covariance matrix if any (otherwise build it temporarily)
  /// - Return kFALSE if computation failed
  
  // Check MCS covariance matrix and recompute it if need
  Bool_t deleteMCSCov = kFALSE;
  if (!mcsCovariances) {
    mcsCovariances = new TMatrixD(fNTrackHits,fNTrackHits);
    deleteMCSCov = kTRUE;
    ComputeMCSCovariances(*mcsCovariances);
  }
  
  // Resize the weights matrices; alocate memory
  if (discardedHit) {
    hitWeightsNB.ResizeTo(fNTrackHits-1,fNTrackHits-1);
    hitWeightsB.ResizeTo(fNTrackHits-1,fNTrackHits-1);
  } else {
    hitWeightsNB.ResizeTo(fNTrackHits,fNTrackHits);
    hitWeightsB.ResizeTo(fNTrackHits,fNTrackHits);
  }
  
  // Define variables
  AliMUONHitForRec *hitForRec1, *hitForRec2;
  Int_t currentHitNumber1, currentHitNumber2;
  
  // Compute the covariance matrices
  currentHitNumber1 = 0;
  for (Int_t hitNumber1 = 0; hitNumber1 < fNTrackHits; hitNumber1++) { 
    hitForRec1 = ((AliMUONTrackParam*) fTrackParamAtHit->UncheckedAt(hitNumber1))->GetHitForRecPtr();
    
    if (hitForRec1 == discardedHit) continue;
    
    // Loop over next hits
    currentHitNumber2 = currentHitNumber1;
    for (Int_t hitNumber2 = hitNumber1; hitNumber2 < fNTrackHits; hitNumber2++) {
      hitForRec2 = ((AliMUONTrackParam*) fTrackParamAtHit->UncheckedAt(hitNumber2))->GetHitForRecPtr();
      
      if (hitForRec2 == discardedHit) continue;
      
      // Fill with MCS covariances
      hitWeightsNB(currentHitNumber1, currentHitNumber2) = (*mcsCovariances)(hitNumber1,hitNumber2);
      
      // Equal contribution from multiple scattering in non bending and bending directions
      hitWeightsB(currentHitNumber1, currentHitNumber2) = hitWeightsNB(currentHitNumber1, currentHitNumber2);
      
      // Add contribution from hit resolution to diagonal element and symmetrize the matrix
      if (currentHitNumber1 == currentHitNumber2) {
	
	// In non bending plane
        hitWeightsNB(currentHitNumber1, currentHitNumber1) += hitForRec1->GetNonBendingReso2();
	// In bending plane
	hitWeightsB(currentHitNumber1, currentHitNumber1) += hitForRec1->GetBendingReso2();
	
      } else {
	
	// In non bending plane
	hitWeightsNB(currentHitNumber2, currentHitNumber1) = hitWeightsNB(currentHitNumber1, currentHitNumber2);
	// In bending plane
	hitWeightsB(currentHitNumber2, currentHitNumber1) = hitWeightsB(currentHitNumber1, currentHitNumber2);
	
      }
      
      currentHitNumber2++;
    }
    
    currentHitNumber1++;
  }
    
  // Inversion of covariance matrices to get the weights
  if (hitWeightsNB.Determinant() != 0 && hitWeightsB.Determinant() != 0) {
    hitWeightsNB.Invert();
    hitWeightsB.Invert();
  } else {
    AliWarning(" Determinant = 0");
    hitWeightsNB.ResizeTo(0,0);
    hitWeightsB.ResizeTo(0,0);
    if(deleteMCSCov) delete mcsCovariances;
    return kFALSE;
  }
  
  if(deleteMCSCov) delete mcsCovariances;
  
  return kTRUE;
  
}

  //__________________________________________________________________________
void AliMUONTrack::ComputeMCSCovariances(TMatrixD& mcsCovariances) const
{
  /// Compute the multiple scattering covariance matrix
  /// (assume that track parameters at each hit are corrects)
  
  // Reset the size of the covariance matrix if needed
  if (mcsCovariances.GetNrows() != fNTrackHits) mcsCovariances.ResizeTo(fNTrackHits,fNTrackHits);
  
  // Define variables
  Int_t nChambers = AliMUONConstants::NTrackingCh();
  AliMUONTrackParam* trackParamAtHit;
  AliMUONHitForRec *hitForRec;
  AliMUONTrackParam extrapTrackParam;
  Int_t currentChamber = 0, expectedChamber = 0, size = 0;
  Double_t *mcsAngle2 = new Double_t[2*nChambers];
  Double_t *zMCS = new Double_t[2*nChambers];
  Int_t *indices = new Int_t[2*fNTrackHits];
  
  // Compute multiple scattering dispersion angle at each chamber
  // and save the z position where it is calculated
  for (Int_t hitNumber = 0; hitNumber < fNTrackHits; hitNumber++) {
    trackParamAtHit = (AliMUONTrackParam*) fTrackParamAtHit->UncheckedAt(hitNumber);
    hitForRec = trackParamAtHit->GetHitForRecPtr();
    
    // look for missing chambers if any
    currentChamber = hitForRec->GetChamberNumber();
    while (currentChamber > expectedChamber) {
      
      // Save the z position where MCS dispersion is calculated
      zMCS[size] = AliMUONConstants::DefaultChamberZ(expectedChamber);
      
      // Do not take into account MCS in chambers prior the first hit
      if (hitNumber > 0) {
        
        // Get track parameters at missing chamber z
        extrapTrackParam = *trackParamAtHit;
        AliMUONTrackExtrap::ExtrapToZ(&extrapTrackParam, zMCS[expectedChamber]);
        
        // Save multiple scattering dispersion angle in missing chamber
        mcsAngle2[size] = AliMUONTrackExtrap::GetMCSAngle2(extrapTrackParam,AliMUONConstants::ChamberThicknessInX0(),1.);
        
      } else mcsAngle2[size] = 0.;
      
      expectedChamber++;
      size++;
    }
    
    // Save z position where MCS dispersion is calculated
    zMCS[size] = trackParamAtHit->GetZ();
    
    // Save multiple scattering dispersion angle in current chamber
    mcsAngle2[size] = AliMUONTrackExtrap::GetMCSAngle2(*trackParamAtHit,AliMUONConstants::ChamberThicknessInX0(),1.);
    
    // Save indice in zMCS array corresponding to the current cluster
    indices[hitNumber] = size;
    
    expectedChamber = currentChamber + 1;
    size++;
  }
  
  // complete array of z if last hit is on the last but one chamber
  if (currentChamber != nChambers-1) zMCS[size++] = AliMUONConstants::DefaultChamberZ(nChambers-1);
  
  // Compute the covariance matrix
  for (Int_t hitNumber1 = 0; hitNumber1 < fNTrackHits; hitNumber1++) { 
    
    for (Int_t hitNumber2 = hitNumber1; hitNumber2 < fNTrackHits; hitNumber2++) {
      
      // Initialization to 0 (diagonal plus upper triangular part)
      mcsCovariances(hitNumber1,hitNumber2) = 0.;
      
      // Compute contribution from multiple scattering in upstream chambers
      for (Int_t k = 0; k < indices[hitNumber1]; k++) { 	
	mcsCovariances(hitNumber1,hitNumber2) += (zMCS[indices[hitNumber1]] - zMCS[k]) * (zMCS[indices[hitNumber2]] - zMCS[k]) * mcsAngle2[k];
      }
      
      // Symetrize the matrix
      mcsCovariances(hitNumber2,hitNumber1) = mcsCovariances(hitNumber1,hitNumber2);
    }
    
  }
    
  delete [] mcsAngle2;
  delete [] zMCS;
  delete [] indices;
  
}

  //__________________________________________________________________________
Int_t AliMUONTrack::HitsInCommon(AliMUONTrack* track) const
{
  /// Returns the number of hits in common between the current track ("this")
  /// and the track pointed to by "track".
  Int_t hitsInCommon = 0;
  AliMUONTrackParam *trackParamAtHit1, *trackParamAtHit2;
  // Loop over hits of first track
  trackParamAtHit1 = (AliMUONTrackParam*) this->fTrackParamAtHit->First();
  while (trackParamAtHit1) {
    // Loop over hits of second track
    trackParamAtHit2 = (AliMUONTrackParam*) track->fTrackParamAtHit->First();
    while (trackParamAtHit2) {
      // Increment "hitsInCommon" if both TrackParamAtHits point to the same HitForRec
      if ((trackParamAtHit1->GetHitForRecPtr()) == (trackParamAtHit2->GetHitForRecPtr())) {
        hitsInCommon++;
	break;
      }
      trackParamAtHit2 = (AliMUONTrackParam*) track->fTrackParamAtHit->After(trackParamAtHit2);
    } // trackParamAtHit2
    trackParamAtHit1 = (AliMUONTrackParam*) this->fTrackParamAtHit->After(trackParamAtHit1);
  } // trackParamAtHit1
  return hitsInCommon;
}

  //__________________________________________________________________________
Double_t AliMUONTrack::GetNormalizedChi2() const
{
  /// return the chi2 value divided by the number of degrees of freedom (or 1.e10 if ndf < 0)
  
  Double_t numberOfDegFree = (2. * fNTrackHits - 5.);
  if (numberOfDegFree > 0.) return fGlobalChi2 / numberOfDegFree;
  else return 1.e10;
}

  //__________________________________________________________________________
Bool_t* AliMUONTrack::CompatibleTrack(AliMUONTrack * track, Double_t sigma2Cut) const
{
  /// Return kTRUE/kFALSE for each chamber if hit is compatible or not 
  TClonesArray *hitArray, *thisHitArray;
  AliMUONHitForRec *hit, *thisHit;
  Int_t chamberNumber;
  Float_t deltaZ;
  Float_t deltaZMax = 1.; // 1 cm
  Float_t chi2 = 0;
  Bool_t *nCompHit = new Bool_t[AliMUONConstants::NTrackingCh()]; 

  for ( Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    nCompHit[ch] = kFALSE;
  }

  thisHitArray = this->GetHitForRecAtHit();

  hitArray =  track->GetHitForRecAtHit();

  for (Int_t iHthis = 0; iHthis < thisHitArray->GetEntriesFast(); iHthis++) {
    thisHit = (AliMUONHitForRec*) thisHitArray->At(iHthis);
    chamberNumber = thisHit->GetChamberNumber();
    if (chamberNumber < 0 || chamberNumber > AliMUONConstants::NTrackingCh()) continue; 
    nCompHit[chamberNumber] = kFALSE;
    for (Int_t iH = 0; iH < hitArray->GetEntriesFast(); iH++) {
      hit = (AliMUONHitForRec*) hitArray->At(iH);
      deltaZ = TMath::Abs(thisHit->GetZ() - hit->GetZ());
      chi2 = thisHit->NormalizedChi2WithHitForRec(hit,sigma2Cut); // set cut to 4 sigmas
      if (chi2 < 3. && deltaZ < deltaZMax) {
	nCompHit[chamberNumber] = kTRUE;
	break;
      }
    }  
  }
  
  return nCompHit;
}

  //__________________________________________________________________________
void AliMUONTrack::RecursiveDump(void) const
{
  /// Recursive dump of AliMUONTrack, i.e. with dump of TrackParamAtHit's and attached HitForRec's
  AliMUONTrackParam *trackParamAtHit;
  AliMUONHitForRec *hitForRec;
  cout << "Recursive dump of Track: " << this << endl;
  // Track
  this->Dump();
  for (Int_t trackHitIndex = 0; trackHitIndex < fNTrackHits; trackHitIndex++) {
    trackParamAtHit = (AliMUONTrackParam*) ((*fTrackParamAtHit)[trackHitIndex]);
    // TrackHit
    cout << "TrackParamAtHit: " << trackParamAtHit << " (index: " << trackHitIndex << ")" << endl;
    trackParamAtHit->Dump();
    hitForRec = trackParamAtHit->GetHitForRecPtr();
    // HitForRec
    cout << "HitForRec: " << hitForRec << endl;
    hitForRec->Dump();
  }
  return;
}
  
//_____________________________________________-
void AliMUONTrack::Print(Option_t*) const
{
  /// Printing Track information 

  cout << "<AliMUONTrack> No.Clusters=" << setw(2)   << GetNTrackHits() << 
      ", Match2Trig=" << setw(1) << GetMatchTrigger()  << 
      ", LoTrgNum=" << setw(3) << GetLoTrgNum()  << 
    ", Chi2-tracking-trigger=" << setw(8) << setprecision(5) <<  GetChi2MatchTrigger();
  cout << Form(" HitTriggerPattern %x",fHitsPatternInTrigCh) << endl;
  GetTrackParamAtHit()->First()->Print("FULL");
}

//__________________________________________________________________________
void AliMUONTrack::SetLocalTrigger(Int_t loCirc, Int_t loStripX, Int_t loStripY, Int_t loDev, Int_t loLpt, Int_t loHpt)
{
  /// pack the local trigger information and store

  if (loCirc < 0) return;

  fLocalTrigger = 0;
  fLocalTrigger += loCirc;
  fLocalTrigger += loStripX << 8;
  fLocalTrigger += loStripY << 13;
  fLocalTrigger += loDev    << 17;
  fLocalTrigger += loLpt    << 22;
  fLocalTrigger += loHpt    << 24;

}

