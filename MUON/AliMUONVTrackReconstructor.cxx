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
/// \class AliMUONVTrackReconstructor
/// Virtual MUON track reconstructor in ALICE (class renamed from AliMUONEventReconstructor)
///
/// This class contains as data:
/// * a pointer to the array of hits to be reconstructed (the event)
/// * a pointer to the array of reconstructed tracks
///
/// It contains as methods, among others:
/// * EventReconstruct to build the muon tracks
/// * EventReconstructTrigger to build the trigger tracks
///
/// Several options and adjustable parameters are available for both Kalman and Original
/// tracking algorithms (hard coded for the moment in AliMUONVTrackReconstructor.cxx):
/// - *fgkSigmaToCutForTracking* : quality cut used to select new clusters to be
///   attached to the track candidate and to select good tracks.
/// - *fgkMakeTrackCandidatesFast* : if this flag is set to 'true', the track candidates
///   are made assuming linear propagation between stations 4 and 5.
/// - *fgkTrackAllTracks* : according to the value of this flag, in case that several
///   new clusters pass the quality cut, either we consider all the possibilities
///   (duplicating tracks) or we attach only the best cluster.
/// - *fgkRecoverTracks* : if this flag is set to 'true', we try to recover the tracks
///   lost during the tracking by removing the worst of the 2 clusters attached in the
///   previous station.
/// - *fgkImproveTracks* : if this flag is set to 'true', we try to improve the quality
///   of the tracks at the end of the tracking by removing clusters that do not pass
///   new quality cut (within the limit that we must keep at least one cluster per
///   the station).
/// - *fgkSigmaToCutForImprovement* : quality cut used when we try to improve the
///   quality of the tracks.
///
///  \author Philippe Pillot
//-----------------------------------------------------------------------------

#include "AliMUONVTrackReconstructor.h"

#include "AliMUONConstants.h"
#include "AliMUONHitForRec.h"
#include "AliMUONObjectPair.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackHitPattern.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMpDEManager.h"

#include "AliLog.h"
#include "AliCodeTimer.h"
#include "AliTracker.h"

#include <TClonesArray.h>
#include <TMath.h>
#include <TMatrixD.h>

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONVTrackReconstructor) // Class implementation in ROOT context
/// \endcond


//************* Defaults parameters for reconstruction
const Double_t AliMUONVTrackReconstructor::fgkDefaultMinBendingMomentum = 3.0;
const Double_t AliMUONVTrackReconstructor::fgkDefaultMaxBendingMomentum = 3000.0;
const Double_t AliMUONVTrackReconstructor::fgkDefaultMaxNormChi2MatchTrigger = 16.0;
const Double_t AliMUONVTrackReconstructor::fgkMaxTrackingDistanceBending    = 2.;
const Double_t AliMUONVTrackReconstructor::fgkMaxTrackingDistanceNonBending = 2.;

const Double_t AliMUONVTrackReconstructor::fgkSigmaToCutForTracking = 6.0;
const Bool_t   AliMUONVTrackReconstructor::fgkMakeTrackCandidatesFast = kFALSE;
const Bool_t   AliMUONVTrackReconstructor::fgkTrackAllTracks = kTRUE;
const Bool_t   AliMUONVTrackReconstructor::fgkRecoverTracks = kTRUE;
const Bool_t   AliMUONVTrackReconstructor::fgkComplementTracks = kTRUE;
const Bool_t   AliMUONVTrackReconstructor::fgkImproveTracks = kTRUE;
const Double_t AliMUONVTrackReconstructor::fgkSigmaToCutForImprovement = 5.0;


  //__________________________________________________________________________
AliMUONVTrackReconstructor::AliMUONVTrackReconstructor()
  : TObject(),
    fMinBendingMomentum(fgkDefaultMinBendingMomentum),
    fMaxBendingMomentum(fgkDefaultMaxBendingMomentum),
    fMaxNormChi2MatchTrigger(fgkDefaultMaxNormChi2MatchTrigger),
    fHitsForRecPtr(0x0),
    fNHitsForRec(0),
    fNHitsForRecPerChamber(0x0),
    fIndexOfFirstHitForRecPerChamber(0x0),
    fRecTracksPtr(0x0),
    fNRecTracks(0)
{
  /// Constructor for class AliMUONVTrackReconstructor
  fNHitsForRecPerChamber = new Int_t[AliMUONConstants::NTrackingCh()];
  fIndexOfFirstHitForRecPerChamber = new Int_t[AliMUONConstants::NTrackingCh()];

  // Memory allocation for the TClonesArray of hits for reconstruction
  // Is 10000 the right size ????
  fHitsForRecPtr = new TClonesArray("AliMUONHitForRec", 10000);
  
  // Memory allocation for the TClonesArray of reconstructed tracks
  fRecTracksPtr = new TClonesArray("AliMUONTrack", 10);
  
  // set the magnetic field for track extrapolations
  const AliMagF* kField = AliTracker::GetFieldMap();
  if (!kField) AliFatal("No field available");
  AliMUONTrackExtrap::SetField(kField);
}

  //__________________________________________________________________________
AliMUONVTrackReconstructor::~AliMUONVTrackReconstructor()
{
  /// Destructor for class AliMUONVTrackReconstructor
  delete [] fNHitsForRecPerChamber;
  delete [] fIndexOfFirstHitForRecPerChamber;
  delete fHitsForRecPtr;
  delete fRecTracksPtr;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::ResetHitsForRec()
{
  /// To reset the TClonesArray of HitsForRec,
  /// and the number of HitForRec and the index of the first HitForRec per chamber
  if (fHitsForRecPtr) fHitsForRecPtr->Clear("C");
  fNHitsForRec = 0;
  for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    fNHitsForRecPerChamber[ch] = 0;
    fIndexOfFirstHitForRecPerChamber[ch] = 0;
  }
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::ResetTracks()
{
  /// To reset the TClonesArray of reconstructed tracks
  if (fRecTracksPtr) fRecTracksPtr->Clear("C");
  fNRecTracks = 0;
  return;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::EventReconstruct(const AliMUONVClusterStore& clusterStore, AliMUONVTrackStore& trackStore)
{
  /// To reconstruct one event
  AliDebug(1,"");
  AliCodeTimerAuto("");
  
  ResetTracks();
  ResetHitsForRec();
  AddHitsForRecFromRawClusters(clusterStore);
  MakeTracks();

  // Add tracks to MUON data container 
  for (Int_t i=0; i<fNRecTracks; ++i) 
  {
    AliMUONTrack * track = (AliMUONTrack*) fRecTracksPtr->At(i);
    trackStore.Add(*track);
  }
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::AddHitsForRecFromRawClusters(const AliMUONVClusterStore& clusterStore)
{
  /// Build internal array of hit for rec from clusterStore
  
  AliMUONVCluster* clus(0x0);
  Int_t iclus(0);
  
  TIter next(clusterStore.CreateIterator());
  
  while ( ( clus = static_cast<AliMUONVCluster*>(next()) ) )
  {
    // new AliMUONHitForRec from raw cluster
    // and increment number of AliMUONHitForRec's (total and in chamber)
    AliMUONHitForRec* hitForRec = new ((*fHitsForRecPtr)[fNHitsForRec]) AliMUONHitForRec(clus);
    fNHitsForRec++;
    // more information into HitForRec
    hitForRec->SetNonBendingReso2(clus->GetErrX2());
    hitForRec->SetBendingReso2(clus->GetErrY2());
    //  original raw cluster
    Int_t ch = AliMpDEManager::GetChamberId(clus->GetDetElemId());
    hitForRec->SetChamberNumber(ch);
    hitForRec->SetHitNumber(iclus);
    // Z coordinate of the raw cluster (cm)
    hitForRec->SetZ(clus->GetZ());
    if (AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 3) {
      cout << "Chamber " << ch <<" raw cluster  " << iclus << " : " << endl;
      clus->Print("full");
      cout << "AliMUONHitForRec number (1...): " << fNHitsForRec << endl;
      hitForRec->Print("full");
    }
    ++iclus;
  } // end of chamber loop
  
  SortHitsForRecWithIncreasingChamber(); 
  
  AliDebug(1,"End of AddHitsForRecFromRawClusters");
  
  if (AliLog::GetGlobalDebugLevel() > 0) 
  {
    AliDebug(1, Form("NHitsForRec: %d",fNHitsForRec));
    for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) 
    {
      AliDebug(1, Form("Chamber(0...): %d",ch));
      AliDebug(1, Form("NHitsForRec: %d", fNHitsForRecPerChamber[ch]));
      AliDebug(1, Form("Index(first HitForRec): %d", fIndexOfFirstHitForRecPerChamber[ch]));
      for (Int_t hit = fIndexOfFirstHitForRecPerChamber[ch];
           hit < fIndexOfFirstHitForRecPerChamber[ch] + fNHitsForRecPerChamber[ch];
           hit++) {
        AliDebug(1, Form("HitForRec index(0...): %d",hit));
        ((*fHitsForRecPtr)[hit])->Dump();
      }
    }
  }
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::SortHitsForRecWithIncreasingChamber()
{
  /// Sort HitsForRec's in increasing order with respect to chamber number.
  /// Uses the function "Compare".
  /// Update the information for HitsForRec per chamber too.
  Int_t ch, nhits, prevch;
  fHitsForRecPtr->Sort();
  for (ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    fNHitsForRecPerChamber[ch] = 0;
    fIndexOfFirstHitForRecPerChamber[ch] = 0;
  }
  prevch = 0; // previous chamber
  nhits = 0; // number of hits in current chamber
  // Loop over HitsForRec
  for (Int_t hit = 0; hit < fNHitsForRec; hit++) {
    // chamber number (0...)
    ch = ((AliMUONHitForRec*)  ((*fHitsForRecPtr)[hit]))->GetChamberNumber();
    // increment number of hits in current chamber
    (fNHitsForRecPerChamber[ch])++;
    // update index of first HitForRec in current chamber
    // if chamber number different from previous one
    if (ch != prevch) {
      fIndexOfFirstHitForRecPerChamber[ch] = hit;
      prevch = ch;
    }
  }
  return;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::MakeTracks()
{
  /// To make the tracks from the list of segments and points in all stations
  AliDebug(1,"Enter MakeTracks");
  // Look for candidates from at least 3 aligned points in stations(1..) 4 and 5
  MakeTrackCandidates();
  if (fRecTracksPtr->GetEntriesFast() == 0) return;
  // Follow tracks in stations(1..) 3, 2 and 1
  FollowTracks();
  // Complement the reconstructed tracks
  if (fgkComplementTracks) ComplementTracks();
  // Improve the reconstructed tracks
  if (fgkImproveTracks) ImproveTracks();
  // Remove double tracks
  RemoveDoubleTracks();
  // Fill AliMUONTrack data members
  Finalize();
}

  //__________________________________________________________________________
TClonesArray* AliMUONVTrackReconstructor::MakeSegmentsInStation(Int_t station)
{
  /// To make the list of segments in station(0..) "Station" from the list of hits to be reconstructed.
  /// Return a new TClonesArray of segments.
  /// It is the responsibility of the user to delete it afterward.
  AliDebug(1,Form("Enter MakeSegmentsPerStation (0...) %d",station));
  
  AliMUONHitForRec *hit1Ptr, *hit2Ptr;
  AliMUONObjectPair *segment;
  Double_t bendingSlope = 0, impactParam = 0., bendingMomentum = 0.; // to avoid compilation warning
                                                                     // first and second chambers (0...) in the station
  Int_t ch1 = 2 * station;
  Int_t ch2 = ch1 + 1;
  
  // list of segments
  TClonesArray *segments = new TClonesArray("AliMUONObjectPair", fNHitsForRecPerChamber[ch2]);
  segments->SetOwner(kTRUE);
  
  // Loop over HitForRec's in the first chamber of the station
  for (Int_t hit1 = fIndexOfFirstHitForRecPerChamber[ch1];
       hit1 < fIndexOfFirstHitForRecPerChamber[ch1] + fNHitsForRecPerChamber[ch1];
       hit1++) 
  {
    // pointer to the HitForRec
    hit1Ptr = (AliMUONHitForRec*) ((*fHitsForRecPtr)[hit1]);
    // Loop over HitsForRec's in the second chamber of the station
    for (Int_t hit2 = fIndexOfFirstHitForRecPerChamber[ch2];
         hit2 < fIndexOfFirstHitForRecPerChamber[ch2] + fNHitsForRecPerChamber[ch2];
         hit2++) 
    {
      // pointer to the HitForRec
      hit2Ptr = (AliMUONHitForRec*) ((*fHitsForRecPtr)[hit2]);
      if ( hit1Ptr->GetZ() - hit2Ptr->GetZ() != 0. ) 
      {
        // bending slope
        bendingSlope = (hit1Ptr->GetBendingCoor() - hit2Ptr->GetBendingCoor()) / (hit1Ptr->GetZ() - hit2Ptr->GetZ());
        // impact parameter
        impactParam = hit1Ptr->GetBendingCoor() - hit1Ptr->GetZ() * bendingSlope;
        // absolute value of bending momentum
        bendingMomentum = TMath::Abs(AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(impactParam));
      } else 
      {
        AliWarning("hit1Ptr->GetZ() = hit2Ptr->GetZ(): no segment created");
        continue;
      }   
      // check for bending momentum within tolerances
      if ((bendingMomentum < fMaxBendingMomentum) && (bendingMomentum > fMinBendingMomentum)) 
      {
        // make new segment
        segment = new ((*segments)[segments->GetLast()+1]) AliMUONObjectPair(hit1Ptr, hit2Ptr, kFALSE, kFALSE);
        if (AliLog::GetGlobalDebugLevel() > 1) {
          cout << "segmentIndex(0...): " << segments->GetLast() << endl;
          segment->Dump();
          cout << "HitForRec in first chamber" << endl;
          hit1Ptr->Dump();
          cout << "HitForRec in second chamber" << endl;
          hit2Ptr->Dump();
        }
      }
    } //for (Int_t hit2
  } // for (Int_t hit1...
  AliDebug(1,Form("Station: %d  NSegments:  %d ", station, segments->GetEntriesFast()));
  return segments;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::RemoveIdenticalTracks()
{
  /// To remove identical tracks:
  /// Tracks are considered identical if they have all their hits in common.
  /// One keeps the track with the larger number of hits if need be
  AliMUONTrack *track1, *track2, *trackToRemove;
  Int_t hitsInCommon, nHits1, nHits2;
  Bool_t removedTrack1;
  // Loop over first track of the pair
  track1 = (AliMUONTrack*) fRecTracksPtr->First();
  while (track1) {
    removedTrack1 = kFALSE;
    nHits1 = track1->GetNTrackHits();
    // Loop over second track of the pair
    track2 = (AliMUONTrack*) fRecTracksPtr->After(track1);
    while (track2) {
      nHits2 = track2->GetNTrackHits();
      // number of hits in common between two tracks
      hitsInCommon = track1->HitsInCommon(track2);
      // check for identical tracks
      if ((hitsInCommon == nHits1) || (hitsInCommon == nHits2)) {
        // decide which track to remove
        if (nHits2 > nHits1) {
	  // remove track1 and continue the first loop with the track next to track1
	  trackToRemove = track1;
	  track1 = (AliMUONTrack*) fRecTracksPtr->After(track1);
          fRecTracksPtr->Remove(trackToRemove);
	  fRecTracksPtr->Compress(); // this is essential to retrieve the TClonesArray afterwards
	  fNRecTracks--;
	  removedTrack1 = kTRUE;
	  break;
	} else {
	  // remove track2 and continue the second loop with the track next to track2
	  trackToRemove = track2;
	  track2 = (AliMUONTrack*) fRecTracksPtr->After(track2);
	  fRecTracksPtr->Remove(trackToRemove);
	  fRecTracksPtr->Compress(); // this is essential to retrieve the TClonesArray afterwards
	  fNRecTracks--;
        }
      } else track2 = (AliMUONTrack*) fRecTracksPtr->After(track2);
    } // track2
    if (removedTrack1) continue;
    track1 = (AliMUONTrack*) fRecTracksPtr->After(track1);
  } // track1
  return;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::RemoveDoubleTracks()
{
  /// To remove double tracks:
  /// Tracks are considered identical if more than half of the hits of the track
  /// which has the smaller number of hits are in common with the other track.
  /// Among two identical tracks, one keeps the track with the larger number of hits
  /// or, if these numbers are equal, the track with the minimum chi2.
  AliMUONTrack *track1, *track2, *trackToRemove;
  Int_t hitsInCommon, nHits1, nHits2;
  Bool_t removedTrack1;
  // Loop over first track of the pair
  track1 = (AliMUONTrack*) fRecTracksPtr->First();
  while (track1) {
    removedTrack1 = kFALSE;
    nHits1 = track1->GetNTrackHits();
    // Loop over second track of the pair
    track2 = (AliMUONTrack*) fRecTracksPtr->After(track1);
    while (track2) {
      nHits2 = track2->GetNTrackHits();
      // number of hits in common between two tracks
      hitsInCommon = track1->HitsInCommon(track2);
      // check for identical tracks
      if (((nHits1 < nHits2) && (2 * hitsInCommon > nHits1)) || (2 * hitsInCommon > nHits2)) {
        // decide which track to remove
        if ((nHits1 > nHits2) || ((nHits1 == nHits2) && (track1->GetFitFMin() <= track2->GetFitFMin()))) {
	  // remove track2 and continue the second loop with the track next to track2
	  trackToRemove = track2;
	  track2 = (AliMUONTrack*) fRecTracksPtr->After(track2);
	  fRecTracksPtr->Remove(trackToRemove);
	  fRecTracksPtr->Compress(); // this is essential to retrieve the TClonesArray afterwards
	  fNRecTracks--;
        } else {
	  // else remove track1 and continue the first loop with the track next to track1
	  trackToRemove = track1;
	  track1 = (AliMUONTrack*) fRecTracksPtr->After(track1);
          fRecTracksPtr->Remove(trackToRemove);
	  fRecTracksPtr->Compress(); // this is essential to retrieve the TClonesArray afterwards
	  fNRecTracks--;
	  removedTrack1 = kTRUE;
	  break;
        }
      } else track2 = (AliMUONTrack*) fRecTracksPtr->After(track2);
    } // track2
    if (removedTrack1) continue;
    track1 = (AliMUONTrack*) fRecTracksPtr->After(track1);
  } // track1
  return;
}

  //__________________________________________________________________________
Double_t AliMUONVTrackReconstructor::TryOneHitForRec(const AliMUONTrackParam &trackParam, AliMUONHitForRec* hitForRec,
						     AliMUONTrackParam &trackParamAtHit, Bool_t updatePropagator)
{
/// Test the compatibility between the track and the hitForRec (using trackParam's covariance matrix):
/// return the corresponding Chi2
/// return trackParamAtHit
  
  // extrapolate track parameters and covariances at the z position of the tested hit
  trackParamAtHit = trackParam;
  AliMUONTrackExtrap::ExtrapToZCov(&trackParamAtHit, hitForRec->GetZ(), updatePropagator);
  
  // set pointer to hit into trackParamAtHit
  trackParamAtHit.SetHitForRecPtr(hitForRec);
  
  // Set differences between trackParam and hitForRec in the bending and non bending directions
  Double_t dX = hitForRec->GetNonBendingCoor() - trackParamAtHit.GetNonBendingCoor();
  Double_t dY = hitForRec->GetBendingCoor() - trackParamAtHit.GetBendingCoor();
  
  // Calculate errors and covariances
  const TMatrixD& kParamCov = trackParamAtHit.GetCovariances();
  Double_t sigma2X = kParamCov(0,0) + hitForRec->GetNonBendingReso2();
  Double_t sigma2Y = kParamCov(2,2) + hitForRec->GetBendingReso2();
  
  // Compute chi2
  return dX * dX / sigma2X + dY * dY / sigma2Y;
  
}

  //__________________________________________________________________________
Bool_t AliMUONVTrackReconstructor::TryOneHitForRecFast(const AliMUONTrackParam &trackParam, AliMUONHitForRec* hitForRec)
{
/// Test the compatibility between the track and the hitForRec within a wide fix window
/// assuming linear propagation of the track:
/// return kTRUE if they are compatibles
  
  Double_t dZ = hitForRec->GetZ() - trackParam.GetZ();
  Double_t dX = hitForRec->GetNonBendingCoor() - (trackParam.GetNonBendingCoor() + trackParam.GetNonBendingSlope() * dZ);
  Double_t dY = hitForRec->GetBendingCoor() - (trackParam.GetBendingCoor() + trackParam.GetBendingSlope() * dZ);
  
  if (TMath::Abs(dX) > fgkMaxTrackingDistanceNonBending || TMath::Abs(dY) > fgkMaxTrackingDistanceBending) return kFALSE;
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Double_t AliMUONVTrackReconstructor::TryTwoHitForRecFast(const AliMUONTrackParam &trackParamAtHit1, AliMUONHitForRec* hitForRec2,
							 AliMUONTrackParam &trackParamAtHit2)
{
/// Test the compatibility between the track and the 2 hitForRec together (using trackParam's covariance matrix)
/// assuming linear propagation between the two hits:
/// return the corresponding Chi2 accounting for covariances between the 2 hitForRec
/// return trackParamAtHit1 & 2
  
  // extrapolate linearly track parameters at the z position of the second hit
  trackParamAtHit2 = trackParamAtHit1;
  AliMUONTrackExtrap::LinearExtrapToZ(&trackParamAtHit2, hitForRec2->GetZ());
  
  // set pointer to hit2 into trackParamAtHit2
  trackParamAtHit2.SetHitForRecPtr(hitForRec2);
  
  // Set differences between track and hitForRec in the bending and non bending directions
  AliMUONHitForRec* hitForRec1 = trackParamAtHit1.GetHitForRecPtr();
  Double_t dX1 = hitForRec1->GetNonBendingCoor() - trackParamAtHit1.GetNonBendingCoor();
  Double_t dX2 = hitForRec2->GetNonBendingCoor() - trackParamAtHit2.GetNonBendingCoor();
  Double_t dY1 = hitForRec1->GetBendingCoor() - trackParamAtHit1.GetBendingCoor();
  Double_t dY2 = hitForRec2->GetBendingCoor() - trackParamAtHit2.GetBendingCoor();
  
  // Calculate errors and covariances
  const TMatrixD& kParamCov1 = trackParamAtHit1.GetCovariances();
  const TMatrixD& kParamCov2 = trackParamAtHit2.GetCovariances();
  Double_t dZ = trackParamAtHit2.GetZ() - trackParamAtHit1.GetZ();
  Double_t sigma2X1 = kParamCov1(0,0) + hitForRec1->GetNonBendingReso2();
  Double_t sigma2X2 = kParamCov2(0,0) + hitForRec2->GetNonBendingReso2();
  Double_t covX1X2  = kParamCov1(0,0) + dZ * kParamCov1(0,1);
  Double_t sigma2Y1 = kParamCov1(2,2) + hitForRec1->GetBendingReso2();
  Double_t sigma2Y2 = kParamCov2(2,2) + hitForRec2->GetBendingReso2();
  Double_t covY1Y2  = kParamCov1(2,2) + dZ * kParamCov1(2,3);
  
  // Compute chi2
  Double_t detX = sigma2X1 * sigma2X2 - covX1X2 * covX1X2;
  Double_t detY = sigma2Y1 * sigma2Y2 - covY1Y2 * covY1Y2;
  if (detX == 0. || detY == 0.) return 1.e10;
  return   (dX1 * dX1 * sigma2X2 + dX2 * dX2 * sigma2X1 - 2. * dX1 * dX2 * covX1X2) / detX
  	 + (dY1 * dY1 * sigma2Y2 + dY2 * dY2 * sigma2Y1 - 2. * dY1 * dY2 * covY1Y2) / detY;
  
}

  //__________________________________________________________________________
Bool_t AliMUONVTrackReconstructor::FollowLinearTrackInStation(AliMUONTrack &trackCandidate, Int_t nextStation)
{
  /// Follow trackCandidate in station(0..) nextStation assuming linear propagation, and search for compatible HitForRec(s)
  /// Keep all possibilities or only the best one(s) according to the flag fgkTrackAllTracks:
  /// kTRUE:  duplicate "trackCandidate" if there are several possibilities and add the new tracks at the end of
  ///         fRecTracksPtr to avoid conficts with other track candidates at this current stage of the tracking procedure.
  ///         Remove the obsolete "trackCandidate" at the end.
  /// kFALSE: add only the best hit(s) to the "trackCandidate". Try to add a couple of hits in priority.
  /// return kTRUE if new hits have been found (otherwise return kFALSE)
  AliDebug(1,Form("Enter FollowLinearTrackInStation(1..) %d", nextStation+1));
  
  // Order the chamber according to the propagation direction (tracking starts with chamber 2):
  // - nextStation == station(1...) 5 => forward propagation
  // - nextStation < station(1...) 5 => backward propagation
  Int_t ch1, ch2;
  if (nextStation==4) {
    ch1 = 2*nextStation+1;
    ch2 = 2*nextStation;
  } else {
    ch1 = 2*nextStation;
    ch2 = 2*nextStation+1;
  }
  
  Double_t chi2WithOneHitForRec = 1.e10;
  Double_t chi2WithTwoHitForRec = 1.e10;
  Double_t maxChi2WithOneHitForRec = 2. * fgkSigmaToCutForTracking * fgkSigmaToCutForTracking; // 2 because 2 quantities in chi2
  Double_t maxChi2WithTwoHitForRec = 4. * fgkSigmaToCutForTracking * fgkSigmaToCutForTracking; // 4 because 4 quantities in chi2
  Double_t bestChi2WithOneHitForRec = maxChi2WithOneHitForRec;
  Double_t bestChi2WithTwoHitForRec = maxChi2WithTwoHitForRec;
  Bool_t foundOneHit = kFALSE;
  Bool_t foundTwoHits = kFALSE;
  AliMUONTrack *newTrack = 0x0;
  AliMUONHitForRec *hitForRecCh1, *hitForRecCh2;
  AliMUONTrackParam extrapTrackParamAtHit1;
  AliMUONTrackParam extrapTrackParamAtHit2;
  AliMUONTrackParam bestTrackParamAtHit1;
  AliMUONTrackParam bestTrackParamAtHit2;
  Bool_t *hitForRecCh1Used = new Bool_t[fNHitsForRecPerChamber[ch1]];
  for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) hitForRecCh1Used[hit1] = kFALSE;
  
  // Get track parameters
  AliMUONTrackParam trackParam(*(AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->First());
  
  // Add MCS effect
  AliMUONTrackExtrap::AddMCSEffect(&trackParam,AliMUONConstants::ChamberThicknessInX0(),1.);
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "FollowLinearTrackInStation: look for hits in chamber(1..): " << ch2+1 << endl;
  }
  
  // look for candidates in chamber 2 
  for (Int_t hit2 = 0; hit2 < fNHitsForRecPerChamber[ch2]; hit2++) {
    
    hitForRecCh2 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch2]+hit2);
    
    // try to add the current hit fast
    if (!TryOneHitForRecFast(trackParam, hitForRecCh2)) continue;
    
    // try to add the current hit accuratly
    extrapTrackParamAtHit2 = trackParam;
    AliMUONTrackExtrap::LinearExtrapToZ(&extrapTrackParamAtHit2, hitForRecCh2->GetZ());
    chi2WithOneHitForRec = TryOneHitForRec(extrapTrackParamAtHit2, hitForRecCh2, extrapTrackParamAtHit2);
    
    // if good chi2 then try to attach a hitForRec in the other chamber too
    if (chi2WithOneHitForRec < maxChi2WithOneHitForRec) {
      Bool_t foundSecondHit = kFALSE;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowLinearTrackInStation: found one hit in chamber(1..): " << ch2+1
	     << " (Chi2 = " << chi2WithOneHitForRec << ")" << endl;
        cout << "                      look for second hits in chamber(1..): " << ch1+1 << " ..." << endl;
      }
      
      // add MCS effect
      AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtHit2,AliMUONConstants::ChamberThicknessInX0(),1.);
      
      for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) {
        
	hitForRecCh1 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch1]+hit1);
	
        // try to add the current hit fast
        if (!TryOneHitForRecFast(extrapTrackParamAtHit2, hitForRecCh1)) continue;
	
	// try to add the current hit in addition to the one found on the previous chamber
	chi2WithTwoHitForRec = TryTwoHitForRecFast(extrapTrackParamAtHit2, hitForRecCh1, extrapTrackParamAtHit1);
        
	// if good chi2 then consider to add the 2 hitForRec to the "trackCandidate"
	if (chi2WithTwoHitForRec < maxChi2WithTwoHitForRec) {
	  foundSecondHit = kTRUE;
	  foundTwoHits = kTRUE;
          
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowLinearTrackInStation: found one hit in chamber(1..): " << ch1+1
	  	 << " (Chi2 = " << chi2WithTwoHitForRec << ")" << endl;
	  }
	  
	  if (fgkTrackAllTracks) {
	    // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
            newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	    extrapTrackParamAtHit1.SetRemovable(kTRUE);
	    newTrack->AddTrackParamAtHit(&extrapTrackParamAtHit1,hitForRecCh1);
	    extrapTrackParamAtHit2.SetRemovable(kTRUE);
	    newTrack->AddTrackParamAtHit(&extrapTrackParamAtHit2,hitForRecCh2);
	    newTrack->GetTrackParamAtHit()->Sort();
	    fNRecTracks++;
	    
	    // Tag hitForRecCh1 as used
	    hitForRecCh1Used[hit1] = kTRUE;
	    
	    // Printout for debuging
	    if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	      cout << "FollowLinearTrackInStation: added two hits in station(1..): " << nextStation+1 << endl;
	      if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	    }
	    
          } else if (chi2WithTwoHitForRec < bestChi2WithTwoHitForRec) {
	    // keep track of the best couple of hits
	    bestChi2WithTwoHitForRec = chi2WithTwoHitForRec;
	    bestTrackParamAtHit1 = extrapTrackParamAtHit1;
	    bestTrackParamAtHit2 = extrapTrackParamAtHit2;
          }
	  
	}
	
      }
      
      // if no hitForRecCh1 found then consider to add hitForRecCh2 only
      if (!foundSecondHit) {
	foundOneHit = kTRUE;
        
	if (fgkTrackAllTracks) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
          newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  extrapTrackParamAtHit2.SetRemovable(kFALSE);
	  newTrack->AddTrackParamAtHit(&extrapTrackParamAtHit2,hitForRecCh2);
	  newTrack->GetTrackParamAtHit()->Sort();
	  fNRecTracks++;
	  
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowLinearTrackInStation: added one hit in chamber(1..): " << ch2+1 << endl;
	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	  }
	  
	} else if (!foundTwoHits && chi2WithOneHitForRec < bestChi2WithOneHitForRec) {
	  // keep track of the best single hitForRec except if a couple of hits has already been found
	  bestChi2WithOneHitForRec = chi2WithOneHitForRec;
	  bestTrackParamAtHit1 = extrapTrackParamAtHit2;
        }
	
      }
      
    }
    
  }
  
  // look for candidates in chamber 1 not already attached to a track
  // if we want to keep all possible tracks or if no good couple of hitForRec has been found
  if (fgkTrackAllTracks || !foundTwoHits) {
    
    // Printout for debuging
    if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
      cout << "FollowLinearTrackInStation: look for single hits in chamber(1..): " << ch1+1 << endl;
    }
    
    //Extrapolate trackCandidate to chamber "ch2"
    AliMUONTrackExtrap::LinearExtrapToZ(&trackParam, AliMUONConstants::DefaultChamberZ(ch2));
    
    // add MCS effect for next step
    AliMUONTrackExtrap::AddMCSEffect(&trackParam,AliMUONConstants::ChamberThicknessInX0(),1.);
      
    for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) {
      
      hitForRecCh1 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch1]+hit1);
      
      if (hitForRecCh1Used[hit1]) continue; // Skip hitForRec already used
      
      // try to add the current hit fast
      if (!TryOneHitForRecFast(trackParam, hitForRecCh1)) continue;
	
      // try to add the current hit accuratly
      extrapTrackParamAtHit1 = trackParam;
      AliMUONTrackExtrap::LinearExtrapToZ(&extrapTrackParamAtHit1, hitForRecCh1->GetZ());
      chi2WithOneHitForRec = TryOneHitForRec(extrapTrackParamAtHit1, hitForRecCh1, extrapTrackParamAtHit1);
    
      // if good chi2 then consider to add hitForRecCh1
      // We do not try to attach a hitForRec in the other chamber too since it has already been done above
      if (chi2WithOneHitForRec < maxChi2WithOneHitForRec) {
	foundOneHit = kTRUE;
  	
	// Printout for debuging
  	if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	  cout << "FollowLinearTrackInStation: found one hit in chamber(1..): " << ch1+1
  	       << " (Chi2 = " << chi2WithOneHitForRec << ")" << endl;
  	}
	
	if (fgkTrackAllTracks) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
  	  newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  extrapTrackParamAtHit1.SetRemovable(kFALSE);
	  newTrack->AddTrackParamAtHit(&extrapTrackParamAtHit1,hitForRecCh1);
	  newTrack->GetTrackParamAtHit()->Sort();
	  fNRecTracks++;
  	  
	  // Printout for debuging
  	  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	    cout << "FollowLinearTrackInStation: added one hit in chamber(1..): " << ch1+1 << endl;
  	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
  	  }
	  
	} else if (chi2WithOneHitForRec < bestChi2WithOneHitForRec) {
	  // keep track of the best single hitForRec except if a couple of hits has already been found
	  bestChi2WithOneHitForRec = chi2WithOneHitForRec;
	  bestTrackParamAtHit1 = extrapTrackParamAtHit1;
  	}
	
      }
      
    }
    
  }
  
  // fill out the best track if required else clean up the fRecTracksPtr array
  if (!fgkTrackAllTracks) {
    if (foundTwoHits) {
      bestTrackParamAtHit1.SetRemovable(kTRUE);
      trackCandidate.AddTrackParamAtHit(&bestTrackParamAtHit1,bestTrackParamAtHit1.GetHitForRecPtr());
      bestTrackParamAtHit2.SetRemovable(kTRUE);
      trackCandidate.AddTrackParamAtHit(&bestTrackParamAtHit2,bestTrackParamAtHit2.GetHitForRecPtr());
      trackCandidate.GetTrackParamAtHit()->Sort();
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowLinearTrackInStation: added the two best hits in station(1..): " << nextStation+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else if (foundOneHit) {
      bestTrackParamAtHit1.SetRemovable(kFALSE);
      trackCandidate.AddTrackParamAtHit(&bestTrackParamAtHit1,bestTrackParamAtHit1.GetHitForRecPtr());
      trackCandidate.GetTrackParamAtHit()->Sort();
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowLinearTrackInStation: added the best hit in chamber(1..): " << bestTrackParamAtHit1.GetHitForRecPtr()->GetChamberNumber()+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else {
      delete [] hitForRecCh1Used;
      return kFALSE;
    }
    
  } else if (foundOneHit || foundTwoHits) {
    
    // remove obsolete track
    fRecTracksPtr->Remove(&trackCandidate);
    fNRecTracks--;
    
  } else {
    delete [] hitForRecCh1Used;  
    return kFALSE;
  }
  
  delete [] hitForRecCh1Used;
  return kTRUE;
  
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::ValidateTracksWithTrigger(AliMUONVTrackStore& trackStore,
                                                           const AliMUONVTriggerTrackStore& triggerTrackStore,
                                                           const AliMUONVTriggerStore& triggerStore,
                                                           const AliMUONTrackHitPattern& trackHitPattern)
{
  /// Try to match track from tracking system with trigger track
  AliCodeTimerAuto("");
  
  static const Double_t kDistSigma[3]={1,1,0.02}; // sigma of distributions (trigger-track) X,Y,slopeY
  
  Int_t matchTrigger;
  Int_t loTrgNum(-1);
  Double_t distTriggerTrack[3];
  Double_t xTrack, yTrack, ySlopeTrack, chi2MatchTrigger, minChi2MatchTrigger, chi2;

  TIter itTrack(trackStore.CreateIterator());
  AliMUONTrack* track;
  
  while ( ( track = static_cast<AliMUONTrack*>(itTrack()) ) )
  {
    matchTrigger = 0;
    chi2MatchTrigger = 0.;
    loTrgNum = -1;
    Int_t doubleMatch=-1; // Check if track matches 2 trigger tracks
    Double_t doubleChi2 = -1.;
    
    AliMUONTrackParam trackParam(*((AliMUONTrackParam*) (track->GetTrackParamAtHit()->Last())));
    AliMUONTrackExtrap::ExtrapToZ(&trackParam, AliMUONConstants::DefaultChamberZ(10)); // extrap to 1st trigger chamber
    
    xTrack = trackParam.GetNonBendingCoor();
    yTrack = trackParam.GetBendingCoor();
    ySlopeTrack = trackParam.GetBendingSlope();
    minChi2MatchTrigger = 999.;
    
    AliMUONTriggerTrack *triggerTrack;
    TIter itTriggerTrack(triggerTrackStore.CreateIterator());
    while ( ( triggerTrack = static_cast<AliMUONTriggerTrack*>(itTriggerTrack() ) ) )
    {
      distTriggerTrack[0] = (triggerTrack->GetX11()-xTrack)/kDistSigma[0];
      distTriggerTrack[1] = (triggerTrack->GetY11()-yTrack)/kDistSigma[1];
      distTriggerTrack[2] = (TMath::Tan(triggerTrack->GetThetay())-ySlopeTrack)/kDistSigma[2];
      chi2 = 0.;
      for (Int_t iVar = 0; iVar < 3; iVar++) chi2 += distTriggerTrack[iVar]*distTriggerTrack[iVar];
      chi2 /= 3.; // Normalized Chi2: 3 degrees of freedom (X,Y,slopeY)
      if (chi2 < fMaxNormChi2MatchTrigger) 
      {
        Bool_t isDoubleTrack = (TMath::Abs(chi2 - minChi2MatchTrigger)<1.);
        if (chi2 < minChi2MatchTrigger && chi2 < fMaxNormChi2MatchTrigger) 
        {
          if(isDoubleTrack)
          {
            doubleMatch = loTrgNum;
            doubleChi2 = chi2MatchTrigger;
          }
          minChi2MatchTrigger = chi2;
          chi2MatchTrigger = chi2;
          loTrgNum = triggerTrack->GetLoTrgNum();
          AliMUONLocalTrigger* locTrg = triggerStore.FindLocal(loTrgNum);
          matchTrigger=1;
          if(locTrg->LoLpt()>0)matchTrigger=2;
          if(locTrg->LoHpt()>0)matchTrigger=3;
        }
        else if(isDoubleTrack) 
        {
          doubleMatch = triggerTrack->GetLoTrgNum();
          doubleChi2 = chi2;
        }
      }
    }
    if(doubleMatch>=0)
    { // If two trigger tracks match, select the one passing more trigger cuts
      AliDebug(1, Form("Two candidates found: %i and %i",loTrgNum,doubleMatch));
      AliMUONLocalTrigger* locTrg1 = triggerStore.FindLocal(doubleMatch);
      if((locTrg1->LoLpt()>0 && matchTrigger<2) || (locTrg1->LoHpt() && matchTrigger<3))
      {
        if(locTrg1->LoHpt()>0)matchTrigger=3;
        else matchTrigger=2;
        loTrgNum = doubleMatch;
        chi2MatchTrigger=doubleChi2;
      }
    }
    
    track->SetMatchTrigger(matchTrigger);
    track->SetLoTrgNum(loTrgNum);
    track->SetChi2MatchTrigger(chi2MatchTrigger);
    
    AliMUONLocalTrigger* locTrg = static_cast<AliMUONLocalTrigger*>(triggerStore.FindLocal(loTrgNum));
    
    if (locTrg)
    {    
      track->SetLocalTrigger(locTrg->LoCircuit(),
                             locTrg->LoStripX(),
                             locTrg->LoStripY(),
                             locTrg->LoDev(),
                             locTrg->LoLpt(),
                             locTrg->LoHpt());
    }    
  }  
  
  trackHitPattern.GetHitPattern(trackStore,triggerStore);
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::EventReconstructTrigger(const AliMUONTriggerCircuit& circuit,
                                                         const AliMUONVTriggerStore& triggerStore,
                                                         AliMUONVTriggerTrackStore& triggerTrackStore)
{
  /// To make the trigger tracks from Local Trigger
  AliDebug(1, "");
  AliCodeTimerAuto("");
  
  AliMUONGlobalTrigger* globalTrigger = triggerStore.Global();
  
  UChar_t gloTrigPat = 0;

  if (globalTrigger)
  {
    gloTrigPat = globalTrigger->GetGlobalResponse();
  }
  
  TIter next(triggerStore.CreateIterator());
  AliMUONLocalTrigger* locTrg(0x0);

  Float_t z11 = AliMUONConstants::DefaultChamberZ(10);
  Float_t z21 = AliMUONConstants::DefaultChamberZ(12);
      
  AliMUONTriggerTrack triggerTrack;
  
  while ( ( locTrg = static_cast<AliMUONLocalTrigger*>(next()) ) )
  {
    Bool_t xTrig=kFALSE;
    Bool_t yTrig=kFALSE;
    
    Int_t localBoardId = locTrg->LoCircuit();
    if ( locTrg->LoSdev()==1 && locTrg->LoDev()==0 && 
         locTrg->LoStripX()==0) xTrig=kFALSE; // no trigger in X
    else xTrig=kTRUE;                         // trigger in X
    if (locTrg->LoTrigY()==1 && 
        locTrg->LoStripY()==15 ) yTrig = kFALSE; // no trigger in Y
    else yTrig = kTRUE;                          // trigger in Y
    
    if (xTrig && yTrig) 
    { // make Trigger Track if trigger in X and Y
      
      Float_t y11 = circuit.GetY11Pos(localBoardId, locTrg->LoStripX()); 
      // need first to convert deviation to [0-30] 
      // (see AliMUONLocalTriggerBoard::LocalTrigger)
      Int_t deviation = locTrg->LoDev(); 
      Int_t sign = 0;
      if ( !locTrg->LoSdev() &&  deviation ) sign=-1;
      if ( !locTrg->LoSdev() && !deviation ) sign= 0;
      if (  locTrg->LoSdev() == 1 )          sign=+1;
      deviation *= sign;
      deviation += 15;
      Int_t stripX21 = locTrg->LoStripX()+deviation+1;
      Float_t y21 = circuit.GetY21Pos(localBoardId, stripX21);       
      Float_t x11 = circuit.GetX11Pos(localBoardId, locTrg->LoStripY());
      
      AliDebug(1, Form(" MakeTriggerTrack %d %d %d %d %f %f %f \n",locTrg->LoCircuit(),
                       locTrg->LoStripX(),locTrg->LoStripX()+locTrg->LoDev()+1,locTrg->LoStripY(),y11, y21, x11));
      
      Float_t thetax = TMath::ATan2( x11 , z11 );
      Float_t thetay = TMath::ATan2( (y21-y11) , (z21-z11) );
      
      triggerTrack.SetX11(x11);
      triggerTrack.SetY11(y11);
      triggerTrack.SetThetax(thetax);
      triggerTrack.SetThetay(thetay);
      triggerTrack.SetGTPattern(gloTrigPat);
      triggerTrack.SetLoTrgNum(localBoardId);
      
      triggerTrackStore.Add(triggerTrack);
    } // board is fired 
  } // end of loop on Local Trigger
}

