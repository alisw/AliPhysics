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

////////////////////////////////////
///
/// \class AliMUONVTrackReconstructor
/// Virtual MUON track reconstructor in ALICE (class renamed from AliMUONEventReconstructor)
///
/// This class contains as data:
/// * a pointer to the array of hits to be reconstructed (the event)
/// * a pointer to the array of segments made with these hits inside each station
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
///
////////////////////////////////////

#include "AliMUONVTrackReconstructor.h"

#include "AliMUONConstants.h"
#include "AliMUONRawCluster.h"
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
#include "AliMUONRawCluster.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMpDEManager.h"

#include "AliLog.h"
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

const Double_t AliMUONVTrackReconstructor::fgkSigmaToCutForTracking = 6.0;
const Double_t AliMUONVTrackReconstructor::fgkSigmaToCutForImprovement = 4.0;
const Bool_t   AliMUONVTrackReconstructor::fgkTrackAllTracks = kTRUE;
const Double_t AliMUONVTrackReconstructor::fgkMaxTrackingDistanceBending    = 2.;
const Double_t AliMUONVTrackReconstructor::fgkMaxTrackingDistanceNonBending = 2.;
const Bool_t   AliMUONVTrackReconstructor::fgkRecoverTracks = kTRUE;
const Bool_t   AliMUONVTrackReconstructor::fgkImproveTracks = kTRUE;


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
  
  TIter next(clusterStore.CreateIterator());
  AliMUONRawCluster* clus(0x0);
  Int_t iclus(0);
  
  while ( ( clus = static_cast<AliMUONRawCluster*>(next()) ) )
  {
    // new AliMUONHitForRec from raw cluster
    // and increment number of AliMUONHitForRec's (total and in chamber)
    AliMUONHitForRec* hitForRec = new ((*fHitsForRecPtr)[fNHitsForRec]) AliMUONHitForRec(clus);
    fNHitsForRec++;
    // more information into HitForRec
    hitForRec->SetBendingReso2(clus->GetErrY() * clus->GetErrY());
    hitForRec->SetNonBendingReso2(clus->GetErrX() * clus->GetErrX());
    //  original raw cluster
    Int_t ch = AliMpDEManager::GetChamberId(clus->GetDetElemId());
    hitForRec->SetChamberNumber(ch);
    hitForRec->SetHitNumber(iclus);
    // Z coordinate of the raw cluster (cm)
    hitForRec->SetZ(clus->GetZ(0));
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
  
  if (!trackParam.CovariancesExist()) AliWarning(" track parameter covariance matrix does not exist");
  
  // extrapolate track parameters and covariances at the z position of the tested hit
  if (&trackParam != &trackParamAtHit) {
    trackParamAtHit = trackParam;
    AliMUONTrackExtrap::ExtrapToZCov(&trackParamAtHit, hitForRec->GetZ(), updatePropagator);
  }
  
  // set pointer to hit into trackParamAtHit
  trackParamAtHit.SetHitForRecPtr(hitForRec);
  
  // Set differences between trackParam and hitForRec in the bending and non bending directions
  TMatrixD dPos(2,1);
  dPos(0,0) = hitForRec->GetNonBendingCoor() - trackParamAtHit.GetNonBendingCoor();
  dPos(1,0) = hitForRec->GetBendingCoor() - trackParamAtHit.GetBendingCoor();
  
  // quick test of hitForRec compatibility within a wide road of x*y = 10*1 cm2 to save computing time
  if (TMath::Abs(dPos(0,0)) > fgkMaxTrackingDistanceNonBending ||
      TMath::Abs(dPos(1,0)) > fgkMaxTrackingDistanceBending) return 1.e10;
  
  // Set the error matrix from trackParam covariances and hitForRec resolution
  const TMatrixD& kParamCov = trackParamAtHit.GetCovariances();
  TMatrixD error(2,2);
  error(0,0) = kParamCov(0,0) + hitForRec->GetNonBendingReso2();
  error(0,1) = kParamCov(0,2);
  error(1,0) = kParamCov(2,0);
  error(1,1) = kParamCov(2,2) + hitForRec->GetBendingReso2();
  
  // Invert the error matrix for Chi2 calculation
  if (error.Determinant() != 0) {
    error.Invert();
  } else {
    AliWarning(" Determinant error=0");
    return 1.e10;
  }
  
  // Compute the Chi2 value
  TMatrixD tmp(dPos,TMatrixD::kTransposeMult,error);
  TMatrixD result(tmp,TMatrixD::kMult,dPos);
  
  return result(0,0);
  
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::ValidateTracksWithTrigger(AliMUONVTrackStore& trackStore,
                                                           const AliMUONVTriggerTrackStore& triggerTrackStore,
                                                           const AliMUONVTriggerStore& triggerStore,
                                                           const AliMUONTrackHitPattern& trackHitPattern)
{
  /// Try to match track from tracking system with trigger track
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

