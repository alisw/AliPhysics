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

/// \class AliMUONTrackReconstructorK
/////////////////////////////////////
///
/// MUON track reconstructor using the kalman method
///
/// This class contains as data:
/// - the parameters for the track reconstruction
///
/// It contains as methods, among others:
/// - MakeTracks to build the tracks
///
////////////////////////////////////

#include "AliMUONTrackReconstructorK.h"
#include "AliMUONConstants.h"
#include "AliMUONHitForRec.h"
#include "AliMUONObjectPair.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrackK.h" 

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrackReconstructorK) // Class implementation in ROOT context
ClassImp(AliMUONConstants)
/// \endcond

//__________________________________________________________________________
AliMUONTrackReconstructorK::AliMUONTrackReconstructorK(const Option_t* TrackMethod)
  : AliMUONVTrackReconstructor(),
    fTrackMethod(2) //tracking method (2-Kalman 3-Combination-Kalman/Clustering)
{
  /// Constructor for class AliMUONTrackReconstructorK

  if (strstr(TrackMethod,"Kalman")) {
    AliInfo(" *** Tracking with the Kalman filter *** ");
    fTrackMethod = 2;
  } else if (strstr(TrackMethod,"Combi")) {
    AliInfo(" *** Combined cluster / track finder ***");
    fTrackMethod = 3;
  } else AliFatal(Form("Tracking method %s not available",TrackMethod));
  
  // Memory allocation for the TClonesArray of reconstructed tracks
  fRecTracksPtr = new TClonesArray("AliMUONTrackK", 10);
}

//__________________________________________________________________________
AliMUONTrackReconstructorK::~AliMUONTrackReconstructorK(void)
{
  /// Destructor for class AliMUONTrackReconstructorK
  delete fRecTracksPtr;
}

//__________________________________________________________________________
void AliMUONTrackReconstructorK::MakeTracks(void)
{
  /// To make the tracks from the list of segments and points in all stations
  AliDebug(1,"Enter MakeTracks");
  // The order may be important for the following Reset's
  //AZ ResetTracks();
  MakeTrackCandidates();
  if (fRecTracksPtr->GetEntriesFast() == 0) return;
  // Follow tracks in stations(1..) 3, 2 and 1
  FollowTracks();
  // Remove double tracks
  RemoveDoubleTracks();
  // Print out some track info (if necessary)
  EventDump();
  // Fill AliMUONTrack data members
  FillMUONTrack();
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::MakeTrackCandidates(void)
{
  /// To make initial tracks for Kalman filter from segments in stations(1..)  4 and 5
  Int_t istat, iseg;
  TClonesArray *segments;
  AliMUONObjectPair *segment;
  AliMUONTrackK *trackK;

  AliDebug(1,"Enter MakeTrackCandidatesK");

  AliMUONTrackK a(this, fHitsForRecPtr);
  // Loop over stations(1...) 5 and 4
  for (istat=4; istat>=3; istat--) 
  {
    // Make segments in the station
    segments = MakeSegmentsInStation(istat);
    // Loop over segments in the station
    for (iseg=0; iseg<segments->GetEntriesFast(); iseg++) 
    {
      // Transform segments to tracks
      segment = (AliMUONObjectPair*) ((*segments)[iseg]);
      trackK = new ((*fRecTracksPtr)[fNRecTracks++]) AliMUONTrackK(*segment);
    } // for (iseg=0;...)
    delete segments;
  } // for (istat=4;...)
}

//__________________________________________________________________________
void AliMUONTrackReconstructorK::FollowTracks(void)
{
  /// Follow tracks using Kalman filter
  Bool_t ok;
  Int_t icand, ichamBeg = 0, ichamEnd, chamBits;
  AliMUONTrackK *trackK;
  AliMUONHitForRec *hit;
  AliMUONRawCluster *clus;
  TClonesArray *rawclusters;
  clus = 0; rawclusters = 0;

  Double_t simpleBPosition = 0.5 * (AliMUONConstants::CoilZ() + AliMUONConstants::YokeZ());
  Double_t simpleBLength = 0.5 * (AliMUONConstants::CoilL() + AliMUONConstants::YokeL());
  Double_t zDipole1 = simpleBPosition + 0.5 * simpleBLength;
  Double_t zDipole2 = zDipole1 - simpleBLength;

  // Print hits
  trackK = (AliMUONTrackK*) ((*fRecTracksPtr)[0]);

  icand = -1;
  Int_t nSeeds;
  nSeeds = fNRecTracks; // starting number of seeds
  // Loop over track candidates
  while (icand < fNRecTracks-1) {
    icand ++;
    if (trackK->DebugLevel()>0) cout << " *** Kalman track candidate No. " << icand << endl;
    trackK = (AliMUONTrackK*) ((*fRecTracksPtr)[icand]);
    if (trackK->GetRecover() < 0) continue; // failed track

    // Discard candidate which will produce the double track
    /*
    if (icand > 0) {
      ok = CheckCandidate(icand,nSeeds);
      if (!ok) {
        trackK->SetRecover(-1); // mark candidate to be removed
        continue;
      }
    }
    */

    ok = kTRUE;
    if (trackK->GetRecover() == 0) 
      hit = (AliMUONHitForRec*) trackK->GetTrackHits()->Last(); // last hit
    else 
      hit = trackK->GetHitLastOk(); // hit where track stopped

    if (hit) ichamBeg = hit->GetChamberNumber();
    if (ichamBeg == 8 && trackK->GetStation0() == 4) ichamBeg = 10;

    ichamEnd = 0;
    // Check propagation direction
    if (!hit) { ichamBeg = ichamEnd; AliFatal(" ??? "); }
    else if (trackK->GetTrackDir() < 0) {
      ichamEnd = 9; // forward propagation
      ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kFALSE,zDipole1,zDipole2);
      if (ok) {
        ichamBeg = ichamEnd + 1;
        ichamEnd = 6; // backward propagation
	// Change weight matrix and zero fChi2 for backpropagation
        trackK->StartBack();
	trackK->SetTrackDir(1);
        ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kTRUE,zDipole1,zDipole2);
        ichamBeg = ichamEnd;
        ichamEnd = 0;
      }
    } else {
      if (trackK->GetBPFlag()) {
	// backpropagation
        if (ichamBeg == 10) ichamEnd = 8;
        else {
          if (ichamBeg == 9) ichamBeg++;
	  ichamEnd = 6; 
	  // Change weight matrix and zero fChi2 for backpropagation
	  trackK->StartBack();
	}
        ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kTRUE,zDipole1,zDipole2);
        ichamBeg = ichamEnd;
        ichamEnd = 0;
      }
    }

    if (ok) {
      trackK->SetTrackDir(1);
      trackK->SetBPFlag(kFALSE);
      ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kFALSE,zDipole1,zDipole2);
    }
    if (!ok) { trackK->SetRecover(-1); continue; } // mark candidate to be removed

    // Apply smoother
    if (trackK->GetRecover() >= 0) {
      ok = trackK->Smooth();
      if (!ok) trackK->SetRecover(-1); // mark candidate to be removed
    }

    // Majority 3 of 4 in first 2 stations
    if (!ok) continue;
    chamBits = 0;
    Double_t chi2max = 0;
    for (Int_t i=0; i<trackK->GetNTrackHits(); i++) {
      hit = (AliMUONHitForRec*) (*trackK->GetTrackHits())[i];
      chamBits |= BIT(hit->GetChamberNumber());
      if (trackK->GetChi2PerPoint(i) > chi2max) chi2max = trackK->GetChi2PerPoint(i);
    }
    if (!((chamBits&3)==3 || (chamBits>>2&3)==3) && chi2max > 25) {
      //trackK->Recover();
      trackK->SetRecover(-1); //mark candidate to be removed
      continue;
    }
    if (ok) trackK->SetTrackQuality(0); // compute "track quality"
  } // while

  for (Int_t i=0; i<fNRecTracks; i++) {
    trackK = (AliMUONTrackK*) ((*fRecTracksPtr)[i]);
    if (trackK->GetRecover() < 0) fRecTracksPtr->RemoveAt(i);
  }

  // Compress TClonesArray
  fRecTracksPtr->Compress();
  fNRecTracks = fRecTracksPtr->GetEntriesFast();
  return;
}

//__________________________________________________________________________
Bool_t AliMUONTrackReconstructorK::CheckCandidate(Int_t icand, Int_t nSeeds) const
{
  /// Discards track candidate if it will produce the double track (having
  /// the same seed segment hits as hits of a good track found before)
  AliMUONTrackK *track1, *track2;
  AliMUONHitForRec *hit1, *hit2, *hit;
  AliMUONObjectPair *segment1, *segment2;

  track1 = (AliMUONTrackK*) ((*fRecTracksPtr)[icand]);
  hit1 = (AliMUONHitForRec*) (*track1->GetTrackHits())[0]; // 1'st hit
  hit2 = (AliMUONHitForRec*) (*track1->GetTrackHits())[1]; // 2'nd hit

  for (Int_t i=0; i<icand; i++) {
    track2 = (AliMUONTrackK*) ((*fRecTracksPtr)[i]);
    //if (track2->GetRecover() < 0) continue;
    if (track2->GetRecover() < 0 && icand >= nSeeds) continue;
    segment1 = track1->GetStartSegment();
    segment2 = track2->GetStartSegment();
    if (segment1->First()  == segment2->First() &&
        segment1->Second() == segment2->Second()) {
      return kFALSE;
    } else {
      Int_t nSame = 0;
      for (Int_t j=0; j<track2->GetNTrackHits(); j++) {
        hit = (AliMUONHitForRec*) (*track2->GetTrackHits())[j];
        if (hit == hit1 || hit == hit2) {
          nSame++;
          if (nSame == 2) return kFALSE;
        }
      } // for (Int_t j=0;
    }
  } // for (Int_t i=0;
  return kTRUE;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::RemoveDoubleTracks(void)
{
  /// Removes double tracks (sharing more than half of their hits).
  /// Keeps the track with higher quality
  AliMUONTrackK *track1, *track2, *trackToKill;

  // Sort tracks according to their quality
  fRecTracksPtr->Sort();

  // Loop over first track of the pair
  track1 = (AliMUONTrackK*) fRecTracksPtr->First();
  Int_t debug = track1->DebugLevel();
  while (track1) {
    // Loop over second track of the pair
    track2 = (AliMUONTrackK*) fRecTracksPtr->After(track1);
    while (track2) {
      // Check whether or not to keep track2
      if (!track2->KeepTrack(track1)) {
        if (debug >= 0) cout << " Killed track: " << 1/(*track2->GetTrackParameters())(4,0) <<
	  " " << track2->GetTrackQuality() << endl;
        trackToKill = track2;
        track2 = (AliMUONTrackK*) fRecTracksPtr->After(track2);
        trackToKill->Kill();
        fRecTracksPtr->Compress();
      } else track2 = (AliMUONTrackK*) fRecTracksPtr->After(track2);
    } // track2
    track1 = (AliMUONTrackK*) fRecTracksPtr->After(track1);
  } // track1

  fNRecTracks = fRecTracksPtr->GetEntriesFast();
  if (debug >= 0) cout << " Number of Kalman tracks: " << fNRecTracks << endl;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::FillMUONTrack()
{
  // Set track parameters at hits for Kalman track. Fill fTrackParamAtHit of AliMUONTrack's
  AliMUONTrackK *track;
  track = (AliMUONTrackK*) fRecTracksPtr->First();
  while (track) {
    track->FillMUONTrack();
    track = (AliMUONTrackK*) fRecTracksPtr->After(track);
  } 
  return;
}


  //__________________________________________________________________________
void AliMUONTrackReconstructorK::EventDump(void)
{
  /// Dump reconstructed event (track parameters at vertex and at first hit),
  /// and the particle parameters
  AliDebug(1,"****** enter EventDump ******");
  AliDebug(1, Form("Number of Reconstructed tracks : %d", fNRecTracks));
                                                                                
  Int_t debug =   ((AliMUONTrackK*) fRecTracksPtr->First())->DebugLevel();
  if (debug < 0) return;
                                                                                
  for (Int_t i = 0; i < fNRecTracks; ++i) {
    AliMUONTrackK *track = (AliMUONTrackK*) fRecTracksPtr->UncheckedAt(i);
    track->Print("");
  }
}
