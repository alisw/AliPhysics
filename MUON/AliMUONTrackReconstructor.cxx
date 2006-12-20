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

/// \class AliMUONTrackReconstructor
/// MUON track reconstructor using the original method
///
/// This class contains as data:
/// - the parameters for the track reconstruction
///
/// It contains as methods, among others:
/// - MakeTracks to build the tracks
///

#include <stdlib.h>
#include <Riostream.h>
#include <TMatrixD.h>

#include "AliMUONVTrackReconstructor.h"
#include "AliMUONTrackReconstructor.h"
#include "AliMUONData.h"
#include "AliMUONConstants.h"
#include "AliMUONRawCluster.h"
#include "AliMUONHitForRec.h"
#include "AliMUONObjectPair.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliLog.h"
#include <TMinuit.h>

// Functions to be minimized with Minuit
void TrackChi2(Int_t &NParam, Double_t *Gradient, Double_t &Chi2, Double_t *Param, Int_t Flag);
void TrackChi2MCS(Int_t &NParam, Double_t *Gradient, Double_t &Chi2, Double_t *Param, Int_t Flag);

Double_t MultipleScatteringAngle2(AliMUONTrackParam *param);

/// \cond CLASSIMP
ClassImp(AliMUONTrackReconstructor) // Class implementation in ROOT context
/// \endcond

//************* Defaults parameters for reconstruction
const Double_t AliMUONTrackReconstructor::fgkMaxNormChi2 = 100.0;
const Bool_t AliMUONTrackReconstructor::fgkTrackAllTracks = kFALSE;

//__________________________________________________________________________
AliMUONTrackReconstructor::AliMUONTrackReconstructor(AliMUONData* data)
  : AliMUONVTrackReconstructor(data)
{
  /// Constructor for class AliMUONTrackReconstructor
  
  // Memory allocation for the TClonesArray of reconstructed tracks
  fRecTracksPtr = new TClonesArray("AliMUONTrack", 10);
}

  //__________________________________________________________________________
AliMUONTrackReconstructor::~AliMUONTrackReconstructor(void)
{
  /// Destructor for class AliMUONTrackReconstructor
  delete fRecTracksPtr;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::AddHitsForRecFromRawClusters()
{
  /// To add to the list of hits for reconstruction all the raw clusters
  TTree *treeR;
  AliMUONHitForRec *hitForRec;
  AliMUONRawCluster *clus;
  Int_t iclus, nclus;
  TClonesArray *rawclusters;
  AliDebug(1,"Enter AddHitsForRecFromRawClusters");
  
  treeR = fMUONData->TreeR();
  if (!treeR) {
    AliError("TreeR must be loaded");
    exit(0);
  }
  
  fMUONData->SetTreeAddress("RC");
  fMUONData->GetRawClusters(); // only one entry  
  
  // Loop over tracking chambers
  for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    rawclusters =fMUONData->RawClusters(ch);
    nclus = (Int_t) (rawclusters->GetEntries());
    // Loop over (cathode correlated) raw clusters
    for (iclus = 0; iclus < nclus; iclus++) {
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(iclus);
      // new AliMUONHitForRec from raw cluster
      // and increment number of AliMUONHitForRec's (total and in chamber)
      hitForRec = new ((*fHitsForRecPtr)[fNHitsForRec]) AliMUONHitForRec(clus);
      fNHitsForRec++;
      // more information into HitForRec
      hitForRec->SetBendingReso2(clus->GetErrY() * clus->GetErrY());
      hitForRec->SetNonBendingReso2(clus->GetErrX() * clus->GetErrX());
      //  original raw cluster
      hitForRec->SetChamberNumber(ch);
      hitForRec->SetHitNumber(iclus);
      // Z coordinate of the raw cluster (cm)
      hitForRec->SetZ(clus->GetZ(0));
      StdoutToAliDebug(3,
                       cout << "Chamber " << ch <<
                       " raw cluster  " << iclus << " : " << endl;
                       clus->Print("full");
                       cout << "AliMUONHitForRec number (1...): " << fNHitsForRec << endl;
                       hitForRec->Print("full");
                       );
    } // end of cluster loop
  } // end of chamber loop
  SortHitsForRecWithIncreasingChamber(); 
  
  AliDebug(1,"End of AddHitsForRecFromRawClusters");
    if (AliLog::GetGlobalDebugLevel() > 0) {
      AliDebug(1, Form("NHitsForRec: %d",fNHitsForRec));
      for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
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
  
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::MakeTracks(void)
{
  /// To make the tracks from the list of segments and points in all stations
  AliDebug(1,"Enter MakeTracks");
  // Look for candidates from at least 3 aligned points in stations(1..) 4 and 5
  MakeTrackCandidates();
  // Follow tracks in stations(1..) 3, 2 and 1
  FollowTracks();
  // Remove double tracks
  RemoveDoubleTracks();
  // Propagate tracks to the vertex through absorber
  ExtrapTracksToVertex();
  // Fill out the AliMUONTrack's
  FillMUONTrack();
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::MakeTrackCandidates(void)
{
  /// To make track candidates:
  /// Start with segments station(1..) 4 or 5 then follow track in station 5 or 4.
  /// Good candidates are made of at least three hitForRec's.
  /// Keep only best candidates or all of them according to the flag fgkTrackAllTracks.
  TClonesArray *segments;
  AliMUONObjectPair *segment;
  AliMUONHitForRec *hitForRec1, *hitForRec2;
  AliMUONTrack *track;
  AliMUONTrackParam *trackParamAtFirstHit;
  Int_t iCandidate = 0;

  AliDebug(1,"Enter MakeTrackCandidates");

  // Loop over stations(1..) 5 and 4 and make track candidates
  for (Int_t istat=4; istat>=3; istat--) {
    // Make segments in the station
    segments = MakeSegmentsInStation(istat);
    // Loop over segments
    for (Int_t iseg=0; iseg<segments->GetEntriesFast(); iseg++) {
      AliDebug(1,Form("Making primary candidate(1..) %d",++iCandidate));
      // Transform segments to tracks and put them at the end of fRecTracksPtr
      segment = (AliMUONObjectPair*) ((*segments)[iseg]);
      hitForRec1 = (AliMUONHitForRec*) segment->First();
      hitForRec2 = (AliMUONHitForRec*) segment->Second();
      track = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(hitForRec1, hitForRec2);
      fNRecTracks++;
      // Add MCS effects in parameter covariances
      trackParamAtFirstHit = (AliMUONTrackParam*) (track->GetTrackParamAtHit()->First());
      AliMUONTrackExtrap::AddMCSEffectInTrackParamCov(trackParamAtFirstHit,AliMUONConstants::ChamberThicknessInX0(),1.);
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
        cout<<endl<<"Track parameter covariances at first hit with multiple Coulomb scattering effects:"<<endl;
        trackParamAtFirstHit->GetCovariances()->Print();
      }
      // Look for compatible hitForRec(s) in the other station
      FollowTrackInStation(track,7-istat);
    }
    // delete the array of segments
    delete segments;
  }
  fRecTracksPtr->Compress(); // this is essential before checking tracks
  
  // Keep all different tracks or only the best ones as required
  if (fgkTrackAllTracks) RemoveIdenticalTracks();
  else RemoveDoubleTracks();
  
  AliDebug(1,Form("Number of good candidates = %d",fNRecTracks));
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::RemoveIdenticalTracks(void)
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
void AliMUONTrackReconstructor::RemoveDoubleTracks(void)
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
void AliMUONTrackReconstructor::FollowTracks(void)
{
  /// Follow tracks in stations(1..) 3, 2 and 1
  AliDebug(1,"Enter FollowTracks");
  
  AliMUONTrack *track, *nextTrack;
  AliMUONTrackParam *trackParamAtFirstHit;
  Double_t numberOfDegFree, chi2Norm;
  Int_t CurrentNRecTracks;
  
  for (Int_t station = 2; station >= 0; station--) {
    // Save the actual number of reconstructed track in case of
    // tracks are added or suppressed during the tracking procedure
    // !! Do not compress fRecTracksPtr until the end of the loop over tracks !!
    CurrentNRecTracks = fNRecTracks;
    for (Int_t iRecTrack = 0; iRecTrack <CurrentNRecTracks; iRecTrack++) {
      AliDebug(1,Form("FollowTracks: track candidate(1..) %d", iRecTrack+1));
      track = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iRecTrack);
      // Fit the track:
      // Do not take into account the multiple scattering to speed up the fit
      // Calculate the track parameter covariance matrix
      // If "station" is station(1..) 3 then use the vertex to better constrain the fit
      if (station==2) {
        SetVertexForFit(track);
        track->SetFitWithVertex(kTRUE);
      } else track->SetFitWithVertex(kFALSE);
      Fit(track,kFALSE, kTRUE);
      // Remove the track if the normalized chi2 is too high
      numberOfDegFree = (2. * track->GetNTrackHits() - 5.);
      if (numberOfDegFree > 0) chi2Norm = track->GetFitFMin() / numberOfDegFree;
      else chi2Norm = 1.e10;
      if (chi2Norm > fgkMaxNormChi2) {
  	fRecTracksPtr->Remove(track);
  	fNRecTracks--;
  	continue;
      }
      // Add MCS effects in parameter covariances
      trackParamAtFirstHit = (AliMUONTrackParam*) (track->GetTrackParamAtHit()->First());
      AliMUONTrackExtrap::AddMCSEffectInTrackParamCov(trackParamAtFirstHit,AliMUONConstants::ChamberThicknessInX0(),1.);
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
        cout<<endl<<"Track parameter covariances at first hit with multiple Coulomb scattering effects:"<<endl;
        trackParamAtFirstHit->GetCovariances()->Print();
      }
      // Look for compatible hitForRec in station(0..) "station"
      FollowTrackInStation(track,station);
    }
    // Compress fRecTracksPtr for the next step
    fRecTracksPtr->Compress();
    // Keep only the best tracks if required
    if (!fgkTrackAllTracks) RemoveDoubleTracks();
  }
  
  // Last fit of track candidates with all station
  // Take into account the multiple scattering and remove bad tracks
  Int_t trackIndex = -1;
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    trackIndex++;
    nextTrack = (AliMUONTrack*) fRecTracksPtr->After(track); // prepare next track
    track->SetFitWithVertex(kFALSE); // just to be sure
    Fit(track,kTRUE, kFALSE);
    // Printout for debuging
    if (AliLog::GetGlobalDebugLevel() >= 3) {
      cout << "FollowTracks: track candidate(0..) " << trackIndex << " after final fit" << endl;
      track->RecursiveDump();
    } 
    // Remove the track if the normalized chi2 is too high
    numberOfDegFree = (2.0 * track->GetNTrackHits() - 5);
    if (numberOfDegFree > 0) chi2Norm = track->GetFitFMin() / numberOfDegFree;
    else chi2Norm = 1.e10;
    if (chi2Norm > fgkMaxNormChi2) {
      fRecTracksPtr->Remove(track);
      fNRecTracks--;
    }
    track = nextTrack;
  }
  fRecTracksPtr->Compress();
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::FollowTrackInStation(AliMUONTrack* trackCandidate, Int_t nextStation)
{
  /// Follow trackCandidate in station(0..) nextStation and search for compatible HitForRec(s)
  /// Keep all possibilities or only the best one(s) according to the flag fgkTrackAllTracks:
  /// kTRUE:  duplicate "trackCandidate" if there are several possibilities and add the new tracks at the end of
  ///         fRecTracksPtr to avoid conficts with other track candidates at this current stage of the tracking procedure.
  ///         Remove the obsolete "trackCandidate" at the end.
  /// kFALSE: add only the best hit(s) to the "trackCandidate". Try to add a couple of hits in priority.
  AliDebug(1,Form("Enter FollowTrackInStation(1..) %d", nextStation+1));
  
  Int_t ch1 = 2*nextStation;
  Int_t ch2 = 2*nextStation+1;
  Double_t zCh2 = AliMUONConstants::DefaultChamberZ(ch2);
  Double_t chi2WithOneHitForRec = 1.e10;
  Double_t chi2WithTwoHitForRec = 1.e10;
  Double_t maxChi2WithOneHitForRec = 2.*fgkMaxNormChi2; // 2 because 2 quantities in chi2
  Double_t maxChi2WithTwoHitForRec = 4.*fgkMaxNormChi2; // 4 because 4 quantities in chi2
  Double_t bestChi2WithOneHitForRec = maxChi2WithOneHitForRec;
  Double_t bestChi2WithTwoHitForRec = maxChi2WithTwoHitForRec;
  AliMUONTrack *newTrack = 0x0;
  AliMUONHitForRec *hitForRecCh1, *hitForRecCh2;
  AliMUONHitForRec *bestHitForRec1 = 0x0, *bestHitForRec2 = 0x0;
  Bool_t *hitForRecCh1Used = new Bool_t[fNHitsForRecPerChamber[ch1]];
  for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) hitForRecCh1Used[hit1] = kFALSE;
  //
  //Extrapolate trackCandidate to chamber "ch2" to save computing time in the next steps
  AliMUONTrackParam *extrapTrackParamPtr = trackCandidate->GetExtrapTrackParam();
  *extrapTrackParamPtr = *((AliMUONTrackParam*)(trackCandidate->GetTrackParamAtHit()->First()));
  AliMUONTrackExtrap::ExtrapToZCov(extrapTrackParamPtr, zCh2);
  AliMUONTrackParam extrapTrackParamSave(*extrapTrackParamPtr);
  //
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    TMatrixD* paramCovForDebug = extrapTrackParamPtr->GetCovariances();
    cout<<endl<<"Track parameter covariances at first hit extrapolated to z = "<<zCh2<<":"<<endl;
    paramCovForDebug->Print();
  }
  //
  // look for candidates in chamber 2 
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "FollowTrackInStation: look for hits in chamber(1..): " << ch2+1 << endl;
  }
  for (Int_t hit2 = 0; hit2 < fNHitsForRecPerChamber[ch2]; hit2++) {
    hitForRecCh2 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch2]+hit2);
    // extrapolate track parameters and covariances only once for this hit
    AliMUONTrackExtrap::ExtrapToZCov(extrapTrackParamPtr, hitForRecCh2->GetZ());
    chi2WithOneHitForRec = trackCandidate->TryOneHitForRec(hitForRecCh2);
    // if good chi2 then try to attach a hitForRec in the other chamber too
    if (chi2WithOneHitForRec < maxChi2WithOneHitForRec) {
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: look for second hits in chamber(1..): " << ch1+1 << endl;
      }
      Bool_t foundSecondHit = kFALSE;
      for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) {
        hitForRecCh1 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch1]+hit1);
	chi2WithTwoHitForRec = trackCandidate->TryTwoHitForRec(hitForRecCh2, hitForRecCh1); // order hits like that to save computing time
        // if good chi2 then create a new track by adding the 2 hitForRec to the "trackCandidate"
	if (chi2WithTwoHitForRec < maxChi2WithTwoHitForRec) {
	  foundSecondHit = kTRUE;
          if (fgkTrackAllTracks) {
	    // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
            newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(*trackCandidate);
	    fNRecTracks++;
            AliMUONTrackParam TrackParam1(extrapTrackParamSave);
            AliMUONTrackExtrap::ExtrapToZ(&TrackParam1, hitForRecCh1->GetZ());
	    newTrack->AddTrackParamAtHit(&TrackParam1,hitForRecCh1);
            AliMUONTrackParam TrackParam2(extrapTrackParamSave);
            AliMUONTrackExtrap::ExtrapToZ(&TrackParam2, hitForRecCh2->GetZ());
	    newTrack->AddTrackParamAtHit(&TrackParam2,hitForRecCh2);
            // Sort TrackParamAtHit according to increasing Z
            newTrack->GetTrackParamAtHit()->Sort();
	    // Update the chi2 of the new track
	    if (newTrack->GetFitFMin()<0) newTrack->SetFitFMin(chi2WithTwoHitForRec);
	    else newTrack->SetFitFMin(newTrack->GetFitFMin() + chi2WithTwoHitForRec);
	    // Tag hitForRecCh1 as used
	    hitForRecCh1Used[hit1] = kTRUE;
	    // Printout for debuging
	    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	      cout << "FollowTrackInStation: added two hits in station(1..): " << nextStation+1
	           << " (Chi2 = " << chi2WithTwoHitForRec << ")" << endl;
	      if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	    }
          } else if (chi2WithTwoHitForRec < bestChi2WithTwoHitForRec) {
	    // keep track of the best couple of hits
	    bestChi2WithTwoHitForRec = chi2WithTwoHitForRec;
	    bestHitForRec1 = hitForRecCh1;
	    bestHitForRec2 = hitForRecCh2;
          }
	}
      }
      // if no hitForRecCh1 found then consider to add hitForRecCh2 only
      if (!foundSecondHit) {
        if (fgkTrackAllTracks) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
          newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(*trackCandidate);
	  fNRecTracks++;
          AliMUONTrackParam TrackParam1(extrapTrackParamSave);
          AliMUONTrackExtrap::ExtrapToZ(&TrackParam1, hitForRecCh2->GetZ());
          newTrack->AddTrackParamAtHit(&TrackParam1,hitForRecCh2);
          // Sort TrackParamAtHit according to increasing Z
          newTrack->GetTrackParamAtHit()->Sort();
	  // Update the chi2 of the new track
	  if (newTrack->GetFitFMin()<0) newTrack->SetFitFMin(chi2WithOneHitForRec);
	  else newTrack->SetFitFMin(newTrack->GetFitFMin() + chi2WithOneHitForRec);
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowTrackInStation: added one hit in chamber(1..): " << ch2+1
	  	 << " (Chi2 = " << chi2WithOneHitForRec << ")" << endl;
	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	  }
	} if (!bestHitForRec2 && chi2WithOneHitForRec < bestChi2WithOneHitForRec) {
	  // keep track of the best single hitForRec except if a couple
          // of hits has already been found (i.e. bestHitForRec2!=0x0)
	  bestChi2WithOneHitForRec = chi2WithOneHitForRec;
	  bestHitForRec1 = hitForRecCh2;
        }
      }
    }
    // reset the extrapolated track parameter for next step
    trackCandidate->SetExtrapTrackParam(&extrapTrackParamSave);
  }
  //
  // look for candidates in chamber 1 not already attached to a track
  // if we want to keep all possible tracks or if no good couple of hitForRec has been found
  if (fgkTrackAllTracks || !bestHitForRec2) {
    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
      cout << "FollowTrackInStation: look for single hits in chamber(1..): " << ch1+1 << endl;
    }
    for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) {
      hitForRecCh1 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch1]+hit1);
      if (hitForRecCh1Used[hit1]) continue; // Skip hitForRec already used
      chi2WithOneHitForRec = trackCandidate->TryOneHitForRec(hitForRecCh1);
      // if good chi2 then create a new track by adding the good hitForRec in "ch1" to the "trackCandidate"
      // We do not try to attach a hitForRec in the other chamber too since it has already been done above
      if (chi2WithOneHitForRec < maxChi2WithOneHitForRec) {
  	if (fgkTrackAllTracks) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
  	  newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(*trackCandidate);
	  fNRecTracks++;
  	  AliMUONTrackParam TrackParam1(extrapTrackParamSave);
  	  AliMUONTrackExtrap::ExtrapToZ(&TrackParam1, hitForRecCh1->GetZ());
  	  newTrack->AddTrackParamAtHit(&TrackParam1,hitForRecCh1);
  	  // Sort TrackParamAtHit according to increasing Z
  	  newTrack->GetTrackParamAtHit()->Sort();
	  // Update the chi2 of the new track
	  if (newTrack->GetFitFMin()<0) newTrack->SetFitFMin(chi2WithOneHitForRec);
	  else newTrack->SetFitFMin(newTrack->GetFitFMin() + chi2WithOneHitForRec);
  	  // Printout for debuging
  	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	    cout << "FollowTrackInStation: added one hit in chamber(1..): " << ch1+1
  		 << " (Chi2 = " << chi2WithOneHitForRec << ")" << endl;
  	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
  	  }
  	} if (!bestHitForRec2 && chi2WithOneHitForRec < bestChi2WithOneHitForRec) {
  	  // keep track of the best single hitForRec except if a couple
  	  // of hits has already been found (i.e. bestHitForRec1!=0x0)
  	  bestChi2WithOneHitForRec = chi2WithOneHitForRec;
  	  bestHitForRec1 = hitForRecCh1;
  	}
      }
    }
  }
  //
  // fill out the best track if required else clean up the fRecTracksPtr array
  if (!fgkTrackAllTracks && bestHitForRec1) {
    AliMUONTrackParam TrackParam1(extrapTrackParamSave);
    AliMUONTrackExtrap::ExtrapToZ(&TrackParam1, bestHitForRec1->GetZ());
    trackCandidate->AddTrackParamAtHit(&TrackParam1,bestHitForRec1);
    if (bestHitForRec2) {
      AliMUONTrackParam TrackParam2(extrapTrackParamSave);
      AliMUONTrackExtrap::ExtrapToZ(&TrackParam2, bestHitForRec2->GetZ());
      trackCandidate->AddTrackParamAtHit(&TrackParam2,bestHitForRec2);
      // Sort TrackParamAtHit according to increasing Z
      trackCandidate->GetTrackParamAtHit()->Sort();
      // Update the chi2 of the new track
      if (trackCandidate->GetFitFMin()<0) trackCandidate->SetFitFMin(bestChi2WithTwoHitForRec);
      else trackCandidate->SetFitFMin(trackCandidate->GetFitFMin() + bestChi2WithTwoHitForRec);
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the two best hits in station(1..): " << nextStation+1
             << " (Chi2 = " << bestChi2WithTwoHitForRec << ")" << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
    } else {
      // Sort TrackParamAtHit according to increasing Z
      trackCandidate->GetTrackParamAtHit()->Sort();
      // Update the chi2 of the new track
      if (trackCandidate->GetFitFMin()<0) trackCandidate->SetFitFMin(bestChi2WithOneHitForRec);
      else trackCandidate->SetFitFMin(trackCandidate->GetFitFMin() + bestChi2WithOneHitForRec);
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the best hit in station(1..): " << nextStation+1
             << " (Chi2 = " << bestChi2WithOneHitForRec << ")" << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
    }
  } else {
    fRecTracksPtr->Remove(trackCandidate); // obsolete track
    fNRecTracks--;
  }
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::SetVertexForFit(AliMUONTrack* trackCandidate)
{
  /// Add the vertex as a measured hit to constrain the fit of the "trackCandidate"
  /// Compute the vertex resolution from natural vertex dispersion and
  /// multiple scattering effets according to trackCandidate path in absorber
  /// It is necessary to account for multiple scattering effects here instead of during the fit of
  /// the "trackCandidate" to do not influence the result by changing track resolution at vertex
  AliDebug(1,"Enter SetVertexForFit");
  
  Double_t nonBendingReso2 = fNonBendingVertexDispersion * fNonBendingVertexDispersion;
  Double_t bendingReso2 = fBendingVertexDispersion * fBendingVertexDispersion;
  // add multiple scattering effets
  AliMUONTrackParam paramAtVertex(*((AliMUONTrackParam*)(trackCandidate->GetTrackParamAtHit()->First())));
  paramAtVertex.DeleteCovariances(); // to be sure to account only for multiple scattering
  AliMUONTrackExtrap::ExtrapToVertexUncorrected(&paramAtVertex,0.);
  TMatrixD* paramCov = paramAtVertex.GetCovariances();
  nonBendingReso2 += (*paramCov)(0,0);
  bendingReso2 += (*paramCov)(2,2);
  // Set the vertex
  AliMUONHitForRec vertex; // Coordinates set to (0.,0.,0.) by default
  vertex.SetNonBendingReso2(nonBendingReso2);
  vertex.SetBendingReso2(bendingReso2);
  trackCandidate->SetVertex(&vertex);
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::Fit(AliMUONTrack *track, Bool_t includeMCS, Bool_t calcCov)
{
  /// Fit the track "track" with or without multiple Coulomb scattering according to "includeMCS".
  
  Double_t benC, errorParam, invBenP, nonBenC, x, y;
  AliMUONTrackParam *trackParam;
  Double_t arg[1], fedm, errdef, fitFMin;
  Int_t npari, nparx;
  Int_t status, covStatus;
  
  // Clear MINUIT parameters
  gMinuit->mncler();
  // Give the fitted track to MINUIT
  gMinuit->SetObjectFit(track);
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    // Define print level
    arg[0] = 1;
    gMinuit->mnexcm("SET PRI", arg, 1, status);
    // Print covariance matrix
    gMinuit->mnexcm("SHO COV", arg, 0, status);
  } else {
    arg[0] = -1;
    gMinuit->mnexcm("SET PRI", arg, 1, status);
  }
  // No warnings
  gMinuit->mnexcm("SET NOW", arg, 0, status);
  // Define strategy
  //arg[0] = 2;
  //gMinuit->mnexcm("SET STR", arg, 1, status);
  
  // Switch between available FCN according to "includeMCS"
  if (includeMCS) gMinuit->SetFCN(TrackChi2MCS);
  else gMinuit->SetFCN(TrackChi2);
  
  // Set fitted parameters (!! The order is very important for the covariance matrix !!)
  trackParam = (AliMUONTrackParam*) (track->GetTrackParamAtHit()->First());
  // could be tried with no limits for the search (min=max=0) ????
  // mandatory limits in non Bending to avoid NaN values of parameters
  gMinuit->mnparm(0, "X", trackParam->GetNonBendingCoor(), 0.03, -500.0, 500.0, status);
  gMinuit->mnparm(1, "NonBenS", trackParam->GetNonBendingSlope(), 0.001, -0.5, 0.5, status);
  // mandatory limits in Bending to avoid NaN values of parameters
  gMinuit->mnparm(2, "Y", trackParam->GetBendingCoor(), 0.10, -500.0, 500.0, status);
  gMinuit->mnparm(3, "BenS", trackParam->GetBendingSlope(), 0.001, -0.5, 0.5, status);
  gMinuit->mnparm(4, "InvBenP", trackParam->GetInverseBendingMomentum(), 0.003, -0.4, 0.4, status);
  
  // minimization
  gMinuit->mnexcm("MIGRAD", arg, 0, status);
  
  // Calculate the covariance matrix more accurately if required
  if (calcCov) gMinuit->mnexcm("HESSE", arg, 0, status);
  
  // get results into "invBenP", "benC", "nonBenC" ("x", "y")
  gMinuit->GetParameter(0, x, errorParam);
  trackParam->SetNonBendingCoor(x);
  gMinuit->GetParameter(1, nonBenC, errorParam);
  trackParam->SetNonBendingSlope(nonBenC);
  gMinuit->GetParameter(2, y, errorParam);
  trackParam->SetBendingCoor(y);
  gMinuit->GetParameter(3, benC, errorParam);
  trackParam->SetBendingSlope(benC);
  gMinuit->GetParameter(4, invBenP, errorParam);
  trackParam->SetInverseBendingMomentum(invBenP);
  
  // global result of the fit
  gMinuit->mnstat(fitFMin, fedm, errdef, npari, nparx, covStatus);
  track->SetFitFMin(fitFMin);
  
  // Get the covariance matrix if required
  if (calcCov) {
    // Covariance matrix according to HESSE status
    // If problem then keep only the diagonal terms (variances)
    Double_t matrix[5][5];
    gMinuit->mnemat(&matrix[0][0],5);
    if (covStatus == 3) trackParam->SetCovariances(matrix);
    else trackParam->SetVariances(matrix);
  } else *(trackParam->GetCovariances()) = 0;
  
}

  //__________________________________________________________________________
void TrackChi2(Int_t & /*NParam*/, Double_t * /*Gradient*/, Double_t &Chi2, Double_t *Param, Int_t /*Flag*/)
{
  /// Return the "Chi2" to be minimized with Minuit for track fitting,
  /// with "NParam" parameters
  /// and their current values in array pointed to by "Param".
  /// Assumes that the track hits are sorted according to increasing Z.
  /// Track parameters at each TrackHit are updated accordingly.
  /// Multiple Coulomb scattering is not taken into account
  
  AliMUONTrack *trackBeingFitted = (AliMUONTrack*) gMinuit->GetObjectFit();
//  AliMUONTrack *trackBeingFitted = (AliMUONTrack*) AliMUONTrackReconstructor::Fitter()->GetObjectFit();
  AliMUONTrackParam param1;
  AliMUONTrackParam* trackParamAtHit;
  AliMUONHitForRec* hitForRec;
  Chi2 = 0.0; // initialize Chi2
  Double_t dX, dY;
  
  // copy of track parameters to be fitted
  param1 = *((AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->First()));
  param1.SetNonBendingCoor(Param[0]);
  param1.SetNonBendingSlope(Param[1]);
  param1.SetBendingCoor(Param[2]);
  param1.SetBendingSlope(Param[3]);
  param1.SetInverseBendingMomentum(Param[4]);
  
  // Take the vertex into account in the fit if required
  if (trackBeingFitted->GetFitWithVertex()) {
    AliMUONTrackParam paramAtVertex(param1);
    AliMUONTrackExtrap::ExtrapToZ(&paramAtVertex, 0.);
    AliMUONHitForRec *vertex = trackBeingFitted->GetVertex();
    if (!vertex) {
      cout<<"Error in TrackChi2: Want to use the vertex in tracking but it has not been created!!"<<endl;
      exit(-1);
    }
    dX = vertex->GetNonBendingCoor() - paramAtVertex.GetNonBendingCoor();
    dY = vertex->GetBendingCoor() - paramAtVertex.GetBendingCoor();
    Chi2 += dX * dX / vertex->GetNonBendingReso2() + dY * dY / vertex->GetBendingReso2();
  }
  
  // Follow track through all planes of track hits
  trackParamAtHit = (AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->First());
  while (trackParamAtHit) {
    hitForRec = trackParamAtHit->GetHitForRecPtr();
    // extrapolation to the plane of the hitForRec attached to the current trackParamAtHit
    AliMUONTrackExtrap::ExtrapToZ(&param1, hitForRec->GetZ());
    // update track parameters of the current hit
    trackParamAtHit->SetTrackParam(param1);
    // Increment Chi2
    // done hit per hit, with hit resolution,
    // and not with point and angle like in "reco_muon.F" !!!!
    dX = hitForRec->GetNonBendingCoor() - param1.GetNonBendingCoor();
    dY = hitForRec->GetBendingCoor() - param1.GetBendingCoor();
    Chi2 = Chi2 + dX * dX / hitForRec->GetNonBendingReso2() + dY * dY / hitForRec->GetBendingReso2();
    trackParamAtHit = (AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->After(trackParamAtHit));
  }
}

  //__________________________________________________________________________
void TrackChi2MCS(Int_t & /*NParam*/, Double_t * /*Gradient*/, Double_t &Chi2, Double_t *Param, Int_t /*Flag*/)
{
  /// Return the "Chi2" to be minimized with Minuit for track fitting,
  /// with "NParam" parameters
  /// and their current values in array pointed to by "Param".
  /// Assumes that the track hits are sorted according to increasing Z.
  /// Track parameters at each TrackHit are updated accordingly.
  /// Multiple Coulomb scattering is taken into account with covariance matrix.
  
  AliMUONTrack *trackBeingFitted = (AliMUONTrack*) gMinuit->GetObjectFit();
//  AliMUONTrack *trackBeingFitted = (AliMUONTrack*) AliMUONTrackReconstructor::Fitter()->GetObjectFit();
  AliMUONTrackParam param1;
  AliMUONTrackParam* trackParamAtHit;
  AliMUONHitForRec* hitForRec;
  Chi2 = 0.0; // initialize Chi2
  Int_t chCurrent, chPrev = 0, hitNumber, hitNumber1, hitNumber2, hitNumber3;
  Double_t z1, z2, z3;
  AliMUONTrackParam *trackParamAtHit1, *trackParamAtHit2, *trackParamAtHit3;
  AliMUONHitForRec *hitForRec1, *hitForRec2;
  Double_t hbc1, hbc2, pbc1, pbc2;
  Double_t hnbc1, hnbc2, pnbc1, pnbc2;
  Int_t numberOfHit = trackBeingFitted->GetNTrackHits();
  TMatrixD *covBending = new TMatrixD(numberOfHit, numberOfHit);
  TMatrixD *covNonBending = new TMatrixD(numberOfHit, numberOfHit);
  Double_t *msa2 = new Double_t[numberOfHit];
  
  // copy of track parameters to be fitted
  param1 = *((AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->First()));
  param1.SetNonBendingCoor(Param[0]);
  param1.SetNonBendingSlope(Param[1]);
  param1.SetBendingCoor(Param[2]);
  param1.SetBendingSlope(Param[3]);
  param1.SetInverseBendingMomentum(Param[4]);

  // Take the vertex into account in the fit if required
  if (trackBeingFitted->GetFitWithVertex()) {
    AliMUONTrackParam paramAtVertex(param1);
    AliMUONTrackExtrap::ExtrapToZ(&paramAtVertex, 0.);
    AliMUONHitForRec *vertex = trackBeingFitted->GetVertex();
    if (!vertex) {
      cout<<"Error in TrackChi2MCS: Want to use the vertex in tracking but it has not been created!!"<<endl;
      exit(-1);
    }
    Double_t dX = vertex->GetNonBendingCoor() - paramAtVertex.GetNonBendingCoor();
    Double_t dY = vertex->GetBendingCoor() - paramAtVertex.GetBendingCoor();
    Chi2 += dX * dX / vertex->GetNonBendingReso2() + dY * dY / vertex->GetBendingReso2();
  }
  
  // Predicted coordinates and multiple scattering angles are first calculated
  for (hitNumber = 0; hitNumber < numberOfHit; hitNumber++) {
    trackParamAtHit = (AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->UncheckedAt(hitNumber));
    hitForRec = trackParamAtHit->GetHitForRecPtr();
    // extrapolation to the plane of the hitForRec attached to the current trackParamAtHit
    AliMUONTrackExtrap::ExtrapToZ(&param1, hitForRec->GetZ());
    // update track parameters of the current hit
    trackParamAtHit->SetTrackParam(param1);
    // square of multiple scattering angle at current hit, with one chamber
    msa2[hitNumber] = MultipleScatteringAngle2(&param1);
    // correction for eventual missing hits or multiple hits in a chamber,
    // according to the number of chambers
    // between the current hit and the previous one
    chCurrent = hitForRec->GetChamberNumber();
    if (hitNumber > 0) msa2[hitNumber] = msa2[hitNumber] * (chCurrent - chPrev);
    chPrev = chCurrent;
  }

  // Calculates the covariance matrix
  for (hitNumber1 = 0; hitNumber1 < numberOfHit; hitNumber1++) { 
    trackParamAtHit1 = (AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->UncheckedAt(hitNumber1));
    hitForRec1 = trackParamAtHit1->GetHitForRecPtr();
    z1 = hitForRec1->GetZ();
    for (hitNumber2 = hitNumber1; hitNumber2 < numberOfHit; hitNumber2++) {
      trackParamAtHit2 = (AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->UncheckedAt(hitNumber2));
      z2 = trackParamAtHit2->GetHitForRecPtr()->GetZ();
      // initialization to 0 (diagonal plus upper triangular part)
      (*covBending)(hitNumber2, hitNumber1) = 0.0;
      // contribution from multiple scattering in bending plane:
      // loop over upstream hits
      for (hitNumber3 = 0; hitNumber3 < hitNumber1; hitNumber3++) { 	
        trackParamAtHit3 = (AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->UncheckedAt(hitNumber3));
	z3 = trackParamAtHit3->GetHitForRecPtr()->GetZ();
	(*covBending)(hitNumber2, hitNumber1) = (*covBending)(hitNumber2, hitNumber1) + ((z1 - z3) * (z2 - z3) * msa2[hitNumber3]); 
      }
      // equal contribution from multiple scattering in non bending plane
      (*covNonBending)(hitNumber2, hitNumber1) = (*covBending)(hitNumber2, hitNumber1);
      if (hitNumber1 == hitNumber2) {
	// Diagonal elements: add contribution from position measurements
	// in bending plane
	(*covBending)(hitNumber2, hitNumber1) = (*covBending)(hitNumber2, hitNumber1) + hitForRec1->GetBendingReso2();
	// and in non bending plane
	(*covNonBending)(hitNumber2, hitNumber1) = (*covNonBending)(hitNumber2, hitNumber1) + hitForRec1->GetNonBendingReso2();
      } else {
	// Non diagonal elements: symmetrization
	// for bending plane
	(*covBending)(hitNumber1, hitNumber2) = (*covBending)(hitNumber2, hitNumber1);
	// and non bending plane
	(*covNonBending)(hitNumber1, hitNumber2) = (*covNonBending)(hitNumber2, hitNumber1);
      }
    } // for (hitNumber2 = hitNumber1;...
  } // for (hitNumber1 = 0;...
    
  // Inversion of covariance matrices
  Int_t ifailBending;
  gMinuit->mnvert(&((*covBending)(0,0)), numberOfHit, numberOfHit, numberOfHit, ifailBending);
  Int_t ifailNonBending;
  gMinuit->mnvert(&((*covNonBending)(0,0)), numberOfHit, numberOfHit, numberOfHit, ifailNonBending);

  // It would be worth trying to calculate the inverse of the covariance matrix
  // only once per fit, since it cannot change much in principle,
  // and it would save a lot of computing time !!!!
  
  // Calculates Chi2
  if ((ifailBending == 0) && (ifailNonBending == 0)) {
    // with Multiple Scattering if inversion correct
    for (hitNumber1 = 0; hitNumber1 < numberOfHit ; hitNumber1++) { 
      trackParamAtHit1 = (AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->UncheckedAt(hitNumber1));
      hitForRec1 = trackParamAtHit1->GetHitForRecPtr();
      hbc1 = hitForRec1->GetBendingCoor();
      pbc1 = trackParamAtHit1->GetBendingCoor();
      hnbc1 = hitForRec1->GetNonBendingCoor();
      pnbc1 = trackParamAtHit1->GetNonBendingCoor();
      for (hitNumber2 = 0; hitNumber2 < numberOfHit; hitNumber2++) {
	trackParamAtHit2 = (AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->UncheckedAt(hitNumber2));
        hitForRec2 = trackParamAtHit2->GetHitForRecPtr();
	hbc2 = hitForRec2->GetBendingCoor();
	pbc2 = trackParamAtHit2->GetBendingCoor();
	hnbc2 = hitForRec2->GetNonBendingCoor();
	pnbc2 = trackParamAtHit2->GetNonBendingCoor();
	Chi2 += ((*covBending)(hitNumber2, hitNumber1) * (hbc1 - pbc1) * (hbc2 - pbc2)) +
		((*covNonBending)(hitNumber2, hitNumber1) * (hnbc1 - pnbc1) * (hnbc2 - pnbc2));
      }
    }
  } else {
    // without Multiple Scattering if inversion impossible
    for (hitNumber1 = 0; hitNumber1 < numberOfHit ; hitNumber1++) { 
      trackParamAtHit1 = (AliMUONTrackParam*) (trackBeingFitted->GetTrackParamAtHit()->UncheckedAt(hitNumber1));
      hitForRec1 = trackParamAtHit1->GetHitForRecPtr();
      hbc1 = hitForRec1->GetBendingCoor();
      pbc1 = trackParamAtHit1->GetBendingCoor();
      hnbc1 = hitForRec1->GetNonBendingCoor();
      pnbc1 = trackParamAtHit1->GetNonBendingCoor();
      Chi2 += ((hbc1 - pbc1) * (hbc1 - pbc1) / hitForRec1->GetBendingReso2()) +
	      ((hnbc1 - pnbc1) * (hnbc1 - pnbc1) / hitForRec1->GetNonBendingReso2());
    }
  }
  
  delete covBending;
  delete covNonBending;
  delete [] msa2;
}

  //__________________________________________________________________________
Double_t MultipleScatteringAngle2(AliMUONTrackParam *param)
{
  /// Returns square of multiple Coulomb scattering angle
  /// from TrackParamAtHit pointed to by "param"
  Double_t slopeBending, slopeNonBending, radiationLength, inverseBendingMomentum2, inverseTotalMomentum2;
  Double_t varMultipleScatteringAngle;
  slopeBending = param->GetBendingSlope();
  slopeNonBending = param->GetNonBendingSlope();
  // thickness in radiation length for the current track,
  // taking local angle into account
  radiationLength = AliMUONConstants::ChamberThicknessInX0() *
		    TMath::Sqrt(1.0 + slopeBending*slopeBending + slopeNonBending*slopeNonBending);
  inverseBendingMomentum2 =  param->GetInverseBendingMomentum() * param->GetInverseBendingMomentum();
  inverseTotalMomentum2 = inverseBendingMomentum2 * (1.0 + slopeBending * slopeBending) /
			  (1.0 + slopeBending *slopeBending + slopeNonBending * slopeNonBending); 
  varMultipleScatteringAngle = 0.0136 * (1.0 + 0.038 * TMath::Log(radiationLength));
  // The velocity is assumed to be 1 !!!!
  varMultipleScatteringAngle = inverseTotalMomentum2 * radiationLength * varMultipleScatteringAngle * varMultipleScatteringAngle;
  return varMultipleScatteringAngle;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::ExtrapTracksToVertex()
{
  /// propagate tracks to the vertex through the absorber (Branson)
  AliMUONTrack *track;
  AliMUONTrackParam trackParamVertex;
  AliMUONHitForRec *vertex;
  Double_t vX, vY, vZ;
  
  for (Int_t iRecTrack = 0; iRecTrack <fNRecTracks; iRecTrack++) {
    track = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iRecTrack);
    trackParamVertex = *((AliMUONTrackParam*)(track->GetTrackParamAtHit()->First()));
    vertex = track->GetVertex();
    if (vertex) { // if track as a vertex defined, use it
      vX = vertex->GetNonBendingCoor();
      vY = vertex->GetBendingCoor();
      vZ = vertex->GetZ();
    } else { // else vertex = (0.,0.,0.)
      vX = 0.;
      vY = 0.;
      vZ = 0.;
    }
    AliMUONTrackExtrap::ExtrapToVertex(&trackParamVertex, vX, vY, vZ);
    track->SetTrackParamAtVertex(&trackParamVertex);
    
    if (AliLog::GetGlobalDebugLevel() > 0) {
      cout << "FollowTracks: track candidate(0..): " << iRecTrack
  	   << " after extrapolation to vertex" << endl;
      track->RecursiveDump();
    }
  }
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::FillMUONTrack()
{
  /// Fill fHitForRecAtHit of AliMUONTrack's
  AliMUONTrack *track;
  AliMUONTrackParam *trackParamAtHit;
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    trackParamAtHit = (AliMUONTrackParam*) (track->GetTrackParamAtHit()->First());
    while (trackParamAtHit) {
      track->AddHitForRecAtHit(trackParamAtHit->GetHitForRecPtr());
      trackParamAtHit = (AliMUONTrackParam*) (track->GetTrackParamAtHit()->After(trackParamAtHit)); 
    }
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
  }
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::EventDump(void)
{
  /// Dump reconstructed event (track parameters at vertex and at first hit),
  /// and the particle parameters
  AliMUONTrack *track;
  AliMUONTrackParam *trackParam, *trackParam1;
  Double_t bendingSlope, nonBendingSlope, pYZ;
  Double_t pX, pY, pZ, x, y, z, c;
  Int_t trackIndex, nTrackHits;
 
  AliDebug(1,"****** enter EventDump ******");
  AliDebug(1, Form("Number of Reconstructed tracks : %d", fNRecTracks)); 
  
  fRecTracksPtr->Compress(); // for simple loop without "Next" since no hole
  // Loop over reconstructed tracks
  for (trackIndex = 0; trackIndex < fNRecTracks; trackIndex++) {
    AliDebug(1, Form("track number: %d", trackIndex));
    // function for each track for modularity ????
    track = (AliMUONTrack*) ((*fRecTracksPtr)[trackIndex]);
    nTrackHits = track->GetNTrackHits();
    AliDebug(1, Form("Number of track hits: %d ", nTrackHits));
    // track parameters at Vertex
    trackParam = track->GetTrackParamAtVertex();
    x = trackParam->GetNonBendingCoor();
    y = trackParam->GetBendingCoor();
    z = trackParam->GetZ();
    bendingSlope = trackParam->GetBendingSlope();
    nonBendingSlope = trackParam->GetNonBendingSlope();
    pYZ = 1/TMath::Abs(trackParam->GetInverseBendingMomentum());
    pZ = pYZ/TMath::Sqrt(1+bendingSlope*bendingSlope);
    pX = pZ * nonBendingSlope;
    pY = pZ * bendingSlope;
    c = TMath::Sign(1.0, trackParam->GetInverseBendingMomentum());
    AliDebug(1, Form("Track parameters at Vertex z= %f: X= %f Y= %f pX= %f pY= %f pZ= %f c= %f\n",
		     z, x, y, pX, pY, pZ, c));

    // track parameters at first hit
    trackParam1 = (AliMUONTrackParam*) track->GetTrackParamAtHit()->First();
    x = trackParam1->GetNonBendingCoor();
    y = trackParam1->GetBendingCoor();
    z = trackParam1->GetZ();
    bendingSlope = trackParam1->GetBendingSlope();
    nonBendingSlope = trackParam1->GetNonBendingSlope();
    pYZ = 1/TMath::Abs(trackParam1->GetInverseBendingMomentum());
    pZ = pYZ/TMath::Sqrt(1.0 + bendingSlope * bendingSlope);
    pX = pZ * nonBendingSlope;
    pY = pZ * bendingSlope;
    c = TMath::Sign(1.0, trackParam1->GetInverseBendingMomentum());
    AliDebug(1, Form("track parameters at z= %f: X= %f Y= %f pX= %f pY= %f pZ= %f c= %f\n",
		     z, x, y, pX, pY, pZ, c));
  }
  // informations about generated particles NO !!!!!!!!
  
//    for (Int_t iPart = 0; iPart < np; iPart++) {
//      p = gAlice->Particle(iPart);
//      printf(" particle %d: type= %d px= %f py= %f pz= %f pdg= %d\n",
//  	   iPart, p->GetPdgCode(), p->Px(), p->Py(), p->Pz(), p->GetPdgCode());    
//    }
  return;
}


