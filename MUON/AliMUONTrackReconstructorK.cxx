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
/// \class AliMUONTrackReconstructorK
///
/// MUON track reconstructor using the kalman method
///
/// This class contains as data:
/// - the parameters for the track reconstruction
///
/// It contains as methods, among others:
/// - MakeTracks to build the tracks
///
//-----------------------------------------------------------------------------

#include "AliMUONTrackReconstructorK.h"
#include "AliMUONConstants.h"
#include "AliMUONHitForRec.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TMath.h>
#include <TMatrixD.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrackReconstructorK) // Class implementation in ROOT context
/// \endcond

  //__________________________________________________________________________
AliMUONTrackReconstructorK::AliMUONTrackReconstructorK()
  : AliMUONVTrackReconstructor()
{
  /// Constructor for class AliMUONTrackReconstructorK
  AliInfo("*** Tracking with Kalman Filter ***");
}

  //__________________________________________________________________________
AliMUONTrackReconstructorK::~AliMUONTrackReconstructorK()
{
/// Destructor
} 

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::MakeTrackCandidates()
{
  /// To make track candidates (assuming linear propagation if the flag fgkMakeTrackCandidatesFast is set to kTRUE):
  /// Start with segments station(1..) 4 or 5 then follow track in station 5 or 4.
  /// Good candidates are made of at least three hitForRec's.
  /// Keep only best candidates or all of them according to the flag fgkTrackAllTracks.
  
  TClonesArray *segments;
  AliMUONTrack *track;
  Int_t iCandidate = 0;
  Bool_t hitFound;

  AliDebug(1,"Enter MakeTrackCandidates");

  // Loop over stations(1..) 5 and 4 and make track candidates
  for (Int_t istat=4; istat>=3; istat--) {
    
    // Make segments in the station
    segments = MakeSegmentsInStation(istat);
    
    // Loop over segments
    for (Int_t iseg=0; iseg<segments->GetEntriesFast(); iseg++) {
      AliDebug(1,Form("Making primary candidate(1..) %d",++iCandidate));
      
      // Transform segments to tracks and put them at the end of fRecTracksPtr
      track = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack((AliMUONObjectPair*)((*segments)[iseg]));
      fNRecTracks++;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
        cout<<endl<<"Track parameter covariances at first hit:"<<endl;
        ((AliMUONTrackParam*) (track->GetTrackParamAtHit()->First()))->GetCovariances().Print();
      }
      
      // Look for compatible hitForRec(s) in the other station
      if (AliMUONReconstructor::GetRecoParam()->MakeTrackCandidatesFast()) hitFound = FollowLinearTrackInStation(*track,7-istat);
      else {
        // First recompute track parameters and covariances on station(1..) 5 using Kalman filter
        // (to make sure all tracks are treated in the same way)
        if (istat == 4) RetraceTrack(*track,kFALSE);
	hitFound = FollowTrackInStation(*track,7-istat);
      }
      
      // Remove track if no hit found
      if (!hitFound) {
        fRecTracksPtr->Remove(track);
  	fNRecTracks--;
      }
      
    }
    // delete the array of segments
    delete segments;
  }
  
  
  // Retrace tracks using Kalman filter and remove bad ones
  if (AliMUONReconstructor::GetRecoParam()->MakeTrackCandidatesFast()) {
    fRecTracksPtr->Compress(); // this is essential before checking tracks
    
    Int_t currentNRecTracks = fNRecTracks;
    for (Int_t iRecTrack = 0; iRecTrack < currentNRecTracks; iRecTrack++) {
      track = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iRecTrack);
      
      // Recompute track parameters and covariances using Kalman filter
      RetraceTrack(*track,kTRUE);
      
      // Remove the track if the normalized chi2 is too high
      if (track->GetNormalizedChi2() > AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking() *
	                               AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking()) {
        fRecTracksPtr->Remove(track);
        fNRecTracks--;
      }
      
    }
    
  }
  
  fRecTracksPtr->Compress(); // this is essential before checking tracks
  
  // Keep all different tracks or only the best ones as required
  if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) RemoveIdenticalTracks();
  else RemoveDoubleTracks();
  
  AliDebug(1,Form("Number of good candidates = %d",fNRecTracks));
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::RetraceTrack(AliMUONTrack &trackCandidate, Bool_t resetSeed)
{
  /// Re-run the kalman filter from the most downstream hit to the most uptream one
  AliDebug(1,"Enter RetraceTrack");
  
  AliMUONTrackParam* startingTrackParam = (AliMUONTrackParam*) trackCandidate.GetTrackParamAtHit()->Last();
  
  // Reset the "seed" (= track parameters and their covariances at last hit) if required
  if (resetSeed) {
    // => Shift track parameters at the position of the last hit
    AliMUONHitForRec* hitForRecAtHit = startingTrackParam->GetHitForRecPtr();
    startingTrackParam->SetNonBendingCoor(hitForRecAtHit->GetNonBendingCoor());
    startingTrackParam->SetBendingCoor(hitForRecAtHit->GetBendingCoor());
    
    // => Re-compute and reset track parameters covariances at last hit (as if the other hits did not exist)
    const TMatrixD& kParamCov = startingTrackParam->GetCovariances();
    TMatrixD newParamCov(5,5);
    newParamCov.Zero();
    // Non bending plane
    newParamCov(0,0) = hitForRecAtHit->GetNonBendingReso2();
    newParamCov(1,1) = 100.*kParamCov(1,1);
    // Bending plane
    newParamCov(2,2) = hitForRecAtHit->GetBendingReso2();
    newParamCov(3,3) = 100.*kParamCov(3,3);
    // Inverse bending momentum
    newParamCov(4,4) = 0.5*startingTrackParam->GetInverseBendingMomentum() * 0.5*startingTrackParam->GetInverseBendingMomentum();
    startingTrackParam->SetCovariances(newParamCov);
    
    // Reset the track chi2
    startingTrackParam->SetTrackChi2(0.);
  
  }
  
  // Redo the tracking
  RetracePartialTrack(trackCandidate, startingTrackParam);
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::RetracePartialTrack(AliMUONTrack &trackCandidate, const AliMUONTrackParam* startingTrackParam)
{
  /// Re-run the kalman filter from the hit attached to startingTrackParam to the most uptream hit
  AliDebug(1,"Enter RetracePartialTrack");
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "RetracePartialTrack: track chi2 before re-tracking: " << trackCandidate.GetFitFMin() << endl;
  }
  
  // Reset the track chi2
  trackCandidate.SetFitFMin(startingTrackParam->GetTrackChi2());
  
  // loop over attached hits until the first one and recompute track parameters and covariances using kalman filter
  Int_t expectedChamber = startingTrackParam->GetHitForRecPtr()->GetChamberNumber() - 1;
  Int_t currentChamber;
  Double_t addChi2TrackAtHit;
  AliMUONTrackParam* trackParamAtHit = (AliMUONTrackParam*) trackCandidate.GetTrackParamAtHit()->Before(startingTrackParam); 
  while (trackParamAtHit) {
    
    // reset track parameters and their covariances
    trackParamAtHit->SetParameters(startingTrackParam->GetParameters());
    trackParamAtHit->SetZ(startingTrackParam->GetZ());
    trackParamAtHit->SetCovariances(startingTrackParam->GetCovariances());
    
    // add MCS effect
    AliMUONTrackExtrap::AddMCSEffect(trackParamAtHit,AliMUONConstants::ChamberThicknessInX0(),1.);
    
    // reset propagator for smoother
    if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) trackParamAtHit->ResetPropagator();
    
    // add MCS in missing chambers if any (at most 2 chambers can be missing according to tracking criteria)
    currentChamber = trackParamAtHit->GetHitForRecPtr()->GetChamberNumber();
    while (currentChamber < expectedChamber) {
      // extrapolation to the missing chamber (update the propagator)
      AliMUONTrackExtrap::ExtrapToZCov(trackParamAtHit, AliMUONConstants::DefaultChamberZ(expectedChamber),
				       AliMUONReconstructor::GetRecoParam()->UseSmoother());
      // add MCS effect
      AliMUONTrackExtrap::AddMCSEffect(trackParamAtHit,AliMUONConstants::ChamberThicknessInX0(),1.);
      expectedChamber--;
    }
    
    // extrapolation to the plane of the hitForRec attached to the current trackParamAtHit (update the propagator)
    AliMUONTrackExtrap::ExtrapToZCov(trackParamAtHit, trackParamAtHit->GetHitForRecPtr()->GetZ(),
				     AliMUONReconstructor::GetRecoParam()->UseSmoother());
    
    if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) {
      // save extrapolated parameters for smoother
      trackParamAtHit->SetExtrapParameters(trackParamAtHit->GetParameters());
      
      // save extrapolated covariance matrix for smoother
      trackParamAtHit->SetExtrapCovariances(trackParamAtHit->GetCovariances());
    }
    
    // Compute new track parameters including "hitForRecCh2" using kalman filter
    addChi2TrackAtHit = RunKalmanFilter(*trackParamAtHit);
    
    // Update the track chi2
    trackCandidate.SetFitFMin(trackCandidate.GetFitFMin() + addChi2TrackAtHit);
    trackParamAtHit->SetTrackChi2(trackCandidate.GetFitFMin());
    
    // prepare next step, add MCS effects in parameter covariances
    expectedChamber--;
    startingTrackParam = trackParamAtHit;
    trackParamAtHit = (AliMUONTrackParam*) (trackCandidate.GetTrackParamAtHit()->Before(startingTrackParam)); 
  }
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "RetracePartialTrack: track chi2 after re-tracking: " << trackCandidate.GetFitFMin() << endl;
  }
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::FollowTracks()
{
  /// Follow tracks in stations(1..) 3, 2 and 1
  AliDebug(1,"Enter FollowTracks");
  
  AliMUONTrack *track;
  Int_t currentNRecTracks;
  Bool_t hitFound;
  
  for (Int_t station = 2; station >= 0; station--) {
    
    // Save the actual number of reconstructed track in case of
    // tracks are added or suppressed during the tracking procedure
    // !! Do not compress fRecTracksPtr until the end of the loop over tracks !!
    currentNRecTracks = fNRecTracks;
    
    for (Int_t iRecTrack = 0; iRecTrack <currentNRecTracks; iRecTrack++) {
      AliDebug(1,Form("FollowTracks: track candidate(1..) %d", iRecTrack+1));
      
      track = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iRecTrack);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
        cout<<endl<<"Track parameter covariances at first hit:"<<endl;
        ((AliMUONTrackParam*) (track->GetTrackParamAtHit()->First()))->GetCovariances().Print();
      }
      
      // Look for compatible hitForRec in station(0..) "station"
      hitFound = FollowTrackInStation(*track,station);
      
      // Try to recover track if required
      if (!hitFound && AliMUONReconstructor::GetRecoParam()->RecoverTracks()) hitFound = RecoverTrack(*track,station);
      
      // remove track if no hit found
      if (!hitFound) {
	fRecTracksPtr->Remove(track);
  	fNRecTracks--;
      }
      
    }
    
    fRecTracksPtr->Compress(); // this is essential before checking tracks
    
    // Keep only the best tracks if required
    if (!AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) RemoveDoubleTracks();
    
  }
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructorK::FollowTrackInStation(AliMUONTrack &trackCandidate, Int_t nextStation)
{
  /// Follow trackCandidate in station(0..) nextStation and search for compatible HitForRec(s)
  /// Keep all possibilities or only the best one(s) according to the flag fgkTrackAllTracks:
  /// kTRUE:  duplicate "trackCandidate" if there are several possibilities and add the new tracks at the end of
  ///         fRecTracksPtr to avoid conficts with other track candidates at this current stage of the tracking procedure.
  ///         Remove the obsolete "trackCandidate" at the end.
  /// kFALSE: add only the best hit(s) to the "trackCandidate". Try to add a couple of hits in priority.
  /// return kTRUE if new hits have been found (otherwise return kFALSE)
  AliDebug(1,Form("Enter FollowTrackInStation(1..) %d", nextStation+1));
  
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
  
  Double_t chi2OfHitForRec;
  Double_t maxChi2OfHitForRec = 2. * AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking() *
				     AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking(); // 2 because 2 quantities in chi2
  Double_t addChi2TrackAtHit1;
  Double_t addChi2TrackAtHit2;
  Double_t bestAddChi2TrackAtHit1 = 1.e10;
  Double_t bestAddChi2TrackAtHit2 = 1.e10;
  Bool_t foundOneHit = kFALSE;
  Bool_t foundTwoHits = kFALSE;
  AliMUONTrack *newTrack = 0x0;
  AliMUONHitForRec *hitForRecCh1, *hitForRecCh2;
  AliMUONTrackParam extrapTrackParam;
  AliMUONTrackParam extrapTrackParamAtHit1;
  AliMUONTrackParam extrapTrackParamAtHit2;
  AliMUONTrackParam bestTrackParamAtHit1;
  AliMUONTrackParam bestTrackParamAtHit2;
  Bool_t *hitForRecCh1Used = new Bool_t[fNHitsForRecPerChamber[ch1]];
  for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) hitForRecCh1Used[hit1] = kFALSE;

  // Get track parameters
  AliMUONTrackParam extrapTrackParamAtCh(*(AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->First());
  
  // Add MCS effect
  AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(),1.);
  
  // reset propagator for smoother
  if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) extrapTrackParamAtCh.ResetPropagator();
  
  // Add MCS in the missing chamber if any (only 1 chamber can be missing according to tracking criteria)
  if (ch1 < ch2 && extrapTrackParamAtCh.GetHitForRecPtr()->GetChamberNumber() > ch2 + 1) {
    // extrapolation to the missing chamber
    AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(ch2 + 1),
				     AliMUONReconstructor::GetRecoParam()->UseSmoother());
    // add MCS effect
    AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(),1.);
  }
  
  //Extrapolate trackCandidate to chamber "ch2"
  AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(ch2),
				   AliMUONReconstructor::GetRecoParam()->UseSmoother());
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    cout<<endl<<"Track parameter covariances at first hit extrapolated to z = "<<AliMUONConstants::DefaultChamberZ(ch2)<<":"<<endl;
    extrapTrackParamAtCh.GetCovariances().Print();
  }
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "FollowTrackInStation: look for hits in chamber(1..): " << ch2+1 << endl;
  }
  
  // look for candidates in chamber 2 
  for (Int_t hit2 = 0; hit2 < fNHitsForRecPerChamber[ch2]; hit2++) {
    
    hitForRecCh2 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch2]+hit2);
    
    // try to add the current hit fast
    if (!TryOneHitForRecFast(extrapTrackParamAtCh, hitForRecCh2)) continue;
    
    // try to add the current hit accuratly
    chi2OfHitForRec = TryOneHitForRec(extrapTrackParamAtCh, hitForRecCh2, extrapTrackParamAtHit2,
				      AliMUONReconstructor::GetRecoParam()->UseSmoother());
    
    // if good chi2 then try to attach a hitForRec in the other chamber too
    if (chi2OfHitForRec < maxChi2OfHitForRec) {
      Bool_t foundSecondHit = kFALSE;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: found one hit in chamber(1..): " << ch2+1
	     << " (Chi2 = " << chi2OfHitForRec << ")" << endl;
        cout << "                      look for second hits in chamber(1..): " << ch1+1 << " ..." << endl;
      }
      
      if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) {
        // save extrapolated parameters for smoother
        extrapTrackParamAtHit2.SetExtrapParameters(extrapTrackParamAtHit2.GetParameters());
        
        // save extrapolated covariance matrix for smoother
        extrapTrackParamAtHit2.SetExtrapCovariances(extrapTrackParamAtHit2.GetCovariances());
      }
      
      // Compute new track parameters including "hitForRecCh2" using kalman filter
      addChi2TrackAtHit2 = RunKalmanFilter(extrapTrackParamAtHit2);
      
      // copy new track parameters for next step
      extrapTrackParam = extrapTrackParamAtHit2;
      
      // add MCS effect
      AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParam,AliMUONConstants::ChamberThicknessInX0(),1.);
      
      // reset propagator for smoother
      if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) extrapTrackParam.ResetPropagator();
      
      //Extrapolate track parameters to chamber "ch1"
      AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParam, AliMUONConstants::DefaultChamberZ(ch1),
				       AliMUONReconstructor::GetRecoParam()->UseSmoother());
      
      for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) {
        
	hitForRecCh1 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch1]+hit1);
	
    	// try to add the current hit fast
    	if (!TryOneHitForRecFast(extrapTrackParam, hitForRecCh1)) continue;
    	
    	// try to add the current hit accuratly
    	chi2OfHitForRec = TryOneHitForRec(extrapTrackParam, hitForRecCh1, extrapTrackParamAtHit1,
					  AliMUONReconstructor::GetRecoParam()->UseSmoother());
    	
	// if good chi2 then consider to add the 2 hitForRec to the "trackCandidate"
	if (chi2OfHitForRec < maxChi2OfHitForRec) {
	  foundSecondHit = kTRUE;
	  foundTwoHits = kTRUE;
          
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowTrackInStation: found one hit in chamber(1..): " << ch1+1
	  	 << " (Chi2 = " << chi2OfHitForRec << ")" << endl;
	  }
          
          if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) {
            // save extrapolated parameters for smoother
            extrapTrackParamAtHit1.SetExtrapParameters(extrapTrackParamAtHit1.GetParameters());
            
            // save extrapolated covariance matrix for smoother
            extrapTrackParamAtHit1.SetExtrapCovariances(extrapTrackParamAtHit1.GetCovariances());
          }
          
          // Compute new track parameters including "hitForRecCh1" using kalman filter
          addChi2TrackAtHit1 = RunKalmanFilter(extrapTrackParamAtHit1);
          
	  if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
	    // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
            newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	    UpdateTrack(*newTrack,extrapTrackParamAtHit1,extrapTrackParamAtHit2,addChi2TrackAtHit1,addChi2TrackAtHit2);
	    fNRecTracks++;
	    
	    // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
	    // (going in the right direction)
	    if (nextStation == 4) RetraceTrack(*newTrack,kTRUE);
	    
	    // Tag hitForRecCh1 as used
	    hitForRecCh1Used[hit1] = kTRUE;
	    
	    // Printout for debuging
	    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	      cout << "FollowTrackInStation: added two hits in station(1..): " << nextStation+1 << endl;
	      if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	    }
	    
          } else if (addChi2TrackAtHit1+addChi2TrackAtHit2 < bestAddChi2TrackAtHit1+bestAddChi2TrackAtHit2) {
	    // keep track of the best couple of hits
	    bestAddChi2TrackAtHit1 = addChi2TrackAtHit1;
	    bestAddChi2TrackAtHit2 = addChi2TrackAtHit2;
	    bestTrackParamAtHit1 = extrapTrackParamAtHit1;
	    bestTrackParamAtHit2 = extrapTrackParamAtHit2;
          }
	  
	}
	
      }
      
      // if no hitForRecCh1 found then consider to add hitForRecCh2 only
      if (!foundSecondHit) {
	foundOneHit = kTRUE;
        
	if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
          newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  UpdateTrack(*newTrack,extrapTrackParamAtHit2,addChi2TrackAtHit2);
	  fNRecTracks++;
	  
	  // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
	  // (going in the right direction)
	  if (nextStation == 4) RetraceTrack(*newTrack,kTRUE);
	  
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowTrackInStation: added one hit in chamber(1..): " << ch2+1 << endl;
	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	  }
	  
	} else if (!foundTwoHits && addChi2TrackAtHit2 < bestAddChi2TrackAtHit1) {
	  // keep track of the best single hitForRec except if a couple of hits has already been found
	  bestAddChi2TrackAtHit1 = addChi2TrackAtHit2;
	  bestTrackParamAtHit1 = extrapTrackParamAtHit2;
        }
	
      }
      
    }
    
  }
  
  // look for candidates in chamber 1 not already attached to a track
  // if we want to keep all possible tracks or if no good couple of hitForRec has been found
  if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks() || !foundTwoHits) {
    
    // Printout for debuging
    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
      cout << "FollowTrackInStation: look for single hits in chamber(1..): " << ch1+1 << endl;
    }
    
    // add MCS effect for next step
    AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(),1.);
    
    //Extrapolate trackCandidate to chamber "ch1"
    AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(ch1),
				     AliMUONReconstructor::GetRecoParam()->UseSmoother());
    
    for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) {
      
      hitForRecCh1 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch1]+hit1);
      
      if (hitForRecCh1Used[hit1]) continue; // Skip hitForRec already used
      
      // try to add the current hit fast
      if (!TryOneHitForRecFast(extrapTrackParamAtCh, hitForRecCh1)) continue;
      
      // try to add the current hit accuratly
      chi2OfHitForRec = TryOneHitForRec(extrapTrackParamAtCh, hitForRecCh1, extrapTrackParamAtHit1,
					AliMUONReconstructor::GetRecoParam()->UseSmoother());
      
      // if good chi2 then consider to add hitForRecCh1
      // We do not try to attach a hitForRec in the other chamber too since it has already been done above
      if (chi2OfHitForRec < maxChi2OfHitForRec) {
	foundOneHit = kTRUE;
  	
	// Printout for debuging
  	if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	  cout << "FollowTrackInStation: found one hit in chamber(1..): " << ch1+1
  	       << " (Chi2 = " << chi2OfHitForRec << ")" << endl;
  	}
        
	if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) {
          // save extrapolated parameters for smoother
          extrapTrackParamAtHit1.SetExtrapParameters(extrapTrackParamAtHit1.GetParameters());
          
          // save extrapolated covariance matrix for smoother
          extrapTrackParamAtHit1.SetExtrapCovariances(extrapTrackParamAtHit1.GetCovariances());
        }
        
        // Compute new track parameters including "hitForRecCh1" using kalman filter
        addChi2TrackAtHit1 = RunKalmanFilter(extrapTrackParamAtHit1);
        
	if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
  	  newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  UpdateTrack(*newTrack,extrapTrackParamAtHit1,addChi2TrackAtHit1);
	  fNRecTracks++;
  	  
	  // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
	  // (going in the right direction)
	  if (nextStation == 4) RetraceTrack(*newTrack,kTRUE);
	  
	  // Printout for debuging
  	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	    cout << "FollowTrackInStation: added one hit in chamber(1..): " << ch1+1 << endl;
  	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
  	  }
	  
  	} else if (addChi2TrackAtHit1 < bestAddChi2TrackAtHit1) {
	  // keep track of the best single hitForRec except if a couple of hits has already been found
	  bestAddChi2TrackAtHit1 = addChi2TrackAtHit1;
	  bestTrackParamAtHit1 = extrapTrackParamAtHit1;
  	}
	
      }
      
    }
    
  }
  
  // fill out the best track if required else clean up the fRecTracksPtr array
  if (!AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
    if (foundTwoHits) {
      UpdateTrack(trackCandidate,bestTrackParamAtHit1,bestTrackParamAtHit2,bestAddChi2TrackAtHit1,bestAddChi2TrackAtHit2);
      
      // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
      // (going in the right direction)
      if (nextStation == 4) RetraceTrack(trackCandidate,kTRUE);

      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the two best hits in station(1..): " << nextStation+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else if (foundOneHit) {
      UpdateTrack(trackCandidate,bestTrackParamAtHit1,bestAddChi2TrackAtHit1);
      
      // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
      // (going in the right direction)
      if (nextStation == 4) RetraceTrack(trackCandidate,kTRUE);

      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the best hit in chamber(1..): " << bestTrackParamAtHit1.GetHitForRecPtr()->GetChamberNumber()+1 << endl;
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
Double_t AliMUONTrackReconstructorK::RunKalmanFilter(AliMUONTrackParam &trackParamAtHit)
{
  /// Compute new track parameters and their covariances including new hit using kalman filter
  /// return the additional track chi2
  AliDebug(1,"Enter RunKalmanFilter");
  
  // Get actual track parameters (p)
  TMatrixD param(trackParamAtHit.GetParameters());
  
  // Get new hit parameters (m)
  AliMUONHitForRec *hitForRecAtHit = trackParamAtHit.GetHitForRecPtr();
  TMatrixD hit(5,1);
  hit.Zero();
  hit(0,0) = hitForRecAtHit->GetNonBendingCoor();
  hit(2,0) = hitForRecAtHit->GetBendingCoor();
  
  // Compute the actual parameter weight (W)
  TMatrixD paramWeight(trackParamAtHit.GetCovariances());
  if (paramWeight.Determinant() != 0) {
    paramWeight.Invert();
  } else {
    AliWarning(" Determinant = 0");
    return 1.e10;
  }
  
  // Compute the new hit weight (U)
  TMatrixD hitWeight(5,5);
  hitWeight.Zero();
  hitWeight(0,0) = 1. / hitForRecAtHit->GetNonBendingReso2();
  hitWeight(2,2) = 1. / hitForRecAtHit->GetBendingReso2();

  // Compute the new parameters covariance matrix ( (W+U)^-1 )
  TMatrixD newParamCov(paramWeight,TMatrixD::kPlus,hitWeight);
  if (newParamCov.Determinant() != 0) {
    newParamCov.Invert();
  } else {
    AliWarning(" Determinant = 0");
    return 1.e10;
  }
  
  // Save the new parameters covariance matrix
  trackParamAtHit.SetCovariances(newParamCov);
  
  // Compute the new parameters (p' = ((W+U)^-1)U(m-p) + p)
  TMatrixD tmp(hit,TMatrixD::kMinus,param);
  TMatrixD tmp2(hitWeight,TMatrixD::kMult,tmp); // U(m-p)
  TMatrixD newParam(newParamCov,TMatrixD::kMult,tmp2); // ((W+U)^-1)U(m-p)
  newParam += param; // ((W+U)^-1)U(m-p) + p
  
  // Save the new parameters
  trackParamAtHit.SetParameters(newParam);
  
  // Compute the additional chi2 (= ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m))
  tmp = newParam; // p'
  tmp -= param; // (p'-p)
  TMatrixD tmp3(paramWeight,TMatrixD::kMult,tmp); // W(p'-p)
  TMatrixD addChi2Track(tmp,TMatrixD::kTransposeMult,tmp3); // ((p'-p)^-1)W(p'-p)
  tmp = newParam; // p'
  tmp -= hit; // (p'-m)
  TMatrixD tmp4(hitWeight,TMatrixD::kMult,tmp); // U(p'-m)
  addChi2Track += TMatrixD(tmp,TMatrixD::kTransposeMult,tmp4); // ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m)
  
  return addChi2Track(0,0);
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtHit, Double_t addChi2)
{
  /// Add 1 hit to the track candidate
  /// Update chi2 of the track 
  
  // Flag hit as being not removable
  trackParamAtHit.SetRemovable(kFALSE);
  trackParamAtHit.SetLocalChi2(0.); // --> Local chi2 not used
  
  // Update the track chi2 into TrackParamAtHit
  trackParamAtHit.SetTrackChi2(track.GetFitFMin() + addChi2);
  
  // Update the chi2 of the new track
  track.SetFitFMin(trackParamAtHit.GetTrackChi2());
  
  // Update array of TrackParamAtHit
  track.AddTrackParamAtHit(&trackParamAtHit,trackParamAtHit.GetHitForRecPtr());
  track.GetTrackParamAtHit()->Sort();
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtHit1, AliMUONTrackParam &trackParamAtHit2,
					     Double_t addChi2AtHit1, Double_t addChi2AtHit2)
{
  /// Add 2 hits to the track candidate
  /// Update track and local chi2
  
  // Update local chi2 at first hit
  AliMUONHitForRec* hit1 = trackParamAtHit1.GetHitForRecPtr();
  Double_t deltaX = trackParamAtHit1.GetNonBendingCoor() - hit1->GetNonBendingCoor();
  Double_t deltaY = trackParamAtHit1.GetBendingCoor() - hit1->GetBendingCoor();
  Double_t localChi2AtHit1 = deltaX*deltaX / hit1->GetNonBendingReso2() +
  			     deltaY*deltaY / hit1->GetBendingReso2();
  trackParamAtHit1.SetLocalChi2(localChi2AtHit1);
  
  // Flag first hit as being removable
  trackParamAtHit1.SetRemovable(kTRUE);
  
  // Update local chi2 at second hit
  AliMUONHitForRec* hit2 = trackParamAtHit2.GetHitForRecPtr();
  AliMUONTrackParam extrapTrackParamAtHit2(trackParamAtHit1);
  AliMUONTrackExtrap::ExtrapToZ(&extrapTrackParamAtHit2, trackParamAtHit2.GetZ());
  deltaX = extrapTrackParamAtHit2.GetNonBendingCoor() - hit2->GetNonBendingCoor();
  deltaY = extrapTrackParamAtHit2.GetBendingCoor() - hit2->GetBendingCoor();
  Double_t localChi2AtHit2 = deltaX*deltaX / hit2->GetNonBendingReso2() +
  			     deltaY*deltaY / hit2->GetBendingReso2();
  trackParamAtHit2.SetLocalChi2(localChi2AtHit2);
  
  // Flag second hit as being removable
  trackParamAtHit2.SetRemovable(kTRUE);
  
  // Update the track chi2 into TrackParamAtHit1
  trackParamAtHit1.SetTrackChi2(track.GetFitFMin() + addChi2AtHit1);
  
  // Update the track chi2 into TrackParamAtHit2
  trackParamAtHit2.SetTrackChi2(trackParamAtHit1.GetTrackChi2() + addChi2AtHit2);
  
  // Update the chi2 of the new track
  track.SetFitFMin(trackParamAtHit2.GetTrackChi2());
  
  // Update array of TrackParamAtHit
  track.AddTrackParamAtHit(&trackParamAtHit1,trackParamAtHit1.GetHitForRecPtr());
  track.AddTrackParamAtHit(&trackParamAtHit2,trackParamAtHit2.GetHitForRecPtr());
  track.GetTrackParamAtHit()->Sort();
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructorK::RecoverTrack(AliMUONTrack &trackCandidate, Int_t nextStation)
{
  /// Try to recover the track candidate in the next station
  /// by removing the worst of the two hits attached in the current station
  /// Return kTRUE if recovering succeeds
  AliDebug(1,"Enter RecoverTrack");
  
  // Do not try to recover track until we have attached hit(s) on station(1..) 3
  if (nextStation > 1) return kFALSE;
  
  Int_t worstHitNumber = -1;
  Double_t localChi2, worstChi2 = 0.;
  
  // Look for the hit to remove
  for (Int_t hitNumber = 0; hitNumber < 2; hitNumber++) {
    AliMUONTrackParam *trackParamAtHit = (AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->UncheckedAt(hitNumber);
    
    // check if current hit is removable
    if (!trackParamAtHit->IsRemovable()) return kFALSE;
    
    // Pick up hit with the worst chi2
    localChi2 = trackParamAtHit->GetLocalChi2();
    if (localChi2 > worstChi2) {
      worstChi2 = localChi2;
      worstHitNumber = hitNumber;
    }
  }
  
  // Reset best hit as being NOT removable
  ((AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->UncheckedAt((worstHitNumber+1)%2))->SetRemovable(kFALSE);
  
  // Remove the worst hit
  trackCandidate.RemoveTrackParamAtHit((AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->UncheckedAt(worstHitNumber));
  
  // Re-calculate track parameters at the (new) first hit
  RetracePartialTrack(trackCandidate,(AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->UncheckedAt(1));
  
  // Look for new hit(s) in next station
  return FollowTrackInStation(trackCandidate,nextStation);
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructorK::RunSmoother(AliMUONTrack &track)
{
  /// Compute new track parameters and their covariances using smoother
  AliDebug(1,"Enter UseSmoother");
  
  AliMUONTrackParam *previousTrackParam = (AliMUONTrackParam*) track.GetTrackParamAtHit()->First();
  
  // Smoothed parameters and covariances at first hit = filtered parameters and covariances
  previousTrackParam->SetSmoothParameters(previousTrackParam->GetParameters());
  previousTrackParam->SetSmoothCovariances(previousTrackParam->GetCovariances());
  
  // Compute local chi2 at first hit
  AliMUONHitForRec *hitForRecAtHit = previousTrackParam->GetHitForRecPtr();
  Double_t dX = hitForRecAtHit->GetNonBendingCoor() - previousTrackParam->GetNonBendingCoor();
  Double_t dY = hitForRecAtHit->GetBendingCoor() - previousTrackParam->GetBendingCoor();
  Double_t localChi2 = dX * dX / hitForRecAtHit->GetNonBendingReso2() + dY * dY / hitForRecAtHit->GetBendingReso2();
  
  // Save local chi2 at first hit
  previousTrackParam->SetLocalChi2(localChi2);
  
  AliMUONTrackParam *currentTrackParam = (AliMUONTrackParam*) track.GetTrackParamAtHit()->After(previousTrackParam);
  while (currentTrackParam) {
    
    // Get variables
    const TMatrixD &extrapParameters          = previousTrackParam->GetExtrapParameters();  // X(k+1 k)
    const TMatrixD &filteredParameters        = currentTrackParam->GetParameters();         // X(k k)
    const TMatrixD &previousSmoothParameters  = previousTrackParam->GetSmoothParameters();  // X(k+1 n)
    const TMatrixD &propagator                = previousTrackParam->GetPropagator();        // F(k)
    const TMatrixD &extrapCovariances         = previousTrackParam->GetExtrapCovariances(); // C(k+1 k)
    const TMatrixD &filteredCovariances       = currentTrackParam->GetCovariances();        // C(k k)
    const TMatrixD &previousSmoothCovariances = previousTrackParam->GetSmoothCovariances(); // C(k+1 n)
    
    // Compute smoother gain: A(k) = C(kk) * F(f)^t * (C(k+1 k))^-1
    TMatrixD extrapWeight(extrapCovariances);
    if (extrapWeight.Determinant() != 0) {
      extrapWeight.Invert(); // (C(k+1 k))^-1
    } else {
      AliWarning(" Determinant = 0");
      return kFALSE;
    }
    TMatrixD smootherGain(filteredCovariances,TMatrixD::kMultTranspose,propagator); // C(kk) * F(f)^t
    smootherGain *= extrapWeight; // C(kk) * F(f)^t * (C(k+1 k))^-1
    
    // Compute smoothed parameters: X(k n) = X(k k) + A(k) * (X(k+1 n) - X(k+1 k))
    TMatrixD tmpParam(previousSmoothParameters,TMatrixD::kMinus,extrapParameters); // X(k+1 n) - X(k+1 k)
    TMatrixD smoothParameters(smootherGain,TMatrixD::kMult,tmpParam); // A(k) * (X(k+1 n) - X(k+1 k))
    smoothParameters += filteredParameters; // X(k k) + A(k) * (X(k+1 n) - X(k+1 k))
    
    // Save smoothed parameters
    currentTrackParam->SetSmoothParameters(smoothParameters);
    
    // Compute smoothed covariances: C(k n) = C(k k) + A(k) * (C(k+1 n) - C(k+1 k)) * (A(k))^t
    TMatrixD tmpCov(previousSmoothCovariances,TMatrixD::kMinus,extrapCovariances); // C(k+1 n) - C(k+1 k)
    TMatrixD tmpCov2(tmpCov,TMatrixD::kMultTranspose,smootherGain); // (C(k+1 n) - C(k+1 k)) * (A(k))^t
    TMatrixD smoothCovariances(smootherGain,TMatrixD::kMult,tmpCov2); // A(k) * (C(k+1 n) - C(k+1 k)) * (A(k))^t
    smoothCovariances += filteredCovariances; // C(k k) + A(k) * (C(k+1 n) - C(k+1 k)) * (A(k))^t
    
    // Save smoothed covariances
    currentTrackParam->SetSmoothCovariances(smoothCovariances);
    
    // Compute smoothed residual: r(k n) = hit - X(k n)
    hitForRecAtHit = currentTrackParam->GetHitForRecPtr();
    TMatrixD smoothResidual(2,1);
    smoothResidual.Zero();
    smoothResidual(0,0) = hitForRecAtHit->GetNonBendingCoor() - smoothParameters(0,0);
    smoothResidual(1,0) = hitForRecAtHit->GetBendingCoor() - smoothParameters(2,0);
    
    // Compute weight of smoothed residual: W(k n) = (hitCov - C(k n))^-1
    TMatrixD smoothResidualWeight(2,2);
    smoothResidualWeight(0,0) = hitForRecAtHit->GetNonBendingReso2() - smoothCovariances(0,0);
    smoothResidualWeight(0,1) = - smoothCovariances(0,2);
    smoothResidualWeight(1,0) = - smoothCovariances(2,0);
    smoothResidualWeight(1,1) = hitForRecAtHit->GetBendingReso2() - smoothCovariances(2,2);
    if (smoothResidualWeight.Determinant() != 0) {
      smoothResidualWeight.Invert();
    } else {
      AliWarning(" Determinant = 0");
      return kFALSE;
    }
    
    // Compute local chi2 = (r(k n))^t * W(k n) * r(k n)
    TMatrixD tmpChi2(smoothResidual,TMatrixD::kTransposeMult,smoothResidualWeight); // (r(k n))^t * W(k n)
    TMatrixD localChi2(tmpChi2,TMatrixD::kMult,smoothResidual); // (r(k n))^t * W(k n) * r(k n)
    
    // Save local chi2
    currentTrackParam->SetLocalChi2(localChi2(0,0));
    
    previousTrackParam = currentTrackParam;
    currentTrackParam = (AliMUONTrackParam*) track.GetTrackParamAtHit()->After(previousTrackParam);
  }
  
  return kTRUE;
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::ComplementTracks()
{
  /// Complete tracks by adding missing clusters (if there is an overlap between
  /// two detection elements, the track may have two clusters in the same chamber)
  /// Recompute track parameters and covariances at each clusters
  AliDebug(1,"Enter ComplementTracks");
  
  Int_t chamberId, detElemId;
  Double_t chi2OfHitForRec, addChi2TrackAtHit, bestAddChi2TrackAtHit;
  Double_t maxChi2OfHitForRec = 2. * AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking() *
                                     AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking(); // 2 because 2 quantities in chi2
  Bool_t foundOneHit, trackModified;
  AliMUONHitForRec *hitForRec;
  AliMUONTrackParam *trackParam, *previousTrackParam, *nextTrackParam;
  AliMUONTrackParam trackParamAtHit;
  AliMUONTrackParam bestTrackParamAtHit;
  
  // Remove double track to complete only "good" tracks
  RemoveDoubleTracks();
  
  AliMUONTrack *track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    trackModified = kFALSE;
    
    trackParam = (AliMUONTrackParam*)track->GetTrackParamAtHit()->First();
    previousTrackParam = trackParam;
    while (trackParam) {
      foundOneHit = kFALSE;
      bestAddChi2TrackAtHit = 1.e10;
      
      // prepare nextTrackParam before adding new cluster because of the sorting
      nextTrackParam = (AliMUONTrackParam*)track->GetTrackParamAtHit()->After(trackParam);
      
      chamberId = trackParam->GetHitForRecPtr()->GetChamberNumber();
      detElemId = trackParam->GetHitForRecPtr()->GetDetElemId();
      
      // look for one second candidate in the same chamber
      for (Int_t hit = 0; hit < fNHitsForRecPerChamber[chamberId]; hit++) {
	
	hitForRec = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[chamberId]+hit);
        
	// look for a cluster in another detection element
	if (hitForRec->GetDetElemId() == detElemId) continue;
	
	// try to add the current hit fast
	if (!TryOneHitForRecFast(*trackParam, hitForRec)) continue;
	
	// try to add the current hit accurately
	// never use track parameters at last cluster because the covariance matrix is meaningless
	if (nextTrackParam) chi2OfHitForRec = TryOneHitForRec(*trackParam, hitForRec, trackParamAtHit);
	else chi2OfHitForRec = TryOneHitForRec(*previousTrackParam, hitForRec, trackParamAtHit);
	
	// if good chi2 then consider to add this cluster to the track
	if (chi2OfHitForRec < maxChi2OfHitForRec) {
          
	  // Compute local track parameters including "hitForRec" using kalman filter
          addChi2TrackAtHit = RunKalmanFilter(trackParamAtHit);
          
	  // keep track of the best cluster
	  if (addChi2TrackAtHit < bestAddChi2TrackAtHit) {
	    bestAddChi2TrackAtHit = addChi2TrackAtHit;
	    bestTrackParamAtHit = trackParamAtHit;
	    foundOneHit = kTRUE;
	  }
	  
	}
	
      }
      
      // add new cluster if any
      if (foundOneHit) {
	UpdateTrack(*track,bestTrackParamAtHit,bestAddChi2TrackAtHit);
	bestTrackParamAtHit.SetAloneInChamber(kFALSE);
	trackParam->SetAloneInChamber(kFALSE);
	trackModified = kTRUE;
      }
      
      previousTrackParam = trackParam;
      trackParam = nextTrackParam;
    }
    
    // re-compute track parameters using kalman filter if needed
    if (trackModified) RetraceTrack(*track,kTRUE);
    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
  }
  
}

//__________________________________________________________________________
void AliMUONTrackReconstructorK::ImproveTracks()
{
  /// Improve tracks by removing clusters with local chi2 highter than the defined cut
  /// Recompute track parameters and covariances at the remaining clusters
  AliDebug(1,"Enter ImproveTracks");
  
  Double_t localChi2, worstLocalChi2;
  Int_t worstChamber, previousChamber;
  AliMUONTrack *track, *nextTrack;
  AliMUONTrackParam *trackParamAtHit, *worstTrackParamAtHit, *previousTrackParam, *nextTrackParam;
  Bool_t smoothed;
  Double_t sigmaCut2 = AliMUONReconstructor::GetRecoParam()->GetSigmaCutForImprovement() *
                       AliMUONReconstructor::GetRecoParam()->GetSigmaCutForImprovement();
  
  // Remove double track to improve only "good" tracks
  RemoveDoubleTracks();
  
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    
    // prepare next track in case the actual track is suppressed
    nextTrack = (AliMUONTrack*) fRecTracksPtr->After(track);
    
    while (!track->IsImproved()) {
      
      // Run smoother if required
      smoothed = kFALSE;
      if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) smoothed = RunSmoother(*track);
      
      // Use standard procedure to compute local chi2 if smoother not required or not working
      if (!smoothed) {
        
        // Update track parameters and covariances
        track->UpdateCovTrackParamAtHit();
        
        // Compute local chi2 of each hits
        track->ComputeLocalChi2(kTRUE);
      }
      
      // Look for the hit to remove
      worstTrackParamAtHit = NULL;
      worstLocalChi2 = 0.;
      trackParamAtHit = (AliMUONTrackParam*)track->GetTrackParamAtHit()->First();
      while (trackParamAtHit) {
        
        // Pick up hit with the worst chi2
        localChi2 = trackParamAtHit->GetLocalChi2();
        if (localChi2 > worstLocalChi2) {
          worstLocalChi2 = localChi2;
          worstTrackParamAtHit = trackParamAtHit;
        }
        
	trackParamAtHit = (AliMUONTrackParam*)track->GetTrackParamAtHit()->After(trackParamAtHit);
      }
      
      // Check if bad removable hit found
      if (!worstTrackParamAtHit) {
        track->SetImproved(kTRUE);
        break;
      }
      
      // Check whether the worst chi2 is under requirement or not
      if (worstLocalChi2 < 2. * sigmaCut2) { // 2 because 2 quantities in chi2
        track->SetImproved(kTRUE);
        break;
      }
      
      // if the worst hit is not removable then remove the entire track
      if (!worstTrackParamAtHit->IsRemovable() && worstTrackParamAtHit->IsAloneInChamber()) {
	fRecTracksPtr->Remove(track);
  	fNRecTracks--;
        break;
      }
      
      // Reset the second hit in the same station as being not removable
      // or reset the second hit in the same chamber as being alone
      worstChamber = worstTrackParamAtHit->GetHitForRecPtr()->GetChamberNumber();
      previousTrackParam = (AliMUONTrackParam*) track->GetTrackParamAtHit()->Before(worstTrackParamAtHit);
      nextTrackParam = (AliMUONTrackParam*) track->GetTrackParamAtHit()->After(worstTrackParamAtHit);
      if (worstTrackParamAtHit->IsAloneInChamber()) { // Worst hit removable and alone in chamber
	
	if (worstChamber%2 == 0) { // Modify flags in next chamber
	  
	  nextTrackParam->SetRemovable(kFALSE);
	  if (!nextTrackParam->IsAloneInChamber()) // Make sure both hits in second chamber are not removable anymore
	    ((AliMUONTrackParam*) track->GetTrackParamAtHit()->After(nextTrackParam))->SetRemovable(kFALSE);
	  
	} else { // Modify flags in previous chamber
	  
	  previousTrackParam->SetRemovable(kFALSE);
	  if (!previousTrackParam->IsAloneInChamber()) // Make sure both hits in second chamber are not removable anymore
	    ((AliMUONTrackParam*) track->GetTrackParamAtHit()->Before(previousTrackParam))->SetRemovable(kFALSE);
	  
	}
	
      } else { // Worst hit not alone in its chamber
        
	if (previousTrackParam) previousChamber = previousTrackParam->GetHitForRecPtr()->GetChamberNumber();
	else previousChamber = -1;
	
	if (previousChamber == worstChamber) { // the second hit on the same chamber is the previous one
	  
	  previousTrackParam->SetAloneInChamber(kTRUE);
	  // transfert the removability to the second hit
	  if (worstTrackParamAtHit->IsRemovable()) previousTrackParam->SetRemovable(kTRUE);
	  
	} else { // the second hit on the same chamber is the next one
	  
	  nextTrackParam->SetAloneInChamber(kTRUE);
	  // transfert the removability to the second hit
	  if (worstTrackParamAtHit->IsRemovable()) nextTrackParam->SetRemovable(kTRUE);
	  
	}
	
      }
      
      // Remove the worst hit
      track->RemoveTrackParamAtHit(worstTrackParamAtHit);
      
      // Re-calculate track parameters
      // - from the hit immediately downstream the one suppressed
      // - or from the begining - if parameters have been re-computed using the standard method (kalman parameters have been lost)
      //			- or if the removed hit was the last one
      if (smoothed && nextTrackParam) RetracePartialTrack(*track,nextTrackParam);
      else RetraceTrack(*track,kTRUE);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "ImproveTracks: track " << fRecTracksPtr->IndexOf(track)+1 << " improved " << endl;
      }
      
    }
    
    track = nextTrack;
  }
  
  // compress the array in case of some tracks have been removed
  fRecTracksPtr->Compress();
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::Finalize()
{
  /// Fill AliMUONTrack's fHitForRecAtHit array
  
  AliMUONTrack *track;
  AliMUONTrackParam *trackParamAtHit;
  Bool_t smoothed = kFALSE;
  
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    
    // update track parameters (using smoother if required) if not already done
    if (!track->IsImproved()) {
      smoothed = kFALSE;
      if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) smoothed = RunSmoother(*track);
      if (!smoothed) track->UpdateCovTrackParamAtHit();
    }
    
    trackParamAtHit = (AliMUONTrackParam*) (track->GetTrackParamAtHit()->First());
    while (trackParamAtHit) {
      
      // copy smoothed parameters and covariances if any
      if (smoothed) {
        trackParamAtHit->SetParameters(trackParamAtHit->GetSmoothParameters());
        trackParamAtHit->SetCovariances(trackParamAtHit->GetSmoothCovariances());
      }
      
      // update array of track hits
      track->AddHitForRecAtHit(trackParamAtHit->GetHitForRecPtr());
      
      trackParamAtHit = (AliMUONTrackParam*) (track->GetTrackParamAtHit()->After(trackParamAtHit));
    }
    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
    
  }
    
}

