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
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
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
  /// Constructor
}

  //__________________________________________________________________________
AliMUONTrackReconstructorK::~AliMUONTrackReconstructorK()
{
/// Destructor
} 

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::MakeTrackCandidates(const AliMUONVClusterStore& clusterStore)
{
  /// To make track candidates (assuming linear propagation if the flag fgkMakeTrackCandidatesFast is set to kTRUE):
  /// Start with segments station(1..) 4 or 5 then follow track in station 5 or 4.
  /// Good candidates are made of at least three clusters.
  /// Keep only best candidates or all of them according to the flag fgkTrackAllTracks.
  
  TClonesArray *segments;
  AliMUONTrack *track;
  Int_t iCandidate = 0;
  Bool_t clusterFound;

  AliDebug(1,"Enter MakeTrackCandidates");

  // Loop over stations(1..) 5 and 4 and make track candidates
  for (Int_t istat=4; istat>=3; istat--) {
    
    // Make segments in the station
    segments = MakeSegmentsInStation(clusterStore, istat);
    
    // Loop over segments
    for (Int_t iseg=0; iseg<segments->GetEntriesFast(); iseg++) {
      AliDebug(1,Form("Making primary candidate(1..) %d",++iCandidate));
      
      // Transform segments to tracks and put them at the end of fRecTracksPtr
      track = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack((AliMUONObjectPair*)((*segments)[iseg]));
      fNRecTracks++;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
        cout<<endl<<"Track parameter covariances at first cluster:"<<endl;
        ((AliMUONTrackParam*) (track->GetTrackParamAtCluster()->First()))->GetCovariances().Print();
      }
      
      // Look for compatible cluster(s) in the other station
      if (AliMUONReconstructor::GetRecoParam()->MakeTrackCandidatesFast())
	clusterFound = FollowLinearTrackInStation(*track, clusterStore, 7-istat);
      else {
        // First recompute track parameters and covariances on station(1..) 5 using Kalman filter
        // (to make sure all tracks are treated in the same way)
        if (istat == 4) RetraceTrack(*track, kFALSE);
	clusterFound = FollowTrackInStation(*track, clusterStore, 7-istat);
      }
      
      // Remove track if no cluster found
      if (!clusterFound) {
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
  /// Re-run the kalman filter from the most downstream cluster to the most uptream one
  AliDebug(1,"Enter RetraceTrack");
  
  AliMUONTrackParam* startingTrackParam = (AliMUONTrackParam*) trackCandidate.GetTrackParamAtCluster()->Last();
  
  // Reset the "seed" (= track parameters and their covariances at last cluster) if required
  if (resetSeed) {
    // => Shift track parameters at the position of the last cluster
    AliMUONVCluster* cluster = startingTrackParam->GetClusterPtr();
    startingTrackParam->SetNonBendingCoor(cluster->GetX());
    startingTrackParam->SetBendingCoor(cluster->GetY());
    
    // => Re-compute and reset track parameters covariances at last cluster (as if the other clusters did not exist)
    const TMatrixD& kParamCov = startingTrackParam->GetCovariances();
    TMatrixD newParamCov(5,5);
    newParamCov.Zero();
    // Non bending plane
    newParamCov(0,0) = cluster->GetErrX2();
    newParamCov(1,1) = 100.*kParamCov(1,1);
    // Bending plane
    newParamCov(2,2) = cluster->GetErrY2();
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
  /// Re-run the kalman filter from the cluster attached to startingTrackParam to the most uptream cluster
  AliDebug(1,"Enter RetracePartialTrack");
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "RetracePartialTrack: track chi2 before re-tracking: " << trackCandidate.GetGlobalChi2() << endl;
  }
  
  // Reset the track chi2
  trackCandidate.SetGlobalChi2(startingTrackParam->GetTrackChi2());
  
  // loop over attached clusters until the first one and recompute track parameters and covariances using kalman filter
  Int_t expectedChamber = startingTrackParam->GetClusterPtr()->GetChamberId() - 1;
  Int_t currentChamber;
  Double_t addChi2TrackAtCluster;
  AliMUONTrackParam* trackParamAtCluster = (AliMUONTrackParam*) trackCandidate.GetTrackParamAtCluster()->Before(startingTrackParam); 
  while (trackParamAtCluster) {
    
    // reset track parameters and their covariances
    trackParamAtCluster->SetParameters(startingTrackParam->GetParameters());
    trackParamAtCluster->SetZ(startingTrackParam->GetZ());
    trackParamAtCluster->SetCovariances(startingTrackParam->GetCovariances());
    
    // add MCS effect
    AliMUONTrackExtrap::AddMCSEffect(trackParamAtCluster,AliMUONConstants::ChamberThicknessInX0(),1.);
    
    // reset propagator for smoother
    if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) trackParamAtCluster->ResetPropagator();
    
    // add MCS in missing chambers if any (at most 2 chambers can be missing according to tracking criteria)
    currentChamber = trackParamAtCluster->GetClusterPtr()->GetChamberId();
    while (currentChamber < expectedChamber) {
      // extrapolation to the missing chamber (update the propagator)
      AliMUONTrackExtrap::ExtrapToZCov(trackParamAtCluster, AliMUONConstants::DefaultChamberZ(expectedChamber),
				       AliMUONReconstructor::GetRecoParam()->UseSmoother());
      // add MCS effect
      AliMUONTrackExtrap::AddMCSEffect(trackParamAtCluster,AliMUONConstants::ChamberThicknessInX0(),1.);
      expectedChamber--;
    }
    
    // extrapolation to the plane of the cluster attached to the current trackParamAtCluster (update the propagator)
    AliMUONTrackExtrap::ExtrapToZCov(trackParamAtCluster, trackParamAtCluster->GetClusterPtr()->GetZ(),
				     AliMUONReconstructor::GetRecoParam()->UseSmoother());
    
    if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) {
      // save extrapolated parameters for smoother
      trackParamAtCluster->SetExtrapParameters(trackParamAtCluster->GetParameters());
      
      // save extrapolated covariance matrix for smoother
      trackParamAtCluster->SetExtrapCovariances(trackParamAtCluster->GetCovariances());
    }
    
    // Compute new track parameters including "clusterCh2" using kalman filter
    addChi2TrackAtCluster = RunKalmanFilter(*trackParamAtCluster);
    
    // Update the track chi2
    trackCandidate.SetGlobalChi2(trackCandidate.GetGlobalChi2() + addChi2TrackAtCluster);
    trackParamAtCluster->SetTrackChi2(trackCandidate.GetGlobalChi2());
    
    // prepare next step
    expectedChamber = currentChamber - 1;
    startingTrackParam = trackParamAtCluster;
    trackParamAtCluster = (AliMUONTrackParam*) (trackCandidate.GetTrackParamAtCluster()->Before(startingTrackParam)); 
  }
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "RetracePartialTrack: track chi2 after re-tracking: " << trackCandidate.GetGlobalChi2() << endl;
  }
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::FollowTracks(const AliMUONVClusterStore& clusterStore)
{
  /// Follow tracks in stations(1..) 3, 2 and 1
  AliDebug(1,"Enter FollowTracks");
  
  AliMUONTrack *track;
  Int_t currentNRecTracks;
  Bool_t clusterFound;
  
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
        cout<<endl<<"Track parameter covariances at first cluster:"<<endl;
        ((AliMUONTrackParam*) (track->GetTrackParamAtCluster()->First()))->GetCovariances().Print();
      }
      
      // Look for compatible cluster(s) in station(0..) "station"
      clusterFound = FollowTrackInStation(*track, clusterStore, station);
      
      // Try to recover track if required
      if (!clusterFound && AliMUONReconstructor::GetRecoParam()->RecoverTracks())
	clusterFound = RecoverTrack(*track, clusterStore, station);
      
      // remove track if no cluster found
      if (!clusterFound) {
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
Bool_t AliMUONTrackReconstructorK::FollowTrackInStation(AliMUONTrack &trackCandidate, const AliMUONVClusterStore& clusterStore, Int_t nextStation)
{
  /// Follow trackCandidate in station(0..) nextStation and search for compatible cluster(s)
  /// Keep all possibilities or only the best one(s) according to the flag fgkTrackAllTracks:
  /// kTRUE:  duplicate "trackCandidate" if there are several possibilities and add the new tracks at the end of
  ///         fRecTracksPtr to avoid conficts with other track candidates at this current stage of the tracking procedure.
  ///         Remove the obsolete "trackCandidate" at the end.
  /// kFALSE: add only the best cluster(s) to the "trackCandidate". Try to add a couple of clusters in priority.
  /// return kTRUE if new cluster(s) have been found (otherwise return kFALSE)
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
  
  Double_t chi2OfCluster;
  Double_t maxChi2OfCluster = 2. * AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking() *
				   AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking(); // 2 because 2 quantities in chi2
  Double_t addChi2TrackAtCluster1;
  Double_t addChi2TrackAtCluster2;
  Double_t bestAddChi2TrackAtCluster1 = 1.e10;
  Double_t bestAddChi2TrackAtCluster2 = 1.e10;
  Bool_t foundOneCluster = kFALSE;
  Bool_t foundTwoClusters = kFALSE;
  AliMUONTrack *newTrack = 0x0;
  AliMUONVCluster *clusterCh1, *clusterCh2;
  AliMUONTrackParam extrapTrackParam;
  AliMUONTrackParam extrapTrackParamAtCluster1;
  AliMUONTrackParam extrapTrackParamAtCluster2;
  AliMUONTrackParam bestTrackParamAtCluster1;
  AliMUONTrackParam bestTrackParamAtCluster2;
  
  Int_t nClusters = clusterStore.GetSize();
  Bool_t *clusterCh1Used = new Bool_t[nClusters];
  for (Int_t i = 0; i < nClusters; i++) clusterCh1Used[i] = kFALSE;
  Int_t iCluster1;
  
  // Get track parameters
  AliMUONTrackParam extrapTrackParamAtCh(*(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->First());
  
  // Add MCS effect
  AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(),1.);
  
  // reset propagator for smoother
  if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) extrapTrackParamAtCh.ResetPropagator();
  
  // Add MCS in the missing chamber if any (only 1 chamber can be missing according to tracking criteria)
  if (ch1 < ch2 && extrapTrackParamAtCh.GetClusterPtr()->GetChamberId() > ch2 + 1) {
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
    cout<<endl<<"Track parameter covariances at first cluster extrapolated to z = "<<AliMUONConstants::DefaultChamberZ(ch2)<<":"<<endl;
    extrapTrackParamAtCh.GetCovariances().Print();
  }
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "FollowTrackInStation: look for clusters in chamber(1..): " << ch2+1 << endl;
  }
  
  // Create iterators to loop over clusters in both chambers
  TIter nextInCh1(clusterStore.CreateChamberIterator(ch1,ch1));
  TIter nextInCh2(clusterStore.CreateChamberIterator(ch2,ch2));
  
  // look for candidates in chamber 2
  while ( ( clusterCh2 = static_cast<AliMUONVCluster*>(nextInCh2()) ) ) {
    
    // try to add the current cluster fast
    if (!TryOneClusterFast(extrapTrackParamAtCh, clusterCh2)) continue;
    
    // try to add the current cluster accuratly
    chi2OfCluster = TryOneCluster(extrapTrackParamAtCh, clusterCh2, extrapTrackParamAtCluster2,
				  AliMUONReconstructor::GetRecoParam()->UseSmoother());
    
    // if good chi2 then try to attach a cluster in the other chamber too
    if (chi2OfCluster < maxChi2OfCluster) {
      Bool_t foundSecondCluster = kFALSE;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: found one cluster in chamber(1..): " << ch2+1
	     << " (Chi2 = " << chi2OfCluster << ")" << endl;
        cout << "                      look for second clusters in chamber(1..): " << ch1+1 << " ..." << endl;
      }
      
      if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) {
        // save extrapolated parameters for smoother
        extrapTrackParamAtCluster2.SetExtrapParameters(extrapTrackParamAtCluster2.GetParameters());
        
        // save extrapolated covariance matrix for smoother
        extrapTrackParamAtCluster2.SetExtrapCovariances(extrapTrackParamAtCluster2.GetCovariances());
      }
      
      // Compute new track parameters including "clusterCh2" using kalman filter
      addChi2TrackAtCluster2 = RunKalmanFilter(extrapTrackParamAtCluster2);
      
      // copy new track parameters for next step
      extrapTrackParam = extrapTrackParamAtCluster2;
      
      // add MCS effect
      AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParam,AliMUONConstants::ChamberThicknessInX0(),1.);
      
      // reset propagator for smoother
      if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) extrapTrackParam.ResetPropagator();
      
      //Extrapolate track parameters to chamber "ch1"
      AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParam, AliMUONConstants::DefaultChamberZ(ch1),
				       AliMUONReconstructor::GetRecoParam()->UseSmoother());
      
      // reset cluster iterator of chamber 1
      nextInCh1.Reset();
      iCluster1 = -1;
      
      // look for second candidates in chamber 1
      while ( ( clusterCh1 = static_cast<AliMUONVCluster*>(nextInCh1()) ) ) {
	iCluster1++;
	
    	// try to add the current cluster fast
    	if (!TryOneClusterFast(extrapTrackParam, clusterCh1)) continue;
    	
    	// try to add the current cluster accuratly
    	chi2OfCluster = TryOneCluster(extrapTrackParam, clusterCh1, extrapTrackParamAtCluster1,
					  AliMUONReconstructor::GetRecoParam()->UseSmoother());
    	
	// if good chi2 then consider to add the 2 clusters to the "trackCandidate"
	if (chi2OfCluster < maxChi2OfCluster) {
	  foundSecondCluster = kTRUE;
	  foundTwoClusters = kTRUE;
          
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowTrackInStation: found one cluster in chamber(1..): " << ch1+1
	  	 << " (Chi2 = " << chi2OfCluster << ")" << endl;
	  }
          
          if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) {
            // save extrapolated parameters for smoother
            extrapTrackParamAtCluster1.SetExtrapParameters(extrapTrackParamAtCluster1.GetParameters());
            
            // save extrapolated covariance matrix for smoother
            extrapTrackParamAtCluster1.SetExtrapCovariances(extrapTrackParamAtCluster1.GetCovariances());
          }
          
          // Compute new track parameters including "clusterCh1" using kalman filter
          addChi2TrackAtCluster1 = RunKalmanFilter(extrapTrackParamAtCluster1);
          
	  if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
	    // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new clusters
            newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	    UpdateTrack(*newTrack,extrapTrackParamAtCluster1,extrapTrackParamAtCluster2,addChi2TrackAtCluster1,addChi2TrackAtCluster2);
	    fNRecTracks++;
	    
	    // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
	    // (going in the right direction)
	    if (nextStation == 4) RetraceTrack(*newTrack,kTRUE);
	    
	    // Tag clusterCh1 as used
	    clusterCh1Used[iCluster1] = kTRUE;
	    
	    // Printout for debuging
	    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	      cout << "FollowTrackInStation: added two clusters in station(1..): " << nextStation+1 << endl;
	      if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	    }
	    
          } else if (addChi2TrackAtCluster1+addChi2TrackAtCluster2 < bestAddChi2TrackAtCluster1+bestAddChi2TrackAtCluster2) {
	    // keep track of the best couple of clusters
	    bestAddChi2TrackAtCluster1 = addChi2TrackAtCluster1;
	    bestAddChi2TrackAtCluster2 = addChi2TrackAtCluster2;
	    bestTrackParamAtCluster1 = extrapTrackParamAtCluster1;
	    bestTrackParamAtCluster2 = extrapTrackParamAtCluster2;
          }
	  
	}
	
      }
      
      // if no clusterCh1 found then consider to add clusterCh2 only
      if (!foundSecondCluster) {
	foundOneCluster = kTRUE;
        
	if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
          newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  UpdateTrack(*newTrack,extrapTrackParamAtCluster2,addChi2TrackAtCluster2);
	  fNRecTracks++;
	  
	  // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
	  // (going in the right direction)
	  if (nextStation == 4) RetraceTrack(*newTrack,kTRUE);
	  
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowTrackInStation: added one cluster in chamber(1..): " << ch2+1 << endl;
	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	  }
	  
	} else if (!foundTwoClusters && addChi2TrackAtCluster2 < bestAddChi2TrackAtCluster1) {
	  // keep track of the best single cluster except if a couple of clusters has already been found
	  bestAddChi2TrackAtCluster1 = addChi2TrackAtCluster2;
	  bestTrackParamAtCluster1 = extrapTrackParamAtCluster2;
        }
	
      }
      
    }
    
  }
  
  // look for candidates in chamber 1 not already attached to a track
  // if we want to keep all possible tracks or if no good couple of clusters has been found
  if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks() || !foundTwoClusters) {
    
    // Printout for debuging
    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
      cout << "FollowTrackInStation: look for single clusters in chamber(1..): " << ch1+1 << endl;
    }
    
    // add MCS effect for next step
    AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(),1.);
    
    //Extrapolate trackCandidate to chamber "ch1"
    AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(ch1),
				     AliMUONReconstructor::GetRecoParam()->UseSmoother());
    
    // reset cluster iterator of chamber 1
    nextInCh1.Reset();
    iCluster1 = -1;
    
    // look for second candidates in chamber 1
    while ( ( clusterCh1 = static_cast<AliMUONVCluster*>(nextInCh1()) ) ) {
      iCluster1++;
      
      if (clusterCh1Used[iCluster1]) continue; // Skip clusters already used
      
      // try to add the current cluster fast
      if (!TryOneClusterFast(extrapTrackParamAtCh, clusterCh1)) continue;
      
      // try to add the current cluster accuratly
      chi2OfCluster = TryOneCluster(extrapTrackParamAtCh, clusterCh1, extrapTrackParamAtCluster1,
				    AliMUONReconstructor::GetRecoParam()->UseSmoother());
      
      // if good chi2 then consider to add clusterCh1
      // We do not try to attach a cluster in the other chamber too since it has already been done above
      if (chi2OfCluster < maxChi2OfCluster) {
	foundOneCluster = kTRUE;
  	
	// Printout for debuging
  	if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	  cout << "FollowTrackInStation: found one cluster in chamber(1..): " << ch1+1
  	       << " (Chi2 = " << chi2OfCluster << ")" << endl;
  	}
        
	if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) {
          // save extrapolated parameters for smoother
          extrapTrackParamAtCluster1.SetExtrapParameters(extrapTrackParamAtCluster1.GetParameters());
          
          // save extrapolated covariance matrix for smoother
          extrapTrackParamAtCluster1.SetExtrapCovariances(extrapTrackParamAtCluster1.GetCovariances());
        }
        
        // Compute new track parameters including "clusterCh1" using kalman filter
        addChi2TrackAtCluster1 = RunKalmanFilter(extrapTrackParamAtCluster1);
        
	if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
  	  newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  UpdateTrack(*newTrack,extrapTrackParamAtCluster1,addChi2TrackAtCluster1);
	  fNRecTracks++;
  	  
	  // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
	  // (going in the right direction)
	  if (nextStation == 4) RetraceTrack(*newTrack,kTRUE);
	  
	  // Printout for debuging
  	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	    cout << "FollowTrackInStation: added one cluster in chamber(1..): " << ch1+1 << endl;
  	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
  	  }
	  
  	} else if (addChi2TrackAtCluster1 < bestAddChi2TrackAtCluster1) {
	  // keep track of the best single cluster except if a couple of clusters has already been found
	  bestAddChi2TrackAtCluster1 = addChi2TrackAtCluster1;
	  bestTrackParamAtCluster1 = extrapTrackParamAtCluster1;
  	}
	
      }
      
    }
    
  }
  
  // fill out the best track if required else clean up the fRecTracksPtr array
  if (!AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
    if (foundTwoClusters) {
      UpdateTrack(trackCandidate,bestTrackParamAtCluster1,bestTrackParamAtCluster2,bestAddChi2TrackAtCluster1,bestAddChi2TrackAtCluster2);
      
      // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
      // (going in the right direction)
      if (nextStation == 4) RetraceTrack(trackCandidate,kTRUE);

      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the two best clusters in station(1..): " << nextStation+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else if (foundOneCluster) {
      UpdateTrack(trackCandidate,bestTrackParamAtCluster1,bestAddChi2TrackAtCluster1);
      
      // if we are arrived on station(1..) 5, recompute track parameters and covariances starting from this station
      // (going in the right direction)
      if (nextStation == 4) RetraceTrack(trackCandidate,kTRUE);

      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructorK") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the best cluster in chamber(1..): " << bestTrackParamAtCluster1.GetClusterPtr()->GetChamberId()+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else {
      delete [] clusterCh1Used;
      return kFALSE;
    }
  } else if (foundOneCluster || foundTwoClusters) {
    
    // remove obsolete track
    fRecTracksPtr->Remove(&trackCandidate);
    fNRecTracks--;
    
  } else {
    delete [] clusterCh1Used;
    return kFALSE;
  }  
  delete [] clusterCh1Used;
  return kTRUE;
  
}

  //__________________________________________________________________________
Double_t AliMUONTrackReconstructorK::RunKalmanFilter(AliMUONTrackParam &trackParamAtCluster)
{
  /// Compute new track parameters and their covariances including new cluster using kalman filter
  /// return the additional track chi2
  AliDebug(1,"Enter RunKalmanFilter");
  
  // Get actual track parameters (p)
  TMatrixD param(trackParamAtCluster.GetParameters());
  
  // Get new cluster parameters (m)
  AliMUONVCluster *cluster = trackParamAtCluster.GetClusterPtr();
  TMatrixD clusterParam(5,1);
  clusterParam.Zero();
  clusterParam(0,0) = cluster->GetX();
  clusterParam(2,0) = cluster->GetY();
  
  // Compute the actual parameter weight (W)
  TMatrixD paramWeight(trackParamAtCluster.GetCovariances());
  if (paramWeight.Determinant() != 0) {
    paramWeight.Invert();
  } else {
    AliWarning(" Determinant = 0");
    return 1.e10;
  }
  
  // Compute the new cluster weight (U)
  TMatrixD clusterWeight(5,5);
  clusterWeight.Zero();
  clusterWeight(0,0) = 1. / cluster->GetErrX2();
  clusterWeight(2,2) = 1. / cluster->GetErrY2();

  // Compute the new parameters covariance matrix ( (W+U)^-1 )
  TMatrixD newParamCov(paramWeight,TMatrixD::kPlus,clusterWeight);
  if (newParamCov.Determinant() != 0) {
    newParamCov.Invert();
  } else {
    AliWarning(" Determinant = 0");
    return 1.e10;
  }
  
  // Save the new parameters covariance matrix
  trackParamAtCluster.SetCovariances(newParamCov);
  
  // Compute the new parameters (p' = ((W+U)^-1)U(m-p) + p)
  TMatrixD tmp(clusterParam,TMatrixD::kMinus,param);
  TMatrixD tmp2(clusterWeight,TMatrixD::kMult,tmp); // U(m-p)
  TMatrixD newParam(newParamCov,TMatrixD::kMult,tmp2); // ((W+U)^-1)U(m-p)
  newParam += param; // ((W+U)^-1)U(m-p) + p
  
  // Save the new parameters
  trackParamAtCluster.SetParameters(newParam);
  
  // Compute the additional chi2 (= ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m))
  tmp = newParam; // p'
  tmp -= param; // (p'-p)
  TMatrixD tmp3(paramWeight,TMatrixD::kMult,tmp); // W(p'-p)
  TMatrixD addChi2Track(tmp,TMatrixD::kTransposeMult,tmp3); // ((p'-p)^-1)W(p'-p)
  tmp = newParam; // p'
  tmp -= clusterParam; // (p'-m)
  TMatrixD tmp4(clusterWeight,TMatrixD::kMult,tmp); // U(p'-m)
  addChi2Track += TMatrixD(tmp,TMatrixD::kTransposeMult,tmp4); // ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m)
  
  return addChi2Track(0,0);
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtCluster, Double_t addChi2)
{
  /// Add 1 cluster to the track candidate
  /// Update chi2 of the track 
  
  // Flag cluster as being not removable
  trackParamAtCluster.SetRemovable(kFALSE);
  trackParamAtCluster.SetLocalChi2(0.); // --> Local chi2 not used
  
  // Update the track chi2 into trackParamAtCluster
  trackParamAtCluster.SetTrackChi2(track.GetGlobalChi2() + addChi2);
  
  // Update the chi2 of the new track
  track.SetGlobalChi2(trackParamAtCluster.GetTrackChi2());
  
  // Update array of TrackParamAtCluster
  track.AddTrackParamAtCluster(trackParamAtCluster,*(trackParamAtCluster.GetClusterPtr()));
  track.GetTrackParamAtCluster()->Sort();
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtCluster1, AliMUONTrackParam &trackParamAtCluster2,
					     Double_t addChi2AtCluster1, Double_t addChi2AtCluster2)
{
  /// Add 2 clusters to the track candidate (order is important)
  /// Update track and local chi2
  
  // Update local chi2 at first cluster
  AliMUONVCluster* cluster1 = trackParamAtCluster1.GetClusterPtr();
  Double_t deltaX = trackParamAtCluster1.GetNonBendingCoor() - cluster1->GetX();
  Double_t deltaY = trackParamAtCluster1.GetBendingCoor() - cluster1->GetY();
  Double_t localChi2AtCluster1 = deltaX*deltaX / cluster1->GetErrX2() +
  			         deltaY*deltaY / cluster1->GetErrY2();
  trackParamAtCluster1.SetLocalChi2(localChi2AtCluster1);
  
  // Flag first cluster as being removable
  trackParamAtCluster1.SetRemovable(kTRUE);
  
  // Update local chi2 at second cluster
  AliMUONVCluster* cluster2 = trackParamAtCluster2.GetClusterPtr();
  AliMUONTrackParam extrapTrackParamAtCluster2(trackParamAtCluster1);
  AliMUONTrackExtrap::ExtrapToZ(&extrapTrackParamAtCluster2, trackParamAtCluster2.GetZ());
  deltaX = extrapTrackParamAtCluster2.GetNonBendingCoor() - cluster2->GetX();
  deltaY = extrapTrackParamAtCluster2.GetBendingCoor() - cluster2->GetY();
  Double_t localChi2AtCluster2 = deltaX*deltaX / cluster2->GetErrX2() +
  			         deltaY*deltaY / cluster2->GetErrY2();
  trackParamAtCluster2.SetLocalChi2(localChi2AtCluster2);
  
  // Flag second cluster as being removable
  trackParamAtCluster2.SetRemovable(kTRUE);
  
  // Update the track chi2 into trackParamAtCluster2
  trackParamAtCluster2.SetTrackChi2(track.GetGlobalChi2() + addChi2AtCluster2);
  
  // Update the track chi2 into trackParamAtCluster1
  trackParamAtCluster1.SetTrackChi2(trackParamAtCluster2.GetTrackChi2() + addChi2AtCluster1);
  
  // Update the chi2 of the new track
  track.SetGlobalChi2(trackParamAtCluster1.GetTrackChi2());
  
  // Update array of trackParamAtCluster
  track.AddTrackParamAtCluster(trackParamAtCluster1,*cluster1);
  track.AddTrackParamAtCluster(trackParamAtCluster2,*cluster2);
  track.GetTrackParamAtCluster()->Sort();
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructorK::RecoverTrack(AliMUONTrack &trackCandidate, const AliMUONVClusterStore& clusterStore, Int_t nextStation)
{
  /// Try to recover the track candidate in the next station
  /// by removing the worst of the two clusters attached in the current station
  /// Return kTRUE if recovering succeeds
  AliDebug(1,"Enter RecoverTrack");
  
  // Do not try to recover track until we have attached cluster(s) on station(1..) 3
  if (nextStation > 1) return kFALSE;
  
  Int_t worstClusterNumber = -1;
  Double_t localChi2, worstLocalChi2 = 0.;
  
  // Look for the cluster to remove
  for (Int_t clusterNumber = 0; clusterNumber < 2; clusterNumber++) {
    AliMUONTrackParam *trackParamAtCluster = (AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->UncheckedAt(clusterNumber);
    
    // check if current cluster is removable
    if (!trackParamAtCluster->IsRemovable()) return kFALSE;
    
    // Pick up cluster with the worst chi2
    localChi2 = trackParamAtCluster->GetLocalChi2();
    if (localChi2 > worstLocalChi2) {
      worstLocalChi2 = localChi2;
      worstClusterNumber = clusterNumber;
    }
  }
  
  // Reset best cluster as being NOT removable
  ((AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->UncheckedAt((worstClusterNumber+1)%2))->SetRemovable(kFALSE);
  
  // Remove the worst cluster
  trackCandidate.RemoveTrackParamAtCluster((AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->UncheckedAt(worstClusterNumber));
  
  // Re-calculate track parameters at the (new) first cluster
  RetracePartialTrack(trackCandidate,(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->UncheckedAt(1));
  
  // Look for new cluster(s) in next station
  return FollowTrackInStation(trackCandidate, clusterStore, nextStation);
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructorK::RunSmoother(AliMUONTrack &track)
{
  /// Compute new track parameters and their covariances using smoother
  AliDebug(1,"Enter UseSmoother");
  
  AliMUONTrackParam *previousTrackParam = (AliMUONTrackParam*) track.GetTrackParamAtCluster()->First();
  
  // Smoothed parameters and covariances at first cluster = filtered parameters and covariances
  previousTrackParam->SetSmoothParameters(previousTrackParam->GetParameters());
  previousTrackParam->SetSmoothCovariances(previousTrackParam->GetCovariances());
  
  // Compute local chi2 at first cluster
  AliMUONVCluster *cluster = previousTrackParam->GetClusterPtr();
  Double_t dX = cluster->GetX() - previousTrackParam->GetNonBendingCoor();
  Double_t dY = cluster->GetY() - previousTrackParam->GetBendingCoor();
  Double_t localChi2 = dX * dX / cluster->GetErrX2() + dY * dY / cluster->GetErrY2();
  
  // Save local chi2 at first cluster
  previousTrackParam->SetLocalChi2(localChi2);
  
  AliMUONTrackParam *currentTrackParam = (AliMUONTrackParam*) track.GetTrackParamAtCluster()->After(previousTrackParam);
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
    
    // Compute smoothed residual: r(k n) = cluster - X(k n)
    cluster = currentTrackParam->GetClusterPtr();
    TMatrixD smoothResidual(2,1);
    smoothResidual.Zero();
    smoothResidual(0,0) = cluster->GetX() - smoothParameters(0,0);
    smoothResidual(1,0) = cluster->GetY() - smoothParameters(2,0);
    
    // Compute weight of smoothed residual: W(k n) = (clusterCov - C(k n))^-1
    TMatrixD smoothResidualWeight(2,2);
    smoothResidualWeight(0,0) = cluster->GetErrX2() - smoothCovariances(0,0);
    smoothResidualWeight(0,1) = - smoothCovariances(0,2);
    smoothResidualWeight(1,0) = - smoothCovariances(2,0);
    smoothResidualWeight(1,1) = cluster->GetErrY2() - smoothCovariances(2,2);
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
    currentTrackParam = (AliMUONTrackParam*) track.GetTrackParamAtCluster()->After(previousTrackParam);
  }
  
  return kTRUE;
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::ComplementTracks(const AliMUONVClusterStore& clusterStore)
{
  /// Complete tracks by adding missing clusters (if there is an overlap between
  /// two detection elements, the track may have two clusters in the same chamber)
  /// Recompute track parameters and covariances at each clusters
  AliDebug(1,"Enter ComplementTracks");
  
  Int_t chamberId, detElemId;
  Double_t chi2OfCluster, addChi2TrackAtCluster, bestAddChi2TrackAtCluster;
  Double_t maxChi2OfCluster = 2. * AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking() *
				   AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking(); // 2 because 2 quantities in chi2
  Bool_t foundOneCluster, trackModified;
  AliMUONVCluster *cluster;
  AliMUONTrackParam *trackParam, *previousTrackParam, *nextTrackParam, trackParamAtCluster, bestTrackParamAtCluster;
  
  // Remove double track to complete only "good" tracks
  RemoveDoubleTracks();
  
  AliMUONTrack *track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    trackModified = kFALSE;
    
    trackParam = (AliMUONTrackParam*)track->GetTrackParamAtCluster()->First();
    previousTrackParam = trackParam;
    while (trackParam) {
      foundOneCluster = kFALSE;
      bestAddChi2TrackAtCluster = 1.e10;
      chamberId = trackParam->GetClusterPtr()->GetChamberId();
      detElemId = trackParam->GetClusterPtr()->GetDetElemId();
      
      // prepare nextTrackParam before adding new cluster because of the sorting
      nextTrackParam = (AliMUONTrackParam*)track->GetTrackParamAtCluster()->After(trackParam);
      
      // Create iterators to loop over clusters in current chamber
      TIter nextInCh(clusterStore.CreateChamberIterator(chamberId,chamberId));
      
      // look for one second candidate in the same chamber
      while ( ( cluster = static_cast<AliMUONVCluster*>(nextInCh()) ) ) {
        
	// look for a cluster in another detection element
	if (cluster->GetDetElemId() == detElemId) continue;
	
	// try to add the current cluster fast
	if (!TryOneClusterFast(*trackParam, cluster)) continue;
	
	// try to add the current cluster accurately
	// never use track parameters at last cluster because the covariance matrix is meaningless
	if (nextTrackParam) chi2OfCluster = TryOneCluster(*trackParam, cluster, trackParamAtCluster);
	else chi2OfCluster = TryOneCluster(*previousTrackParam, cluster, trackParamAtCluster);
	
	// if good chi2 then consider to add this cluster to the track
	if (chi2OfCluster < maxChi2OfCluster) {
          
	  // Compute local track parameters including current cluster using kalman filter
          addChi2TrackAtCluster = RunKalmanFilter(trackParamAtCluster);
          
	  // keep track of the best cluster
	  if (addChi2TrackAtCluster < bestAddChi2TrackAtCluster) {
	    bestAddChi2TrackAtCluster = addChi2TrackAtCluster;
	    bestTrackParamAtCluster = trackParamAtCluster;
	    foundOneCluster = kTRUE;
	  }
	  
	}
	
      }
      
      // add new cluster if any
      if (foundOneCluster) {
	UpdateTrack(*track,bestTrackParamAtCluster,bestAddChi2TrackAtCluster);
	bestTrackParamAtCluster.SetAloneInChamber(kFALSE);
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
  AliMUONTrackParam *trackParamAtCluster, *worstTrackParamAtCluster, *previousTrackParam, *nextTrackParam;
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
        track->UpdateCovTrackParamAtCluster();
        
        // Compute local chi2 of each clusters
        track->ComputeLocalChi2(kTRUE);
      }
      
      // Look for the cluster to remove
      worstTrackParamAtCluster = NULL;
      worstLocalChi2 = 0.;
      trackParamAtCluster = (AliMUONTrackParam*)track->GetTrackParamAtCluster()->First();
      while (trackParamAtCluster) {
        
        // Pick up cluster with the worst chi2
        localChi2 = trackParamAtCluster->GetLocalChi2();
        if (localChi2 > worstLocalChi2) {
          worstLocalChi2 = localChi2;
          worstTrackParamAtCluster = trackParamAtCluster;
        }
        
	trackParamAtCluster = (AliMUONTrackParam*)track->GetTrackParamAtCluster()->After(trackParamAtCluster);
      }
      
      // Check if bad removable cluster found
      if (!worstTrackParamAtCluster) {
        track->SetImproved(kTRUE);
        break;
      }
      
      // Check whether the worst chi2 is under requirement or not
      if (worstLocalChi2 < 2. * sigmaCut2) { // 2 because 2 quantities in chi2
        track->SetImproved(kTRUE);
        break;
      }
      
      // if the worst cluster is not removable then remove the entire track
      if (!worstTrackParamAtCluster->IsRemovable() && worstTrackParamAtCluster->IsAloneInChamber()) {
	fRecTracksPtr->Remove(track);
  	fNRecTracks--;
        break;
      }
      
      // Reset the second cluster in the same station as being not removable
      // or reset the second cluster in the same chamber as being alone
      worstChamber = worstTrackParamAtCluster->GetClusterPtr()->GetChamberId();
      previousTrackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->Before(worstTrackParamAtCluster);
      nextTrackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->After(worstTrackParamAtCluster);
      if (worstTrackParamAtCluster->IsAloneInChamber()) { // Worst cluster removable and alone in chamber
	
	if (worstChamber%2 == 0) { // Modify flags in next chamber
	  
	  nextTrackParam->SetRemovable(kFALSE);
	  if (!nextTrackParam->IsAloneInChamber()) // Make sure both clusters in second chamber are not removable anymore
	    ((AliMUONTrackParam*) track->GetTrackParamAtCluster()->After(nextTrackParam))->SetRemovable(kFALSE);
	  
	} else { // Modify flags in previous chamber
	  
	  previousTrackParam->SetRemovable(kFALSE);
	  if (!previousTrackParam->IsAloneInChamber()) // Make sure both clusters in second chamber are not removable anymore
	    ((AliMUONTrackParam*) track->GetTrackParamAtCluster()->Before(previousTrackParam))->SetRemovable(kFALSE);
	  
	}
	
      } else { // Worst cluster not alone in its chamber
        
	if (previousTrackParam) previousChamber = previousTrackParam->GetClusterPtr()->GetChamberId();
	else previousChamber = -1;
	
	if (previousChamber == worstChamber) { // the second cluster on the same chamber is the previous one
	  
	  previousTrackParam->SetAloneInChamber(kTRUE);
	  // transfert the removability to the second cluster
	  if (worstTrackParamAtCluster->IsRemovable()) previousTrackParam->SetRemovable(kTRUE);
	  
	} else { // the second cluster on the same chamber is the next one
	  
	  nextTrackParam->SetAloneInChamber(kTRUE);
	  // transfert the removability to the second cluster
	  if (worstTrackParamAtCluster->IsRemovable()) nextTrackParam->SetRemovable(kTRUE);
	  
	}
	
      }
      
      // Remove the worst cluster
      track->RemoveTrackParamAtCluster(worstTrackParamAtCluster);
      
      // Re-calculate track parameters
      // - from the cluster immediately downstream the one suppressed
      // - or from the begining - if parameters have been re-computed using the standard method (kalman parameters have been lost)
      //			- or if the removed cluster was the last one
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
  /// Update track parameters and covariances at each attached cluster
  
  AliMUONTrack *track;
  AliMUONTrackParam *trackParamAtCluster;
  Bool_t smoothed = kFALSE;
  
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    
    // update track parameters (using smoother if required) if not already done
    if (!track->IsImproved()) {
      smoothed = kFALSE;
      if (AliMUONReconstructor::GetRecoParam()->UseSmoother()) smoothed = RunSmoother(*track);
      if (!smoothed) track->UpdateCovTrackParamAtCluster();
    }
    
    // copy smoothed parameters and covariances if any
    if (smoothed) {
      
      trackParamAtCluster = (AliMUONTrackParam*) (track->GetTrackParamAtCluster()->First());
      while (trackParamAtCluster) {
	
	trackParamAtCluster->SetParameters(trackParamAtCluster->GetSmoothParameters());
	trackParamAtCluster->SetCovariances(trackParamAtCluster->GetSmoothCovariances());
	
	trackParamAtCluster = (AliMUONTrackParam*) (track->GetTrackParamAtCluster()->After(trackParamAtCluster));
      }
      
    }
    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
  }
    
}

