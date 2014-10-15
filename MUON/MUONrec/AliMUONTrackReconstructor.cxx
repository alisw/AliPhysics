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
/// \class AliMUONTrackReconstructor
/// MUON track reconstructor using the original method
///
/// This class contains as data:
/// - the parameters for the track reconstruction
///
/// It contains as methods, among others:
/// - MakeTracks to build the tracks
//-----------------------------------------------------------------------------

#include "AliMUONTrackReconstructor.h"

#include "AliMUONConstants.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterServer.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONRecoParam.h"
#include "AliMUONGeometryTransformer.h"

#include "AliMpArea.h"

#include "AliLog.h"

#include <TMinuit.h>
#include <Riostream.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TClonesArray.h>

// Functions to be minimized with Minuit
void TrackChi2(Int_t &nParam, Double_t *gradient, Double_t &chi2, Double_t *param, Int_t flag);

using std::endl;
using std::cout;
/// \cond CLASSIMP
ClassImp(AliMUONTrackReconstructor) // Class implementation in ROOT context
/// \endcond

  //__________________________________________________________________________
AliMUONTrackReconstructor::AliMUONTrackReconstructor(const AliMUONRecoParam* recoParam, AliMUONVClusterServer* clusterServer,
						     const AliMUONGeometryTransformer* transformer)
  : AliMUONVTrackReconstructor(recoParam, clusterServer, transformer)
{
  /// Constructor
}

  //__________________________________________________________________________
AliMUONTrackReconstructor::~AliMUONTrackReconstructor()
{
/// Destructor
} 

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::MakeTrackCandidates(AliMUONVClusterStore& clusterStore)
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

  // Unless we're doing combined tracking, we'll clusterize all stations at once
  Int_t firstChamber(0);
  Int_t lastChamber(9);
  
  if (GetRecoParam()->CombineClusterTrackReco()) {
    // ... Here's the exception : ask the clustering to reconstruct
    // clusters *only* in station 4 and 5 for combined tracking
    firstChamber = 6;
  }

  for (Int_t i = firstChamber; i <= lastChamber; ++i ) 
  {
    if (fClusterServer && GetRecoParam()->UseChamber(i)) fClusterServer->Clusterize(i, clusterStore, AliMpArea(), GetRecoParam());
  }
  
  // Loop over stations(1..) 5 and 4 and make track candidates
  for (Int_t istat=4; istat>=3; istat--) {
    
    // Make segments in the station
    segments = MakeSegmentsBetweenChambers(clusterStore, 2*istat, 2*istat+1);
    
    // Loop over segments
    for (Int_t iseg=0; iseg<segments->GetEntriesFast(); iseg++) 
    {
      AliDebug(1,Form("Making primary candidate(1..) %d",++iCandidate));
      
      // Transform segments to tracks and put them at the end of fRecTracksPtr
      track = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack((AliMUONObjectPair*)((*segments)[iseg]),GetRecoParam()->GetBendingVertexDispersion());
      fNRecTracks++;
      
      // Look for compatible cluster(s) in the other station
      if (GetRecoParam()->MakeTrackCandidatesFast())
	clusterFound = FollowLinearTrackInStation(*track, clusterStore, 7-istat);
      else clusterFound = FollowTrackInStation(*track, clusterStore, 7-istat);
      
      // Remove track if no cluster found on a requested station
      // or abort tracking if there are too many candidates
      if (GetRecoParam()->RequestStation(7-istat)) {
	if (!clusterFound) {
	  fRecTracksPtr->Remove(track);
	  fNRecTracks--;
	} else if (fNRecTracks > GetRecoParam()->GetMaxTrackCandidates()) {
	  AliError(Form("Too many track candidates (%d tracks). Stop tracking.", fNRecTracks));
	  return kFALSE;
	}
      } else {
	if ((fNRecTracks + segments->GetEntriesFast() - iseg - 1) > GetRecoParam()->GetMaxTrackCandidates()) {
	  AliError(Form("Too many track candidates (%d tracks). Stop tracking.", fNRecTracks + segments->GetEntriesFast() - iseg - 1));
	  return kFALSE;
	}
      }
      
    }
    
  }
  
  // Keep all different tracks or only the best ones as required
  if (GetRecoParam()->TrackAllTracks()) RemoveIdenticalTracks();
  else RemoveDoubleTracks();
  
  AliDebug(1,Form("Number of good candidates = %d",fNRecTracks));
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::MakeMoreTrackCandidates(AliMUONVClusterStore& clusterStore)
{
  /// To make extra track candidates (assuming linear propagation if the flag fgkMakeTrackCandidatesFast is set to kTRUE):
  /// clustering is supposed to be already done
  /// Start with segments made of 1 cluster in each of the stations 4 and 5 then follow track in remaining chambers.
  /// Good candidates are made of at least three clusters if both station are requested (two otherwise).
  /// Keep only best candidates or all of them according to the flag fgkTrackAllTracks.
  
  TClonesArray *segments;
  AliMUONObjectPair *segment;
  AliMUONTrack *track;
  Int_t iCandidate = 0, iCurrentTrack, nCurrentTracks;
  Bool_t clusterFound;
  
  AliDebug(1,"Enter MakeMoreTrackCandidates");
  
  // Double loop over chambers in stations(1..) 4 and 5 to make track candidates
  for (Int_t ich1 = 6; ich1 <= 7; ich1++) {
    for (Int_t ich2 = 8; ich2 <= 9; ich2++) {
      
      // Make segments in the station
      segments = MakeSegmentsBetweenChambers(clusterStore, ich1, ich2);
      
      /// Remove segments already attached to a track
      RemoveUsedSegments(*segments);
      
      // Loop over segments
      for (Int_t iSegment=0; iSegment<segments->GetEntriesFast(); iSegment++)
      {
	AliDebug(1,Form("Making primary candidate(1..) %d",++iCandidate));
	segment = (AliMUONObjectPair*) segments->UncheckedAt(iSegment);
	
	// Transform segments to tracks and put them at the end of fRecTracksPtr
	iCurrentTrack = fRecTracksPtr->GetLast()+1;
	track = new ((*fRecTracksPtr)[iCurrentTrack]) AliMUONTrack(segment,GetRecoParam()->GetBendingVertexDispersion());
	fNRecTracks++;
	
	// Look for compatible cluster(s) in the second chamber of station 5
	clusterFound = FollowLinearTrackInChamber(*track, clusterStore, 17-ich2);
	
	// skip the original track in case it has been removed
	if (GetRecoParam()->TrackAllTracks() && clusterFound) iCurrentTrack++;
	
	// loop over every new tracks
	nCurrentTracks = fRecTracksPtr->GetLast()+1;
	while (iCurrentTrack < nCurrentTracks) {
	  track = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iCurrentTrack);
	  
	  // Look for compatible cluster(s) in the second chamber of station 4
	  FollowLinearTrackInChamber(*track, clusterStore, 13-ich1);
	  
	  iCurrentTrack++;
	}
	
	// abort tracking if there are too many candidates
	if ((fNRecTracks + segments->GetEntriesFast() - iSegment - 1) > GetRecoParam()->GetMaxTrackCandidates()) {
	  AliError(Form("Too many track candidates (%d tracks). Stop tracking.", fNRecTracks + segments->GetEntriesFast() - iSegment - 1));
	  return kFALSE;
	}
	
      }
      
    }
  }
  
  // Keep only the best tracks if required
  if (!GetRecoParam()->TrackAllTracks()) RemoveDoubleTracks();
  else fRecTracksPtr->Compress();
  
  AliDebug(1,Form("Number of good candidates = %d",fNRecTracks));
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::FollowTracks(AliMUONVClusterStore& clusterStore)
{
  /// Follow tracks in stations(1..) 3, 2 and 1
  AliDebug(1,"Enter FollowTracks");
  
  AliMUONTrack *track, *nextTrack;
  AliMUONTrackParam *trackParam, *nextTrackParam;
  Int_t currentNRecTracks;
  
  Double_t sigmaCut2 = GetRecoParam()->GetSigmaCutForTracking() *
                       GetRecoParam()->GetSigmaCutForTracking();
  
  for (Int_t station = 2; station >= 0; station--) {
    
    // Save the actual number of reconstructed track in case of
    // tracks are added or suppressed during the tracking procedure
    // !! Do not compress fRecTracksPtr until the end of the loop over tracks !!
    currentNRecTracks = fNRecTracks;
    
    for (Int_t iRecTrack = 0; iRecTrack <currentNRecTracks; iRecTrack++) {
      AliDebug(1,Form("FollowTracks: track candidate(1..) %d", iRecTrack+1));
      
      track = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iRecTrack);
      
      // Fit the track:
      // Do not take into account the multiple scattering to speed up the fit
      // Calculate the track parameter covariance matrix
      // If there is no cluster out of station 4 or 5 then use the vertex to better constrain the fit
      if (((AliMUONTrackParam*) track->GetTrackParamAtCluster()->First())->GetClusterPtr()->GetChamberId() > 5)
	Fit(*track, kFALSE, kTRUE, kTRUE);
      else Fit(*track, kFALSE, kFALSE, kTRUE);
      
      // remove tracks out of limits
      if (!IsAcceptable(*((AliMUONTrackParam*)track->GetTrackParamAtCluster()->First()))) {
	fRecTracksPtr->Remove(track);
	fNRecTracks--;
	continue;
      }
      
      // remove track if the normalized chi2 is too high
      if (track->GetNormalizedChi2() > sigmaCut2) {
	fRecTracksPtr->Remove(track);
	fNRecTracks--;
	continue;
      }
      
      // save parameters from fit into smoothed parameters to complete track afterward
      if (GetRecoParam()->ComplementTracks()) {
	
	if (station==2) { // save track parameters on stations 4 and 5
	  
	  // extrapolate track parameters and covariances at each cluster
	  // remove the track in case of failure
	  if (!track->UpdateCovTrackParamAtCluster()) {
	    fRecTracksPtr->Remove(track);
	    fNRecTracks--;
	    continue;
	  }
	  
	  // save them
	  trackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->First();
	  while (trackParam) {
	    trackParam->SetSmoothParameters(trackParam->GetParameters());
	    trackParam->SetSmoothCovariances(trackParam->GetCovariances());
	    trackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->After(trackParam);
	  }
	  
	} else { // or save track parameters on last station only
	  
	  trackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->First();
	  if (trackParam->GetClusterPtr()->GetChamberId() < 2*(station+2)) {
	    
	    // save parameters from fit
	    trackParam->SetSmoothParameters(trackParam->GetParameters());
	    trackParam->SetSmoothCovariances(trackParam->GetCovariances());
	    
	    // save parameters extrapolated to the second chamber of the same station if it has been hit
	    nextTrackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->After(trackParam);
	    if (nextTrackParam->GetClusterPtr()->GetChamberId() < 2*(station+2)) {
	      
	      // reset parameters and covariances
	      nextTrackParam->SetParameters(trackParam->GetParameters());
	      nextTrackParam->SetZ(trackParam->GetZ());
	      nextTrackParam->SetCovariances(trackParam->GetCovariances());
	      
	      // extrapolate them to the z of the corresponding cluster
	      // remove the track in case of failure
	      if (!AliMUONTrackExtrap::ExtrapToZCov(nextTrackParam, nextTrackParam->GetClusterPtr()->GetZ())) {
		fRecTracksPtr->Remove(track);
		fNRecTracks--;
		continue;
	      }
	      
	      // save them
	      nextTrackParam->SetSmoothParameters(nextTrackParam->GetParameters());
	      nextTrackParam->SetSmoothCovariances(nextTrackParam->GetCovariances());
	      
	    }
	    
	  }
	  
	}
	
      }
      
      // Look for compatible cluster(s) in station(0..) "station"
      if (!FollowTrackInStation(*track, clusterStore, station)) {
	
	// Try to recover track if required
	if (GetRecoParam()->RecoverTracks()) {
	  
	  // work on a copy of the track if this station is not required
	  // to keep the case where no cluster is reconstructed as a possible candidate
	  if (!GetRecoParam()->RequestStation(station)) {
	    track = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(*track);
	    fNRecTracks++;
	  }
	  
	  // try to recover
	  if (!RecoverTrack(*track, clusterStore, station)) {
	    // remove track if no cluster found
	    fRecTracksPtr->Remove(track);
	    fNRecTracks--;
	  }
	  
	} else if (GetRecoParam()->RequestStation(station)) {
	  // remove track if no cluster found
	  fRecTracksPtr->Remove(track);
	  fNRecTracks--;
	} 
	
      }
      
      // abort tracking if there are too many candidates
      if (fNRecTracks > GetRecoParam()->GetMaxTrackCandidates()) {
	AliError(Form("Too many track candidates (%d tracks). Stop tracking.", fNRecTracks));
	return kFALSE;
      }
      
    }
    
    // Compress fRecTracksPtr for the next step
    fRecTracksPtr->Compress();
    
    // Keep only the best tracks if required
    if (!GetRecoParam()->TrackAllTracks()) RemoveDoubleTracks();
    
  }
  
  // Last fit of track candidates with all stations
  // Take into account the multiple scattering and remove bad tracks
  Int_t trackIndex = -1;
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    
    trackIndex++;
    nextTrack = (AliMUONTrack*) fRecTracksPtr->After(track); // prepare next track
    
    Fit(*track, kTRUE, kFALSE, kTRUE);
    
    // Printout for debuging
    if (AliLog::GetGlobalDebugLevel() >= 3) {
      cout << "FollowTracks: track candidate(0..) " << trackIndex << " after final fit" << endl;
      track->RecursiveDump();
    } 
    
    // Remove the track if the normalized chi2 is too high
    if (track->GetNormalizedChi2() > sigmaCut2) {
      fRecTracksPtr->Remove(track);
      fNRecTracks--;
      track = nextTrack;
      continue;
    }
    
    // save parameters from fit into smoothed parameters to complete track afterward
    if (GetRecoParam()->ComplementTracks()) {
      
      trackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->First();
      if (trackParam->GetClusterPtr()->GetChamberId() < 2) {
	
	// save parameters from fit
	trackParam->SetSmoothParameters(trackParam->GetParameters());
	trackParam->SetSmoothCovariances(trackParam->GetCovariances());
	
	// save parameters extrapolated to the second chamber of the same station if it has been hit
	nextTrackParam = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->After(trackParam);
	if (nextTrackParam->GetClusterPtr()->GetChamberId() < 2) {
	  
	  // reset parameters and covariances
	  nextTrackParam->SetParameters(trackParam->GetParameters());
	  nextTrackParam->SetZ(trackParam->GetZ());
	  nextTrackParam->SetCovariances(trackParam->GetCovariances());
	  
	  // extrapolate them to the z of the corresponding cluster
	  // remove the track in case of failure
	  if (!AliMUONTrackExtrap::ExtrapToZCov(nextTrackParam, nextTrackParam->GetClusterPtr()->GetZ())) {
	    fRecTracksPtr->Remove(track);
	    fNRecTracks--;
	    track = nextTrack;
	    continue;
	  }
	  
	  // save them
	  nextTrackParam->SetSmoothParameters(nextTrackParam->GetParameters());
	  nextTrackParam->SetSmoothCovariances(nextTrackParam->GetCovariances());
	  
	}
	
      }
      
    }
    
    track = nextTrack;
    
  }
  
  fRecTracksPtr->Compress();
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::FollowTrackInChamber(AliMUONTrack &trackCandidate, AliMUONVClusterStore& clusterStore, Int_t nextChamber)
{
  /// Follow trackCandidate in chamber(0..) nextChamber and search for compatible cluster(s)
  /// Keep all possibilities or only the best one(s) according to the flag fgkTrackAllTracks:
  /// kTRUE:  duplicate "trackCandidate" if there are several possibilities and add the new tracks at the end of
  ///         fRecTracksPtr to avoid conficts with other track candidates at this current stage of the tracking procedure.
  ///         Remove the obsolete "trackCandidate" at the end.
  /// kFALSE: add only the best cluster(s) to the "trackCandidate". Try to add a couple of clusters in priority.
  AliDebug(1,Form("Enter FollowTrackInChamber(1..) %d", nextChamber+1));
  
  Double_t chi2WithOneCluster = AliMUONTrack::MaxChi2();
  Double_t maxChi2WithOneCluster = 2. * GetRecoParam()->GetSigmaCutForTracking() *
					GetRecoParam()->GetSigmaCutForTracking(); // 2 because 2 quantities in chi2
  Double_t bestChi2WithOneCluster = maxChi2WithOneCluster;
  Bool_t foundOneCluster = kFALSE;
  AliMUONTrack *newTrack = 0x0;
  AliMUONVCluster *cluster;
  AliMUONTrackParam extrapTrackParam;
  AliMUONTrackParam extrapTrackParamAtCh;
  AliMUONTrackParam extrapTrackParamAtCluster;
  AliMUONTrackParam bestTrackParamAtCluster;
  
  // Get track parameters according to the propagation direction
  if (nextChamber > 7) extrapTrackParamAtCh = *(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->Last();
  else extrapTrackParamAtCh = *(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->First();
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    cout<<endl<<"Track parameters and covariances at first cluster:"<<endl;
    extrapTrackParamAtCh.GetParameters().Print();
    extrapTrackParamAtCh.GetCovariances().Print();
  }
  
  // Add MCS effect
  Int_t currentChamber = extrapTrackParamAtCh.GetClusterPtr()->GetChamberId();
  AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(currentChamber),-1.);
  
  // Add MCS in the missing chamber(s) if any
  while (currentChamber > nextChamber + 1) {
    // extrapolation to the missing chamber
    currentChamber--;
    if (!AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(currentChamber))) return kFALSE;
    // add MCS effect
    AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(currentChamber),-1.);
  }
  
  //Extrapolate trackCandidate to chamber
  if (!AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(nextChamber))) return kFALSE;
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    cout<<endl<<"Track parameters and covariances at first cluster extrapolated to z = "<<AliMUONConstants::DefaultChamberZ(nextChamber)<<":"<<endl;
    extrapTrackParamAtCh.GetParameters().Print();
    extrapTrackParamAtCh.GetCovariances().Print();
  }
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "FollowTrackInStation: look for clusters in chamber(1..): " << nextChamber+1 << endl;
  }
  
  // Ask the clustering to reconstruct new clusters around the track position in the current chamber
  // except for station 4 and 5 that are already entirely clusterized
  if (GetRecoParam()->CombineClusterTrackReco()) {
    if (nextChamber < 6) AskForNewClustersInChamber(extrapTrackParamAtCh, clusterStore, nextChamber);
  }
  
  // Create iterators to loop over clusters in both chambers
  TIter next(clusterStore.CreateChamberIterator(nextChamber,nextChamber));
  
  // look for cluster in chamber
  while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) ) {
    
    // try to add the current cluster fast
    if (!TryOneClusterFast(extrapTrackParamAtCh, cluster)) continue;
    
    // try to add the current cluster accuratly
    chi2WithOneCluster = TryOneCluster(extrapTrackParamAtCh, cluster, extrapTrackParamAtCluster);
    
    // if good chi2 then consider to add cluster
    if (chi2WithOneCluster < maxChi2WithOneCluster) {
      foundOneCluster = kTRUE;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	cout << "FollowTrackInStation: found one cluster in chamber(1..): " << nextChamber+1
	<< " (Chi2 = " << chi2WithOneCluster << ")" << endl;
	cluster->Print();
      }
      
      if (GetRecoParam()->TrackAllTracks()) {
	// copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
	newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	UpdateTrack(*newTrack,extrapTrackParamAtCluster);
	fNRecTracks++;
	
	// Printout for debuging
	if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	  cout << "FollowTrackInStation: added one cluster in chamber(1..): " << nextChamber+1 << endl;
	  if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	}
	
      } else if (chi2WithOneCluster < bestChi2WithOneCluster) {
	// keep track of the best single cluster except if a couple of clusters has already been found
	bestChi2WithOneCluster = chi2WithOneCluster;
	bestTrackParamAtCluster = extrapTrackParamAtCluster;
      }
      
    }
    
  }
  
  // fill out the best track if required else clean up the fRecTracksPtr array
  if (!GetRecoParam()->TrackAllTracks()) {
    if (foundOneCluster) {
      UpdateTrack(trackCandidate,bestTrackParamAtCluster);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the best cluster in chamber(1..): " << bestTrackParamAtCluster.GetClusterPtr()->GetChamberId()+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else return kFALSE;
    
  } else if (foundOneCluster) {
    
    // remove obsolete track
    fRecTracksPtr->Remove(&trackCandidate);
    fNRecTracks--;
    
  } else return kFALSE;
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::FollowTrackInStation(AliMUONTrack &trackCandidate, AliMUONVClusterStore& clusterStore, Int_t nextStation)
{
  /// Follow trackCandidate in station(0..) nextStation and search for compatible cluster(s)
  /// Keep all possibilities or only the best one(s) according to the flag fgkTrackAllTracks:
  /// kTRUE:  duplicate "trackCandidate" if there are several possibilities and add the new tracks at the end of
  ///         fRecTracksPtr to avoid conficts with other track candidates at this current stage of the tracking procedure.
  ///         Remove the obsolete "trackCandidate" at the end.
  /// kFALSE: add only the best cluster(s) to the "trackCandidate". Try to add a couple of clusters in priority.
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
  
  Double_t chi2WithOneCluster = AliMUONTrack::MaxChi2();
  Double_t chi2WithTwoClusters = AliMUONTrack::MaxChi2();
  Double_t maxChi2WithOneCluster = 2. * GetRecoParam()->GetSigmaCutForTracking() *
					GetRecoParam()->GetSigmaCutForTracking(); // 2 because 2 quantities in chi2
  Double_t maxChi2WithTwoClusters = 4. * GetRecoParam()->GetSigmaCutForTracking() *
					 GetRecoParam()->GetSigmaCutForTracking(); // 4 because 4 quantities in chi2
  Double_t bestChi2WithOneCluster = maxChi2WithOneCluster;
  Double_t bestChi2WithTwoClusters = maxChi2WithTwoClusters;
  Bool_t foundOneCluster = kFALSE;
  Bool_t foundTwoClusters = kFALSE;
  AliMUONTrack *newTrack = 0x0;
  AliMUONVCluster *clusterCh1, *clusterCh2;
  AliMUONTrackParam extrapTrackParam;
  AliMUONTrackParam extrapTrackParamAtCh;
  AliMUONTrackParam extrapTrackParamAtCluster1;
  AliMUONTrackParam extrapTrackParamAtCluster2;
  AliMUONTrackParam bestTrackParamAtCluster1;
  AliMUONTrackParam bestTrackParamAtCluster2;
  
  // Get track parameters according to the propagation direction
  if (nextStation==4) extrapTrackParamAtCh = *(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->Last();
  else extrapTrackParamAtCh = *(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->First();
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    cout<<endl<<"Track parameters and covariances at first cluster:"<<endl;
    extrapTrackParamAtCh.GetParameters().Print();
    extrapTrackParamAtCh.GetCovariances().Print();
  }
  
  // Add MCS effect
  Int_t currentChamber = extrapTrackParamAtCh.GetClusterPtr()->GetChamberId();
  AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(currentChamber),-1.);
  
  // Add MCS in the missing chamber(s) if any
  while (ch1 < ch2 && currentChamber > ch2 + 1) {
    // extrapolation to the missing chamber
    currentChamber--;
    if (!AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(currentChamber))) return kFALSE;
    // add MCS effect
    AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(currentChamber),-1.);
  }
  
  //Extrapolate trackCandidate to chamber "ch2"
  if (!AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(ch2))) return kFALSE;
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    cout<<endl<<"Track parameters and covariances at first cluster extrapolated to z = "<<AliMUONConstants::DefaultChamberZ(ch2)<<":"<<endl;
    extrapTrackParamAtCh.GetParameters().Print();
    extrapTrackParamAtCh.GetCovariances().Print();
  }
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "FollowTrackInStation: look for clusters in chamber(1..): " << ch2+1 << endl;
  }
  
  // Ask the clustering to reconstruct new clusters around the track position in the current station
  // except for station 4 and 5 that are already entirely clusterized
  if (GetRecoParam()->CombineClusterTrackReco()) {
    if (nextStation < 3) AskForNewClustersInStation(extrapTrackParamAtCh, clusterStore, nextStation);
  }
  
  Int_t nClusters = clusterStore.GetSize();
  Bool_t *clusterCh1Used = new Bool_t[nClusters];
  for (Int_t i = 0; i < nClusters; i++) clusterCh1Used[i] = kFALSE;
  Int_t iCluster1;
  
  // Create iterators to loop over clusters in both chambers
  TIter nextInCh1(clusterStore.CreateChamberIterator(ch1,ch1));
  TIter nextInCh2(clusterStore.CreateChamberIterator(ch2,ch2));
  
  // look for candidates in chamber 2
  while ( ( clusterCh2 = static_cast<AliMUONVCluster*>(nextInCh2()) ) ) {
    
    // try to add the current cluster fast
    if (!TryOneClusterFast(extrapTrackParamAtCh, clusterCh2)) continue;
    
    // try to add the current cluster accuratly
    chi2WithOneCluster = TryOneCluster(extrapTrackParamAtCh, clusterCh2, extrapTrackParamAtCluster2);
    
    // if good chi2 then try to attach a cluster in the other chamber too
    if (chi2WithOneCluster < maxChi2WithOneCluster) {
      Bool_t foundSecondCluster = kFALSE;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: found one cluster in chamber(1..): " << ch2+1
	     << " (Chi2 = " << chi2WithOneCluster << ")" << endl;
	clusterCh2->Print();
        cout << "                      look for second clusters in chamber(1..): " << ch1+1 << " ..." << endl;
      }
      
      // add MCS effect for next step
      AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCluster2,AliMUONConstants::ChamberThicknessInX0(ch2),-1.);
      
      // copy new track parameters for next step
      extrapTrackParam = extrapTrackParamAtCluster2;
      
      //Extrapolate track parameters to chamber "ch1"
      Bool_t normalExtrap = AliMUONTrackExtrap::ExtrapToZ(&extrapTrackParam, AliMUONConstants::DefaultChamberZ(ch1));
      
      // reset cluster iterator of chamber 1
      nextInCh1.Reset();
      iCluster1 = -1;
      
      // look for second candidates in chamber 1
      if (normalExtrap) while ( ( clusterCh1 = static_cast<AliMUONVCluster*>(nextInCh1()) ) ) {
        iCluster1++;
	
    	// try to add the current cluster fast
    	if (!TryOneClusterFast(extrapTrackParam, clusterCh1)) continue;
    	
    	// try to add the current cluster accurately
	chi2WithTwoClusters = TryTwoClusters(extrapTrackParamAtCluster2, clusterCh1, extrapTrackParamAtCluster1);
        
	// if good chi2 then create a new track by adding the 2 clusters to the "trackCandidate"
	if (chi2WithTwoClusters < maxChi2WithTwoClusters) {
	  foundSecondCluster = kTRUE;
          foundTwoClusters = kTRUE;
          
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowTrackInStation: found second cluster in chamber(1..): " << ch1+1
	         << " (Global Chi2 = " << chi2WithTwoClusters << ")" << endl;
	    clusterCh1->Print();
	  }
	  
	  if (GetRecoParam()->TrackAllTracks()) {
	    // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new clusters
            newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	    UpdateTrack(*newTrack,extrapTrackParamAtCluster1,extrapTrackParamAtCluster2);
	    fNRecTracks++;
	    
	    // Tag clusterCh1 as used
	    clusterCh1Used[iCluster1] = kTRUE;
	    
	    // Printout for debuging
	    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	      cout << "FollowTrackInStation: added two clusters in station(1..): " << nextStation+1 << endl;
	      if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	    }
	    
          } else if (chi2WithTwoClusters < bestChi2WithTwoClusters) {
	    // keep track of the best couple of clusters
	    bestChi2WithTwoClusters = chi2WithTwoClusters;
	    bestTrackParamAtCluster1 = extrapTrackParamAtCluster1;
	    bestTrackParamAtCluster2 = extrapTrackParamAtCluster2;
          }
	  
	}
	
      }
      
      // if no clusterCh1 found then consider to add clusterCh2 only
      if (!foundSecondCluster) {
        foundOneCluster = kTRUE;
        
	if (GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
          newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  UpdateTrack(*newTrack,extrapTrackParamAtCluster2);
	  fNRecTracks++;
	  
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowTrackInStation: added one cluster in chamber(1..): " << ch2+1 << endl;
	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	  }
	  
	} else if (!foundTwoClusters && chi2WithOneCluster < bestChi2WithOneCluster) {
	  // keep track of the best single cluster except if a couple of clusters has already been found
	  bestChi2WithOneCluster = chi2WithOneCluster;
	  bestTrackParamAtCluster1 = extrapTrackParamAtCluster2;
        }
	
      }
      
    }
    
  }
  
  // look for candidates in chamber 1 not already attached to a track
  // if we want to keep all possible tracks or if no good couple of clusters has been found
  if (GetRecoParam()->TrackAllTracks() || !foundTwoClusters) {
    
    // add MCS effect for next step
    AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(ch2),-1.);
    
    //Extrapolate trackCandidate to chamber "ch1"
    Bool_t normalExtrap = AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(ch1));
    
    // Printout for debuging
    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
      cout<<endl<<"Track parameters and covariances at first cluster extrapolated to z = "<<AliMUONConstants::DefaultChamberZ(ch1)<<":"<<endl;
      extrapTrackParamAtCh.GetParameters().Print();
      extrapTrackParamAtCh.GetCovariances().Print();
    }
    
    // Printout for debuging
    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
      cout << "FollowTrackInStation: look for single clusters in chamber(1..): " << ch1+1 << endl;
    }
    
    // reset cluster iterator of chamber 1
    nextInCh1.Reset();
    iCluster1 = -1;
    
    // look for second candidates in chamber 1
    if (normalExtrap) while ( ( clusterCh1 = static_cast<AliMUONVCluster*>(nextInCh1()) ) ) {
      iCluster1++;
      
      if (clusterCh1Used[iCluster1]) continue; // Skip cluster already used
      
      // try to add the current cluster fast
      if (!TryOneClusterFast(extrapTrackParamAtCh, clusterCh1)) continue;
      
      // try to add the current cluster accuratly
      chi2WithOneCluster = TryOneCluster(extrapTrackParamAtCh, clusterCh1, extrapTrackParamAtCluster1);
      
      // if good chi2 then consider to add clusterCh1
      // We do not try to attach a cluster in the other chamber too since it has already been done above
      if (chi2WithOneCluster < maxChi2WithOneCluster) {
        foundOneCluster = kTRUE;
	
	// Printout for debuging
  	if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	  cout << "FollowTrackInStation: found one cluster in chamber(1..): " << ch1+1
  	       << " (Chi2 = " << chi2WithOneCluster << ")" << endl;
	  clusterCh1->Print();
  	}
	
  	if (GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
  	  newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  UpdateTrack(*newTrack,extrapTrackParamAtCluster1);
	  fNRecTracks++;
  	  
	  // Printout for debuging
  	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	    cout << "FollowTrackInStation: added one cluster in chamber(1..): " << ch1+1 << endl;
  	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
  	  }
	  
  	} else if (chi2WithOneCluster < bestChi2WithOneCluster) {
	  // keep track of the best single cluster except if a couple of clusters has already been found
  	  bestChi2WithOneCluster = chi2WithOneCluster;
	  bestTrackParamAtCluster1 = extrapTrackParamAtCluster1;
  	}
	
      }
      
    }
    
  }
  
  // fill out the best track if required else clean up the fRecTracksPtr array
  if (!GetRecoParam()->TrackAllTracks()) {
    if (foundTwoClusters) {
      UpdateTrack(trackCandidate,bestTrackParamAtCluster1,bestTrackParamAtCluster2);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the two best clusters in station(1..): " << nextStation+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else if (foundOneCluster) {
      UpdateTrack(trackCandidate,bestTrackParamAtCluster1);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
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
Double_t AliMUONTrackReconstructor::TryTwoClusters(const AliMUONTrackParam &trackParamAtCluster1, AliMUONVCluster* cluster2,
						   AliMUONTrackParam &trackParamAtCluster2)
{
/// Test the compatibility between the track and the 2 clusters together (using trackParam's covariance matrix):
/// return the corresponding Chi2 accounting for covariances between the 2 clusters
/// return trackParamAtCluster1 & 2
  
  // extrapolate track parameters at the z position of the second cluster (no need to extrapolate the covariances)
  // and set pointer to cluster into trackParamAtCluster
  trackParamAtCluster2.SetParameters(trackParamAtCluster1.GetParameters());
  trackParamAtCluster2.SetZ(trackParamAtCluster1.GetZ());
  trackParamAtCluster2.SetClusterPtr(cluster2);
  if (!AliMUONTrackExtrap::ExtrapToZ(&trackParamAtCluster2, cluster2->GetZ())) return 2.*AliMUONTrack::MaxChi2();
  
  // Set differences between track and the 2 clusters in the bending and non bending directions
  AliMUONVCluster* cluster1 = trackParamAtCluster1.GetClusterPtr();
  TMatrixD dPos(4,1);
  dPos(0,0) = cluster1->GetX() - trackParamAtCluster1.GetNonBendingCoor();
  dPos(1,0) = cluster1->GetY() - trackParamAtCluster1.GetBendingCoor();
  dPos(2,0) = cluster2->GetX() - trackParamAtCluster2.GetNonBendingCoor();
  dPos(3,0) = cluster2->GetY() - trackParamAtCluster2.GetBendingCoor();
  
  // Calculate the error matrix from the track parameter covariances at first cluster
  TMatrixD error(4,4);
  error.Zero();
  if (trackParamAtCluster1.CovariancesExist()) {
    // Save track parameters at first cluster
    TMatrixD paramAtCluster1Save(trackParamAtCluster1.GetParameters());
    
    // Save track coordinates at second cluster
    Double_t nonBendingCoor2 = trackParamAtCluster2.GetNonBendingCoor();
    Double_t bendingCoor2    = trackParamAtCluster2.GetBendingCoor();
    
    // copy track parameters at first cluster for jacobian calculation
    AliMUONTrackParam trackParam(trackParamAtCluster1);
    
    // Get the pointer to the parameter covariance matrix at first cluster
    const TMatrixD& kParamCov = trackParam.GetCovariances();
    
    // Calculate the jacobian related to the transformation between track parameters
    // at first cluster and track coordinates at the 2 cluster z-positions
    TMatrixD jacob(4,5);
    jacob.Zero();
    // first derivative at the first cluster:
    jacob(0,0) = 1.; // dx1/dx
    jacob(1,2) = 1.; // dy1/dy
    // first derivative at the second cluster:
    TMatrixD dParam(5,1);
    Double_t direction[5] = {-1.,-1.,1.,1.,-1.};
    for (Int_t i=0; i<5; i++) {
      // Skip jacobian calculation for parameters with no associated error
      if (kParamCov(i,i) == 0.) continue;
      // Small variation of parameter i only
      for (Int_t j=0; j<5; j++) {
        if (j==i) {
          dParam(j,0) = TMath::Sqrt(kParamCov(i,i));
	  dParam(j,0) *= TMath::Sign(1.,direction[j]*paramAtCluster1Save(j,0)); // variation always in the same direction
        } else dParam(j,0) = 0.;
      }
      
      // Set new track parameters at first cluster
      trackParam.SetParameters(paramAtCluster1Save);
      trackParam.AddParameters(dParam);
      trackParam.SetZ(cluster1->GetZ());
      
      // Extrapolate new track parameters to the z position of the second cluster
      if (!AliMUONTrackExtrap::ExtrapToZ(&trackParam,cluster2->GetZ())) return 2.*AliMUONTrack::MaxChi2();
      
      // Calculate the jacobian
      jacob(2,i) = (trackParam.GetNonBendingCoor() - nonBendingCoor2) / dParam(i,0); // dx2/dParami
      jacob(3,i) = (trackParam.GetBendingCoor()    - bendingCoor2   ) / dParam(i,0); // dy2/dParami
    }
    
    // Calculate the error matrix
    TMatrixD tmp(jacob,TMatrixD::kMult,kParamCov);
    error = TMatrixD(tmp,TMatrixD::kMultTranspose,jacob);
  }
  
  // Add cluster resolution to the error matrix
  error(0,0) += cluster1->GetErrX2();
  error(1,1) += cluster1->GetErrY2();
  error(2,2) += cluster2->GetErrX2();
  error(3,3) += cluster2->GetErrY2();
  
  // invert the error matrix for Chi2 calculation
  if (error.Determinant() != 0) {
    error.Invert();
  } else {
    AliWarning(" Determinant error=0");
    return 2.*AliMUONTrack::MaxChi2();
  }
  
  // Compute the Chi2 value
  TMatrixD tmp2(dPos,TMatrixD::kTransposeMult,error);
  TMatrixD result(tmp2,TMatrixD::kMult,dPos);
  
  return result(0,0);
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtCluster)
{
  /// Add 1 cluster to the track candidate
  /// Update chi2 of the track 
  
  // Compute local chi2
  AliMUONVCluster* cluster = trackParamAtCluster.GetClusterPtr();
  Double_t deltaX = trackParamAtCluster.GetNonBendingCoor() - cluster->GetX();
  Double_t deltaY = trackParamAtCluster.GetBendingCoor() - cluster->GetY();
  Double_t localChi2 = deltaX*deltaX / cluster->GetErrX2() +
  		       deltaY*deltaY / cluster->GetErrY2();
  
  // Flag cluster as being not removable
  if (GetRecoParam()->RequestStation(cluster->GetChamberId()/2))
    trackParamAtCluster.SetRemovable(kFALSE);
  else trackParamAtCluster.SetRemovable(kTRUE);
  trackParamAtCluster.SetLocalChi2(0.); // --> Local chi2 not used
  
  // Update the chi2 of the new track
  track.SetGlobalChi2(track.GetGlobalChi2() + localChi2);
  
  // Update TrackParamAtCluster
  track.AddTrackParamAtCluster(trackParamAtCluster,*cluster);
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtCluster1, AliMUONTrackParam &trackParamAtCluster2)
{
  /// Add 2 clusters to the track candidate
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
  deltaX = trackParamAtCluster2.GetNonBendingCoor() - cluster2->GetX();
  deltaY = trackParamAtCluster2.GetBendingCoor() - cluster2->GetY();
  Double_t localChi2AtCluster2 = deltaX*deltaX / cluster2->GetErrX2() +
  			         deltaY*deltaY / cluster2->GetErrY2();
  trackParamAtCluster2.SetLocalChi2(localChi2AtCluster2);
  
  // Flag first cluster as being removable
  trackParamAtCluster2.SetRemovable(kTRUE);
  
  // Update the chi2 of the new track
  track.SetGlobalChi2(track.GetGlobalChi2() + localChi2AtCluster1 + localChi2AtCluster2);
  
  // Update TrackParamAtCluster
  track.AddTrackParamAtCluster(trackParamAtCluster1,*cluster1);
  track.AddTrackParamAtCluster(trackParamAtCluster2,*cluster2);
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::RecoverTrack(AliMUONTrack &trackCandidate, AliMUONVClusterStore& clusterStore, Int_t nextStation)
{
  /// Try to recover the track candidate in the next station
  /// by removing the worst of the two clusters attached in the current station
  /// Return kTRUE if recovering succeeds
  AliDebug(1,"Enter RecoverTrack");
  
  // Do not try to recover track until we have attached cluster(s) on station(1..) 3
  if (nextStation > 1) return kFALSE;
  
  Int_t worstClusterNumber = -1;
  Double_t localChi2, worstLocalChi2 = -1.;
  
  // Look for the cluster to remove
  for (Int_t clusterNumber = 0; clusterNumber < 2; clusterNumber++) {
    AliMUONTrackParam *trackParamAtCluster = (AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->UncheckedAt(clusterNumber);
    
    // check if current cluster is in the previous station
    if (trackParamAtCluster->GetClusterPtr()->GetChamberId()/2 != nextStation+1) break;
    
    // check if current cluster is removable
    if (!trackParamAtCluster->IsRemovable()) return kFALSE;
    
    // reset the current cluster as beig not removable if it is on a required station
    if (GetRecoParam()->RequestStation(nextStation+1)) trackParamAtCluster->SetRemovable(kFALSE);
    
    // Pick up cluster with the worst chi2
    localChi2 = trackParamAtCluster->GetLocalChi2();
    if (localChi2 > worstLocalChi2) {
      worstLocalChi2 = localChi2;
      worstClusterNumber = clusterNumber;
    }
    
  }
  
  // check if worst cluster found
  if (worstClusterNumber < 0) return kFALSE;
  
  // Remove the worst cluster
  trackCandidate.RemoveTrackParamAtCluster((AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->UncheckedAt(worstClusterNumber));
  
  // Re-fit the track:
  // Do not take into account the multiple scattering to speed up the fit
  // Calculate the track parameter covariance matrix
  Fit(trackCandidate, kFALSE, kFALSE, kTRUE);
  
  // skip track if the normalized chi2 is too high
  if (trackCandidate.GetNormalizedChi2() > GetRecoParam()->GetSigmaCutForTracking() * GetRecoParam()->GetSigmaCutForTracking()) return kFALSE;
  
  // skip track out of limits
  if (!IsAcceptable(*((AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->First()))) return kFALSE;
  
  // Look for new cluster(s) in next station
  return FollowTrackInStation(trackCandidate,clusterStore,nextStation);
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::SetVertexErrXY2ForFit(AliMUONTrack &trackCandidate)
{
  /// Compute the vertex resolution square from natural vertex dispersion and
  /// multiple scattering effets according to trackCandidate path in absorber
  /// It is necessary to account for multiple scattering effects here instead of during the fit of
  /// the "trackCandidate" to do not influence the result by changing track resolution at vertex
  AliDebug(1,"Enter SetVertexForFit");
  
  Double_t nonBendingReso2 = GetRecoParam()->GetNonBendingVertexDispersion() *
                             GetRecoParam()->GetNonBendingVertexDispersion();
  Double_t bendingReso2 = GetRecoParam()->GetBendingVertexDispersion() *
			  GetRecoParam()->GetBendingVertexDispersion();
  
  // add multiple scattering effets
  AliMUONTrackParam paramAtVertex(*((AliMUONTrackParam*)(trackCandidate.GetTrackParamAtCluster()->First())));
  paramAtVertex.DeleteCovariances(); // to be sure to account only for multiple scattering
  if (!AliMUONTrackExtrap::ExtrapToZ(&paramAtVertex,AliMUONConstants::AbsZEnd())) {
    nonBendingReso2 = 0.;
    bendingReso2 = 0.;
  } else {
    AliMUONTrackExtrap::ExtrapToVertexUncorrected(&paramAtVertex,0.);
    const TMatrixD& kParamCov = paramAtVertex.GetCovariances();
    nonBendingReso2 += kParamCov(0,0);
    bendingReso2 += kParamCov(2,2);
  }
  
  // Set the vertex resolution square
  trackCandidate.SetVertexErrXY2(nonBendingReso2,bendingReso2);
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::Fit(AliMUONTrack &track, Bool_t includeMCS, Bool_t fitWithVertex, Bool_t calcCov)
{
  /// Fit the track
  /// w/wo multiple Coulomb scattering according to "includeMCS".
  /// w/wo constraining the vertex according to "fitWithVertex".
  /// calculating or not the covariance matrix according to "calcCov".
  AliDebug(1,"Enter Fit");
  
  Double_t benC, errorParam, invBenP, nonBenC, x, y;
  AliMUONTrackParam *trackParam;
  Double_t arg[1], fedm, errdef, globalChi2;
  Int_t npari, nparx;
  Int_t status, covStatus;
  
  // Instantiate gMinuit if not already done
  if (!gMinuit) gMinuit = new TMinuit(6);
  // Clear MINUIT parameters
  gMinuit->mncler();
  // Give the fitted track to MINUIT
  gMinuit->SetObjectFit(&track);
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
  
  // set flag w/wo multiple scattering according to "includeMCS"
  track.FitWithMCS(includeMCS);
  if (includeMCS) {
    // compute cluster weights only once
    if (!track.UpdateTrackParamAtCluster() || !track.ComputeClusterWeights()) {
      AliWarning("cannot take into account the multiple scattering effects");
      track.FitWithMCS(kFALSE);
    }
  }
  
  track.FitWithVertex(fitWithVertex);
  if (fitWithVertex) SetVertexErrXY2ForFit(track);
  
  // Set fitting function
  gMinuit->SetFCN(TrackChi2);
  
  // Set fitted parameters (!! The order is very important for the covariance matrix !!)
  // Mandatory limits to avoid NaN values of parameters
  trackParam = (AliMUONTrackParam*) (track.GetTrackParamAtCluster()->First());
  Double_t maxIBM = 1. / GetRecoParam()->GetMinBendingMomentum();
  gMinuit->mnparm(0, "X", trackParam->GetNonBendingCoor(), 0.03, -500.0, 500.0, status);
  gMinuit->mnparm(1, "NonBenS", trackParam->GetNonBendingSlope(), 0.001, -1., 1., status);
  gMinuit->mnparm(2, "Y", trackParam->GetBendingCoor(), 0.10, -500.0, 500.0, status);
  gMinuit->mnparm(3, "BenS", trackParam->GetBendingSlope(), 0.001, -1.5, 1.5, status);
  gMinuit->mnparm(4, "InvBenP", trackParam->GetInverseBendingMomentum(), 0.003, -maxIBM, maxIBM, status);
  
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
  gMinuit->mnstat(globalChi2, fedm, errdef, npari, nparx, covStatus);
  track.SetGlobalChi2(globalChi2);
  
  // Get the covariance matrix if required
  if (calcCov) {
    // Covariance matrix according to HESSE status
    // If problem then keep only the diagonal terms (variances)
    Double_t matrix[5][5];
    for (Int_t i=0; i<5; i++) for (Int_t j=0; j<5; j++) matrix[i][j] = 0.;
    gMinuit->mnemat(&matrix[0][0],5);
    if (covStatus == 3) trackParam->SetCovariances(matrix);
    else trackParam->SetVariances(matrix);
  } else trackParam->DeleteCovariances();
  
}

  //__________________________________________________________________________
void TrackChi2(Int_t & /*nParam*/, Double_t * /*gradient*/, Double_t &chi2, Double_t *param, Int_t /*flag*/)
{
  /// Return the "Chi2" to be minimized with Minuit for track fitting.
  /// Assumes that the trackParamAtCluster are sorted according to increasing Z.
  /// Track parameters at each cluster are updated accordingly.
  /// Vertex is used according to the flag "trackBeingFitted->GetFitWithVertex()".
  /// Multiple Coulomb scattering is taken into account according to the flag "trackBeingFitted->GetFitWithMCS()".
  
  AliMUONTrack *trackBeingFitted = (AliMUONTrack*) gMinuit->GetObjectFit();
  AliMUONTrackParam* trackParamAtCluster = (AliMUONTrackParam*) trackBeingFitted->GetTrackParamAtCluster()->First();
  Double_t dX, dY;
  chi2 = 0.; // initialize chi2
  
  // update track parameters
  trackParamAtCluster->SetNonBendingCoor(param[0]);
  trackParamAtCluster->SetNonBendingSlope(param[1]);
  trackParamAtCluster->SetBendingCoor(param[2]);
  trackParamAtCluster->SetBendingSlope(param[3]);
  trackParamAtCluster->SetInverseBendingMomentum(param[4]);
  if (!trackBeingFitted->UpdateTrackParamAtCluster()) {
    chi2 = 2.*AliMUONTrack::MaxChi2();
    return;
  }
  
  // Take the vertex into account in the fit if required
  if (trackBeingFitted->FitWithVertex()) {
    Double_t nonBendingReso2,bendingReso2;
    trackBeingFitted->GetVertexErrXY2(nonBendingReso2,bendingReso2);
    AliMUONTrackParam paramAtVertex(*trackParamAtCluster);
    if (nonBendingReso2 != 0. && bendingReso2 != 0. && AliMUONTrackExtrap::ExtrapToZ(&paramAtVertex, 0.)) { // vextex position = (0,0,0)
      dX = paramAtVertex.GetNonBendingCoor();
      dY = paramAtVertex.GetBendingCoor();
      chi2 += dX * dX / nonBendingReso2 + dY * dY / bendingReso2;
    } else {
      chi2 = 2.*AliMUONTrack::MaxChi2();
      return;
    }
  }
  
  // compute chi2 w/wo multiple scattering
  chi2 += trackBeingFitted->ComputeGlobalChi2(trackBeingFitted->FitWithMCS());
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::ComplementTracks(const AliMUONVClusterStore& clusterStore)
{
  /// Complete tracks by adding missing clusters (if there is an overlap between
  /// two detection elements, the track may have two clusters in the same chamber).
  /// Re-fit track parameters and covariances.
  /// Return kTRUE if one or more tracks have been complemented.
  AliDebug(1,"Enter ComplementTracks");
  
  Int_t chamberId, detElemId;
  Double_t chi2OfCluster, bestChi2OfCluster;
  Double_t sigmaCut2 = GetRecoParam()->GetSigmaCutForTracking() *
                       GetRecoParam()->GetSigmaCutForTracking();
  Bool_t foundOneCluster, trackModified, hasChanged = kFALSE;
  AliMUONVCluster* cluster;
  AliMUONTrackParam *trackParam, *nextTrackParam, copyOfTrackParam, trackParamAtCluster, bestTrackParamAtCluster;
  
  AliMUONTrack *track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    trackModified = kFALSE;
    
    trackParam = (AliMUONTrackParam*)track->GetTrackParamAtCluster()->First();
    while (trackParam) {
      foundOneCluster = kFALSE;
      bestChi2OfCluster = 2. * sigmaCut2; // 2 because 2 quantities in chi2
      chamberId = trackParam->GetClusterPtr()->GetChamberId();
      detElemId = trackParam->GetClusterPtr()->GetDetElemId();
      
      // prepare nextTrackParam before adding new cluster because of the sorting
      nextTrackParam = (AliMUONTrackParam*)track->GetTrackParamAtCluster()->After(trackParam);
      
      // recover track parameters from local fit and put them into a copy of trackParam
      copyOfTrackParam.SetZ(trackParam->GetZ());
      copyOfTrackParam.SetParameters(trackParam->GetSmoothParameters());
      copyOfTrackParam.SetCovariances(trackParam->GetSmoothCovariances());
      
      // Create iterators to loop over clusters in current chamber
      TIter nextInCh(clusterStore.CreateChamberIterator(chamberId,chamberId));
      
      // look for one second candidate in the same chamber
      while ( ( cluster = static_cast<AliMUONVCluster*>(nextInCh()) ) ) {
        
	// look for a cluster in another detection element
	if (cluster->GetDetElemId() == detElemId) continue;
	
	// try to add the current cluster fast
	if (!TryOneClusterFast(copyOfTrackParam, cluster)) continue;
	
	// try to add the current cluster accurately
	chi2OfCluster = TryOneCluster(copyOfTrackParam, cluster, trackParamAtCluster);
	
	// if better chi2 then prepare to add this cluster to the track
	if (chi2OfCluster < bestChi2OfCluster) {
	  bestChi2OfCluster = chi2OfCluster;
	  bestTrackParamAtCluster = trackParamAtCluster;
	  foundOneCluster = kTRUE;
	}
	
      }
      
      // add new cluster if any
      if (foundOneCluster) {
	
	// Printout for debuging
	if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	  cout << "ComplementTracks: found one cluster in chamber(1..): " << chamberId+1 << endl;
	  bestTrackParamAtCluster.GetClusterPtr()->Print();
	  cout<<endl<<"Track parameters and covariances at cluster:"<<endl;
	  bestTrackParamAtCluster.GetParameters().Print();
	  bestTrackParamAtCluster.GetCovariances().Print();
	}
	
	trackParam->SetRemovable(kTRUE);
	bestTrackParamAtCluster.SetRemovable(kTRUE);
	track->AddTrackParamAtCluster(bestTrackParamAtCluster,*(bestTrackParamAtCluster.GetClusterPtr()));
	trackModified = kTRUE;
	hasChanged = kTRUE;
      }
      
      trackParam = nextTrackParam;
    }
    
    // re-fit track parameters if needed
    if (trackModified) Fit(*track, kTRUE, kFALSE, kTRUE);
    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
  }
  
  return hasChanged;
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::ImproveTrack(AliMUONTrack &track)
{
  /// Improve the given track by removing clusters with local chi2 highter than the defined cut
  /// Recompute track parameters and covariances at the remaining clusters
  AliDebug(1,"Enter ImproveTrack");
  
  Double_t localChi2, worstLocalChi2;
  AliMUONTrackParam *trackParamAtCluster, *worstTrackParamAtCluster;
  Double_t sigmaCut2 = GetRecoParam()->GetSigmaCutForImprovement() *
		       GetRecoParam()->GetSigmaCutForImprovement();
  
  while (!track.IsImproved()) {
    
    // identify removable clusters
    track.TagRemovableClusters(GetRecoParam()->RequestedStationMask());
    
    // Update track parameters and covariances
    if (!track.UpdateCovTrackParamAtCluster()) {
      AliWarning("unable to update track parameters and covariances --> stop improvement");
      break;
    }
    
    // Compute local chi2 of each clusters
    track.ComputeLocalChi2(kTRUE);
    
    // Look for the cluster to remove
    worstTrackParamAtCluster = NULL;
    worstLocalChi2 = 0.;
    trackParamAtCluster = (AliMUONTrackParam*) track.GetTrackParamAtCluster()->First();
    while (trackParamAtCluster) {
      
      // Pick up cluster with the worst chi2
      localChi2 = trackParamAtCluster->GetLocalChi2();
      if (localChi2 > worstLocalChi2) {
	worstLocalChi2 = localChi2;
	worstTrackParamAtCluster = trackParamAtCluster;
      }
      
    trackParamAtCluster = (AliMUONTrackParam*) track.GetTrackParamAtCluster()->After(trackParamAtCluster);
    }
    
    // Check if worst cluster found
    if (!worstTrackParamAtCluster) {
      AliWarning("Bad local chi2 values?");
      break;
    }
    
    // Check whether the worst chi2 is under requirement or not
    if (worstLocalChi2 < 2. * sigmaCut2) { // 2 because 2 quantities in chi2
      track.SetImproved(kTRUE);
      break;
    }
    
    // if the worst cluster is not removable then stop improvement
    if (!worstTrackParamAtCluster->IsRemovable()) break;
    
    // Remove the worst cluster
    track.RemoveTrackParamAtCluster(worstTrackParamAtCluster);
    
    // Re-fit the track:
    // Take into account the multiple scattering
    // Calculate the track parameter covariance matrix
    Fit(track, kTRUE, kFALSE, kTRUE);
    
    // Printout for debuging
    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
      cout << "ImproveTracks: track " << fRecTracksPtr->IndexOf(&track)+1 << " improved " << endl;
    }
    
  }
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::FinalizeTrack(AliMUONTrack &track)
{
  /// Recompute track parameters and covariances at each attached cluster
  /// from those at the first one, if not already done
  AliDebug(1,"Enter FinalizeTrack");
  if (!track.IsImproved() && !track.UpdateCovTrackParamAtCluster()) {
    AliWarning("finalization failed due to extrapolation problem");
    return kFALSE;
  }
  return kTRUE;
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::RefitTrack(AliMUONTrack &track, Bool_t enableImprovement)
{
  /// re-fit the given track
  
  // check validity of the track
  if (track.GetNClusters() < 3) {
    AliWarning("the track does not contain enough clusters --> unable to refit");
    return kFALSE;
  }
  
  // reset the seed (i.e. parameters at first cluster) before fitting
  AliMUONTrackParam* firstTrackParam = (AliMUONTrackParam*) track.GetTrackParamAtCluster()->First();
  if (firstTrackParam->GetInverseBendingMomentum() == 0.) {
    AliWarning("track parameters at first chamber are not initialized --> unable to refit");
    return kFALSE;
  }
  
  // compute track parameters at each cluster from parameters at the first one
  // necessary to compute multiple scattering effect during refitting
  if (!track.UpdateTrackParamAtCluster()) {
    AliWarning("bad track refitting due to extrapolation failure");
    return kFALSE;
  }
  
  // Re-fit the track:
  // Take into account the multiple scattering
  // Calculate the track parameter covariance matrix
  Fit(track, kTRUE, kFALSE, kTRUE);
  
  // Improve the reconstructed tracks if required
  track.SetImproved(kFALSE);
  if (enableImprovement && GetRecoParam()->ImproveTracks()) ImproveTrack(track);
  
  // Fill AliMUONTrack data members
  if (track.GetGlobalChi2() < AliMUONTrack::MaxChi2()) return FinalizeTrack(track);
  else {
    AliWarning("track not finalized due to extrapolation failure");
    return kFALSE;
  }
  
}

