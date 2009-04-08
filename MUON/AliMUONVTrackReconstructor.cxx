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
/// This class contains as data a pointer to the array of reconstructed tracks
///
/// It contains as methods, among others:
/// * EventReconstruct to build the muon tracks
/// * EventReconstructTrigger to build the trigger tracks
/// * ValidateTracksWithTrigger to match tracker/trigger tracks
///
/// Several options and adjustable parameters are available for both KALMAN and ORIGINAL
/// tracking algorithms. They can be changed through the AliMUONRecoParam object
/// set in the reconstruction macro or read from the CDB
/// (see methods in AliMUONRecoParam.h file for details)
///
/// Main parameters and options are:
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
///   new quality cut (the track is removed is it does not contain enough cluster anymore).
/// - *fgkComplementTracks* : if this flag is set to 'true', we try to improve the quality
///   of the tracks at the end of the tracking by adding potentially missing clusters
///   (we may have 2 clusters in the same chamber because of the overlapping of detection
///   elements, which is not handle by the tracking algorithm).
/// - *fgkSigmaToCutForImprovement* : quality cut used when we try to improve the
///   quality of the tracks.
///
///  \author Philippe Pillot
//-----------------------------------------------------------------------------

#include "AliMUONVTrackReconstructor.h"

#include "AliMUONConstants.h"
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
#include "AliMUONVClusterServer.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONRecoParam.h"

#include "AliMpDEManager.h"
#include "AliMpArea.h"

#include "AliLog.h"
#include "AliCodeTimer.h"
#include "AliTracker.h"

#include <TClonesArray.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVector2.h>

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONVTrackReconstructor) // Class implementation in ROOT context
/// \endcond

  //__________________________________________________________________________
AliMUONVTrackReconstructor::AliMUONVTrackReconstructor(const AliMUONRecoParam* recoParam,
                                                       AliMUONVClusterServer* clusterServer)
: TObject(),
fRecTracksPtr(0x0),
fNRecTracks(0),
fClusterServer(clusterServer),
fkRecoParam(recoParam)
{
  /// Constructor for class AliMUONVTrackReconstructor
  /// WARNING: if clusterServer=0x0, no clusterization will be possible at this level
  
  // Memory allocation for the TClonesArray of reconstructed tracks
  fRecTracksPtr = new TClonesArray("AliMUONTrack", 100);
  
  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField();
}

  //__________________________________________________________________________
AliMUONVTrackReconstructor::~AliMUONVTrackReconstructor()
{
  /// Destructor for class AliMUONVTrackReconstructor
  delete fRecTracksPtr;
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
void AliMUONVTrackReconstructor::EventReconstruct(AliMUONVClusterStore& clusterStore, AliMUONVTrackStore& trackStore)
{
  /// To reconstruct one event
  AliDebug(1,"");
  AliCodeTimerAuto("");
  
  // Reset array of tracks
  ResetTracks();
  
  // Look for candidates from clusters in stations(1..) 4 and 5 (abort in case of failure)
  if (!MakeTrackCandidates(clusterStore)) return;
  
  // Look for extra candidates from clusters in stations(1..) 4 and 5 (abort in case of failure)
  if (GetRecoParam()->MakeMoreTrackCandidates()) {
    if (!MakeMoreTrackCandidates(clusterStore)) return;
  }
  
  // Stop tracking if no candidate found
  if (fRecTracksPtr->GetEntriesFast() == 0) return;
  
  // Follow tracks in stations(1..) 3, 2 and 1 (abort in case of failure)
  if (!FollowTracks(clusterStore)) return;
  
  // Complement the reconstructed tracks
  if (GetRecoParam()->ComplementTracks()) {
    if (ComplementTracks(clusterStore)) RemoveIdenticalTracks();
  }
  
  // Improve the reconstructed tracks
  if (GetRecoParam()->ImproveTracks()) ImproveTracks();
  
  // Remove connected tracks
  if (GetRecoParam()->RemoveConnectedTracksInSt12()) RemoveConnectedTracks(kFALSE);
  else RemoveConnectedTracks(kTRUE);
  
  // Fill AliMUONTrack data members
  Finalize();
  
  // Add tracks to MUON data container 
  for (Int_t i=0; i<fNRecTracks; ++i) 
  {
    AliMUONTrack * track = (AliMUONTrack*) fRecTracksPtr->At(i);
    track->SetUniqueID(i+1);
    trackStore.Add(*track);
  }
}

  //__________________________________________________________________________
TClonesArray* AliMUONVTrackReconstructor::MakeSegmentsBetweenChambers(const AliMUONVClusterStore& clusterStore, Int_t ch1, Int_t ch2)
{
  /// To make the list of segments from the list of clusters in the 2 given chambers.
  /// Return a new TClonesArray of segments.
  /// It is the responsibility of the user to delete it afterward.
  AliDebug(1,Form("Enter MakeSegmentsBetweenChambers (1..) %d-%d", ch1+1, ch2+1));
  
  AliMUONVCluster *cluster1, *cluster2;
  AliMUONObjectPair *segment;
  Double_t nonBendingSlope = 0, bendingSlope = 0, impactParam = 0., bendingMomentum = 0.; // to avoid compilation warning
  
  // Create iterators to loop over clusters in both chambers
  TIter nextInCh1(clusterStore.CreateChamberIterator(ch1,ch1));
  TIter nextInCh2(clusterStore.CreateChamberIterator(ch2,ch2));
  
  // list of segments
  TClonesArray *segments = new TClonesArray("AliMUONObjectPair", 100);
  
  // Loop over clusters in the first chamber of the station
  while ( ( cluster1 = static_cast<AliMUONVCluster*>(nextInCh1()) ) ) {
    
    // reset cluster iterator of chamber 2
    nextInCh2.Reset();
    
    // Loop over clusters in the second chamber of the station
    while ( ( cluster2 = static_cast<AliMUONVCluster*>(nextInCh2()) ) ) {
      
      // non bending slope
      nonBendingSlope = (cluster1->GetX() - cluster2->GetX()) / (cluster1->GetZ() - cluster2->GetZ());
      
      // check if non bending slope is within tolerances
      if (TMath::Abs(nonBendingSlope) > GetRecoParam()->GetMaxNonBendingSlope()) continue;
      
      // bending slope
      bendingSlope = (cluster1->GetY() - cluster2->GetY()) / (cluster1->GetZ() - cluster2->GetZ());
      
      // check the bending momentum of the bending slope depending if the field is ON or OFF
      if (AliMUONTrackExtrap::IsFieldON()) {
	
	// impact parameter
	impactParam = cluster1->GetY() - cluster1->GetZ() * bendingSlope;
	
	// absolute value of bending momentum
	bendingMomentum = TMath::Abs(AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(impactParam));
	
	// check if bending momentum is within tolerances
	if (bendingMomentum < GetRecoParam()->GetMinBendingMomentum() ||
	    bendingMomentum > GetRecoParam()->GetMaxBendingMomentum()) continue;
	
      } else {
	
	// check if non bending slope is within tolerances
	if (TMath::Abs(bendingSlope) > GetRecoParam()->GetMaxBendingSlope()) continue;
      
      }
      // make new segment
      segment = new ((*segments)[segments->GetLast()+1]) AliMUONObjectPair(cluster1, cluster2, kFALSE, kFALSE);
      
      // Printout for debug
      if (AliLog::GetGlobalDebugLevel() > 1) {
	cout << "segmentIndex(0...): " << segments->GetLast() << endl;
	segment->Dump();
	cout << "Cluster in first chamber" << endl;
	cluster1->Print();
	cout << "Cluster in second chamber" << endl;
	cluster2->Print();
      }
      
    }
    
  }
  
  // Printout for debug
  AliDebug(1,Form("chambers%d-%d: NSegments =  %d ", ch1+1, ch2+1, segments->GetEntriesFast()));
  
  return segments;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::RemoveUsedSegments(TClonesArray& segments)
{
  /// To remove pairs of clusters already attached to a track
  AliDebug(1,"Enter RemoveUsedSegments");
  Int_t nSegments = segments.GetEntriesFast();
  Int_t nTracks = fRecTracksPtr->GetEntriesFast();
  AliMUONObjectPair *segment;
  AliMUONTrack *track;
  AliMUONVCluster *cluster, *cluster1, *cluster2;
  Bool_t foundCluster1, foundCluster2, removeSegment;
  
  // Loop over segments
  for (Int_t iSegment=0; iSegment<nSegments; iSegment++) {
    segment = (AliMUONObjectPair*) segments.UncheckedAt(iSegment);
    
    cluster1 = (AliMUONVCluster*) segment->First();
    cluster2 = (AliMUONVCluster*) segment->Second();
    removeSegment = kFALSE;
    
    // Loop over tracks
    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      track = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iTrack);
      
      // skip empty slot
      if (!track) continue;
      
      foundCluster1 = kFALSE;
      foundCluster2 = kFALSE;
      
      // Loop over clusters
      Int_t nClusters = track->GetNClusters();
      for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
        cluster = ((AliMUONTrackParam*) track->GetTrackParamAtCluster()->UncheckedAt(iCluster))->GetClusterPtr();
	
	// check if both clusters are in that track
	if (cluster == cluster1) foundCluster1 = kTRUE;
	else if (cluster == cluster2) foundCluster2 = kTRUE;
	
	if (foundCluster1 && foundCluster2) {
	  removeSegment = kTRUE;
	  break;
	}
	
      }
      
      if (removeSegment) break;
      
    }
    
    if (removeSegment) segments.RemoveAt(iSegment);
      
  }
  
  segments.Compress();
  
  // Printout for debug
  AliDebug(1,Form("NSegments =  %d ", segments.GetEntriesFast()));
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::RemoveIdenticalTracks()
{
  /// To remove identical tracks:
  /// Tracks are considered identical if they have all their clusters in common.
  /// One keeps the track with the larger number of clusters if need be
  AliMUONTrack *track1, *track2;
  Int_t nTracks = fRecTracksPtr->GetEntriesFast();
  Int_t clustersInCommon, nClusters1, nClusters2;
  // Loop over first track of the pair
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
    track1 = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iTrack1);
    // skip empty slot
    if (!track1) continue;
    nClusters1 = track1->GetNClusters();
    // Loop over second track of the pair
    for (Int_t iTrack2 = iTrack1+1; iTrack2 < nTracks; iTrack2++) {
      track2 = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iTrack2);
      // skip empty slot
      if (!track2) continue;
      nClusters2 = track2->GetNClusters();
      // number of clusters in common between two tracks
      clustersInCommon = track1->ClustersInCommon(track2);
      // check for identical tracks
      if ((clustersInCommon == nClusters1) || (clustersInCommon == nClusters2)) {
        // decide which track to remove
        if (nClusters2 > nClusters1) {
	  // remove track1 and continue the first loop with the track next to track1
          fRecTracksPtr->RemoveAt(iTrack1);
	  fNRecTracks--;
	  break;
	} else {
	  // remove track2 and continue the second loop with the track next to track2
	  fRecTracksPtr->RemoveAt(iTrack2);
	  fNRecTracks--;
        }
      }
    } // track2
  } // track1
  fRecTracksPtr->Compress(); // this is essential to retrieve the TClonesArray afterwards
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::RemoveDoubleTracks()
{
  /// To remove double tracks:
  /// Tracks are considered identical if more than half of the clusters of the track
  /// which has the smaller number of clusters are in common with the other track.
  /// Among two identical tracks, one keeps the track with the larger number of clusters
  /// or, if these numbers are equal, the track with the minimum chi2.
  AliMUONTrack *track1, *track2;
  Int_t nTracks = fRecTracksPtr->GetEntriesFast();
  Int_t clustersInCommon2, nClusters1, nClusters2;
  // Loop over first track of the pair
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
    track1 = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iTrack1);
    // skip empty slot
    if (!track1) continue;
    nClusters1 = track1->GetNClusters();
    // Loop over second track of the pair
    for (Int_t iTrack2 = iTrack1+1; iTrack2 < nTracks; iTrack2++) {
      track2 = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iTrack2);
      // skip empty slot
      if (!track2) continue;
      nClusters2 = track2->GetNClusters();
      // number of clusters in common between two tracks
      clustersInCommon2 = 2 * track1->ClustersInCommon(track2);
      // check for identical tracks
      if (clustersInCommon2 > nClusters1 || clustersInCommon2 > nClusters2) {
        // decide which track to remove
        if ((nClusters1 > nClusters2) || ((nClusters1 == nClusters2) && (track1->GetGlobalChi2() <= track2->GetGlobalChi2()))) {
	  // remove track2 and continue the second loop with the track next to track2
	  fRecTracksPtr->RemoveAt(iTrack2);
	  fNRecTracks--;
        } else {
	  // else remove track1 and continue the first loop with the track next to track1
          fRecTracksPtr->RemoveAt(iTrack1);
	  fNRecTracks--;
	  break;
        }
      }
    } // track2
  } // track1
  fRecTracksPtr->Compress(); // this is essential to retrieve the TClonesArray afterwards
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::RemoveConnectedTracks(Bool_t inSt345)
{
  /// To remove double tracks:
  /// Tracks are considered identical if they share 1 cluster or more.
  /// If inSt345=kTRUE only stations 3, 4 and 5 are considered.
  /// Among two identical tracks, one keeps the track with the larger number of clusters
  /// or, if these numbers are equal, the track with the minimum chi2.
  AliMUONTrack *track1, *track2, *trackToRemove;
  Int_t clustersInCommon, nClusters1, nClusters2;
  Bool_t removedTrack1;
  // Loop over first track of the pair
  track1 = (AliMUONTrack*) fRecTracksPtr->First();
  while (track1) {
    removedTrack1 = kFALSE;
    nClusters1 = track1->GetNClusters();
    // Loop over second track of the pair
    track2 = (AliMUONTrack*) fRecTracksPtr->After(track1);
    while (track2) {
      nClusters2 = track2->GetNClusters();
      // number of clusters in common between two tracks
      if (inSt345) clustersInCommon = track1->ClustersInCommonInSt345(track2);
      else clustersInCommon = track1->ClustersInCommon(track2);
      // check for identical tracks
      if (clustersInCommon > 0) {
        // decide which track to remove
        if ((nClusters1 > nClusters2) || ((nClusters1 == nClusters2) && (track1->GetGlobalChi2() <= track2->GetGlobalChi2()))) {
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
void AliMUONVTrackReconstructor::AskForNewClustersInChamber(const AliMUONTrackParam &trackParam,
							    AliMUONVClusterStore& clusterStore, Int_t chamber)
{
  /// Ask the clustering to reconstruct new clusters around the track candidate position
  
  // check if the current chamber is useable
  if (!fClusterServer || !GetRecoParam()->UseChamber(chamber)) return;
  
  // maximum distance between the center of the chamber and a detection element
  // (accounting for the inclination of the chamber)
  static const Double_t kMaxDZ = 15.; // 15 cm
  
  // extrapolate track parameters to the chamber
  AliMUONTrackParam extrapTrackParam(trackParam);
  AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParam, AliMUONConstants::DefaultChamberZ(chamber));
  
  // build the searching area using the track resolution and the maximum-distance-to-track value
  const TMatrixD& kParamCov = extrapTrackParam.GetCovariances();
  Double_t errX2 = kParamCov(0,0) + kMaxDZ * kMaxDZ * kParamCov(1,1) + 2. * kMaxDZ * TMath::Abs(kParamCov(0,1));
  Double_t errY2 = kParamCov(2,2) + kMaxDZ * kMaxDZ * kParamCov(3,3) + 2. * kMaxDZ * TMath::Abs(kParamCov(2,3));
  Double_t dX = TMath::Abs(trackParam.GetNonBendingSlope()) * kMaxDZ +
		GetRecoParam()->GetMaxNonBendingDistanceToTrack() +
		GetRecoParam()->GetSigmaCutForTracking() * TMath::Sqrt(errX2);
  Double_t dY = TMath::Abs(trackParam.GetBendingSlope()) * kMaxDZ +
		GetRecoParam()->GetMaxBendingDistanceToTrack() +
		GetRecoParam()->GetSigmaCutForTracking() * TMath::Sqrt(errY2);
  AliMpArea area(extrapTrackParam.GetNonBendingCoor(), 
                 extrapTrackParam.GetBendingCoor(),
                 dX, dY);
  
  // ask to cluterize in the given area of the given chamber
  fClusterServer->Clusterize(chamber, clusterStore, area, GetRecoParam());
  
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::AskForNewClustersInStation(const AliMUONTrackParam &trackParam,
							    AliMUONVClusterStore& clusterStore, Int_t station)
{
  /// Ask the clustering to reconstruct new clusters around the track candidate position
  /// in the 2 chambers of the given station
  AskForNewClustersInChamber(trackParam, clusterStore, 2*station+1);
  AskForNewClustersInChamber(trackParam, clusterStore, 2*station);
}

  //__________________________________________________________________________
Double_t AliMUONVTrackReconstructor::TryOneCluster(const AliMUONTrackParam &trackParam, AliMUONVCluster* cluster,
						   AliMUONTrackParam &trackParamAtCluster, Bool_t updatePropagator)
{
/// Test the compatibility between the track and the cluster (using trackParam's covariance matrix):
/// return the corresponding Chi2
/// return trackParamAtCluster
  
  // extrapolate track parameters and covariances at the z position of the tested cluster
  trackParamAtCluster = trackParam;
  AliMUONTrackExtrap::ExtrapToZCov(&trackParamAtCluster, cluster->GetZ(), updatePropagator);
  
  // set pointer to cluster into trackParamAtCluster
  trackParamAtCluster.SetClusterPtr(cluster);
  
  // Set differences between trackParam and cluster in the bending and non bending directions
  Double_t dX = cluster->GetX() - trackParamAtCluster.GetNonBendingCoor();
  Double_t dY = cluster->GetY() - trackParamAtCluster.GetBendingCoor();
  
  // Calculate errors and covariances
  const TMatrixD& kParamCov = trackParamAtCluster.GetCovariances();
  Double_t sigmaX2 = kParamCov(0,0) + cluster->GetErrX2();
  Double_t sigmaY2 = kParamCov(2,2) + cluster->GetErrY2();
  
  // Compute chi2
  return dX * dX / sigmaX2 + dY * dY / sigmaY2;
  
}

  //__________________________________________________________________________
Bool_t AliMUONVTrackReconstructor::TryOneClusterFast(const AliMUONTrackParam &trackParam, const AliMUONVCluster* cluster)
{
/// Test the compatibility between the track and the cluster
/// given the track resolution + the maximum-distance-to-track value
/// and assuming linear propagation of the track:
/// return kTRUE if they are compatibles
  
  Double_t dZ = cluster->GetZ() - trackParam.GetZ();
  Double_t dX = cluster->GetX() - (trackParam.GetNonBendingCoor() + trackParam.GetNonBendingSlope() * dZ);
  Double_t dY = cluster->GetY() - (trackParam.GetBendingCoor() + trackParam.GetBendingSlope() * dZ);
  const TMatrixD& kParamCov = trackParam.GetCovariances();
  Double_t errX2 = kParamCov(0,0) + dZ * dZ * kParamCov(1,1) + 2. * dZ * kParamCov(0,1);
  Double_t errY2 = kParamCov(2,2) + dZ * dZ * kParamCov(3,3) + 2. * dZ * kParamCov(2,3);

  Double_t dXmax = GetRecoParam()->GetSigmaCutForTracking() * TMath::Sqrt(errX2) +
                   GetRecoParam()->GetMaxNonBendingDistanceToTrack();
  Double_t dYmax = GetRecoParam()->GetSigmaCutForTracking() * TMath::Sqrt(errY2) +
		   GetRecoParam()->GetMaxBendingDistanceToTrack();
  
  if (TMath::Abs(dX) > dXmax || TMath::Abs(dY) > dYmax) return kFALSE;
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Double_t AliMUONVTrackReconstructor::TryTwoClustersFast(const AliMUONTrackParam &trackParamAtCluster1, AliMUONVCluster* cluster2,
						        AliMUONTrackParam &trackParamAtCluster2)
{
/// Test the compatibility between the track and the 2 clusters together (using trackParam's covariance matrix)
/// assuming linear propagation between the two clusters:
/// return the corresponding Chi2 accounting for covariances between the 2 clusters
/// return trackParamAtCluster2
  
  // extrapolate linearly track parameters and covariances at the z position of the second cluster
  trackParamAtCluster2 = trackParamAtCluster1;
  AliMUONTrackExtrap::LinearExtrapToZ(&trackParamAtCluster2, cluster2->GetZ());
  
  // set pointer to cluster2 into trackParamAtCluster2
  trackParamAtCluster2.SetClusterPtr(cluster2);
  
  // Set differences between track and clusters in the bending and non bending directions
  AliMUONVCluster* cluster1 = trackParamAtCluster1.GetClusterPtr();
  Double_t dX1 = cluster1->GetX() - trackParamAtCluster1.GetNonBendingCoor();
  Double_t dX2 = cluster2->GetX() - trackParamAtCluster2.GetNonBendingCoor();
  Double_t dY1 = cluster1->GetY() - trackParamAtCluster1.GetBendingCoor();
  Double_t dY2 = cluster2->GetY() - trackParamAtCluster2.GetBendingCoor();
  
  // Calculate errors and covariances
  const TMatrixD& kParamCov1 = trackParamAtCluster1.GetCovariances();
  const TMatrixD& kParamCov2 = trackParamAtCluster2.GetCovariances();
  Double_t dZ = trackParamAtCluster2.GetZ() - trackParamAtCluster1.GetZ();
  Double_t sigma2X1 = kParamCov1(0,0) + cluster1->GetErrX2();
  Double_t sigma2X2 = kParamCov2(0,0) + cluster2->GetErrX2();
  Double_t covX1X2  = kParamCov1(0,0) + dZ * kParamCov1(0,1);
  Double_t sigma2Y1 = kParamCov1(2,2) + cluster1->GetErrY2();
  Double_t sigma2Y2 = kParamCov2(2,2) + cluster2->GetErrY2();
  Double_t covY1Y2  = kParamCov1(2,2) + dZ * kParamCov1(2,3);
  
  // Compute chi2
  Double_t detX = sigma2X1 * sigma2X2 - covX1X2 * covX1X2;
  Double_t detY = sigma2Y1 * sigma2Y2 - covY1Y2 * covY1Y2;
  if (detX == 0. || detY == 0.) return 1.e10;
  return   (dX1 * dX1 * sigma2X2 + dX2 * dX2 * sigma2X1 - 2. * dX1 * dX2 * covX1X2) / detX
  	 + (dY1 * dY1 * sigma2Y2 + dY2 * dY2 * sigma2Y1 - 2. * dY1 * dY2 * covY1Y2) / detY;
  
}

  //__________________________________________________________________________
Bool_t AliMUONVTrackReconstructor::FollowLinearTrackInChamber(AliMUONTrack &trackCandidate, const AliMUONVClusterStore& clusterStore,
							      Int_t nextChamber)
{
  /// Follow trackCandidate in chamber(0..) nextChamber assuming linear propagation, and search for compatible cluster(s)
  /// Keep all possibilities or only the best one(s) according to the flag fgkTrackAllTracks:
  /// kTRUE:  duplicate "trackCandidate" if there are several possibilities and add the new tracks at the end of
  ///         fRecTracksPtr to avoid conficts with other track candidates at this current stage of the tracking procedure.
  ///         Remove the obsolete "trackCandidate" at the end.
  /// kFALSE: add only the best cluster(s) to the "trackCandidate". Try to add a couple of clusters in priority.
  /// return kTRUE if new cluster(s) have been found (otherwise return kFALSE)
  AliDebug(1,Form("Enter FollowLinearTrackInChamber(1..) %d", nextChamber+1));
  
  Double_t chi2WithOneCluster = 1.e10;
  Double_t maxChi2WithOneCluster = 2. * GetRecoParam()->GetSigmaCutForTracking() *
					GetRecoParam()->GetSigmaCutForTracking(); // 2 because 2 quantities in chi2
  Double_t bestChi2WithOneCluster = maxChi2WithOneCluster;
  Bool_t foundOneCluster = kFALSE;
  AliMUONTrack *newTrack = 0x0;
  AliMUONVCluster *cluster;
  AliMUONTrackParam trackParam;
  AliMUONTrackParam extrapTrackParamAtCluster;
  AliMUONTrackParam bestTrackParamAtCluster;
  
  // Get track parameters according to the propagation direction
  if (nextChamber > 7) trackParam = *(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->Last();
  else trackParam = *(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->First();
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    cout<<endl<<"Track parameters and covariances at first cluster:"<<endl;
    trackParam.GetParameters().Print();
    trackParam.GetCovariances().Print();
  }
  
  // Add MCS effect
  AliMUONTrackExtrap::AddMCSEffect(&trackParam,AliMUONConstants::ChamberThicknessInX0(),1.);
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "FollowLinearTrackInChamber: look for cluster in chamber(1..): " << nextChamber+1 << endl;
  }
  
  // Create iterators to loop over clusters in chamber
  TIter next(clusterStore.CreateChamberIterator(nextChamber,nextChamber));
  
  // look for candidates in chamber
  while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) ) {
    
    // try to add the current cluster fast
    if (!TryOneClusterFast(trackParam, cluster)) continue;
    
    // try to add the current cluster accuratly
    extrapTrackParamAtCluster = trackParam;
    AliMUONTrackExtrap::LinearExtrapToZ(&extrapTrackParamAtCluster, cluster->GetZ());
    chi2WithOneCluster = TryOneCluster(extrapTrackParamAtCluster, cluster, extrapTrackParamAtCluster);
    
    // if good chi2 then consider to add cluster
    if (chi2WithOneCluster < maxChi2WithOneCluster) {
      foundOneCluster = kTRUE;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	cout << "FollowLinearTrackInChamber: found one cluster in chamber(1..): " << nextChamber+1
	<< " (Chi2 = " << chi2WithOneCluster << ")" << endl;
	cluster->Print();
      }
      
      if (GetRecoParam()->TrackAllTracks()) {
	// copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
	newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	if (GetRecoParam()->RequestStation(nextChamber/2))
	  extrapTrackParamAtCluster.SetRemovable(kFALSE);
	else extrapTrackParamAtCluster.SetRemovable(kTRUE);
	newTrack->AddTrackParamAtCluster(extrapTrackParamAtCluster,*cluster);
	newTrack->SetGlobalChi2(trackCandidate.GetGlobalChi2()+chi2WithOneCluster);
	fNRecTracks++;
	
	// Printout for debuging
	if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	  cout << "FollowLinearTrackInChamber: added one cluster in chamber(1..): " << nextChamber+1 << endl;
	  if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	}
	
      } else if (chi2WithOneCluster < bestChi2WithOneCluster) {
	// keep track of the best cluster
	bestChi2WithOneCluster = chi2WithOneCluster;
	bestTrackParamAtCluster = extrapTrackParamAtCluster;
      }
      
    }
    
  }
  
  // fill out the best track if required else clean up the fRecTracksPtr array
  if (!GetRecoParam()->TrackAllTracks()) {
    if (foundOneCluster) {
      if (GetRecoParam()->RequestStation(nextChamber/2))
	bestTrackParamAtCluster.SetRemovable(kFALSE);
      else bestTrackParamAtCluster.SetRemovable(kTRUE);
      trackCandidate.AddTrackParamAtCluster(bestTrackParamAtCluster,*(bestTrackParamAtCluster.GetClusterPtr()));
      trackCandidate.SetGlobalChi2(trackCandidate.GetGlobalChi2()+bestChi2WithOneCluster);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowLinearTrackInChamber: added the best cluster in chamber(1..): " << bestTrackParamAtCluster.GetClusterPtr()->GetChamberId()+1 << endl;
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
Bool_t AliMUONVTrackReconstructor::FollowLinearTrackInStation(AliMUONTrack &trackCandidate, const AliMUONVClusterStore& clusterStore,
							      Int_t nextStation)
{
  /// Follow trackCandidate in station(0..) nextStation assuming linear propagation, and search for compatible cluster(s)
  /// Keep all possibilities or only the best one(s) according to the flag fgkTrackAllTracks:
  /// kTRUE:  duplicate "trackCandidate" if there are several possibilities and add the new tracks at the end of
  ///         fRecTracksPtr to avoid conficts with other track candidates at this current stage of the tracking procedure.
  ///         Remove the obsolete "trackCandidate" at the end.
  /// kFALSE: add only the best cluster(s) to the "trackCandidate". Try to add a couple of clusters in priority.
  /// return kTRUE if new cluster(s) have been found (otherwise return kFALSE)
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
  
  Double_t chi2WithOneCluster = 1.e10;
  Double_t chi2WithTwoClusters = 1.e10;
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
  AliMUONTrackParam trackParam;
  AliMUONTrackParam extrapTrackParamAtCluster1;
  AliMUONTrackParam extrapTrackParamAtCluster2;
  AliMUONTrackParam bestTrackParamAtCluster1;
  AliMUONTrackParam bestTrackParamAtCluster2;
  
  Int_t nClusters = clusterStore.GetSize();
  Bool_t *clusterCh1Used = new Bool_t[nClusters];
  for (Int_t i = 0; i < nClusters; i++) clusterCh1Used[i] = kFALSE;
  Int_t iCluster1;
  
  // Get track parameters according to the propagation direction
  if (nextStation==4) trackParam = *(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->Last();
  else trackParam = *(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->First();
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    cout<<endl<<"Track parameters and covariances at first cluster:"<<endl;
    trackParam.GetParameters().Print();
    trackParam.GetCovariances().Print();
  }
  
  // Add MCS effect
  AliMUONTrackExtrap::AddMCSEffect(&trackParam,AliMUONConstants::ChamberThicknessInX0(),1.);
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "FollowLinearTrackInStation: look for clusters in chamber(1..): " << ch2+1 << endl;
  }
  
  // Create iterators to loop over clusters in both chambers
  TIter nextInCh1(clusterStore.CreateChamberIterator(ch1,ch1));
  TIter nextInCh2(clusterStore.CreateChamberIterator(ch2,ch2));
  
  // look for candidates in chamber 2
  while ( ( clusterCh2 = static_cast<AliMUONVCluster*>(nextInCh2()) ) ) {
    
    // try to add the current cluster fast
    if (!TryOneClusterFast(trackParam, clusterCh2)) continue;
    
    // try to add the current cluster accuratly
    extrapTrackParamAtCluster2 = trackParam;
    AliMUONTrackExtrap::LinearExtrapToZ(&extrapTrackParamAtCluster2, clusterCh2->GetZ());
    chi2WithOneCluster = TryOneCluster(extrapTrackParamAtCluster2, clusterCh2, extrapTrackParamAtCluster2);
    
    // if good chi2 then try to attach a cluster in the other chamber too
    if (chi2WithOneCluster < maxChi2WithOneCluster) {
      Bool_t foundSecondCluster = kFALSE;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowLinearTrackInStation: found one cluster in chamber(1..): " << ch2+1
	     << " (Chi2 = " << chi2WithOneCluster << ")" << endl;
	clusterCh2->Print();
        cout << "                      look for second clusters in chamber(1..): " << ch1+1 << " ..." << endl;
      }
      
      // add MCS effect
      AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCluster2,AliMUONConstants::ChamberThicknessInX0(),1.);
      
      // reset cluster iterator of chamber 1
      nextInCh1.Reset();
      iCluster1 = -1;
      
      // look for second candidates in chamber 1
      while ( ( clusterCh1 = static_cast<AliMUONVCluster*>(nextInCh1()) ) ) {
        iCluster1++;
	
        // try to add the current cluster fast
        if (!TryOneClusterFast(extrapTrackParamAtCluster2, clusterCh1)) continue;
	
	// try to add the current cluster in addition to the one found in the previous chamber
	chi2WithTwoClusters = TryTwoClustersFast(extrapTrackParamAtCluster2, clusterCh1, extrapTrackParamAtCluster1);
        
	// if good chi2 then consider to add the 2 clusters to the "trackCandidate"
	if (chi2WithTwoClusters < maxChi2WithTwoClusters) {
	  foundSecondCluster = kTRUE;
	  foundTwoClusters = kTRUE;
          
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowLinearTrackInStation: found one cluster in chamber(1..): " << ch1+1
		 << " (Chi2 = " << chi2WithTwoClusters << ")" << endl;
	    clusterCh1->Print();
	  }
	  
	  if (GetRecoParam()->TrackAllTracks()) {
	    // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new clusters
            newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	    extrapTrackParamAtCluster1.SetRemovable(kTRUE);
	    newTrack->AddTrackParamAtCluster(extrapTrackParamAtCluster1,*clusterCh1);
	    extrapTrackParamAtCluster2.SetRemovable(kTRUE);
	    newTrack->AddTrackParamAtCluster(extrapTrackParamAtCluster2,*clusterCh2);
	    newTrack->SetGlobalChi2(newTrack->GetGlobalChi2()+chi2WithTwoClusters);
	    fNRecTracks++;
	    
	    // Tag clusterCh1 as used
	    clusterCh1Used[iCluster1] = kTRUE;
	    
	    // Printout for debuging
	    if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	      cout << "FollowLinearTrackInStation: added two clusters in station(1..): " << nextStation+1 << endl;
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
      
      // if no cluster found in chamber1 then consider to add clusterCh2 only
      if (!foundSecondCluster) {
	foundOneCluster = kTRUE;
        
	if (GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
          newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  if (GetRecoParam()->RequestStation(nextStation))
	    extrapTrackParamAtCluster2.SetRemovable(kFALSE);
	  else extrapTrackParamAtCluster2.SetRemovable(kTRUE);
	  newTrack->AddTrackParamAtCluster(extrapTrackParamAtCluster2,*clusterCh2);
	  newTrack->SetGlobalChi2(newTrack->GetGlobalChi2()+chi2WithOneCluster);
	  fNRecTracks++;
	  
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowLinearTrackInStation: added one cluster in chamber(1..): " << ch2+1 << endl;
	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	  }
	  
	} else if (!foundTwoClusters && chi2WithOneCluster < bestChi2WithOneCluster) {
	  // keep track of the best cluster except if a couple of clusters has already been found
	  bestChi2WithOneCluster = chi2WithOneCluster;
	  bestTrackParamAtCluster1 = extrapTrackParamAtCluster2;
        }
	
      }
      
    }
    
  }
  
  // look for candidates in chamber 1 not already attached to a track
  // if we want to keep all possible tracks or if no good couple of clusters has been found
  if (GetRecoParam()->TrackAllTracks() || !foundTwoClusters) {
    
    // Printout for debuging
    if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
      cout << "FollowLinearTrackInStation: look for single cluster in chamber(1..): " << ch1+1 << endl;
    }
    
    //Extrapolate trackCandidate to chamber "ch2"
    AliMUONTrackExtrap::LinearExtrapToZ(&trackParam, AliMUONConstants::DefaultChamberZ(ch2));
    
    // add MCS effect for next step
    AliMUONTrackExtrap::AddMCSEffect(&trackParam,AliMUONConstants::ChamberThicknessInX0(),1.);
    
    // reset cluster iterator of chamber 1
    nextInCh1.Reset();
    iCluster1 = -1;
    
    // look for second candidates in chamber 1
    while ( ( clusterCh1 = static_cast<AliMUONVCluster*>(nextInCh1()) ) ) {
      iCluster1++;
      
      if (clusterCh1Used[iCluster1]) continue; // Skip clusters already used
      
      // try to add the current cluster fast
      if (!TryOneClusterFast(trackParam, clusterCh1)) continue;
      
      // try to add the current cluster accuratly
      extrapTrackParamAtCluster1 = trackParam;
      AliMUONTrackExtrap::LinearExtrapToZ(&extrapTrackParamAtCluster1, clusterCh1->GetZ());
      chi2WithOneCluster = TryOneCluster(extrapTrackParamAtCluster1, clusterCh1, extrapTrackParamAtCluster1);
      
      // if good chi2 then consider to add clusterCh1
      // We do not try to attach a cluster in the other chamber too since it has already been done above
      if (chi2WithOneCluster < maxChi2WithOneCluster) {
	foundOneCluster = kTRUE;
  	
	// Printout for debuging
  	if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	  cout << "FollowLinearTrackInStation: found one cluster in chamber(1..): " << ch1+1
  	       << " (Chi2 = " << chi2WithOneCluster << ")" << endl;
	  clusterCh1->Print();
  	}
	
	if (GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
  	  newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  if (GetRecoParam()->RequestStation(nextStation))
	    extrapTrackParamAtCluster1.SetRemovable(kFALSE);
	  else extrapTrackParamAtCluster1.SetRemovable(kTRUE);
	  newTrack->AddTrackParamAtCluster(extrapTrackParamAtCluster1,*clusterCh1);
	  newTrack->SetGlobalChi2(newTrack->GetGlobalChi2()+chi2WithOneCluster);
	  fNRecTracks++;
  	  
	  // Printout for debuging
  	  if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	    cout << "FollowLinearTrackInStation: added one cluster in chamber(1..): " << ch1+1 << endl;
  	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
  	  }
	  
	} else if (chi2WithOneCluster < bestChi2WithOneCluster) {
	  // keep track of the best cluster except if a couple of clusters has already been found
	  bestChi2WithOneCluster = chi2WithOneCluster;
	  bestTrackParamAtCluster1 = extrapTrackParamAtCluster1;
  	}
	
      }
      
    }
    
  }
  
  // fill out the best track if required else clean up the fRecTracksPtr array
  if (!GetRecoParam()->TrackAllTracks()) {
    if (foundTwoClusters) {
      bestTrackParamAtCluster1.SetRemovable(kTRUE);
      trackCandidate.AddTrackParamAtCluster(bestTrackParamAtCluster1,*(bestTrackParamAtCluster1.GetClusterPtr()));
      bestTrackParamAtCluster2.SetRemovable(kTRUE);
      trackCandidate.AddTrackParamAtCluster(bestTrackParamAtCluster2,*(bestTrackParamAtCluster2.GetClusterPtr()));
      trackCandidate.SetGlobalChi2(trackCandidate.GetGlobalChi2()+bestChi2WithTwoClusters);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowLinearTrackInStation: added the two best clusters in station(1..): " << nextStation+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else if (foundOneCluster) {
      if (GetRecoParam()->RequestStation(nextStation))
	bestTrackParamAtCluster1.SetRemovable(kFALSE);
      else bestTrackParamAtCluster1.SetRemovable(kTRUE);
      trackCandidate.AddTrackParamAtCluster(bestTrackParamAtCluster1,*(bestTrackParamAtCluster1.GetClusterPtr()));
      trackCandidate.SetGlobalChi2(trackCandidate.GetGlobalChi2()+bestChi2WithOneCluster);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowLinearTrackInStation: added the best cluster in chamber(1..): " << bestTrackParamAtCluster1.GetClusterPtr()->GetChamberId()+1 << endl;
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
void AliMUONVTrackReconstructor::ImproveTracks()
{
  /// Improve tracks by removing clusters with local chi2 highter than the defined cut
  /// Recompute track parameters and covariances at the remaining clusters
  AliDebug(1,"Enter ImproveTracks");
  
  AliMUONTrack *track, *nextTrack;
  
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    
    // prepare next track in case the actual track is suppressed
    nextTrack = (AliMUONTrack*) fRecTracksPtr->After(track);
    
    ImproveTrack(*track);
    
    // remove track if improvement failed
    if (!track->IsImproved()) {
      fRecTracksPtr->Remove(track);
      fNRecTracks--;
    }
    
    track = nextTrack;
  }
  
  // compress the array in case of some tracks have been removed
  fRecTracksPtr->Compress();
  
}

//__________________________________________________________________________
void AliMUONVTrackReconstructor::Finalize()
{
  /// Recompute track parameters and covariances at each attached cluster from those at the first one
  /// Set the label pointing to the corresponding MC track
  
  AliMUONTrack *track;
  
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    
    FinalizeTrack(*track);
    
    track->FindMCLabel();
    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
    
  }
  
}

//__________________________________________________________________________
void AliMUONVTrackReconstructor::ValidateTracksWithTrigger(AliMUONVTrackStore& trackStore,
                                                           const AliMUONVTriggerTrackStore& triggerTrackStore,
                                                           const AliMUONVTriggerStore& triggerStore,
                                                           const AliMUONTrackHitPattern& trackHitPattern)
{
  /// Try to match track from tracking system with trigger track
  AliCodeTimerAuto("");

  trackHitPattern.ExecuteValidation(trackStore, triggerTrackStore, triggerStore);
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
    Bool_t xTrig=locTrg->IsTrigX();
    Bool_t yTrig=locTrg->IsTrigY();
    
    Int_t localBoardId = locTrg->LoCircuit();
    
    if (xTrig && yTrig) 
    { // make Trigger Track if trigger in X and Y
      
      Float_t y11 = circuit.GetY11Pos(localBoardId, locTrg->LoStripX()); 
      // need first to convert deviation to [0-30] 
      // (see AliMUONLocalTriggerBoard::LocalTrigger)
      Int_t deviation = locTrg->GetDeviation(); 
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

