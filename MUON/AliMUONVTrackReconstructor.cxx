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
/// tracking algorithms. They can be changed by using:
/// AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLow(High)FluxParam();
/// muonRecoParam->Set...(); // see methods in AliMUONRecoParam.h for details
/// AliMUONReconstructor::SetRecoParam(muonRecoParam);
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

  //__________________________________________________________________________
AliMUONVTrackReconstructor::AliMUONVTrackReconstructor()
  : TObject(),
    fRecTracksPtr(0x0),
    fNRecTracks(0)
{
  /// Constructor for class AliMUONVTrackReconstructor
  
  // Memory allocation for the TClonesArray of reconstructed tracks
  fRecTracksPtr = new TClonesArray("AliMUONTrack", 100);
  
  // set the magnetic field for track extrapolations
  const AliMagF* kField = AliTracker::GetFieldMap();
  if (!kField) AliFatal("No field available");
  AliMUONTrackExtrap::SetField(kField);
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
void AliMUONVTrackReconstructor::EventReconstruct(const AliMUONVClusterStore& clusterStore, AliMUONVTrackStore& trackStore)
{
  /// To reconstruct one event
  AliDebug(1,"");
  AliCodeTimerAuto("");
  
  // Reset array of tracks
  ResetTracks();
  
  // Look for candidates from at least 3 aligned points in stations(1..) 4 and 5
  MakeTrackCandidates(clusterStore);
  
  // Stop tracking if no candidate found
  if (fRecTracksPtr->GetEntriesFast() == 0) return;
  
  // Follow tracks in stations(1..) 3, 2 and 1
  FollowTracks(clusterStore);
  
  // Complement the reconstructed tracks
  if (AliMUONReconstructor::GetRecoParam()->ComplementTracks()) ComplementTracks(clusterStore);
  
  // Improve the reconstructed tracks
  if (AliMUONReconstructor::GetRecoParam()->ImproveTracks()) ImproveTracks();
  
  // Remove double tracks
  RemoveDoubleTracks();
  
  // Fill AliMUONTrack data members
  Finalize();
  
  // Add tracks to MUON data container 
  for (Int_t i=0; i<fNRecTracks; ++i) 
  {
    AliMUONTrack * track = (AliMUONTrack*) fRecTracksPtr->At(i);
    trackStore.Add(*track);
  }
}

  //__________________________________________________________________________
TClonesArray* AliMUONVTrackReconstructor::MakeSegmentsInStation(const AliMUONVClusterStore& clusterStore, Int_t station)
{
  /// To make the list of segments in station(0..) "Station" from the list of clusters to be reconstructed.
  /// Return a new TClonesArray of segments.
  /// It is the responsibility of the user to delete it afterward.
  AliDebug(1,Form("Enter MakeSegmentsPerStation (1..) %d",station+1));
  
  AliMUONVCluster *cluster1, *cluster2;
  AliMUONObjectPair *segment;
  Double_t bendingSlope = 0, impactParam = 0., bendingMomentum = 0.; // to avoid compilation warning
  Int_t ch1 = 2 * station;
  Int_t ch2 = ch1 + 1;
  
  // Create iterators to loop over clusters in both chambers
  TIter nextInCh1(clusterStore.CreateChamberIterator(ch1,ch1));
  TIter nextInCh2(clusterStore.CreateChamberIterator(ch2,ch2));
  
  // list of segments
  TClonesArray *segments = new TClonesArray("AliMUONObjectPair", 100);
  segments->SetOwner(kTRUE);
  
  // Loop over clusters in the first chamber of the station
  while ( ( cluster1 = static_cast<AliMUONVCluster*>(nextInCh1()) ) ) {
    
    // reset cluster iterator of chamber 2
    nextInCh2.Reset();
    
    // Loop over clusters in the second chamber of the station
    while ( ( cluster2 = static_cast<AliMUONVCluster*>(nextInCh2()) ) ) {
      
      // bending slope
      bendingSlope = (cluster1->GetY() - cluster2->GetY()) / (cluster1->GetZ() - cluster2->GetZ());
      
      // impact parameter
      impactParam = cluster1->GetY() - cluster1->GetZ() * bendingSlope;
     
      // absolute value of bending momentum
      bendingMomentum = TMath::Abs(AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(impactParam));
      
      // check for bending momentum within tolerances
      if ((bendingMomentum < AliMUONReconstructor::GetRecoParam()->GetMaxBendingMomentum()) &&
	  (bendingMomentum > AliMUONReconstructor::GetRecoParam()->GetMinBendingMomentum())) {
        
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
    
  }
  
  // Printout for debug
  AliDebug(1,Form("Station: %d  NSegments:  %d ", station+1, segments->GetEntriesFast()));
  
  return segments;
}

  //__________________________________________________________________________
void AliMUONVTrackReconstructor::RemoveIdenticalTracks()
{
  /// To remove identical tracks:
  /// Tracks are considered identical if they have all their clusters in common.
  /// One keeps the track with the larger number of clusters if need be
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
      clustersInCommon = track1->ClustersInCommon(track2);
      // check for identical tracks
      if ((clustersInCommon == nClusters1) || (clustersInCommon == nClusters2)) {
        // decide which track to remove
        if (nClusters2 > nClusters1) {
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
  /// Tracks are considered identical if more than half of the clusters of the track
  /// which has the smaller number of clusters are in common with the other track.
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
      clustersInCommon = track1->ClustersInCommon(track2);
      // check for identical tracks
      if (((nClusters1 < nClusters2) && (2 * clustersInCommon > nClusters1)) || (2 * clustersInCommon > nClusters2)) {
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
Bool_t AliMUONVTrackReconstructor::TryOneClusterFast(const AliMUONTrackParam &trackParam, AliMUONVCluster* cluster)
{
/// Test the compatibility between the track and the cluster within a wide fix window
/// assuming linear propagation of the track:
/// return kTRUE if they are compatibles
  
  Double_t dZ = cluster->GetZ() - trackParam.GetZ();
  Double_t dX = cluster->GetX() - (trackParam.GetNonBendingCoor() + trackParam.GetNonBendingSlope() * dZ);
  Double_t dY = cluster->GetY() - (trackParam.GetBendingCoor() + trackParam.GetBendingSlope() * dZ);
  
  if (TMath::Abs(dX) > AliMUONReconstructor::GetRecoParam()->GetMaxNonBendingDistanceToTrack() ||
      TMath::Abs(dY) > AliMUONReconstructor::GetRecoParam()->GetMaxBendingDistanceToTrack()) return kFALSE;
  
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
  Double_t maxChi2WithOneCluster = 2. * AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking() *
					AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking(); // 2 because 2 quantities in chi2
  Double_t maxChi2WithTwoClusters = 4. * AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking() *
					 AliMUONReconstructor::GetRecoParam()->GetSigmaCutForTracking(); // 4 because 4 quantities in chi2
  Double_t bestChi2WithOneCluster = maxChi2WithOneCluster;
  Double_t bestChi2WithTwoClusters = maxChi2WithTwoClusters;
  Bool_t foundOneCluster = kFALSE;
  Bool_t foundTwoClusters = kFALSE;
  AliMUONTrack *newTrack = 0x0;
  AliMUONVCluster *clusterCh1, *clusterCh2;
  AliMUONTrackParam extrapTrackParamAtCluster1;
  AliMUONTrackParam extrapTrackParamAtCluster2;
  AliMUONTrackParam bestTrackParamAtCluster1;
  AliMUONTrackParam bestTrackParamAtCluster2;
  
  Int_t nClusters = clusterStore.GetSize();
  Bool_t *clusterCh1Used = new Bool_t[nClusters];
  for (Int_t i = 0; i < nClusters; i++) clusterCh1Used[i] = kFALSE;
  Int_t iCluster1;
  
  // Get track parameters
  AliMUONTrackParam trackParam(*(AliMUONTrackParam*)trackCandidate.GetTrackParamAtCluster()->First());
  
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
	  }
	  
	  if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
	    // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new clusters
            newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	    extrapTrackParamAtCluster1.SetRemovable(kTRUE);
	    newTrack->AddTrackParamAtCluster(extrapTrackParamAtCluster1,*clusterCh1);
	    extrapTrackParamAtCluster2.SetRemovable(kTRUE);
	    newTrack->AddTrackParamAtCluster(extrapTrackParamAtCluster2,*clusterCh2);
	    newTrack->GetTrackParamAtCluster()->Sort();
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
        
	if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
          newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  extrapTrackParamAtCluster2.SetRemovable(kFALSE);
	  newTrack->AddTrackParamAtCluster(extrapTrackParamAtCluster2,*clusterCh2);
	  newTrack->GetTrackParamAtCluster()->Sort();
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
  if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks() || !foundTwoClusters) {
    
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
  	}
	
	if (AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new cluster
  	  newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  extrapTrackParamAtCluster1.SetRemovable(kFALSE);
	  newTrack->AddTrackParamAtCluster(extrapTrackParamAtCluster1,*clusterCh1);
	  newTrack->GetTrackParamAtCluster()->Sort();
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
  if (!AliMUONReconstructor::GetRecoParam()->TrackAllTracks()) {
    if (foundTwoClusters) {
      bestTrackParamAtCluster1.SetRemovable(kTRUE);
      trackCandidate.AddTrackParamAtCluster(bestTrackParamAtCluster1,*(bestTrackParamAtCluster1.GetClusterPtr()));
      bestTrackParamAtCluster2.SetRemovable(kTRUE);
      trackCandidate.AddTrackParamAtCluster(bestTrackParamAtCluster2,*(bestTrackParamAtCluster2.GetClusterPtr()));
      trackCandidate.GetTrackParamAtCluster()->Sort();
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONVTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowLinearTrackInStation: added the two best clusters in station(1..): " << nextStation+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else if (foundOneCluster) {
      bestTrackParamAtCluster1.SetRemovable(kFALSE);
      trackCandidate.AddTrackParamAtCluster(bestTrackParamAtCluster1,*(bestTrackParamAtCluster1.GetClusterPtr()));
      trackCandidate.GetTrackParamAtCluster()->Sort();
      
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
void AliMUONVTrackReconstructor::ValidateTracksWithTrigger(AliMUONVTrackStore& trackStore,
                                                           const AliMUONVTriggerTrackStore& triggerTrackStore,
                                                           const AliMUONVTriggerStore& triggerStore,
                                                           const AliMUONTrackHitPattern& trackHitPattern)
{
  /// Try to match track from tracking system with trigger track
  AliCodeTimerAuto("");

  const Double_t kDeltaZ = TMath::Abs(AliMUONConstants::DefaultChamberZ(12) - AliMUONConstants::DefaultChamberZ(10));

  //static const Double_t kDistSigma[3]={1,1,0.02}; // sigma of distributions (trigger-track) X,Y,slopeY
  // sigma of distributions (trigger-track) X,Y,slopeY
  const Double_t kDistSigma[3]={AliMUONConstants::TriggerNonBendingReso(),
				AliMUONConstants::TriggerBendingReso(),
				1.414 * AliMUONConstants::TriggerBendingReso()/kDeltaZ};

  const Double_t kTrigNonBendReso = AliMUONConstants::TriggerNonBendingReso();
  const Double_t kTrigBendReso = AliMUONConstants::TriggerBendingReso();
  const Double_t kTrigSlopeBendReso = 1.414 * AliMUONConstants::TriggerBendingReso()/kDeltaZ;
  const Double_t kTrigCovSlopeBend = - kTrigBendReso * kTrigBendReso / kDeltaZ;

  // Covariance matrix 3x3 (X,Y,slopeY) for trigger tracks
  TMatrixD trigCov(3,3);
  trigCov.Zero();
  trigCov(0,0) = kTrigNonBendReso * kTrigNonBendReso;
  trigCov(1,1) = kTrigBendReso * kTrigBendReso;
  trigCov(2,2) = kTrigSlopeBendReso * kTrigSlopeBendReso;
  trigCov(1,2) = trigCov(2,1) = kTrigCovSlopeBend;

  Int_t matchTrigger;
  Int_t loTrgNum(-1);
  Double_t distTriggerTrack[3], sigma2[3];
  Double_t xTrack, yTrack, ySlopeTrack, chi2MatchTrigger, minChi2MatchTrigger, chi2;

  TIter itTrack(trackStore.CreateIterator());
  AliMUONTrack* track;

  const Float_t kZFilterOut = AliMUONConstants::MuonFilterZEnd();
  const Float_t kFilterThickness = TMath::Abs(kZFilterOut-AliMUONConstants::MuonFilterZBeg()); // cm
  const Int_t kFirstTrigCh = AliMUONConstants::NTrackingCh();
  
  while ( ( track = static_cast<AliMUONTrack*>(itTrack()) ) )
  {
    matchTrigger = 0;
    chi2MatchTrigger = 0.;
    loTrgNum = -1;
    Int_t doubleMatch=-1; // Check if track matches 2 trigger tracks
    Double_t doubleChi2 = -1.;
    
    AliMUONTrackParam trackParam(*((AliMUONTrackParam*) (track->GetTrackParamAtCluster()->Last())));

    AliMUONTrackExtrap::ExtrapToZCov(&trackParam, kZFilterOut); // Extrap to muon filter end
    AliMUONTrackExtrap::AddMCSEffect(&trackParam, kFilterThickness, AliMUONConstants::MuonFilterX0()); // Add MCS effects
    AliMUONTrackExtrap::ExtrapToZCov(&trackParam, AliMUONConstants::DefaultChamberZ(kFirstTrigCh)); // extrap to 1st trigger chamber

    const TMatrixD& kParamCov = trackParam.GetCovariances();
    
    xTrack = trackParam.GetNonBendingCoor();
    yTrack = trackParam.GetBendingCoor();
    ySlopeTrack = trackParam.GetBendingSlope();

    // Covariance matrix 3x3 (X,Y,slopeY) for tracker tracks
    TMatrixD trackCov(3,3);
    trackCov.Zero();
    trackCov(0,0) = kParamCov(0,0);
    trackCov(1,1) = kParamCov(2,2);
    trackCov(2,2) = kParamCov(3,3);
    trackCov(1,2) = kParamCov(2,3);
    trackCov(2,1) = kParamCov(3,2);

    TMatrixD sumCov(trackCov,TMatrixD::kPlus,trigCov);

    Bool_t isCovOK = kTRUE;

    if (sumCov.Determinant() != 0) {
      sumCov.Invert();
    } else {
      AliWarning(" Determinant = 0");
      isCovOK = kFALSE;
      sigma2[0] = kParamCov(0,0);
      sigma2[1] = kParamCov(2,2);
      sigma2[2] = kParamCov(3,3);
      for (Int_t iVar = 0; iVar < 3; iVar++) sigma2[iVar] += kDistSigma[iVar] * kDistSigma[iVar];
    }

    minChi2MatchTrigger = 999.;

    AliMUONTriggerTrack *triggerTrack;
    TIter itTriggerTrack(triggerTrackStore.CreateIterator());
    while ( ( triggerTrack = static_cast<AliMUONTriggerTrack*>(itTriggerTrack() ) ) )
    {
      distTriggerTrack[0] = triggerTrack->GetX11()-xTrack;
      distTriggerTrack[1] = triggerTrack->GetY11()-yTrack;
      distTriggerTrack[2] = TMath::Tan(triggerTrack->GetThetay())-ySlopeTrack;

      if(isCovOK){
	TMatrixD paramDiff(3,1);
	for(Int_t iVar = 0; iVar < 3; iVar++)
	  paramDiff(iVar,0) = distTriggerTrack[iVar];
	
      	TMatrixD tmp(sumCov,TMatrixD::kMult,paramDiff);
	TMatrixD chi2M(paramDiff,TMatrixD::kTransposeMult,tmp);
	chi2 = chi2M(0,0);
      }
      else {
	chi2 = 0.;
	for (Int_t iVar = 0; iVar < 3; iVar++) chi2 += distTriggerTrack[iVar]*distTriggerTrack[iVar]/sigma2[iVar];
      }

      chi2 /= 3.; // Normalized Chi2: 3 degrees of freedom (X,Y,slopeY)
      if (chi2 < AliMUONReconstructor::GetRecoParam()->GetMaxNormChi2MatchTrigger()) 
      {
        Bool_t isDoubleTrack = (TMath::Abs(chi2 - minChi2MatchTrigger)<1.);
        if (chi2 < minChi2MatchTrigger && chi2 < AliMUONReconstructor::GetRecoParam()->GetMaxNormChi2MatchTrigger()) 
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

