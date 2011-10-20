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

#include "AliMUONReconstructor.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONObjectPair.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONConstants.h"
#include "AliMUONTrackParam.h"

#include "AliLog.h"

#include <TMath.h>

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrack) // Class implementation in ROOT context
/// \endcond


const Double_t AliMUONTrack::fgkMaxChi2 = 1.e10; ///< maximum chi2 above which the track can be considered as abnormal


//__________________________________________________________________________
AliMUONTrack::AliMUONTrack()
  : TObject(),
    fTrackParamAtCluster(0x0),
    fFitWithVertex(kFALSE),
    fVertexErrXY2(),
    fFitWithMCS(kFALSE),
    fClusterWeightsNonBending(0x0),
    fClusterWeightsBending(0x0),
    fGlobalChi2(-1.),
    fImproved(kFALSE),
    fMatchTrigger(-1),
    fChi2MatchTrigger(0.),
    fTrackID(-1),
    fTrackParamAtVertex(0x0),
    fHitsPatternInTrigCh(0),
    fLocalTrigger(0),
    fConnected(kFALSE)
{
  /// Default constructor
  fVertexErrXY2[0] = 0.;
  fVertexErrXY2[1] = 0.;
}

  //__________________________________________________________________________
AliMUONTrack::AliMUONTrack(AliMUONObjectPair *segment, Double_t bendingVertexDispersion)
  : TObject(),
    fTrackParamAtCluster(new TObjArray(20)),
    fFitWithVertex(kFALSE),
    fVertexErrXY2(),
    fFitWithMCS(kFALSE),
    fClusterWeightsNonBending(0x0),
    fClusterWeightsBending(0x0),
    fGlobalChi2(0.),
    fImproved(kFALSE),
    fMatchTrigger(-1),
    fChi2MatchTrigger(0.),
    fTrackID(-1),
    fTrackParamAtVertex(0x0),
    fHitsPatternInTrigCh(0),
    fLocalTrigger(0),
    fConnected(kFALSE)
{
  /// Constructor from two clusters
  
  fTrackParamAtCluster->SetOwner(kTRUE);
  
  fVertexErrXY2[0] = 0.;
  fVertexErrXY2[1] = 0.;
  
  // Pointers to clusters from the segment
  AliMUONVCluster* firstCluster = (AliMUONVCluster*) segment->First();
  AliMUONVCluster* lastCluster = (AliMUONVCluster*) segment->Second();
  
  // Compute track parameters
  Double_t z1 = firstCluster->GetZ();
  Double_t z2 = lastCluster->GetZ();
  Double_t dZ = z1 - z2;
  // Non bending plane
  Double_t nonBendingCoor1 = firstCluster->GetX();
  Double_t nonBendingCoor2 = lastCluster->GetX();
  Double_t nonBendingSlope = (nonBendingCoor1 - nonBendingCoor2) / dZ;
  // Bending plane
  Double_t bendingCoor1 = firstCluster->GetY();
  Double_t bendingCoor2 = lastCluster->GetY();
  Double_t bendingSlope = (bendingCoor1 - bendingCoor2) / dZ;
  // Inverse bending momentum
  Double_t bendingImpact = bendingCoor1 - z1 * bendingSlope;
  Double_t inverseBendingMomentum = 1. / AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(bendingImpact);
  
  // Set track parameters at first cluster
  AliMUONTrackParam trackParamAtFirstCluster;
  trackParamAtFirstCluster.SetZ(z1);
  trackParamAtFirstCluster.SetNonBendingCoor(nonBendingCoor1);
  trackParamAtFirstCluster.SetNonBendingSlope(nonBendingSlope);
  trackParamAtFirstCluster.SetBendingCoor(bendingCoor1);
  trackParamAtFirstCluster.SetBendingSlope(bendingSlope);
  trackParamAtFirstCluster.SetInverseBendingMomentum(inverseBendingMomentum);
  
  // Set track parameters at last cluster
  AliMUONTrackParam trackParamAtLastCluster;
  trackParamAtLastCluster.SetZ(z2);
  trackParamAtLastCluster.SetNonBendingCoor(nonBendingCoor2);
  trackParamAtLastCluster.SetNonBendingSlope(nonBendingSlope);
  trackParamAtLastCluster.SetBendingCoor(bendingCoor2);
  trackParamAtLastCluster.SetBendingSlope(bendingSlope);
  trackParamAtLastCluster.SetInverseBendingMomentum(inverseBendingMomentum);
  
  // Compute and set track parameters covariances at first cluster
  TMatrixD paramCov(5,5);
  paramCov.Zero();
  // Non bending plane
  paramCov(0,0) = firstCluster->GetErrX2();
  paramCov(0,1) = firstCluster->GetErrX2() / dZ;
  paramCov(1,0) = paramCov(0,1);
  paramCov(1,1) = ( firstCluster->GetErrX2() + lastCluster->GetErrX2() ) / dZ / dZ;
  // Bending plane
  paramCov(2,2) = firstCluster->GetErrY2();
  paramCov(2,3) = firstCluster->GetErrY2() / dZ;
  paramCov(3,2) = paramCov(2,3);
  paramCov(3,3) = ( firstCluster->GetErrY2() + lastCluster->GetErrY2() ) / dZ / dZ;
  // Inverse bending momentum (vertex resolution + bending slope resolution + 10% error on dipole parameters+field)
  if (AliMUONTrackExtrap::IsFieldON()) {
    paramCov(4,4) = ( ( bendingVertexDispersion*bendingVertexDispersion +
		       (z1 * z1 * lastCluster->GetErrY2() + z2 * z2 * firstCluster->GetErrY2()) / dZ / dZ) /
		     bendingImpact / bendingImpact + 0.1 * 0.1) * inverseBendingMomentum * inverseBendingMomentum ;
    paramCov(2,4) = - z2 * firstCluster->GetErrY2() * inverseBendingMomentum / bendingImpact / dZ;
    paramCov(4,2) = paramCov(2,4);
    paramCov(3,4) = - (z1 * lastCluster->GetErrY2() + z2 * firstCluster->GetErrY2()) * inverseBendingMomentum / bendingImpact / dZ / dZ;
    paramCov(4,3) = paramCov(3,4);
  } else paramCov(4,4) = inverseBendingMomentum*inverseBendingMomentum;
  trackParamAtFirstCluster.SetCovariances(paramCov);
  
  // Compute and set track parameters covariances at last cluster
  // Non bending plane
  paramCov(0,0) = lastCluster->GetErrX2();
  paramCov(0,1) = - lastCluster->GetErrX2() / dZ;
  paramCov(1,0) = paramCov(0,1);
  // Bending plane
  paramCov(2,2) = lastCluster->GetErrY2();
  paramCov(2,3) = - lastCluster->GetErrY2() / dZ;
  paramCov(3,2) = paramCov(2,3);
  // Inverse bending momentum (vertex resolution + bending slope resolution + 10% error on dipole parameters+field)
  if (AliMUONTrackExtrap::IsFieldON()) {
    paramCov(2,4) = z1 * lastCluster->GetErrY2() * inverseBendingMomentum / bendingImpact / dZ;
    paramCov(4,2) = paramCov(2,4);
  }
  trackParamAtLastCluster.SetCovariances(paramCov);
  
  // Add track parameters at clusters
  AddTrackParamAtCluster(trackParamAtFirstCluster,*firstCluster);
  AddTrackParamAtCluster(trackParamAtLastCluster,*lastCluster);
  
}

//__________________________________________________________________________
AliMUONTrack::AliMUONTrack(const AliMUONTrack& track)
  : TObject(track),
    fTrackParamAtCluster(0x0),
    fFitWithVertex(track.fFitWithVertex),
    fVertexErrXY2(),
    fFitWithMCS(track.fFitWithMCS),
    fClusterWeightsNonBending(0x0),
    fClusterWeightsBending(0x0),
    fGlobalChi2(track.fGlobalChi2),
    fImproved(track.fImproved),
    fMatchTrigger(track.fMatchTrigger),
    fChi2MatchTrigger(track.fChi2MatchTrigger),
    fTrackID(track.fTrackID),
    fTrackParamAtVertex(0x0),
    fHitsPatternInTrigCh(track.fHitsPatternInTrigCh),
    fLocalTrigger(track.fLocalTrigger),
    fConnected(track.fConnected)
{
  ///copy constructor
  
  // necessary to make a copy of the objects and not only the pointers in TObjArray.
  if (track.fTrackParamAtCluster) {
    fTrackParamAtCluster = new TObjArray(track.fTrackParamAtCluster->GetSize());
    fTrackParamAtCluster->SetOwner(kTRUE);
    for (Int_t i = 0; i < track.GetNClusters(); i++)
      fTrackParamAtCluster->AddLast(new AliMUONTrackParam(*static_cast<AliMUONTrackParam*>(track.fTrackParamAtCluster->UncheckedAt(i))));
  }
  
  // copy vertex resolution square used during the tracking procedure
  fVertexErrXY2[0] = track.fVertexErrXY2[0];
  fVertexErrXY2[1] = track.fVertexErrXY2[1];
  
  // copy cluster weights matrices if any
  if (track.fClusterWeightsNonBending) fClusterWeightsNonBending = new TMatrixD(*(track.fClusterWeightsNonBending));
  if (track.fClusterWeightsBending) fClusterWeightsBending = new TMatrixD(*(track.fClusterWeightsBending));
  
  // copy track parameters at vertex if any
  if (track.fTrackParamAtVertex) fTrackParamAtVertex = new AliMUONTrackParam(*(track.fTrackParamAtVertex));
  
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
  
  // clear memory
  Clear();
  
  // necessary to make a copy of the objects and not only the pointers in TObjArray
  if (track.fTrackParamAtCluster) {
    fTrackParamAtCluster = new TObjArray(track.fTrackParamAtCluster->GetSize());
    fTrackParamAtCluster->SetOwner(kTRUE);
    for (Int_t i = 0; i < track.GetNClusters(); i++)
      fTrackParamAtCluster->AddLast(new AliMUONTrackParam(*static_cast<AliMUONTrackParam*>(track.fTrackParamAtCluster->UncheckedAt(i))));
  }
  
  // copy cluster weights matrix if any
  if (track.fClusterWeightsNonBending) {
    if (fClusterWeightsNonBending) {
      fClusterWeightsNonBending->ResizeTo(*(track.fClusterWeightsNonBending));
      *fClusterWeightsNonBending = *(track.fClusterWeightsNonBending);
    } else fClusterWeightsNonBending = new TMatrixD(*(track.fClusterWeightsNonBending));
  }
  
  // copy cluster weights matrix if any
  if (track.fClusterWeightsBending) {
    if (fClusterWeightsBending) {
      fClusterWeightsBending->ResizeTo(*(track.fClusterWeightsBending));
      *fClusterWeightsBending = *(track.fClusterWeightsBending);
    } else fClusterWeightsBending = new TMatrixD(*(track.fClusterWeightsBending));
  }
  
  // copy track parameters at vertex if any
  if (track.fTrackParamAtVertex) {
    if (fTrackParamAtVertex) *fTrackParamAtVertex = *(track.fTrackParamAtVertex);
    else fTrackParamAtVertex = new AliMUONTrackParam(*(track.fTrackParamAtVertex));
  }
  
  fFitWithVertex      =  track.fFitWithVertex;
  fVertexErrXY2[0]    =  track.fVertexErrXY2[0];
  fVertexErrXY2[1]    =  track.fVertexErrXY2[1];
  fFitWithMCS         =  track.fFitWithMCS;
  fGlobalChi2         =  track.fGlobalChi2;
  fImproved           =  track.fImproved;
  fMatchTrigger       =  track.fMatchTrigger;
  fChi2MatchTrigger   =  track.fChi2MatchTrigger;
  fTrackID            =  track.fTrackID; 
  fHitsPatternInTrigCh = track.fHitsPatternInTrigCh;
  fLocalTrigger        = track.fLocalTrigger;
  fConnected          =  track.fConnected;

  return *this;
}

  //__________________________________________________________________________
AliMUONTrack::~AliMUONTrack()
{
  /// Destructor
  delete fTrackParamAtCluster;
  delete fClusterWeightsNonBending;
  delete fClusterWeightsBending;
  delete fTrackParamAtVertex;
}

  //__________________________________________________________________________
void AliMUONTrack::Clear(Option_t* /*opt*/)
{
  /// Clear arrays
  delete fTrackParamAtCluster; fTrackParamAtCluster = 0x0;
  delete fClusterWeightsNonBending; fClusterWeightsNonBending = 0x0;
  delete fClusterWeightsBending; fClusterWeightsBending = 0x0;
  delete fTrackParamAtVertex; fTrackParamAtVertex = 0x0;
}

  //__________________________________________________________________________
void AliMUONTrack::Reset()
{
  /// Reset to default values
  SetUniqueID(0);
  fFitWithVertex = kFALSE;
  fVertexErrXY2[0] = 0.;
  fVertexErrXY2[1] = 0.;
  fFitWithMCS = kFALSE;
  fGlobalChi2 = -1.;
  fImproved = kFALSE;
  fMatchTrigger = -1;
  fChi2MatchTrigger = 0.;
  fTrackID = -1;
  fHitsPatternInTrigCh = 0;
  fLocalTrigger = 0;
  fConnected = kFALSE;
  delete fTrackParamAtCluster; fTrackParamAtCluster = 0x0;
  delete fClusterWeightsNonBending; fClusterWeightsNonBending = 0x0;
  delete fClusterWeightsBending; fClusterWeightsBending = 0x0;
  delete fTrackParamAtVertex; fTrackParamAtVertex = 0x0;
}

  //__________________________________________________________________________
TObjArray* AliMUONTrack::GetTrackParamAtCluster() const
{
  /// return array of track parameters at cluster (create it if needed)
  if (!fTrackParamAtCluster) {
    fTrackParamAtCluster = new TObjArray(20);
    fTrackParamAtCluster->SetOwner(kTRUE);
  }
  return fTrackParamAtCluster;
}

  //__________________________________________________________________________
void AliMUONTrack::AddTrackParamAtCluster(const AliMUONTrackParam &trackParam, AliMUONVCluster &cluster, Bool_t copy)
{
  /// Copy given track parameters into a new TrackParamAtCluster
  /// Link parameters with the associated cluster
  /// If copy=kTRUE: the cluster is copied then passed the trackParam which become its owner 
  ///     otherwise: make sure to do not delete the cluster until it is used by the track
  
  // check chamber ID of the associated cluster
  if (cluster.GetChamberId() < 0 || cluster.GetChamberId() > AliMUONConstants::NTrackingCh()) {
    AliError(Form("Chamber ID of the associated cluster is not valid (ChamberId=%d)",cluster.GetChamberId()));
    return;
  }
  
  // check whether track parameters are given at the correct cluster z position
  if (TMath::Abs(cluster.GetZ() - trackParam.GetZ())>1.e-5) {   // AU
    AliError("track parameters are given at a different z position than the one of the associated cluster");
    return;
  }
  
  // add parameters to the array of track parameters
  if (!fTrackParamAtCluster) {
    fTrackParamAtCluster = new TObjArray(20);
    fTrackParamAtCluster->SetOwner(kTRUE);
  }
  AliMUONTrackParam* trackParamAtCluster = new AliMUONTrackParam(trackParam);
  fTrackParamAtCluster->AddLast(trackParamAtCluster);
  
  // link parameters with the associated cluster or its copy
  if (copy) {
    AliMUONVCluster *clusterCopy = static_cast<AliMUONVCluster*>(cluster.Clone());
    trackParamAtCluster->SetClusterPtr(clusterCopy, kTRUE);
  } else trackParamAtCluster->SetClusterPtr(&cluster);
  
  // sort the array of track parameters
  fTrackParamAtCluster->Sort();
}

  //__________________________________________________________________________
void AliMUONTrack::RemoveTrackParamAtCluster(AliMUONTrackParam *trackParam)
{
  /// Remove trackParam from the array of TrackParamAtCluster and delete it since the array is owner
  
  if (fTrackParamAtCluster) {
    
    AliMUONTrackParam* trackParamAtCluster = static_cast<AliMUONTrackParam*>(fTrackParamAtCluster->Remove(trackParam));
    
    if (trackParamAtCluster) {
      
      // clean memory
      delete trackParamAtCluster;
      
      // remove hole
      fTrackParamAtCluster->Compress();
      
    } else AliWarning("object to remove does not exist in array fTrackParamAtCluster");
    
  } else AliWarning("array fTrackParamAtCluster does not exist");
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrack::UpdateTrackParamAtCluster()
{
  /// Update track parameters at each attached cluster
  /// Return kFALSE in case of failure (i.e. extrapolation problem)
  
  Int_t nClusters = GetNClusters();
  if (nClusters == 0) {
    AliWarning("no cluster attached to the track");
    return kFALSE;
  }
  
  Bool_t extrapStatus = kTRUE;
  AliMUONTrackParam* startingTrackParam = static_cast<AliMUONTrackParam*>(fTrackParamAtCluster->UncheckedAt(0));
  
  for (Int_t i = 1; i < nClusters; i++) {
    AliMUONTrackParam* trackParamAtCluster = static_cast<AliMUONTrackParam*>(fTrackParamAtCluster->UncheckedAt(i));
    
    // reset track parameters and their covariances
    trackParamAtCluster->SetParameters(startingTrackParam->GetParameters());
    trackParamAtCluster->SetZ(startingTrackParam->GetZ());
    
    // extrapolation to the given z
    if (!AliMUONTrackExtrap::ExtrapToZ(trackParamAtCluster, trackParamAtCluster->GetClusterPtr()->GetZ())) extrapStatus = kFALSE;
    
    // prepare next step
    startingTrackParam = trackParamAtCluster;
  }

  // set global chi2 to max value in case of problem during track extrapolation
  if (!extrapStatus) SetGlobalChi2(2.*MaxChi2());
  return extrapStatus;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrack::UpdateCovTrackParamAtCluster()
{
  /// Update track parameters and their covariances at each attached cluster
  /// Include effects of multiple scattering in chambers
  /// Return kFALSE in case of failure (i.e. extrapolation problem)
  
  Int_t nClusters = GetNClusters();
  if (nClusters == 0) {
    AliWarning("no cluster attached to the track");
    return kFALSE;
  }
  
  Bool_t extrapStatus = kTRUE;
  AliMUONTrackParam* startingTrackParam = static_cast<AliMUONTrackParam*>(fTrackParamAtCluster->UncheckedAt(0));
  Int_t expectedChamber = startingTrackParam->GetClusterPtr()->GetChamberId() + 1;
  Int_t currentChamber;
  
  for (Int_t i = 1; i < nClusters; i++) {
    AliMUONTrackParam* trackParamAtCluster = static_cast<AliMUONTrackParam*>(fTrackParamAtCluster->UncheckedAt(i));
    
    // reset track parameters and their covariances
    trackParamAtCluster->SetParameters(startingTrackParam->GetParameters());
    trackParamAtCluster->SetZ(startingTrackParam->GetZ());
    trackParamAtCluster->SetCovariances(startingTrackParam->GetCovariances());
    
    // add MCS effect
    AliMUONTrackExtrap::AddMCSEffect(trackParamAtCluster,AliMUONConstants::ChamberThicknessInX0(expectedChamber-1),-1.);
    
    // add MCS in missing chambers if any
    currentChamber = trackParamAtCluster->GetClusterPtr()->GetChamberId();
    while (currentChamber > expectedChamber) {
      // extrapolation to the missing chamber
      if (!AliMUONTrackExtrap::ExtrapToZCov(trackParamAtCluster, AliMUONConstants::DefaultChamberZ(expectedChamber))) extrapStatus = kFALSE;
      // add MCS effect
      AliMUONTrackExtrap::AddMCSEffect(trackParamAtCluster,AliMUONConstants::ChamberThicknessInX0(expectedChamber),-1.);
      expectedChamber++;
    }
    
    // extrapolation to the z of the current cluster
    if (!AliMUONTrackExtrap::ExtrapToZCov(trackParamAtCluster, trackParamAtCluster->GetClusterPtr()->GetZ())) extrapStatus = kFALSE;
    
    // prepare next step
    expectedChamber = currentChamber + 1;
    startingTrackParam = trackParamAtCluster;
  }
  
  // set global chi2 to max value in case of problem during track extrapolation
  if (!extrapStatus) SetGlobalChi2(2.*MaxChi2());
  return extrapStatus;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrack::IsValid(UInt_t requestedStationMask, Bool_t request2ChInSameSt45)
{
  /// check the validity of the current track:
  /// at least one cluster per requested station
  /// and at least 2 chambers in stations 4 & 5 that contain cluster(s)
  /// + if request2ChInSameSt45 = kTRUE: 2 chambers hit in the same station (4 or 5)
  
  Int_t nClusters = GetNClusters();
  AliMUONTrackParam *trackParam;
  Int_t currentCh, currentSt, previousCh = -1, nChHitInSt4 = 0, nChHitInSt5 = 0;
  UInt_t presentStationMask(0);
  
  // first loop over clusters
  for (Int_t i = 0; i < nClusters; i++) {
    trackParam = (AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(i);
    
    currentCh = trackParam->GetClusterPtr()->GetChamberId();
    currentSt = currentCh/2;
    
    // build present station mask
    presentStationMask |= ( 1 << currentSt );
    
    // count the number of chambers hit in station 4 that contain cluster(s)
    if (currentSt == 3 && currentCh != previousCh) {
      nChHitInSt4++;
      previousCh = currentCh;
    }
    
    // count the number of chambers hit in station 5 that contain cluster(s)
    if (currentSt == 4 && currentCh != previousCh) {
      nChHitInSt5++;
      previousCh = currentCh;
    }
    
  }
  
  // at least one cluster per requested station
  if ((requestedStationMask & presentStationMask) != requestedStationMask) return kFALSE;
  
  // 2 chambers hit in the same station (4 or 5)
  if (request2ChInSameSt45) return (nChHitInSt4 == 2 || nChHitInSt5 == 2);
  // or 2 chambers hit in station 4 & 5 together
  else return (nChHitInSt4+nChHitInSt5 >= 2);
  
}

  //__________________________________________________________________________
void AliMUONTrack::TagRemovableClusters(UInt_t requestedStationMask) {
  /// Identify clusters that can be removed from the track,
  /// with the only requirements to have at least 1 cluster per requested station
  /// and at least 2 chambers over 4 in stations 4 & 5 that contain cluster(s)
  
  Int_t nClusters = GetNClusters();
  AliMUONTrackParam *trackParam, *nextTrackParam;
  Int_t currentCh, nextCh, currentSt, nextSt, previousCh = -1, nChHitInSt45 = 0;
  
  // first loop over clusters
  for (Int_t i = 0; i < nClusters; i++) {
    trackParam = (AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(i);
    
    currentCh = trackParam->GetClusterPtr()->GetChamberId();
    currentSt = currentCh/2;
    
    // reset flags to kFALSE for all clusters in required station
    if ((1 << currentSt) & requestedStationMask) trackParam->SetRemovable(kFALSE);
    else trackParam->SetRemovable(kTRUE);
    
    // count the number of chambers in station 4 & 5 that contain cluster(s)
    if (currentCh > 5 && currentCh != previousCh) {
      nChHitInSt45++;
      previousCh = currentCh;
    }
    
  }
  
  // second loop over clusters
  for (Int_t i = 0; i < nClusters; i++) {
    trackParam = (AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(i);
    
    currentCh = trackParam->GetClusterPtr()->GetChamberId();
    currentSt = currentCh/2;
    
    // make sure they are more than 2 clusters in 2 different chambers of stations 4 & 5
    // but 2 clusters in he same chamber will still be flagged as removable
    if (nChHitInSt45 < 3 && currentSt > 2) {
      
      if (i == nClusters-1) {
	
	trackParam->SetRemovable(kFALSE);
      
      } else {
	
	nextTrackParam = (AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(i+1);
	nextCh = nextTrackParam->GetClusterPtr()->GetChamberId();
	
	// set clusters in the same chamber as being removable
	if (nextCh == currentCh) {
	  trackParam->SetRemovable(kTRUE);
	  nextTrackParam->SetRemovable(kTRUE);
	  i++; // skip cluster already checked
	} else {
	  trackParam->SetRemovable(kFALSE);
	}
	
      }
      
    } else {
      
      // skip clusters already flag as removable
      if (trackParam->IsRemovable()) continue;
      
      // loop over next track parameters
      for (Int_t j = i+1; j < nClusters; j++) {
	nextTrackParam = (AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(j);
	
	nextCh = nextTrackParam->GetClusterPtr()->GetChamberId();
	nextSt = nextCh/2;
	
	// set clusters in the same station as being removable
	if (nextSt == currentSt) {
	  trackParam->SetRemovable(kTRUE);
	  nextTrackParam->SetRemovable(kTRUE);
	  i++; // skip cluster already checked
	}
	
      }
      
    }
      
  }
    
}

  //__________________________________________________________________________
Bool_t AliMUONTrack::ComputeLocalChi2(Bool_t accountForMCS)
{
  /// Compute each cluster contribution to the chi2 of the track
  /// accounting for multiple scattering or not according to the flag
  /// - Also recompute the weight matrices of the attached clusters if accountForMCS=kTRUE
  /// - Assume that track parameters at each cluster are corrects
  /// - Return kFALSE if computation failed
  AliDebug(1,"Enter ComputeLocalChi2");
  
  if (!fTrackParamAtCluster) {
    AliWarning("no cluster attached to this track");
    return kFALSE;
  }
  
  if (accountForMCS) { // Compute local chi2 taking into account multiple scattering effects
      
    // Compute MCS covariance matrix only once
    Int_t nClusters = GetNClusters();
    TMatrixD mcsCovariances(nClusters,nClusters);
    ComputeMCSCovariances(mcsCovariances);
    
    // Make sure cluster weights are consistent with following calculations
    if (!ComputeClusterWeights(&mcsCovariances)) {
      AliWarning("cannot take into account the multiple scattering effects");
      return ComputeLocalChi2(kFALSE);
    }
    
    // Compute chi2 of the track
    Double_t globalChi2 = ComputeGlobalChi2(kTRUE);
    if (globalChi2 < 0.) return kFALSE;
    
    // Loop over removable clusters and compute their local chi2
    AliMUONTrackParam* trackParamAtCluster;
    AliMUONTrackParam* trackParamAtCluster1;
    AliMUONVCluster *cluster, *discardedCluster;
    Int_t iCluster1, iCluster2, iCurrentCluster1, iCurrentCluster2;
    TMatrixD clusterWeightsNB(nClusters-1,nClusters-1);
    TMatrixD clusterWeightsB(nClusters-1,nClusters-1);
    Double_t *dX = new Double_t[nClusters-1];
    Double_t *dY = new Double_t[nClusters-1];
    Double_t globalChi2b;
    for (Int_t iCluster = 0; iCluster < nClusters ; iCluster++) { 
      trackParamAtCluster = static_cast<AliMUONTrackParam*>(fTrackParamAtCluster->UncheckedAt(iCluster));
      
      discardedCluster = trackParamAtCluster->GetClusterPtr();
      
      // Recompute cluster weights without the current cluster
      if (!ComputeClusterWeights(clusterWeightsNB, clusterWeightsB, &mcsCovariances, discardedCluster)) {
  	AliWarning("cannot take into account the multiple scattering effects");
	delete [] dX;
	delete [] dY;
	return ComputeLocalChi2(kFALSE);
      }
      
      // Compute track chi2 without the current cluster
      globalChi2b = 0.;
      iCurrentCluster1 = 0;
      for (iCluster1 = 0; iCluster1 < nClusters ; iCluster1++) { 
    	trackParamAtCluster1 = (AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(iCluster1);
    	cluster = trackParamAtCluster1->GetClusterPtr();
        
        if (cluster == discardedCluster) continue;
        
        // Compute and save residuals
    	dX[iCurrentCluster1] = cluster->GetX() - trackParamAtCluster1->GetNonBendingCoor();
    	dY[iCurrentCluster1] = cluster->GetY() - trackParamAtCluster1->GetBendingCoor();
        
        iCurrentCluster2 = 0;
    	for (iCluster2 = 0; iCluster2 < iCluster1; iCluster2++) {
    	  cluster = ((AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(iCluster2))->GetClusterPtr();
          
          if (cluster == discardedCluster) continue;
          
          // Add contribution from covariances
          globalChi2b += (clusterWeightsNB(iCurrentCluster1, iCurrentCluster2) +
        		  clusterWeightsNB(iCurrentCluster2, iCurrentCluster1)) * dX[iCurrentCluster1] * dX[iCurrentCluster2] +
        		 (clusterWeightsB(iCurrentCluster1, iCurrentCluster2) +
        		  clusterWeightsB(iCurrentCluster2, iCurrentCluster1)) * dY[iCurrentCluster1] * dY[iCurrentCluster2];
          
          iCurrentCluster2++;
    	}
        
        // Add contribution from variances
    	globalChi2b += clusterWeightsNB(iCurrentCluster1, iCurrentCluster1) * dX[iCurrentCluster1] * dX[iCurrentCluster1] +
        	       clusterWeightsB(iCurrentCluster1, iCurrentCluster1) * dY[iCurrentCluster1] * dY[iCurrentCluster1];
    	
        iCurrentCluster1++;
      }

      // Set local chi2
      trackParamAtCluster->SetLocalChi2(globalChi2 - globalChi2b);
    }
    
    delete [] dX;
    delete [] dY;
    
  } else { // without multiple scattering effects
    
    Int_t nClusters = GetNClusters();
    AliMUONTrackParam* trackParamAtCluster;
    AliMUONVCluster *discardedCluster;
    Double_t dX, dY;
    for (Int_t iCluster = 0; iCluster < nClusters ; iCluster++) { 
      trackParamAtCluster = static_cast<AliMUONTrackParam*>(fTrackParamAtCluster->UncheckedAt(iCluster));
      
      discardedCluster = trackParamAtCluster->GetClusterPtr();
      
      // Compute residuals
      dX = discardedCluster->GetX() - trackParamAtCluster->GetNonBendingCoor();
      dY = discardedCluster->GetY() - trackParamAtCluster->GetBendingCoor();
      
      // Set local chi2
      trackParamAtCluster->SetLocalChi2(dX * dX / discardedCluster->GetErrX2() + dY * dY / discardedCluster->GetErrY2());
    }
  
  }
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Double_t AliMUONTrack::ComputeGlobalChi2(Bool_t accountForMCS)
{
  /// Compute the chi2 of the track accounting for multiple scattering or not according to the flag
  /// - Assume that track parameters at each cluster are corrects
  /// - Assume the cluster weights matrices are corrects
  /// - Return a value of chi2 higher than the maximum allowed if computation failed
  AliDebug(1,"Enter ComputeGlobalChi2");
  
  if (!fTrackParamAtCluster) {
    AliWarning("no cluster attached to this track");
    return 2.*MaxChi2();
  }
  
  Double_t chi2 = 0.;
  
  if (accountForMCS) {
    
    // Check the weight matrices. If weight matrices are not available compute chi2 without MCS
    if (!fClusterWeightsNonBending || !fClusterWeightsBending) {
      AliWarning("cluster weights including multiple scattering effects are not available\n\t\t --> compute chi2 WITHOUT multiple scattering");
      return ComputeGlobalChi2(kFALSE);
    }
    Int_t nClusters = GetNClusters();
    if (fClusterWeightsNonBending->GetNrows() != nClusters || fClusterWeightsBending->GetNcols() != nClusters) {
      AliWarning("cluster weights including multiple scattering effects are not available\n\t\t --> compute chi2 WITHOUT multiple scattering");
      return ComputeGlobalChi2(kFALSE);
    }
    
    // Compute chi2
    AliMUONVCluster *cluster;
    Double_t *dX = new Double_t[nClusters];
    Double_t *dY = new Double_t[nClusters];
    AliMUONTrackParam* trackParamAtCluster;
    for (Int_t iCluster1 = 0; iCluster1 < nClusters; iCluster1++) { 
      trackParamAtCluster = (AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(iCluster1);
      cluster = trackParamAtCluster->GetClusterPtr();
      dX[iCluster1] = cluster->GetX() - trackParamAtCluster->GetNonBendingCoor();
      dY[iCluster1] = cluster->GetY() - trackParamAtCluster->GetBendingCoor();
      for (Int_t iCluster2 = 0; iCluster2 < iCluster1; iCluster2++) {
        chi2 += ((*fClusterWeightsNonBending)(iCluster1, iCluster2) + (*fClusterWeightsNonBending)(iCluster2, iCluster1)) * dX[iCluster1] * dX[iCluster2] +
		((*fClusterWeightsBending)(iCluster1, iCluster2) + (*fClusterWeightsBending)(iCluster2, iCluster1)) * dY[iCluster1] * dY[iCluster2];
      }
      chi2 += ((*fClusterWeightsNonBending)(iCluster1, iCluster1) * dX[iCluster1] * dX[iCluster1]) +
	      ((*fClusterWeightsBending)(iCluster1, iCluster1) * dY[iCluster1] * dY[iCluster1]);
    }
    delete [] dX;
    delete [] dY;
    
  } else {
    
    AliMUONVCluster *cluster;
    Double_t dX, dY;
    AliMUONTrackParam* trackParamAtCluster;
    Int_t nClusters = GetNClusters();
    for (Int_t iCluster = 0; iCluster < nClusters ; iCluster++) { 
      trackParamAtCluster = (AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(iCluster);
      cluster = trackParamAtCluster->GetClusterPtr();
      dX = cluster->GetX() - trackParamAtCluster->GetNonBendingCoor();
      dY = cluster->GetY() - trackParamAtCluster->GetBendingCoor();
      chi2 += dX * dX / cluster->GetErrX2() + dY * dY / cluster->GetErrY2();
    }
    
  }
  
  return chi2;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrack::ComputeClusterWeights(TMatrixD* mcsCovariances)
{
  /// Compute the weight matrices of the attached clusters, in non bending and bending direction,
  /// accounting for multiple scattering correlations and cluster resolution
  /// - Use the provided MCS covariance matrix if any (otherwise build it temporarily)
  /// - Assume that track parameters at each cluster are corrects
  /// - Return kFALSE if computation failed
  AliDebug(1,"Enter ComputeClusterWeights1");
  
  if (!fTrackParamAtCluster) {
    AliWarning("no cluster attached to this track");
    return kFALSE;
  }
  
  // Alocate memory
  Int_t nClusters = GetNClusters();
  if (!fClusterWeightsNonBending) fClusterWeightsNonBending = new TMatrixD(nClusters,nClusters);
  if (!fClusterWeightsBending) fClusterWeightsBending = new TMatrixD(nClusters,nClusters);
  
  // Compute weights matrices
  if (!ComputeClusterWeights(*fClusterWeightsNonBending, *fClusterWeightsBending, mcsCovariances)) return kFALSE;
  
  return kTRUE;
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrack::ComputeClusterWeights(TMatrixD& clusterWeightsNB, TMatrixD& clusterWeightsB,
					   TMatrixD* mcsCovariances, const AliMUONVCluster* discardedCluster) const
{
  /// Compute the weight matrices, in non bending and bending direction,
  /// of the other attached clusters assuming the discarded one does not exist
  /// accounting for multiple scattering correlations and cluster resolution
  /// - Use the provided MCS covariance matrix if any (otherwise build it temporarily)
  /// - Return kFALSE if computation failed
  AliDebug(1,"Enter ComputeClusterWeights2");
  
  // Check MCS covariance matrix and recompute it if need
  Int_t nClusters = GetNClusters();
  Bool_t deleteMCSCov = kFALSE;
  if (!mcsCovariances) {
    mcsCovariances = new TMatrixD(nClusters,nClusters);
    deleteMCSCov = kTRUE;
    ComputeMCSCovariances(*mcsCovariances);
  }
  
  // Resize the weights matrices; alocate memory
  if (discardedCluster) {
    clusterWeightsNB.ResizeTo(nClusters-1,nClusters-1);
    clusterWeightsB.ResizeTo(nClusters-1,nClusters-1);
  } else {
    clusterWeightsNB.ResizeTo(nClusters,nClusters);
    clusterWeightsB.ResizeTo(nClusters,nClusters);
  }
  
  // Define variables
  AliMUONVCluster *cluster1, *cluster2;
  Int_t iCurrentCluster1, iCurrentCluster2;
  
  // Compute the covariance matrices
  iCurrentCluster1 = 0;
  for (Int_t iCluster1 = 0; iCluster1 < nClusters; iCluster1++) { 
    cluster1 = ((AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(iCluster1))->GetClusterPtr();
    
    if (cluster1 == discardedCluster) continue;
    
    // Loop over next clusters
    iCurrentCluster2 = iCurrentCluster1;
    for (Int_t iCluster2 = iCluster1; iCluster2 < nClusters; iCluster2++) {
      cluster2 = ((AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(iCluster2))->GetClusterPtr();
      
      if (cluster2 == discardedCluster) continue;
      
      // Fill with MCS covariances
      clusterWeightsNB(iCurrentCluster1, iCurrentCluster2) = (*mcsCovariances)(iCluster1,iCluster2);
      
      // Equal contribution from multiple scattering in non bending and bending directions
      clusterWeightsB(iCurrentCluster1, iCurrentCluster2) = clusterWeightsNB(iCurrentCluster1, iCurrentCluster2);
      
      // Add contribution from cluster resolution to diagonal element and symmetrize the matrix
      if (iCurrentCluster1 == iCurrentCluster2) {
	
	// In non bending plane
        clusterWeightsNB(iCurrentCluster1, iCurrentCluster1) += cluster1->GetErrX2();
	// In bending plane
	clusterWeightsB(iCurrentCluster1, iCurrentCluster1) += cluster1->GetErrY2();
	
      } else {
	
	// In non bending plane
	clusterWeightsNB(iCurrentCluster2, iCurrentCluster1) = clusterWeightsNB(iCurrentCluster1, iCurrentCluster2);
	// In bending plane
	clusterWeightsB(iCurrentCluster2, iCurrentCluster1) = clusterWeightsB(iCurrentCluster1, iCurrentCluster2);
	
      }
      
      iCurrentCluster2++;
    }
    
    iCurrentCluster1++;
  }
    
  // Inversion of covariance matrices to get the weights
  if (clusterWeightsNB.Determinant() != 0 && clusterWeightsB.Determinant() != 0) {
    clusterWeightsNB.Invert();
    clusterWeightsB.Invert();
  } else {
    AliWarning(" Determinant = 0");
    clusterWeightsNB.ResizeTo(0,0);
    clusterWeightsB.ResizeTo(0,0);
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
  /// (assume that track parameters at each cluster are corrects)
  AliDebug(1,"Enter ComputeMCSCovariances");
  
  // Reset the size of the covariance matrix if needed
  Int_t nClusters = GetNClusters();
  if (mcsCovariances.GetNrows() != nClusters) mcsCovariances.ResizeTo(nClusters,nClusters);
  
  // Define variables
  Int_t nChambers = AliMUONConstants::NTrackingCh();
  AliMUONTrackParam* trackParamAtCluster;
  AliMUONTrackParam extrapTrackParam;
  Int_t currentChamber = 0, expectedChamber = 0, size = 0;
  Double_t *mcsAngle2 = new Double_t[2*nChambers];
  Double_t *zMCS = new Double_t[2*nChambers];
  Int_t *indices = new Int_t[2*nClusters];
  
  // Compute multiple scattering dispersion angle at each chamber
  // and save the z position where it is calculated
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
    trackParamAtCluster = (AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(iCluster);
    
    // look for missing chambers if any
    currentChamber = trackParamAtCluster->GetClusterPtr()->GetChamberId();
    while (currentChamber > expectedChamber) {
      
      // Save the z position where MCS dispersion is calculated
      zMCS[size] = AliMUONConstants::DefaultChamberZ(expectedChamber);
      
      // Do not take into account MCS in chambers prior the first cluster
      if (iCluster > 0) {
        
        // Get track parameters at missing chamber z
        extrapTrackParam = *trackParamAtCluster;
        AliMUONTrackExtrap::ExtrapToZ(&extrapTrackParam, zMCS[size]);
        
        // Save multiple scattering dispersion angle in missing chamber
        mcsAngle2[size] = AliMUONTrackExtrap::GetMCSAngle2(extrapTrackParam,AliMUONConstants::ChamberThicknessInX0(expectedChamber),1.);
        
      } else mcsAngle2[size] = 0.;
      
      expectedChamber++;
      size++;
    }
    
    // Save z position where MCS dispersion is calculated
    zMCS[size] = trackParamAtCluster->GetZ();
    
    // Save multiple scattering dispersion angle in current chamber
    mcsAngle2[size] = AliMUONTrackExtrap::GetMCSAngle2(*trackParamAtCluster,AliMUONConstants::ChamberThicknessInX0(currentChamber),1.);
    
    // Save indice in zMCS array corresponding to the current cluster
    indices[iCluster] = size;
    
    expectedChamber = currentChamber + 1;
    size++;
  }
  
  // complete array of z if last cluster is on the last but one chamber
  if (currentChamber != nChambers-1) zMCS[size++] = AliMUONConstants::DefaultChamberZ(nChambers-1);
  
  // Compute the covariance matrix
  for (Int_t iCluster1 = 0; iCluster1 < nClusters; iCluster1++) { 
    
    for (Int_t iCluster2 = iCluster1; iCluster2 < nClusters; iCluster2++) {
      
      // Initialization to 0 (diagonal plus upper triangular part)
      mcsCovariances(iCluster1,iCluster2) = 0.;
      
      // Compute contribution from multiple scattering in upstream chambers
      for (Int_t k = 0; k < indices[iCluster1]; k++) { 	
	mcsCovariances(iCluster1,iCluster2) += (zMCS[indices[iCluster1]] - zMCS[k]) * (zMCS[indices[iCluster2]] - zMCS[k]) * mcsAngle2[k];
      }
      
      // Symetrize the matrix
      mcsCovariances(iCluster2,iCluster1) = mcsCovariances(iCluster1,iCluster2);
    }
    
  }
    
  delete [] mcsAngle2;
  delete [] zMCS;
  delete [] indices;
  
}

  //__________________________________________________________________________
Int_t AliMUONTrack::ClustersInCommon(AliMUONTrack* track, Int_t stMin, Int_t stMax) const
{
  /// Returns the number of clusters in common in stations [stMin, stMax]
  /// between the current track ("this") and the track pointed to by "track".
  
  if (!track || !track->fTrackParamAtCluster || !this->fTrackParamAtCluster) return 0;
  
  Int_t chMin = 2 * stMin;
  Int_t chMax = 2 * stMax + 1;
  Int_t clustersInCommon = 0;
  
  // Loop over clusters of first track
  Int_t nCl1 = this->GetNClusters();
  for(Int_t iCl1 = 0; iCl1 < nCl1; iCl1++) {
    AliMUONVCluster* cl1 = ((AliMUONTrackParam*) this->fTrackParamAtCluster->UncheckedAt(iCl1))->GetClusterPtr();
    
    Int_t chCl1 = cl1->GetChamberId();
    if (chCl1 < chMin || chCl1 > chMax) continue;
    
    // Loop over clusters of second track
    Int_t nCl2 = track->GetNClusters();
    for(Int_t iCl2 = 0; iCl2 < nCl2; iCl2++) {
      AliMUONVCluster* cl2 = ((AliMUONTrackParam*) track->fTrackParamAtCluster->UncheckedAt(iCl2))->GetClusterPtr();
      
      Int_t chCl2 = cl2->GetChamberId();
      if (chCl2 < chMin || chCl2 > chMax) continue;
      
      // Increment "clustersInCommon" if both clusters have the same ID
      if (cl1->GetUniqueID() == cl2->GetUniqueID()) {
	clustersInCommon++;
	break;
      }
      
    }
    
  }
  
  return clustersInCommon;
}

  //__________________________________________________________________________
Int_t AliMUONTrack::GetNDF() const
{
  /// return the number of degrees of freedom
  
  Int_t ndf = 2 * GetNClusters() - 5;
  return (ndf > 0) ? ndf : 0;
}

  //__________________________________________________________________________
Double_t AliMUONTrack::GetNormalizedChi2() const
{
  /// return the chi2 value divided by the number of degrees of freedom (or FLT_MAX if ndf <= 0)
  
  Double_t ndf = (Double_t) GetNDF();
  return (ndf > 0.) ? fGlobalChi2 / ndf : 2.*MaxChi2();
}

  //__________________________________________________________________________
Int_t AliMUONTrack::FindCompatibleClusters(const AliMUONTrack &track, Double_t sigmaCut, Bool_t compatibleCluster[10]) const
{
  /// Try to match clusters from this track with clusters from the given track within the provided sigma cut:
  /// - Fill the array compatibleCluster[iCh] with kTRUE if a compatible cluster has been found in chamber iCh.
  /// - Return the number of clusters of "this" track matched with one cluster of the given track.
  AliMUONVCluster *cluster1, *cluster2;
  Double_t chi2, dX, dY;
  Double_t chi2Max = sigmaCut * sigmaCut;
  
  // initialization
  Int_t nMatchClusters = 0;
  for ( Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) compatibleCluster[ch] = kFALSE;

  if (!track.fTrackParamAtCluster || !this->fTrackParamAtCluster) return nMatchClusters;
  
  // Loop over clusters of first track
  Int_t nCl1 = this->GetNClusters();
  for(Int_t iCl1 = 0; iCl1 < nCl1; iCl1++) {
    cluster1 = static_cast<AliMUONTrackParam*>(this->fTrackParamAtCluster->UncheckedAt(iCl1))->GetClusterPtr();
    
    // Loop over clusters of second track
    Int_t nCl2 = track.GetNClusters();
    for(Int_t iCl2 = 0; iCl2 < nCl2; iCl2++) {
      cluster2 = static_cast<AliMUONTrackParam*>(track.fTrackParamAtCluster->UncheckedAt(iCl2))->GetClusterPtr();
      
      // check DE Id
      if (cluster1->GetDetElemId() != cluster2->GetDetElemId()) continue;
      
      // check local chi2
      dX = cluster1->GetX() - cluster2->GetX();
      dY = cluster1->GetY() - cluster2->GetY();
      chi2 = dX * dX / (cluster1->GetErrX2() + cluster2->GetErrX2()) + dY * dY / (cluster1->GetErrY2() + cluster2->GetErrY2());
      if (chi2 > 2. * chi2Max) continue; // 2 because 2 quantities in chi2
      
      compatibleCluster[cluster1->GetChamberId()] = kTRUE;
      nMatchClusters++;
      break;
    }
    
  }
  
  return nMatchClusters;
}

//__________________________________________________________________________
Bool_t AliMUONTrack::Match(AliMUONTrack &track, Double_t sigmaCut, Int_t &nMatchClusters) const
{
  /// Try to match this track with the given track. Matching conditions:
  /// - more than 50% of clusters from this track matched with clusters from the given track
  /// - at least 1 cluster matched before and 1 cluster matched after the dipole
  
  Bool_t compTrack[10];
  nMatchClusters = FindCompatibleClusters(track, sigmaCut, compTrack);
  
  if ((compTrack[0] || compTrack[1] || compTrack[2] || compTrack[3]) && // at least 1 cluster matched in st 1 & 2
      (compTrack[6] || compTrack[7] || compTrack[8] || compTrack[9]) && // at least 1 cluster matched in st 4 & 5
      2 * nMatchClusters > GetNClusters()) return kTRUE;                // more than 50% of clusters matched
  else return kFALSE;
  
}

//__________________________________________________________________________
void AliMUONTrack::SetTrackParamAtVertex(const AliMUONTrackParam* trackParam)
{
  /// set track parameters at vertex
  if (trackParam == 0x0) return;
  if (fTrackParamAtVertex) *fTrackParamAtVertex = *trackParam;
  else fTrackParamAtVertex = new AliMUONTrackParam(*trackParam);
}

//__________________________________________________________________________
void AliMUONTrack::RecursiveDump() const
{
  /// Recursive dump of AliMUONTrack, i.e. with dump of trackParamAtCluster and attached clusters
  AliMUONTrackParam *trackParamAtCluster;
  AliMUONVCluster *cluster;
  cout << "Recursive dump of Track: " << this << endl;
  // Track
  this->Dump();
  for (Int_t iCluster = 0; iCluster < GetNClusters(); iCluster++) {
    trackParamAtCluster = (AliMUONTrackParam*) ((*fTrackParamAtCluster)[iCluster]);
    // trackParamAtCluster
    cout << "trackParamAtCluster: " << trackParamAtCluster << " (index: " << iCluster << ")" << endl;
    trackParamAtCluster->Dump();
    cluster = trackParamAtCluster->GetClusterPtr();
    // cluster
    cout << "cluster: " << cluster << endl;
    cluster->Print();
  }
  return;
}
  
//_____________________________________________-
void AliMUONTrack::Print(Option_t*) const
{
  /// Printing Track information 

  streamsize curW = cout.width();
  streamsize curPrecision = cout.precision();
  cout << "<AliMUONTrack> No.Clusters=" << setw(2)   << GetNClusters() << 
      ", Match2Trig=" << setw(1) << GetMatchTrigger()  << 
      ", LoTrgNum=" << setw(3) << LoCircuit()  << 
    ", Chi2-tracking-trigger=" << setw(8) << setprecision(5) <<  GetChi2MatchTrigger();
  cout << Form(" HitTriggerPattern %x",fHitsPatternInTrigCh);
  cout << Form(" MClabel=%d",fTrackID) << endl;
  if (fTrackParamAtCluster) fTrackParamAtCluster->First()->Print("FULL");
  cout.width(curW);
  cout.precision(curPrecision);
}

//__________________________________________________________________________
void AliMUONTrack::SetLocalTrigger(Int_t loCirc, Int_t loStripX, Int_t loStripY, Int_t loDev, Int_t loLpt, Int_t loHpt, UChar_t respWithoutChamber)
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
  fLocalTrigger += respWithoutChamber << 26;

}

//__________________________________________________________________________
void AliMUONTrack::FindMCLabel()
{
  /// Determine the MC label from the label of the attached clusters and fill fMCLabel data member:
  /// More than 50% of clusters, including 1 before and 1 after the dipole, must share the same label
  
  Int_t nClusters = GetNClusters();
  Int_t halfCluster = nClusters/2;
  
  // reset MC label
  fTrackID = -1;
  
  // loop over first clusters (if nClusters left < (nClusters-halfCluster) the conditions cannot be fulfilled)
  for (Int_t iCluster1 = 0; iCluster1 < nClusters-halfCluster; iCluster1++) {
    AliMUONVCluster* cluster1 = ((AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(iCluster1))->GetClusterPtr();
    
    // if the first cluster is not on station 1 or 2 the conditions cannot be fulfilled
    if (cluster1->GetChamberId() > 3) return;
    
    Int_t label1 = cluster1->GetMCLabel();
    if (label1 < 0) continue;
    
    Int_t nIdenticalLabel = 1;
    
    // Loop over next clusters
    for (Int_t iCluster2 = iCluster1+1; iCluster2 < nClusters; iCluster2++) {
      AliMUONVCluster* cluster2 = ((AliMUONTrackParam*) fTrackParamAtCluster->UncheckedAt(iCluster2))->GetClusterPtr();
      
      if (cluster2->GetMCLabel() != label1) continue;
      
      nIdenticalLabel++;
      
      // stop as soon as conditions are fulfilled
      if (nIdenticalLabel > halfCluster && cluster2->GetChamberId() > 5) {
        fTrackID = label1;
	return;
      }
      
    }
    
  }
  
}

