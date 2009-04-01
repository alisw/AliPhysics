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

// $Id$

#include "AliMUONTriggerTrackToTrackerClusters.h"

///\class AliMUONTriggerTrackToTrackerClusters
/// 
/// Class to convert trigger tracks into "fake" clusters in stations 4 and 5
///
/// Only intent is to be able to reconstruct data where stations 4 and 5 were
/// not functionning, typically early cosmic runs
///
///\author Laurent Aphecetche, Subatech

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONTrackParam.h"
#include "AliMpArea.h"
#include "AliMpDEManager.h"
#include <TMath.h>

///\cond CLASSIMP
ClassImp(AliMUONTriggerTrackToTrackerClusters)
///\endcond

//_____________________________________________________________________________
AliMUONTriggerTrackToTrackerClusters::AliMUONTriggerTrackToTrackerClusters(const AliMUONGeometryTransformer& transformer,
                                                                           AliMUONVTriggerTrackStore* trackStore)
: TObject(), fkTransformer(transformer), fTriggerTrackStore(trackStore)
{
  /// ctor. We do not take ownership of the trigger track store.
}

//_____________________________________________________________________________
AliMUONTriggerTrackToTrackerClusters::~AliMUONTriggerTrackToTrackerClusters()
{
  /// dtor
}

//_____________________________________________________________________________
Int_t 
AliMUONTriggerTrackToTrackerClusters::DetElemId(Int_t chamber, Double_t x, Double_t y,
                                                Double_t ex, Double_t ey,
                                                Double_t& z) const
{
  /// Find in which detection element (x,y) (global) position is.
  
  AliMpDEIterator it;
  
  AliMpArea a(x, y, ex, ey);
  
  it.First(chamber);
  
  while ( !it.IsDone() )
  {
    Int_t detElemId = it.CurrentDEId();
    
    AliMpArea* area = fkTransformer.GetDEArea(detElemId);
    
    if ( area->Overlap(a) ) 
    {
      // get z of the center of that DE.
      Double_t dummyx, dummyy;
      fkTransformer.Local2Global(detElemId,0,0,0,dummyx,dummyy,z);
      return detElemId;
    }
    it.Next();
  }
  
  return -1;
}

//_____________________________________________________________________________
Int_t 
AliMUONTriggerTrackToTrackerClusters::GenerateClusters(Int_t iChamber,
                                                       AliMUONVClusterStore& clusterStore) const
{
  /// Generate clusters in given chamber
  /// Return the number of clusters added to the clusterStore
  
  AliCodeTimerAuto(Form("Chamber %d",iChamber));
  
  TIter next(fTriggerTrackStore->CreateIterator());
  
  AliMUONTriggerTrack* track;
  Int_t nadded(0);
  
  while ( ( track = static_cast<AliMUONTriggerTrack*>(next()) ) )
  {
    nadded += GenerateClusters(iChamber,*track,clusterStore);
  }
  return nadded;
}

//_____________________________________________________________________________
Int_t 
AliMUONTriggerTrackToTrackerClusters::GenerateClusters(Int_t iChamber,
                                                       const AliMUONTriggerTrack& track,
                                                       AliMUONVClusterStore& clusterStore) const
{
  /// From a trigger track, generate 1 cluster in given chamber
  
  /// Get a (rough) guestimate of the track momentum
  
  Int_t nadded(0);
  
  Double_t z = AliMUONConstants::DefaultChamberZ(10);
  
  Double_t bendingCoord = track.GetY11();
  Double_t bendingSlope = TMath::Tan(track.GetThetay());
  
  Double_t bendingImpact = bendingCoord - z * bendingSlope;
  
  AliDebug(1,Form("TriggerTrack impact parameter=%e",bendingImpact));
  
  //  StdoutToAliDebug(1,track.Print());
  
  Double_t inverseBendingMomentum = 1. / AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(bendingImpact);
  
  // Construct an AliMUONTrackParam from the trigger track, in order to be able to extrapolate it
  // to chambers 6..9 planes.
  
  AliMUONTrackParam trackParam;
  
  trackParam.SetZ(z);
  trackParam.SetNonBendingCoor(track.GetX11());
  trackParam.SetNonBendingSlope(TMath::Tan(track.GetThetax()));
  trackParam.SetBendingCoor(bendingCoord);
  trackParam.SetBendingSlope(bendingSlope);
  trackParam.SetInverseBendingMomentum(inverseBendingMomentum);
  
  Double_t dZ = TMath::Abs(AliMUONConstants::DefaultChamberZ(12) - AliMUONConstants::DefaultChamberZ(10));
  
  Double_t sigmaX = AliMUONConstants::TriggerNonBendingReso();
  
  Double_t sigmaY = AliMUONConstants::TriggerBendingReso();
  
  // Compute and set track parameters covariances
  TMatrixD paramCov(5,5);
  paramCov.Zero();
  
  // Non bending plane
  paramCov(0,0) = sigmaX*sigmaX;
  paramCov(0,1) = -sigmaX/dZ;
  paramCov(1,0) = paramCov(0,1);
  
  paramCov(1,1) = 2.0*sigmaX/dZ/dZ;
  
  // Bending plane
  paramCov(2,2) = sigmaY*sigmaY;
  paramCov(2,3) = -sigmaY/dZ;
  paramCov(3,2) = paramCov(2,3);
  
  paramCov(3,3) = 2.0*sigmaY/dZ/dZ;
  
  // Inverse bending momentum (50% error)
  paramCov(4,4) = 0.5*inverseBendingMomentum * 0.5*inverseBendingMomentum;
  
  // Set covariances
  trackParam.SetCovariances(paramCov);
  
  // Now we extrapolate this trackParam to chambers 6 -> 9
  
  const Float_t kFilterThickness = TMath::Abs(AliMUONConstants::MuonFilterZEnd()-AliMUONConstants::MuonFilterZBeg()); // cm
  
  Int_t nclusters = clusterStore.GetSize();
  
  AliMUONTrackParam tp(trackParam);
  
  Double_t zg = AliMUONConstants::DefaultChamberZ(iChamber);
  AliMUONTrackExtrap::ExtrapToZCov(&tp, zg); // Extrap to iChamber
  
  AliMUONTrackExtrap::AddMCSEffect(&tp, kFilterThickness, AliMUONConstants::MuonFilterX0()); // Add MCS effects
  
  AliDebug(1,Form("iChamber=%d",iChamber));
  
  StdoutToAliDebug(1,tp.Print("FULLCOV"););
  
  Double_t x = tp.GetNonBendingCoor();
  Double_t y = tp.GetBendingCoor();
  const TMatrixD& cov = tp.GetCovariances();
  Double_t ex = TMath::Sqrt(cov(0,0));
  Double_t ey = TMath::Sqrt(cov(2,2));
  
  Double_t zde;
  
  Int_t detElemId = DetElemId(iChamber,x,y,ex,ey,zde);
  
  AliDebug(1,Form("zg = %e zde = %e",zg,zde));
  
  if ( AliMpDEManager::IsValidDetElemId(detElemId) ) 
  {
    AliMUONVCluster* rawCluster = clusterStore.Add(AliMpDEManager::GetChamberId(detElemId), detElemId, nclusters);
    
    ++nclusters;
    ++nadded;
    
    rawCluster->SetCharge(100.0);
    rawCluster->SetXYZ(x, y, zg);
    rawCluster->SetErrXY(ex,ey);
  }
  else
  {
    AliWarning(Form("No DE found at xg=%e yg=%e",detElemId,x,y));
  }
  
  return nadded;
}
