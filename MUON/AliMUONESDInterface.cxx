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

#include "AliMUONESDInterface.h"
#include "AliMUONTrack.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUON2DMapIterator.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONConstants.h"
#include "AliMUONTracker.h"
#include "AliMUONRecoParam.h"
#include "AliMUONVTrackReconstructor.h"

#include "AliMpExMapIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpPad.h"

#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDMuonCluster.h"
#include "AliESDMuonPad.h"
#include "AliLog.h"

#include <TClass.h>
#include <TIterator.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <Riostream.h>

//-----------------------------------------------------------------------------
/// \class AliMUONESDInterface
///
/// There are 2 way of using thid converter between MUON track/cluster/digit
/// and ESDMuon track/cluster/pad:
/// 
/// 1) using the static methods converting the objects one by one
///
/// 2) loading a whole ESDEvent and using the finders and/or the iterators
///    to access the corresponding MUON objects
///
/// note: You can set the recoParam used to refit the MUON track with ResetTracker(...);
///       By default we use Kalman filter + Smoother
///
/// \author Philippe Pillot
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONESDInterface)
/// \endcond

AliMUONRecoParam* AliMUONESDInterface::fgRecoParam = 0x0;
AliMUONVTrackReconstructor* AliMUONESDInterface::fgTracker = 0x0;

TString AliMUONESDInterface::fgTrackStoreName = "AliMUONTrackStoreV1";
TString AliMUONESDInterface::fgClusterStoreName = "AliMUONClusterStoreV2";
TString AliMUONESDInterface::fgDigitStoreName = "AliMUONDigitStoreV2R";
TString AliMUONESDInterface::fgTriggerStoreName = "AliMUONTriggerStoreV1";

//_____________________________________________________________________________
AliMUONESDInterface::AliMUONESDInterface()
: TObject(),
fTracks(0x0),
fDigits(0x0),
fTriggers(0x0),
fClusterMap(0x0),
fDigitMap(0x0)
{
  /// Default constructor
}

//_____________________________________________________________________________
AliMUONESDInterface::~AliMUONESDInterface()
{
  /// Destructor
  delete fTracks;
  delete fDigits;
  delete fTriggers;
  delete fClusterMap;
  delete fDigitMap;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                    methods to play with internal objects                    //
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

//_____________________________________________________________________________
void AliMUONESDInterface::Clear(Option_t*)
{
  /// clear memory
  delete fTracks; fTracks = 0x0;
  delete fDigits; fDigits = 0x0;
  delete fTriggers; fTriggers = 0x0;
  delete fClusterMap; fClusterMap = 0x0;
  delete fDigitMap; fDigitMap = 0x0;
}

//_____________________________________________________________________________
void AliMUONESDInterface::Reset()
{
  /// reset stores and maps
  
  if (fTracks) fTracks->Clear("C");
  else fTracks = NewTrackStore();
  
  if (fDigits) fDigits->Clear("C");
  else fDigits = NewDigitStore();
  
  if (fTriggers) fTriggers->Clear("C");
  else fTriggers = NewTriggerStore();
  
  if (fClusterMap) fClusterMap->Clear();
  else fClusterMap = new AliMpExMap;
  fClusterMap->SetOwner(kTRUE);
  
  if (fDigitMap) fDigitMap->Clear();
  else fDigitMap = new AliMpExMap;
  fDigitMap->SetOwner(kTRUE);
}

//_____________________________________________________________________________
void AliMUONESDInterface::LoadEvent(AliESDEvent& esdEvent)
{
  /// Extract MUON data from the given ESD event
  
  // reset data members
  Reset();
  
  // loop over ESD tracks and fill the stores
  Int_t nTracks = (Int_t) esdEvent.GetNumberOfMuonTracks(); 
  for (Int_t iTrack = 0; iTrack <  nTracks; iTrack++) {
    
    // get ESD track
    AliESDMuonTrack* esdTrack = esdEvent.GetMuonTrack(iTrack);
    
    // fill trigger store if related info are availables
    if (esdTrack->ContainTriggerData()) Add(*esdTrack, *fTriggers);
    
    // fill tracker data if availables
    if (!esdTrack->ContainTrackerData()) continue;
    
    // add it to track store
    AliMUONTrack* track = Add(*esdTrack, *fTracks);
    
    // prepare cluster map
    AliMpExMap* cMap = new AliMpExMap;
    cMap->SetOwner(kFALSE);
    fClusterMap->Add(esdTrack->GetUniqueID(), cMap);
    
    // prepare digit maps
    AliMpExMap* dMaps = new AliMpExMap;
    dMaps->SetOwner(kTRUE);
    fDigitMap->Add(esdTrack->GetUniqueID(), dMaps);
    
    // loop over ESD clusters
    Int_t nClusters = esdTrack->GetNClusters();
    for (Int_t iCluster = 0; iCluster <  nClusters; iCluster++) {
      
      // get ESD cluster
      AliESDMuonCluster *esdCluster = (AliESDMuonCluster*) esdTrack->GetClusters().UncheckedAt(iCluster);
      
      // get the corresponding MUON cluster
      AliMUONVCluster* cluster = FindClusterInTrack(*track, esdCluster->GetUniqueID());
      
      // fill cluster map
      cMap->Add(cluster->GetUniqueID(), cluster);
      
      // prepare digit map
      AliMpExMap* dMap =new AliMpExMap;
      dMap->SetOwner(kFALSE);
      dMaps->Add(esdCluster->GetUniqueID(), dMap);
      
      // loop over ESD pads
      Int_t nPads = esdCluster->GetNPads();
      for (Int_t iPad = 0; iPad < nPads; iPad++) {
	
	// get ESD pad
	AliESDMuonPad *esdPad = (AliESDMuonPad*) esdCluster->GetPads().UncheckedAt(iPad);
	
	// add it to digit store
	AliMUONVDigit* digit = Add(*esdPad, *fDigits);
	
	// fill digit map
	if (digit) dMap->Add(esdPad->GetUniqueID(), digit);
	else dMap->Add(esdPad->GetUniqueID(), fDigits->FindObject(esdPad->GetUniqueID()));
	
      } // end of loop over pads
      
    } // end of loop over clusters
    
  } // end of loop over tracks
  
}

//___________________________________________________________________________
Int_t AliMUONESDInterface::GetNTracks() const
{
  /// return the number of tracks
  return fTracks ? fTracks->GetSize() : 0;
}

//___________________________________________________________________________
Int_t AliMUONESDInterface::GetNClusters() const
{
  /// return the number of clusters
  Int_t nClusters = 0;
  AliMUONTrack *track;
  TIter next(CreateTrackIterator());
  while ((track = static_cast<AliMUONTrack*>(next()))) nClusters += track->GetNClusters();
  return nClusters;
}

//___________________________________________________________________________
Int_t AliMUONESDInterface::GetNClusters(UInt_t trackId) const
{
  /// return the number of clusters in track "trackId"
  AliMUONTrack* track = FindTrack(trackId);
  return track ? track->GetNClusters() : 0;
}

//___________________________________________________________________________
Int_t AliMUONESDInterface::GetNDigits() const
{
  /// return the number of digits
  return fDigits ? fDigits->GetSize() : 0;
}

//___________________________________________________________________________
Int_t AliMUONESDInterface::GetNDigits(UInt_t trackId) const
{
  /// return the number of digits in all clusters of track "trackId"
  Int_t nDigits = 0;
  AliMUONVCluster *cluster;
  TIter next(CreateClusterIterator(trackId));
  while ((cluster = static_cast<AliMUONVCluster*>(next()))) nDigits += cluster->GetNDigits();
  return nDigits;
}

//___________________________________________________________________________
Int_t AliMUONESDInterface::GetNDigits(UInt_t trackId, UInt_t clusterId) const
{
  /// return the number of digits in cluster numbered "iCluster" of track "iTrack"
  AliMUONVCluster* cluster = FindCluster(trackId, clusterId);
  return cluster ? cluster->GetNDigits() : 0;
}

//___________________________________________________________________________
Int_t AliMUONESDInterface::GetNDigitsInCluster(UInt_t clusterId) const
{
  /// return the number of digits in cluster "clusterId"
  AliMUONVCluster* cluster = FindCluster(clusterId);
  return cluster ? cluster->GetNDigits() : 0;
}

//___________________________________________________________________________
Int_t AliMUONESDInterface::GetNTriggers() const
{
  /// return the number of triggers
  return fTriggers ? fTriggers->GetSize() : 0;
}

//___________________________________________________________________________
AliMUONTrack* AliMUONESDInterface::FindTrack(UInt_t trackId) const
{
  /// return track "trackId" (0x0 if not found)
  AliMUONTrack *track = fTracks ? static_cast<AliMUONTrack*>(fTracks->FindObject(trackId)) : 0x0;
  if (!track) AliWarning(Form("track %d does not exist",trackId));
  return track;
}

//___________________________________________________________________________
AliMUONVCluster* AliMUONESDInterface::FindCluster(UInt_t clusterId) const
{
  /// loop over tracks and return the first cluster "clusterId" found (0x0 if not found)
  AliMpExMap *cMap;
  AliMUONVCluster* cluster = 0x0;
  
  if (fClusterMap) {
    
    TIter next(fClusterMap->CreateIterator());
    while ((cMap = static_cast<AliMpExMap*>(next()))) {
      
      cluster = static_cast<AliMUONVCluster*>(cMap->GetValue(clusterId));
      if (cluster) return cluster;
      
    }
    
  }
  
  if (!cluster) AliWarning(Form("cluster %d does not exist",clusterId));
  return 0x0;
}

//___________________________________________________________________________
AliMUONVCluster* AliMUONESDInterface::FindCluster(UInt_t trackId, UInt_t clusterId) const
{
  /// return cluster "clusterId" in track "trackId" (0x0 if not found)
  AliMpExMap *cMap = fClusterMap ? static_cast<AliMpExMap*>(fClusterMap->GetValue(trackId)) : 0x0;
  AliMUONVCluster* cluster = cMap ? static_cast<AliMUONVCluster*>(cMap->GetValue(clusterId)) : 0x0;
  if (!cluster) AliWarning(Form("cluster %d does not exist in track %d", clusterId, trackId));
  return cluster;
}

//___________________________________________________________________________
AliMUONVDigit* AliMUONESDInterface::FindDigit(UInt_t digitId) const
{
  /// return digit "digitId" (0x0 if not found)
  AliMUONVDigit *digit = fDigits ? fDigits->FindObject(digitId) : 0x0;
  if (!digit) AliWarning(Form("digit %d does not exist",digitId));
  return digit;
}

//___________________________________________________________________________
AliMUONLocalTrigger* AliMUONESDInterface::FindLocalTrigger(Int_t boardNumber) const
{
  /// return MUON local trigger "boardNumber"
  return (fTriggers) ? fTriggers->FindLocal(boardNumber) : 0x0;
}

//___________________________________________________________________________
TIterator* AliMUONESDInterface::CreateTrackIterator() const
{
  /// return iterator over all tracks
  return fTracks ? fTracks->CreateIterator() : 0x0;
}

//___________________________________________________________________________
TIterator* AliMUONESDInterface::CreateClusterIterator() const
{
  /// return iterator over all clusters
  return fClusterMap ? new AliMUON2DMapIterator(*fClusterMap) : 0x0;
}

//___________________________________________________________________________
TIterator* AliMUONESDInterface::CreateClusterIterator(UInt_t trackId) const
{
  /// return iterator over clusters of track "trackId"
  AliMpExMap *cMap = fClusterMap ? static_cast<AliMpExMap*>(fClusterMap->GetValue(trackId)) : 0x0;
  return cMap ? cMap->CreateIterator() : 0x0;
}

//___________________________________________________________________________
TIterator* AliMUONESDInterface::CreateDigitIterator() const
{
  /// return iterator over all digits
  return fDigits ? fDigits->CreateIterator() : 0x0;
}

//___________________________________________________________________________
TIterator* AliMUONESDInterface::CreateDigitIterator(UInt_t trackId) const
{
  /// return iterator over all digits of track "trackId"
  AliMpExMap* dMaps = fDigitMap ? static_cast<AliMpExMap*>(fDigitMap->GetValue(trackId)) : 0x0;
  return dMaps ? new AliMUON2DMapIterator(*dMaps) : 0x0;
}

//___________________________________________________________________________
TIterator* AliMUONESDInterface::CreateDigitIterator(UInt_t trackId, UInt_t clusterId) const
{
  /// return iterator over digits of cluster "clusterId" in track "trackId"
  AliMpExMap* dMaps = fDigitMap ? static_cast<AliMpExMap*>(fDigitMap->GetValue(trackId)) : 0x0;
  AliMpExMap* dMap = dMaps ? static_cast<AliMpExMap*>(dMaps->GetValue(clusterId)) : 0x0;
  return dMap ? dMap->CreateIterator() : 0x0;
}

//___________________________________________________________________________
TIterator* AliMUONESDInterface::CreateDigitIteratorInCluster(UInt_t clusterId) const
{
  /// return iterator over digits of the first cluster "clusterId" found by looping over all tracks
  AliMpExMap *dMaps;
  AliMpExMap* dMap = 0x0;
  
  if (fDigitMap) {
    
    TIter next(fDigitMap->CreateIterator());
    while ((dMaps = static_cast<AliMpExMap*>(next()))) {
      
      dMap = static_cast<AliMpExMap*>(dMaps->GetValue(clusterId));
      if (dMap) return dMap->CreateIterator();
      
    }
    
  }
  
  return 0x0;
}

//___________________________________________________________________________
TIterator* AliMUONESDInterface::CreateLocalTriggerIterator() const
{
  /// return iterator over all local trigger
  return fTriggers ? fTriggers->CreateLocalIterator() : 0x0;
}

//___________________________________________________________________________
AliMUONVCluster* AliMUONESDInterface::FindClusterInTrack(const AliMUONTrack& track, UInt_t clusterId) const
{
  /// find the cluster with the given Id into the track
  
  Int_t nClusters = track.GetNClusters();
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
    
    AliMUONVCluster* cluster = ((AliMUONTrackParam*) track.GetTrackParamAtCluster()->UncheckedAt(iCluster))->GetClusterPtr();
    if (cluster->GetUniqueID() == clusterId) return cluster;
    
  }
  
  return 0x0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                                static methods                               //
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

//_____________________________________________________________________________
void AliMUONESDInterface::ResetTracker(const AliMUONRecoParam* recoParam)
{
  /// Reset the MUON tracker using "recoParam" if provided.
  /// If not provided, will use Kalman filter + Smoother
  
  delete fgTracker;
  delete fgRecoParam;
  
  if (recoParam) {
    
    fgRecoParam = new AliMUONRecoParam(*recoParam);
    
  } else {
    
    fgRecoParam = AliMUONRecoParam::GetLowFluxParam();
    fgRecoParam->SetTrackingMode("KALMAN");
    fgRecoParam->UseSmoother(kTRUE);
    fgRecoParam->SetBendingVertexDispersion(10.);
    
  }
  
  fgTracker = AliMUONTracker::CreateTrackReconstructor(fgRecoParam,0x0);
  
}

//_____________________________________________________________________________
AliMUONVTrackStore* AliMUONESDInterface::NewTrackStore()
{
  /// Create an empty track store of type fgTrackStoreName
  TClass* classPtr = TClass::GetClass(fgTrackStoreName);
  if (!classPtr || !classPtr->InheritsFrom("AliMUONVTrackStore")) {
    cout<<"E-AliMUONESDInterface::NewTrackStore: Unable to create store of type "<<fgTrackStoreName.Data()<<endl;
    return 0x0;
  }
  return reinterpret_cast<AliMUONVTrackStore*>(classPtr->New());
}

//_____________________________________________________________________________
AliMUONVClusterStore* AliMUONESDInterface::NewClusterStore()
{
  /// Create an empty cluster store of type fgClusterStoreName
  TClass* classPtr = TClass::GetClass(fgClusterStoreName);
  if (!classPtr || !classPtr->InheritsFrom("AliMUONVClusterStore")) {
    cout<<"E-AliMUONESDInterface::NewClusterStore: Unable to create store of type "<<fgClusterStoreName.Data()<<endl;
    return 0x0;
  }
  return reinterpret_cast<AliMUONVClusterStore*>(classPtr->New());
}

//_____________________________________________________________________________
AliMUONVDigitStore* AliMUONESDInterface::NewDigitStore()
{
  /// Create an empty digit store of type fgDigitStoreName
  TClass* classPtr = TClass::GetClass(fgDigitStoreName);
  if (!classPtr || !classPtr->InheritsFrom("AliMUONVDigitStore")) {
    cout<<"E-AliMUONESDInterface::NewDigitStore: Unable to create store of type "<<fgDigitStoreName.Data()<<endl;
    return 0x0;
  }
  return reinterpret_cast<AliMUONVDigitStore*>(classPtr->New());
}

//_____________________________________________________________________________
AliMUONVTriggerStore* AliMUONESDInterface::NewTriggerStore()
{
  /// Create an empty trigger store of type fgTriggerStoreName
  TClass* classPtr = TClass::GetClass(fgTriggerStoreName);
  if (!classPtr || !classPtr->InheritsFrom("AliMUONVTriggerStore")) {
    cout<<"E-AliMUONESDInterface::NewTriggerStore: Unable to create store of type "<<fgTriggerStoreName.Data()<<endl;
    return 0x0;
  }
  return reinterpret_cast<AliMUONVTriggerStore*>(classPtr->New());
}

//_________________________________________________________________________
void AliMUONESDInterface::GetParamAtVertex(const AliESDMuonTrack& esdTrack, AliMUONTrackParam& trackParam)
{
  /// Get parameters at vertex from ESDMuon track
  trackParam.SetZ(esdTrack.GetZ());
  trackParam.SetNonBendingCoor(esdTrack.GetNonBendingCoor());
  trackParam.SetNonBendingSlope(TMath::Tan(esdTrack.GetThetaX()));
  trackParam.SetBendingCoor(esdTrack.GetBendingCoor());
  trackParam.SetBendingSlope(TMath::Tan(esdTrack.GetThetaY()));
  trackParam.SetInverseBendingMomentum(esdTrack.GetInverseBendingMomentum());
}

//_________________________________________________________________________
void AliMUONESDInterface::GetParamAtDCA(const AliESDMuonTrack& esdTrack, AliMUONTrackParam& trackParam)
{
  /// Get parameters at DCA from ESDMuon track
  trackParam.SetZ(esdTrack.GetZ());
  trackParam.SetNonBendingCoor(esdTrack.GetNonBendingCoorAtDCA());
  trackParam.SetNonBendingSlope(TMath::Tan(esdTrack.GetThetaXAtDCA()));
  trackParam.SetBendingCoor(esdTrack.GetBendingCoorAtDCA());
  trackParam.SetBendingSlope(TMath::Tan(esdTrack.GetThetaYAtDCA()));
  trackParam.SetInverseBendingMomentum(esdTrack.GetInverseBendingMomentumAtDCA());
}

//_________________________________________________________________________
void AliMUONESDInterface::GetParamAtFirstCluster(const AliESDMuonTrack& esdTrack, AliMUONTrackParam& trackParam)
{
  /// Get parameters at first cluster from ESDMuon track
  trackParam.SetZ(esdTrack.GetZUncorrected());
  trackParam.SetNonBendingCoor(esdTrack.GetNonBendingCoorUncorrected());
  trackParam.SetNonBendingSlope(TMath::Tan(esdTrack.GetThetaXUncorrected()));
  trackParam.SetBendingCoor(esdTrack.GetBendingCoorUncorrected());
  trackParam.SetBendingSlope(TMath::Tan(esdTrack.GetThetaYUncorrected()));
  trackParam.SetInverseBendingMomentum(esdTrack.GetInverseBendingMomentumUncorrected());
}

//_________________________________________________________________________
void AliMUONESDInterface::GetParamCov(const AliESDMuonTrack& esdTrack, AliMUONTrackParam& trackParam)
{
  /// Get parameters covariances from ESD track
  
  // get ESD covariance matrix
  TMatrixD covariances(5,5);
  esdTrack.GetCovariances(covariances);
  
  // compute Jacobian to change the coordinate system
  // from (X,thetaX,Y,thetaY,c/pYZ) to (X,slopeX,Y,slopeY,c/pYZ)
  Double_t cosThetaX = TMath::Cos(esdTrack.GetThetaXUncorrected());
  Double_t cosThetaY = TMath::Cos(esdTrack.GetThetaYUncorrected());
  TMatrixD jacob(5,5);
  jacob.Zero();
  jacob(0,0) = 1.;
  jacob(1,1) = 1. / cosThetaX / cosThetaX;
  jacob(2,2) = 1.;
  jacob(3,3) = 1. / cosThetaY / cosThetaY;
  jacob(4,4) = 1.;
  
  // compute covariance matrix in ESD coordinate system
  TMatrixD tmp(covariances,TMatrixD::kMultTranspose,jacob);
  trackParam.SetCovariances(TMatrixD(jacob,TMatrixD::kMult,tmp));
  
}

//_________________________________________________________________________
void AliMUONESDInterface::SetParamAtVertex(const AliMUONTrackParam& trackParam, AliESDMuonTrack& esdTrack)
{
  /// Set parameters in ESD track
  esdTrack.SetZ(trackParam.GetZ());
  esdTrack.SetNonBendingCoor(trackParam.GetNonBendingCoor());
  esdTrack.SetThetaX(TMath::ATan(trackParam.GetNonBendingSlope()));
  esdTrack.SetBendingCoor(trackParam.GetBendingCoor()); 
  esdTrack.SetThetaY(TMath::ATan(trackParam.GetBendingSlope()));
  esdTrack.SetInverseBendingMomentum(trackParam.GetInverseBendingMomentum());
}

//_________________________________________________________________________
void AliMUONESDInterface::SetParamAtDCA(const AliMUONTrackParam& trackParam, AliESDMuonTrack& esdTrack)
{
  /// Set parameters in ESD track
  esdTrack.SetNonBendingCoorAtDCA(trackParam.GetNonBendingCoor());
  esdTrack.SetThetaXAtDCA(TMath::ATan(trackParam.GetNonBendingSlope()));
  esdTrack.SetBendingCoorAtDCA(trackParam.GetBendingCoor()); 
  esdTrack.SetThetaYAtDCA(TMath::ATan(trackParam.GetBendingSlope()));
  esdTrack.SetInverseBendingMomentumAtDCA(trackParam.GetInverseBendingMomentum());
}

//_________________________________________________________________________
void AliMUONESDInterface::SetParamAtFirstCluster(const AliMUONTrackParam& trackParam, AliESDMuonTrack& esdTrack)
{
  /// Set parameters in ESD track
  esdTrack.SetZUncorrected(trackParam.GetZ());
  esdTrack.SetNonBendingCoorUncorrected(trackParam.GetNonBendingCoor());
  esdTrack.SetThetaXUncorrected(TMath::ATan(trackParam.GetNonBendingSlope()));
  esdTrack.SetBendingCoorUncorrected(trackParam.GetBendingCoor()); 
  esdTrack.SetThetaYUncorrected(TMath::ATan(trackParam.GetBendingSlope()));
  esdTrack.SetInverseBendingMomentumUncorrected(trackParam.GetInverseBendingMomentum());
}

//_________________________________________________________________________
void AliMUONESDInterface::SetParamCov(const AliMUONTrackParam& trackParam, AliESDMuonTrack& esdTrack)
{
  /// Set parameters covariances in ESD track
  
  // set null matrix if covariances does not exist
  if (!trackParam.CovariancesExist()) {
    TMatrixD tmp(5,5);
    tmp.Zero();
    esdTrack.SetCovariances(tmp);
    return;
  }
  
  // compute Jacobian to change the coordinate system
  // from (X,slopeX,Y,slopeY,c/pYZ) to (X,thetaX,Y,thetaY,c/pYZ)
  Double_t cosThetaX = TMath::Cos(TMath::ATan(trackParam.GetNonBendingSlope()));
  Double_t cosThetaY = TMath::Cos(TMath::ATan(trackParam.GetBendingSlope()));
  TMatrixD jacob(5,5);
  jacob.Zero();
  jacob(0,0) = 1.;
  jacob(1,1) = cosThetaX * cosThetaX;
  jacob(2,2) = 1.;
  jacob(3,3) = cosThetaY * cosThetaY;
  jacob(4,4) = 1.;
  
  // compute covariance matrix in ESD coordinate system
  TMatrixD tmp(trackParam.GetCovariances(),TMatrixD::kMultTranspose,jacob);
  esdTrack.SetCovariances(TMatrixD(jacob,TMatrixD::kMult,tmp));
  
}

//_____________________________________________________________________________
void AliMUONESDInterface::ESDToMUON(const AliESDMuonTrack& esdTrack, AliMUONTrack& track)
{
  /// Transfert data from ESDMuon track to MUON track.
  /// The track parameters at each cluster are obtained by refitting the track
  /// or by extrapolating the parameters at the first one if the refit failed.
  /// note: You can set the recoParam used to refit the MUON track with ResetTracker(...);
  ///       By default we use Kalman filter + Smoother
  
  // if the ESDMuon track is a ghost then return an empty MUON track
  if (!esdTrack.ContainTrackerData()) {
    track.Reset();
    track.SetUniqueID(esdTrack.GetUniqueID());
    return;
  }
  
  track.Clear("C");
  
  // global info
  track.SetUniqueID(esdTrack.GetUniqueID());
  track.FitWithVertex(kFALSE);
  track.FitWithMCS(kFALSE);
  track.SetImproved(kFALSE);
  track.SetVertexErrXY2(0.,0.);
  track.SetGlobalChi2(esdTrack.GetChi2());
  track.SetMatchTrigger(esdTrack.GetMatchTrigger());
  track.SetLoTrgNum(-1);
  track.SetChi2MatchTrigger(esdTrack.GetChi2MatchTrigger());
  track.SetHitsPatternInTrigCh(esdTrack.GetHitsPatternInTrigCh());
  track.SetLocalTrigger(esdTrack.LoCircuit(), esdTrack.LoStripX(), esdTrack.LoStripY(),
			esdTrack.LoDev(), esdTrack.LoLpt(), esdTrack.LoHpt());
  
  // track parameters at vertex
  AliMUONTrackParam paramAtVertex;
  GetParamAtVertex(esdTrack, paramAtVertex);
  track.SetTrackParamAtVertex(&paramAtVertex);
    
  // track parameters at first cluster
  AliMUONTrackParam param;
  GetParamAtFirstCluster(esdTrack, param);
  GetParamCov(esdTrack, param);
  
  // create empty cluster
  AliMUONVClusterStore* cStore = NewClusterStore();
  if (!cStore) return;
  AliMUONVCluster* cluster = cStore->CreateCluster(0,0,0);
  
  // fill TrackParamAtCluster with track parameters at each cluster if available
  // or with only track parameters at first (fake) cluster if not
  if(esdTrack.ClustersStored()) {
    
    // loop over ESD clusters
    AliESDMuonCluster *esdCluster = (AliESDMuonCluster*) esdTrack.GetClusters().First();
    while (esdCluster) {
      
      // copy cluster information
      ESDToMUON(*esdCluster, *cluster);
      
      // only set the Z parameter to avoid error in the AddTrackParamAtCluster(...) method
      param.SetZ(cluster->GetZ());
      
      // add common track parameters at current cluster
      track.AddTrackParamAtCluster(param, *cluster, kTRUE);
      
      esdCluster = (AliESDMuonCluster*) esdTrack.GetClusters().After(esdCluster);
    }
    
    // recompute parameters at first cluster in case of those stored
    // in ESD are not related to the most upstream cluster
    AliMUONTrackParam *firstTrackParam = (AliMUONTrackParam*) track.GetTrackParamAtCluster()->First();
    firstTrackParam->SetZ(esdTrack.GetZUncorrected()); // reset the z to the one stored in ESD
    AliMUONTrackExtrap::ExtrapToZCov(firstTrackParam,firstTrackParam->GetClusterPtr()->GetZ());
    
    // refit the track to get better parameters and covariances at each cluster (temporary disable track improvement)
    if (!fgTracker) ResetTracker();
    if (!fgTracker->RefitTrack(track, kFALSE)) track.UpdateCovTrackParamAtCluster();
    
  } else {
    
    // get number of the first hit chamber according to the MUONClusterMap
    Int_t firstCh = 0;
    while (firstCh < 10 && !esdTrack.IsInMuonClusterMap(firstCh)) firstCh++;
    
    // produce fake cluster at this chamber
    cluster->SetUniqueID(AliMUONVCluster::BuildUniqueID(firstCh, 0, 0));
    cluster->SetXYZ(param.GetNonBendingCoor(), param.GetBendingCoor(), param.GetZ());
    cluster->SetErrXY(0., 0.);
    
    // add track parameters at first (fake) cluster
    track.AddTrackParamAtCluster(param, *cluster, kTRUE);
    
  }
  
  // set the MC label from ESD track
  track.SetMCLabel(esdTrack.GetLabel());
  
  delete cluster;
  delete cStore;
  
}

//_____________________________________________________________________________
void AliMUONESDInterface::ESDToMUON(const AliESDMuonTrack& esdTrack, AliMUONLocalTrigger& locTrg)
{
  /// Transfert trigger data from ESDMuon track to the MUONLocalTtrigger object
  
  // if the ESDMuon track is a ghost then return an empty MUON track
  if (!esdTrack.ContainTriggerData()) {
    AliMUONLocalTrigger emptyLocTrg;
    locTrg = emptyLocTrg;
    return;
  }
  
  locTrg.SetLoCircuit(esdTrack.LoCircuit());
  locTrg.SetLoStripX(esdTrack.LoStripX());
  locTrg.SetLoStripY(esdTrack.LoStripY());
  locTrg.SetDeviation(esdTrack.LoDev());
  locTrg.SetLoLpt(esdTrack.LoLpt());
  locTrg.SetLoHpt(esdTrack.LoHpt());
  locTrg.SetLoTrigY(1);
  locTrg.SetX1Pattern(esdTrack.GetTriggerX1Pattern());
  locTrg.SetX2Pattern(esdTrack.GetTriggerX2Pattern());
  locTrg.SetX3Pattern(esdTrack.GetTriggerX3Pattern());
  locTrg.SetX4Pattern(esdTrack.GetTriggerX4Pattern());
  locTrg.SetY1Pattern(esdTrack.GetTriggerY1Pattern());
  locTrg.SetY2Pattern(esdTrack.GetTriggerY2Pattern());
  locTrg.SetY3Pattern(esdTrack.GetTriggerY3Pattern());
  locTrg.SetY4Pattern(esdTrack.GetTriggerY4Pattern());
  
}

//_____________________________________________________________________________
void AliMUONESDInterface::ESDToMUON(const AliESDMuonCluster& esdCluster, AliMUONVCluster& cluster)
{
  /// Transfert data from ESDMuon cluster to MUON cluster
  
  cluster.Clear("C");
  
  cluster.SetUniqueID(esdCluster.GetUniqueID());
  cluster.SetXYZ(esdCluster.GetX(), esdCluster.GetY(), esdCluster.GetZ());
  cluster.SetErrXY(esdCluster.GetErrX(),esdCluster.GetErrY());
  cluster.SetCharge(esdCluster.GetCharge());
  cluster.SetChi2(esdCluster.GetChi2());
  cluster.SetMCLabel(esdCluster.GetLabel());
  
  if (esdCluster.PadsStored()) {
    Int_t nPads = esdCluster.GetNPads();
    for (Int_t iPad = 0; iPad < nPads; iPad++)
      cluster.AddDigitId(((AliESDMuonPad*)esdCluster.GetPads().UncheckedAt(iPad))->GetUniqueID());
  }
  
}

//___________________________________________________________________________
void AliMUONESDInterface::ESDToMUON(const AliESDMuonPad& esdPad, AliMUONVDigit& digit)
{
  /// Transfert data from ESDMuon pad to MUON digit
  
  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(esdPad.GetDetElemId(), esdPad.GetManuId());  
  AliMpPad pad = seg->PadByLocation(esdPad.GetManuId(), esdPad.GetManuChannel(), kFALSE);
  
  digit.Saturated(esdPad.IsSaturated());
  digit.Used(kFALSE);
  digit.Calibrated(esdPad.IsCalibrated());
  digit.SetUniqueID(esdPad.GetUniqueID());
  digit.SetCharge(esdPad.GetCharge());
  digit.SetADC(esdPad.GetADC());
  digit.SetPadXY(pad.GetIx(), pad.GetIy());
  
}

//_____________________________________________________________________________
void AliMUONESDInterface::MUONToESD(const AliMUONTrack& track, AliESDMuonTrack& esdTrack, const Double_t vertex[3],
				    const AliMUONVDigitStore* digits, const AliMUONLocalTrigger* locTrg)
{
  /// Transfert data from MUON track to ESDMuon track
  /// Incorporate the ESDPads if the digits are provided
  /// Add trigger info if the MUON track is matched with a trigger track
  
  // empty MUON track -> produce a ghost ESDMuon track if trigger info are available otherwise produce an empty track
  if (track.GetNClusters() == 0) {
    if (locTrg) MUONToESD(*locTrg, esdTrack, track.GetUniqueID());
    else {
      cout<<"W-AliMUONESDInterface::MUONToESD: will produce an empty ESDMuon track"<<endl;
      esdTrack.Reset();
      esdTrack.SetUniqueID(0xFFFFFFFF);
    }
    return;
  }
  
  esdTrack.Clear("C");
  
  // set global info
  esdTrack.SetUniqueID(track.GetUniqueID());
  esdTrack.SetChi2(track.GetGlobalChi2());
  esdTrack.SetNHit(track.GetNClusters());
  esdTrack.SetLabel(track.GetMCLabel());
  
  // set param at first cluster
  AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>((track.GetTrackParamAtCluster())->First());
  SetParamAtFirstCluster(*trackParam, esdTrack);
  SetParamCov(*trackParam, esdTrack);
  
  // set param at vertex
  AliMUONTrackParam trackParamAtVtx(*trackParam);
  AliMUONTrackExtrap::ExtrapToVertex(&trackParamAtVtx, vertex[0], vertex[1], vertex[2], 0., 0.);
  SetParamAtVertex(trackParamAtVtx, esdTrack);
  
  // set param at Distance of Closest Approach
  AliMUONTrackParam trackParamAtDCA(*trackParam);
  AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&trackParamAtDCA, vertex[2]);
  SetParamAtDCA(trackParamAtDCA, esdTrack);
  
  // set muon cluster info
  AliESDMuonCluster esdCluster;
  esdTrack.SetMuonClusterMap(0);
  while (trackParam) {
    MUONToESD(*(trackParam->GetClusterPtr()), esdCluster, digits);
    esdTrack.AddCluster(esdCluster);
    esdTrack.AddInMuonClusterMap(esdCluster.GetChamberId());
    trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
  }
  
  // set trigger info
  esdTrack.SetLocalTrigger(track.GetLocalTrigger());
  esdTrack.SetChi2MatchTrigger(track.GetChi2MatchTrigger());
  esdTrack.SetHitsPatternInTrigCh(track.GetHitsPatternInTrigCh());
  if (locTrg) {
    esdTrack.SetTriggerX1Pattern(locTrg->GetX1Pattern());
    esdTrack.SetTriggerY1Pattern(locTrg->GetY1Pattern());
    esdTrack.SetTriggerX2Pattern(locTrg->GetX2Pattern());
    esdTrack.SetTriggerY2Pattern(locTrg->GetY2Pattern());
    esdTrack.SetTriggerX3Pattern(locTrg->GetX3Pattern());
    esdTrack.SetTriggerY3Pattern(locTrg->GetY3Pattern());
    esdTrack.SetTriggerX4Pattern(locTrg->GetX4Pattern());
    esdTrack.SetTriggerY4Pattern(locTrg->GetY4Pattern());
  } else {
    esdTrack.SetTriggerX1Pattern(0);
    esdTrack.SetTriggerY1Pattern(0);
    esdTrack.SetTriggerX2Pattern(0);
    esdTrack.SetTriggerY2Pattern(0);
    esdTrack.SetTriggerX3Pattern(0);
    esdTrack.SetTriggerY3Pattern(0);
    esdTrack.SetTriggerX4Pattern(0);
    esdTrack.SetTriggerY4Pattern(0);
  }
  
}

//_____________________________________________________________________________
void AliMUONESDInterface::MUONToESD(const AliMUONLocalTrigger& locTrg, AliESDMuonTrack& esdTrack,
				    UInt_t trackId, const AliMUONTriggerTrack* triggerTrack)
{
  /// Build ghost ESDMuon track containing only informations about trigger track
  
  esdTrack.Reset();
  esdTrack.SetUniqueID(trackId);
  
  // set trigger info
  AliMUONTrack muonTrack;
  muonTrack.SetLocalTrigger(locTrg.LoCircuit(),
			    locTrg.LoStripX(),
			    locTrg.LoStripY(),
			    locTrg.GetDeviation(),
			    locTrg.LoLpt(),
			    locTrg.LoHpt());
  esdTrack.SetLocalTrigger(muonTrack.GetLocalTrigger());
  esdTrack.SetChi2MatchTrigger(0.);
  esdTrack.SetTriggerX1Pattern(locTrg.GetX1Pattern());
  esdTrack.SetTriggerY1Pattern(locTrg.GetY1Pattern());
  esdTrack.SetTriggerX2Pattern(locTrg.GetX2Pattern());
  esdTrack.SetTriggerY2Pattern(locTrg.GetY2Pattern());
  esdTrack.SetTriggerX3Pattern(locTrg.GetX3Pattern());
  esdTrack.SetTriggerY3Pattern(locTrg.GetY3Pattern());
  esdTrack.SetTriggerX4Pattern(locTrg.GetX4Pattern());
  esdTrack.SetTriggerY4Pattern(locTrg.GetY4Pattern());
  UShort_t hitPattern = 0;
  if(triggerTrack){
    hitPattern = triggerTrack->GetHitsPatternInTrigCh();
    esdTrack.SetHitsPatternInTrigCh(hitPattern);
    esdTrack.SetThetaXUncorrected(triggerTrack->GetThetax());
    esdTrack.SetThetaYUncorrected(triggerTrack->GetThetay());
    esdTrack.SetNonBendingCoorUncorrected(triggerTrack->GetX11());
    esdTrack.SetBendingCoorUncorrected(triggerTrack->GetY11());
  }
}

//_____________________________________________________________________________
void AliMUONESDInterface::MUONToESD(const AliMUONVCluster& cluster, AliESDMuonCluster& esdCluster, const AliMUONVDigitStore* digits)
{
  /// Transfert data from MUON cluster to ESDMuon cluster
  /// Incorporate the ESDPads if the digits are provided
  
  esdCluster.Clear("C");
  
  esdCluster.SetUniqueID(cluster.GetUniqueID());
  esdCluster.SetXYZ(cluster.GetX(), cluster.GetY(), cluster.GetZ());
  esdCluster.SetErrXY(cluster.GetErrX(), cluster.GetErrY());
  esdCluster.SetCharge(cluster.GetCharge());
  esdCluster.SetChi2(cluster.GetChi2());
  esdCluster.SetLabel(cluster.GetMCLabel());
  
  if (digits) { // transfert all data if required
    
    AliESDMuonPad esdPad;
    for (Int_t i=0; i<cluster.GetNDigits(); i++) {
      AliMUONVDigit* digit = digits->FindObject(cluster.GetDigitId(i));
      if (!digit) {
	cout<<"E-AliMUONESDInterface::MUONToESD: digit "<<cluster.GetDigitId(i)<<" not found"<<endl;
	continue;
      }
      MUONToESD(*digit, esdPad);
      esdCluster.AddPad(esdPad);
    }
    
  }
  
}

//_____________________________________________________________________________
void AliMUONESDInterface::MUONToESD(const AliMUONVDigit& digit, AliESDMuonPad& esdPad)
{
  /// Transfert data from MUON digit to ESDMuon pad
  esdPad.SetUniqueID(digit.GetUniqueID());
  esdPad.SetADC(digit.ADC());
  esdPad.SetCharge(digit.Charge());
  esdPad.SetCalibrated(digit.IsCalibrated());
  esdPad.SetSaturated(digit.IsSaturated());
}

//___________________________________________________________________________
AliMUONTrack* AliMUONESDInterface::Add(const AliESDMuonTrack& esdTrack, AliMUONVTrackStore& trackStore)
{
  /// Create MUON track from ESDMuon track and add it to the store
  /// return a pointer to the track into the store (0x0 if the track already exist)
  if(trackStore.FindObject(esdTrack.GetUniqueID())) return 0x0;
  AliMUONTrack* track = trackStore.Add(AliMUONTrack());
  ESDToMUON(esdTrack, *track);
  return track;
}

//___________________________________________________________________________
void AliMUONESDInterface::Add(const AliESDMuonTrack& esdTrack, AliMUONVTriggerStore& triggerStore)
{
  /// Create MUON local trigger from ESDMuon track and add it to the store if not already there
  if (triggerStore.FindLocal(esdTrack.LoCircuit())) return;
  AliMUONLocalTrigger locTrg;
  ESDToMUON(esdTrack, locTrg);
  triggerStore.Add(locTrg);
}

//___________________________________________________________________________
AliMUONVCluster* AliMUONESDInterface::Add(const AliESDMuonCluster& esdCluster, AliMUONVClusterStore& clusterStore)
{
  /// Create MUON cluster from ESDMuon cluster and add it to the store
  /// return a pointer to the cluster into the store (0x0 if the cluster already exist)
  AliMUONVCluster* cluster = clusterStore.Add(esdCluster.GetChamberId(), esdCluster.GetDetElemId(), esdCluster.GetClusterIndex());
  if (cluster) ESDToMUON(esdCluster, *cluster);
  return cluster;
}

//___________________________________________________________________________
AliMUONVDigit* AliMUONESDInterface::Add(const AliESDMuonPad& esdPad, AliMUONVDigitStore& digitStore)
{
  /// Create MUON digit from ESDMuon digit and add it to the store
  /// return a pointer to the digit into the store (0x0 if the digit already exist)
  AliMUONVDigit* digit = digitStore.Add(esdPad.GetDetElemId(), esdPad.GetManuId(), esdPad.GetManuChannel(), esdPad.GetCathode(), AliMUONVDigitStore::kDeny);
  if (digit) ESDToMUON(esdPad, *digit);
  return digit;
}

