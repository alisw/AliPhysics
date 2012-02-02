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

//====================================================================================================================================================
//
//      Description of an ALICE muon forward track, combining the information of the Muon Spectrometer and the Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliLog.h"
#include "AliMUONTrack.h"
#include "AliMFTCluster.h"
#include "AliMUONVCluster.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "TClonesArray.h"
#include "TMatrixD.h"
#include "TParticle.h"
#include "AliMuonForwardTrack.h"

ClassImp(AliMuonForwardTrack)

//====================================================================================================================================================

AliMuonForwardTrack::AliMuonForwardTrack():
  AliMUONTrack(),
  fMUONTrack(0),
  fMCTrackRef(0),
  fMFTClusters(0)
{

  // default constructor
  for (Int_t iPlane=0; iPlane<fMaxNPlanesMFT; iPlane++) fPlaneExists[iPlane] = kFALSE;
  fMFTClusters = new TClonesArray("AliMFTCluster");

}

//====================================================================================================================================================

AliMuonForwardTrack::AliMuonForwardTrack(AliMUONTrack *MUONTrack):
  AliMUONTrack(),
  fMUONTrack(0),
  fMCTrackRef(0),
  fMFTClusters(0)
{

  SetMUONTrack(MUONTrack);
  for (Int_t iPlane=0; iPlane<fMaxNPlanesMFT; iPlane++) fPlaneExists[iPlane] = kFALSE;
  fMFTClusters = new TClonesArray("AliMFTCluster");

}

//====================================================================================================================================================

AliMuonForwardTrack::AliMuonForwardTrack(const AliMuonForwardTrack& track): 
  AliMUONTrack(track),
  fMUONTrack(track.fMUONTrack),
  fMCTrackRef(track.fMCTrackRef),
  fMFTClusters(track.fMFTClusters)
{

  // copy constructor
  for (Int_t iPlane=0; iPlane<fMaxNPlanesMFT; iPlane++) fPlaneExists[iPlane] = (track.fPlaneExists)[iPlane];
  
}

//====================================================================================================================================================

AliMuonForwardTrack& AliMuonForwardTrack::operator=(const AliMuonForwardTrack& track) {

  // Asignment operator

  // check assignement to self
  if (this == &track) return *this;

  // base class assignement
  AliMUONTrack::operator=(track);
  
  // clear memory
  Clear();
  
  fMUONTrack    = track.fMUONTrack;
  fMCTrackRef   = track.fMCTrackRef;
  fMFTClusters  = track.fMFTClusters;

  for (Int_t iPlane=0; iPlane<fMaxNPlanesMFT; iPlane++) fPlaneExists[iPlane] = (track.fPlaneExists)[iPlane];
  
  return *this;

}

//====================================================================================================================================================

void AliMuonForwardTrack::SetMUONTrack(AliMUONTrack *MUONTrack) {

  if (fMUONTrack) {
    AliInfo("fMUONTrack already exists, nothing will be done");
    return;
  }

  fMUONTrack = MUONTrack;
  
}

//====================================================================================================================================================

void AliMuonForwardTrack::SetMCTrackRef(TParticle *MCTrackRef) {

  if (fMCTrackRef) {
    AliInfo("fMCTrackRef already exists, nothing will be done");
    return;
  }

  fMCTrackRef = MCTrackRef;
  
}

//====================================================================================================================================================

void AliMuonForwardTrack::AddTrackParamAtMFTCluster(AliMUONTrackParam &trackParam, AliMFTCluster &mftCluster) {

  AliDebug(1, Form("fMFTClusters=%p has %d entries", fMFTClusters, fMFTClusters->GetEntries()));
  Int_t iMFTCluster = fMFTClusters->GetEntries();
  AliDebug(1, Form("mftCluster->GetX() = %f  mftCluster->GetY() = %f  mftCluster->GetErrX() = %f  cmftCluster->GetErrY() = %f", 
	   mftCluster.GetX(), mftCluster.GetY(), mftCluster.GetErrX(), mftCluster.GetErrY()));
  AliMUONVCluster *muonCluster = (AliMUONVCluster*) mftCluster.CreateMUONCluster();
  AliDebug(1, Form("Created MUON cluster %p", muonCluster));
  trackParam.SetUniqueID(iMFTCluster);    // we profit of this slot to store the reference to the corresponding MFTCluster
  AliDebug(1, Form("Now adding muonCluster %p and trackParam %p",muonCluster, &trackParam));
  AddTrackParamAtCluster(trackParam, *muonCluster, kTRUE);
  AliDebug(1, Form("GetTrackParamAtCluster() = %p  has %d entries",GetTrackParamAtCluster(), GetTrackParamAtCluster()->GetEntries()));
  // we pass the parameters this->GetTrackParamAtCluster()->First() to the Kalman Filter algorithm: they will be updated!!
  Double_t chi2Kalman = RunKalmanFilter(*(AliMUONTrackParam*)(GetTrackParamAtCluster()->First()));
  Double_t newGlobalChi2 = GetGlobalChi2() + chi2Kalman;
  mftCluster.SetLocalChi2(chi2Kalman);
  mftCluster.SetTrackChi2(newGlobalChi2);
  new ((*fMFTClusters)[iMFTCluster]) AliMFTCluster(mftCluster);
  AliDebug(1, Form("muonCluster->GetZ() = %f, trackParam->GetZ() = %f",muonCluster->GetZ(), trackParam.GetZ()));
  SetGlobalChi2(newGlobalChi2);
  ((AliMUONTrackParam*) GetTrackParamAtCluster()->First())->SetTrackChi2(newGlobalChi2);
  for (Int_t iPar=0; iPar<GetTrackParamAtCluster()->GetEntries(); iPar++) {
    AliDebug(1, Form("GetTrackParamAtCluster()->At(%d)->GetClusterPtr() = %p",
		     iPar, ((AliMUONTrackParam*)(GetTrackParamAtCluster()->At(iPar)))->GetClusterPtr()));
  }

}

//====================================================================================================================================================

AliMUONTrackParam* AliMuonForwardTrack::GetTrackParamAtMUONCluster(Int_t iMUONCluster) {

  if (iMUONCluster<0 || iMUONCluster>=GetNMUONClusters()) {
    AliError("Invalid MUON cluster index. NULL pointer will be returned");
    return NULL;
  }

  AliMUONTrackParam *trackParam = (AliMUONTrackParam*) fMUONTrack->GetTrackParamAtCluster()->At(iMUONCluster); 

  return trackParam;

}

//====================================================================================================================================================

AliMUONTrackParam* AliMuonForwardTrack::GetTrackParamAtMFTCluster(Int_t iMFTCluster) {

  if (iMFTCluster<0 || iMFTCluster>=GetNMFTClusters()) {
    AliError("Invalid MFT cluster index. NULL pointer will be returned");
    return NULL;
  }

  AliMUONTrackParam *trackParam = (AliMUONTrackParam*) GetTrackParamAtCluster()->At(iMFTCluster); 

  return trackParam;

}

//====================================================================================================================================================

AliMUONVCluster* AliMuonForwardTrack::GetMUONCluster(Int_t iMUONCluster) {

  if (iMUONCluster<0 || iMUONCluster>=GetNMUONClusters()) {
    AliError("Invalid MUON cluster index. NULL pointer will be returned");
    return NULL;
  }

  AliMUONTrackParam *trackParam = GetTrackParamAtMUONCluster(iMUONCluster);
  AliMUONVCluster *muonCluster = trackParam->GetClusterPtr();

  return muonCluster;

}

//====================================================================================================================================================

AliMFTCluster* AliMuonForwardTrack::GetMFTCluster(Int_t iMFTCluster) {

  if (iMFTCluster<0 || iMFTCluster>=GetNMFTClusters()) {
    AliError("Invalid MFT cluster index. NULL pointer will be returned");
    return NULL;
  }

  AliMUONTrackParam *trackParam = GetTrackParamAtMFTCluster(iMFTCluster);
  AliMFTCluster *mftCluster = (AliMFTCluster*) fMFTClusters->At(trackParam->GetUniqueID());

  return mftCluster;

}

//====================================================================================================================================================

Double_t AliMuonForwardTrack::RunKalmanFilter(AliMUONTrackParam &trackParamAtCluster) {

  AliDebug(1, Form("Running Kalman filter for parameters %p (z = %f) and cluster %p (z = %f)", 
		   &trackParamAtCluster, trackParamAtCluster.GetZ(), trackParamAtCluster.GetClusterPtr(), trackParamAtCluster.GetClusterPtr()->GetZ()));
  
  // Compute new track parameters and their covariances including new cluster using kalman filter
  // return the *additional* track chi2
  
  // Get actual track parameters (p)
  TMatrixD param(trackParamAtCluster.GetParameters());
  
  // Get new cluster parameters (m)
  AliMUONVCluster *cluster = trackParamAtCluster.GetClusterPtr();
  AliDebug(1, Form("cluster->GetX() = %f  cluster->GetY() = %f  cluster->GetErrX() = %f  cluster->GetErrY() = %f", 
		   cluster->GetX(), cluster->GetY(), cluster->GetErrX(), cluster->GetErrY()));
  TMatrixD clusterParam(5,1);
  clusterParam.Zero();
  clusterParam(0,0) = cluster->GetX();
  clusterParam(2,0) = cluster->GetY();

  // Compute the current parameter weight (W)
  TMatrixD paramWeight(trackParamAtCluster.GetCovariances());
  if (paramWeight.Determinant() != 0) {
    paramWeight.Invert();
  } else {
    Warning("RunKalmanFilter"," Determinant = 0");
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
  } 
  else {
    Warning("RunKalmanFilter"," Determinant = 0");
    return 1.e10;
  }
  
  // Save the new parameters covariance matrix
  trackParamAtCluster.SetCovariances(newParamCov);
  
  // Compute the new parameters (p' = ((W+U)^-1)U(m-p) + p)
  TMatrixD tmp(clusterParam,TMatrixD::kMinus,param);
  TMatrixD tmp2(clusterWeight,TMatrixD::kMult,tmp);    // U(m-p)
  TMatrixD newParam(newParamCov,TMatrixD::kMult,tmp2); // ((W+U)^-1)U(m-p)
  newParam += param;                                   // ((W+U)^-1)U(m-p) + p

  // Save the new parameters
  trackParamAtCluster.SetParameters(newParam);
  
  // Compute the additional chi2 (= ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m))
  tmp = newParam; // p'
  tmp -= param;   // (p'-p)
  TMatrixD tmp3(paramWeight,TMatrixD::kMult,tmp);           // W(p'-p)
  TMatrixD addChi2Track(tmp,TMatrixD::kTransposeMult,tmp3); // ((p'-p)^-1)W(p'-p)
  tmp = newParam;      // p'
  tmp -= clusterParam; // (p'-m)
  TMatrixD tmp4(clusterWeight,TMatrixD::kMult,tmp);            // U(p'-m)
  addChi2Track += TMatrixD(tmp,TMatrixD::kTransposeMult,tmp4); // ((p'-p)^-1)W(p'-p) + ((p'-m)^-1)U(p'-m)
  
  AliDebug(1,Form("Adding Kalman chi2 = %f",addChi2Track(0,0)));

  return addChi2Track(0,0);
  
}

//====================================================================================================================================================

Double_t AliMuonForwardTrack::GetWeightedOffset(Double_t x, Double_t y, Double_t z) {

  AliMUONTrackParam *param = GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZCov(param, z);

  TMatrixD cov(5,5);
  cov = param->GetCovariances();

  TMatrixD covCoordinates(2,2);
  covCoordinates(0,0) = cov(0,0);
  covCoordinates(0,1) = cov(0,2);
  covCoordinates(1,0) = cov(2,0);
  covCoordinates(1,1) = cov(2,2);
  
  TMatrixD covCoordinatesInverse = covCoordinates.Invert();

  Double_t dX = param->GetNonBendingCoor() - x;
  Double_t dY = param->GetBendingCoor()    - y;
  
  Double_t weightedOffset = TMath::Sqrt(0.5*(dX*dX*covCoordinatesInverse(0,0) + 
					     dY*dY*covCoordinatesInverse(1,1) + 
					     2.*dX*dY*covCoordinatesInverse(0,1)));

  return weightedOffset;

}

//====================================================================================================================================================

Double_t AliMuonForwardTrack::GetOffsetX(Double_t x, Double_t z) {

  AliMUONTrackParam *param = GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZCov(param, z);
  Double_t dX = param->GetNonBendingCoor() - x;
  return dX;

}

//====================================================================================================================================================

Double_t AliMuonForwardTrack::GetOffsetY(Double_t y, Double_t z) {

  AliMUONTrackParam *param = GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZCov(param, z);
  Double_t dY = param->GetBendingCoor() - y;
  return dY;

}

//====================================================================================================================================================

Double_t AliMuonForwardTrack::GetOffset(Double_t x, Double_t y, Double_t z) {

  AliMUONTrackParam *param = GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZCov(param, z);
  Double_t dX = param->GetNonBendingCoor() - x;
  Double_t dY = param->GetBendingCoor()    - y;
  return TMath::Sqrt(dX*dX + dY*dY);

}

//====================================================================================================================================================

