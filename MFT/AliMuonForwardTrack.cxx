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
#include "AliMFTConstants.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "AliMUONConstants.h"

ClassImp(AliMuonForwardTrack)

//====================================================================================================================================================

AliMuonForwardTrack::AliMuonForwardTrack():
  AliMUONTrack(),
  fMUONTrack(0),
  fMCTrackRef(0),
  fMFTClusters(0),
  fNWrongClustersMC(-1),
  fTrackMCId(-1),
  fKinem(0,0,0,0),
  fParamCovMatrix(5,5)
{

  // default constructor
  
  for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) fPlaneExists[iPlane] = kFALSE;
  for (Int_t iParent=0; iParent<fgkNParentsMax; iParent++) {
    fParentMCLabel[iParent] = -1;
    fParentPDGCode[iParent] =  0;
  }
  fMFTClusters = new TClonesArray("AliMFTCluster");
  fMFTClusters -> SetOwner(kTRUE);

}

//====================================================================================================================================================

AliMuonForwardTrack::AliMuonForwardTrack(AliMUONTrack *MUONTrack):
  AliMUONTrack(),
  fMUONTrack(0),
  fMCTrackRef(0),
  fMFTClusters(0),
  fNWrongClustersMC(-1),
  fTrackMCId(-1),
  fKinem(0,0,0,0),
  fParamCovMatrix(5,5)
{

  SetMUONTrack(MUONTrack);
  for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) fPlaneExists[iPlane] = kFALSE;
  for (Int_t iParent=0; iParent<fgkNParentsMax; iParent++) {
    fParentMCLabel[iParent] = -1;
    fParentPDGCode[iParent] =  0;
  }
  fMFTClusters = new TClonesArray("AliMFTCluster");
  fMFTClusters -> SetOwner(kTRUE);

}

//====================================================================================================================================================

AliMuonForwardTrack::AliMuonForwardTrack(const AliMuonForwardTrack& track): 
  AliMUONTrack(track),
  fMUONTrack(0x0),
  fMCTrackRef(0x0),
  fMFTClusters(0x0),
  fNWrongClustersMC(track.fNWrongClustersMC),
  fTrackMCId(track.fTrackMCId),
  fKinem(track.fKinem),
  fParamCovMatrix(track.fParamCovMatrix)
{

  // copy constructor
  fMUONTrack        = new AliMUONTrack(*(track.fMUONTrack));
  if (track.fMCTrackRef) fMCTrackRef = new TParticle(*(track.fMCTrackRef));
  fMFTClusters      = new TClonesArray(*(track.fMFTClusters));
  fMFTClusters->SetOwner(kTRUE);
  for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) fPlaneExists[iPlane] = (track.fPlaneExists)[iPlane];
  for (Int_t iParent=0; iParent<fgkNParentsMax; iParent++) {
    fParentMCLabel[iParent] = (track.fParentMCLabel)[iParent];
    fParentPDGCode[iParent] = (track.fParentPDGCode)[iParent];
  }
  
}

//====================================================================================================================================================

AliMuonForwardTrack& AliMuonForwardTrack::operator=(const AliMuonForwardTrack& track) {

  // Asignment operator

  // check assignement to self
  if (this == &track) return *this;

  // base class assignement
  AliMUONTrack::operator=(track);
  
  // clear memory
  Clear("");
  
  fMUONTrack        = new AliMUONTrack(*(track.fMUONTrack));
  if (track.fMCTrackRef) fMCTrackRef = new TParticle(*(track.fMCTrackRef));
  fMFTClusters      = new TClonesArray(*(track.fMFTClusters));
  fMFTClusters->SetOwner(kTRUE);
  fNWrongClustersMC = track.fNWrongClustersMC;
  fTrackMCId        = track.fTrackMCId;
  fKinem            = track.fKinem;
  fParamCovMatrix   = track.fParamCovMatrix;

  for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) fPlaneExists[iPlane] = (track.fPlaneExists)[iPlane];
  for (Int_t iParent=0; iParent<fgkNParentsMax; iParent++) {
    fParentMCLabel[iParent] = (track.fParentMCLabel)[iParent];
    fParentPDGCode[iParent] = (track.fParentPDGCode)[iParent];
  }
  
  return *this;

}

//====================================================================================================================================================

void AliMuonForwardTrack::Clear(const Option_t* /*opt*/) {

  // Clear arrays
  fMFTClusters -> Delete(); 
  delete fMFTClusters; fMFTClusters = 0x0;
  delete fMUONTrack;   fMUONTrack   = 0x0;
  delete fMCTrackRef;  fMCTrackRef  = 0x0;
  
}

//====================================================================================================================================================

AliMuonForwardTrack::~AliMuonForwardTrack() {

  delete fMUONTrack;
  delete fMCTrackRef;
  fMFTClusters -> Delete();
  delete fMFTClusters;
  
}

//====================================================================================================================================================

void AliMuonForwardTrack::SetMUONTrack(AliMUONTrack *MUONTrack) {

  if (fMUONTrack) {
    AliInfo("fMUONTrack already exists, nothing will be done");
    return;
  }

  fMUONTrack = MUONTrack;
  SetGlobalChi2(fMUONTrack->GetGlobalChi2());
  
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

  AliDebug(1, Form("Before adding: this->fMFTClusters=%p has %d entries", this->fMFTClusters, this->fMFTClusters->GetEntries()));
  Int_t iMFTCluster = this->fMFTClusters->GetEntries();
  AliDebug(1, Form("mftCluster->GetX() = %f  mftCluster->GetY() = %f  mftCluster->GetErrX() = %f  cmftCluster->GetErrY() = %f", 
	   mftCluster.GetX(), mftCluster.GetY(), mftCluster.GetErrX(), mftCluster.GetErrY()));
  AliMUONVCluster *muonCluster = (AliMUONVCluster*) mftCluster.CreateMUONCluster();
  AliDebug(1, Form("Created MUON cluster %p", muonCluster));
  trackParam.SetUniqueID(iMFTCluster);    // we profit of this slot to store the reference to the corresponding MFTCluster
  AliDebug(1, Form("Now adding muonCluster %p and trackParam %p",muonCluster, &trackParam));
  AddTrackParamAtCluster(trackParam, *muonCluster, kTRUE);
  // we pass the parameters this->GetTrackParamAtCluster()->First() to the Kalman Filter algorithm: they will be updated!!
  Double_t chi2Kalman = RunKalmanFilter(*(AliMUONTrackParam*)(GetTrackParamAtCluster()->First()));
  AliDebug(1, Form("Adding Kalman chi2 = %f to global chi2 = %f", chi2Kalman, GetGlobalChi2()));
  Double_t newGlobalChi2 = GetGlobalChi2() + chi2Kalman;
  mftCluster.SetLocalChi2(chi2Kalman);
  mftCluster.SetTrackChi2(newGlobalChi2);
  new ((*(this->fMFTClusters))[iMFTCluster]) AliMFTCluster(mftCluster);
  AliDebug(1, Form("GetTrackParamAtCluster() = %p  has %d entries while this->fMFTClusters=%p has %d entries",
		   GetTrackParamAtCluster(), GetTrackParamAtCluster()->GetEntries(), this->fMFTClusters, this->fMFTClusters->GetEntries()));
  AliDebug(1, Form("muonCluster->GetZ() = %f, trackParam->GetZ() = %f",muonCluster->GetZ(), trackParam.GetZ()));
  SetGlobalChi2(newGlobalChi2);
  AliDebug(1, Form("New global chi2 = %f", GetGlobalChi2()));
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
    AliError(Form("Invalid MFT cluster index (%d). GetNMFTClusters()=%d. NULL pointer will be returned", iMFTCluster, GetNMFTClusters()));
    return NULL;
  }

  AliMUONTrackParam *trackParam = GetTrackParamAtMFTCluster(iMFTCluster);
  AliMFTCluster *mftCluster = (AliMFTCluster*) this->fMFTClusters->At(trackParam->GetUniqueID());

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

Double_t AliMuonForwardTrack::GetDCA(Double_t x, Double_t y, Double_t z) {

  // Distance of Closest Approach, according to the standard MUON terminology. Actually, the offset of the track w.r.t. the primary vertex,
  // where the extrapolation of the track DOES NOT include the MFT information

  AliMUONTrackParam *param = GetTrackParamAtMUONCluster(0);
  AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(param, z);
  Double_t dX = param->GetNonBendingCoor() - x;
  Double_t dY = param->GetBendingCoor()    - y;
  return TMath::Sqrt(dX*dX + dY*dY);

}

//====================================================================================================================================================

Double_t AliMuonForwardTrack::GetMomentumSpectrometer(Double_t z) {

  // Momentum of the track at the primary vertex plane, where the extrapolation of the track DOES NOT include the MFT information

  AliMUONTrackParam *param = GetTrackParamAtMUONCluster(0);
  AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(param, z);
  return param->P();

}

//====================================================================================================================================================

Double_t AliMuonForwardTrack::GetThetaAbs() {

  // it is the angle defined by the imaginary line goingo from the vertex to the exit point of the track at the end of the hadron absorber
  
  Double_t z = AliMUONConstants::AbsZEnd();
  AliMUONTrackParam *param = GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZ(param, z);
  Double_t x = param->GetNonBendingCoor();
  Double_t y = param->GetBendingCoor();

  return 180. +TMath::ATan(TMath::Sqrt(x*x + y*y)/z)*TMath::RadToDeg();

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrack::IsFromResonance() {

  Bool_t result = kFALSE;

  if ( GetParentPDGCode(0) ==    113 ||
       GetParentPDGCode(0) ==    221 ||
       GetParentPDGCode(0) ==    223 ||
       GetParentPDGCode(0) ==    331 ||
       GetParentPDGCode(0) ==    333 ||
       GetParentPDGCode(0) ==    443 ||
       GetParentPDGCode(0) == 100443 ||
       GetParentPDGCode(0) ==    553 ||
       GetParentPDGCode(0) == 100553 ) result = kTRUE;
  
  if (result) AliDebug(1, Form("Muon comes from a resonance %d", GetParentPDGCode(0)));
  
  return result; 
  
}

//====================================================================================================================================================

Bool_t AliMuonForwardTrack::IsDirectCharm() {

  Bool_t result = kFALSE;

  if (IsPDGCharm(GetParentPDGCode(0)) && !IsPDGBeauty(GetParentPDGCode(1))) result = kTRUE;
  
  if (result) AliDebug(1, Form("Muon comes from a charmed hadron %d", GetParentPDGCode(0)));
  
  return result; 
  
}

//====================================================================================================================================================

Bool_t AliMuonForwardTrack::IsDirectBeauty() {

  Bool_t result = kFALSE;

  if (IsPDGBeauty(GetParentPDGCode(0))) result = kTRUE;
  
  if (result) AliDebug(1, Form("Muon comes from a beauty hadron %d", GetParentPDGCode(0)));
  
  return result; 
  
}

//====================================================================================================================================================

Bool_t AliMuonForwardTrack::IsChainBeauty() {

  Bool_t result = kFALSE;

  if (IsPDGCharm(GetParentPDGCode(0)) && IsPDGBeauty(GetParentPDGCode(1))) result = kTRUE;
  
  if (result) AliDebug(1, Form("Muon comes from a charmed hadron %d which comes from a beauty hadron %d", GetParentPDGCode(0), GetParentPDGCode(1)));
  
  return result; 
  
}

//====================================================================================================================================================

Bool_t AliMuonForwardTrack::IsPDGCharm(Int_t pdg) {

  Bool_t result = kFALSE;

  if ( TMath::Abs(pdg) ==   411 ||
       TMath::Abs(pdg) ==   421 ||
       TMath::Abs(pdg) == 10411 ||
       TMath::Abs(pdg) == 10421 ||
       TMath::Abs(pdg) ==   413 ||
       TMath::Abs(pdg) ==   423 ||
       TMath::Abs(pdg) == 10413 ||
       TMath::Abs(pdg) == 10423 ||
       TMath::Abs(pdg) == 20413 ||
       TMath::Abs(pdg) == 20423 ||
       TMath::Abs(pdg) ==   415 ||
       TMath::Abs(pdg) ==   425 ||
       TMath::Abs(pdg) ==   431 ||
       TMath::Abs(pdg) == 10431 ||
       TMath::Abs(pdg) ==   433 ||
       TMath::Abs(pdg) == 10433 ||
       TMath::Abs(pdg) == 20433 ||
       TMath::Abs(pdg) ==   435 ) result = kTRUE;
  
  return result; 
  
}

//====================================================================================================================================================

Bool_t AliMuonForwardTrack::IsPDGBeauty(Int_t pdg) {

  Bool_t result = kFALSE;

  if ( TMath::Abs(pdg) ==   511 ||
       TMath::Abs(pdg) ==   521 ||
       TMath::Abs(pdg) == 10511 ||
       TMath::Abs(pdg) == 10521 ||
       TMath::Abs(pdg) ==   513 ||
       TMath::Abs(pdg) ==   523 ||
       TMath::Abs(pdg) == 10513 ||
       TMath::Abs(pdg) == 10523 ||
       TMath::Abs(pdg) == 20513 ||
       TMath::Abs(pdg) == 20523 ||
       TMath::Abs(pdg) ==   515 ||
       TMath::Abs(pdg) ==   525 ||
       TMath::Abs(pdg) ==   531 ||
       TMath::Abs(pdg) == 10531 ||
       TMath::Abs(pdg) ==   533 ||
       TMath::Abs(pdg) == 10533 ||
       TMath::Abs(pdg) == 20533 ||
       TMath::Abs(pdg) ==   535 ||
       TMath::Abs(pdg) ==   541 ||
       TMath::Abs(pdg) == 10541 ||
       TMath::Abs(pdg) ==   543 ||
       TMath::Abs(pdg) == 10543 ||
       TMath::Abs(pdg) == 20543 ||
       TMath::Abs(pdg) ==   545 ) result = kTRUE;
  
  return result; 
  
}

//====================================================================================================================================================

Bool_t AliMuonForwardTrack::IsFromBackground() {

  Bool_t result = kFALSE;

  if (!IsFromResonance() && !IsDirectCharm() && !IsDirectBeauty() && !IsChainBeauty()) result = kTRUE;

  if (result) AliDebug(1, Form("Muon comes from a background source %d", GetParentPDGCode(0)));

  return result;

}

//====================================================================================================================================================

void AliMuonForwardTrack::EvalKinem(Double_t z) {

  AliMUONTrackParam *param = GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZCov(param, z);

  Double_t mMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t energy = TMath::Sqrt(param->P()*param->P() + mMu*mMu);
  fKinem.SetPxPyPzE(param->Px(), param->Py(), param->Pz(), energy);

  TMatrixD cov(5,5);
  fParamCovMatrix = param->GetCovariances();

}

//====================================================================================================================================================

