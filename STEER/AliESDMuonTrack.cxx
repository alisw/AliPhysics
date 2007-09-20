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

///////////////////////////////////////////////////////////////////////////////
//
//  Class to describe the MUON tracks
//  in the Event Summary Data class
//  This is where the results of reconstruction
//  are stored for the muons
//  Author: G.Martinez
//
///////////////////////////////////////////////////////////////////////////////

#include "AliESDMuonTrack.h"

#include <TLorentzVector.h>
#include <TMath.h>

ClassImp(AliESDMuonTrack)

//_____________________________________________________________________________
AliESDMuonTrack::AliESDMuonTrack ():
  AliVParticle(),
  fInverseBendingMomentum(0),
  fThetaX(0),
  fThetaY(0),
  fZ(0),
  fBendingCoor(0),
  fNonBendingCoor(0),
  fInverseBendingMomentumUncorrected(0),
  fThetaXUncorrected(0),
  fThetaYUncorrected(0),
  fZUncorrected(0),
  fBendingCoorUncorrected(0),
  fNonBendingCoorUncorrected(0),
  fChi2(0),
  fNHit(0),
  fLocalTrigger(234),
  fChi2MatchTrigger(0),
  fHitsPatternInTrigCh(0),
  fMuonClusterMap(0)
{
  //
  // Default constructor
  //
  for (Int_t i = 0; i < 15; i++) fCovariances[i] = 0;
}


//_____________________________________________________________________________
AliESDMuonTrack::AliESDMuonTrack (const AliESDMuonTrack& MUONTrack):
  AliVParticle(MUONTrack),
  fInverseBendingMomentum(MUONTrack.fInverseBendingMomentum),
  fThetaX(MUONTrack.fThetaX),
  fThetaY(MUONTrack.fThetaY),
  fZ(MUONTrack.fZ),
  fBendingCoor(MUONTrack.fBendingCoor),
  fNonBendingCoor(MUONTrack.fNonBendingCoor),
  fInverseBendingMomentumUncorrected(MUONTrack.fInverseBendingMomentumUncorrected),
  fThetaXUncorrected(MUONTrack.fThetaXUncorrected),
  fThetaYUncorrected(MUONTrack.fThetaYUncorrected),
  fZUncorrected(MUONTrack.fZUncorrected),
  fBendingCoorUncorrected(MUONTrack.fBendingCoorUncorrected),
  fNonBendingCoorUncorrected(MUONTrack.fNonBendingCoorUncorrected),
  fChi2(MUONTrack.fChi2),
  fNHit(MUONTrack.fNHit),
  fLocalTrigger(MUONTrack.fLocalTrigger),
  fChi2MatchTrigger(MUONTrack.fChi2MatchTrigger),
  fHitsPatternInTrigCh(MUONTrack.fHitsPatternInTrigCh),
  fMuonClusterMap(MUONTrack.fMuonClusterMap)
{
  //
  // Copy constructor
  // Deep copy implemented
  //
  for (Int_t i = 0; i < 15; i++) fCovariances[i] = MUONTrack.fCovariances[i];
}

//_____________________________________________________________________________
AliESDMuonTrack& AliESDMuonTrack::operator=(const AliESDMuonTrack& MUONTrack)
{
  // 
  // Equal operator for a deep copy
  //
  if (this == &MUONTrack)
    return *this;

  AliVParticle::operator=(MUONTrack); // don't forget to invoke the base class' assignment operator
  
  fInverseBendingMomentum = MUONTrack.fInverseBendingMomentum; 
  fThetaX                 = MUONTrack.fThetaX;           
  fThetaY                 = MUONTrack.fThetaY;           
  fZ                      = MUONTrack.fZ;                
  fBendingCoor            = MUONTrack.fBendingCoor;      
  fNonBendingCoor         = MUONTrack.fNonBendingCoor;   
  
  fInverseBendingMomentumUncorrected = MUONTrack.fInverseBendingMomentumUncorrected; 
  fThetaXUncorrected                 = MUONTrack.fThetaXUncorrected;           
  fThetaYUncorrected                 = MUONTrack.fThetaYUncorrected;           
  fZUncorrected                      = MUONTrack.fZUncorrected;                
  fBendingCoorUncorrected            = MUONTrack.fBendingCoorUncorrected;      
  fNonBendingCoorUncorrected         = MUONTrack.fNonBendingCoorUncorrected;   
  
  for (Int_t i = 0; i < 15; i++) fCovariances[i] = MUONTrack.fCovariances[i];
  
  fChi2                   = MUONTrack.fChi2;             
  fNHit                   = MUONTrack.fNHit; 

  fLocalTrigger           = MUONTrack.fLocalTrigger;  
  fChi2MatchTrigger       = MUONTrack.fChi2MatchTrigger; 

  fHitsPatternInTrigCh    = MUONTrack.fHitsPatternInTrigCh;
 
  fMuonClusterMap	  = MUONTrack.fMuonClusterMap;
  
  return *this;
}

//_____________________________________________________________________________
void AliESDMuonTrack::GetCovariances(TMatrixD& cov) const
{
  // return covariance matrix of uncorrected parameters
  cov.ResizeTo(5,5);
  for (Int_t i = 0; i < 5; i++)
    for (Int_t j = 0; j <= i; j++)
      cov(i,j) = cov (j,i) = fCovariances[i*(i+1)/2 + j];
}

//_____________________________________________________________________________
void AliESDMuonTrack::SetCovariances(const TMatrixD& cov)
{
  // set reduced covariance matrix of uncorrected parameters
  for (Int_t i = 0; i < 5; i++)
    for (Int_t j = 0; j <= i; j++)
      fCovariances[i*(i+1)/2 + j] = cov(i,j);

}

//_____________________________________________________________________________
void AliESDMuonTrack::GetCovarianceXYZPxPyPz(Double_t cov[21]) const
{
  // return reduced covariance matrix of uncorrected parameters in (X,Y,Z,Px,Py,Pz) coordinate system
  // 
  // Cov(x,x) ... :   cov[0]
  // Cov(y,x) ... :   cov[1]  cov[2]
  // Cov(z,x) ... :   cov[3]  cov[4]  cov[5]
  // Cov(px,x)... :   cov[6]  cov[7]  cov[8]  cov[9]
  // Cov(py,x)... :   cov[10] cov[11] cov[12] cov[13] cov[14]
  // Cov(pz,x)... :   cov[15] cov[16] cov[17] cov[18] cov[19] cov[20]
  //
  // Get ESD covariance matrix into a TMatrixD
  TMatrixD covESD(5,5);
  GetCovariances(covESD);

  // compute Jacobian to change the coordinate system
  // from (X,thetaX,Y,thetaY,c/pYZ) to (X,Y,Z,pX,pY,pZ)
  Double_t tanThetaX = TMath::Tan(fThetaXUncorrected);
  Double_t tanThetaY = TMath::Tan(fThetaYUncorrected);
  Double_t cosThetaX2 = TMath::Cos(fThetaXUncorrected) * TMath::Cos(fThetaXUncorrected);
  Double_t cosThetaY2 = TMath::Cos(fThetaYUncorrected) * TMath::Cos(fThetaYUncorrected);
  Double_t pZ = PzUncorrected();
  Double_t dpZdthetaY = - fInverseBendingMomentumUncorrected * fInverseBendingMomentumUncorrected *
			  pZ * pZ * pZ * tanThetaY / cosThetaY2;
  Double_t dpZdinvpYZ = - pZ / fInverseBendingMomentumUncorrected;
  TMatrixD jacob(6,5);
  jacob.Zero();
  jacob(0,0) = 1.;
  jacob(1,2) = 1.;
  jacob(3,1) = pZ / cosThetaX2;
  jacob(3,3) = dpZdthetaY * tanThetaX;
  jacob(3,4) = dpZdinvpYZ * tanThetaX;
  jacob(4,3) = dpZdthetaY * tanThetaY + pZ / cosThetaY2;
  jacob(4,4) = dpZdinvpYZ * tanThetaY;
  jacob(5,3) = dpZdthetaY;
  jacob(5,4) = dpZdinvpYZ;
  
  // compute covariance matrix in AOD coordinate system
  TMatrixD tmp(covESD,TMatrixD::kMultTranspose,jacob);
  TMatrixD covAOD(jacob,TMatrixD::kMult,tmp);
  
  // Get AOD covariance matrix into co[21]
  for (Int_t i = 0; i < 6; i++)
    for (Int_t j = 0; j <= i; j++)
      cov[i*(i+1)/2 + j] = covAOD(i,j);
  
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::Px() const
{
  // return p_x from track parameters
  Double_t nonBendingSlope = TMath::Tan(fThetaX);
  Double_t bendingSlope    = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return pZ * nonBendingSlope;
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::Py() const
{
  // return p_y from track parameters
  Double_t bendingSlope = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return pZ * bendingSlope;
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::Pz() const
{
  // return p_z from track parameters
  Double_t bendingSlope = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  return -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::P() const
{
  // return p from track parameters
  Double_t nonBendingSlope = TMath::Tan(fThetaX);
  Double_t bendingSlope    = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return -pZ * TMath::Sqrt(1.0 + bendingSlope*bendingSlope + nonBendingSlope*nonBendingSlope);
}

//_____________________________________________________________________________
void AliESDMuonTrack::LorentzP(TLorentzVector& vP) const
{
  // return Lorentz momentum vector from track parameters
  Double_t muonMass = M();
  Double_t nonBendingSlope = TMath::Tan(fThetaX);
  Double_t bendingSlope    = TMath::Tan(fThetaY);
  Double_t pYZ = (fInverseBendingMomentum != 0.) ? TMath::Abs(1. / fInverseBendingMomentum) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  Double_t pX  = pZ * nonBendingSlope;
  Double_t pY  = pZ * bendingSlope;
  Double_t e   = TMath::Sqrt(muonMass*muonMass + pX*pX + pY*pY + pZ*pZ);
  vP.SetPxPyPzE(pX, pY, pZ, e);
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::PxUncorrected() const
{
  // return p_x from track parameters
  Double_t nonBendingSlope = TMath::Tan(fThetaXUncorrected);
  Double_t bendingSlope    = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return pZ * nonBendingSlope;
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::PyUncorrected() const
{
  // return p_y from track parameters
  Double_t bendingSlope = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return pZ * bendingSlope;
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::PzUncorrected() const
{
  // return p_z from track parameters
  Double_t bendingSlope = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  return -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
}

//_____________________________________________________________________________
Double_t AliESDMuonTrack::PUncorrected() const
{
  // return p from track parameters
  Double_t nonBendingSlope = TMath::Tan(fThetaXUncorrected);
  Double_t bendingSlope    = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  return -pZ * TMath::Sqrt(1.0 + bendingSlope*bendingSlope + nonBendingSlope*nonBendingSlope);
}

//_____________________________________________________________________________
void AliESDMuonTrack::LorentzPUncorrected(TLorentzVector& vP) const
{
  // return Lorentz momentum vector from track parameters
  Double_t muonMass = M();
  Double_t nonBendingSlope = TMath::Tan(fThetaXUncorrected);
  Double_t bendingSlope    = TMath::Tan(fThetaYUncorrected);
  Double_t pYZ = (fInverseBendingMomentumUncorrected != 0.) ? TMath::Abs(1. / fInverseBendingMomentumUncorrected) : 0.;
  Double_t pZ  = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);  // spectro. (z<0)
  Double_t pX  = pZ * nonBendingSlope;
  Double_t pY  = pZ * bendingSlope;
  Double_t e   = TMath::Sqrt(muonMass*muonMass + pX*pX + pY*pY + pZ*pZ);
  vP.SetPxPyPzE(pX, pY, pZ, e);
}

//_____________________________________________________________________________
Int_t AliESDMuonTrack::GetMatchTrigger() const
{
  //  backward compatibility after replacing fMatchTrigger by fLocalTrigger
  //  0 track does not match trigger
  //  1 track match but does not pass pt cut
  //  2 track match Low pt cut
  //  3 track match High pt cut

  if (LoCircuit() == -1) {
    return 0;
  } else if (LoLpt() == 0 && LoHpt() == 0) {
    return 1;
  } else if (LoLpt() >  0 && LoHpt() == 0) {
    return 2;
  } else {
    return 3;
  }

}

//_____________________________________________________________________________
void AliESDMuonTrack::AddInMuonClusterMap(Int_t chamber)
{
  // Update the muon cluster map by adding this chamber(0..)
  
  static const UInt_t kMask[10] = {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200};
  
  fMuonClusterMap |= kMask[chamber];
  
}

//_____________________________________________________________________________
Bool_t AliESDMuonTrack::IsInMuonClusterMap(Int_t chamber) const
{
  // return kTRUE if this chamber(0..) is in the muon cluster map
  
  static const UInt_t kMask[10] = {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200};
  
  return ((fMuonClusterMap | kMask[chamber]) == fMuonClusterMap) ? kTRUE : kFALSE;
  
}

