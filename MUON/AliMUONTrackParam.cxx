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
// Class AliMUONTrackParam
//-------------------------
// Track parameters in ALICE dimuon spectrometer
//-----------------------------------------------------------------------------

#include "AliMUONTrackParam.h"
#include "AliMUONVCluster.h"

#include "AliESDMuonTrack.h"
#include "AliLog.h"

#include <TMath.h>

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrackParam) // Class implementation in ROOT context
/// \endcond

  //_________________________________________________________________________
AliMUONTrackParam::AliMUONTrackParam()
  : TObject(),
    fZ(0.),
    fParameters(5,1),
    fCovariances(0x0),
    fPropagator(0x0),
    fExtrapParameters(0x0),
    fExtrapCovariances(0x0),
    fSmoothParameters(0x0),
    fSmoothCovariances(0x0),
    fClusterPtr(0x0),
    fOwnCluster(kFALSE),
    fRemovable(kFALSE),
    fAloneInChamber(kTRUE),
    fTrackChi2(0.),
    fLocalChi2(0.)
{
  /// Constructor
}

  //_________________________________________________________________________
AliMUONTrackParam::AliMUONTrackParam(const AliMUONTrackParam& theMUONTrackParam)
  : TObject(theMUONTrackParam),
    fZ(theMUONTrackParam.fZ),
    fParameters(theMUONTrackParam.fParameters),
    fCovariances(0x0),
    fPropagator(0x0),
    fExtrapParameters(0x0),
    fExtrapCovariances(0x0),
    fSmoothParameters(0x0),
    fSmoothCovariances(0x0),
    fClusterPtr(0x0),
    fOwnCluster(theMUONTrackParam.fOwnCluster),
    fRemovable(theMUONTrackParam.fRemovable),
    fAloneInChamber(theMUONTrackParam.fAloneInChamber),
    fTrackChi2(theMUONTrackParam.fTrackChi2),
    fLocalChi2(theMUONTrackParam.fLocalChi2)
{
  /// Copy constructor
  if (theMUONTrackParam.fCovariances) fCovariances = new TMatrixD(*(theMUONTrackParam.fCovariances));
  if (theMUONTrackParam.fPropagator) fPropagator = new TMatrixD(*(theMUONTrackParam.fPropagator));
  if (theMUONTrackParam.fExtrapParameters) fExtrapParameters = new TMatrixD(*(theMUONTrackParam.fExtrapParameters));
  if (theMUONTrackParam.fExtrapCovariances) fExtrapCovariances = new TMatrixD(*(theMUONTrackParam.fExtrapCovariances));
  if (theMUONTrackParam.fSmoothParameters) fSmoothParameters = new TMatrixD(*(theMUONTrackParam.fSmoothParameters));
  if (theMUONTrackParam.fSmoothCovariances) fSmoothCovariances = new TMatrixD(*(theMUONTrackParam.fSmoothCovariances));
  
  if(fOwnCluster) fClusterPtr = theMUONTrackParam.fClusterPtr->CreateCopy();
  else fClusterPtr = theMUONTrackParam.fClusterPtr;
}

  //_________________________________________________________________________
AliMUONTrackParam& AliMUONTrackParam::operator=(const AliMUONTrackParam& theMUONTrackParam)
{
  /// Asignment operator
  if (this == &theMUONTrackParam)
    return *this;

  // base class assignement
  TObject::operator=(theMUONTrackParam);

  fZ = theMUONTrackParam.fZ; 
  
  fParameters = theMUONTrackParam.fParameters;
  
  if (theMUONTrackParam.fCovariances) {
    if (fCovariances) *fCovariances = *(theMUONTrackParam.fCovariances);
    else fCovariances = new TMatrixD(*(theMUONTrackParam.fCovariances));
  } else {
    delete fCovariances;
    fCovariances = 0x0;
  }
  
  if (theMUONTrackParam.fPropagator) {
    if (fPropagator) *fPropagator = *(theMUONTrackParam.fPropagator);
    else fPropagator = new TMatrixD(*(theMUONTrackParam.fPropagator));
  } else {
    delete fPropagator;
    fPropagator = 0x0;
  }
  
  if (theMUONTrackParam.fExtrapParameters) {
    if (fExtrapParameters) *fExtrapParameters = *(theMUONTrackParam.fExtrapParameters);
    else fExtrapParameters = new TMatrixD(*(theMUONTrackParam.fExtrapParameters));
  } else {
    delete fExtrapParameters;
    fExtrapParameters = 0x0;
  }
  
  if (theMUONTrackParam.fExtrapCovariances) {
    if (fExtrapCovariances) *fExtrapCovariances = *(theMUONTrackParam.fExtrapCovariances);
    else fExtrapCovariances = new TMatrixD(*(theMUONTrackParam.fExtrapCovariances));
  } else {
    delete fExtrapCovariances;
    fExtrapCovariances = 0x0;
  }
  
  if (theMUONTrackParam.fSmoothParameters) {
    if (fSmoothParameters) *fSmoothParameters = *(theMUONTrackParam.fSmoothParameters);
    else fSmoothParameters = new TMatrixD(*(theMUONTrackParam.fSmoothParameters));
  } else {
    delete fSmoothParameters;
    fSmoothParameters = 0x0;
  }
  
  if (theMUONTrackParam.fSmoothCovariances) {
    if (fSmoothCovariances) *fSmoothCovariances = *(theMUONTrackParam.fSmoothCovariances);
    else fSmoothCovariances = new TMatrixD(*(theMUONTrackParam.fSmoothCovariances));
  } else {
    delete fSmoothCovariances;
    fSmoothCovariances = 0x0;
  }
  
  fOwnCluster = theMUONTrackParam.fOwnCluster;
  if(fOwnCluster) fClusterPtr = theMUONTrackParam.fClusterPtr->CreateCopy();
  else fClusterPtr = theMUONTrackParam.fClusterPtr;
  
  fRemovable = theMUONTrackParam.fRemovable;
  
  fAloneInChamber = theMUONTrackParam.fAloneInChamber;
  
  fTrackChi2 = theMUONTrackParam.fTrackChi2;
  fLocalChi2 = theMUONTrackParam.fLocalChi2;
  
  return *this;
}

  //__________________________________________________________________________
AliMUONTrackParam::~AliMUONTrackParam()
{
/// Destructor
  DeleteCovariances();
  delete fPropagator;
  delete fExtrapParameters;
  delete fExtrapCovariances;
  delete fSmoothParameters;
  delete fSmoothCovariances;
  if(fOwnCluster) delete fClusterPtr;
}

  //__________________________________________________________________________
void
AliMUONTrackParam::Clear(Option_t* /*opt*/)
{
  /// clear memory
  DeleteCovariances();
  delete fPropagator; fPropagator = 0x0;
  delete fExtrapParameters; fExtrapParameters = 0x0;
  delete fExtrapCovariances; fExtrapCovariances = 0x0;
  delete fSmoothParameters; fSmoothParameters = 0x0;
  delete fSmoothCovariances; fSmoothCovariances = 0x0;
  if(fOwnCluster) {
    delete fClusterPtr; fClusterPtr = 0x0;
  }
}

  //_________________________________________________________________________
void AliMUONTrackParam::GetParamFrom(const AliESDMuonTrack& esdMuonTrack)
{
  /// Get parameters from ESD track
  fZ = esdMuonTrack.GetZ(); 
  fParameters(0,0) = esdMuonTrack.GetNonBendingCoor();
  fParameters(1,0) = TMath::Tan(esdMuonTrack.GetThetaX());
  fParameters(2,0) = esdMuonTrack.GetBendingCoor(); 
  fParameters(3,0) = TMath::Tan(esdMuonTrack.GetThetaY());
  fParameters(4,0) = esdMuonTrack.GetInverseBendingMomentum();
}

  //_________________________________________________________________________
void AliMUONTrackParam::SetParamFor(AliESDMuonTrack& esdMuonTrack) const
{
  /// Set parameters in ESD track
  esdMuonTrack.SetZ(fZ);
  esdMuonTrack.SetNonBendingCoor(fParameters(0,0));
  esdMuonTrack.SetThetaX(TMath::ATan(fParameters(1,0)));
  esdMuonTrack.SetBendingCoor(fParameters(2,0)); 
  esdMuonTrack.SetThetaY(TMath::ATan(fParameters(3,0)));
  esdMuonTrack.SetInverseBendingMomentum(fParameters(4,0));
}

  //_________________________________________________________________________
void AliMUONTrackParam::GetParamFromUncorrected(const AliESDMuonTrack& esdMuonTrack)
{
  /// Get parameters from ESD track
  fZ = esdMuonTrack.GetZUncorrected(); 
  fParameters(0,0) = esdMuonTrack.GetNonBendingCoorUncorrected();
  fParameters(1,0) = TMath::Tan(esdMuonTrack.GetThetaXUncorrected());
  fParameters(2,0) = esdMuonTrack.GetBendingCoorUncorrected(); 
  fParameters(3,0) = TMath::Tan(esdMuonTrack.GetThetaYUncorrected());
  fParameters(4,0) = esdMuonTrack.GetInverseBendingMomentumUncorrected();
}

  //_________________________________________________________________________
void AliMUONTrackParam::SetParamForUncorrected(AliESDMuonTrack& esdMuonTrack) const
{
  /// Set parameters in ESD track
  esdMuonTrack.SetZUncorrected(fZ);
  esdMuonTrack.SetNonBendingCoorUncorrected(fParameters(0,0));
  esdMuonTrack.SetThetaXUncorrected(TMath::ATan(fParameters(1,0)));
  esdMuonTrack.SetBendingCoorUncorrected(fParameters(2,0)); 
  esdMuonTrack.SetThetaYUncorrected(TMath::ATan(fParameters(3,0)));
  esdMuonTrack.SetInverseBendingMomentumUncorrected(fParameters(4,0));
}

  //_________________________________________________________________________
void AliMUONTrackParam::GetCovFrom(const AliESDMuonTrack& esdMuonTrack)
{
  /// Get parameters covariances from ESD track
  
  // Get ESD covariance matrix
  if (!fCovariances) fCovariances = new TMatrixD(5,5);
  esdMuonTrack.GetCovariances(*fCovariances);

  // compute Jacobian to change the coordinate system
  // from (X,thetaX,Y,thetaY,c/pYZ) to (X,slopeX,Y,slopeY,c/pYZ)
  Double_t cosThetaX = TMath::Cos(TMath::ATan(fParameters(1,0)));
  Double_t cosThetaY = TMath::Cos(TMath::ATan(fParameters(3,0)));
  TMatrixD jacob(5,5);
  jacob.Zero();
  jacob(0,0) = 1.;
  jacob(1,1) = 1. / cosThetaX / cosThetaX;
  jacob(2,2) = 1.;
  jacob(3,3) = 1. / cosThetaY / cosThetaY;
  jacob(4,4) = 1.;
  
  // compute covariance matrix in ESD coordinate system
  TMatrixD tmp(*fCovariances,TMatrixD::kMultTranspose,jacob);
  *fCovariances = TMatrixD(jacob,TMatrixD::kMult,tmp);
  
}

  //_________________________________________________________________________
void AliMUONTrackParam::SetCovFor(AliESDMuonTrack& esdMuonTrack) const
{
  /// Set parameters covariances in ESD track
  
  // set null matrix if covariances does not exist
  if (!fCovariances) {
    TMatrixD tmp(5,5);
    tmp.Zero();
    esdMuonTrack.SetCovariances(tmp);
    return;
  }
  
  // compute Jacobian to change the coordinate system
  // from (X,slopeX,Y,slopeY,c/pYZ) to (X,thetaX,Y,thetaY,c/pYZ)
  Double_t cosThetaX = TMath::Cos(TMath::ATan(fParameters(1,0)));
  Double_t cosThetaY = TMath::Cos(TMath::ATan(fParameters(3,0)));
  TMatrixD jacob(5,5);
  jacob.Zero();
  jacob(0,0) = 1.;
  jacob(1,1) = cosThetaX * cosThetaX;
  jacob(2,2) = 1.;
  jacob(3,3) = cosThetaY * cosThetaY;
  jacob(4,4) = 1.;
  
  // compute covariance matrix in ESD coordinate system
  TMatrixD tmp(*fCovariances,TMatrixD::kMultTranspose,jacob);
  TMatrixD tmp2(jacob,TMatrixD::kMult,tmp);
  esdMuonTrack.SetCovariances(tmp2);

}

  //__________________________________________________________________________
Double_t AliMUONTrackParam::Px() const
{
  /// return p_x from track parameters
  Double_t pYZ, pZ, pX;
  pYZ = 0;
  if (  TMath::Abs(fParameters(4,0)) > 0 )
    pYZ = TMath::Abs(1.0 / fParameters(4,0));
  pZ = -pYZ / (TMath::Sqrt(1.0 + fParameters(3,0) * fParameters(3,0)));  // spectro. (z<0)
  pX = pZ * fParameters(1,0); 
  return pX;
}

  //__________________________________________________________________________
Double_t AliMUONTrackParam::Py() const
{
  /// return p_y from track parameters
  Double_t pYZ, pZ, pY;
  pYZ = 0;
  if (  TMath::Abs(fParameters(4,0)) > 0 )
    pYZ = TMath::Abs(1.0 / fParameters(4,0));
  pZ = -pYZ / (TMath::Sqrt(1.0 + fParameters(3,0) * fParameters(3,0)));  // spectro. (z<0)
  pY = pZ * fParameters(3,0); 
  return pY;
}

  //__________________________________________________________________________
Double_t AliMUONTrackParam::Pz() const
{
  /// return p_z from track parameters
  Double_t pYZ, pZ;
  pYZ = 0;
  if (  TMath::Abs(fParameters(4,0)) > 0 )
    pYZ = TMath::Abs(1.0 / fParameters(4,0));
  pZ = -pYZ / (TMath::Sqrt(1.0 + fParameters(3,0) * fParameters(3,0)));  // spectro. (z<0)
  return pZ;
}

  //__________________________________________________________________________
Double_t AliMUONTrackParam::P() const
{
  /// return p from track parameters
  Double_t  pYZ, pZ, p;
  pYZ = 0;
  if (  TMath::Abs(fParameters(4,0)) > 0 )
    pYZ = TMath::Abs(1.0 / fParameters(4,0));
  pZ = -pYZ / (TMath::Sqrt(1.0 + fParameters(3,0) * fParameters(3,0)));  // spectro. (z<0)
  p = TMath::Abs(pZ) * 
    TMath::Sqrt(1.0 + fParameters(3,0) * fParameters(3,0) + fParameters(1,0) * fParameters(1,0));
  return p;
  
}

  //__________________________________________________________________________
const TMatrixD& AliMUONTrackParam::GetCovariances() const
{
  /// Return the covariance matrix (create it before if needed)
  if (!fCovariances) {
    fCovariances = new TMatrixD(5,5);
    fCovariances->Zero();
  }
  return *fCovariances;
}

  //__________________________________________________________________________
void AliMUONTrackParam::SetCovariances(const TMatrixD& covariances)
{
  /// Set the covariance matrix
  if (fCovariances) *fCovariances = covariances;
  else fCovariances = new TMatrixD(covariances);
}

  //__________________________________________________________________________
void AliMUONTrackParam::SetCovariances(const Double_t matrix[5][5])
{
  /// Set the covariance matrix
  if (fCovariances) fCovariances->SetMatrixArray(&(matrix[0][0]));
  else fCovariances = new TMatrixD(5,5,&(matrix[0][0]));
}

  //__________________________________________________________________________
void AliMUONTrackParam::SetVariances(const Double_t matrix[5][5])
{
  /// Set the diagonal terms of the covariance matrix (variances)
  if (!fCovariances) fCovariances = new TMatrixD(5,5);
  fCovariances->Zero();
  for (Int_t i=0; i<5; i++) (*fCovariances)(i,i) = matrix[i][i];
}

  //__________________________________________________________________________
void AliMUONTrackParam::DeleteCovariances()
{
  /// Delete the covariance matrix
  delete fCovariances;
  fCovariances = 0x0;
}

  //__________________________________________________________________________
const TMatrixD& AliMUONTrackParam::GetPropagator() const
{
  /// Return the propagator (create it before if needed)
  if (!fPropagator) {
    fPropagator = new TMatrixD(5,5);
    fPropagator->UnitMatrix();
  }
  return *fPropagator;
}

  //__________________________________________________________________________
void AliMUONTrackParam::ResetPropagator()
{
  /// Reset the propagator
  if (fPropagator) fPropagator->UnitMatrix();
}

  //__________________________________________________________________________
void AliMUONTrackParam::UpdatePropagator(const TMatrixD& propagator)
{
  /// Update the propagator
  if (fPropagator) *fPropagator = TMatrixD(propagator,TMatrixD::kMult,*fPropagator);
  else fPropagator = new TMatrixD(propagator);
}

  //__________________________________________________________________________
const TMatrixD& AliMUONTrackParam::GetExtrapParameters() const
{
  /// Return extrapolated parameters (create it before if needed)
  if (!fExtrapParameters) {
    fExtrapParameters = new TMatrixD(5,1);
    fExtrapParameters->Zero();
  }
  return *fExtrapParameters;
  }

  //__________________________________________________________________________
void AliMUONTrackParam::SetExtrapParameters(const TMatrixD& extrapParameters)
{
  /// Set extrapolated parameters
  if (fExtrapParameters) *fExtrapParameters = extrapParameters;
  else fExtrapParameters = new TMatrixD(extrapParameters);
}

  //__________________________________________________________________________
const TMatrixD& AliMUONTrackParam::GetExtrapCovariances() const
{
  /// Return the extrapolated covariance matrix (create it before if needed)
  if (!fExtrapCovariances) {
    fExtrapCovariances = new TMatrixD(5,5);
    fExtrapCovariances->Zero();
  }
  return *fExtrapCovariances;
  }

  //__________________________________________________________________________
void AliMUONTrackParam::SetExtrapCovariances(const TMatrixD& extrapCovariances)
{
  /// Set the extrapolated covariance matrix
  if (fExtrapCovariances) *fExtrapCovariances = extrapCovariances;
  else fExtrapCovariances = new TMatrixD(extrapCovariances);
}

  //__________________________________________________________________________
const TMatrixD& AliMUONTrackParam::GetSmoothParameters() const
{
  /// Return the smoothed parameters (create it before if needed)
  if (!fSmoothParameters) {
    fSmoothParameters = new TMatrixD(5,1);
    fSmoothParameters->Zero();
  }
  return *fSmoothParameters;
  }

  //__________________________________________________________________________
void AliMUONTrackParam::SetSmoothParameters(const TMatrixD& smoothParameters)
{
  /// Set the smoothed parameters
  if (fSmoothParameters) *fSmoothParameters = smoothParameters;
  else fSmoothParameters = new TMatrixD(smoothParameters);
}

  //__________________________________________________________________________
const TMatrixD& AliMUONTrackParam::GetSmoothCovariances() const
{
  /// Return the smoothed covariance matrix (create it before if needed)
  if (!fSmoothCovariances) {
    fSmoothCovariances = new TMatrixD(5,5);
    fSmoothCovariances->Zero();
  }
  return *fSmoothCovariances;
  }

  //__________________________________________________________________________
void AliMUONTrackParam::SetSmoothCovariances(const TMatrixD& smoothCovariances)
{
  /// Set the smoothed covariance matrix
  if (fSmoothCovariances) *fSmoothCovariances = smoothCovariances;
  else fSmoothCovariances = new TMatrixD(smoothCovariances);
}

  //__________________________________________________________________________
Int_t AliMUONTrackParam::Compare(const TObject* trackParam) const
{
  /// "Compare" function to sort with decreasing Z (spectro. muon Z <0).
  /// Returns 1 (0, -1) if the current Z
  /// is smaller than (equal to, larger than) Z of trackParam
  if (fZ < ((AliMUONTrackParam*)trackParam)->GetZ()) return(1);
  else if (fZ == ((AliMUONTrackParam*)trackParam)->GetZ()) return(0);
  else return(-1);
}

//_____________________________________________-
void AliMUONTrackParam::Print(Option_t* opt) const
{
  /// Printing TrackParam information 
  /// "full" option for printing all the information about the TrackParam
  TString sopt(opt);
  sopt.ToUpper();
 
  if ( sopt.Contains("FULL") ) { 
    cout << "<AliMUONTrackParam> Bending P=" << setw(5) << setprecision(3)  << 1./fParameters(4,0) << 
      ", NonBendSlope=" << setw(5) << setprecision(3)  << fParameters(1,0)*180./TMath::Pi() <<
      ", BendSlope=" << setw(5) << setprecision(3)     << fParameters(3,0)*180./TMath::Pi()  << 
      ", (x,y,z)_IP=(" <<  setw(5) << setprecision(3) << fParameters(0,0) <<
      "," <<  setw(5) << setprecision(3) << fParameters(2,0) <<
      "," <<  setw(5) << setprecision(3) << fZ <<
      ") cm, (px,py,pz)=(" << setw(5) << setprecision(3) << Px() <<
      "," << setw(5) << setprecision(3) << Py() <<
      "," << setw(5) << setprecision(3) << Pz() << ") GeV/c" << endl;
  }
  else {
    cout << "<AliMUONTrackParam>"  << endl;
  }
    
}
