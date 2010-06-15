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
// Class AliMUONTrackExtrap
// ------------------------
// Tools for track extrapolation in ALICE dimuon spectrometer
// Author: Philippe Pillot
//-----------------------------------------------------------------------------

#include "AliMUONTrackExtrap.h" 
#include "AliMUONTrackParam.h"
#include "AliMUONConstants.h"
#include "AliMUONReconstructor.h"

#include "AliMagF.h"
#include "AliExternalTrackParam.h"

#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TMath.h>
#include <TDatabasePDG.h>

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrackExtrap) // Class implementation in ROOT context
/// \endcond

const Double_t AliMUONTrackExtrap::fgkSimpleBPosition = 0.5 * (AliMUONConstants::CoilZ() + AliMUONConstants::YokeZ());
const Double_t AliMUONTrackExtrap::fgkSimpleBLength = 0.5 * (AliMUONConstants::CoilL() + AliMUONConstants::YokeL());
      Double_t AliMUONTrackExtrap::fgSimpleBValue = 0.;
      Bool_t   AliMUONTrackExtrap::fgFieldON = kFALSE;
const Bool_t   AliMUONTrackExtrap::fgkUseHelix = kFALSE;
const Int_t    AliMUONTrackExtrap::fgkMaxStepNumber = 5000;
const Double_t AliMUONTrackExtrap::fgkHelixStepLength = 6.;
const Double_t AliMUONTrackExtrap::fgkRungeKuttaMaxResidue = 0.002;

//__________________________________________________________________________
void AliMUONTrackExtrap::SetField()
{
  /// set field on/off flag;  
  /// set field at the centre of the dipole
  const Double_t x[3] = {50.,50.,fgkSimpleBPosition};
  Double_t b[3] = {0.,0.,0.};
  TGeoGlobalMagField::Instance()->Field(x,b);
  fgSimpleBValue = b[0];
  fgFieldON = (TMath::Abs(fgSimpleBValue) > 1.e-10) ? kTRUE : kFALSE;
  
}

//__________________________________________________________________________
Double_t AliMUONTrackExtrap::GetImpactParamFromBendingMomentum(Double_t bendingMomentum)
{
  /// Returns impact parameter at vertex in bending plane (cm),
  /// from the signed bending momentum "BendingMomentum" in bending plane (GeV/c),
  /// using simple values for dipole magnetic field.
  /// The sign of "BendingMomentum" is the sign of the charge.
  
  if (bendingMomentum == 0.) return 1.e10;
  
  const Double_t kCorrectionFactor = 1.1; // impact parameter is 10% underestimated
  
  return kCorrectionFactor * (-0.0003 * fgSimpleBValue * fgkSimpleBLength * fgkSimpleBPosition / bendingMomentum);
}

//__________________________________________________________________________
Double_t 
AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(Double_t impactParam)
{
  /// Returns signed bending momentum in bending plane (GeV/c),
  /// the sign being the sign of the charge for particles moving forward in Z,
  /// from the impact parameter "ImpactParam" at vertex in bending plane (cm),
  /// using simple values for dipole magnetic field.
  
  if (impactParam == 0.) return 1.e10;
  
  const Double_t kCorrectionFactor = 1.1; // bending momentum is 10% underestimated
  
  if (fgFieldON) 
  {
    return kCorrectionFactor * (-0.0003 * fgSimpleBValue * fgkSimpleBLength * fgkSimpleBPosition / impactParam);
  }
  else 
  {
    return AliMUONConstants::GetMostProbBendingMomentum();
  }
}

//__________________________________________________________________________
void AliMUONTrackExtrap::LinearExtrapToZ(AliMUONTrackParam* trackParam, Double_t zEnd)
{
  /// Track parameters linearly extrapolated to the plane at "zEnd".
  /// On return, results from the extrapolation are updated in trackParam.
  
  if (trackParam->GetZ() == zEnd) return; // nothing to be done if same z
  
  // Compute track parameters
  Double_t dZ = zEnd - trackParam->GetZ();
  trackParam->SetNonBendingCoor(trackParam->GetNonBendingCoor() + trackParam->GetNonBendingSlope() * dZ);
  trackParam->SetBendingCoor(trackParam->GetBendingCoor() + trackParam->GetBendingSlope() * dZ);
  trackParam->SetZ(zEnd);
}

//__________________________________________________________________________
void AliMUONTrackExtrap::LinearExtrapToZCov(AliMUONTrackParam* trackParam, Double_t zEnd, Bool_t updatePropagator)
{
  /// Track parameters and their covariances linearly extrapolated to the plane at "zEnd".
  /// On return, results from the extrapolation are updated in trackParam.
  
  if (trackParam->GetZ() == zEnd) return; // nothing to be done if same z
  
  // No need to propagate the covariance matrix if it does not exist
  if (!trackParam->CovariancesExist()) {
    cout<<"W-AliMUONTrackExtrap::LinearExtrapToZCov: Covariance matrix does not exist"<<endl;
    // Extrapolate linearly track parameters to "zEnd"
    LinearExtrapToZ(trackParam,zEnd);
    return;
  }
  
  // Compute track parameters
  Double_t dZ = zEnd - trackParam->GetZ();
  trackParam->SetNonBendingCoor(trackParam->GetNonBendingCoor() + trackParam->GetNonBendingSlope() * dZ);
  trackParam->SetBendingCoor(trackParam->GetBendingCoor() + trackParam->GetBendingSlope() * dZ);
  trackParam->SetZ(zEnd);
  
  // Calculate the jacobian related to the track parameters linear extrapolation to "zEnd"
  TMatrixD jacob(5,5);
  jacob.UnitMatrix();
  jacob(0,1) = dZ;
  jacob(2,3) = dZ;
  
  // Extrapolate track parameter covariances to "zEnd"
  TMatrixD tmp(trackParam->GetCovariances(),TMatrixD::kMultTranspose,jacob);
  TMatrixD tmp2(jacob,TMatrixD::kMult,tmp);
  trackParam->SetCovariances(tmp2);
  
  // Update the propagator if required
  if (updatePropagator) trackParam->UpdatePropagator(jacob);
}

//__________________________________________________________________________
Bool_t AliMUONTrackExtrap::ExtrapToZ(AliMUONTrackParam* trackParam, Double_t zEnd)
{
  /// Interface to track parameter extrapolation to the plane at "Z" using Helix or Rungekutta algorithm.
  /// On return, the track parameters resulting from the extrapolation are updated in trackParam.
  if (!fgFieldON) {
    AliMUONTrackExtrap::LinearExtrapToZ(trackParam,zEnd);
    return kTRUE;
  }
  else if (fgkUseHelix) return AliMUONTrackExtrap::ExtrapToZHelix(trackParam,zEnd);
  else return AliMUONTrackExtrap::ExtrapToZRungekutta(trackParam,zEnd);
}

//__________________________________________________________________________
Bool_t AliMUONTrackExtrap::ExtrapToZHelix(AliMUONTrackParam* trackParam, Double_t zEnd)
{
  /// Track parameter extrapolation to the plane at "Z" using Helix algorithm.
  /// On return, the track parameters resulting from the extrapolation are updated in trackParam.
  if (trackParam->GetZ() == zEnd) return kTRUE; // nothing to be done if same Z
  Double_t forwardBackward; // +1 if forward, -1 if backward
  if (zEnd < trackParam->GetZ()) forwardBackward = 1.0; // spectro. z<0 
  else forwardBackward = -1.0;
  Double_t v3[7], v3New[7]; // 7 in parameter ????
  Int_t i3, stepNumber;
  // For safety: return kTRUE or kFALSE ????
  // Parameter vector for calling EXTRAP_ONESTEP
  ConvertTrackParamForExtrap(trackParam, forwardBackward, v3);
  // sign of charge (sign of fInverseBendingMomentum if forward motion)
  // must be changed if backward extrapolation
  Double_t chargeExtrap = forwardBackward * TMath::Sign(Double_t(1.0), trackParam->GetInverseBendingMomentum());
  // Extrapolation loop
  stepNumber = 0;
  while (((-forwardBackward * (v3[2] - zEnd)) <= 0.0) && (stepNumber < fgkMaxStepNumber)) { // spectro. z<0
    stepNumber++;
    ExtrapOneStepHelix(chargeExtrap, fgkHelixStepLength, v3, v3New);
    if ((-forwardBackward * (v3New[2] - zEnd)) > 0.0) break; // one is beyond Z spectro. z<0
							     // better use TArray ????
    for (i3 = 0; i3 < 7; i3++) {v3[i3] = v3New[i3];}
  }
  // check fgkMaxStepNumber ????
  // Interpolation back to exact Z (2nd order)
  // should be in function ???? using TArray ????
  Double_t dZ12 = v3New[2] - v3[2]; // 1->2
  if (TMath::Abs(dZ12) > 0) {
    Double_t dZ1i = zEnd - v3[2]; // 1-i
    Double_t dZi2 = v3New[2] - zEnd; // i->2
    Double_t xPrime = (v3New[0] - v3[0]) / dZ12;
    Double_t xSecond = ((v3New[3] / v3New[5]) - (v3[3] / v3[5])) / dZ12;
    Double_t yPrime = (v3New[1] - v3[1]) / dZ12;
    Double_t ySecond = ((v3New[4] / v3New[5]) - (v3[4] / v3[5])) / dZ12;
    v3[0] = v3[0] + xPrime * dZ1i - 0.5 * xSecond * dZ1i * dZi2; // X
    v3[1] = v3[1] + yPrime * dZ1i - 0.5 * ySecond * dZ1i * dZi2; // Y
    v3[2] = zEnd; // Z
    Double_t xPrimeI = xPrime - 0.5 * xSecond * (dZi2 - dZ1i);
    Double_t yPrimeI = yPrime - 0.5 * ySecond * (dZi2 - dZ1i);
    // (PX, PY, PZ)/PTOT assuming forward motion
    v3[5] = 1.0 / TMath::Sqrt(1.0 + xPrimeI * xPrimeI + yPrimeI * yPrimeI); // PZ/PTOT
    v3[3] = xPrimeI * v3[5]; // PX/PTOT
    v3[4] = yPrimeI * v3[5]; // PY/PTOT
  } else {
    cout<<"W-AliMUONTrackExtrap::ExtrapToZHelix: Extrap. to Z not reached, Z = "<<zEnd<<endl;
  }
  // Recover track parameters (charge back for forward motion)
  RecoverTrackParam(v3, chargeExtrap * forwardBackward, trackParam);
  return kTRUE;
}

//__________________________________________________________________________
Bool_t AliMUONTrackExtrap::ExtrapToZRungekutta(AliMUONTrackParam* trackParam, Double_t zEnd)
{
  /// Track parameter extrapolation to the plane at "Z" using Rungekutta algorithm.
  /// On return, the track parameters resulting from the extrapolation are updated in trackParam.
  if (trackParam->GetZ() == zEnd) return kTRUE; // nothing to be done if same Z
  Double_t forwardBackward; // +1 if forward, -1 if backward
  if (zEnd < trackParam->GetZ()) forwardBackward = 1.0; // spectro. z<0 
  else forwardBackward = -1.0;
  // sign of charge (sign of fInverseBendingMomentum if forward motion)
  // must be changed if backward extrapolation
  Double_t chargeExtrap = forwardBackward * TMath::Sign(Double_t(1.0), trackParam->GetInverseBendingMomentum());
  Double_t v3[7], v3New[7];
  Double_t dZ, step;
  Int_t stepNumber = 0;
  
  // Extrapolation loop (until within tolerance or the track turn around)
  Double_t residue = zEnd - trackParam->GetZ();
  Bool_t uturn = kFALSE;
  Bool_t trackingFailed = kFALSE;
  Bool_t tooManyStep = kFALSE;
  while (TMath::Abs(residue) > fgkRungeKuttaMaxResidue && stepNumber <= fgkMaxStepNumber) {
    
    dZ = zEnd - trackParam->GetZ();
    // step lenght assuming linear trajectory
    step = dZ * TMath::Sqrt(1.0 + trackParam->GetBendingSlope()*trackParam->GetBendingSlope() +
			    trackParam->GetNonBendingSlope()*trackParam->GetNonBendingSlope());
    ConvertTrackParamForExtrap(trackParam, forwardBackward, v3);
    
    do { // reduce step lenght while zEnd oversteped
      if (stepNumber > fgkMaxStepNumber) {
        cout<<"W-AliMUONTrackExtrap::ExtrapToZRungekutta: Too many trials: "<<stepNumber<<endl;
	tooManyStep = kTRUE;
	break;
      }
      stepNumber ++;
      step = TMath::Abs(step);
      if (!AliMUONTrackExtrap::ExtrapOneStepRungekutta(chargeExtrap,step,v3,v3New)) {
	trackingFailed = kTRUE;
	break;
      }
      residue = zEnd - v3New[2];
      step *= dZ/(v3New[2]-trackParam->GetZ());
    } while (residue*dZ < 0 && TMath::Abs(residue) > fgkRungeKuttaMaxResidue);
    
    if (trackingFailed) break;
    else if (v3New[5]*v3[5] < 0) { // the track turned around
      cout<<"W-AliMUONTrackExtrap::ExtrapToZRungekutta: The track turned around"<<endl;
      uturn = kTRUE;
      break;
    } else RecoverTrackParam(v3New, chargeExtrap * forwardBackward, trackParam);
    
  }
  
  // terminate the extropolation with a straight line up to the exact "zEnd" value
  if (trackingFailed || uturn) {
    
    // track ends +-100 meters away in the bending direction
    dZ = zEnd - v3[2];
    Double_t bendingSlope = TMath::Sign(1.e4,-fgSimpleBValue*trackParam->GetInverseBendingMomentum()) / dZ;
    Double_t pZ = TMath::Abs(1. / trackParam->GetInverseBendingMomentum()) / TMath::Sqrt(1.0 + bendingSlope * bendingSlope);
    Double_t nonBendingSlope = TMath::Sign(TMath::Abs(v3[3]) * v3[6] / pZ, trackParam->GetNonBendingSlope());
    trackParam->SetNonBendingCoor(trackParam->GetNonBendingCoor() + dZ * nonBendingSlope);
    trackParam->SetNonBendingSlope(nonBendingSlope);
    trackParam->SetBendingCoor(trackParam->GetBendingCoor() + dZ * bendingSlope);
    trackParam->SetBendingSlope(bendingSlope);
    trackParam->SetZ(zEnd);
    
    return kFALSE;
    
  } else {
    
    // track extrapolated normally
    trackParam->SetNonBendingCoor(trackParam->GetNonBendingCoor() + residue * trackParam->GetNonBendingSlope());
    trackParam->SetBendingCoor(trackParam->GetBendingCoor() + residue * trackParam->GetBendingSlope());
    trackParam->SetZ(zEnd);
    
    return !tooManyStep;
    
  }
  
}

//__________________________________________________________________________
void AliMUONTrackExtrap::ConvertTrackParamForExtrap(AliMUONTrackParam* trackParam, Double_t forwardBackward, Double_t *v3)
{
  /// Set vector of Geant3 parameters pointed to by "v3" from track parameters in trackParam.
  /// Since AliMUONTrackParam is only geometry, one uses "forwardBackward"
  /// to know whether the particle is going forward (+1) or backward (-1).
  v3[0] = trackParam->GetNonBendingCoor(); // X
  v3[1] = trackParam->GetBendingCoor(); // Y
  v3[2] = trackParam->GetZ(); // Z
  Double_t pYZ = TMath::Abs(1.0 / trackParam->GetInverseBendingMomentum());
  Double_t pZ = pYZ / TMath::Sqrt(1.0 + trackParam->GetBendingSlope() * trackParam->GetBendingSlope());
  v3[6] = TMath::Sqrt(pYZ * pYZ + pZ * pZ * trackParam->GetNonBendingSlope() * trackParam->GetNonBendingSlope()); // PTOT
  v3[5] = -forwardBackward * pZ / v3[6]; // PZ/PTOT spectro. z<0
  v3[3] = trackParam->GetNonBendingSlope() * v3[5]; // PX/PTOT
  v3[4] = trackParam->GetBendingSlope() * v3[5]; // PY/PTOT
}

//__________________________________________________________________________
void AliMUONTrackExtrap::RecoverTrackParam(Double_t *v3, Double_t charge, AliMUONTrackParam* trackParam)
{
  /// Set track parameters in trackParam from Geant3 parameters pointed to by "v3",
  /// assumed to be calculated for forward motion in Z.
  /// "InverseBendingMomentum" is signed with "charge".
  trackParam->SetNonBendingCoor(v3[0]); // X
  trackParam->SetBendingCoor(v3[1]); // Y
  trackParam->SetZ(v3[2]); // Z
  Double_t pYZ = v3[6] * TMath::Sqrt((1.-v3[3])*(1.+v3[3]));
  trackParam->SetInverseBendingMomentum(charge/pYZ);
  trackParam->SetBendingSlope(v3[4]/v3[5]);
  trackParam->SetNonBendingSlope(v3[3]/v3[5]);
}

//__________________________________________________________________________
Bool_t AliMUONTrackExtrap::ExtrapToZCov(AliMUONTrackParam* trackParam, Double_t zEnd, Bool_t updatePropagator)
{
  /// Track parameters and their covariances extrapolated to the plane at "zEnd".
  /// On return, results from the extrapolation are updated in trackParam.
  
  if (trackParam->GetZ() == zEnd) return kTRUE; // nothing to be done if same z
  
  if (!fgFieldON) { // linear extrapolation if no magnetic field
    AliMUONTrackExtrap::LinearExtrapToZCov(trackParam,zEnd,updatePropagator);
    return kTRUE;
  }
  
  // No need to propagate the covariance matrix if it does not exist
  if (!trackParam->CovariancesExist()) {
    cout<<"W-AliMUONTrackExtrap::ExtrapToZCov: Covariance matrix does not exist"<<endl;
    // Extrapolate track parameters to "zEnd"
    return ExtrapToZ(trackParam,zEnd);
  }
  
  // Save the actual track parameters
  AliMUONTrackParam trackParamSave(*trackParam);
  TMatrixD paramSave(trackParamSave.GetParameters());
  Double_t zBegin = trackParamSave.GetZ();
  
  // Get reference to the parameter covariance matrix
  const TMatrixD& kParamCov = trackParam->GetCovariances();
	
  // Extrapolate track parameters to "zEnd"
  // Do not update the covariance matrix if the extrapolation failed
  if (!ExtrapToZ(trackParam,zEnd)) return kFALSE;
  
  // Get reference to the extrapolated parameters
  const TMatrixD& extrapParam = trackParam->GetParameters();
  
  // Calculate the jacobian related to the track parameters extrapolation to "zEnd"
  Bool_t extrapStatus = kTRUE;
  TMatrixD jacob(5,5);
  jacob.Zero();
  TMatrixD dParam(5,1);
  Double_t direction[5] = {-1.,-1.,1.,1.,-1.};
  for (Int_t i=0; i<5; i++) {
    // Skip jacobian calculation for parameters with no associated error
    if (kParamCov(i,i) <= 0.) continue;
    
    // Small variation of parameter i only
    for (Int_t j=0; j<5; j++) {
      if (j==i) {
        dParam(j,0) = TMath::Sqrt(kParamCov(i,i));
	dParam(j,0) *= TMath::Sign(1.,direction[j]*paramSave(j,0)); // variation always in the same direction
      } else dParam(j,0) = 0.;
    }
    
    // Set new parameters
    trackParamSave.SetParameters(paramSave);
    trackParamSave.AddParameters(dParam);
    trackParamSave.SetZ(zBegin);
    
    // Extrapolate new track parameters to "zEnd"
    if (!ExtrapToZ(&trackParamSave,zEnd)) {
      cout<<"W-AliMUONTrackExtrap::ExtrapToZCov: Bad covariance matrix"<<endl;
      extrapStatus = kFALSE;
    }
    
    // Calculate the jacobian
    TMatrixD jacobji(trackParamSave.GetParameters(),TMatrixD::kMinus,extrapParam);
    jacobji *= 1. / dParam(i,0);
    jacob.SetSub(0,i,jacobji);
  }
  
  // Extrapolate track parameter covariances to "zEnd"
  TMatrixD tmp(kParamCov,TMatrixD::kMultTranspose,jacob);
  TMatrixD tmp2(jacob,TMatrixD::kMult,tmp);
  trackParam->SetCovariances(tmp2);
  
  // Update the propagator if required
  if (updatePropagator) trackParam->UpdatePropagator(jacob);
  
  return extrapStatus;
}

//__________________________________________________________________________
void AliMUONTrackExtrap::AddMCSEffectInAbsorber(AliMUONTrackParam* param, Double_t signedPathLength, Double_t f0, Double_t f1, Double_t f2)
{
  /// Add to the track parameter covariances the effects of multiple Coulomb scattering
  /// signedPathLength must have the sign of (zOut - zIn) where all other parameters are assumed to be given at zOut.
  
  // absorber related covariance parameters
  Double_t bendingSlope = param->GetBendingSlope();
  Double_t nonBendingSlope = param->GetNonBendingSlope();
  Double_t inverseBendingMomentum = param->GetInverseBendingMomentum();
  Double_t alpha2 = 0.0136 * 0.0136 * inverseBendingMomentum * inverseBendingMomentum * (1.0 + bendingSlope * bendingSlope) /
                    (1.0 + bendingSlope *bendingSlope + nonBendingSlope * nonBendingSlope); // velocity = 1
  Double_t pathLength = TMath::Abs(signedPathLength);
  Double_t varCoor = alpha2 * (pathLength * pathLength * f0 - 2. * pathLength * f1 + f2);
  Double_t covCorrSlope = TMath::Sign(1.,signedPathLength) * alpha2 * (pathLength * f0 - f1);
  Double_t varSlop = alpha2 * f0;
  
  // Set MCS covariance matrix
  TMatrixD newParamCov(param->GetCovariances());
  // Non bending plane
  newParamCov(0,0) += varCoor;       newParamCov(0,1) += covCorrSlope;
  newParamCov(1,0) += covCorrSlope;  newParamCov(1,1) += varSlop;
  // Bending plane
  newParamCov(2,2) += varCoor;       newParamCov(2,3) += covCorrSlope;
  newParamCov(3,2) += covCorrSlope;  newParamCov(3,3) += varSlop;
  
  // Set momentum related covariances if B!=0
  if (fgFieldON) {
    // compute derivative d(q/Pxy) / dSlopeX and d(q/Pxy) / dSlopeY
    Double_t dqPxydSlopeX = inverseBendingMomentum * nonBendingSlope / (1. + nonBendingSlope*nonBendingSlope + bendingSlope*bendingSlope);
    Double_t dqPxydSlopeY = - inverseBendingMomentum * nonBendingSlope*nonBendingSlope * bendingSlope /
                              (1. + bendingSlope*bendingSlope) / (1. + nonBendingSlope*nonBendingSlope + bendingSlope*bendingSlope);
    // Inverse bending momentum (due to dependences with bending and non bending slopes)
    newParamCov(4,0) += dqPxydSlopeX * covCorrSlope; newParamCov(0,4) += dqPxydSlopeX * covCorrSlope;
    newParamCov(4,1) += dqPxydSlopeX * varSlop;      newParamCov(1,4) += dqPxydSlopeX * varSlop;
    newParamCov(4,2) += dqPxydSlopeY * covCorrSlope; newParamCov(2,4) += dqPxydSlopeY * covCorrSlope;
    newParamCov(4,3) += dqPxydSlopeY * varSlop;      newParamCov(3,4) += dqPxydSlopeY * varSlop;
    newParamCov(4,4) += (dqPxydSlopeX*dqPxydSlopeX + dqPxydSlopeY*dqPxydSlopeY) * varSlop;
  }
  
  // Set new covariances
  param->SetCovariances(newParamCov);
}

//__________________________________________________________________________
void AliMUONTrackExtrap::CorrectMCSEffectInAbsorber(AliMUONTrackParam* param,
						    Double_t xVtx, Double_t yVtx, Double_t zVtx,
						    Double_t errXVtx, Double_t errYVtx,
						    Double_t absZBeg, Double_t pathLength, Double_t f0, Double_t f1, Double_t f2)
{
  /// Correct parameters and corresponding covariances using Branson correction
  /// - input param are parameters and covariances at the end of absorber
  /// - output param are parameters and covariances at vertex
  /// Absorber correction parameters are supposed to be calculated at the current track z-position
  
  // Position of the Branson plane (spectro. (z<0))
  Double_t zB = (f1>0.) ? absZBeg - f2/f1 : 0.;
  
  // Add MCS effects to current parameter covariances (spectro. (z<0))
  AddMCSEffectInAbsorber(param, -pathLength, f0, f1, f2);
  
  // Get track parameters and covariances in the Branson plane corrected for magnetic field effect
  ExtrapToZCov(param,zVtx);
  LinearExtrapToZCov(param,zB);
  
  // compute track parameters at vertex
  TMatrixD newParam(5,1);
  newParam(0,0) = xVtx;
  newParam(1,0) = (param->GetNonBendingCoor() - xVtx) / (zB - zVtx);
  newParam(2,0) = yVtx;
  newParam(3,0) = (param->GetBendingCoor() - yVtx) / (zB - zVtx);
  newParam(4,0) = param->GetCharge() / param->P() *
                  TMath::Sqrt(1.0 + newParam(1,0)*newParam(1,0) + newParam(3,0)*newParam(3,0)) /
		  TMath::Sqrt(1.0 + newParam(3,0)*newParam(3,0));
  
  // Get covariances in (X, SlopeX, Y, SlopeY, q*PTot) coordinate system
  TMatrixD paramCovP(param->GetCovariances());
  Cov2CovP(param->GetParameters(),paramCovP);
  
  // Get the covariance matrix in the (XVtx, X, YVtx, Y, q*PTot) coordinate system
  TMatrixD paramCovVtx(5,5);
  paramCovVtx.Zero();
  paramCovVtx(0,0) = errXVtx * errXVtx;
  paramCovVtx(1,1) = paramCovP(0,0);
  paramCovVtx(2,2) = errYVtx * errYVtx;
  paramCovVtx(3,3) = paramCovP(2,2);
  paramCovVtx(4,4) = paramCovP(4,4);
  paramCovVtx(1,3) = paramCovP(0,2);
  paramCovVtx(3,1) = paramCovP(2,0);
  paramCovVtx(1,4) = paramCovP(0,4);
  paramCovVtx(4,1) = paramCovP(4,0);
  paramCovVtx(3,4) = paramCovP(2,4);
  paramCovVtx(4,3) = paramCovP(4,2);
  
  // Jacobian of the transformation (XVtx, X, YVtx, Y, q*PTot) -> (XVtx, SlopeXVtx, YVtx, SlopeYVtx, q*PTotVtx)
  TMatrixD jacob(5,5);
  jacob.UnitMatrix();
  jacob(1,0) = - 1. / (zB - zVtx);
  jacob(1,1) = 1. / (zB - zVtx);
  jacob(3,2) = - 1. / (zB - zVtx);
  jacob(3,3) = 1. / (zB - zVtx);
  
  // Compute covariances at vertex in the (XVtx, SlopeXVtx, YVtx, SlopeYVtx, q*PTotVtx) coordinate system
  TMatrixD tmp(paramCovVtx,TMatrixD::kMultTranspose,jacob);
  TMatrixD newParamCov(jacob,TMatrixD::kMult,tmp);
  
  // Compute covariances at vertex in the (XVtx, SlopeXVtx, YVtx, SlopeYVtx, q/PyzVtx) coordinate system
  CovP2Cov(newParam,newParamCov);
  
  // Set parameters and covariances at vertex
  param->SetParameters(newParam);
  param->SetZ(zVtx);
  param->SetCovariances(newParamCov);
}

//__________________________________________________________________________
void AliMUONTrackExtrap::CorrectELossEffectInAbsorber(AliMUONTrackParam* param, Double_t eLoss, Double_t sigmaELoss2)
{
  /// Correct parameters for energy loss and add energy loss fluctuation effect to covariances
  
  // Get parameter covariances in (X, SlopeX, Y, SlopeY, q*PTot) coordinate system
  TMatrixD newParamCov(param->GetCovariances());
  Cov2CovP(param->GetParameters(),newParamCov);
  
  // Compute new parameters corrected for energy loss
  Double_t muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); // GeV
  Double_t p = param->P();
  Double_t e = TMath::Sqrt(p*p + muMass*muMass);
  Double_t eCorr = e + eLoss;
  Double_t pCorr = TMath::Sqrt(eCorr*eCorr - muMass*muMass);
  Double_t nonBendingSlope = param->GetNonBendingSlope();
  Double_t bendingSlope = param->GetBendingSlope();
  param->SetInverseBendingMomentum(param->GetCharge() / pCorr *
				   TMath::Sqrt(1.0 + nonBendingSlope*nonBendingSlope + bendingSlope*bendingSlope) /
				   TMath::Sqrt(1.0 + bendingSlope*bendingSlope));
  
  // Add effects of energy loss fluctuation to covariances
  newParamCov(4,4) += eCorr * eCorr / pCorr / pCorr * sigmaELoss2;
  
  // Get new parameter covariances in (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
  CovP2Cov(param->GetParameters(),newParamCov);
  
  // Set new parameter covariances
  param->SetCovariances(newParamCov);
}

//__________________________________________________________________________
Bool_t AliMUONTrackExtrap::GetAbsorberCorrectionParam(Double_t trackXYZIn[3], Double_t trackXYZOut[3], Double_t pTotal,
						      Double_t &pathLength, Double_t &f0, Double_t &f1, Double_t &f2,
						      Double_t &meanRho, Double_t &totalELoss, Double_t &sigmaELoss2)
{
  /// Parameters used to correct for Multiple Coulomb Scattering and energy loss in absorber
  /// Calculated assuming a linear propagation from trackXYZIn to trackXYZOut (order is important)
  // pathLength: path length between trackXYZIn and trackXYZOut (cm)
  // f0:         0th moment of z calculated with the inverse radiation-length distribution
  // f1:         1st moment of z calculated with the inverse radiation-length distribution
  // f2:         2nd moment of z calculated with the inverse radiation-length distribution
  // meanRho:    average density of crossed material (g/cm3)
  // totalELoss: total energy loss in absorber
  
  // Reset absorber's parameters
  pathLength = 0.;
  f0 = 0.;
  f1 = 0.;
  f2 = 0.;
  meanRho = 0.;
  totalELoss = 0.;
  sigmaELoss2 = 0.;
  
  // Check whether the geometry is available
  if (!gGeoManager) {
    cout<<"E-AliMUONTrackExtrap::GetAbsorberCorrectionParam: no TGeo"<<endl;
    return kFALSE;
  }
  
  // Initialize starting point and direction
  pathLength = TMath::Sqrt((trackXYZOut[0] - trackXYZIn[0])*(trackXYZOut[0] - trackXYZIn[0])+
			   (trackXYZOut[1] - trackXYZIn[1])*(trackXYZOut[1] - trackXYZIn[1])+
			   (trackXYZOut[2] - trackXYZIn[2])*(trackXYZOut[2] - trackXYZIn[2]));
  if (pathLength < TGeoShape::Tolerance()) return kFALSE;
  Double_t b[3];
  b[0] = (trackXYZOut[0] - trackXYZIn[0]) / pathLength;
  b[1] = (trackXYZOut[1] - trackXYZIn[1]) / pathLength;
  b[2] = (trackXYZOut[2] - trackXYZIn[2]) / pathLength;
  TGeoNode *currentnode = gGeoManager->InitTrack(trackXYZIn, b);
  if (!currentnode) {
    cout<<"E-AliMUONTrackExtrap::GetAbsorberCorrectionParam: start point out of geometry"<<endl;
    return kFALSE;
  }
  
  // loop over absorber slices and calculate absorber's parameters
  Double_t rho = 0.; // material density (g/cm3)
  Double_t x0 = 0.;  // radiation-length (cm-1)
  Double_t atomicA = 0.; // A of material
  Double_t atomicZ = 0.; // Z of material
  Double_t atomicZoverA = 0.; // Z/A of material
  Double_t localPathLength = 0;
  Double_t remainingPathLength = pathLength;
  Double_t zB = trackXYZIn[2];
  Double_t zE, dzB, dzE;
  do {
    // Get material properties
    TGeoMaterial *material = currentnode->GetVolume()->GetMedium()->GetMaterial();
    rho = material->GetDensity();
    x0 = material->GetRadLen();
    atomicA = material->GetA();
    atomicZ = material->GetZ();
    if(material->IsMixture()){
      TGeoMixture * mixture = (TGeoMixture*)material;
      atomicZoverA = 0.;
      Double_t sum = 0.;
      for (Int_t iel=0;iel<mixture->GetNelements();iel++){
	sum  += mixture->GetWmixt()[iel];
	atomicZoverA += mixture->GetZmixt()[iel]*mixture->GetWmixt()[iel]/mixture->GetAmixt()[iel];
      }
      atomicZoverA/=sum;
    }
    else atomicZoverA = atomicZ/atomicA;
    
    // Get path length within this material
    gGeoManager->FindNextBoundary(remainingPathLength);
    localPathLength = gGeoManager->GetStep() + 1.e-6;
    // Check if boundary within remaining path length. If so, make sure to cross the boundary to prepare the next step
    if (localPathLength >= remainingPathLength) localPathLength = remainingPathLength;
    else {
      currentnode = gGeoManager->Step();
      if (!currentnode) {
        cout<<"E-AliMUONTrackExtrap::GetAbsorberCorrectionParam: navigation failed"<<endl;
	f0 = f1 = f2 = meanRho = totalELoss = sigmaELoss2 = 0.;
	return kFALSE;
      }
      if (!gGeoManager->IsEntering()) {
        // make another small step to try to enter in new absorber slice
        gGeoManager->SetStep(0.001);
	currentnode = gGeoManager->Step();
	if (!gGeoManager->IsEntering() || !currentnode) {
          cout<<"E-AliMUONTrackExtrap::GetAbsorberCorrectionParam: navigation failed"<<endl;
	  f0 = f1 = f2 = meanRho = totalELoss = sigmaELoss2 = 0.;
	  return kFALSE;
	}
        localPathLength += 0.001;
      }
    }
    
    // calculate absorber's parameters
    zE = b[2] * localPathLength + zB;
    dzB = zB - trackXYZIn[2];
    dzE = zE - trackXYZIn[2];
    f0 += localPathLength / x0;
    f1 += (dzE*dzE - dzB*dzB) / b[2] / b[2] / x0 / 2.;
    f2 += (dzE*dzE*dzE - dzB*dzB*dzB) / b[2] / b[2] / b[2] / x0 / 3.;
    meanRho += localPathLength * rho;
    totalELoss += BetheBloch(pTotal, localPathLength, rho, atomicZ, atomicZoverA);
    sigmaELoss2 += EnergyLossFluctuation2(pTotal, localPathLength, rho, atomicZoverA);
    
    // prepare next step
    zB = zE;
    remainingPathLength -= localPathLength;
  } while (remainingPathLength > TGeoShape::Tolerance());
  
  meanRho /= pathLength;
  
  return kTRUE;
}

//__________________________________________________________________________
Double_t AliMUONTrackExtrap::GetMCSAngle2(const AliMUONTrackParam& param, Double_t dZ, Double_t x0)
{
  /// Return the angular dispersion square due to multiple Coulomb scattering
  /// through a material of thickness "dZ" and of radiation length "x0"
  /// assuming linear propagation and using the small angle approximation.
  
  Double_t bendingSlope = param.GetBendingSlope();
  Double_t nonBendingSlope = param.GetNonBendingSlope();
  Double_t inverseTotalMomentum2 = param.GetInverseBendingMomentum() * param.GetInverseBendingMomentum() *
                                   (1.0 + bendingSlope * bendingSlope) /
                                   (1.0 + bendingSlope *bendingSlope + nonBendingSlope * nonBendingSlope); 
  // Path length in the material
  Double_t pathLength = TMath::Abs(dZ) * TMath::Sqrt(1.0 + bendingSlope*bendingSlope + nonBendingSlope*nonBendingSlope);
  // relativistic velocity
  Double_t velo = 1.;
  // Angular dispersion square of the track (variance) in a plane perpendicular to the trajectory
  Double_t theta02 = 0.0136 / velo * (1 + 0.038 * TMath::Log(pathLength/x0));
  
  return theta02 * theta02 * inverseTotalMomentum2 * pathLength / x0;
}

//__________________________________________________________________________
void AliMUONTrackExtrap::AddMCSEffect(AliMUONTrackParam *param, Double_t dZ, Double_t x0)
{
  /// Add to the track parameter covariances the effects of multiple Coulomb scattering
  /// through a material of thickness "Abs(dZ)" and of radiation length "x0"
  /// assuming linear propagation and using the small angle approximation.
  /// dZ = zOut - zIn (sign is important) and "param" is assumed to be given zOut.
  /// If x0 <= 0., assume dZ = pathLength/x0 and consider the material thickness as negligible.
  
  Double_t bendingSlope = param->GetBendingSlope();
  Double_t nonBendingSlope = param->GetNonBendingSlope();
  Double_t inverseBendingMomentum = param->GetInverseBendingMomentum();
  Double_t inverseTotalMomentum2 = inverseBendingMomentum * inverseBendingMomentum *
                                   (1.0 + bendingSlope * bendingSlope) /
                                   (1.0 + bendingSlope *bendingSlope + nonBendingSlope * nonBendingSlope); 
  // Path length in the material
  Double_t signedPathLength = dZ * TMath::Sqrt(1.0 + bendingSlope*bendingSlope + nonBendingSlope*nonBendingSlope);
  Double_t pathLengthOverX0 = (x0 > 0.) ? TMath::Abs(signedPathLength) / x0 : TMath::Abs(signedPathLength);
  // relativistic velocity
  Double_t velo = 1.;
  // Angular dispersion square of the track (variance) in a plane perpendicular to the trajectory
  Double_t theta02 = 0.0136 / velo * (1 + 0.038 * TMath::Log(pathLengthOverX0));
  theta02 *= theta02 * inverseTotalMomentum2 * pathLengthOverX0;
  
  Double_t varCoor 	= (x0 > 0.) ? signedPathLength * signedPathLength * theta02 / 3. : 0.;
  Double_t varSlop 	= theta02;
  Double_t covCorrSlope = (x0 > 0.) ? signedPathLength * theta02 / 2. : 0.;
  
  // Set MCS covariance matrix
  TMatrixD newParamCov(param->GetCovariances());
  // Non bending plane
  newParamCov(0,0) += varCoor;       newParamCov(0,1) += covCorrSlope;
  newParamCov(1,0) += covCorrSlope;  newParamCov(1,1) += varSlop;
  // Bending plane
  newParamCov(2,2) += varCoor;       newParamCov(2,3) += covCorrSlope;
  newParamCov(3,2) += covCorrSlope;  newParamCov(3,3) += varSlop;
  
  // Set momentum related covariances if B!=0
  if (fgFieldON) {
    // compute derivative d(q/Pxy) / dSlopeX and d(q/Pxy) / dSlopeY
    Double_t dqPxydSlopeX = inverseBendingMomentum * nonBendingSlope / (1. + nonBendingSlope*nonBendingSlope + bendingSlope*bendingSlope);
    Double_t dqPxydSlopeY = - inverseBendingMomentum * nonBendingSlope*nonBendingSlope * bendingSlope /
                              (1. + bendingSlope*bendingSlope) / (1. + nonBendingSlope*nonBendingSlope + bendingSlope*bendingSlope);
    // Inverse bending momentum (due to dependences with bending and non bending slopes)
    newParamCov(4,0) += dqPxydSlopeX * covCorrSlope; newParamCov(0,4) += dqPxydSlopeX * covCorrSlope;
    newParamCov(4,1) += dqPxydSlopeX * varSlop;      newParamCov(1,4) += dqPxydSlopeX * varSlop;
    newParamCov(4,2) += dqPxydSlopeY * covCorrSlope; newParamCov(2,4) += dqPxydSlopeY * covCorrSlope;
    newParamCov(4,3) += dqPxydSlopeY * varSlop;      newParamCov(3,4) += dqPxydSlopeY * varSlop;
    newParamCov(4,4) += (dqPxydSlopeX*dqPxydSlopeX + dqPxydSlopeY*dqPxydSlopeY) * varSlop;
  }
  
  // Set new covariances
  param->SetCovariances(newParamCov);
}

//__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToVertex(AliMUONTrackParam* trackParam,
					Double_t xVtx, Double_t yVtx, Double_t zVtx,
					Double_t errXVtx, Double_t errYVtx,
					Bool_t correctForMCS, Bool_t correctForEnergyLoss)
{
  /// Main method for extrapolation to the vertex:
  /// Returns the track parameters and covariances resulting from the extrapolation of the current trackParam
  /// Changes parameters and covariances according to multiple scattering and energy loss corrections:
  /// if correctForMCS=kTRUE:  compute parameters using Branson correction and add correction resolution to covariances
  /// if correctForMCS=kFALSE: add parameter dispersion due to MCS in parameter covariances
  /// if correctForEnergyLoss=kTRUE:  correct parameters for energy loss and add energy loss fluctuation to covariances
  /// if correctForEnergyLoss=kFALSE: do nothing about energy loss
  
  if (trackParam->GetZ() == zVtx) return; // nothing to be done if already at vertex
  
  if (trackParam->GetZ() > zVtx) { // spectro. (z<0)
    cout<<"E-AliMUONTrackExtrap::ExtrapToVertex: Starting Z ("<<trackParam->GetZ()
        <<") upstream the vertex (zVtx = "<<zVtx<<")"<<endl;
    return;
  }
  
  // Check the vertex position relatively to the absorber
  if (zVtx < AliMUONConstants::AbsZBeg() && zVtx > AliMUONConstants::AbsZEnd()) { // spectro. (z<0)
    cout<<"W-AliMUONTrackExtrap::ExtrapToVertex: Ending Z ("<<zVtx
        <<") inside the front absorber ("<<AliMUONConstants::AbsZBeg()<<","<<AliMUONConstants::AbsZEnd()<<")"<<endl;
  } else if (zVtx < AliMUONConstants::AbsZEnd() ) { // spectro. (z<0)
    cout<<"W-AliMUONTrackExtrap::ExtrapToVertex: Ending Z ("<<zVtx
        <<") downstream the front absorber (zAbsorberEnd = "<<AliMUONConstants::AbsZEnd()<<")"<<endl;
    if (trackParam->CovariancesExist()) ExtrapToZCov(trackParam,zVtx);
    else ExtrapToZ(trackParam,zVtx);
    return;
  }
  
  // Check the track position relatively to the absorber and extrapolate track parameters to the end of the absorber if needed
  if (trackParam->GetZ() > AliMUONConstants::AbsZBeg()) { // spectro. (z<0)
    cout<<"W-AliMUONTrackExtrap::ExtrapToVertex: Starting Z ("<<trackParam->GetZ()
        <<") upstream the front absorber (zAbsorberBegin = "<<AliMUONConstants::AbsZBeg()<<")"<<endl;
    if (trackParam->CovariancesExist()) ExtrapToZCov(trackParam,zVtx);
    else ExtrapToZ(trackParam,zVtx);
    return;
  } else if (trackParam->GetZ() > AliMUONConstants::AbsZEnd()) { // spectro. (z<0)
    cout<<"W-AliMUONTrackExtrap::ExtrapToVertex: Starting Z ("<<trackParam->GetZ()
        <<") inside the front absorber ("<<AliMUONConstants::AbsZBeg()<<","<<AliMUONConstants::AbsZEnd()<<")"<<endl;
  } else {
    if (trackParam->CovariancesExist()) ExtrapToZCov(trackParam,AliMUONConstants::AbsZEnd());
    else ExtrapToZ(trackParam,AliMUONConstants::AbsZEnd());
  }
  
  // Get absorber correction parameters assuming linear propagation in absorber
  Double_t trackXYZOut[3];
  trackXYZOut[0] = trackParam->GetNonBendingCoor();
  trackXYZOut[1] = trackParam->GetBendingCoor();
  trackXYZOut[2] = trackParam->GetZ();
  Double_t trackXYZIn[3];
  if (correctForMCS) { // assume linear propagation until the vertex
    trackXYZIn[2] = TMath::Min(zVtx, AliMUONConstants::AbsZBeg()); // spectro. (z<0)
    trackXYZIn[0] = trackXYZOut[0] + (xVtx - trackXYZOut[0]) / (zVtx - trackXYZOut[2]) * (trackXYZIn[2] - trackXYZOut[2]);
    trackXYZIn[1] = trackXYZOut[1] + (yVtx - trackXYZOut[1]) / (zVtx - trackXYZOut[2]) * (trackXYZIn[2] - trackXYZOut[2]);
  } else {
    AliMUONTrackParam trackParamIn(*trackParam);
    ExtrapToZ(&trackParamIn, TMath::Min(zVtx, AliMUONConstants::AbsZBeg()));
    trackXYZIn[0] = trackParamIn.GetNonBendingCoor();
    trackXYZIn[1] = trackParamIn.GetBendingCoor();
    trackXYZIn[2] = trackParamIn.GetZ();
  }
  Double_t pTot = trackParam->P();
  Double_t pathLength, f0, f1, f2, meanRho, totalELoss, sigmaELoss2;
  if (!GetAbsorberCorrectionParam(trackXYZIn,trackXYZOut,pTot,pathLength,f0,f1,f2,meanRho,totalELoss,sigmaELoss2)) {
    cout<<"E-AliMUONTrackExtrap::ExtrapToVertex: Unable to take into account the absorber effects"<<endl;
    if (trackParam->CovariancesExist()) ExtrapToZCov(trackParam,zVtx);
    else ExtrapToZ(trackParam,zVtx);
    return;
  }
  
  // Compute track parameters and covariances at vertex according to correctForMCS and correctForEnergyLoss flags
  if (correctForMCS) {
    
    if (correctForEnergyLoss) {
      
      // Correct for multiple scattering and energy loss
      CorrectELossEffectInAbsorber(trackParam, 0.5*totalELoss, 0.5*sigmaELoss2);
      CorrectMCSEffectInAbsorber(trackParam, xVtx, yVtx, zVtx, errXVtx, errYVtx,
				 trackXYZIn[2], pathLength, f0, f1, f2);
      CorrectELossEffectInAbsorber(trackParam, 0.5*totalELoss, 0.5*sigmaELoss2);
      
    } else {
      
      // Correct for multiple scattering
      CorrectMCSEffectInAbsorber(trackParam, xVtx, yVtx, zVtx, errXVtx, errYVtx,
				 trackXYZIn[2], pathLength, f0, f1, f2);
    }
    
  } else {
    
    if (correctForEnergyLoss) {
      
      // Correct for energy loss add multiple scattering dispersion in covariance matrix
      CorrectELossEffectInAbsorber(trackParam, 0.5*totalELoss, 0.5*sigmaELoss2);
      AddMCSEffectInAbsorber(trackParam, -pathLength, f0, f1, f2); // (spectro. (z<0))
      ExtrapToZCov(trackParam, trackXYZIn[2]);
      CorrectELossEffectInAbsorber(trackParam, 0.5*totalELoss, 0.5*sigmaELoss2);
      ExtrapToZCov(trackParam, zVtx);
      
    } else {
      
      // add multiple scattering dispersion in covariance matrix
      AddMCSEffectInAbsorber(trackParam, -pathLength, f0, f1, f2); // (spectro. (z<0))
      ExtrapToZCov(trackParam, zVtx);
      
    }
    
  }
  
}

//__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToVertex(AliMUONTrackParam* trackParam,
					Double_t xVtx, Double_t yVtx, Double_t zVtx,
					Double_t errXVtx, Double_t errYVtx)
{
  /// Extrapolate track parameters to vertex, corrected for multiple scattering and energy loss effects
  /// Add branson correction resolution and energy loss fluctuation to parameter covariances
  ExtrapToVertex(trackParam, xVtx, yVtx, zVtx, errXVtx, errYVtx, kTRUE, kTRUE);
}

//__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToVertexWithoutELoss(AliMUONTrackParam* trackParam,
						    Double_t xVtx, Double_t yVtx, Double_t zVtx,
						    Double_t errXVtx, Double_t errYVtx)
{
  /// Extrapolate track parameters to vertex, corrected for multiple scattering effects only
  /// Add branson correction resolution to parameter covariances
  ExtrapToVertex(trackParam, xVtx, yVtx, zVtx, errXVtx, errYVtx, kTRUE, kFALSE);
}

//__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(AliMUONTrackParam* trackParam, Double_t zVtx)
{
  /// Extrapolate track parameters to vertex, corrected for energy loss effects only
  /// Add dispersion due to multiple scattering and energy loss fluctuation to parameter covariances
  ExtrapToVertex(trackParam, 0., 0., zVtx, 0., 0., kFALSE, kTRUE);
}

//__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToVertexUncorrected(AliMUONTrackParam* trackParam, Double_t zVtx)
{
  /// Extrapolate track parameters to vertex without multiple scattering and energy loss corrections
  /// Add dispersion due to multiple scattering to parameter covariances
  ExtrapToVertex(trackParam, 0., 0., zVtx, 0., 0., kFALSE, kFALSE);
}

//__________________________________________________________________________
Double_t AliMUONTrackExtrap::TotalMomentumEnergyLoss(AliMUONTrackParam* trackParam, Double_t xVtx, Double_t yVtx, Double_t zVtx)
{
  /// Calculate the total momentum energy loss in-between the track position and the vertex assuming a linear propagation
  
  if (trackParam->GetZ() == zVtx) return 0.; // nothing to be done if already at vertex
  
  // Check whether the geometry is available
  if (!gGeoManager) {
    cout<<"E-AliMUONTrackExtrap::TotalMomentumEnergyLoss: no TGeo"<<endl;
    return 0.;
  }
  
  // Get encountered material correction parameters assuming linear propagation from vertex to the track position
  Double_t trackXYZOut[3];
  trackXYZOut[0] = trackParam->GetNonBendingCoor();
  trackXYZOut[1] = trackParam->GetBendingCoor();
  trackXYZOut[2] = trackParam->GetZ();
  Double_t trackXYZIn[3];
  trackXYZIn[0] = xVtx;
  trackXYZIn[1] = yVtx;
  trackXYZIn[2] = zVtx;
  Double_t pTot = trackParam->P();
  Double_t pathLength, f0, f1, f2, meanRho, totalELoss, sigmaELoss2;
  GetAbsorberCorrectionParam(trackXYZIn,trackXYZOut,pTot,pathLength,f0,f1,f2,meanRho,totalELoss,sigmaELoss2);
  
  // total momentum corrected for energy loss
  Double_t muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); // GeV
  Double_t e = TMath::Sqrt(pTot*pTot + muMass*muMass);
  Double_t eCorr = e + totalELoss;
  Double_t pTotCorr = TMath::Sqrt(eCorr*eCorr - muMass*muMass);
  
  return pTotCorr - pTot;
}

//__________________________________________________________________________
Double_t AliMUONTrackExtrap::BetheBloch(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicZ, Double_t atomicZoverA)
{
  /// Returns the mean total momentum energy loss of muon with total momentum='pTotal'
  /// in the absorber layer of lenght='pathLength', density='rho', A='atomicA' and Z='atomicZ'
  Double_t muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); // GeV
  
  // mean exitation energy (GeV)
  Double_t i;
  if (atomicZ < 13) i = (12. * atomicZ + 7.) * 1.e-9;
  else i = (9.76 * atomicZ + 58.8 * TMath::Power(atomicZ,-0.19)) * 1.e-9;
  
  return pathLength * rho * AliExternalTrackParam::BetheBlochGeant(pTotal/muMass, rho, 0.20, 3.00, i, atomicZoverA);
}

//__________________________________________________________________________
Double_t AliMUONTrackExtrap::EnergyLossFluctuation2(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicZoverA)
{
  /// Returns the total momentum energy loss fluctuation of muon with total momentum='pTotal'
  /// in the absorber layer of lenght='pathLength', density='rho', A='atomicA' and Z='atomicZ'
  Double_t muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass(); // GeV
  //Double_t eMass = 0.510998918e-3; // GeV
  Double_t k = 0.307075e-3; // GeV.g^-1.cm^2
  Double_t p2=pTotal*pTotal;
  Double_t beta2=p2/(p2 + muMass*muMass);
  
  Double_t fwhm = 2. * k * rho * pathLength * atomicZoverA / beta2; // FWHM of the energy loss Landau distribution
  Double_t sigma2 = fwhm * fwhm / (8.*log(2.)); // gaussian: fwmh = 2 * srqt(2*ln(2)) * sigma (i.e. fwmh = 2.35 * sigma)
  
  //sigma2 = k * rho * pathLength * atomicZ / atomicA * eMass; // sigma2 of the energy loss gaussian distribution
  
  return sigma2;
}

//__________________________________________________________________________
void AliMUONTrackExtrap::Cov2CovP(const TMatrixD &param, TMatrixD &cov)
{
  /// change coordinate system: (X, SlopeX, Y, SlopeY, q/Pyz) -> (X, SlopeX, Y, SlopeY, q*PTot)
  /// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
  
  // charge * total momentum
  Double_t qPTot = TMath::Sqrt(1. + param(1,0)*param(1,0) + param(3,0)*param(3,0)) /
                   TMath::Sqrt(1. + param(3,0)*param(3,0)) / param(4,0);
  
  // Jacobian of the opposite transformation
  TMatrixD jacob(5,5);
  jacob.UnitMatrix();
  jacob(4,1) = qPTot * param(1,0) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,3) = - qPTot * param(1,0) * param(1,0) * param(3,0) /
                 (1. + param(3,0)*param(3,0)) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,4) = - qPTot / param(4,0);
  
  // compute covariances in new coordinate system
  TMatrixD tmp(cov,TMatrixD::kMultTranspose,jacob);
  cov.Mult(jacob,tmp);
}

//__________________________________________________________________________
void AliMUONTrackExtrap::CovP2Cov(const TMatrixD &param, TMatrixD &covP)
{
  /// change coordinate system: (X, SlopeX, Y, SlopeY, q*PTot) -> (X, SlopeX, Y, SlopeY, q/Pyz)
  /// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
  
  // charge * total momentum
  Double_t qPTot = TMath::Sqrt(1. + param(1,0)*param(1,0) + param(3,0)*param(3,0)) /
                   TMath::Sqrt(1. + param(3,0)*param(3,0)) / param(4,0);
  
  // Jacobian of the transformation
  TMatrixD jacob(5,5);
  jacob.UnitMatrix();
  jacob(4,1) = param(4,0) * param(1,0) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,3) = - param(4,0) * param(1,0) * param(1,0) * param(3,0) /
                 (1. + param(3,0)*param(3,0)) / (1. + param(1,0)*param(1,0) + param(3,0)*param(3,0));
  jacob(4,4) = - param(4,0) / qPTot;
  
  // compute covariances in new coordinate system
  TMatrixD tmp(covP,TMatrixD::kMultTranspose,jacob);
  covP.Mult(jacob,tmp);
}

 //__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapOneStepHelix(Double_t charge, Double_t step, const Double_t *vect, Double_t *vout)
{
/// <pre>
///    ******************************************************************
///    *                                                                *
///    *  Performs the tracking of one step in a magnetic field         *
///    *  The trajectory is assumed to be a helix in a constant field   *
///    *  taken at the mid point of the step.                           *
///    *  Parameters:                                                   *
///    *   input                                                        *
///    *     STEP =arc length of the step asked                         *
///    *     VECT =input vector (position,direction cos and momentum)   *
///    *     CHARGE=  electric charge of the particle                   *
///    *   output                                                       *
///    *     VOUT = same as VECT after completion of the step           *
///    *                                                                *
///    *    ==>Called by : USER, GUSWIM                               *
///    *       Author    m.hansroul  *********                          *
///    *       modified  s.egli, s.v.levonian                           *
///    *       modified  v.perevoztchikov
///    *                                                                *
///    ******************************************************************
/// </pre>

// modif: everything in double precision

    Double_t xyz[3], h[4], hxp[3];
    Double_t h2xy, hp, rho, tet;
    Double_t sint, sintt, tsint, cos1t;
    Double_t f1, f2, f3, f4, f5, f6;

    const Int_t kix  = 0;
    const Int_t kiy  = 1;
    const Int_t kiz  = 2;
    const Int_t kipx = 3;
    const Int_t kipy = 4;
    const Int_t kipz = 5;
    const Int_t kipp = 6;

    const Double_t kec = 2.9979251e-4;
    //
    //    ------------------------------------------------------------------
    //
    //       units are kgauss,centimeters,gev/c
    //
    vout[kipp] = vect[kipp];
    if (TMath::Abs(charge) < 0.00001) {
      for (Int_t i = 0; i < 3; i++) {
	vout[i] = vect[i] + step * vect[i+3];
	vout[i+3] = vect[i+3];
      }
      return;
    }
    xyz[0]    = vect[kix] + 0.5 * step * vect[kipx];
    xyz[1]    = vect[kiy] + 0.5 * step * vect[kipy];
    xyz[2]    = vect[kiz] + 0.5 * step * vect[kipz];

    //cmodif: call gufld (xyz, h) changed into:
    TGeoGlobalMagField::Instance()->Field(xyz,h);
 
    h2xy = h[0]*h[0] + h[1]*h[1];
    h[3] = h[2]*h[2]+ h2xy;
    if (h[3] < 1.e-12) {
      for (Int_t i = 0; i < 3; i++) {
	vout[i] = vect[i] + step * vect[i+3];
	vout[i+3] = vect[i+3];
      }
      return;
    }
    if (h2xy < 1.e-12*h[3]) {
      ExtrapOneStepHelix3(charge*h[2], step, vect, vout);
      return;
    }
    h[3] = TMath::Sqrt(h[3]);
    h[0] /= h[3];
    h[1] /= h[3];
    h[2] /= h[3];
    h[3] *= kec;

    hxp[0] = h[1]*vect[kipz] - h[2]*vect[kipy];
    hxp[1] = h[2]*vect[kipx] - h[0]*vect[kipz];
    hxp[2] = h[0]*vect[kipy] - h[1]*vect[kipx];
 
    hp = h[0]*vect[kipx] + h[1]*vect[kipy] + h[2]*vect[kipz];

    rho = -charge*h[3]/vect[kipp];
    tet = rho * step;

    if (TMath::Abs(tet) > 0.15) {
      sint = TMath::Sin(tet);
      sintt = (sint/tet);
      tsint = (tet-sint)/tet;
      cos1t = 2.*(TMath::Sin(0.5*tet))*(TMath::Sin(0.5*tet))/tet;
    } else {
      tsint = tet*tet/36.;
      sintt = (1. - tsint);
      sint = tet*sintt;
      cos1t = 0.5*tet;
    }

    f1 = step * sintt;
    f2 = step * cos1t;
    f3 = step * tsint * hp;
    f4 = -tet*cos1t;
    f5 = sint;
    f6 = tet * cos1t * hp;
 
    vout[kix] = vect[kix] + f1*vect[kipx] + f2*hxp[0] + f3*h[0];
    vout[kiy] = vect[kiy] + f1*vect[kipy] + f2*hxp[1] + f3*h[1];
    vout[kiz] = vect[kiz] + f1*vect[kipz] + f2*hxp[2] + f3*h[2];
 
    vout[kipx] = vect[kipx] + f4*vect[kipx] + f5*hxp[0] + f6*h[0];
    vout[kipy] = vect[kipy] + f4*vect[kipy] + f5*hxp[1] + f6*h[1];
    vout[kipz] = vect[kipz] + f4*vect[kipz] + f5*hxp[2] + f6*h[2];
 
    return;
}

 //__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapOneStepHelix3(Double_t field, Double_t step, const Double_t *vect, Double_t *vout)
{
/// <pre>
///	******************************************************************
///	*								 *
///	*	Tracking routine in a constant field oriented		 *
///	*	along axis 3						 *
///	*	Tracking is performed with a conventional		 *
///	*	helix step method					 *
///	*								 *
///	*    ==>Called by : USER, GUSWIM				 *
///	*	Authors    R.Brun, M.Hansroul  *********		 *
///	*	Rewritten  V.Perevoztchikov
///	*								 *
///	******************************************************************
/// </pre>

    Double_t hxp[3];
    Double_t h4, hp, rho, tet;
    Double_t sint, sintt, tsint, cos1t;
    Double_t f1, f2, f3, f4, f5, f6;

    const Int_t kix  = 0;
    const Int_t kiy  = 1;
    const Int_t kiz  = 2;
    const Int_t kipx = 3;
    const Int_t kipy = 4;
    const Int_t kipz = 5;
    const Int_t kipp = 6;

    const Double_t kec = 2.9979251e-4;

// 
//     ------------------------------------------------------------------
// 
//       units are kgauss,centimeters,gev/c
// 
    vout[kipp] = vect[kipp];
    h4 = field * kec;

    hxp[0] = - vect[kipy];
    hxp[1] = + vect[kipx];
 
    hp = vect[kipz];

    rho = -h4/vect[kipp];
    tet = rho * step;
    if (TMath::Abs(tet) > 0.15) {
      sint = TMath::Sin(tet);
      sintt = (sint/tet);
      tsint = (tet-sint)/tet;
      cos1t = 2.* TMath::Sin(0.5*tet) * TMath::Sin(0.5*tet)/tet;
    } else {
      tsint = tet*tet/36.;
      sintt = (1. - tsint);
      sint = tet*sintt;
      cos1t = 0.5*tet;
    }

    f1 = step * sintt;
    f2 = step * cos1t;
    f3 = step * tsint * hp;
    f4 = -tet*cos1t;
    f5 = sint;
    f6 = tet * cos1t * hp;
 
    vout[kix] = vect[kix] + f1*vect[kipx] + f2*hxp[0];
    vout[kiy] = vect[kiy] + f1*vect[kipy] + f2*hxp[1];
    vout[kiz] = vect[kiz] + f1*vect[kipz] + f3;
 
    vout[kipx] = vect[kipx] + f4*vect[kipx] + f5*hxp[0];
    vout[kipy] = vect[kipy] + f4*vect[kipy] + f5*hxp[1];
    vout[kipz] = vect[kipz] + f4*vect[kipz] + f6;

    return;
}

 //__________________________________________________________________________
Bool_t AliMUONTrackExtrap::ExtrapOneStepRungekutta(Double_t charge, Double_t step, const Double_t* vect, Double_t* vout)
{
/// <pre>
///	******************************************************************
///	*								 *
///	*  Runge-Kutta method for tracking a particle through a magnetic *
///	*  field. Uses Nystroem algorithm (See Handbook Nat. Bur. of	 *
///	*  Standards, procedure 25.5.20)				 *
///	*								 *
///	*  Input parameters						 *
///	*	CHARGE    Particle charge				 *
///	*	STEP	  Step size					 *
///	*	VECT	  Initial co-ords,direction cosines,momentum	 *
///	*  Output parameters						 *
///	*	VOUT	  Output co-ords,direction cosines,momentum	 *
///	*  User routine called  					 *
///	*	CALL GUFLD(X,F) 					 *
///	*								 *
///	*    ==>Called by : USER, GUSWIM				 *
///	*	Authors    R.Brun, M.Hansroul  *********		 *
///	*		   V.Perevoztchikov (CUT STEP implementation)	 *
///	*								 *
///	*								 *
///	******************************************************************
/// </pre>

    Double_t h2, h4, f[4];
    Double_t xyzt[3], a, b, c, ph,ph2;
    Double_t secxs[4],secys[4],seczs[4],hxp[3];
    Double_t g1, g2, g3, g4, g5, g6, ang2, dxt, dyt, dzt;
    Double_t est, at, bt, ct, cba;
    Double_t f1, f2, f3, f4, rho, tet, hnorm, hp, rho1, sint, cost;
    
    Double_t x;
    Double_t y;
    Double_t z;
    
    Double_t xt;
    Double_t yt;
    Double_t zt;

    Double_t maxit = 1992;
    Double_t maxcut = 11;

    const Double_t kdlt   = 1e-4;
    const Double_t kdlt32 = kdlt/32.;
    const Double_t kthird = 1./3.;
    const Double_t khalf  = 0.5;
    const Double_t kec = 2.9979251e-4;

    const Double_t kpisqua = 9.86960440109;
    const Int_t kix  = 0;
    const Int_t kiy  = 1;
    const Int_t kiz  = 2;
    const Int_t kipx = 3;
    const Int_t kipy = 4;
    const Int_t kipz = 5;
  
    // *.
    // *.    ------------------------------------------------------------------
    // *.
    // *             this constant is for units cm,gev/c and kgauss
    // *
    Int_t iter = 0;
    Int_t ncut = 0;
    for(Int_t j = 0; j < 7; j++)
      vout[j] = vect[j];

    Double_t  pinv   = kec * charge / vect[6];
    Double_t tl = 0.;
    Double_t h = step;
    Double_t rest;

 
    do {
      rest  = step - tl;
      if (TMath::Abs(h) > TMath::Abs(rest)) h = rest;
      //cmodif: call gufld(vout,f) changed into:
      TGeoGlobalMagField::Instance()->Field(vout,f);

      // *
      // *             start of integration
      // *
      x      = vout[0];
      y      = vout[1];
      z      = vout[2];
      a      = vout[3];
      b      = vout[4];
      c      = vout[5];

      h2     = khalf * h;
      h4     = khalf * h2;
      ph     = pinv * h;
      ph2    = khalf * ph;
      secxs[0] = (b * f[2] - c * f[1]) * ph2;
      secys[0] = (c * f[0] - a * f[2]) * ph2;
      seczs[0] = (a * f[1] - b * f[0]) * ph2;
      ang2 = (secxs[0]*secxs[0] + secys[0]*secys[0] + seczs[0]*seczs[0]);
      if (ang2 > kpisqua) break;

      dxt    = h2 * a + h4 * secxs[0];
      dyt    = h2 * b + h4 * secys[0];
      dzt    = h2 * c + h4 * seczs[0];
      xt     = x + dxt;
      yt     = y + dyt;
      zt     = z + dzt;
      // *
      // *              second intermediate point
      // *

      est = TMath::Abs(dxt) + TMath::Abs(dyt) + TMath::Abs(dzt);
      if (est > h) {
	if (ncut++ > maxcut) break;
	h *= khalf;
	continue;
      }
 
      xyzt[0] = xt;
      xyzt[1] = yt;
      xyzt[2] = zt;

      //cmodif: call gufld(xyzt,f) changed into:
      TGeoGlobalMagField::Instance()->Field(xyzt,f);

      at     = a + secxs[0];
      bt     = b + secys[0];
      ct     = c + seczs[0];

      secxs[1] = (bt * f[2] - ct * f[1]) * ph2;
      secys[1] = (ct * f[0] - at * f[2]) * ph2;
      seczs[1] = (at * f[1] - bt * f[0]) * ph2;
      at     = a + secxs[1];
      bt     = b + secys[1];
      ct     = c + seczs[1];
      secxs[2] = (bt * f[2] - ct * f[1]) * ph2;
      secys[2] = (ct * f[0] - at * f[2]) * ph2;
      seczs[2] = (at * f[1] - bt * f[0]) * ph2;
      dxt    = h * (a + secxs[2]);
      dyt    = h * (b + secys[2]);
      dzt    = h * (c + seczs[2]);
      xt     = x + dxt;
      yt     = y + dyt;
      zt     = z + dzt;
      at     = a + 2.*secxs[2];
      bt     = b + 2.*secys[2];
      ct     = c + 2.*seczs[2];

      est = TMath::Abs(dxt)+TMath::Abs(dyt)+TMath::Abs(dzt);
      if (est > 2.*TMath::Abs(h)) {
	if (ncut++ > maxcut) break;
	h *= khalf;
	continue;
      }
 
      xyzt[0] = xt;
      xyzt[1] = yt;
      xyzt[2] = zt;

      //cmodif: call gufld(xyzt,f) changed into:
      TGeoGlobalMagField::Instance()->Field(xyzt,f);

      z      = z + (c + (seczs[0] + seczs[1] + seczs[2]) * kthird) * h;
      y      = y + (b + (secys[0] + secys[1] + secys[2]) * kthird) * h;
      x      = x + (a + (secxs[0] + secxs[1] + secxs[2]) * kthird) * h;

      secxs[3] = (bt*f[2] - ct*f[1])* ph2;
      secys[3] = (ct*f[0] - at*f[2])* ph2;
      seczs[3] = (at*f[1] - bt*f[0])* ph2;
      a      = a+(secxs[0]+secxs[3]+2. * (secxs[1]+secxs[2])) * kthird;
      b      = b+(secys[0]+secys[3]+2. * (secys[1]+secys[2])) * kthird;
      c      = c+(seczs[0]+seczs[3]+2. * (seczs[1]+seczs[2])) * kthird;

      est    = TMath::Abs(secxs[0]+secxs[3] - (secxs[1]+secxs[2]))
	+ TMath::Abs(secys[0]+secys[3] - (secys[1]+secys[2]))
	+ TMath::Abs(seczs[0]+seczs[3] - (seczs[1]+seczs[2]));

      if (est > kdlt && TMath::Abs(h) > 1.e-4) {
	if (ncut++ > maxcut) break;
	h *= khalf;
	continue;
      }

      ncut = 0;
      // *               if too many iterations, go to helix
      if (iter++ > maxit) break;

      tl += h;
      if (est < kdlt32) 
	h *= 2.;
      cba    = 1./ TMath::Sqrt(a*a + b*b + c*c);
      vout[0] = x;
      vout[1] = y;
      vout[2] = z;
      vout[3] = cba*a;
      vout[4] = cba*b;
      vout[5] = cba*c;
      rest = step - tl;
      if (step < 0.) rest = -rest;
      if (rest < 1.e-5*TMath::Abs(step)) return kTRUE;

    } while(1);

    // angle too big, use helix
    cout<<"W-AliMUONTrackExtrap::ExtrapOneStepRungekutta: Ruge-Kutta failed: switch to helix"<<endl;

    f1  = f[0];
    f2  = f[1];
    f3  = f[2];
    f4  = TMath::Sqrt(f1*f1+f2*f2+f3*f3);
    if (f4 < 1.e-10) {
      cout<<"E-AliMUONTrackExtrap::ExtrapOneStepRungekutta: magnetic field at (";
      cout<<xyzt[0]<<", "<<xyzt[1]<<", "<<xyzt[2]<<") = "<<f4<<": giving up"<<endl;
      return kFALSE;
    }
    rho = -f4*pinv;
    tet = rho * step;
 
    hnorm = 1./f4;
    f1 = f1*hnorm;
    f2 = f2*hnorm;
    f3 = f3*hnorm;

    hxp[0] = f2*vect[kipz] - f3*vect[kipy];
    hxp[1] = f3*vect[kipx] - f1*vect[kipz];
    hxp[2] = f1*vect[kipy] - f2*vect[kipx];
 
    hp = f1*vect[kipx] + f2*vect[kipy] + f3*vect[kipz];

    rho1 = 1./rho;
    sint = TMath::Sin(tet);
    cost = 2.*TMath::Sin(khalf*tet)*TMath::Sin(khalf*tet);

    g1 = sint*rho1;
    g2 = cost*rho1;
    g3 = (tet-sint) * hp*rho1;
    g4 = -cost;
    g5 = sint;
    g6 = cost * hp;
 
    vout[kix] = vect[kix] + g1*vect[kipx] + g2*hxp[0] + g3*f1;
    vout[kiy] = vect[kiy] + g1*vect[kipy] + g2*hxp[1] + g3*f2;
    vout[kiz] = vect[kiz] + g1*vect[kipz] + g2*hxp[2] + g3*f3;
 
    vout[kipx] = vect[kipx] + g4*vect[kipx] + g5*hxp[0] + g6*f1;
    vout[kipy] = vect[kipy] + g4*vect[kipy] + g5*hxp[1] + g6*f2;
    vout[kipz] = vect[kipz] + g4*vect[kipz] + g5*hxp[2] + g6*f3;

    return kTRUE;
}

