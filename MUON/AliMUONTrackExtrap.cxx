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

///////////////////////////////////////////////////
//
// Tools
// for
// track
// extrapolation
// in
// ALICE
// dimuon
// spectrometer
//
///////////////////////////////////////////////////

#include <Riostream.h>
#include <TMatrixD.h>

#include "AliMUONTrackExtrap.h" 
#include "AliMUONTrackParam.h"
#include "AliMUONConstants.h"
#include "AliMagF.h" 
#include "AliLog.h" 
#include "AliTracker.h"

ClassImp(AliMUONTrackExtrap) // Class implementation in ROOT context

const AliMagF* AliMUONTrackExtrap::fgkField = 0x0;
const Int_t    AliMUONTrackExtrap::fgkMaxStepNumber = 5000;
const Double_t AliMUONTrackExtrap::fgkStepLength = 6.;

  //__________________________________________________________________________
Double_t AliMUONTrackExtrap::GetImpactParamFromBendingMomentum(Double_t bendingMomentum)
{
  /// Returns impact parameter at vertex in bending plane (cm),
  /// from the signed bending momentum "BendingMomentum" in bending plane (GeV/c),
  /// using simple values for dipole magnetic field.
  /// The sign of "BendingMomentum" is the sign of the charge.
  
  if (bendingMomentum == 0.) return 1.e10;
  
  Double_t simpleBPosition = 0.5 * (AliMUONConstants::CoilZ() + AliMUONConstants::YokeZ());
  Double_t simpleBLength = 0.5 * (AliMUONConstants::CoilL() + AliMUONConstants::YokeL());
  Float_t b[3], x[3] = {0.,0.,(Float_t) simpleBPosition};
  if (fgkField) fgkField->Field(x,b);
  else {
    cout<<"F-AliMUONTrackExtrap::GetField: fgkField = 0x0"<<endl;
    exit(-1);
  }
  Double_t simpleBValue = (Double_t) b[0];
  
  return (-0.0003 * simpleBValue * simpleBLength * simpleBPosition / bendingMomentum);
}

  //__________________________________________________________________________
Double_t AliMUONTrackExtrap::GetBendingMomentumFromImpactParam(Double_t impactParam)
{
  /// Returns signed bending momentum in bending plane (GeV/c),
  /// the sign being the sign of the charge for particles moving forward in Z,
  /// from the impact parameter "ImpactParam" at vertex in bending plane (cm),
  /// using simple values for dipole magnetic field.
  
  if (impactParam == 0.) return 1.e10;
  
  Double_t simpleBPosition = 0.5 * (AliMUONConstants::CoilZ() + AliMUONConstants::YokeZ());
  Double_t simpleBLength = 0.5 * (AliMUONConstants::CoilL() + AliMUONConstants::YokeL());
  Float_t b[3], x[3] = {0.,0.,(Float_t) simpleBPosition};
  if (fgkField) fgkField->Field(x,b);
  else {
    cout<<"F-AliMUONTrackExtrap::GetField: fgkField = 0x0"<<endl;
    exit(-1);
  }
  Double_t simpleBValue = (Double_t) b[0];
  
  return (-0.0003 * simpleBValue * simpleBLength * simpleBPosition / impactParam);
}

  //__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToZ(AliMUONTrackParam* trackParam, Double_t zEnd)
{
  /// Track parameter extrapolation to the plane at "Z".
  /// On return, the track parameters resulting from the extrapolation are updated in trackParam.
  if (trackParam->GetZ() == zEnd) return; // nothing to be done if same Z
  Double_t forwardBackward; // +1 if forward, -1 if backward
  if (zEnd < trackParam->GetZ()) forwardBackward = 1.0; // spectro. z<0 
  else forwardBackward = -1.0;
  Double_t v3[7], v3New[7]; // 7 in parameter ????
  Int_t i3, stepNumber;
  // For safety: return kTRUE or kFALSE ????
  // Parameter vector for calling EXTRAP_ONESTEP
  ConvertTrackParamForExtrap(trackParam, v3, forwardBackward);
  // sign of charge (sign of fInverseBendingMomentum if forward motion)
  // must be changed if backward extrapolation
  Double_t chargeExtrap = forwardBackward * TMath::Sign(Double_t(1.0), trackParam->GetInverseBendingMomentum());
  // Extrapolation loop
  stepNumber = 0;
  while (((-forwardBackward * (v3[2] - zEnd)) <= 0.0) && (stepNumber < fgkMaxStepNumber)) { // spectro. z<0
    stepNumber++;
    // Option for switching between helix and Runge-Kutta ???? 
    //ExtrapOneStepRungekutta(chargeExtrap, fgkStepLength, v3, v3New);
    ExtrapOneStepHelix(chargeExtrap, fgkStepLength, v3, v3New);
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
    cout<<"W-AliMUONTrackExtrap::ExtrapToZ: Extrap. to Z not reached, Z = "<<zEnd<<endl;
  }
  // Track parameters from 3 parameters,
  // with charge back for forward motion
  RecoverTrackParam(v3, chargeExtrap * forwardBackward, trackParam);
}

  //__________________________________________________________________________
void AliMUONTrackExtrap::ConvertTrackParamForExtrap(AliMUONTrackParam* trackParam, Double_t *v3, Double_t forwardBackward)
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
  Double_t pYZ = v3[6] * TMath::Sqrt(1.0 - v3[3] * v3[3]);
  trackParam->SetInverseBendingMomentum(charge/pYZ);
  trackParam->SetBendingSlope(v3[4]/v3[5]);
  trackParam->SetNonBendingSlope(v3[3]/v3[5]);
}

  //__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToZCov(AliMUONTrackParam* trackParam, Double_t zEnd)
{
  /// Track parameters and their covariances extrapolated to the plane at "zEnd".
  /// On return, results from the extrapolation are updated in trackParam.
  
  if (trackParam->GetZ() == zEnd) return; // nothing to be done if same z
  
  // Save the actual track parameters
  AliMUONTrackParam trackParamSave(*trackParam);
  Double_t nonBendingCoor 	  = trackParamSave.GetNonBendingCoor();
  Double_t nonBendingSlope 	  = trackParamSave.GetNonBendingSlope();
  Double_t bendingCoor 		  = trackParamSave.GetBendingCoor();
  Double_t bendingSlope 	  = trackParamSave.GetBendingSlope();
  Double_t inverseBendingMomentum = trackParamSave.GetInverseBendingMomentum();
  Double_t zBegin		  = trackParamSave.GetZ();
  
  // Extrapolate track parameters to "zEnd"
  ExtrapToZ(trackParam,zEnd);
  Double_t extrapNonBendingCoor 	= trackParam->GetNonBendingCoor();
  Double_t extrapNonBendingSlope 	= trackParam->GetNonBendingSlope();
  Double_t extrapBendingCoor 		= trackParam->GetBendingCoor();
  Double_t extrapBendingSlope		= trackParam->GetBendingSlope();
  Double_t extrapInverseBendingMomentum = trackParam->GetInverseBendingMomentum();
  
  // Get the pointer to the parameter covariance matrix
  if (!trackParam->CovariancesExist()) {
    //cout<<"W-AliMUONTrackExtrap::ExtrapToZCov: track parameter covariance matrix does not exist"<<endl;
    //cout<<"                                    -> nothing to extrapolate !!"<<endl;
    return;
  }
  TMatrixD* paramCov = trackParam->GetCovariances();
  
  // Calculate the jacobian related to the track parameters extrapolation to "zEnd"
  TMatrixD jacob(5,5);
  jacob = 0.;
  Double_t dParam[5];
  for (Int_t i=0; i<5; i++) {
    // Skip jacobian calculation for parameters with no associated error
    if ((*paramCov)(i,i) == 0.) continue;
    // Small variation of parameter i only
    for (Int_t j=0; j<5; j++) {
      if (j==i) {
        dParam[j] = TMath::Sqrt((*paramCov)(i,i));
	if (j == 4) dParam[j] *= TMath::Sign(1.,-inverseBendingMomentum); // variation always in the same direction
      } else dParam[j] = 0.;
    }
    // Set new parameters
    trackParamSave.SetNonBendingCoor	    (nonBendingCoor	    + dParam[0]);
    trackParamSave.SetNonBendingSlope	    (nonBendingSlope	    + dParam[1]);
    trackParamSave.SetBendingCoor	    (bendingCoor	    + dParam[2]);
    trackParamSave.SetBendingSlope	    (bendingSlope	    + dParam[3]);
    trackParamSave.SetInverseBendingMomentum(inverseBendingMomentum + dParam[4]);
    trackParamSave.SetZ			    (zBegin);
    // Extrapolate new track parameters to "zEnd"
    ExtrapToZ(&trackParamSave,zEnd);
    // Calculate the jacobian
    jacob(0,i) = (trackParamSave.GetNonBendingCoor()	     - extrapNonBendingCoor 	   ) / dParam[i];
    jacob(1,i) = (trackParamSave.GetNonBendingSlope()	     - extrapNonBendingSlope 	   ) / dParam[i];
    jacob(2,i) = (trackParamSave.GetBendingCoor()	     - extrapBendingCoor    	   ) / dParam[i];
    jacob(3,i) = (trackParamSave.GetBendingSlope()	     - extrapBendingSlope   	   ) / dParam[i];
    jacob(4,i) = (trackParamSave.GetInverseBendingMomentum() - extrapInverseBendingMomentum) / dParam[i];
  }
  
  // Extrapolate track parameter covariances to "zEnd"
  TMatrixD tmp((*paramCov),TMatrixD::kMultTranspose,jacob);
  (*paramCov) = TMatrixD(jacob,TMatrixD::kMult,tmp);
  
}

  //__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToStation(AliMUONTrackParam* trackParamIn, Int_t station, AliMUONTrackParam *trackParamOut)
{
  /// Track parameters extrapolated from "trackParamIn" to both chambers of the station(0..) "station"
  /// are returned in the array (dimension 2) of track parameters pointed to by "TrackParamOut"
  /// (index 0 and 1 for first and second chambers).
  Double_t extZ[2], z1, z2;
  Int_t i1 = -1, i2 = -1; // = -1 to avoid compilation warnings
  // range of station to be checked ????
  z1 = AliMUONConstants::DefaultChamberZ(2 * station);
  z2 = AliMUONConstants::DefaultChamberZ(2 * station + 1);
  // First and second Z to extrapolate at
  if ((z1 > trackParamIn->GetZ()) && (z2 > trackParamIn->GetZ())) {i1 = 0; i2 = 1;}
  else if ((z1 < trackParamIn->GetZ()) && (z2 < trackParamIn->GetZ())) {i1 = 1; i2 = 0;}
  else {
    cout<<"E-AliMUONTrackExtrap::ExtrapToStation: Starting Z ("<<trackParamIn->GetZ()
    	<<") in between z1 ("<<z1<<") and z2 ("<<z2<<") of station(0..)"<<station<<endl;
    exit(-1);
  }
  extZ[i1] = z1;
  extZ[i2] = z2;
  // copy of track parameters
  trackParamOut[i1] = *trackParamIn;
  // first extrapolation
  ExtrapToZ(&(trackParamOut[i1]),extZ[0]);
  trackParamOut[i2] = trackParamOut[i1];
  // second extrapolation
  ExtrapToZ(&(trackParamOut[i2]),extZ[1]);
  return;
}

  //__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToVertexUncorrected(AliMUONTrackParam* trackParam, Double_t zVtx)
{
  /// Extrapolation to the vertex (at the z position "zVtx") without Branson and Field correction.
  /// Returns the track parameters resulting from the extrapolation in the current TrackParam.
  /// Include multiple Coulomb scattering effects in trackParam covariances.
  
  if (trackParam->GetZ() > zVtx) { // spectro. (z<0)
    cout<<"W-AliMUONTrackExtrap::ExtrapToVertexUncorrected: Starting Z ("<<trackParam->GetZ()
    	<<") upstream the vertex (zVtx = "<<zVtx<<")"<<endl;
    exit(-1);
  }
  
  if (zVtx < AliMUONConstants::ZAbsorberEnd()) { // spectro. (z<0)
    cout<<"W-AliMUONTrackExtrap::ExtrapToVertexUncorrected: Ending Z ("<<zVtx
    	<<") downstream the front absorber (zAbsorberEnd = "<<AliMUONConstants::ZAbsorberEnd()<<")"<<endl;
    
    ExtrapToZCov(trackParam,zVtx);
    return;
  }
  
  // First Extrapolates track parameters upstream to the "Z" end of the front absorber
  if (trackParam->GetZ() < AliMUONConstants::ZAbsorberEnd()) { // spectro. (z<0)
    ExtrapToZCov(trackParam,AliMUONConstants::ZAbsorberEnd());
  } else {
    cout<<"W-AliMUONTrackExtrap::ExtrapToVertexUncorrected: Starting Z ("<<trackParam->GetZ()
    	<<") upstream or inside the front absorber (zAbsorberEnd = "<<AliMUONConstants::ZAbsorberEnd()<<")"<<endl;
  }
  
  // Then go through all the absorber layers
  Double_t tan3 = TMath::Tan(3./180.*TMath::Pi());
  Double_t r0Norm, x0, z, zElement, dZ, nonBendingCoor, bendingCoor;
  for (Int_t iElement=AliMUONConstants::NAbsorberElements()-1; iElement>=0; iElement--) {
    zElement = AliMUONConstants::ZAbsorberElement(iElement);
    z = trackParam->GetZ();
    if (z > zElement) continue; // spectro. (z<0)
    nonBendingCoor = trackParam->GetNonBendingCoor();
    bendingCoor = trackParam->GetBendingCoor();
    r0Norm = nonBendingCoor * nonBendingCoor + bendingCoor * bendingCoor;
    r0Norm  = TMath::Sqrt(r0Norm) / TMath::Abs(trackParam->GetZ()) / tan3;
    if (r0Norm > 1.) x0 = AliMUONConstants::X0AbsorberOut(iElement); // outer part of the absorber
    else x0 = AliMUONConstants::X0AbsorberIn(iElement); // inner part of the absorber
    
    if (zVtx > zElement) {
      ExtrapToZCov(trackParam,zElement); // extrapolate to absorber element "iElement"
      dZ = zElement - z;
      AddMCSEffectInTrackParamCov(trackParam,dZ,x0); // include MCS effect in covariances
    } else {
      ExtrapToZCov(trackParam,zVtx); // extrapolate to zVtx
      dZ = zVtx - z;
      AddMCSEffectInTrackParamCov(trackParam,dZ,x0); // include MCS effect in covariances
      break;
    }
  }
  
  // finally go to the vertex
  ExtrapToZCov(trackParam,zVtx);
  
}

  //__________________________________________________________________________
void AliMUONTrackExtrap::AddMCSEffectInTrackParamCov(AliMUONTrackParam *param, Double_t dZ, Double_t x0)
{
  /// Add to the track parameter covariances the effects of multiple Coulomb scattering
  /// through a material of thickness "dZ" and of radiation length "x0"
  /// assuming linear propagation and using the small angle approximation.
  
  Double_t bendingSlope = param->GetBendingSlope();
  Double_t nonBendingSlope = param->GetNonBendingSlope();
  Double_t inverseTotalMomentum2 = param->GetInverseBendingMomentum() * param->GetInverseBendingMomentum() *
  				   (1.0 + bendingSlope * bendingSlope) /
			  	   (1.0 + bendingSlope *bendingSlope + nonBendingSlope * nonBendingSlope); 
  // Path length in the material
  Double_t pathLength = TMath::Abs(dZ) * TMath::Sqrt(1.0 + bendingSlope*bendingSlope + nonBendingSlope*nonBendingSlope);
  Double_t pathLength2 = pathLength * pathLength;
  // relativistic velocity
  Double_t velo = 1.;
  // Angular dispersion square of the track (variance) in a plane perpendicular to the trajectory
  Double_t theta02 = 0.0136 / velo * (1 + 0.038 * TMath::Log(pathLength/x0));
  theta02 *= theta02 * inverseTotalMomentum2 * pathLength / x0;
  
  // Add effects of multiple Coulomb scattering in track parameter covariances
  TMatrixD* paramCov = param->GetCovariances();
  Double_t varCoor 	= pathLength2 * theta02 / 3.;
  Double_t varSlop 	= theta02;
  Double_t covCorrSlope = pathLength * theta02 / 2.;
  // Non bending plane
  (*paramCov)(0,0) += varCoor;		(*paramCov)(0,1) += covCorrSlope;
  (*paramCov)(1,0) += covCorrSlope;	(*paramCov)(1,1) += varSlop;
  // Bending plane
  (*paramCov)(2,2) += varCoor;		(*paramCov)(2,3) += covCorrSlope;
  (*paramCov)(3,2) += covCorrSlope;	(*paramCov)(3,3) += varSlop;
  
}

  //__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapToVertex(AliMUONTrackParam* trackParam, Double_t xVtx, Double_t yVtx, Double_t zVtx)
{
  /// Extrapolation to the vertex.
  /// Returns the track parameters resulting from the extrapolation in the current TrackParam.
  /// Changes parameters according to Branson correction through the absorber 
  
  // Extrapolates track parameters upstream to the "Z" end of the front absorber
  ExtrapToZ(trackParam,AliMUONConstants::ZAbsorberEnd()); // !!!
  // Makes Branson correction (multiple scattering + energy loss)
  BransonCorrection(trackParam,xVtx,yVtx,zVtx);
  // Makes a simple magnetic field correction through the absorber
  FieldCorrection(trackParam,AliMUONConstants::ZAbsorberEnd());
}


//  Keep this version for future developments
  //__________________________________________________________________________
// void AliMUONTrackExtrap::BransonCorrection(AliMUONTrackParam* trackParam)
// {
//   /// Branson correction of track parameters
//   // the entry parameters have to be calculated at the end of the absorber
//   Double_t zEndAbsorber, zBP, xBP, yBP;
//   Double_t  pYZ, pX, pY, pZ, pTotal, xEndAbsorber, yEndAbsorber, radiusEndAbsorber2, pT, theta;
//   Int_t sign;
//   // Would it be possible to calculate all that from Geant configuration ????
//   // and to get the Branson parameters from a function in ABSO module ????
//   // with an eventual contribution from other detectors like START ????
//   // Radiation lengths outer part theta > 3 degres
//   static Double_t x01[9] = { 18.8,    // C (cm)
// 			     10.397,   // Concrete (cm)
// 			     0.56,    // Plomb (cm)
// 			     47.26,   // Polyethylene (cm)
// 			     0.56,   // Plomb (cm)
// 			     47.26,   // Polyethylene (cm)
// 			     0.56,   // Plomb (cm)
// 			     47.26,   // Polyethylene (cm)
// 			     0.56 };   // Plomb (cm)
//   // inner part theta < 3 degres
//   static Double_t x02[3] = { 18.8,    // C (cm)
// 			     10.397,   // Concrete (cm)
// 			     0.35 };    // W (cm) 
//   // z positions of the materials inside the absober outer part theta > 3 degres
//   static Double_t z1[10] = { 90, 315, 467, 472, 477, 482, 487, 492, 497, 502 };
//   // inner part theta < 3 degres
//   static Double_t z2[4] = { 90, 315, 467, 503 };
//   static Bool_t first = kTRUE;
//   static Double_t zBP1, zBP2, rLimit;
//   // Calculates z positions of the Branson's planes: zBP1 for outer part and zBP2 for inner part (only at the first call)
//   if (first) {
//     first = kFALSE;
//     Double_t aNBP = 0.0;
//     Double_t aDBP = 0.0;
//     Int_t iBound;
//     
//     for (iBound = 0; iBound < 9; iBound++) {
//       aNBP = aNBP +
// 	(z1[iBound+1] * z1[iBound+1] * z1[iBound+1] -
// 	 z1[iBound]   * z1[iBound]   * z1[iBound]    ) / x01[iBound];
//       aDBP = aDBP +
// 	(z1[iBound+1] * z1[iBound+1] - z1[iBound]   * z1[iBound]    ) / x01[iBound];
//     }
//     zBP1 = (2.0 * aNBP) / (3.0 * aDBP);
//     aNBP = 0.0;
//     aDBP = 0.0;
//     for (iBound = 0; iBound < 3; iBound++) {
//       aNBP = aNBP +
// 	(z2[iBound+1] * z2[iBound+1] * z2[iBound+1] -
// 	 z2[iBound]   * z2[iBound ]  * z2[iBound]    ) / x02[iBound];
//       aDBP = aDBP +
// 	(z2[iBound+1] * z2[iBound+1] - z2[iBound] * z2[iBound]) / x02[iBound];
//     }
//     zBP2 = (2.0 * aNBP) / (3.0 * aDBP);
//     rLimit = z2[3] * TMath::Tan(3.0 * (TMath::Pi()) / 180.);
//   }
// 
//   pYZ = TMath::Abs(1.0 / trackParam->GetInverseBendingMomentum());
//   sign = 1;      
//   if (trackParam->GetInverseBendingMomentum() < 0) sign = -1;     
//   pZ = pYZ / (TMath::Sqrt(1.0 + trackParam->GetBendingSlope() * trackParam->GetBendingSlope())); 
//   pX = pZ * trackParam->GetNonBendingSlope(); 
//   pY = pZ * trackParam->GetBendingSlope(); 
//   pTotal = TMath::Sqrt(pYZ *pYZ + pX * pX);
//   xEndAbsorber = trackParam->GetNonBendingCoor(); 
//   yEndAbsorber = trackParam->GetBendingCoor(); 
//   radiusEndAbsorber2 = xEndAbsorber * xEndAbsorber + yEndAbsorber * yEndAbsorber;
// 
//   if (radiusEndAbsorber2 > rLimit*rLimit) {
//     zEndAbsorber = z1[9];
//     zBP = zBP1;
//   } else {
//     zEndAbsorber = z2[3];
//     zBP = zBP2;
//   }
// 
//   xBP = xEndAbsorber - (pX / pZ) * (zEndAbsorber - zBP);
//   yBP = yEndAbsorber - (pY / pZ) * (zEndAbsorber - zBP);
// 
//   // new parameters after Branson and energy loss corrections
//   pZ = pTotal * zBP / TMath::Sqrt(xBP * xBP + yBP * yBP + zBP * zBP);
//   pX = pZ * xBP / zBP;
//   pY = pZ * yBP / zBP;
//   trackParam->SetBendingSlope(pY/pZ);
//   trackParam->SetNonBendingSlope(pX/pZ);
//   
//   pT = TMath::Sqrt(pX * pX + pY * pY);      
//   theta = TMath::ATan2(pT, pZ); 
//   pTotal = TotalMomentumEnergyLoss(rLimit, pTotal, theta, xEndAbsorber, yEndAbsorber);
// 
//   trackParam->SetInverseBendingMomentum((sign / pTotal) *
//     TMath::Sqrt(1.0 +
// 		trackParam->GetBendingSlope() * trackParam->GetBendingSlope() +
// 		trackParam->GetNonBendingSlope() * trackParam->GetNonBendingSlope()) /
//     TMath::Sqrt(1.0 + trackParam->GetBendingSlope() * trackParam->GetBendingSlope()));
// 
//   // vertex position at (0,0,0)
//   // should be taken from vertex measurement ???
//   trackParam->SetBendingCoor(0.);
//   trackParam->SetNonBendingCoor(0.);
//   trackParam->SetZ(0.);
// }

void AliMUONTrackExtrap::BransonCorrection(AliMUONTrackParam* trackParam, Double_t xVtx, Double_t yVtx, Double_t zVtx)
{
  /// Branson correction of track parameters
  // the entry parameters have to be calculated at the end of the absorber
  // simplified version: the z positions of Branson's planes are no longer calculated
  // but are given as inputs. One can use the macros MUONTestAbso.C and DrawTestAbso.C
  // to test this correction. 
  // Would it be possible to calculate all that from Geant configuration ????
  // and to get the Branson parameters from a function in ABSO module ????
  // with an eventual contribution from other detectors like START ????
  // change to take into account the vertex postition (real, reconstruct,....)

  Double_t  zBP, xBP, yBP;
  Double_t  pYZ, pX, pY, pZ, pTotal, xEndAbsorber, yEndAbsorber, radiusEndAbsorber2, pT, theta;
  Int_t sign;
  static Bool_t first = kTRUE;
  static Double_t zBP1, zBP2, rLimit, thetaLimit;
  // zBP1 for outer part and zBP2 for inner part (only at the first call)
  if (first) {
    first = kFALSE;
  
    thetaLimit = 3.0 * (TMath::Pi()) / 180.;
    rLimit = TMath::Abs(AliMUONConstants::ZAbsorberEnd()) * TMath::Tan(thetaLimit);
    zBP1 = -450; // values close to those calculated with EvalAbso.C
    zBP2 = -480;
  }

  pYZ = TMath::Abs(1.0 / trackParam->GetInverseBendingMomentum());
  sign = 1;      
  if (trackParam->GetInverseBendingMomentum() < 0) sign = -1;  
  pZ = trackParam->Pz();
  pX = trackParam->Px(); 
  pY = trackParam->Py(); 
  pTotal = TMath::Sqrt(pYZ *pYZ + pX * pX);
  xEndAbsorber = trackParam->GetNonBendingCoor(); 
  yEndAbsorber = trackParam->GetBendingCoor(); 
  radiusEndAbsorber2 = xEndAbsorber * xEndAbsorber + yEndAbsorber * yEndAbsorber;

  if (radiusEndAbsorber2 > rLimit*rLimit) {
    zBP = zBP1;
  } else {
    zBP = zBP2;
  }

  xBP = xEndAbsorber - (pX / pZ) * (AliMUONConstants::ZAbsorberEnd() - zBP);
  yBP = yEndAbsorber - (pY / pZ) * (AliMUONConstants::ZAbsorberEnd() - zBP);

  // new parameters after Branson and energy loss corrections
//   Float_t zSmear = zBP - gRandom->Gaus(0.,2.);  // !!! possible smearing of Z vertex position

  Float_t zSmear = zBP ;
  
  pZ = pTotal * (zSmear-zVtx) / TMath::Sqrt((xBP-xVtx) * (xBP-xVtx) + (yBP-yVtx) * (yBP-yVtx) +( zSmear-zVtx) * (zSmear-zVtx) );
  pX = pZ * (xBP - xVtx)/ (zSmear-zVtx);
  pY = pZ * (yBP - yVtx) / (zSmear-zVtx);
  trackParam->SetBendingSlope(pY/pZ);
  trackParam->SetNonBendingSlope(pX/pZ);

  
  pT = TMath::Sqrt(pX * pX + pY * pY);      
  theta = TMath::ATan2(pT, TMath::Abs(pZ)); 
  pTotal = TotalMomentumEnergyLoss(thetaLimit, pTotal, theta);

  trackParam->SetInverseBendingMomentum((sign / pTotal) *
    TMath::Sqrt(1.0 +
		trackParam->GetBendingSlope() * trackParam->GetBendingSlope() +
		trackParam->GetNonBendingSlope() * trackParam->GetNonBendingSlope()) /
    TMath::Sqrt(1.0 + trackParam->GetBendingSlope() * trackParam->GetBendingSlope()));

  // vertex position at (0,0,0)
  // should be taken from vertex measurement ???

  trackParam->SetBendingCoor(xVtx);
  trackParam->SetNonBendingCoor(yVtx);
  trackParam->SetZ(zVtx);

}

  //__________________________________________________________________________
Double_t AliMUONTrackExtrap::TotalMomentumEnergyLoss(Double_t thetaLimit, Double_t pTotal, Double_t theta)
{
  /// Returns the total momentum corrected from energy loss in the front absorber
  // One can use the macros MUONTestAbso.C and DrawTestAbso.C
  // to test this correction. 
  // Momentum energy loss behaviour evaluated with the simulation of single muons (april 2002)
  Double_t deltaP, pTotalCorrected;

   // Parametrization to be redone according to change of absorber material ????
  // See remark in function BransonCorrection !!!!
  // The name is not so good, and there are many arguments !!!!
  if (theta  < thetaLimit ) {
    if (pTotal < 20) {
      deltaP = 2.5938 + 0.0570 * pTotal - 0.001151 * pTotal * pTotal;
    } else {
      deltaP = 3.0714 + 0.011767 *pTotal;
    }
    deltaP *= 0.75; // AZ
  } else {
    if (pTotal < 20) {
      deltaP  = 2.1207 + 0.05478 * pTotal - 0.00145079 * pTotal * pTotal;
    } else { 
      deltaP = 2.6069 + 0.0051705 * pTotal;
    }
    deltaP *= 0.9; // AZ
  }
  pTotalCorrected = pTotal + deltaP / TMath::Cos(theta);
  return pTotalCorrected;
}

  //__________________________________________________________________________
void AliMUONTrackExtrap::FieldCorrection(AliMUONTrackParam *trackParam, Double_t zEnd)
{
  /// Correction of the effect of the magnetic field in the absorber
  // Assume a constant field along Z axis.
  Float_t b[3],x[3]; 
  Double_t bZ;
  Double_t pYZ,pX,pY,pZ,pT;
  Double_t pXNew,pYNew;
  Double_t c;

  pYZ = TMath::Abs(1.0 / trackParam->GetInverseBendingMomentum());
  c = TMath::Sign(1.0,trackParam->GetInverseBendingMomentum()); // particle charge 
 
  pZ = trackParam->Pz();
  pX = trackParam->Px(); 
  pY = trackParam->Py();
  pT = TMath::Sqrt(pX*pX+pY*pY);

  if (TMath::Abs(pZ) <= 0) return;
  x[2] = zEnd/2;
  x[0] = x[2]*trackParam->GetNonBendingSlope();  
  x[1] = x[2]*trackParam->GetBendingSlope();

  // Take magn. field value at position x.
  if (fgkField) fgkField->Field(x,b);
  else {
    cout<<"F-AliMUONTrackExtrap::FieldCorrection: fgkField = 0x0"<<endl;
    exit(-1);
  }
  bZ =  b[2];
 
  // Transverse momentum rotation
  // Parameterized with the study of DeltaPhi = phiReco - phiGen as a function of pZ.
  Double_t phiShift = c*0.436*0.0003*bZ*zEnd/pZ; 
 // Rotate momentum around Z axis.
  pXNew = pX*TMath::Cos(phiShift) - pY*TMath::Sin(phiShift);
  pYNew = pX*TMath::Sin(phiShift) + pY*TMath::Cos(phiShift);
 
  trackParam->SetBendingSlope(pYNew/pZ);
  trackParam->SetNonBendingSlope(pXNew/pZ);
  
  trackParam->SetInverseBendingMomentum(c/TMath::Sqrt(pYNew*pYNew+pZ*pZ));
 
}

 //__________________________________________________________________________
void AliMUONTrackExtrap::ExtrapOneStepHelix(Double_t charge, Double_t step, Double_t *vect, Double_t *vout)
{
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
///    *    ==>Called by : <USER>, GUSWIM                               *
///    *       Author    m.hansroul  *********                          *
///    *       modified  s.egli, s.v.levonian                           *
///    *       modified  v.perevoztchikov
///    *                                                                *
///    ******************************************************************

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
    GetField (xyz, h);
 
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
void AliMUONTrackExtrap::ExtrapOneStepHelix3(Double_t field, Double_t step, Double_t *vect, Double_t *vout)
{
///	******************************************************************
///	*								 *
///	*	Tracking routine in a constant field oriented		 *
///	*	along axis 3						 *
///	*	Tracking is performed with a conventional		 *
///	*	helix step method					 *
///	*								 *
///	*    ==>Called by : <USER>, GUSWIM				 *
///	*	Authors    R.Brun, M.Hansroul  *********		 *
///	*	Rewritten  V.Perevoztchikov
///	*								 *
///	******************************************************************

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
void AliMUONTrackExtrap::ExtrapOneStepRungekutta(Double_t charge, Double_t step, Double_t* vect, Double_t* vout)
{
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
///	*    ==>Called by : <USER>, GUSWIM				 *
///	*	Authors    R.Brun, M.Hansroul  *********		 *
///	*		   V.Perevoztchikov (CUT STEP implementation)	 *
///	*								 *
///	*								 *
///	******************************************************************

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

      GetField(vout,f);

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
      GetField(xyzt,f);

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
      GetField(xyzt,f);

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
      if (rest < 1.e-5*TMath::Abs(step)) return;

    } while(1);

    // angle too big, use helix

    f1  = f[0];
    f2  = f[1];
    f3  = f[2];
    f4  = TMath::Sqrt(f1*f1+f2*f2+f3*f3);
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

    return;
}
//___________________________________________________________
 void  AliMUONTrackExtrap::GetField(Double_t *Position, Double_t *Field)
{
  /// interface for arguments in double precision (Why ? ChF)
  Float_t x[3], b[3];

  x[0] = Position[0]; x[1] = Position[1]; x[2] = Position[2];

  if (fgkField) fgkField->Field(x,b);
  else {
    cout<<"F-AliMUONTrackExtrap::GetField: fgkField = 0x0"<<endl;
    exit(-1);
  }
  
  Field[0] = b[0]; Field[1] = b[1]; Field[2] = b[2];

  return;
}

