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
// Track parameters
// in
// ALICE
// dimuon
// spectrometer
//
///////////////////////////////////////////////////

#include <Riostream.h>

#include "AliCallf77.h" 
#include "AliMUON.h"
#include "AliMUONTrackParam.h" 
#include "AliMUONChamber.h"
#include "AliRun.h" 
#include "AliMagF.h" 
#include "AliLog.h" 

ClassImp(AliMUONTrackParam) // Class implementation in ROOT context

  // A few calls in Fortran or from Fortran (extrap.F).
  // Needed, instead of calls to Geant subroutines,
  // because double precision is necessary for track fit converging with Minuit.
  // The "extrap" functions should be translated into C++ ????
#ifndef WIN32 
# define extrap_onestep_helix extrap_onestep_helix_
# define extrap_onestep_helix3 extrap_onestep_helix3_
# define extrap_onestep_rungekutta extrap_onestep_rungekutta_
# define gufld_double gufld_double_
#else 
# define extrap_onestep_helix EXTRAP_ONESTEP_HELIX
# define extrap_onestep_helix3 EXTRAP_ONESTEP_HELIX3
# define extrap_onestep_rungekutta EXTRAP_ONESTEP_RUNGEKUTTA
# define gufld_double GUFLD_DOUBLE
#endif 

extern "C" {
  void type_of_call extrap_onestep_helix
  (Double_t &Charge, Double_t &StepLength, Double_t *VGeant3, Double_t *VGeant3New);

  void type_of_call extrap_onestep_helix3
  (Double_t &Field, Double_t &StepLength, Double_t *VGeant3, Double_t *VGeant3New);

  void type_of_call extrap_onestep_rungekutta
  (Double_t &Charge, Double_t &StepLength, Double_t *VGeant3, Double_t *VGeant3New);

  void type_of_call gufld_double(Double_t *Position, Double_t *Field) {
    // interface to "gAlice->Field()->Field" for arguments in double precision
    Float_t x[3], b[3];
    x[0] = Position[0]; x[1] = Position[1]; x[2] = Position[2];
    gAlice->Field()->Field(x, b);
    Field[0] = b[0]; Field[1] = b[1]; Field[2] = b[2];
  }
}

  //_________________________________________________________________________
AliMUONTrackParam::AliMUONTrackParam()
  : TObject()
{
// Constructor

  fInverseBendingMomentum = 0;
  fBendingSlope = 0;
  fNonBendingSlope = 0;
  fZ = 0;
  fBendingCoor = 0;
  fNonBendingCoor = 0;
}

  //_________________________________________________________________________
AliMUONTrackParam& 
AliMUONTrackParam::operator=(const AliMUONTrackParam& theMUONTrackParam)
{
  if (this == &theMUONTrackParam)
    return *this;

  // base class assignement
  TObject::operator=(theMUONTrackParam);

  fInverseBendingMomentum =  theMUONTrackParam.fInverseBendingMomentum; 
  fBendingSlope           =  theMUONTrackParam.fBendingSlope; 
  fNonBendingSlope        =  theMUONTrackParam.fNonBendingSlope; 
  fZ                      =  theMUONTrackParam.fZ; 
  fBendingCoor            =  theMUONTrackParam.fBendingCoor; 
  fNonBendingCoor         =  theMUONTrackParam.fNonBendingCoor;

  return *this;
}
  //_________________________________________________________________________
AliMUONTrackParam::AliMUONTrackParam(const AliMUONTrackParam& theMUONTrackParam)
  : TObject(theMUONTrackParam)
{
  fInverseBendingMomentum =  theMUONTrackParam.fInverseBendingMomentum; 
  fBendingSlope           =  theMUONTrackParam.fBendingSlope; 
  fNonBendingSlope        =  theMUONTrackParam.fNonBendingSlope; 
  fZ                      =  theMUONTrackParam.fZ; 
  fBendingCoor            =  theMUONTrackParam.fBendingCoor; 
  fNonBendingCoor         =  theMUONTrackParam.fNonBendingCoor;
}

  //__________________________________________________________________________
void AliMUONTrackParam::ExtrapToZ(Double_t Z)
{
  // Track parameter extrapolation to the plane at "Z".
  // On return, the track parameters resulting from the extrapolation
  // replace the current track parameters.
  if (this->fZ == Z) return; // nothing to be done if same Z
  Double_t forwardBackward; // +1 if forward, -1 if backward
  if (Z < this->fZ) forwardBackward = 1.0; // spectro. z<0 
  else forwardBackward = -1.0;
  Double_t vGeant3[7], vGeant3New[7]; // 7 in parameter ????
  Int_t iGeant3, stepNumber;
  Int_t maxStepNumber = 5000; // in parameter ????
  // For safety: return kTRUE or kFALSE ????
  // Parameter vector for calling EXTRAP_ONESTEP
  SetGeant3Parameters(vGeant3, forwardBackward);
  // sign of charge (sign of fInverseBendingMomentum if forward motion)
  // must be changed if backward extrapolation
  Double_t chargeExtrap = forwardBackward *
    TMath::Sign(Double_t(1.0), this->fInverseBendingMomentum);
  Double_t stepLength = 6.0; // in parameter ????
  // Extrapolation loop
  stepNumber = 0;
  while (((-forwardBackward * (vGeant3[2] - Z)) <= 0.0) &&  // spectro. z<0
	 (stepNumber < maxStepNumber)) {
    stepNumber++;
    // Option for switching between helix and Runge-Kutta ???? 
    // extrap_onestep_rungekutta(chargeExtrap, stepLength, vGeant3, vGeant3New);
    extrap_onestep_helix(chargeExtrap, stepLength, vGeant3, vGeant3New);
    if ((-forwardBackward * (vGeant3New[2] - Z)) > 0.0) break; // one is beyond Z spectro. z<0
    // better use TArray ????
    for (iGeant3 = 0; iGeant3 < 7; iGeant3++)
      {vGeant3[iGeant3] = vGeant3New[iGeant3];}
  }
  // check maxStepNumber ????
  // Interpolation back to exact Z (2nd order)
  // should be in function ???? using TArray ????
  Double_t dZ12 = vGeant3New[2] - vGeant3[2]; // 1->2
  Double_t dZ1i = Z - vGeant3[2]; // 1-i
  Double_t dZi2 = vGeant3New[2] - Z; // i->2
  Double_t xPrime = (vGeant3New[0] - vGeant3[0]) / dZ12;
  Double_t xSecond =
    ((vGeant3New[3] / vGeant3New[5]) - (vGeant3[3] / vGeant3[5])) / dZ12;
  Double_t yPrime = (vGeant3New[1] - vGeant3[1]) / dZ12;
  Double_t ySecond =
    ((vGeant3New[4] / vGeant3New[5]) - (vGeant3[4] / vGeant3[5])) / dZ12;
  vGeant3[0] = vGeant3[0] + xPrime * dZ1i - 0.5 * xSecond * dZ1i * dZi2; // X
  vGeant3[1] = vGeant3[1] + yPrime * dZ1i - 0.5 * ySecond * dZ1i * dZi2; // Y
  vGeant3[2] = Z; // Z
  Double_t xPrimeI = xPrime - 0.5 * xSecond * (dZi2 - dZ1i);
  Double_t yPrimeI = yPrime - 0.5 * ySecond * (dZi2 - dZ1i);
  // (PX, PY, PZ)/PTOT assuming forward motion
  vGeant3[5] =
    1.0 / TMath::Sqrt(1.0 + xPrimeI * xPrimeI + yPrimeI * yPrimeI); // PZ/PTOT
  vGeant3[3] = xPrimeI * vGeant3[5]; // PX/PTOT
  vGeant3[4] = yPrimeI * vGeant3[5]; // PY/PTOT
  // Track parameters from Geant3 parameters,
  // with charge back for forward motion
  GetFromGeant3Parameters(vGeant3, chargeExtrap * forwardBackward);
}

  //__________________________________________________________________________
void AliMUONTrackParam::SetGeant3Parameters(Double_t *VGeant3, Double_t ForwardBackward)
{
  // Set vector of Geant3 parameters pointed to by "VGeant3"
  // from track parameters in current AliMUONTrackParam.
  // Since AliMUONTrackParam is only geometry, one uses "ForwardBackward"
  // to know whether the particle is going forward (+1) or backward (-1).
  VGeant3[0] = this->fNonBendingCoor; // X
  VGeant3[1] = this->fBendingCoor; // Y
  VGeant3[2] = this->fZ; // Z
  Double_t pYZ = TMath::Abs(1.0 / this->fInverseBendingMomentum);
  Double_t pZ =
    pYZ / TMath::Sqrt(1.0 + this->fBendingSlope * this->fBendingSlope);
  VGeant3[6] =
    TMath::Sqrt(pYZ * pYZ +
		pZ * pZ * this->fNonBendingSlope * this->fNonBendingSlope); // PTOT
  VGeant3[5] = -ForwardBackward * pZ / VGeant3[6]; // PZ/PTOT spectro. z<0
  VGeant3[3] = this->fNonBendingSlope * VGeant3[5]; // PX/PTOT
  VGeant3[4] = this->fBendingSlope * VGeant3[5]; // PY/PTOT
}

  //__________________________________________________________________________
void AliMUONTrackParam::GetFromGeant3Parameters(Double_t *VGeant3, Double_t Charge)
{
  // Get track parameters in current AliMUONTrackParam
  // from Geant3 parameters pointed to by "VGeant3",
  // assumed to be calculated for forward motion in Z.
  // "InverseBendingMomentum" is signed with "Charge".
  this->fNonBendingCoor = VGeant3[0]; // X
  this->fBendingCoor = VGeant3[1]; // Y
  this->fZ = VGeant3[2]; // Z
  Double_t pYZ = VGeant3[6] * TMath::Sqrt(1.0 - VGeant3[3] * VGeant3[3]);
  this->fInverseBendingMomentum = Charge / pYZ;
  this->fBendingSlope = VGeant3[4] / VGeant3[5];
  this->fNonBendingSlope = VGeant3[3] / VGeant3[5];
}

  //__________________________________________________________________________
void AliMUONTrackParam::ExtrapToStation(Int_t Station, AliMUONTrackParam *TrackParam)
{
  // Track parameters extrapolated from current track parameters ("this")
  // to both chambers of the station(0..) "Station"
  // are returned in the array (dimension 2) of track parameters
  // pointed to by "TrackParam" (index 0 and 1 for first and second chambers).
  Double_t extZ[2], z1, z2;
  Int_t i1 = -1, i2 = -1; // = -1 to avoid compilation warnings
  AliMUON *pMUON = (AliMUON*) gAlice->GetModule("MUON"); // necessary ????
  // range of Station to be checked ????
  z1 = (&(pMUON->Chamber(2 * Station)))->Z(); // Z of first chamber
  z2 = (&(pMUON->Chamber(2 * Station + 1)))->Z(); // Z of second chamber
  // First and second Z to extrapolate at
  if ((z1 > this->fZ) && (z2 > this->fZ)) {i1 = 0; i2 = 1;}
  else if ((z1 < this->fZ) && (z2 < this->fZ)) {i1 = 1; i2 = 0;}
  else {
  	AliError(Form("Starting Z (%f) in between z1 (%f) and z2 (%f) of station(0..)%d",this->fZ,z1,z2,Station));
//     cout << "ERROR in AliMUONTrackParam::CreateExtrapSegmentInStation" << endl;
//     cout << "Starting Z (" << this->fZ << ") in between z1 (" << z1 <<
//       ") and z2 (" << z2 << ") of station(0..) " << Station << endl;
  }
  extZ[i1] = z1;
  extZ[i2] = z2;
  // copy of track parameters
  TrackParam[i1] = *this;
  // first extrapolation
  (&(TrackParam[i1]))->ExtrapToZ(extZ[0]);
  TrackParam[i2] = TrackParam[i1];
  // second extrapolation
  (&(TrackParam[i2]))->ExtrapToZ(extZ[1]);
  return;
}

  //__________________________________________________________________________
void AliMUONTrackParam::ExtrapToVertex()
{
  // Extrapolation to the vertex.
  // Returns the track parameters resulting from the extrapolation,
  // in the current TrackParam.
  // Changes parameters according to Branson correction through the absorber 
  
  Double_t zAbsorber = -503.0; // to be coherent with the Geant absorber geometry !!!!
                               // spectro. (z<0) 
  // Extrapolates track parameters upstream to the "Z" end of the front absorber
  ExtrapToZ(zAbsorber); // !!!
  // Makes Branson correction (multiple scattering + energy loss)
  BransonCorrection();
  // Makes a simple magnetic field correction through the absorber
  FieldCorrection(zAbsorber);
}


//  Keep this version for future developments
  //__________________________________________________________________________
// void AliMUONTrackParam::BransonCorrection()
// {
//   // Branson correction of track parameters
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

//   pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
//   sign = 1;      
//   if (fInverseBendingMomentum < 0) sign = -1;     
//   pZ = pYZ / (TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope)); 
//   pX = pZ * fNonBendingSlope; 
//   pY = pZ * fBendingSlope; 
//   pTotal = TMath::Sqrt(pYZ *pYZ + pX * pX);
//   xEndAbsorber = fNonBendingCoor; 
//   yEndAbsorber = fBendingCoor; 
//   radiusEndAbsorber2 = xEndAbsorber * xEndAbsorber + yEndAbsorber * yEndAbsorber;

//   if (radiusEndAbsorber2 > rLimit*rLimit) {
//     zEndAbsorber = z1[9];
//     zBP = zBP1;
//   } else {
//     zEndAbsorber = z2[3];
//     zBP = zBP2;
//   }

//   xBP = xEndAbsorber - (pX / pZ) * (zEndAbsorber - zBP);
//   yBP = yEndAbsorber - (pY / pZ) * (zEndAbsorber - zBP);

//   // new parameters after Branson and energy loss corrections
//   pZ = pTotal * zBP / TMath::Sqrt(xBP * xBP + yBP * yBP + zBP * zBP);
//   pX = pZ * xBP / zBP;
//   pY = pZ * yBP / zBP;
//   fBendingSlope = pY / pZ;
//   fNonBendingSlope = pX / pZ;
  
//   pT = TMath::Sqrt(pX * pX + pY * pY);      
//   theta = TMath::ATan2(pT, pZ); 
//   pTotal =
//     TotalMomentumEnergyLoss(rLimit, pTotal, theta, xEndAbsorber, yEndAbsorber);

//   fInverseBendingMomentum = (sign / pTotal) *
//     TMath::Sqrt(1.0 +
// 		fBendingSlope * fBendingSlope +
// 		fNonBendingSlope * fNonBendingSlope) /
//     TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope);

//   // vertex position at (0,0,0)
//   // should be taken from vertex measurement ???
//   fBendingCoor = 0.0;
//   fNonBendingCoor = 0;
//   fZ= 0;
// }

void AliMUONTrackParam::BransonCorrection()
{
  // Branson correction of track parameters
  // the entry parameters have to be calculated at the end of the absorber
  // simplified version: the z positions of Branson's planes are no longer calculated
  // but are given as inputs. One can use the macros MUONTestAbso.C and DrawTestAbso.C
  // to test this correction. 
  // Would it be possible to calculate all that from Geant configuration ????
  // and to get the Branson parameters from a function in ABSO module ????
  // with an eventual contribution from other detectors like START ????
  Double_t  zBP, xBP, yBP;
  Double_t  pYZ, pX, pY, pZ, pTotal, xEndAbsorber, yEndAbsorber, radiusEndAbsorber2, pT, theta;
  Int_t sign;
  static Bool_t first = kTRUE;
  static Double_t zBP1, zBP2, rLimit, thetaLimit, zEndAbsorber;
  // zBP1 for outer part and zBP2 for inner part (only at the first call)
  if (first) {
    first = kFALSE;
  
    zEndAbsorber = -503;  // spectro (z<0)
    thetaLimit = 3.0 * (TMath::Pi()) / 180.;
    rLimit = TMath::Abs(zEndAbsorber) * TMath::Tan(thetaLimit);
    zBP1 = -450; // values close to those calculated with EvalAbso.C
    zBP2 = -480;
  }

  pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  sign = 1;      
  if (fInverseBendingMomentum < 0) sign = -1;  
  pZ = Pz();
  pX = Px(); 
  pY = Py(); 
  pTotal = TMath::Sqrt(pYZ *pYZ + pX * pX);
  xEndAbsorber = fNonBendingCoor; 
  yEndAbsorber = fBendingCoor; 
  radiusEndAbsorber2 = xEndAbsorber * xEndAbsorber + yEndAbsorber * yEndAbsorber;

  if (radiusEndAbsorber2 > rLimit*rLimit) {
    zBP = zBP1;
  } else {
    zBP = zBP2;
  }

  xBP = xEndAbsorber - (pX / pZ) * (zEndAbsorber - zBP);
  yBP = yEndAbsorber - (pY / pZ) * (zEndAbsorber - zBP);

  // new parameters after Branson and energy loss corrections
//   Float_t zSmear = zBP - gRandom->Gaus(0.,2.);  // !!! possible smearing of Z vertex position
  Float_t zSmear = zBP;
  
  pZ = pTotal * zSmear / TMath::Sqrt(xBP * xBP + yBP * yBP + zSmear * zSmear);
  pX = pZ * xBP / zSmear;
  pY = pZ * yBP / zSmear;
  fBendingSlope = pY / pZ;
  fNonBendingSlope = pX / pZ;

  
  pT = TMath::Sqrt(pX * pX + pY * pY);      
  theta = TMath::ATan2(pT, TMath::Abs(pZ)); 
  pTotal = TotalMomentumEnergyLoss(thetaLimit, pTotal, theta);

  fInverseBendingMomentum = (sign / pTotal) *
    TMath::Sqrt(1.0 +
		fBendingSlope * fBendingSlope +
		fNonBendingSlope * fNonBendingSlope) /
    TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope);

  // vertex position at (0,0,0)
  // should be taken from vertex measurement ???
  fBendingCoor = 0.0;
  fNonBendingCoor = 0;
  fZ= 0;
}

  //__________________________________________________________________________
Double_t AliMUONTrackParam::TotalMomentumEnergyLoss(Double_t thetaLimit, Double_t pTotal, Double_t theta)
{
  // Returns the total momentum corrected from energy loss in the front absorber
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
  } else {
    if (pTotal < 20) {
      deltaP  = 2.1207 + 0.05478 * pTotal - 0.00145079 * pTotal * pTotal;
    } else { 
      deltaP = 2.6069 + 0.0051705 * pTotal;
    }
  }
  pTotalCorrected = pTotal + deltaP / TMath::Cos(theta);
  return pTotalCorrected;
}

  //__________________________________________________________________________
void AliMUONTrackParam::FieldCorrection(Double_t Z)
{
  // 
  // Correction of the effect of the magnetic field in the absorber
  // Assume a constant field along Z axis.

  Float_t b[3],x[3]; 
  Double_t bZ;
  Double_t pYZ,pX,pY,pZ,pT;
  Double_t pXNew,pYNew;
  Double_t c;

  pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  c = TMath::Sign(1.0,fInverseBendingMomentum); // particle charge 
 
  pZ = Pz();
  pX = Px(); 
  pY = Py();
  pT = TMath::Sqrt(pX*pX+pY*pY);

  if (TMath::Abs(pZ) <= 0) return;
  x[2] = Z/2;
  x[0] = x[2]*fNonBendingSlope;  
  x[1] = x[2]*fBendingSlope;

  // Take magn. field value at position x.
  gAlice->Field()->Field(x, b);
  bZ =  b[2];
 
  // Transverse momentum rotation
  // Parameterized with the study of DeltaPhi = phiReco - phiGen as a function of pZ.
  Double_t phiShift = c*0.436*0.0003*bZ*Z/pZ; 
 // Rotate momentum around Z axis.
  pXNew = pX*TMath::Cos(phiShift) - pY*TMath::Sin(phiShift);
  pYNew = pX*TMath::Sin(phiShift) + pY*TMath::Cos(phiShift);
 
  fBendingSlope = pYNew / pZ;
  fNonBendingSlope = pXNew / pZ;
  
  fInverseBendingMomentum = c / TMath::Sqrt(pYNew*pYNew+pZ*pZ);
 
}
  //__________________________________________________________________________
Double_t AliMUONTrackParam::Px()
{
  // return px from track paramaters
  Double_t pYZ, pZ, pX;
  pYZ = 0;
  if (  TMath::Abs(fInverseBendingMomentum) > 0 )
    pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  pZ = -pYZ / (TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope));  // spectro. (z<0)
  pX = pZ * fNonBendingSlope; 
  return pX;
}
  //__________________________________________________________________________
Double_t AliMUONTrackParam::Py()
{
  // return px from track paramaters
  Double_t pYZ, pZ, pY;
  pYZ = 0;
  if (  TMath::Abs(fInverseBendingMomentum) > 0 )
    pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  pZ = -pYZ / (TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope));  // spectro. (z<0)
  pY = pZ * fBendingSlope; 
  return pY;
}
  //__________________________________________________________________________
Double_t AliMUONTrackParam::Pz()
{
  // return px from track paramaters
  Double_t pYZ, pZ;
  pYZ = 0;
  if (  TMath::Abs(fInverseBendingMomentum) > 0 )
    pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  pZ = -pYZ / (TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope));  // spectro. (z<0)
  return pZ;
}
  //__________________________________________________________________________
Double_t AliMUONTrackParam::P()
{
  // return p from track paramaters
  Double_t  pYZ, pZ, p;
  pYZ = 0;
  if (  TMath::Abs(fInverseBendingMomentum) > 0 )
    pYZ = TMath::Abs(1.0 / fInverseBendingMomentum);
  pZ = -pYZ / (TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope));  // spectro. (z<0)
  p = TMath::Abs(pZ) * 
    TMath::Sqrt(1.0 + fBendingSlope * fBendingSlope + fNonBendingSlope * fNonBendingSlope);
  return p;
  
}
