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

//#include <Riostream.h>
#include "AliMUON.h"
#include "AliMUONTrackParam.h" 
#include "AliMUONConstants.h"
#include "AliESDMuonTrack.h"
#include "AliMagF.h" 
#include "AliLog.h" 
#include "AliTracker.h"
#include "AliMUONHitForRec.h"

ClassImp(AliMUONTrackParam) // Class implementation in ROOT context

  //_________________________________________________________________________
AliMUONTrackParam::AliMUONTrackParam()
  : TObject(),
    fInverseBendingMomentum(0.),
    fBendingSlope(0.),
    fNonBendingSlope(0.),
    fZ(0.),
    fBendingCoor(0.),
    fNonBendingCoor(0.),
    fkField(0x0),
    fHitForRecPtr(0x0)
{
// Constructor
  // get field from outside
  fkField = AliTracker::GetFieldMap();
  if (!fkField) AliWarning("No field available");
}

  //_________________________________________________________________________
AliMUONTrackParam::AliMUONTrackParam(const AliMUONTrackParam& theMUONTrackParam)
  : TObject(theMUONTrackParam),
    fInverseBendingMomentum(theMUONTrackParam.fInverseBendingMomentum), 
    fBendingSlope(theMUONTrackParam.fBendingSlope),
    fNonBendingSlope(theMUONTrackParam.fNonBendingSlope),
    fZ(theMUONTrackParam.fZ),
    fBendingCoor(theMUONTrackParam.fBendingCoor),
    fNonBendingCoor(theMUONTrackParam.fNonBendingCoor),
    fkField(theMUONTrackParam.fkField),
    fHitForRecPtr(theMUONTrackParam.fHitForRecPtr)
{
  // Copy constructor

}

  //_________________________________________________________________________
AliMUONTrackParam& AliMUONTrackParam::operator=(const AliMUONTrackParam& theMUONTrackParam)
{
  // Asignment operator
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
  fkField                 =  theMUONTrackParam.fkField;
  fHitForRecPtr           =  theMUONTrackParam.fHitForRecPtr;

  return *this;
}

  //__________________________________________________________________________
AliMUONTrackParam::~AliMUONTrackParam()
{
/// Destructor
/// Update the number of TrackHit's connected to the attached HitForRec if any
  if (fHitForRecPtr) fHitForRecPtr->SetNTrackHits(fHitForRecPtr->GetNTrackHits() - 1); // decrement NTrackHits of hit
}

  //__________________________________________________________________________
void AliMUONTrackParam::SetTrackParam(AliMUONTrackParam& theMUONTrackParam)
{
  /// Set track parameters from "TrackParam" leaving pointer to fHitForRecPtr unchanged
  fInverseBendingMomentum =  theMUONTrackParam.fInverseBendingMomentum; 
  fBendingSlope           =  theMUONTrackParam.fBendingSlope; 
  fNonBendingSlope        =  theMUONTrackParam.fNonBendingSlope; 
  fZ                      =  theMUONTrackParam.fZ; 
  fBendingCoor            =  theMUONTrackParam.fBendingCoor; 
  fNonBendingCoor         =  theMUONTrackParam.fNonBendingCoor;
  
}

  //__________________________________________________________________________
AliMUONHitForRec* AliMUONTrackParam::GetHitForRecPtr(void) const
{
/// return pointer to HitForRec attached to the current TrackParam
/// this method should not be called when fHitForRecPtr == NULL
  if (!fHitForRecPtr) AliWarning("AliMUONTrackParam::GetHitForRecPtr: fHitForRecPtr == NULL");
  return fHitForRecPtr;
}

  //__________________________________________________________________________
Int_t AliMUONTrackParam::Compare(const TObject* TrackParam) const
{
/// "Compare" function to sort with decreasing Z (spectro. muon Z <0).
/// Returns 1 (0, -1) if Z of current TrackHit
/// is smaller than (equal to, larger than) Z of TrackHit
  if (fHitForRecPtr->GetZ() < ((AliMUONTrackParam*)TrackParam)->fHitForRecPtr->GetZ()) return(1);
  else if (fHitForRecPtr->GetZ() == ((AliMUONTrackParam*)TrackParam)->fHitForRecPtr->GetZ()) return(0);
  else return(-1);
}

  //_________________________________________________________________________
void AliMUONTrackParam::GetParamFrom(const AliESDMuonTrack& esdMuonTrack)
{
  // assigned value form ESD track.
  fInverseBendingMomentum =  esdMuonTrack.GetInverseBendingMomentum();
  fBendingSlope           =  TMath::Tan(esdMuonTrack.GetThetaY());
  fNonBendingSlope        =  TMath::Tan(esdMuonTrack.GetThetaX());
  fZ                      =  esdMuonTrack.GetZ(); 
  fBendingCoor            =  esdMuonTrack.GetBendingCoor(); 
  fNonBendingCoor         =  esdMuonTrack.GetNonBendingCoor();
}

  //_________________________________________________________________________
void AliMUONTrackParam::SetParamFor(AliESDMuonTrack& esdMuonTrack)
{
  // assigned value form ESD track.
  esdMuonTrack.SetInverseBendingMomentum(fInverseBendingMomentum);
  esdMuonTrack.SetThetaX(TMath::ATan(fNonBendingSlope));
  esdMuonTrack.SetThetaY(TMath::ATan(fBendingSlope));
  esdMuonTrack.SetZ(fZ); 
  esdMuonTrack.SetBendingCoor(fBendingCoor); 
  esdMuonTrack.SetNonBendingCoor(fNonBendingCoor);
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
    //ExtrapOneStepRungekutta(chargeExtrap, stepLength, vGeant3, vGeant3New);
    ExtrapOneStepHelix(chargeExtrap, stepLength, vGeant3, vGeant3New);
    if ((-forwardBackward * (vGeant3New[2] - Z)) > 0.0) break; // one is beyond Z spectro. z<0
    // better use TArray ????
    for (iGeant3 = 0; iGeant3 < 7; iGeant3++)
      {vGeant3[iGeant3] = vGeant3New[iGeant3];}
  }
  // check maxStepNumber ????
  // Interpolation back to exact Z (2nd order)
  // should be in function ???? using TArray ????
  Double_t dZ12 = vGeant3New[2] - vGeant3[2]; // 1->2
  if (TMath::Abs(dZ12) > 0) {
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
  } else {
    AliWarning(Form("Extrap. to Z not reached, Z = %f",Z));    
  }
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
  // range of Station to be checked ????
  z1 = AliMUONConstants::DefaultChamberZ(2 * Station);
  z2 = AliMUONConstants::DefaultChamberZ(2 * Station + 1);
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
void AliMUONTrackParam::ExtrapToVertex(Double_t xVtx, Double_t yVtx, Double_t zVtx)
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
  BransonCorrection(xVtx,yVtx,zVtx);
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

void AliMUONTrackParam::BransonCorrection(Double_t xVtx,Double_t yVtx,Double_t zVtx)
{
  // Branson correction of track parameters
  // the entry parameters have to be calculated at the end of the absorber
  // simplified version: the z positions of Branson's planes are no longer calculated
  // but are given as inputs. One can use the macros MUONTestAbso.C and DrawTestAbso.C
  // to test this correction. 
  // Would it be possible to calculate all that from Geant configuration ????
  // and to get the Branson parameters from a function in ABSO module ????
  // with an eventual contribution from other detectors like START ????
  //change to take into account the vertex postition (real, reconstruct,....)

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

  Float_t zSmear = zBP ;
  
   pZ = pTotal * (zSmear-zVtx) / TMath::Sqrt((xBP-xVtx) * (xBP-xVtx) + (yBP-yVtx) * (yBP-yVtx) +( zSmear-zVtx) * (zSmear-zVtx) );
   pX = pZ * (xBP - xVtx)/ (zSmear-zVtx);
   pY = pZ * (yBP - yVtx) / (zSmear-zVtx);
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

  fBendingCoor = xVtx;
  fNonBendingCoor = yVtx;
  fZ= zVtx;

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
  fkField->Field(x, b);
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
Double_t AliMUONTrackParam::Px() const
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
Double_t AliMUONTrackParam::Py() const
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
Double_t AliMUONTrackParam::Pz() const
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
Double_t AliMUONTrackParam::P() const
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
 //__________________________________________________________________________
void AliMUONTrackParam::ExtrapOneStepHelix(Double_t charge, Double_t step, 
					 Double_t *vect, Double_t *vout) const
{
//    ******************************************************************
//    *                                                                *
//    *  Performs the tracking of one step in a magnetic field         *
//    *  The trajectory is assumed to be a helix in a constant field   *
//    *  taken at the mid point of the step.                           *
//    *  Parameters:                                                   *
//    *   input                                                        *
//    *     STEP =arc length of the step asked                         *
//    *     VECT =input vector (position,direction cos and momentum)   *
//    *     CHARGE=  electric charge of the particle                   *
//    *   output                                                       *
//    *     VOUT = same as VECT after completion of the step           *
//    *                                                                *
//    *    ==>Called by : <USER>, GUSWIM                               *
//    *       Author    m.hansroul  *********                          *
//    *       modified  s.egli, s.v.levonian                           *
//    *       modified  v.perevoztchikov
//    *                                                                *
//    ******************************************************************
//

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
void AliMUONTrackParam::ExtrapOneStepHelix3(Double_t field, Double_t step, 
					       Double_t *vect, Double_t *vout) const
{
// 
//     ******************************************************************
//     *                                                                *
//     *       Tracking routine in a constant field oriented            *
//     *       along axis 3                                             *
//     *       Tracking is performed with a conventional                *
//     *       helix step method                                        *
//     *                                                                *
//     *    ==>Called by : <USER>, GUSWIM                               *
//     *       Authors    R.Brun, M.Hansroul  *********                 *
//     *       Rewritten  V.Perevoztchikov
//     *                                                                *
//     ******************************************************************
// 

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
void AliMUONTrackParam::ExtrapOneStepRungekutta(Double_t charge, Double_t step, 
						     Double_t* vect, Double_t* vout) const
{
// 
//     ******************************************************************
//     *                                                                *
//     *  Runge-Kutta method for tracking a particle through a magnetic *
//     *  field. Uses Nystroem algorithm (See Handbook Nat. Bur. of     *
//     *  Standards, procedure 25.5.20)                                 *
//     *                                                                *
//     *  Input parameters                                              *
//     *       CHARGE    Particle charge                                *
//     *       STEP      Step size                                      *
//     *       VECT      Initial co-ords,direction cosines,momentum     *
//     *  Output parameters                                             *
//     *       VOUT      Output co-ords,direction cosines,momentum      *
//     *  User routine called                                           *
//     *       CALL GUFLD(X,F)                                          *
//     *                                                                *
//     *    ==>Called by : <USER>, GUSWIM                               *
//     *       Authors    R.Brun, M.Hansroul  *********                 *
//     *                  V.Perevoztchikov (CUT STEP implementation)    *
//     *                                                                *
//     *                                                                *
//     ******************************************************************
// 

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
 void  AliMUONTrackParam::GetField(Double_t *Position, Double_t *Field) const
{
    // interface for arguments in double precision (Why ? ChF)

    Float_t x[3], b[3];

    x[0] = Position[0]; x[1] = Position[1]; x[2] = Position[2];

    fkField->Field(x, b);
    Field[0] = b[0]; Field[1] = b[1]; Field[2] = b[2];

    return;
  }
//_____________________________________________-
void AliMUONTrackParam::Print(Option_t* opt) const
{
//
  // Printing TrackParam information 
  // "full" option for printing all the information about the TrackParam
  //
  TString sopt(opt);
  sopt.ToUpper();
 
  if ( sopt.Contains("FULL") ) { 
    cout << "<AliMUONTrackParam> Bending P=" << setw(5) << setprecision(3)  << 1./GetInverseBendingMomentum() << 
      ", NonBendSlope=" << setw(5) << setprecision(3)  << GetNonBendingSlope()*180./TMath::Pi() <<
      ", BendSlope=" << setw(5) << setprecision(3)     << GetBendingSlope()*180./TMath::Pi()  << 
      ", (x,y,z)_IP=(" <<  setw(5) << setprecision(3) << GetNonBendingCoor() <<
      "," <<  setw(5) << setprecision(3) << GetBendingCoor() <<
      "," <<  setw(5) << setprecision(3) << GetZ() <<
      ") cm, (px,py,pz)=(" << setw(5) << setprecision(3) << Px() <<
      "," << setw(5) << setprecision(3) << Py() <<
      "," << setw(5) << setprecision(3) << Pz() << ") GeV/c" << endl;
  }
  else {
    cout << "<AliMUONTrackParam>"  << endl;
  }
    
}
