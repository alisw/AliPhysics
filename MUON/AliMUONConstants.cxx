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

#include "AliMUONConstants.h"

#include "AliLog.h"

#include "TMath.h"
#include "TClass.h"
#include "AliMpConstants.h"

//-----------------------------------------------------------------------------
/// \class AliMUONConstants
/// This class holds various constants to be used in many places,
/// such as the number of tracking and trigger chambers, 
/// some geometrical constants (to build the initial geometry for simulation)
/// and mathieson distribution default values.
/// Those constants should as much as possible replace hard-coded values
/// which are to be considered strictly illegal in the MUON code (or any code,
/// by the way).
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONConstants)
/// \endcond

Int_t   AliMUONConstants::fgNTrackingSt = 5;
Int_t   AliMUONConstants::fgNTriggerCircuit = 234;
Int_t   AliMUONConstants::fgNDetElem = 228;
Int_t   AliMUONConstants::fgNGeomModules = 20;
Float_t AliMUONConstants::fgkTriggerTofLimit = 75E-9;

Float_t AliMUONConstants::fgDefaultChamberZ[14] = 
  {-526.16, -545.24, -676.4, -695.4, // St12
   -967.5, -998.5, -1276.5, -1307.5, -1406.6, -1437.6,// updated 08/05, EDMS id 335328 (A. Tournaire)
   -1603.5, -1620.5, -1703.5, -1720.5}; // M1 & M2


// These are used by AliMUONConstants::ChamberNumber and must be calculated once
// by that method from fgDzCh, fgDzSlat, fgDefaultChamberZ and fgSt345inclination,
// so for now we set everything to zero.
Float_t AliMUONConstants::fgDefaultChamberMinZ[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Float_t AliMUONConstants::fgDefaultChamberMaxZ[14] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

Float_t AliMUONConstants::fgDefaultRatioTriggerChamber[4] =
{1., 1.01060, 1.06236, 1.07296};


Float_t  AliMUONConstants::fgSt345inclination = 0.794; // in degrees, rotation axis is X axis 

Float_t  AliMUONConstants::fgDzCh   = 15.5/2.;
Float_t  AliMUONConstants::fgDzSlat = 8.5/2.;

Float_t  AliMUONConstants::fgSqrtKx3 = 0.7131;
Float_t  AliMUONConstants::fgSqrtKy3 = 0.7642;

Float_t  AliMUONConstants::fgSqrtKx3St1 = 0.7000;
Float_t  AliMUONConstants::fgSqrtKy3St1 = 0.7550;

Float_t  AliMUONConstants::fgChargeCorrel    = 0.11;
Float_t  AliMUONConstants::fgChargeCorrelSt1 = 1.0; //??? 
Float_t  AliMUONConstants::fgPitch     = 0.25;
Float_t  AliMUONConstants::fgPitchSt1  = 0.21; 

// From Alain TOURNAIRE    
// ALICE / ALICE Engineering baseline / Dimuonspectrometer (DIS) v7-1
// EDMS Id 335328 for "search in EDMS 
// These are the diameter (Dmin == innner and Dmax - outner) values of the active surface
// In the case of Dmax, the value corresponds to the maximum diameter of the active surface with 2pi coverture in phi
Float_t  AliMUONConstants::fgDmin[7]  = {   36.4,  46.2,  63.0,   79.0,   79.0,  98.8,  100.0};  // cm
Float_t  AliMUONConstants::fgDmax[7]  = {  176.6, 229.0, 308.84, 418.2,  522.0, 850.0, 900.0};   // cm
 
Int_t    AliMUONConstants::fgMaxZoom = 20;

// Defaults parameters for dipole magnet
// From ALICE Dimuon - parameters / geometry table,
// V7-3 (version 7 created 24/03/2004 updated 25/10/2005)
Double_t AliMUONConstants::fgCoilZ = -994.05;
Double_t AliMUONConstants::fgCoilL = 502.1;
Double_t AliMUONConstants::fgYokeZ = -986.6;
Double_t AliMUONConstants::fgYokeL = 309.4;

// Defaults parameters for absorber (27/06/07)
const Double_t AliMUONConstants::fgkAbsZBeg = -90.;
const Double_t AliMUONConstants::fgkAbsZEnd = -505.;
    
// Default trigger chamber resolution (cm)
// Warning: the resolution refers only to ALIGNMENT
// For the total resolution the strip width should be taken into account!
const Double_t AliMUONConstants::fgkTriggerNonBendingReso = 0.2;
const Double_t AliMUONConstants::fgkTriggerBendingReso = 0.2;

// Defaults parameters for muon filter (19/11/07)
const Double_t AliMUONConstants::fgkMuonFilterZBeg = -1471.;
const Double_t AliMUONConstants::fgkMuonFilterZEnd = -1471.-120.;
const Double_t AliMUONConstants::fgkMuonFilterX0 = 1.76;

// Defaults parameters for track reconstruction
Double_t AliMUONConstants::fgChamberThicknessInX0[10] = {0.065, 0.065, 0.075, 0.075, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035};

const Double_t AliMUONConstants::fgkMostProbBendingMomentum = 2.0;

Float_t AliMUONConstants::fgAverageChamberT[14]=
  {17.64*1E-9, 18.28*1E-9, 22.68*1E-9, 23.33*1E-9, 32.42*1E-9, 33.48*1E-9, 42.76*1E-9,
   43.81*1E-9, 47.13*1E-9, 48.17*1E-9, 53.75*1E-9, 54.32*1E-9, 57.12*1E-9, 57.67*1E-9};

//______________________________________________________________________________
Int_t AliMUONConstants::NCh()
{
  return AliMpConstants::NofChambers();
}

//______________________________________________________________________________
Int_t AliMUONConstants::NTrackingCh()
{
  return AliMpConstants::NofTrackingChambers();
}

//______________________________________________________________________________
Int_t AliMUONConstants::ChamberNumber(Float_t z, bool warn)
{
  // return chamber number according z position of hit. Should be taken from geometry ?

  if (fgDefaultChamberMinZ[0] == 0) // Are the min/max Z arrays initialised?
  {
    // The min and max Z arrays need to be calculated.
    for (Int_t i = 0; i < NCh(); i++)
    {
      Float_t a = 0, b = 0;
      if (4 <= i and i < 10)
      {
        Float_t dzAngle = TMath::Tan(TMath::Pi()*St345Inclination()/180.) * Rmax(i/2);
        // We add 2.5mm since Rmax is an under-estimate.
        a = DefaultChamberZ(i) + DzSlat() + DzCh() + dzAngle + 0.25;
        b = DefaultChamberZ(i) - DzSlat() - DzCh() - dzAngle - 0.25;
      }
      else
      {
        a = DefaultChamberZ(i) + DzSlat();
        b = DefaultChamberZ(i) - DzSlat();
      }
      fgDefaultChamberMinZ[i] = TMath::Min(a, b);
      fgDefaultChamberMaxZ[i] = TMath::Max(a, b);
    }
  }

  // We can apply a binary search for the chamber since the fgDefaultChamberMinZ and
  // fgDefaultChamberMaxZ arrays are ordered.
  Int_t mini = 0, maxi = NCh()-1;
  while (mini <= maxi)
  {
    Int_t iChamber = (maxi + mini) / 2;
    if (z < fgDefaultChamberMinZ[iChamber])
      mini = iChamber+1;
    else if (z > fgDefaultChamberMaxZ[iChamber])
      maxi = iChamber-1;
    else
      // We are between min and max Z of chamber number iChamber so we found our chamber.
      return iChamber;
  }

  if (warn) AliWarningClass(Form("No chamber number found for z = %f",z));
  return -1;
}

//______________________________________________________________________________
Float_t AliMUONConstants::ReducedQTot(Float_t qtot, Float_t timeDif)
{
  // return a reduced charge if the hit belongs to a track from a pileup event
  Float_t q = qtot*1.19*(1.24-timeDif*1E6)*TMath::Exp(-(0.97-timeDif*1E6)*(0.97-timeDif*1E6)/2.42);
  return q;
}
