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

///
/// This class holds various constants to be used in many places,
/// such as the number of tracking and trigger chambers, 
/// some geometrical constants (to build the initial geometry for simulation)
/// and mathieson distribution default values.
/// Those constants should as much as possible replace hard-coded values
/// which are to be considered strictly illegal in the MUON code (or any code,
/// by the way).
///

/// \cond CLASSIMP
ClassImp(AliMUONConstants)
/// \endcond

Int_t   AliMUONConstants::fgNCh = 14;
Int_t   AliMUONConstants::fgNTrackingCh = 10;
Int_t   AliMUONConstants::fgNTriggerCh = 4;
Int_t   AliMUONConstants::fgNTriggerCircuit = 234;
Int_t   AliMUONConstants::fgNDetElem = 228;
Int_t   AliMUONConstants::fgNGeomModules = 20;



Float_t AliMUONConstants::fgDefaultChamberZ[14] = 
  {-526.16, -545.24, -676.4, -695.4, // St12
   -967.5, -998.5, -1276.5, -1307.5, -1406.6, -1437.6,// updated 08/05, EDMS id 335328 (A. Tournaire)
   -1603.5, -1620.5, -1703.5, -1720.5}; // M1 & M2

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
Float_t  AliMUONConstants::fgDmin[7]  = {   36.4,  46.2,  63.0,   79.0,   79.0,  99.0,  100.0};  // cm
Float_t  AliMUONConstants::fgDmax[7]  = {  176.6, 229.0, 308.84, 418.2,  522.0, 850.0, 900.0};   // cm
 
Int_t    AliMUONConstants::fgMaxZoom = 20;

//______________________________________________________________________________
Int_t AliMUONConstants::ChamberNumber(Float_t z) 
{
  // return chamber number according z position of hit. Should be taken from geometry ?
 
  Float_t dMaxChamber = DzSlat() + DzCh() + 0.25; // cm st 3 &4 & 5
  dMaxChamber += 3.00;                            // factor for inclination of chamber  
  // dMaxChamber += Rmax(4) * TMath::Sin(fgSt345inclination*TMath::Pi()/360.); 
                                                  // factor for inclination of chamber 
  if ( z >  (DefaultChamberZ(4)+50.)) dMaxChamber = 7.; // cm stations 1 & 2
  Int_t iChamber;

  for (iChamber = 0; iChamber < 10; iChamber++) {
    if (TMath::Abs(z-DefaultChamberZ(iChamber)) < dMaxChamber) {
      return iChamber;
    }
  }
  
  if ( z > DefaultChamberZ(NTrackingCh()-1) ) {
    AliWarningClass(Form("No chamber number found for z = %f",z));
    // for (iChamber = 0; iChamber < 10; iChamber++) {
    //   cout << iChamber << " zpos: " << DefaultChamberZ(iChamber)
    //        << "  from " << DefaultChamberZ(iChamber) + dMaxChamber
    // 	      << "  to " << DefaultChamberZ(iChamber) - dMaxChamber 
    //        << "  delta " << dMaxChamber << endl;
    //}
  }
  return -1;
}
