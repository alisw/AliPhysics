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

#include <TMath.h>

Int_t   AliMUONConstants::fgNCh = 14;
Int_t   AliMUONConstants::fgNTrackingCh = 10;
Int_t   AliMUONConstants::fgNTriggerCh = 4;
Int_t   AliMUONConstants::fgNTriggerCircuit = 234;
Int_t   AliMUONConstants::fgNDetElem = 228;



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
// These are the diameter (innner and ounner) values of the active surface
Float_t  AliMUONConstants::fgDmin[7]  = {   36.4,  46.2,  63.0,   79.0,   79.0,  99.0,  100.0};  
Float_t  AliMUONConstants::fgDmax[7]  = {  176.6, 229.0, 308.84, 418.2,  522.0, 850.0, 900.0}; 
 
Int_t    AliMUONConstants::fgMaxZoom = 20;

ClassImp(AliMUONConstants)

//______________________________________________________________________________
Int_t AliMUONConstants::ChamberNumber(Float_t z) 
{
  // return chamber number according z position of hit. Should be taken from geometry ?
 
  Float_t dMaxChamber = DzSlat() + DzCh() + 0.25; // cm st 3 &4 & 5
  if ( z >  (DefaultChamberZ(4)+50.)) dMaxChamber = 7.; // cm stations 1 & 2
  Int_t iChamber;

  for (iChamber = 0; iChamber < 10; iChamber++) {
    
    if (TMath::Abs(z-DefaultChamberZ(iChamber)) < dMaxChamber) {
      return iChamber;
    }
  }
  return -1;
}
