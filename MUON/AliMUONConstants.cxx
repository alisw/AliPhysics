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

#include <TMath.h>
#include "AliMUONConstants.h"

ClassImp(AliMUONConstants)

Int_t   AliMUONConstants::fgNCh = 14;
Int_t   AliMUONConstants::fgNTrackingCh = 10;
Int_t   AliMUONConstants::fgNTriggerCh = 4;
Int_t   AliMUONConstants::fgNTriggerCircuit = 234;

Float_t AliMUONConstants::fgDefaultChamberZ[14] = 
  {-533.5, -546.5, -678.5, -693.5, // St12
   -966.9, -997.9, -1274.5, -1305.5, -1408.6, -1439.6, // St345  update sept04 Ch. Finck 
   -1603.5, -1620.5, -1703.5, -1720.5}; // M1 & M2

Float_t  AliMUONConstants::fgDzCh   = 15.5/2.;
Float_t  AliMUONConstants::fgDzSlat = 8.5/2.;

Float_t  AliMUONConstants::fgSqrtKx3Slat = 0.7131;
Float_t  AliMUONConstants::fgSqrtKy3Slat = 0.7642;

Float_t  AliMUONConstants::fgSqrtKx3St12 = 0.7000;
Float_t  AliMUONConstants::fgSqrtKy3St12 = 0.7550;

Float_t  AliMUONConstants::fgChargeCorrelSlat = 0.11;
Float_t  AliMUONConstants::fgChargeCorrelSt12 = 0.0; //???

Float_t  AliMUONConstants::fgPitchSlat = 0.25;
Float_t  AliMUONConstants::fgPitchSt12 = 0.20; 

Float_t  AliMUONConstants::fgDmin[7] = {  36.4,  46.2,  66.0,   80.,   80., 100., 100.};    
Float_t  AliMUONConstants::fgDmax[7]  = {183., 245., 395.,  560.,  563., 850., 900.};  
Int_t    AliMUONConstants::fgMaxZoom = 20;

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
