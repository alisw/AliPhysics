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

ClassImp(AliMUONConstants)

Int_t   AliMUONConstants::fgNCh = 14;
Int_t   AliMUONConstants::fgNTrackingCh = 10;
Int_t   AliMUONConstants::fgNTriggerCh = 4;
Int_t   AliMUONConstants::fgNTriggerCircuit = 234;
Int_t   AliMUONConstants::fgNofDetElements[14] =
{ 4, 4, 4, 4, 18, 18, 26, 26, 26, 26, 18, 18, 18, 18 };

Float_t AliMUONConstants::fgDefaultChamberZ[14] = 
  {-533.5, -546.5, -678.5, -693.5, // St12
   -966.9, -997.9, -1274.5, -1305.5, -1408.6, -1439.6, // St345  update sept04 Ch. Finck 
   -1603.5, -1620.5, -1703.5, -1720.5}; // M1 & M2

Float_t  AliMUONConstants::fgDzCh   = 15.5/2.;
Float_t  AliMUONConstants::fgDzSlat = 8.5/2.;


Float_t  AliMUONConstants::fgDmin[7] = {  36.4,  46.2,  66.0,   80.,   80., 100., 100.};    
Float_t  AliMUONConstants::fgDmax[7]  = {183., 245., 395.,  560.,  563., 850., 900.};  
Int_t   AliMUONConstants::fgMaxZoom = 20;

//______________________________________________________________________________
Int_t AliMUONConstants::GetChamberId(Int_t detElemId)
{ 
// Get chamber Id from detection element Id
// ---

  return detElemId/100 - 1;
}  

//______________________________________________________________________________
Int_t AliMUONConstants::GetFirstDetElemId(Int_t chamberId)
{
// Get first detection element Id for chamber specified by chamber Id
// ---

  return (chamberId+1)*100;
}  
