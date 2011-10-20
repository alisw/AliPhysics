/**************************************************************************
 * Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
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

//====================================================================================================================================================
//
//      Digit description for the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliDigit.h"
#include "AliMFTDigit.h"

ClassImp(AliMFTDigit)

//====================================================================================================================================================

AliMFTDigit::AliMFTDigit():
  AliDigit(),
  fNMCTracks(0),
  fPixelX(-1),
  fPixelY(-1),
  fPixelZ(0),
  fPixelCenterX(0),
  fPixelCenterY(0),  
  fPixelCenterZ(0),  
  fPixelWidthX(0),
  fPixelWidthY(0),  
  fPixelWidthZ(0),  
  fPlane(-1),
  fDetElemID(-1),
  fEloss(0),
  fNElectrons(0)
{

  // default cosntructor

  for (Int_t iTr=0; iTr<fNMaxMCTracks; iTr++) fMCLabel[iTr] = -1; 

}

//====================================================================================================================================================
