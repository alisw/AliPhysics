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

/*
$Log
*/

////////////////////////////////////////////////
//  Hits classes for set ZDC                  //
////////////////////////////////////////////////


#include "AliZDCHit.h"

ClassImp(AliZDCHit)
  
//_____________________________________________________________________________
AliZDCHit::AliZDCHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Add a ZDC hit
  //
  Int_t i;
  for(i=0; i<2; i++) {
     fVolume[i] = vol[i];
  }
  fX 		= hits[0];
  fY 		= hits[1];
  fZ 		= hits[2];
  fPrimKinEn 	= hits[3];
  fXImpact 	= hits[4];
  fYImpact 	= hits[5];
  fSFlag 	= hits[6];
  fLightPMQ 	= hits[7];
  fLightPMC 	= hits[8];
  fEnergy 	= hits[9]; 
  
}
