/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *;
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
////////////////////////////////////////////////
//  MergedHits classes for set ZDC            //
////////////////////////////////////////////////
//

#include "AliZDCMergedHit.h"
#include "AliRun.h"

ClassImp(AliZDCMergedHit)
  
//_____________________________________________________________________________
AliZDCMergedHit::AliZDCMergedHit(Int_t *sector, Float_t *mhits)
{
  //
  // Add a ZDC hit for merging
  //
  Int_t i;
  for(i=0; i<2; i++) {
     fSector[i] = sector[i];
  }
  fPrimKinEn 	= mhits[0];
  fXImpact 	= mhits[1];
  fYImpact 	= mhits[2];
  fSFlag 	= mhits[3];
  fLightPMQ 	= mhits[4];
  fLightPMC 	= mhits[5];
  fEnergy 	= mhits[6]; 
  
}
