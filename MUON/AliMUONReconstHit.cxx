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

#include "AliMUONReconstHit.h"

ClassImp(AliMUONReconstHit)

//_____________________________________________________________________________
AliMUONReconstHit::AliMUONReconstHit() 
  : TObject()
{
// Default constructor

  fCorrelIndex[0]=fCorrelIndex[1]=fCorrelIndex[2]=fCorrelIndex[3]=0;
  fX[0]=fX[1]=fX[2]=fX[3]=0; fY[0]=fY[1]=fY[2]=fY[3]=0; 
}

//_____________________________________________________________________________
AliMUONReconstHit::AliMUONReconstHit(Int_t *idx, Float_t *x, Float_t *y)
  : TObject()
{
    //
    // Creates a MUON correlation object
    //
    for(Int_t i=0; i<4; i++) {
	fCorrelIndex[i]  = idx[i];
	fX[i]    = x[i];
	fY[i]    = y[i];
    }
}
