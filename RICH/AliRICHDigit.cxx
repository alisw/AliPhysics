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
  $Log$
*/


#include "AliRICHDigit.h"

ClassImp(AliRICHDigit)
//_____________________________________________________________________________
AliRICHDigit::AliRICHDigit(Int_t *digits)
{
    //
    // Creates a RICH digit object to be updated
    //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fSignal      = digits[2];
    
}
//_____________________________________________________________________________
AliRICHDigit::AliRICHDigit(Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Creates a RICH digit object
    //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fSignal      = digits[2];
    for(Int_t i=0; i<100; i++) {
	fTcharges[i]  = charges[i];
	fTracks[i]    = tracks[i];
    }
}
