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
    fPhysics     = digits[3];
    fHit         = digits[4];
    
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
    fPhysics     = digits[3];
    fHit         = digits[4];

    for(Int_t i=0; i<kMAXTRACKSPERRICHDIGIT; i++) {
	fTcharges[i]  = charges[i];
	fTracks[i]    = tracks[i];
    }
}

void AliRICHDigit::Print(Option_t*)const
{
  Info("","PadX=%3i, PadY=%3i, Signal=%4i, Phys=%4i, Hit=%5i TID1=%5i, TID2=%5i, TID3=%5i TC1=%4i TC2=%4i",
           fPadX,    fPadY,    fSignal,    fPhysics, fHit,   fTracks[0],fTracks[1],fTracks[2],fTcharges[0],fTcharges[1]);
}//void AliRICHdigit::Print(Option_t *option)const
