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

#include "AliMUONDigit.h"


ClassImp(AliMUONDigit)
//_____________________________________________________________________________
 AliMUONDigit::AliMUONDigit(const AliMUONDigit& digits):TObject(digits)
{
// copy constructor
  
    fPadX        = digits.fPadX;
    fPadY        = digits.fPadY;
    fCathode     = digits.fCathode;
    fSignal      = digits.fSignal;
    fPhysics     = digits.fPhysics;
    fHit         = digits.fHit;
    fDetElemId   = digits.fDetElemId;

    for(Int_t i=0; i<kMAXTRACKS; i++) {
	fTcharges[i]  = digits.fTcharges[i];
	fTracks[i]    = digits.fTracks[i];
    }
}

//_____________________________________________________________________________
 AliMUONDigit::AliMUONDigit()
   : TObject(),
     fPadX(0),
     fPadY(0),
     fCathode(0),
     fSignal(0),
     fPhysics(0),
     fHit(0),
     fDetElemId(0)
{
// Default constructor

    for(Int_t i=0; i<kMAXTRACKS; i++) {
	fTcharges[i]  = 0;
	fTracks[i]    = 0;
    }
}

//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit(Int_t *digits)
{
  //
  // Creates a MUON digit object to be updated
  //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fCathode     = digits[2];
    fSignal      = digits[3];
    fPhysics     = digits[4];
    fHit         = digits[5];
    fDetElemId   = digits[6];
}
//_____________________________________________________________________________
AliMUONDigit::AliMUONDigit(Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Creates a MUON digit object
    //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fCathode     = digits[2];
    fSignal      = digits[3];
    fPhysics     = digits[4];
    fHit         = digits[5];
    fDetElemId   = digits[6];

    for(Int_t i=0; i<kMAXTRACKS; i++) {
	fTcharges[i]  = charges[i];
	fTracks[i]    = tracks[i];
    }
}

AliMUONDigit::~AliMUONDigit()
{
    // Destructor 
}

//_____________________________________________________________________________
AliMUONDigit& AliMUONDigit::operator=(const AliMUONDigit& digits)
{
    if (this == &digits)
      return *this;
 
    fPadX        = digits.fPadX;
    fPadY        = digits.fPadY;
    fCathode     = digits.fCathode;
    fSignal      = digits.fSignal;
    fPhysics     = digits.fPhysics;
    fHit         = digits.fHit;
    fDetElemId   = digits.fDetElemId;

    for(Int_t i=0; i<kMAXTRACKS; i++) {
	fTcharges[i]  = digits.fTcharges[i];
	fTracks[i]    = digits.fTracks[i];
    }

    return *this;
}
//_____________________________________________________________________________
Int_t AliMUONDigit::Compare(const TObject *obj) const
{
// sort by idDE

 AliMUONDigit* d = (AliMUONDigit*) obj;

 return ( fDetElemId > d->DetElemId()) ? 1 : -1;

}
