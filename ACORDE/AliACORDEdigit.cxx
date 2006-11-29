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

////////////////////////////////////////////////////////////////////////////
//  ACORDE digit: Id
//
// The digits are made in FinishEvent() by summing all the hits in a
// counter.
//   The main parts of the code need to be written.
//
////////////////////////////////////////////////////////////////////////////

#include "AliACORDEdigit.h"

#include <TArrayF.h>
#include <TArrayI.h>

ClassImp(AliACORDEdigit)

//_____________________________________________________________________________
AliACORDEdigit::AliACORDEdigit()
  : AliDigit(),
    fSector(0),
    fPlate(0),
    fStrip(0),
    fPadx(0),
    fPadz(0),
    fNDigits(0),
    fTdc(0),
    fAdc(0)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliACORDEdigit::AliACORDEdigit(Int_t* tracks, Int_t *vol, Float_t *digit)
  : AliDigit(tracks),
    fSector(vol[0]),
    fPlate(vol[1]),
    fStrip(vol[2]),
    fPadx(vol[3]),
    fPadz(vol[4]),
    fNDigits(1),
    fTdc(new TArrayF(fNDigits)),
    fAdc(new TArrayF(fNDigits))
{
  
  //
  // Creates ACORDE digit
  // The creator for the AliACORDEdigit class. This routine fills the
  // AliACORDEdigit data members from the array digits. 
  //
  (*fTdc)[0] = digit[0];
  (*fAdc)[0] = digit[1];
}

//_____________________________________________________________________________
AliACORDEdigit::AliACORDEdigit(const AliACORDEdigit& digit)
  : AliDigit(digit),
    fSector(digit.fSector),
    fPlate(digit.fPlate),
    fStrip(digit.fStrip),
    fPadx(digit.fPadx),
    fPadz(digit.fPadz),
    fNDigits(digit.fNDigits),
    fTdc(digit.fTdc),  
    fAdc(digit.fAdc)
{
  //
  //-- Copy constructor
  //
}

//_____________________________________________________________________________
AliACORDEdigit::~AliACORDEdigit()
{
  //
  //
  //
  if ( fAdc ) { delete fAdc; fAdc = 0; }
  if ( fTdc ) { delete fTdc; fTdc = 0; }
}

//_____________________________________________________________________________
AliACORDEdigit& AliACORDEdigit::operator=(const AliACORDEdigit& digit)
{
  //
  //-- Asingment operator.
  //
  fSector = digit.fSector;
  fPlate  = digit.fPlate;
  fStrip  = digit.fStrip;
  fPadx   = digit.fPadx;
  fPadz   = digit.fPadz;
  fNDigits = digit.fNDigits;
  fTdc = digit.fTdc;
  fAdc = digit.fAdc;
  return *this;
}
