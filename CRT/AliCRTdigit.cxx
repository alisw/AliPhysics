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

#include "AliCRTdigit.h"

#include <TArrayF.h>
#include <TArrayI.h>

ClassImp(AliCRTdigit)

//_____________________________________________________________________________
AliCRTdigit::AliCRTdigit()
  : AliDigit(),
    fSector(0),
    fPlate(0),
    fStrip(0),
    fPadx(0),
    fPadz(0),
    fNDigits(0),
    fTdc(0),
    fAdc(0),
    fTracks(0)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliCRTdigit::AliCRTdigit(Int_t tracknum, Int_t *vol, Float_t *digit)
  : AliDigit(),
    fSector(vol[0]),
    fPlate(vol[1]),
    fStrip(vol[2]),
    fPadx(vol[3]),
    fPadz(vol[4]),
    fNDigits(1),
    fTdc(new TArrayF(fNDigits)),
    fAdc(new TArrayF(fNDigits)),
    fTracks(new TArrayI(fNDigits))
{
  
  //
  // Creates CRT digit
  // The creator for the AliCRTdigit class. This routine fills the
  // AliCRTdigit data members from the array digits. 
  //
  (*fTdc)[0] = digit[0];
  (*fAdc)[0] = digit[1];
  //fTracks = new TArrayI(kMAXDIGITS*fNDigits);
  (*fTracks)[0] = tracknum;
  //for (Int_t i = 1; i <kMAXDIGITS*fNDigits; i++) {
  for (Int_t i = 1; i <fNDigits; i++) {
    (*fTracks)[i] = -1;
  }

}

//_____________________________________________________________________________
AliCRTdigit::AliCRTdigit(const AliCRTdigit& digit)
  : AliDigit(digit),
    fSector(digit.fSector),
    fPlate(digit.fPlate),
    fStrip(digit.fStrip),
    fPadx(digit.fPadx),
    fPadz(digit.fPadz),
    fNDigits(digit.fNDigits),
    fTdc(digit.fTdc),  
    fAdc(digit.fAdc),
    fTracks(digit.fTracks)
{
  //
  //-- Copy constructor
  //
}

//_____________________________________________________________________________
AliCRTdigit::~AliCRTdigit()
{
  //
  //
  //
  if ( fTracks ) { delete fTracks; fTracks = 0; }
  if ( fAdc ) { delete fAdc; fAdc = 0; }
  if ( fTdc ) { delete fTdc; fTdc = 0; }
}

//_____________________________________________________________________________
AliCRTdigit& AliCRTdigit::operator=(const AliCRTdigit& digit)
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
  fTracks = digit.fTracks;
  return *this;
}
