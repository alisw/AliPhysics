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

#include "AliCRT.h"
#include "AliCRTdigit.h"
#include "AliRun.h"

ClassImp(AliCRTdigit)

//_____________________________________________________________________________
AliCRTdigit::AliCRTdigit()
{
  // Default ctor.

  fNDigits = 0;
  fTdc = 0;
  fAdc = 0;
  fTracks = 0;
}

//_____________________________________________________________________________
AliCRTdigit::AliCRTdigit(Int_t tracknum, Int_t *vol,Float_t *digit)
{
  
  //
  // Creates CRT digit
  // The creator for the AliCRTdigit class. This routine fills the
  // AliCRTdigit data members from the array digits. 
  //

  fSector = vol[0];
  fPlate  = vol[1];
  fStrip  = vol[2];
  fPadx   = vol[3];
  fPadz   = vol[4];
  fNDigits = 1;
  fTdc = new TArrayF(fNDigits);
  (*fTdc)[0] = digit[0];
  fAdc = new TArrayF(fNDigits);
  (*fAdc)[0] = digit[1];
  //fTracks = new TArrayI(kMAXDIGITS*fNDigits);
  fTracks = new TArrayI(fNDigits);
  (*fTracks)[0] = tracknum;
  //for (Int_t i = 1; i <kMAXDIGITS*fNDigits; i++) {
  for (Int_t i = 1; i <fNDigits; i++) {
    (*fTracks)[i] = -1;
  }

}

//_____________________________________________________________________________
AliCRTdigit::AliCRTdigit(const AliCRTdigit & digit)
{
  //
  //-- Copy ctor.
  //
  fSector = digit.fSector;
  fPlate  = digit.fPlate;
  fStrip  = digit.fStrip;
  fPadx   = digit.fPadx;
  fPadz   = digit.fPadz;
  fNDigits = digit.fNDigits;
  fTdc = new TArrayF(*digit.fTdc);  
  fAdc = new TArrayF(*digit.fAdc);
  fTracks = new TArrayI(*digit.fTracks);

}

//_____________________________________________________________________________
AliCRTdigit& AliCRTdigit::operator= (const AliCRTdigit & digit)
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
  fTdc = new TArrayF(*digit.fTdc);  
  fAdc = new TArrayF(*digit.fAdc);
  fTracks = new TArrayI(*digit.fTracks);

  return *this;
}
