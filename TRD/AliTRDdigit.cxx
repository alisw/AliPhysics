

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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  The TRD digit                                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDdigit.h"

ClassImp(AliTRDdigit)

//_____________________________________________________________________________
  
  // Marks a raw digit
  const UInt_t AliTRDdigit::fgkRawDigit = 0x00000001; 

//_____________________________________________________________________________
AliTRDdigit::AliTRDdigit()
  :AliDigitNew()
  ,fRow(0)
  ,fCol(0)
  ,fTime(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDdigit::AliTRDdigit(Bool_t isRaw, Int_t *digits, Int_t *amp)
  :AliDigitNew()
  ,fRow(0)
  ,fCol(0)
  ,fTime(0)
{
  //
  // Create a TRD digit
  //

  // Store the volume hierarchy
  fId   = digits[0];

  // Store the row, pad, and time bucket number
  fRow  = digits[1];
  fCol  = digits[2];
  fTime = digits[3];

  // Store the signal amplitude
  fAmp  = amp[0];

  if (isRaw) SetBit(fgkRawDigit);

}

//_____________________________________________________________________________
AliTRDdigit::~AliTRDdigit()
{
  //
  // AliTRDdigit destructor
  //

}

//_____________________________________________________________________________
Int_t AliTRDdigit::DecodeAmp() const
{
  //
  // Decodes the digit amplitude
  //

  return 0;

}
