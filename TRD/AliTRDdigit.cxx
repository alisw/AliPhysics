

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
Revision 1.5  2000/11/01 14:53:20  cblume
Merge with TRD-develop

Revision 1.1.2.4  2000/10/17 02:27:34  cblume
Get rid of global constants

Revision 1.1.2.3  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.2.2  2000/09/22 14:42:05  cblume
Changed data members to UShort_t

Revision 1.4  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.3  2000/06/07 16:25:37  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.2  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.1.2.1  2000/05/08 14:40:29  cblume
Introduce raw digit bit flag and DecodeAmp()

*/

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
AliTRDdigit::AliTRDdigit():AliDigitNew()
{
  //
  // Default constructor
  //

  fRow  = 0;
  fCol  = 0;
  fTime = 0;

}

//_____________________________________________________________________________
AliTRDdigit::AliTRDdigit(Bool_t isRaw, Int_t *digits, Int_t *amp):AliDigitNew()
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
