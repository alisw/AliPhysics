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

///////////////////////////////////////////////////////////////////////////////
//
// This class which defines the trigger bunch-crossing mask
//
//
///////////////////////////////////////////////////////////////////////////////

#include "AliTriggerBCMask.h"

ClassImp(AliTriggerBCMask)

//_____________________________________________________________________________
AliTriggerBCMask::AliTriggerBCMask():
  TNamed()
{
  // Default constructor
  for (Int_t i = 0; i < kNBytesPerBCMask; i++) fBCMask[i] = 0;
}

//_____________________________________________________________________________
AliTriggerBCMask::AliTriggerMask( TString & name, UChar_t *mask ):
  TNamed( name, name )
{
  // Constructor
  for (Int_t i = 0; i < kNBytesPerBCMask; i++) fBCMask[i] = mask[i];
}
//_____________________________________________________________________________
AliTriggerBCMask::~AliTriggerBCMask() 
{ 
  // Destructor
}
//_____________________________________________________________________________
AliTriggerBCMask::AliTriggerBCMask( const AliTriggerBCMask& mask ):
  TNamed( mask ),
{
   // Copy constructor
  for (Int_t i = 0; i < kNBytesPerBCMask; i++) fBCMask[i] = mask.fBCMask[i];
}

//______________________________________________________________________________
AliTriggerBCMask& AliTriggerBCMask::operator=(const AliTriggerBCMask& mask)
{
   // AliTriggerBCMask assignment operator.

   if (this != &mask) {
      TNamed::operator=(mask);
      for (Int_t i = 0; i < kNBytesPerBCMask; i++) fBCMask[i] = mask.fBCMask[i];
   }
   return *this;
}

//_____________________________________________________________________________
Bool_t AliTriggerBCMask::GetMask( UShort_t index)
{
  // Return true or false whenever the mask is active
  // for the bunch-crossing # = index
  UShort_t position = index/8;
  if (position >= kNBytesPerBCMask) return kFALSE;
  UChar_t offset = index%8;
  return (fBCMask[position] & (0x1 << offset));
}
