/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawEquipmentHeader                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <Bytes.h>

#include "AliRawEquipmentHeader.h"

ClassImp(AliRawEquipmentHeader)

//______________________________________________________________________________
AliRawEquipmentHeader::AliRawEquipmentHeader():
  fSize(0),
  fEquipmentType(0),
  fEquipmentID(0),
  fBasicElementSizeType(0)
{
  // Default constructor
  for(Int_t i = 0; i < kAttributeWords; i++)
    fTypeAttribute[i] = 0;
}

//______________________________________________________________________________
void AliRawEquipmentHeader::Swap()
{
   // Swap equipment header data. There is no way to see if the data
   // has already been swapped. This method is only called when the
   // header is read from the DATE event builder (GDC).

   fSize                 = net2host(fSize);
   fEquipmentType        = net2host(fEquipmentType);
   fEquipmentID          = net2host(fEquipmentID);
   fBasicElementSizeType = net2host(fBasicElementSizeType);
   for (int i = 0; i < kAttributeWords; i++)
      fTypeAttribute[i] = net2host(fTypeAttribute[i]);
}

//______________________________________________________________________________
void AliRawEquipmentHeader::Reset()
{
  // Reset the contents of the equipment
  // header data
  fSize = fEquipmentType = fEquipmentID = fBasicElementSizeType = 0;

  for(Int_t i = 0; i < kAttributeWords; i++)
    fTypeAttribute[i] = 0;
}
