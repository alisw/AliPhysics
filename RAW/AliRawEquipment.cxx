// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawEvent                                                          //
//                                                                      //
// Set of classes defining the ALICE RAW event format. The AliRawEvent  //
// class defines a RAW event. It consists of an AliEventHeader object   //
// an AliEquipmentHeader object, an AliRawData object and an array of   //
// sub-events, themselves also being AliRawEvents. The number of        //
// sub-events depends on the number of DATE LDC's.                      //
// The AliRawEvent objects are written to a ROOT file using different   //
// technologies, i.e. to local disk via AliRawDB or via rfiod using     //
// AliRawRFIODB or via rootd using AliRawRootdDB or to CASTOR via       //
// rootd using AliRawCastorDB (and for performance testing there is     //
// also AliRawNullDB).                                                  //
// The AliStats class provides statics information that is added as     //
// a single keyed object to each raw file.                              //
// The AliTagDB provides an interface to a TAG database.                //
// The AliMDC class is usid by the "alimdc" stand-alone program         //
// that reads data directly from DATE.                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliRawEquipmentHeader.h"
#include "AliRawData.h"

#include "AliRawEquipment.h"


ClassImp(AliRawEquipment)


//______________________________________________________________________________
AliRawEquipment::AliRawEquipment():
fEqpHdr(NULL),
fRawData(NULL)
{
   // Create ALICE equipment object.

}

//______________________________________________________________________________
AliRawEquipmentHeader *AliRawEquipment::GetEquipmentHeader()
{
   // Get equipment header part of AliRawEquipment.

   if (!fEqpHdr)
      fEqpHdr = new AliRawEquipmentHeader;

   return fEqpHdr;
}

//______________________________________________________________________________
AliRawData *AliRawEquipment::GetRawData()
{
   // Get raw data part of AliRawEquipment.

   if (!fRawData)
      fRawData = new AliRawData;

   return fRawData;
}

//______________________________________________________________________________
void AliRawEquipment::Reset()
{
   // Reset the equipment in case it needs to be re-used (avoiding costly
   // new/delete cycle). We reset the size marker for the AliRawData
   // object.

   if (fEqpHdr) fEqpHdr->Reset();
   GetRawData()->SetSize(0);
}

//______________________________________________________________________________
AliRawEquipment::~AliRawEquipment()
{
   // Clean up event object. Delete also, possible, private raw data.

   delete fEqpHdr;
   delete fRawData;
}
