// @(#) $Id$
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

#include <TObjArray.h>
#include "TBuffer.h"

#include "AliLog.h"

#include "AliRawEventHeaderBase.h"
#include "AliRawEquipment.h"

#include "AliRawEvent.h"


ClassImp(AliRawEvent)


//______________________________________________________________________________
AliRawEvent::AliRawEvent():
fNEquipments(0),
fNSubEvents(0),
fEvtHdr(NULL),
fEquipments(NULL),
fSubEvents(NULL)
{
   // Create ALICE event object. If ownData is kFALSE we will use a static
   // raw data object, otherwise a private copy will be made.

}

//______________________________________________________________________________
AliRawEventHeaderBase *AliRawEvent::GetHeader()
{
  if (!fEvtHdr) {
    AliFatal("Event header does not exist!");
    return 0x0;
  }

  return fEvtHdr;
}

//______________________________________________________________________________
AliRawVEquipment *AliRawEvent::GetEquipment(Int_t index) const
{
   // Get specified equipment. Returns 0 if equipment does not exist.

   if (!fEquipments)
      return 0;

   return (AliRawEquipment *) fEquipments->At(index);
}

//______________________________________________________________________________
AliRawVEvent *AliRawEvent::GetSubEvent(Int_t index)
{
   // Get specified sub event. Returns 0 if sub event does not exist.

   if (!fSubEvents)
      return 0;

   return (AliRawEvent *) fSubEvents->At(index);
}

//______________________________________________________________________________
AliRawEvent::~AliRawEvent()
{
   // Clean up event object. Delete also, possible, private raw data.

   delete fEvtHdr;
   if (fEquipments)
      fEquipments->Delete();
   delete fEquipments;
   if (fSubEvents)
      fSubEvents->Delete();
   delete fSubEvents;
}

//______________________________________________________________________________
void AliRawEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliRawEvent.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fNEquipments;
      R__b >> fNSubEvents;
      R__b >> fEvtHdr;
      R__b >> fEquipments;
      R__b >> fSubEvents;
      R__b.CheckByteCount(R__s, R__c, AliRawEvent::IsA());
   } else {
      R__c = R__b.WriteVersion(AliRawEvent::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << fNEquipments;
      R__b << fNSubEvents;
      R__b << fEvtHdr;
      R__b << fEquipments;
      R__b << fSubEvents;
      R__b.SetByteCount(R__c, kTRUE);
   }
}
