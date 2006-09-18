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
// The AliRunDB class provides the interface to the run and file        //
// catalogues (AliEn or plain MySQL).                                   //
// The AliStats class provides statics information that is added as     //
// a single keyed object to each raw file.                              //
// The AliTagDB provides an interface to a TAG database.                //
// The AliMDC class is usid by the "alimdc" stand-alone program         //
// that reads data directly from DATE.                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>

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
AliRawEventHeaderBase *AliRawEvent::GetHeader(char*& data)
{
  // Get event header part of AliRawEvent.
  // First the DATE version is identified and then the
  // corresponding event header version object is created
  
  if (!fEvtHdr) {
    fEvtHdr = AliRawEventHeaderBase::Create(data);
  }

  return fEvtHdr;
}

//______________________________________________________________________________
AliRawEventHeaderBase *AliRawEvent::GetHeader()
{
  if (!fEvtHdr) {
      AliFatal("Header version not yet initialized!");
      return 0x0;
    }

  return fEvtHdr;
}

//______________________________________________________________________________
AliRawEquipment *AliRawEvent::NextEquipment()
{
   // Returns next equipment object.

   if (!fEquipments)
      fEquipments = new TObjArray(100); // arbitrary, probably enough to prevent resizing

   if (fEquipments->GetSize() <= fNEquipments) {
      fEquipments->Expand(fNEquipments+10);
      Warning("NextEquipment", "expanded fEquipments by 10 to %d",
              fEquipments->GetSize());
   }

   AliRawEquipment *eq;
   if (!(eq = (AliRawEquipment *)fEquipments->At(fNEquipments))) {
      eq = new AliRawEquipment;
      fEquipments->AddAt(eq, fNEquipments);
   }

   fNEquipments++;

   return eq;
}

//______________________________________________________________________________
AliRawEquipment *AliRawEvent::GetEquipment(Int_t index) const
{
   // Get specified equipment. Returns 0 if equipment does not exist.

   if (!fEquipments)
      return 0;

   return (AliRawEquipment *) fEquipments->At(index);
}

//______________________________________________________________________________
AliRawEvent *AliRawEvent::NextSubEvent()
{
   // Returns next sub-event object.

   if (!fSubEvents)
      fSubEvents = new TObjArray(100); // arbitrary, probably enough to prevent resizing

   if (fSubEvents->GetSize() <= fNSubEvents) {
      fSubEvents->Expand(fNSubEvents+10);
      Warning("NextSubEvent", "expanded fSubEvents by 10 to %d",
              fSubEvents->GetSize());
   }

   AliRawEvent *ev;
   if (!(ev = (AliRawEvent *)fSubEvents->At(fNSubEvents))) {
      ev = new AliRawEvent;
      fSubEvents->AddAt(ev, fNSubEvents);
   }

   fNSubEvents++;

   return ev;
}

//______________________________________________________________________________
AliRawEvent *AliRawEvent::GetSubEvent(Int_t index) const
{
   // Get specified sub event. Returns 0 if sub event does not exist.

   if (!fSubEvents)
      return 0;

   return (AliRawEvent *) fSubEvents->At(index);
}

//______________________________________________________________________________
void AliRawEvent::Reset()
{
   // Reset the event in case it needs to be re-used (avoiding costly
   // new/delete cycle). We reset the size marker for the AliRawData
   // objects and the sub event counter.

   for (int i = 0; i < fNEquipments; i++) {
      AliRawEquipment *eq = (AliRawEquipment *)fEquipments->At(i);
      eq->Reset();
   }
   fNEquipments = 0;
   for (int i = 0; i < fNSubEvents; i++) {
      AliRawEvent *ev = (AliRawEvent *)fSubEvents->At(i);
      ev->Reset();
   }
   fNSubEvents = 0;
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
