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

#include "AliRawEventHeader.h"
#include "AliRawEquipmentHeader.h"
#include "AliRawData.h"

#include "AliRawEvent.h"


ClassImp(AliRawEvent)


//______________________________________________________________________________
AliRawEvent::AliRawEvent()
{
   // Create ALICE event object. If ownData is kFALSE we will use a static
   // raw data object, otherwise a private copy will be made.

   fNSubEvents = 0;
   fEvtHdr     = 0;
   fEqpHdr     = 0;
   fRawData    = 0;
   fSubEvents  = 0;
}

//______________________________________________________________________________
AliRawEvent::AliRawEvent(const AliRawEvent& rawEvent): TObject(rawEvent)
{
// copy constructor

  Fatal("AliRawEvent", "copy constructor not implemented");
}

//______________________________________________________________________________
AliRawEvent& AliRawEvent::operator = (const AliRawEvent& /*rawEvent*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//______________________________________________________________________________
AliRawEventHeader *AliRawEvent::GetHeader()
{
   // Get event header part of AliRawEvent.

   if (!fEvtHdr)
      fEvtHdr = new AliRawEventHeader;

   return fEvtHdr;
}

//______________________________________________________________________________
AliRawEquipmentHeader *AliRawEvent::GetEquipmentHeader()
{
   // Get equipment header part of AliRawEvent.

   if (!fEqpHdr)
      fEqpHdr = new AliRawEquipmentHeader;

   return fEqpHdr;
}

//______________________________________________________________________________
AliRawData *AliRawEvent::GetRawData()
{
   // Get raw data part of AliRawEvent.

   if (!fRawData)
      fRawData = new AliRawData;

   return fRawData;
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

   for (int i = 0; i < fNSubEvents; i++) {
      AliRawEvent *ev = (AliRawEvent *)fSubEvents->At(i);
      ev->GetRawData()->SetSize(0);
   }
   fNSubEvents = 0;
}

//______________________________________________________________________________
AliRawEvent::~AliRawEvent()
{
   // Clean up event object. Delete also, possible, private raw data.

   delete fEvtHdr;
   delete fEqpHdr;
   delete fRawData;
   if (fSubEvents)
      fSubEvents->Delete();
   delete fSubEvents;
}
