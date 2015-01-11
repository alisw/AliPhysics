#ifndef ALIRAWEVENT_H
#define ALIRAWEVENT_H
// @(#) $Id$
// Author: Fons Rademakers  26/11/99
// Updated: Dario Favretto  15/04/2003

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ROOT_TObjArray
#include <TObjArray.h>
#endif


// Forward class declarations
class AliRawEventHeaderBase;

#include "AliRawVEvent.h"

class AliRawEvent : public AliRawVEvent {

public:
   AliRawEvent();
   virtual ~AliRawEvent();

   virtual AliRawEventHeaderBase *GetHeader();
   virtual Int_t                  GetNEquipments() const { return fNEquipments; }
   virtual AliRawVEquipment       *GetEquipment(Int_t index) const;
   virtual Int_t                  GetNSubEvents() const { return fNSubEvents; }
   virtual AliRawVEvent           *GetSubEvent(Int_t index);

private:
   Int_t                  fNEquipments; // number of valid equipments
   Int_t                  fNSubEvents;  // number of valid sub-events
   AliRawEventHeaderBase *fEvtHdr;      // event header object
   TObjArray             *fEquipments;  // AliRawEquipment's
   TObjArray             *fSubEvents;   // sub AliRawEvent's

   AliRawEvent(const AliRawEvent& rawEvent);
   AliRawEvent& operator = (const AliRawEvent& rawEvent);

   ClassDef(AliRawEvent,3)  // ALICE raw event object
};

#endif
