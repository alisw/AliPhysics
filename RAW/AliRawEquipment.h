#ifndef ALIRAWEQUIPMENT_H
#define ALIRAWEQUIPMENT_H
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
#include <TRef.h>

// Forward class declarations
class AliRawDataArray;
class AliRawEquipmentHeader;
class AliRawData;

#include "AliRawVEquipment.h"

class AliRawEquipment : public AliRawVEquipment {

public:
   AliRawEquipment();
   virtual ~AliRawEquipment();

   virtual AliRawEquipmentHeader *GetEquipmentHeader();
   virtual AliRawData            *GetRawData();

private:
   AliRawEquipmentHeader *fEqpHdr;      // equipment header
   AliRawData            *fRawData;     // raw data container
   TRef                   fRawDataRef;  // reference to raw data container

   AliRawEquipment(const AliRawEquipment& rawEvent);
   AliRawEquipment& operator = (const AliRawEquipment& rawEvent);

   ClassDef(AliRawEquipment,3)  // ALICE raw equipment object
};

#endif
