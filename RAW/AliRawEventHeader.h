#ifndef ALIRAWEVENTHEADER_H
#define ALIRAWEVENTHEADER_H
// @(#)alimdc:$Name$:$Id$
// Author: Fons Rademakers  26/11/99

/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRawEventHeader                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#endif


class AliRawEventHeader : public TObject {

public:
   AliRawEventHeader() { fSize = 0; }
   virtual ~AliRawEventHeader() { }

   void         *HeaderBegin() { return (void *) &fSize; }
   Int_t         HeaderSize() const { return (Long_t) &fGDCId - (Long_t) &fSize + sizeof(fGDCId); }
   Bool_t        DataIsSwapped() const;
   Bool_t        IsSwapped() const { return (fMagic == fgkEventMagicNumberSwapped) ? kTRUE : kFALSE; }
   Bool_t        IsValid() const { return IsSwapped() ? kTRUE : ((fMagic == fgkEventMagicNumber) ? kTRUE : kFALSE); }
   void          Swap();

   UInt_t        GetEventSize() const { return fSize; }
   UInt_t        GetMagic() const { return fMagic; }
   UInt_t        GetHeaderLength() const { return fHeadLen; }
   UInt_t        GetVersion() const { return fVersion; }
   UInt_t        GetType() const { return fType; }
   const char   *GetTypeName() const;
   UInt_t        GetRunNumber() const { return fRunNb; }
   UInt_t        GetEventInRun() const;
   const UInt_t *GetId() const { return fId; }
   const UInt_t *GetTriggerPattern() const { return fTriggerPattern; }
   const UInt_t *GetDetectorPattern() const { return fDetectorPattern; }
   const UInt_t *GetTypeAttribute() const { return fTypeAttribute; }
   UInt_t        GetLDCId() const { return fLDCId; }
   UInt_t        GetGDCId() const { return fGDCId; }

   // The following enumeration can be used once the kEventTypeMask has been
   // applied to the raw event type
   enum EAliRawEventType {
     kStartOfRun =       1,    // START_OF_RUN
     kEndOfRun =         2,    // END_OF_RUN
     kStartOfRunFiles =  3,    // START_OF_RUN_FILES
     kEndOfRunFiles =    4,    // END_OF_RUN_FILES
     kStartOfBurst =     5,    // START_OF_BURST
     kEndOfBurst =       6,    // END_OF_BURST
     kPhysicsEvent =     7,    // PHYSICS_EVENT
     kCalibrationEvent = 8,    // CALIBRATION_EVENT
     kFormatError =      9     // EVENT_FORMAT_ERROR
   };

   // Type sizes
   enum {
     kIdWords        = 2,
     kTriggerWords   = 2,
     kDetectorWords  = 1,
     kAttributeWords = 3
   };

private:
   UInt_t fSize;          // size of event in bytes
   UInt_t fMagic;         // magic number used for consistency check
   UInt_t fHeadLen;       // size of header in bytes
   UInt_t fVersion;       // unique version identifier
   UInt_t fType;          // event type
   UInt_t fRunNb;         // run number
   UInt_t fId[kIdWords];  // id field
   UInt_t fTriggerPattern[kTriggerWords];   // trigger pattern
   UInt_t fDetectorPattern[kDetectorWords]; // detector pattern
   UInt_t fTypeAttribute[kAttributeWords];  // system (0,1) and user (2) attributes
   UInt_t fLDCId;         // LDC id
   UInt_t fGDCId;         // GDC id

   static const Int_t fgkEventTypeMin = kStartOfRun;   // minimal event type
   static const Int_t fgkEventTypeMax = kFormatError;  // maximal event type

   static const UInt_t fgkEventMagicNumber        = 0xDA1E5AFE; // magic word
   static const UInt_t fgkEventMagicNumberSwapped = 0xFE5A1EDA; // swapped magic word

   ClassDef(AliRawEventHeader,1)  // Alice raw event header
};

#endif
