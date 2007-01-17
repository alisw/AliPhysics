#ifndef ALIRAWEVENTHEADERVERSIONS_H
#define ALIRAWEVENTHEADERVERSIONS_H

#include "AliRawEventHeaderBase.h"

START_EVENT_HEADER(3,1)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern[1]; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id

END_EVENT_HEADER(3,1)

START_EVENT_HEADER(3,2)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern[1]; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id

END_EVENT_HEADER(3,2)

START_EVENT_HEADER(3,3)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern[1]; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id

END_EVENT_HEADER(3,3)

START_EVENT_HEADER(3,4)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern[1]; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id
   UInt_t fTimestamp;     // event timestamp

END_EVENT_HEADER(3,4)

START_EVENT_HEADER(3,5)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern[1]; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id
   UInt_t fTimestamp;     // event timestamp

END_EVENT_HEADER(3,5)

START_EVENT_HEADER(3,6)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern[1]; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id
   UInt_t fTimestamp;     // event timestamp

END_EVENT_HEADER(3,6)

START_EVENT_HEADER(3,7)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern[1]; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id
   UInt_t fTimestamp;     // event timestamp

END_EVENT_HEADER(3,7)

START_EVENT_HEADER(3,8)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern[1]; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id
   UInt_t fTimestamp;     // event timestamp

END_EVENT_HEADER(3,8)

#endif
