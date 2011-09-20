#ifndef ALIRAWEVENTHEADERVERSIONS_H
#define ALIRAWEVENTHEADERVERSIONS_H

#include "AliRawEventHeaderBase.h"

#undef INIT_HEADER_VARS
#define INIT_HEADER_VARS fType(0), fRunNb(0), fDetectorPattern(0), fLdcId(0), fGdcId(0) 
START_EVENT_HEADER(3,1)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id

END_EVENT_HEADER(3,1)

START_EVENT_HEADER(3,2)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id

END_EVENT_HEADER(3,2)

START_EVENT_HEADER(3,3)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id

END_EVENT_HEADER(3,3)

#undef INIT_HEADER_VARS
#define INIT_HEADER_VARS fType(0), fRunNb(0), fDetectorPattern(0), fLdcId(0), fGdcId(0), fTimestamp(0)
START_EVENT_HEADER(3,4)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern; // detector pattern
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
   UInt_t fDetectorPattern; // detector pattern
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
   UInt_t fDetectorPattern; // detector pattern
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
   UInt_t fDetectorPattern; // detector pattern
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
   UInt_t fDetectorPattern; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id
   UInt_t fTimestamp;     // event timestamp

END_EVENT_HEADER(3,8)

START_EVENT_HEADER(3,9)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id
   UInt_t fTimestamp;     // event timestamp

END_EVENT_HEADER(3,9)

#undef INIT_HEADER_VARS
#define INIT_HEADER_VARS fType(0), fRunNb(0), fDetectorPattern(0), fLdcId(0), fGdcId(0), fTimestamp(0), fTimestampUsec(0)
START_EVENT_HEADER(3,11)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id
   UInt_t fTimestamp;     // event timestamp
   UInt_t fTimestampUsec; // event timestamp (microseconds)

END_EVENT_HEADER(3,11)

START_EVENT_HEADER(3,12)

   UInt_t fType;          // event type
   UInt_t fRunNb;     // run number
   UInt_t fId[2];  // id field
   UInt_t fTriggerPattern[2];   // trigger pattern
   UInt_t fDetectorPattern; // detector pattern
   UInt_t fTypeAttribute[3];  // system (0,1) and user (2) attributes
   UInt_t fLdcId;         // LDC id
   UInt_t fGdcId;         // GDC id
   UInt_t fTimestamp;     // event timestamp
   UInt_t fTimestampUsec; // event timestamp (microseconds)

END_EVENT_HEADER(3,12)

#endif
