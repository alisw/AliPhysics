/****************************************************************************
 *
 * event.h: DATE event data format
 *
 * Revision History:
 *    V01.00  RD PVV    09-Jan-97 Initial version
 *    V01.01  AV        24-Feb-97 Added START_OF_RUN_FILES and triggerNb
 *    V02.02  RD        13-Mar-97 Detector ID mask type added
 *    V02.03  PVV       20-Mar-97 Detector ID on 128 bits
 *    V02.03  RD PVV    20-Mar-97 Added EVENT_DATA_TRUNCATED
 *    V02.04  PVV       06-May-97 Added EVENT_TYPE_MASK
 *    V02.05  RD PVV    09-May-97 Increase EVENT_DATA_MAX_SIZE to 50 K
 *    V02.06  WB MG     22-May-97 Added END_OF_RUN_FILES
 *    V02.07  WB        23-May-97 Added errorCode, deadTime, deadTimeusec
 *                                EVENT_DATA_MAX_SIZE set to 100 * 1024
 *                                MAX_DETECTORS set to 126
 *    V02.08  PVV       02-Jun-97 Modify the encoding of types
 *    V02.09  WB RD PVV 28-Jul-98 Add fileSeqNb in the header.
 *                                Reduce detector mask to 3 long
 *    V02.10  RD        31-Jul-98 (start|end)OfRunFiles added to eventType
 *    V02.11  PVV RD    02-Sep-98 Event type re-numbered
 *                                Equipment bit added to event type
 *    V02.12  AV PVV RD 03-Sep-98 FileSeqNo moved before detectorId
 *    V02.13  RD        08-Oct-98 32 bits fields defined as long32
 *    V02.13  RD        19-Feb-99 Endianness/swap definitions added
 *    V02.14  WB PVV RD 21-Jun-99 typeAttribute added
 *    V02.15  RD        27-Jul-99 Macros for typeAttribute handling added
 *    V02.16  RD        19-Nov-99 Bug in Attributes test/set fixed
 *    V02.17  WB PVV RD 08-May-00 System attributes for SOR and EOR added
 *    V02.18  RD        18-May-00 EVENT_H_ID added
 *
 ***************************************************************************/
/* modified for root convention */


#ifndef __event_h__
#define __event_h__

/* Unique version identifier */
#define EVENT_H_ID 0x00020018

/* Possible Values for type in eventHeaderStruct */

   /* Mask to separate the event type from the event flags */
#define EVENT_TYPE_MASK      ((UInt_t)0x0000FFFF)
#define EVENT_FLAGS_MASK     ((UInt_t)~EVENT_TYPE_MASK)

   /* Event error flags: */
#define EVENT_ERROR          ((UInt_t)0x80000000)
#define EVENT_DATA_TRUNCATED ((UInt_t)0x40000000)

   /* Event flags: */
#define EVENT_EQUIPMENT      ((UInt_t)0x00010000)
#define EVENT_SWAPPED        ((UInt_t)0x00020000)

   /* Event types: */
#define START_OF_RUN         ((UInt_t)1)
#define END_OF_RUN           ((UInt_t)2)
#define START_OF_RUN_FILES   ((UInt_t)3)
#define END_OF_RUN_FILES     ((UInt_t)4)
#define START_OF_BURST       ((UInt_t)5)
#define END_OF_BURST         ((UInt_t)6)
#define PHYSICS_EVENT        ((UInt_t)7)
#define CALIBRATION_EVENT    ((UInt_t)8)
#define END_OF_LINK          ((UInt_t)9)
#define EVENT_FORMAT_ERROR   ((UInt_t)10)

#define EVENT_TYPE_MIN       1
#define EVENT_TYPE_MAX       10

   /* System attributes */
#define START_OF_RUN_START   0
#define START_OF_RUN_END     1
#define END_OF_RUN_START     0
#define END_OF_RUN_END       1
#define NO_EVENT_BUILDING    2

/* The following enumeration can be used once the EVENT_TYPE_MASK has been
 * applied to the raw event type */
enum eventType {
  startOfRun =       START_OF_RUN,
  endOfRun =         END_OF_RUN,
  startOfRunFiles =  START_OF_RUN_FILES,
  endOfRunFiles =    END_OF_RUN_FILES,
  startOfBurst =     START_OF_BURST,
  endOfBurst =       END_OF_BURST,
  physicsEvent =     PHYSICS_EVENT,
  calibrationEvent = CALIBRATION_EVENT,
  endOfLink =        END_OF_LINK,
  formatError =      EVENT_FORMAT_ERROR
};

#define EVENT_MAGIC_NUMBER         ((UInt_t)0xDA1E5AFE)
#define EVENT_MAGIC_NUMBER_SWAPPED ((UInt_t)0xFE5A1EDA)

#define HEADER_LENGTH        sizeof( struct eventHeaderStruct )

/* Maximum size for any event (excluding the header) */
#define EVENT_DATA_MAX_SIZE  ( 100 * 1024 )

/* Maximum size for any event (including the header) */
#define EVENT_MAX_SIZE       ( HEADER_LENGTH + EVENT_DATA_MAX_SIZE )

/* The detectorId mask */
typedef UInt_t detectorIdType;
#define MAX_DETECTOR         ((int)94)
#define SUPER_EVENT_MASK     ((int)(1 << 31))
#define SUPER_EVENT_BIT      ((int)(MAX_DETECTOR+1))
#define MASK_LENGTH          ((int)(MAX_DETECTOR/32)+1)

/* The typeAttribute mask */
#define ATTRIBUTE_WORDS 2
#define ATTRIBUTE_BYTES ( ATTRIBUTE_WORDS * 4 )
#define ATTRIBUTE_BITS  ( ATTRIBUTE_BYTES * 8 )
#define ATTRIBUTE_TO_NUM( bit ) (((bit)&0x20)>>5)
#define ATTRIBUTE_TO_BIT( bit ) (1<<((bit)&0x1f))
#define SET_ATTRIBUTE( mask, bit ) \
  (mask)[ ATTRIBUTE_TO_NUM(bit) ] |= ATTRIBUTE_TO_BIT( bit )
#define CLEAR_ATTRIBUTE( mask, bit ) \
  (mask)[ ATTRIBUTE_TO_NUM(bit) ] &= ~(ATTRIBUTE_TO_BIT( bit ))
#define TEST_ATTRIBUTE( mask, bit ) \
  (((mask)[ ATTRIBUTE_TO_NUM(bit) ] & ATTRIBUTE_TO_BIT( bit )) != 0 )

struct eventHeaderStruct { 
  Int_t eventSize;          /* size of event in Bytes  */
  UInt_t eventMagic;         /* magic number used for consistency check */
  UInt_t eventType;          /* event type */
  UInt_t eventHeadSize;       /* size of header in bytes */
  UInt_t eventRunNb;         /* run number */
  UInt_t burstNb;       /* burst number */
  UInt_t nbInRun;       /* event number in run */
  UInt_t nbInBurst;     /* event number in burst */
  UInt_t triggerNb;     /* trigger number for this detector */
  UInt_t fileSeqNb;     /* file sequence number for multifiles run */
  detectorIdType  detectorId[MASK_LENGTH]; /* detector identification */
  UInt_t time;          /* time in seconds since 0.00 GMT 1.1.1970 */
  UInt_t usec;          /* microseconds */
  UInt_t errorCode;
  UInt_t deadTime;
  UInt_t deadTimeusec;
  UInt_t eventTypeAttribute[ATTRIBUTE_WORDS]; /* event type id mask */
};

struct eventStruct {
    struct eventHeaderStruct eventHeader;
    unsigned short rawData[1];
};
#endif
