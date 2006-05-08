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
 *    V02.19  RD        10-Aug-00 Macros for detectorId handling added
 *    V03.00  RD        23-Nov-00 Version for DATE V4
 *    V03.01  AV KS     05-Apr-02 Introduction of eventLocationDescriptorStruct
 *    V03.02  KS        10-May-02 Added ATTR_KEEP_PAGES for COLE
 *            RD        30-Apr-04 Added definitions for the Common Data Header
 *	      RD	24-Jun-04 Added definitions for HLT DECISION
 *            RD        13-Jul-04 Added macros to OR attributes
 *    V03.03  RD        25-May-05 Added eventTimestamp
 *    V03.04  RD	17-Aug-05 Added VANGUARD and REARGUARD events
 *    V03.05  RD        05-Sep-05 Added SYSTEM and DETECTOR software tri evts
 *    V03.06  RD        14-Sep-05 VANGUARD/REARGUARD changed into START_OF_DATA
 *                                and END_OF_DATA
 *
 * Preprocessor definitions:
 *  NDEBUG  Define BEFORE including this file to disable run-time checks on
 *          various parameters
 *
 * Related facilities:
 *  validateEvent.c  Validation program, should be run after EACH change to
 *                   the definitions given here below
 ***************************************************************************/
#ifndef __event_h__
#define __event_h__

#define EVENT_MAJOR_VERSION_NUMBER 0x0003
#define EVENT_MINOR_VERSION_NUMBER 0x0006

/* ========== System includes ========= */
#include <string.h> /* Needed by: memset, memcpy */
#include <assert.h> /* Needed by: assert */

/* ========== Definitions for the event header ========== */

/* ---------- Header base size ---------- */
/* This value must be updated for each change in the eventHeaderStruct */
#define EVENT_HEAD_BASE_SIZE 68

/* ---------- Event size ---------- */
typedef unsigned long32 eventSizeType;

/* ---------- Magic signature and its byte-swapped version ---------- */
typedef unsigned long32 eventMagicType;
#define EVENT_MAGIC_NUMBER         ((eventMagicType)0xDA1E5AFE)
#define EVENT_MAGIC_NUMBER_SWAPPED ((eventMagicType)0xFE5A1EDA)

/* ---------- Header size ---------- */
typedef unsigned long32 eventHeadSizeType;

/* ---------- Unique version identifier ---------- */
#define EVENT_CURRENT_VERSION \
  (((EVENT_MAJOR_VERSION_NUMBER<<16)&0xffff0000)|\
   (EVENT_MINOR_VERSION_NUMBER&0x0000ffff))
typedef unsigned long32 eventVersionType;

/* ---------- Event type ---------- */
typedef unsigned long32 eventTypeType;
#define START_OF_RUN                    ((eventTypeType) 1)
#define END_OF_RUN                      ((eventTypeType) 2)
#define START_OF_RUN_FILES              ((eventTypeType) 3)
#define END_OF_RUN_FILES                ((eventTypeType) 4)
#define START_OF_BURST                  ((eventTypeType) 5)
#define END_OF_BURST                    ((eventTypeType) 6)
#define PHYSICS_EVENT                   ((eventTypeType) 7)
#define CALIBRATION_EVENT               ((eventTypeType) 8)
#define EVENT_FORMAT_ERROR              ((eventTypeType) 9)
#define START_OF_DATA                   ((eventTypeType)10)
#define END_OF_DATA                     ((eventTypeType)11)
#define SYSTEM_SOFTWARE_TRIGGER_EVENT   ((eventTypeType)12)
#define DETECTOR_SOFTWARE_TRIGGER_EVENT ((eventTypeType)13)
#define EVENT_TYPE_MIN                  1
#define EVENT_TYPE_MAX                  13
enum eventTypeEnum {
  startOfRun                   = START_OF_RUN,
  endOfRun                     = END_OF_RUN,
  startOfRunFiles              = START_OF_RUN_FILES,
  endOfRunFiles                = END_OF_RUN_FILES,
  startOfBurst                 = START_OF_BURST,
  endOfBurst                   = END_OF_BURST,
  physicsEvent                 = PHYSICS_EVENT,
  calibrationEvent             = CALIBRATION_EVENT,
  formatError                  = EVENT_FORMAT_ERROR,
  startOfData                  = START_OF_DATA,
  endOfData                    = END_OF_DATA,
  systemSoftwareTriggerEvent   = SYSTEM_SOFTWARE_TRIGGER_EVENT,
  detectorSoftwareTriggerEvent = DETECTOR_SOFTWARE_TRIGGER_EVENT
};
#define EVENT_TYPE_OK(t) (((t) >= EVENT_TYPE_MIN) && (((t) <= EVENT_TYPE_MAX)))

/* ---------- Run number ---------- */
typedef unsigned long32 eventRunNbType;

/* ---------- The eventId field ---------- */
#define EVENT_ID_BYTES 8
#define EVENT_ID_WORDS ((EVENT_ID_BYTES) >> 2)
typedef unsigned long32 eventIdType[EVENT_ID_WORDS];

   /* PERIOD - ORBIT - BUNCH crossing type events */
#define EVENT_ID_MAX_PERIOD         0x0fffffff
#define EVENT_ID_MAX_ORBIT          0x00ffffff
#define EVENT_ID_MAX_BUNCH_CROSSING 0x00000fff
#define LOAD_EVENT_ID(id,s,o,bc)       \
  (EVENT_ID_SET_PERIOD(id,s),          \
   EVENT_ID_SET_ORBIT(id,o),           \
   EVENT_ID_SET_BUNCH_CROSSING(id,bc))
#define EVENT_ID_GET_BUNCH_CROSSING(id) ((id)[1]&0x00000fff)
#define EVENT_ID_GET_ORBIT(id) \
                     ((((id)[0]<<20)&0xf00000)|(((id)[1]>>12)&0xfffff))
#define EVENT_ID_GET_PERIOD(id) (((id)[0]>>4)&0x0fffffff)

#define EVENT_ID_SET_BUNCH_CROSSING(id,v)                \
  (assert(((v)>=0)&&((v)<=EVENT_ID_MAX_BUNCH_CROSSING)), \
   (id)[1]=((id)[1]&0xfffff000)|((v)&0xfff))
#define EVENT_ID_SET_ORBIT(id,v) \
  (assert(((v)>=0)&&((v)<=EVENT_ID_MAX_ORBIT)),           \
   (id)[0]=(((id)[0])&0xfffffff0)|(((v)&0x00f00000)>>20), \
   (id)[1]=(((id)[1])&0x00000fff)|(((v)&0x000fffff)<<12))
#define EVENT_ID_SET_PERIOD(id,v)                         \
  (assert(((v)>=0)&&((v)<=EVENT_ID_MAX_PERIOD)),          \
   (id)[0]=(((id)[0])&0x0000000f)|(((v)&0x0fffffff)<<4))

   /* RAW type event */
#define EVENT_ID_MAX_NB_IN_RUN   0xffffffff
#define EVENT_ID_MAX_BURST_NB    0x00000fff
#define EVENT_ID_MAX_NB_IN_BURST 0x000fffff
#define LOAD_RAW_EVENT_ID(id,nir,bn,nib)                     \
  (assert(((bn)>=0)  && ((bn)<=EVENT_ID_MAX_BURST_NB) &&     \
	  ((nib)>=0) && ((nib)<=EVENT_ID_MAX_NB_IN_BURST)),  \
   (id)[0]=nir,                                              \
   (id)[1]=((bn<<20)&0xfff00000)|(nib&0x000fffff))
#define EVENT_ID_SET_NB_IN_RUN(id,nir)                  \
  (assert(((nir)>=0)&&((nir)<=EVENT_ID_MAX_NB_IN_RUN)), \
   (id)[0]=(nir))
#define EVENT_ID_SET_BURST_NB(id,bn)                    \
  (assert(((bn)>=0)&&((bn)<=EVENT_ID_MAX_BURST_NB)),    \
   (id)[1]=((id)[1]&0x000fffff)|(((bn)<<20)&0xfff00000))
#define EVENT_ID_SET_NB_IN_BURST(id,nib)                  \
  (assert(((nib)>=0)&&((nib)<=EVENT_ID_MAX_NB_IN_BURST)), \
   (id)[1]=((id)[1]&0xfff00000)|((nib)&0x000fffff))
#define EVENT_ID_GET_NB_IN_RUN(id)   ((id)[0])
#define EVENT_ID_GET_BURST_NB(id)    (((id)[1]>>20)&0x00000fff)
#define EVENT_ID_GET_NB_IN_BURST(id) ((id)[1]&0x000fffff)

   /* All events */
#define EQ_EVENT_ID(a,b) ((((a)[0])==((b)[0]))&&(((a)[1])==((b)[1])))
#define GT_EVENT_ID(a,b) \
    (((((a)[0])>((b)[0])))||((((a)[0])==((b)[0]))&&(((a)[1])>((b)[1]))))
#define LT_EVENT_ID(a,b) \
    ((((a)[0])<((b)[0])) || ((((a)[0])==((b)[0]))&&(((a)[1])<((b)[1]))))
#define GE_EVENT_ID(a,b) (!LT_EVENT_ID(a,b))
#define LE_EVENT_ID(a,b) (!GT_EVENT_ID(a,b))
#define COPY_EVENT_ID(from,to) \
                      memcpy((void*)to,(const void*)from,EVENT_ID_BYTES)
#define ADD_EVENT_ID(a,b) ((a)[1]+=(b)[1],(a)[0]+=(b)[0])
#define SUB_EVENT_ID(a,b) ((a)[1]-=(b)[1],(a)[0]-=(b)[0])
#define ZERO_EVENT_ID(id) memset(id,0,EVENT_ID_BYTES)

/* ---------- Trigger pattern (and relative masks) ---------- */
#define EVENT_TRIGGER_PATTERN_BYTES    8
#define EVENT_TRIGGER_PATTERN_WORDS    ((EVENT_TRIGGER_PATTERN_BYTES)>>2)
typedef unsigned long32 eventTriggerPatternType[EVENT_TRIGGER_PATTERN_WORDS];
#define EVENT_TRIGGER_ID_MIN           1
#define EVENT_TRIGGER_ID_MAX           50
#define CHECK_TRIGGER(t)               (assert(((t)>=EVENT_TRIGGER_ID_MIN) && \
                                               ((t)<=EVENT_TRIGGER_ID_MAX)))
#define TRIGGER_TO_BIT(t)              (1<<((t)&0x1f))
#define TRIGGER_TO_WORD(t)             (CHECK_TRIGGER(t), (t)>>5)
#define ZERO_TRIGGER_PATTERN(p)        memset(p,0,EVENT_TRIGGER_PATTERN_BYTES)
#define SET_TRIGGER_IN_PATTERN(p,id)   (p)[TRIGGER_TO_WORD(id)] |= \
                                                            TRIGGER_TO_BIT(id)
#define CLEAR_TRIGGER_IN_PATTERN(p,id) (p)[TRIGGER_TO_WORD(id)] &= \
                                                         ~(TRIGGER_TO_BIT(id))
#define FLIP_TRIGGER_IN_PATTERN(p,id)  (p)[TRIGGER_TO_WORD(id)] ^= \
                                                            TRIGGER_TO_BIT(id)
#define TEST_TRIGGER_IN_PATTERN(p,id)  (((p)[TRIGGER_TO_WORD(id)] & \
                                                     TRIGGER_TO_BIT(id)) != 0)
#define TRIGGER_PATTERN_INVALID(p)     (((p)[0] & 1) == 0)
#define TRIGGER_PATTERN_VALID(p)       (((p)[0] & 1) != 0)
#define VALIDATE_TRIGGER_PATTERN(p)    ((p)[0] |= 1)
#define INVALIDATE_TRIGGER_PATTERN(p)  ((p)[0] &= 0xfffffffe)
#define COPY_TRIGGER_PATTERN(f,t)      memcpy(t,f,EVENT_TRIGGER_PATTERN_BYTES)
#define TRIGGER_PATTERN_OK(p)          (((p)[1] & 0xfff80000) == 0)

/* ---------- Detectors cluster (and relative masks) ---------- */
#define EVENT_DETECTOR_PATTERN_BYTES 4
#define EVENT_DETECTOR_PATTERN_WORDS (EVENT_DETECTOR_PATTERN_BYTES>>2)
typedef unsigned long32 eventDetectorPatternType[EVENT_DETECTOR_PATTERN_WORDS];
#define EVENT_DETECTOR_ID_MIN  1
#define EVENT_DETECTOR_ID_MAX 24
#define CHECK_DETECTOR(d) (assert(((d) >= EVENT_DETECTOR_ID_MIN) &&\
                                  ((d) <= EVENT_DETECTOR_ID_MAX)))
#define DETECTOR_TO_BIT(d)              (CHECK_DETECTOR(d), 1<<(d))
#define ZERO_DETECTOR_PATTERN(p)        ((p)[0] = 0)
#define SET_DETECTOR_IN_PATTERN(p,d)    ((p)[0] |= DETECTOR_TO_BIT(d))
#define CLEAR_DETECTOR_IN_PATTERN(p,d)  ((p)[0] &= ~(DETECTOR_TO_BIT(d)))
#define FLIP_DETECTOR_IN_PATTERN(p,d)   ((p)[0] ^= DETECTOR_TO_BIT(d))
#define TEST_DETECTOR_IN_PATTERN(p,d)   (((p)[0] & DETECTOR_TO_BIT(d))!=0)
#define DETECTOR_PATTERN_INVALID(p)     (((p)[0] & 1) == 0)
#define DETECTOR_PATTERN_VALID(p)       (((p)[0] & 1) != 0)
#define VALIDATE_DETECTOR_PATTERN(p)    ((p)[0] |= 1)
#define INVALIDATE_DETECTOR_PATTERN(p)  ((p)[0] &= 0xfffffffe)
#define COPY_DETECTOR_PATTERN(f,t)      ((t)[0] = (f)[0])
#define DETECTOR_PATTERN_OK(p)          (((p)[0] & 0xfe000000) == 0)

/* ---------- The sizes and positions of the typeAttribute field ---------- */
#define ALL_ATTRIBUTE_WORDS    3
#define ALL_ATTRIBUTE_BYTES    (ALL_ATTRIBUTE_WORDS * 4)
#define ALL_ATTRIBUTE_BITS     (ALL_ATTRIBUTE_BYTES * 8)
#define USER_ATTRIBUTE_WORDS   2
#define USER_ATTRIBUTE_BYTES   (USER_ATTRIBUTE_WORDS * 4)
#define USER_ATTRIBUTE_BITS    (USER_ATTRIBUTE_BYTES * 8)
#define FIRST_USER_ATTRIBUTE   0
#define LAST_USER_ATTRIBUTE    (USER_ATTRIBUTE_BITS - 1)
#define SYSTEM_ATTRIBUTE_WORDS 1
#define SYSTEM_ATTRIBUTE_BYTES (SYSTEM_ATTRIBUTE_WORDS * 4)
#define SYSTEM_ATTRIBUTE_BITS  (SYSTEM_ATTRIBUTE_BYTES * 8)
#define FIRST_SYSTEM_ATTRIBUTE USER_ATTRIBUTE_BITS
#define LAST_SYSTEM_ATTRIBUTE  (USER_ATTRIBUTE_BITS + \
                                 SYSTEM_ATTRIBUTE_BITS - 1)
typedef unsigned long32 eventTypeAttributeType[ALL_ATTRIBUTE_WORDS];

   /* Word and bit definitions */
#define SYS_ATTR_2_W(b) (assert(((b)>=64)&&((b)<=95)),2)
#define USR_ATTR_2_W(b) (assert(((b)>= 0)&&((b)<=63)),(b)>>5)
#define ATTR_2_W(b)     (assert(((b)>= 0)&&((b)<=95)),(b)>>5)
#define ATTR_2_B(b)     (1<<((b)&0x1f))

   /* Macros to handle all attributes without distinction */
#define RESET_ATTRIBUTES(m)      ((m)[2] = (m)[1] = (m)[0] = 0)
#define SET_ANY_ATTRIBUTE(m,b)   (m)[ATTR_2_W(b)] |=  ATTR_2_B(b)
#define CLEAR_ANY_ATTRIBUTE(m,b) (m)[ATTR_2_W(b)] &= ~(ATTR_2_B(b))
#define FLIP_ANY_ATTRIBUTE(m,b)  (m)[ATTR_2_W(b)] ^=  ATTR_2_B(b)
#define TEST_ANY_ATTRIBUTE(m,b)  (((m)[ATTR_2_W(b)] & ATTR_2_B(b))!= 0)
#define COPY_ALL_ATTRIBUTES( from, to ) \
      memcpy((void *)&to[0], (const void *)&from[0], ALL_ATTRIBUTE_BYTES)
#define OR_ALL_ATTRIBUTES( from, to ) \
      ((to)[2] |= (from)[2], (to)[1] |= (from)[1], (to)[0] |= (from)[0])

   /* Macros to handle SYSTEM attributes */
#define RESET_SYSTEM_ATTRIBUTES(m)  ((m)[2] = 0)
#define SET_SYSTEM_ATTRIBUTE(m,b)   (m)[SYS_ATTR_2_W(b)] |= ATTR_2_B(b)
#define CLEAR_SYSTEM_ATTRIBUTE(m,b) (m)[SYS_ATTR_2_W(b)] &= ~(ATTR_2_B(b))
#define FLIP_SYSTEM_ATTRIBUTE(m,b)  (m)[SYS_ATTR_2_W(b)] ^= ATTR_2_B(b)
#define TEST_SYSTEM_ATTRIBUTE(m,b)  (((m)[SYS_ATTR_2_W(b)] & ATTR_2_B(b)) != 0)
#define COPY_SYSTEM_ATTRIBUTES( from, to ) \
   memcpy((void *)&to[2], (const void *)&from[2], SYSTEM_ATTRIBUTE_BYTES)
#define OR_SYSTEM_ATTRIBUTES( from, to ) ((to)[2] |= (from)[2])

   /* Macros to handle USER attributes */
#define RESET_USER_ATTRIBUTES(m)  ((m)[0] = (m)[1] = 0)
#define SET_USER_ATTRIBUTE(m,b)   (m)[USR_ATTR_2_W(b)] |= ATTR_2_B(b)
#define CLEAR_USER_ATTRIBUTE(m,b) (m)[USR_ATTR_2_W(b)] &= ~(ATTR_2_B(b))
#define FLIP_USER_ATTRIBUTE(m,b)  (m)[USR_ATTR_2_W(b)] ^= ATTR_2_B(b)
#define TEST_USER_ATTRIBUTE(m,b)  (((m)[USR_ATTR_2_W(b)] & ATTR_2_B(b)) != 0)
#define COPY_USER_ATTRIBUTES( from, to ) \
     memcpy((void *)&to[0], (const void *)&from[0], USER_ATTRIBUTE_BYTES)
#define OR_USER_ATTRIBUTES( from, to ) \
      ((to)[1] |= (from)[1], (to)[0] |= (from)[0])

   /* System attributes assignment */
#define ATTR_P_START              64          /* Start of a phase          */
#define ATTR_P_END                65          /* End of a phase            */
#define ATTR_START_OF_RUN_START   ATTR_P_START/* Start of SOR phase        */
#define ATTR_START_OF_RUN_END     ATTR_P_END  /* End of SOR phase          */
#define ATTR_END_OF_RUN_START     ATTR_P_START/* Start of EOR phase        */
#define ATTR_END_OF_RUN_END       ATTR_P_END  /* End of SOR phase          */
#define ATTR_EVENT_SWAPPED        66          /* Swapped event header      */
#define ATTR_EVENT_PAGED     	  67          /* Paged event               */
#define ATTR_SUPER_EVENT          68          /* Super event               */
#define ATTR_ORBIT_BC             69          /* Orbit/bunch crossing in ID*/
#define ATTR_KEEP_PAGES           70          /* Do not deallocate pages   */
#define ATTR_HLT_DECISION	  71	      /* Event contains HLT decis. */

#define ATTR_EVENT_DATA_TRUNCATED 94          /* Truncated payload         */
#define ATTR_EVENT_ERROR          95          /* Invalid event content     */

#define SYSTEM_ATTRIBUTES_OK(m) \
  ((((m)[2]) & ~(ATTR_2_B(ATTR_P_START)              | \
                 ATTR_2_B(ATTR_P_END)                | \
                 ATTR_2_B(ATTR_EVENT_SWAPPED)        | \
                 ATTR_2_B(ATTR_EVENT_PAGED)          | \
                 ATTR_2_B(ATTR_SUPER_EVENT)          | \
                 ATTR_2_B(ATTR_ORBIT_BC)             | \
                 ATTR_2_B(ATTR_KEEP_PAGES)           | \
                 ATTR_2_B(ATTR_HLT_DECISION)         | \
                 ATTR_2_B(ATTR_EVENT_DATA_TRUNCATED) | \
                 ATTR_2_B(ATTR_EVENT_ERROR))) == 0)

/* ---------- LDC and GDC identifier ---------- */
typedef unsigned long32 eventHostIdType;
typedef eventHostIdType eventLdcIdType;
typedef eventHostIdType eventGdcIdType;
#define HOST_ID_MIN ((eventHostIdType)0)         /* The minimum allowed ID */
#define HOST_ID_MAX ((eventHostIdType)511)       /* The maximum allowed ID */
#define VOID_ID     ((eventHostIdType)-1)        /* Unloaded ID            */

/* ---------- timestamp ---------- */
/* The following definition is in common for 32 and 64 bit machines.
   In both architectures, the field must be loaded into a time_t
   variable before being used. Failure to do so may cause undefined
   results up to the early termination of the process.

   The recommended procedure to use this field is the following:

   #include <time.h>

   time_t t;

   t = eventHeaderStruct.eventTimestamp;
   cTime( &t ); (or whatever else can be done with a time_t)

   Please note that the available timestamp will wrap sometime
   around Jan 18, 19:14:07, 2038...
*/
typedef unsigned long32 eventTimestampType;

/* ---------- The event header structure (with + without data) ---------- */
struct eventHeaderStruct { 
  eventSizeType             eventSize;
  eventMagicType            eventMagic;
  eventHeadSizeType         eventHeadSize;
  eventVersionType          eventVersion;
  eventTypeType             eventType;
  eventRunNbType            eventRunNb;
  eventIdType               eventId;
  eventTriggerPatternType   eventTriggerPattern;
  eventDetectorPatternType  eventDetectorPattern;
  eventTypeAttributeType    eventTypeAttribute;
  eventLdcIdType            eventLdcId;
  eventGdcIdType            eventGdcId;
  eventTimestampType        eventTimestamp;
};

struct eventStruct {
    struct eventHeaderStruct eventHeader;
           unsigned short    eventRawData[1];
};

/* ========== Definitions for the Vector ========== */
typedef short        eventVectorBankIdType;
typedef unsigned int eventVectorSizeType;
typedef unsigned int eventVectorOffsetType;

struct eventVectorStruct {
  eventVectorBankIdType eventVectorBankId;
  unsigned              eventVectorPointsToVector : 1;
  eventVectorSizeType   eventVectorSize;
  eventVectorOffsetType eventVectorStartOffset;
};

/* ========== Definitions for the payload descriptor ========== */
typedef        unsigned long32   eventNumEquipmentsType;
typedef struct eventVectorStruct eventExtensionVectorType;

struct vectorPayloadDescriptorStruct {
  eventNumEquipmentsType   eventNumEquipments;
  eventExtensionVectorType eventExtensionVector;
};

/* ========== Definitions for the equipment header ========== */
typedef long32                 equipmentSizeType;
typedef long32                 equipmentTypeType;
typedef long32                 equipmentIdType;
typedef eventTypeAttributeType equipmentTypeAttributeType;
typedef long32                 equipmentBasicElementSizeType;

struct equipmentHeaderStruct {
  equipmentSizeType             equipmentSize;
  equipmentTypeType             equipmentType;
  equipmentIdType               equipmentId;
  equipmentTypeAttributeType    equipmentTypeAttribute;
  equipmentBasicElementSizeType equipmentBasicElementSize;
};

struct equipmentDescriptorStruct {
  struct equipmentHeaderStruct equipmentHeader;
  struct eventVectorStruct     equipmentVector;
};

struct equipmentStruct {
  struct equipmentHeaderStruct equipmentHeader;
         unsigned short        equipmentRawData[1];
};

/* ========== Definition of the event location for the simpleFifo ========== */
struct eventLocationDescriptorStruct {
  eventVectorBankIdType eventBankId;
  eventVectorOffsetType eventOffset;
};

/* ========== Global macros ========== */

/* The macro PAGED_EVENT_SIZE receives in input the ADDRESS of a paged
   event and returns the size (in bytes) of the first page of the
   event */
#define PAGED_EVENT_SIZE( event ) \
  (EVENT_HEAD_BASE_SIZE +sizeof( struct vectorPayloadDescriptorStruct ) + \
   ((*(eventNumEquipmentsType *)((void*)event+EVENT_HEAD_BASE_SIZE))* \
     (sizeof( struct equipmentDescriptorStruct ))))

/* ========== Common data header ========== */
#define CDH_SIZE (8 * 4)
#define CDH_VERSION 1

#define CDH_TRIGGER_OVERLAP_ERROR_BIT            0
#define CDH_TRIGGER_MISSING_ERROR_BIT            1
#define CDH_DATA_PARITY_ERROR_BIT                2
#define CDH_CONTROL_PARITY_ERROR_BIT             3
#define CDH_TRIGGER_INFORMATION_UNAVAILABLE_BIT  4
#define CDH_FEE_ERROR_BIT                        5
#define CDH_HLT_DECISION_BIT                     6
#define CDH_HLT_PAYLOAD_BIT                      7
#define CDH_HLT_DDG_PAYLOAD_BIT                  8

/* Please note how the above data structure has been
   defined for LE systems. Code running on BE systems
   must have all fields reverted to work correctly! */
struct commonDataHeaderStruct {
  unsigned cdhBlockLength               : 32;
  /* ------------------------------------- */
  unsigned cdhEventId1                  : 12;
  unsigned cdhMBZ0                      :  2;
  unsigned cdhL1TriggerMessage          : 10;
  unsigned cdhVersion                   :  8;
  /* ------------------------------------- */
  unsigned cdhEventId2                  : 24;
  unsigned cdhMBZ1                      :  8;
  /* ------------------------------------- */
  unsigned cdhParticipatingSubDetectors : 24;
  unsigned cdhBlockAttributes           :  8;
  /* ------------------------------------- */
  unsigned cdhMiniEventId               : 12;
  unsigned cdhStatusErrorBits           : 16;
  unsigned cdhMBZ2                      :  4;
  /* ------------------------------------- */
  unsigned cdhTriggerClassesLow         : 32;
  /* ------------------------------------- */
  unsigned cdhTriggerClassesHigh        : 18;
  unsigned cdhMBZ3                      : 10;
  unsigned cdhRoiLow                    :  4;
  /* ------------------------------------- */
  unsigned cdhRoiHigh                   : 32;
};

#endif
