// @(#) $Id$

#ifndef ALIHLTDATATYPES_H
#define ALIHLTDATATYPES_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTDataTypes.h
    @author Matthias Richter, Timm Steinbeck, Jochen Thaeder
    @date   
    @brief  Data type declaration for the HLT module.
*/

//////////////////////////////////////////////////////////////////////////
//
// version no of HLT data types
//
//////////////////////////////////////////////////////////////////////////

/* Version   Description
 *   1       first version until June 07; implicite, not tagged
 *   2       introduced June 07, enhanced/cleaned/arranged structure
 *   3       2007-11-15 RAW DDL data type added; some inconsistencies fixed
 *           ('void' and 'any' origins); added signed HLT basic data types
 *           2007-11-23 origin defines have become variables in conjunction
 *           to be used with the operator| (AliHLTComponentDataType)
 *           2007-11-24 added trigger structs and ESD tree data type
 *   4       Component configuration event added
 */
#define ALIHLT_DATA_TYPES_VERSION 4

//////////////////////////////////////////////////////////////////////////
//
// HLT data origin variables.
//
// By converting from defines to variables, the origins can be used with
// the operator|
//
// AliHLTComponentDataType dt;
// dt = kAliHLTDataTypeDDLRaw | gkAliHLTDataOriginTPC;
//
//////////////////////////////////////////////////////////////////////////

/** field size of datat type origin */
const int kAliHLTComponentDataTypefOriginSize=4;


/** invalid data origin */
# define kAliHLTDataOriginVoid "\0\0\0"
/** old invalid data origin, kept for backward compatibility */
# define kAliHLTVoidDataOrigin "\0\0\0"

/** wildcard data type origin */
# define kAliHLTDataOriginAny "***"
/** old wildcard data type origin, kept for backward compatibility */
# define kAliHLTAnyDataOrigin "***"

/** HLT out */
extern const char kAliHLTDataOriginOut[kAliHLTComponentDataTypefOriginSize];

/** HLT/PubSub private internal */
extern const char kAliHLTDataOriginPrivate[kAliHLTComponentDataTypefOriginSize];

/** TPC */
extern const char kAliHLTDataOriginTPC[kAliHLTComponentDataTypefOriginSize];

/** PHOS */
extern const char kAliHLTDataOriginPHOS[kAliHLTComponentDataTypefOriginSize];

/** MUON */
extern const char kAliHLTDataOriginMUON[kAliHLTComponentDataTypefOriginSize];

/** TRD */
extern const char kAliHLTDataOriginTRD[kAliHLTComponentDataTypefOriginSize];

/** ITS */
extern const char kAliHLTDataOriginITS[kAliHLTComponentDataTypefOriginSize];

//////////////////////////////////////////////////////////////////////////
//
// HLT common data type defines
//
//////////////////////////////////////////////////////////////////////////

/** field size of data type id */
const int kAliHLTComponentDataTypefIDsize=8;


/** invalid data type id */
# define kAliHLTVoidDataTypeID "\0\0\0\0\0\0\0"

/** special id for any data type id */
# define kAliHLTAnyDataTypeID "*******"

/** DDL RAW data */
# define kAliHLTDDLRawDataTypeID   {'D','D','L','_','R','A','W',' '}

/** calibration data for file exchange subscriber */
# define kAliHLTFXSCalibDataTypeID {'F','X','S','_','C','A','L',' '}

/** start of run (SOR) event 
 * @ref AliHLTRunDesc
 */
# define kAliHLTSORDataTypeID      {'S','T','A','R','T','O','F','R'}

/** end of run (EOR) event 
 * @ref AliHLTRunDesc
 */
# define kAliHLTEORDataTypeID      {'E','N','D','O','F','R','U','N'}

/** DDL list event 
 * @ref AliHLTEventDDL
 */
# define kAliHLTDDLDataTypeID      {'D','D','L','L','I','S','T',' '}

/** EventType event 
 * - empty payload, specification gives eventType
 */
# define kAliHLTEventDataTypeID    {'E','V','E','N','T','T','Y','P'}

/** ComponentConfiguration event
 * - payload contains the ID of the component in the analysis chain
 */
# define kAliHLTComConfDataTypeID  {'C','O','M','_','C','O','N','F'}

/** ESD data block
 * an AliESD object of varying origin
 * The 'V0' at the end allows a versioning
 */
# define kAliHLTESDObjectDataTypeID    {'A','L','I','E','S','D','V','0'}

/** ESD tree data block
 * TTree with an AliESD object of varying origin
 */
# define kAliHLTESDTreeDataTypeID      {'E','S','D','_','T','R','E','E'}

/** HW Address selection data block
 * - a selection list for 16 bit HW addresses
 * - varying origin
 */
# define kAliHLTHwAddr16DataTypeID     {'H','W','A','D','D','R','1','6'}

using namespace std;

extern "C" {
  //////////////////////////////////////////////////////////////////////////
  //
  // Basic HLT data types
  //
  //////////////////////////////////////////////////////////////////////////

  typedef unsigned char AliHLTUInt8_t;

  typedef signed char AliHLTInt8_t;

  typedef unsigned short AliHLTUInt16_t;

  typedef signed short AliHLTInt16_t;

  typedef unsigned int AliHLTUInt32_t;

  typedef signed int AliHLTInt32_t;

  typedef unsigned long long AliHLTUInt64_t;

  typedef signed long long AliHLTInt64_t;

  typedef float AliHLTFloat32_t;

  typedef double AliHLTFloat64_t;

  typedef AliHLTUInt64_t AliHLTEventID_t;

  //////////////////////////////////////////////////////////////////////////
  //
  // HLT logging levels
  //
  //////////////////////////////////////////////////////////////////////////

  enum AliHLTComponentLogSeverity { 
    kHLTLogNone      = 0,
    kHLTLogBenchmark = 0x1,
    kHLTLogDebug     = 0x2,
    kHLTLogInfo      = 0x4,
    kHLTLogWarning   = 0x8,
    kHLTLogError     = 0x10,
    kHLTLogFatal     = 0x20,
    kHLTLogAll       = 0x3f,
    kHLTLogDefault   = 0x3d 
};

  //////////////////////////////////////////////////////////////////////////
  //
  // HLT data structures for data exchange and external interface
  //
  //////////////////////////////////////////////////////////////////////////

  /**
   * @struct AliHLTComponentEventData
   * Event descriptor
   */
  struct AliHLTComponentEventData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTEventID_t fEventID;
    AliHLTUInt32_t fEventCreation_s;
    AliHLTUInt32_t fEventCreation_us;
    AliHLTUInt32_t fBlockCnt;
  };

  /**
   * @struct AliHLTComponentShmData
   * Shared memory descriptor.
   * Irrelevant for analysis components.
   */
  struct AliHLTComponentShmData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fShmType;
    AliHLTUInt64_t fShmID;
  };

  /**
   * @struct AliHLTComponentDataType
   * Data type descriptor for data blocks transferred through the processing
   * chain.
   */
  struct AliHLTComponentDataType
  {
    AliHLTUInt32_t fStructSize;
    char fID[kAliHLTComponentDataTypefIDsize];                      //!
    char fOrigin[kAliHLTComponentDataTypefOriginSize];              //!
  };

  /**
   * @struct AliHLTComponentBlockData
   * This is the decription of data blocks exchanged between components.
   * \b IMPORTANT: The validity of fPtr and fOffset is different for input and
   * output blocks:
   * - input blocks: The \em fPtr member always points to the beginning of the data
   *                 of size \em fSize. fOffset is ignored and should be in most
   *                 case 0.
   * - output blocks: The \em fPtr member is ignored by the framework. \em fOffset
   *                  must specify the start of the data relative to the output
   *                  buffer. The data block has size \em fSize.
   */
  struct AliHLTComponentBlockData
  {
    /* size and version of the struct */
    AliHLTUInt32_t fStructSize;
    /* shared memory key, ignored by processing components */
    AliHLTComponentShmData fShmKey;
    /* offset of output data relative to the output buffer */
    AliHLTUInt32_t fOffset;
    /* start of the data for input data blocks, fOffset to be ignored*/
    void* fPtr;
    /* size of the data block */
    AliHLTUInt32_t fSize;
    /* data type of the data block */
    AliHLTComponentDataType fDataType;
    /* data specification of the data block */
    AliHLTUInt32_t fSpecification;
  };

  /**
   * @struct AliHLTComponentEventDoneData
   * 
   */
  struct AliHLTComponentEventDoneData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fDataSize;
    void* fData;
  };

  /**
   * @struct AliHLTRunDesc
   * Event descriptor.
   * The struct is send with the SOR and EOR events.
   */
  struct AliHLTRunDesc
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fRunNo;
    AliHLTUInt32_t fRunType;
  };

  //////////////////////////////////////////////////////////////////////////
  //
  // Trigger meta information
  //
  //////////////////////////////////////////////////////////////////////////

  /** field size of fAttribute */
  static const int gkAliHLTBlockDAttributeCount = 8;

  /** field size of fCommonHeader */
  static const int gkAliHLTCommonHeaderCount = 8;

  /** size of the DDL list */
  static const int gkAliHLTDDLListSize = 30;

  /**
   * @struct AliHLTEventDDL
   * DDL list event.
   * The struct is send with the DDLLIST event.
   * Used in the trigger structure for internal apperance of 
   * the DLLs as well as for the HLT readout list send to DAQ 
   * ( as DataType : kAliHLTDataTypeDDL )
   */
  struct AliHLTEventDDL
  {
    AliHLTUInt32_t fCount;
    AliHLTUInt32_t fList[gkAliHLTDDLListSize];
  };

  /**
   * @struct AliHLTEventTriggerData
   */
  struct AliHLTEventTriggerData
  {
    AliHLTUInt8_t  fAttributes[gkAliHLTBlockDAttributeCount]; 
    AliHLTUInt64_t fHLTStatus; // Bit field 
    AliHLTUInt32_t fCommonHeaderWordCnt;
    AliHLTUInt32_t fCommonHeader[gkAliHLTCommonHeaderCount]; 
    AliHLTEventDDL fReadoutList;
  };

  /**
   * @struct AliHLTComponentTriggerData
   * Trigger data
   */
  struct AliHLTComponentTriggerData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fDataSize;
    void* fData;
  };

  //////////////////////////////////////////////////////////////////////////
  //
  // HLT Event Type Specification
  //
  //////////////////////////////////////////////////////////////////////////

  /** Unknown eventType specification */
  static const AliHLTUInt32_t gkAliEventTypeUnknown = ~(AliHLTUInt32_t)0;
  /** SOR eventType specification */ 
  static const AliHLTUInt32_t gkAliEventTypeStartOfRun=1;
  /** Data eventType specification */
  static const AliHLTUInt32_t gkAliEventTypeData=2;
  /** EOR eventType specification */ 
  static const AliHLTUInt32_t gkAliEventTypeEndOfRun=4;
  /** Corrupt eventType specification */
  static const AliHLTUInt32_t gkAliEventTypeCorruptID=8;
  /** Calibration eventType specification */ 
  static const AliHLTUInt32_t gkAliEventTypeCalibration=16;
  /** DataReplay eventType specification */
  static const AliHLTUInt32_t gkAliEventTypeDataReplay=32;
  /** Tick eventType specification */ 
  static const AliHLTUInt32_t gkAliEventTypeTick=64;
  /** Max eventType specification */ 
  static const AliHLTUInt32_t gkAliEventTypeMax=64;

  //////////////////////////////////////////////////////////////////////////
  //
  // HLT defines and defaults
  //
  //////////////////////////////////////////////////////////////////////////

  /** invalid event id */
  const AliHLTEventID_t kAliHLTVoidEventID=~(AliHLTEventID_t)0;

  /** invalid data specification */
  const AliHLTUInt32_t kAliHLTVoidDataSpec = ~(AliHLTUInt32_t)0;

  /** invalid shared memory type */
  const AliHLTUInt32_t gkAliHLTComponentInvalidShmType = 0;

  /** invalid shared memory id */
  const AliHLTUInt64_t gkAliHLTComponentInvalidShmID = ~(AliHLTUInt64_t)0;

  /** invalid data type */
  const AliHLTComponentDataType kAliHLTVoidDataType = {
    sizeof(AliHLTComponentDataType),
    kAliHLTVoidDataTypeID,
    kAliHLTDataOriginVoid
  };

  // there is currently a problem with rootcint if the predefined ids
  // (commented below) are used. rootcint does not find the id if they
  // are char arrays defined with {} and individual chars. If strings
  // are used it works fine
  /** any data type */
  const AliHLTComponentDataType kAliHLTAnyDataType = {
    sizeof(AliHLTComponentDataType),
    kAliHLTAnyDataTypeID,
    kAliHLTDataOriginAny
  };

  /** multiple output data types */
  extern const AliHLTComponentDataType kAliHLTMultipleDataType;

  /** data to file exchange subscriber */
  extern const AliHLTComponentDataType kAliHLTDataTypeFXSCalib;

  /** DDL list data type */
  extern const AliHLTComponentDataType kAliHLTDataTypeDDL;

  /** SOR data type */
  extern const AliHLTComponentDataType kAliHLTDataTypeSOR;

  /** EOR data type */
  extern const AliHLTComponentDataType kAliHLTDataTypeEOR;

  /** Event type specification */
  extern const AliHLTComponentDataType kAliHLTDataTypeEvent;

  /** Configuration event data type */
  extern const AliHLTComponentDataType kAliHLTDataTypeComConf;

  /** RAW DDL data specification, origin is 'any', data publisher origin correctly */
  extern const AliHLTComponentDataType kAliHLTDataTypeDDLRaw;

  /** ESD object data specification, origin is 'any' */
  extern const AliHLTComponentDataType kAliHLTDataTypeESDObject;

  /** ESD Tree data specification, origin is 'any' */
  extern const AliHLTComponentDataType kAliHLTDataTypeESDTree;

  /** 16 bit Hardware address selection data specification, origin is 'any' */
  extern const AliHLTComponentDataType kAliHLTDataTypeHwAddr16;

  //////////////////////////////////////////////////////////////////////////
  //
  // FXS subscriber meta information
  //
  //////////////////////////////////////////////////////////////////////////

  static const int gkAliHLTFXSHeaderfOriginSize = 4;
  static const int gkAliHLTFXSHeaderfFileIDSize = 128;
  static const int gkAliHLTFXSHeaderfDDLNumberSize = 64;

  /** Header in front of the data payload, in order to sent data to the FXS. */
  struct AliHLTFXSHeader
  {
    AliHLTUInt32_t fHeaderVersion;
    AliHLTUInt32_t fRunNumber;
    char fOrigin[gkAliHLTFXSHeaderfOriginSize];
    char fFileID[gkAliHLTFXSHeaderfFileIDSize];
    char fDDLNumber[gkAliHLTFXSHeaderfDDLNumberSize];
  };  

  //////////////////////////////////////////////////////////////////////////
  //
  // Component running environment
  //
  //////////////////////////////////////////////////////////////////////////

  /** logging function */
  typedef int (*AliHLTfctLogging)( void* param, 
				   AliHLTComponentLogSeverity severity,
				   const char* origin,
				   const char* keyword,
				   const char* message);

  /**
   * @struct AliHLTComponentEnvironment
   * Running environment for analysis components.
   * The struct describes function callbacks for 
   */
  struct AliHLTComponentEnvironment
  {
    AliHLTUInt32_t fStructSize;
    void* fParam;
    void* (*fAllocMemoryFunc)( void* param, unsigned long size );
#if 0
    // future addition already foreseen/envisioned
    // IMPORTANT: don not just remove the defines as this breaks the binary
    // compatibility
    int (*fAllocShmMemoryFunc)( void* param, unsigned long size, AliHLTComponentBlockData* blockLocation );
#endif
    int (*fGetEventDoneDataFunc)( void* param, AliHLTEventID_t eventID, unsigned long size, AliHLTComponentEventDoneData** edd );
    AliHLTfctLogging fLoggingFunc;
  };
}

//////////////////////////////////////////////////////////////////////////
//
// Data type helper functions
//
//////////////////////////////////////////////////////////////////////////

inline bool operator==( const AliHLTComponentDataType& dt1, const AliHLTComponentDataType& dt2 )
{
  bool any1=true, any2=true, void1=true, void2=true, match=true;
  for ( int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++ ) {
    any1&=(dt1.fOrigin[i]==kAliHLTDataOriginAny[i]);
    any2&=(dt2.fOrigin[i]==kAliHLTDataOriginAny[i]);
    void1&=(dt1.fOrigin[i]==kAliHLTDataOriginVoid[i]);
    void2&=(dt2.fOrigin[i]==kAliHLTDataOriginVoid[i]);
    match&=dt1.fOrigin[i]==dt2.fOrigin[i];
    if (!(match || (any2 && !void1) || (any1 && !void2)))
      return false;
  }

  any1=true, any2=true, match=true;
  for ( int i = 0; i < kAliHLTComponentDataTypefIDsize; i++ ) {
    any1&=(dt1.fID[i]==kAliHLTAnyDataTypeID[i]);
    any2&=(dt2.fID[i]==kAliHLTAnyDataTypeID[i]);
    void1&=(dt1.fID[i]==kAliHLTVoidDataTypeID[i]);
    void2&=(dt2.fID[i]==kAliHLTVoidDataTypeID[i]);
    match&=dt1.fID[i]==dt2.fID[i];
    if (!(match || (any2 && !void1) || (any1 && !void2)))
      return false;
  }
  return true;
}

inline bool operator!=( const AliHLTComponentDataType& dt1, const AliHLTComponentDataType& dt2 )
{
  return !(dt1==dt2);
}

inline AliHLTComponentDataType operator|(const AliHLTComponentDataType srcdt, const char origin[kAliHLTComponentDataTypefOriginSize])
{
  AliHLTComponentDataType dt=srcdt;
  for ( int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++ )
    dt.fOrigin[i]=origin[i];
  return dt;
}

#endif 
