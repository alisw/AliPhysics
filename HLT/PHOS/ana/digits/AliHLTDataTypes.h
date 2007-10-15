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
 */
#define ALIHLT_DATA_TYPES_VERSION 2

//////////////////////////////////////////////////////////////////////////
//
// HLT data origin defines
//
//////////////////////////////////////////////////////////////////////////

/** field size of datat type origin */
const int kAliHLTComponentDataTypefOriginSize=4;


/** invalid data origin */
# define kAliHLTVoidDataOrigin "\0\0\0"

/** special id for any data type origin */
# define kAliHLTAnyDataOrigin "***"

/** HLT out */
# define kAliHLTDataOriginOut     {'H','L','T',' '}

/** HLT/PubSub private internal */
# define kAliHLTDataOriginPrivate {'P','R','I','V'}

/** TPC */
# define kAliHLTDataOriginTPC     {'T','P','C',' '}

/** PHOS */
# define kAliHLTDataOriginPHOS    {'P','H','O','S'}

/** MUON */
# define kAliHLTDataOriginMUON    {'M','U','O','N'}

/** TRD */
# define kAliHLTDataOriginTRD     {'T','R','D',' '}

/** ITS */
# define kAliHLTDataOriginITS     {'I','T','S',' '}

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

using namespace std;

extern "C" {
  //////////////////////////////////////////////////////////////////////////
  //
  // Basic HLT data types
  //
  //////////////////////////////////////////////////////////////////////////

  typedef unsigned char AliHLTUInt8_t;

  typedef unsigned short AliHLTUInt16_t;

  typedef unsigned int AliHLTUInt32_t;

  typedef unsigned long long AliHLTUInt64_t;

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
   * Descriptor for data blocks.
   */
  struct AliHLTComponentBlockData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTComponentShmData fShmKey;
    AliHLTUInt32_t fOffset;
    void* fPtr;
    AliHLTUInt32_t fSize;
    AliHLTComponentDataType fDataType;
    AliHLTUInt32_t fSpecification;
  };

  /**
   * @struct AliHLTComponentTriggerData
   * Trigger data, not yet defined
   */
  struct AliHLTComponentTriggerData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fDataSize;
    void* fData;
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

  /** size of the DDL list */
  static const int gkAliHLTDDLListSize = 29;

  /**
   * @struct AliHLTEventDDL
   * DDL list event.
   * The struct is send with the DDLLIST event.
   */
  struct AliHLTEventDDL
  {
    AliHLTUInt32_t fCount;
    AliHLTUInt32_t fList[gkAliHLTDDLListSize];
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
    kAliHLTVoidDataOrigin
  };

  // there is currently a problem with rootcint if the predefined ids
  // (commented below) are used. rootcint does not find the id if they
  // are char arrays defined with {} and individual chars. If strings
  // are used it works fine
  /** any data type */
  const AliHLTComponentDataType kAliHLTAnyDataType = {
    sizeof(AliHLTComponentDataType),
    kAliHLTAnyDataTypeID,
    kAliHLTAnyDataOrigin
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
    for ( int i = 0; i < kAliHLTComponentDataTypefIDsize; i++ )
	if ( dt1.fID[i] != dt2.fID[i] )
	    return false;
    for ( int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++ )
	if ( dt1.fOrigin[i] != dt2.fOrigin[i] )
	    return false;
    return true;
    }

inline bool operator!=( const AliHLTComponentDataType& dt1, const AliHLTComponentDataType& dt2 )
    {
    for ( int i = 0; i < kAliHLTComponentDataTypefIDsize; i++ )
	if ( dt1.fID[i] != dt2.fID[i] )
	    return true;
    for ( int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++ )
	if ( dt1.fOrigin[i] != dt2.fOrigin[i] )
	    return true;
    return false;
    }


#endif 
