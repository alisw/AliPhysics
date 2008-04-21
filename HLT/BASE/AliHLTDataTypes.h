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
 *   4       Component configuration and DCS update events added
 *           gkAliHLTDDLListSize set from 29 to 30 according to new PubSub
 *           specs
 *   5       Data types for Run and Event summary, and for monitoring added
 *   6       Common data types for TreeD and TreeR defined
 */
#define ALIHLT_DATA_TYPES_VERSION 6

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

/** field size of datat type origin 
 * @ingroup alihlt_component_datatypes
 */
const int kAliHLTComponentDataTypefOriginSize=4;


/** invalid data origin 
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTDataOriginVoid "\0\0\0"
/** old invalid data origin, kept for backward compatibility */
# define kAliHLTVoidDataOrigin "\0\0\0"

/** wildcard data type origin 
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTDataOriginAny "***"
/** old wildcard data type origin, kept for backward compatibility */
# define kAliHLTAnyDataOrigin "***"

/** Data origin HLT out 
 * @ingroup alihlt_component_datatypes
 */
extern const char kAliHLTDataOriginOut[kAliHLTComponentDataTypefOriginSize];

/** Data origin HLT/PubSub private internal 
 * @ingroup alihlt_component_datatypes
 */
extern const char kAliHLTDataOriginPrivate[kAliHLTComponentDataTypefOriginSize];

/** Data origin TPC 
 * @ingroup alihlt_component_datatypes
 */
extern const char kAliHLTDataOriginTPC[kAliHLTComponentDataTypefOriginSize];

/** Data origin PHOS 
 * @ingroup alihlt_component_datatypes
 */
extern const char kAliHLTDataOriginPHOS[kAliHLTComponentDataTypefOriginSize];

/** Data origin MUON 
 * @ingroup alihlt_component_datatypes
 */
extern const char kAliHLTDataOriginMUON[kAliHLTComponentDataTypefOriginSize];

/** Data origin TRD 
 * @ingroup alihlt_component_datatypes
 */
extern const char kAliHLTDataOriginTRD[kAliHLTComponentDataTypefOriginSize];

/** Data origin ITS 
 * @ingroup alihlt_component_datatypes
 */
extern const char kAliHLTDataOriginITS[kAliHLTComponentDataTypefOriginSize];

//////////////////////////////////////////////////////////////////////////
//
// HLT common data type defines
//
//////////////////////////////////////////////////////////////////////////

/** field size of data type id 
 * @ingroup alihlt_component_datatypes
 */
const int kAliHLTComponentDataTypefIDsize=8;


/** invalid data type id 
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTVoidDataTypeID "\0\0\0\0\0\0\0"

/** special id for any data type id 
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTAnyDataTypeID "*******"

/** DDL RAW data 
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTDDLRawDataTypeID   {'D','D','L','_','R','A','W',' '}

/** calibration data for file exchange subscriber 
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTFXSCalibDataTypeID {'F','X','S','_','C','A','L',' '}

/** start of run (SOR) event 
 * @ref AliHLTRunDesc
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTSORDataTypeID      {'S','T','A','R','T','O','F','R'}

/** end of run (EOR) event 
 * @ref AliHLTRunDesc
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTEORDataTypeID      {'E','N','D','O','F','R','U','N'}

/** DDL list event 
 * @ref AliHLTEventDDL
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTDDLDataTypeID      {'D','D','L','L','I','S','T',' '}

/** EventType event 
 * - empty payload, specification gives eventType
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTEventDataTypeID    {'E','V','E','N','T','T','Y','P'}

/** ComponentConfiguration event
 * - payload contains the CDB path as string
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTComConfDataTypeID  {'C','O','M','_','C','O','N','F'}

/** DCS value update event
 * - payload contains string of relevant detectors
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTUpdtDCSDataTypeID  {'U','P','D','T','_','D','C','S'}

/** ESD data block
 * an AliESD object of varying origin
 * The 'V0' at the end allows a versioning
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTESDObjectDataTypeID    {'A','L','I','E','S','D','V','0'}

/** ESD tree data block
 * TTree with an AliESD object of varying origin
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTESDTreeDataTypeID      {'E','S','D','_','T','R','E','E'}

/** AliRoot TreeD
 * - the digits tree of an AliRoot module
 * @ingroup alihlt_component_datatypes
 */
#define kAliHLTTreeDDataTypeID         {'A','L','I','T','R','E','E','D'}

/** AliRoot TreeR
 * - the rec points tree of an AliRoot module
 * @ingroup alihlt_component_datatypes
 */
#define kAliHLTTreeRDataTypeID         {'A','L','I','T','R','E','E','D'}

/** HW Address selection data block
 * - a selection list for 16 bit HW addresses
 * - varying origin
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTHwAddr16DataTypeID     {'H','W','A','D','D','R','1','6'}

/** Event Statistics
 * - event statistics for given detectors
 * - varying origin
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTEventStatisticsDataTypeID     {'E','V','_','S','T','A','T','I'}

/** Event Summary
 * - event summary
 * - origin : kAliHLTDataOriginOut ( HLT )
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTEventSummaryDataTypeID        {'E','V','_','S','U','M','M','A'}

/** Run Statistics
 * - run statistics for given detectors
 * - varying origin
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTRunStatisticsDataTypeID       {'R','U','N','S','T','A','T','I'}

/** Run Summary
 * - run summary
 * - origin : kAliHLTDataOriginOut ( HLT )
 * @ingroup alihlt_component_datatypes
 */
# define kAliHLTRunSummaryDataTypeID          {'R','U','N','S','U','M','M','A'}

/** general ROOT TObject
 * - a general TObject exported from the HLT analysis
 * - varying origin
 * @ingroup alihlt_component_datatypes
 */
#define kAliHLTTObjectDataTypeID              {'R','O','O','T','T','O','B','J'}

/** ROOT TObjArray
 * - a TObjArray exported from the HLT analysis
 * - varying origin
 * @ingroup alihlt_component_datatypes
 */
#define kAliHLTTObjArrayDataTypeID            {'R','O','O','T','O','B','A','R'}

/** ROOT TTree
 * - a TTree object exported from the HLT analysis
 * - varying origin
 * @ingroup alihlt_component_datatypes
 */
#define kAliHLTTTreeDataTypeID                {'R','O','O','T','T','R','E','E'}

/** ROOT histogram
 * - a histogram object exported from the HLT analysis
 * - class derives from TH1 (directly or indirectly) and inherits all common functionality
 * - varying origin
 * @ingroup alihlt_component_datatypes
 */
#define kAliHLTHistogramDataTypeID            {'R','O','O','T','H','I','S','T'}

/** ROOT TNtuple
 * - a TNtupl object exported from the HLT analysis
 * - varying origin
 * @ingroup alihlt_component_datatypes
 */
#define kAliHLTTNtupleDataTypeID              {'R','O','O','T','T','U','P','L'}

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

  /**
   * Logging severities of the HLT
   */
  enum AliHLTComponentLogSeverity {
    /** no logging */
    kHLTLogNone      = 0,
    /** benchmark messages */
    kHLTLogBenchmark = 0x1,
    /** debug messages */
    kHLTLogDebug     = 0x2,
    /** info messages */
    kHLTLogInfo      = 0x4,
    /** warning messages */
    kHLTLogWarning   = 0x8,
    /** error messages */
    kHLTLogError     = 0x10,
    /** fatal error messages */
    kHLTLogFatal     = 0x20,
    /** few important messages not to be filtered out.
     * redirected to kHLTLogInfo in AliRoot
     */
    kHLTLogImportant = 0x40,
    /** special value to enable all messages */
    kHLTLogAll       = 0x7f,
    /** the default logging filter */
    kHLTLogDefault   = 0x79
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
   * @defgroup alihlt_component_datatypes Common Component Data Types
   * The analysis framework defines a number of common data types for
   * usage in the detector modules, like e.g. ::kAliHLTAnyDataType
   * and ::kAliHLTDataTypeDDLRaw. Those data types always have
   * origin ::kAliHLTDataOriginAny. The correct detector origin can be
   * set by using operator '|'
   * <pre>
   * AliHLTComponentDataType dt=kAliHLTDDLRawDataTypeID|kAliHLTDataOriginTPC
   * </pre>
   * @ingroup alihlt_component
   */

  /**
   * @struct AliHLTComponentDataType
   * Data type descriptor for data blocks transferred through the processing
   * chain.
   * @ingroup alihlt_component_datatypes
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
  const int gkAliHLTBlockDAttributeCount = 8;

  /** field size of fCommonHeader */
   const int gkAliHLTCommonHeaderCount = 8;

  /** size of the DDL list */
  const int gkAliHLTDDLListSize = 30;

  /** Number of Trigger Classes of CTP in CDH */
  const int gkNCTPTriggerClasses = 50;

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
  const AliHLTUInt32_t gkAliEventTypeUnknown = ~(AliHLTUInt32_t)0;
  /** SOR eventType specification */ 
  const AliHLTUInt32_t gkAliEventTypeStartOfRun=1;
  /** Data eventType specification */
  const AliHLTUInt32_t gkAliEventTypeData=2;
  /** EOR eventType specification */ 
  const AliHLTUInt32_t gkAliEventTypeEndOfRun=4;
  /** Corrupt eventType specification */
  const AliHLTUInt32_t gkAliEventTypeCorruptID=8;
  /** Calibration eventType specification */ 
  const AliHLTUInt32_t gkAliEventTypeCalibration=16;
  /** DataReplay eventType specification */
  const AliHLTUInt32_t gkAliEventTypeDataReplay=32;
  /** Configuration eventType specification */
  const AliHLTUInt32_t gkAliEventTypeConfiguration=34;
  /** Update DCS eventType specification */
  const AliHLTUInt32_t gkAliEventTypeReadPreprocessor=35;
  /** Tick eventType specification */ 
  const AliHLTUInt32_t gkAliEventTypeTick=64;
  /** Max eventType specification */ 
  const AliHLTUInt32_t gkAliEventTypeMax=64;

  //////////////////////////////////////////////////////////////////////////
  //
  // HLT defines and defaults
  //
  //////////////////////////////////////////////////////////////////////////

  /** invalid event id 
   * @ingroup alihlt_component_datatypes
   */
  const AliHLTEventID_t kAliHLTVoidEventID=~(AliHLTEventID_t)0;

  /** invalid data specification 
   * @ingroup alihlt_component_datatypes
   */
  const AliHLTUInt32_t kAliHLTVoidDataSpec = ~(AliHLTUInt32_t)0;

  /** invalid shared memory type */
  const AliHLTUInt32_t gkAliHLTComponentInvalidShmType = 0;

  /** invalid shared memory id */
  const AliHLTUInt64_t gkAliHLTComponentInvalidShmID = ~(AliHLTUInt64_t)0;

  /** invalid data type 
   * @ingroup alihlt_component_datatypes
   */
  const AliHLTComponentDataType kAliHLTVoidDataType = {
    sizeof(AliHLTComponentDataType),
    kAliHLTVoidDataTypeID,
    kAliHLTDataOriginVoid
  };

  // there is currently a problem with rootcint if the predefined ids
  // (commented below) are used. rootcint does not find the id if they
  // are char arrays defined with {} and individual chars. If strings
  // are used it works fine
  /** any data type 
   * @ingroup alihlt_component_datatypes
   */
  const AliHLTComponentDataType kAliHLTAnyDataType = {
    sizeof(AliHLTComponentDataType),
    kAliHLTAnyDataTypeID,
    kAliHLTDataOriginAny
  };

  /** multiple output data types 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTMultipleDataType;

  /** data to file exchange subscriber 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeFXSCalib;

  /** DDL list data type 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeDDL;

  /** SOR data type 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeSOR;

  /** EOR data type 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeEOR;

  /** Event type specification 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeEvent;

  /** Configuration event data type 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeComConf;

  /** DCS value update event 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeUpdtDCS;

  /** RAW DDL data specification, origin is 'any', data publisher origin correctly 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeDDLRaw;

  /** ESD object data specification, origin is 'any' 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeESDObject;

  /** ESD Tree data specification, origin is 'any' 
   
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeESDTree;

  /** AliRoot TreeD data specification, origin is 'any' 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeAliTreeD;

  /** AliRoot TreeR data specification, origin is 'any' 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeAliTreeR;

  /** 16 bit Hardware address selection data specification, origin is 'any' 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeHwAddr16;

  /** Event statistics 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeEventStatistics;

  /** Event summary 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeEventSummary;

  /** Event statistics 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeRunStatistics;

  /** Event summary 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeRunSummary;

  //////////////////////////////////////////////////////////////////////////
  //
  // Data Types for Monitoring objects
  //
  //////////////////////////////////////////////////////////////////////////

  /** general ROOT TObject 
   * @ingroup alihlt_component_datatypes
   */
  extern const AliHLTComponentDataType kAliHLTDataTypeTObject;            // {ROOTTOBJ,"***"}
									  		
  /** ROOT TObjArray 
   * @ingroup alihlt_component_datatypes
   */							  		
  extern const AliHLTComponentDataType kAliHLTDataTypeTObjArray;	  // {ROOTOBAR,"***"}
									  		
  /** ROOT TTree 
   * @ingroup alihlt_component_datatypes
   */							  		
  extern const AliHLTComponentDataType kAliHLTDataTypeTTree;		  // {ROOTTREE,"***"}
									  		
  /** ROOT TH1 (can be used for all histograms, they derive from TH1) 
   * @ingroup alihlt_component_datatypes
   */  		
  extern const AliHLTComponentDataType kAliHLTDataTypeHistogram;	  // {ROOTHIST,"***"}
									  		
  /** ROOT TNtuple 
   * @ingroup alihlt_component_datatypes
   */							  		
  extern const AliHLTComponentDataType kAliHLTDataTypeTNtuple;		  // {ROOTTUPL,"***"}

  //////////////////////////////////////////////////////////////////////////
  //
  // FXS subscriber meta information
  //
  //////////////////////////////////////////////////////////////////////////

  const int gkAliHLTFXSHeaderfOriginSize = 4;
  const int gkAliHLTFXSHeaderfFileIDSize = 128;
  const int gkAliHLTFXSHeaderfDDLNumberSize = 64;

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

/** Comparison operator for HLT component data types.
 * The operator takes wildcards into account, i.e. the ::kAliHLTAnyDataType,
 * ::kAliHLTAnyDataTypeID and ::kAliHLTDataOriginAny definitions.
 * @ingroup alihlt_component_datatypes
 */
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

/** Comparison operator for HLT component data types
 * invers of @ref operator==( const AliHLTComponentDataType& dt1, const AliHLTComponentDataType& dt2 )
 * @ingroup alihlt_component_datatypes
 */
inline bool operator!=( const AliHLTComponentDataType& dt1, const AliHLTComponentDataType& dt2 )
{
  return !(dt1==dt2);
}

/** exact comparison of HLT component data types
 * @ingroup alihlt_component_datatypes
 */
inline bool MatchExactly( const AliHLTComponentDataType& dt1, const AliHLTComponentDataType& dt2 )
{
  for ( int i = 0; i < kAliHLTComponentDataTypefIDsize; i++ )
    if ( dt1.fID[i] != dt2.fID[i] )
      return false;
  for ( int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++ )
    if ( dt1.fOrigin[i] != dt2.fOrigin[i] )
      return false;
  return true;
}

/** merge operator for HLT component data types and origins
 * @ingroup alihlt_component_datatypes
 */
inline AliHLTComponentDataType operator|(const AliHLTComponentDataType srcdt, const char origin[kAliHLTComponentDataTypefOriginSize])
{
  AliHLTComponentDataType dt=srcdt;
  for ( int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++ )
    dt.fOrigin[i]=origin[i];
  return dt;
}

#endif 
