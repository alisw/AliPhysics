// @(#) $Id$

#ifndef ALIHLTDATATYPES_H
#define ALIHLTDATATYPES_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTDataTypes.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Data type declaration for the HLT module.
*/

/* #include <vector> */
/* using namespace std; */

extern "C" {

  typedef unsigned char AliHLTUInt8_t;

  typedef unsigned short AliHLTUInt16_t;

  typedef unsigned int AliHLTUInt32_t;

  typedef unsigned long long AliHLTUInt64_t;

  typedef AliHLTUInt64_t AliHLTEventID_t;

  enum AliHLTComponent_LogSeverity { kHLTLogNone=0, kHLTLogBenchmark=1, kHLTLogDebug=2, kHLTLogInfo=4, kHLTLogWarning=8, kHLTLogError=16, kHLTLogFatal=32, kHLTLogAll=0x3f, kHLTLogDefault=0x39 };

  struct AliHLTComponent_EventData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTEventID_t fEventID;
    AliHLTUInt32_t fEventCreation_s;
    AliHLTUInt32_t fEventCreation_us;
    AliHLTUInt32_t fBlockCnt;
  };

  struct AliHLTComponent_ShmData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fShmType;
    AliHLTUInt64_t fShmID;
  };
    const AliHLTUInt32_t gkAliHLTComponent_InvalidShmType = 0;
    const AliHLTUInt64_t gkAliHLTComponent_InvalidShmID = ~(AliHLTUInt64_t)0;

  struct AliHLTComponent_DataType
  {
    AliHLTUInt32_t fStructSize;
    char fID[8];
    char fOrigin[4];
  };
  
  struct AliHLTComponent_BlockData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTComponent_ShmData fShmKey;
    AliHLTUInt32_t fOffset;
    void* fPtr;
    AliHLTUInt32_t fSize;
    AliHLTComponent_DataType fDataType;
    AliHLTUInt32_t fSpecification;
  };

  struct AliHLTComponent_TriggerData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fDataSize;
    void* fData;
  };
  
  struct AliHLTComponent_EventDoneData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fDataSize;
    void* fData;
  };

  typedef int (*AliHLTfctLogging)( void* param, AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* message );

  struct AliHLTComponentEnvironment
  {
    AliHLTUInt32_t fStructSize;
    void* fParam;
    void* (*fAllocMemoryFunc)( void* param, unsigned long size );
#if 0
    int (*fAllocShmMemoryFunc)( void* param, unsigned long size, AliHLTComponent_BlockData* blockLocation ); // future addition already foreseen/envisioned
#endif
/*     int (*fMakeOutputDataBlockListFunc)( void* param, const vector<AliHLTComponent_BlockData>& blocks, AliHLTUInt32_t* blockCount, AliHLTComponent_BlockData** outputBlocks ); */
    int (*fGetEventDoneDataFunc)( void* param, AliHLTEventID_t eventID, unsigned long size, AliHLTComponent_EventDoneData** edd );
    AliHLTfctLogging fLoggingFunc;
  };
}

inline bool operator==( const AliHLTComponent_DataType& dt1, const AliHLTComponent_DataType& dt2 )
    {
    for ( unsigned i = 0; i < 8; i++ )
	if ( dt1.fID[i] != dt2.fID[i] )
	    return false;
    for ( unsigned i = 0; i < 4; i++ )
	if ( dt1.fOrigin[i] != dt2.fOrigin[i] )
	    return false;
    return true;
    }

inline bool operator!=( const AliHLTComponent_DataType& dt1, const AliHLTComponent_DataType& dt2 )
    {
    for ( unsigned i = 0; i < 8; i++ )
	if ( dt1.fID[i] != dt2.fID[i] )
	    return true;
    for ( unsigned i = 0; i < 4; i++ )
	if ( dt1.fOrigin[i] != dt2.fOrigin[i] )
	    return true;
    return false;
    }

#endif 
