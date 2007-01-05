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

using namespace std;

extern "C" {

  typedef unsigned char AliHLTUInt8_t;

  typedef unsigned short AliHLTUInt16_t;

  typedef unsigned int AliHLTUInt32_t;

  typedef unsigned long long AliHLTUInt64_t;

  typedef AliHLTUInt64_t AliHLTEventID_t;

  enum AliHLTComponentLogSeverity { kHLTLogNone=0, kHLTLogBenchmark=1, kHLTLogDebug=2, kHLTLogInfo=4, kHLTLogWarning=8, kHLTLogError=16, kHLTLogFatal=32, kHLTLogAll=0x3f, kHLTLogDefault=0x39 };

  struct AliHLTComponentEventData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTEventID_t fEventID;
    AliHLTUInt32_t fEventCreation_s;
    AliHLTUInt32_t fEventCreation_us;
    AliHLTUInt32_t fBlockCnt;
  };

  struct AliHLTComponentShmData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fShmType;
    AliHLTUInt64_t fShmID;
  };

    const AliHLTUInt32_t gkAliHLTComponentInvalidShmType = 0;
    const AliHLTUInt64_t gkAliHLTComponentInvalidShmID = ~(AliHLTUInt64_t)0;

  const int kAliHLTComponentDataTypefIDsize=8;
  const int kAliHLTComponentDataTypefOriginSize=4;
  struct AliHLTComponentDataType
  {
    AliHLTUInt32_t fStructSize;
    char fID[kAliHLTComponentDataTypefIDsize];
    char fOrigin[kAliHLTComponentDataTypefOriginSize];
  };

# define kAliHLTVoidDataTypeID "\0\0\0\0\0\0\0"
# define kAliHLTVoidDataOrigin "\0\0\0"
  const AliHLTComponentDataType kAliHLTVoidDataType = {sizeof(AliHLTComponentDataType), kAliHLTVoidDataTypeID, kAliHLTVoidDataOrigin};
  
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

  struct AliHLTComponentTriggerData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fDataSize;
    void* fData;
  };
  
  struct AliHLTComponentEventDoneData
  {
    AliHLTUInt32_t fStructSize;
    AliHLTUInt32_t fDataSize;
    void* fData;
  };

  typedef int (*AliHLTfctLogging)( void* param, AliHLTComponentLogSeverity severity, const char* origin, const char* keyword, const char* message );

  struct AliHLTComponentEnvironment
  {
    AliHLTUInt32_t fStructSize;
    void* fParam;
    void* (*fAllocMemoryFunc)( void* param, unsigned long size );
#if 0
    int (*fAllocShmMemoryFunc)( void* param, unsigned long size, AliHLTComponentBlockData* blockLocation ); // future addition already foreseen/envisioned
#endif
    int (*fGetEventDoneDataFunc)( void* param, AliHLTEventID_t eventID, unsigned long size, AliHLTComponentEventDoneData** edd );
    AliHLTfctLogging fLoggingFunc;
  };
}

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
