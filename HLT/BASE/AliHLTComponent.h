// @(#) $Id$

#ifndef ALIHLTCOMPONENT_H
#define ALIHLTCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTComponent
   base class for HLT components
 */

#include <errno.h>
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTDefinitions.h"
#include "TObject.h"

class AliHLTComponentHandler;

class AliHLTComponent : public AliHLTLogging {
 public:
  AliHLTComponent();
  virtual ~AliHLTComponent();

  enum TComponentType { kUnknown=0, kSource=1, kProcessor=2, kSink=3 };
  virtual int Init( AliHLTComponentEnvironment* environ, void* environ_param, int argc, const char** argv );
  virtual int Deinit();
  virtual int ProcessEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
			    AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
			    AliHLTUInt32_t& size, AliHLTUInt32_t& outputBlockCnt, 
			    AliHLTComponent_BlockData*& outputBlocks,
			    AliHLTComponent_EventDoneData*& edd ) = 0;

  // Information member functions for registration.
  virtual TComponentType GetComponentType() = 0; // Source, sink, or processor
  virtual const char* GetComponentID() = 0;
  virtual void GetInputDataTypes( vector<AliHLTComponent_DataType>& ) = 0;
  virtual AliHLTComponent_DataType GetOutputDataType() = 0;
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) = 0;

  // Spawn function, return new class instance
  virtual AliHLTComponent* Spawn() = 0;
 
  static int SetGlobalComponentHandler(AliHLTComponentHandler* pCH, int bOverwrite=0) {
    int iResult=0;
    if (fpComponentHandler==NULL || bOverwrite!=0)
      fpComponentHandler=pCH;
    else
      iResult=-EPERM;
    return iResult;
  }
  static int UnsetGlobalComponentHandler() {
    return SetGlobalComponentHandler(NULL,1);
  }
 protected:

  // Fill various structures with default values.
  void FillBlockData( AliHLTComponent_BlockData& blockData ) {
    blockData.fStructSize = sizeof(blockData);
    FillShmData( blockData.fShmKey );
    blockData.fOffset = ~(AliHLTUInt32_t)0;
    blockData.fPtr = NULL;
    blockData.fSize = 0;
    FillDataType( blockData.fDataType );
    blockData.fSpecification = ~(AliHLTUInt32_t)0;
  }
  void FillShmData( AliHLTComponent_ShmData& shmData ) {
    shmData.fStructSize = sizeof(shmData);
    shmData.fShmType = gkAliHLTComponent_InvalidShmType;
    shmData.fShmID = gkAliHLTComponent_InvalidShmID;
  }
  void FillDataType( AliHLTComponent_DataType& dataType ) {
    dataType.fStructSize = sizeof(dataType);
    memset( dataType.fID, '*', 8 );
    memset( dataType.fOrigin, '*', 4 );
  }
  
  virtual int DoInit( int argc, const char** argv ){
    return 0;
  }

  virtual int DoDeinit(){
    return 0;
  }

  void* AllocMemory( unsigned long size ) {
    if (fEnvironment.fAllocMemoryFunc)
      return (*fEnvironment.fAllocMemoryFunc)(fEnvironment.fParam, size );
    return NULL;
  }

  int MakeOutputDataBlockList( const vector<AliHLTComponent_BlockData>& blocks, AliHLTUInt32_t* blockCount,
			       AliHLTComponent_BlockData** outputBlocks ) {
    if (fEnvironment.fMakeOutputDataBlockListFunc)
      return (*fEnvironment.fMakeOutputDataBlockListFunc)(fEnvironment.fParam, blocks, blockCount, outputBlocks );
    return -ENOSYS;
  }

  int GetEventDoneData( unsigned long size, AliHLTComponent_EventDoneData** edd ) {
    if (fEnvironment.fGetEventDoneDataFunc)
      return (*fEnvironment.fGetEventDoneDataFunc)(fEnvironment.fParam, fCurrentEvent, size, edd );
    return -ENOSYS;
  }


 private:
  static AliHLTComponentHandler* fpComponentHandler;
  AliHLTComponentEnvironment fEnvironment;

  AliHLTEventID_t fCurrentEvent; // Set by ProcessEvent before actual processing starts (e.g. before calling AliHLTProcessor::DoEvent)

  ClassDef(AliHLTComponent, 0)
};
#endif
