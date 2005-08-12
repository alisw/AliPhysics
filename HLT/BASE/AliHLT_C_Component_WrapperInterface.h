// @(#) $Id$

#ifndef ALIHLT_C_COMPONENT_WARAPPERINTERFACE_H
#define ALIHLT_C_COMPONENT_WARAPPERINTERFACE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLT_C_Component_WrapperInterface
   pure C interface to the AliRoot HLT component handler
   utilized by the HLT PubSub Framework
 */

#include <AliHLTDataTypes.h>

#ifdef __cplusplus
extern "C" {
#endif
typedef void* AliHLTComponentHandle;

  const AliHLTComponentHandle kEmptyHLTComponentHandle = 0;

int AliHLT_C_Component_InitSystem( AliHLTComponentEnvironment* environ );
int AliHLT_C_Component_DeinitSystem();
int AliHLT_C_Component_LoadLibrary( const char* libraryPath );
int AliHLT_C_Component_UnloadLibrary( const char* libraryPath );
int AliHLT_C_CreateComponent( const char* componentType, void* environ_param, int argc, const char** argv, AliHLTComponentHandle* handle );
void AliHLT_C_DestroyComponent( AliHLTComponentHandle );
int AliHLT_C_ProcessEvent( AliHLTComponentHandle, const AliHLTComponent_EventData* evtData, const AliHLTComponent_BlockData* blocks, 
                           AliHLTComponent_TriggerData* trigData, AliHLTUInt8_t* outputPtr,
                           AliHLTUInt32_t* size, AliHLTUInt32_t* outputBlockCnt, 
                           AliHLTComponent_BlockData** outputBlocks,
                           AliHLTComponent_EventDoneData** edd );
int AliHLT_C_GetOutputDataType( AliHLTComponentHandle, AliHLTComponent_DataType* dataType );
int AliHLT_C_GetOutputSize( AliHLTComponentHandle, unsigned long* constBase, double* inputMultiplier );

#ifdef __cplusplus
}
#endif

#endif //ALIHLT_C_COMPONENT_WARAPPERINTERFACE_H 
