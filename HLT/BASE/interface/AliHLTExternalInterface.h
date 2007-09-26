// @(#) $Id$

#ifndef ALIHLT_EXTERNALINTERFACE_H
#define ALIHLT_EXTERNALINTERFACE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTExternalInterface.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Pure and dynamic C interface to the AliRoot HLT component handler
    @note   Utilized by the HLT Online (PubSub) framework
*/

#include <AliHLTDataTypes.h>
class AliHLTSystem;
#ifdef __cplusplus
extern "C" {
#endif

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
typedef void* AliHLTComponentHandle;

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
const AliHLTComponentHandle kEmptyHLTComponentHandle = 0;

/* Matthias Dec 2006
 * The names have been changed for Aliroot's coding conventions sake
 * The old names are defined for backward compatibility with the 
 * PublisherSubscriber framework
 */
typedef AliHLTComponentLogSeverity AliHLTComponent_LogSeverity;
typedef AliHLTComponentEventData AliHLTComponent_EventData;
typedef AliHLTComponentShmData AliHLTComponent_ShmData;
typedef AliHLTComponentDataType AliHLTComponent_DataType;
typedef AliHLTComponentBlockData AliHLTComponent_BlockData;
typedef AliHLTComponentTriggerData AliHLTComponent_TriggerData;
typedef AliHLTComponentEventDoneData AliHLTComponent_EventDoneData;
const AliHLTUInt32_t gkAliHLTComponent_InvalidShmType = gkAliHLTComponentInvalidShmType;
const AliHLTUInt64_t gkAliHLTComponent_InvalidShmID = gkAliHLTComponentInvalidShmID;

typedef int (*AliHLTExtFctInitSystem)( AliHLTComponentEnvironment* );

typedef int (*AliHLTExtFctDeinitSystem)();

typedef int (*AliHLTExtFctLoadLibrary)( const char* );

typedef int (*AliHLTExtFctUnloadLibrary)( const char* );

typedef int (*AliHLTExtFctCreateComponent)( const char*, void*, int, const char**, AliHLTComponentHandle* );

typedef void (*AliHLTExtFctDestroyComponent)( AliHLTComponentHandle );

typedef int (*AliHLTExtFctProcessEvent)( AliHLTComponentHandle, const AliHLTComponentEventData*, const AliHLTComponentBlockData*, 
					 AliHLTComponentTriggerData*, AliHLTUInt8_t*,
					 AliHLTUInt32_t*, AliHLTUInt32_t*, 
					 AliHLTComponentBlockData**,
					 AliHLTComponentEventDoneData** );

typedef int (*AliHLTExtFctGetOutputDataType)( AliHLTComponentHandle, AliHLTComponentDataType* );

typedef int (*AliHLTExtFctGetOutputSize)( AliHLTComponentHandle, unsigned long*, double* );

struct AliHLTExternalFuctions_t {
  AliHLTExtFctInitSystem        fctInitSystem;
  AliHLTExtFctDeinitSystem      fctDeinitSystem;
  AliHLTExtFctLoadLibrary       fctLoadLibrary;
  AliHLTExtFctUnloadLibrary     fctUnloadLibrary;
  AliHLTExtFctCreateComponent   fctCreateComponent;
  AliHLTExtFctDestroyComponent  fctDestroyComponent;
  AliHLTExtFctProcessEvent      fctProcessEvent;
  AliHLTExtFctGetOutputDataType fctGetOutputDataType;
  AliHLTExtFctGetOutputSize     fctGetOutputSize;
};

#define ALIHLT_FCT_ENTRY_INITSYSTEM        "AliHLT_C_Component_InitSystem"
#define ALIHLT_FCT_ENTRY_DEINITSYSTEM      "AliHLT_C_Component_DeinitSystem"
#define ALIHLT_FCT_ENTRY_LOADLIBRARY       "AliHLT_C_Component_LoadLibrary"
#define ALIHLT_FCT_ENTRY_UNLOADLIBRARY     "AliHLT_C_Component_UnloadLibrary"
#define ALIHLT_FCT_ENTRY_CREATECOMPONENT   "AliHLT_C_Component_CreateComponent"
#define ALIHLT_FCT_ENTRY_DESTROYCOMPONENT  "AliHLT_C_Component_DestroyComponent"
#define ALIHLT_FCT_ENTRY_PROCESSEVENT      "AliHLT_C_Component_ProcessEvent"
#define ALIHLT_FCT_ENTRY_GETOUTPUTDATATYPE "AliHLT_C_Component_GetOutputDataType"
#define ALIHLT_FCT_ENTRY_GETOUTPUTSIZE     "AliHLT_C_Component_GetOutputSize"

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
int AliHLT_C_Component_InitSystem( AliHLTComponentEnvironment* environ );

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
int AliHLT_C_Component_DeinitSystem();

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
int AliHLT_C_Component_LoadLibrary( const char* libraryPath );

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
int AliHLT_C_Component_UnloadLibrary( const char* libraryPath );

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
int AliHLT_C_CreateComponent( const char* componentType, void* environ_param, int argc, const char** argv, AliHLTComponentHandle* handle );

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
void AliHLT_C_DestroyComponent( AliHLTComponentHandle );

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
int AliHLT_C_ProcessEvent( AliHLTComponentHandle handle, const AliHLTComponent_EventData* evtData, const AliHLTComponent_BlockData* blocks, 
                           AliHLTComponent_TriggerData* trigData, AliHLTUInt8_t* outputPtr,
                           AliHLTUInt32_t* size, AliHLTUInt32_t* outputBlockCnt, 
                           AliHLTComponent_BlockData** outputBlocks,
                           AliHLTComponent_EventDoneData** edd );

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
int AliHLT_C_GetOutputDataType( AliHLTComponentHandle, AliHLTComponent_DataType* dataType );

/**
 * 
 * @ingroup alihlt_wrapper_interface
 */
int AliHLT_C_GetOutputSize( AliHLTComponentHandle, unsigned long* constBase, double* inputMultiplier );

/**
 * Set options for an AliHLTSystem instance.
 * The function is introduced for the sake of backward compatibility.
 * Called from AliHLTReconstructor, which loads the function dynamically.
 * @return neg. error code if failed                                     <br>
 *         -EFAULT    type cast failed                                   <br>
 *         -EINVAL    invalid parameter
 * @ingroup alihlt_system_interface
 */
int AliHLTSystemSetOptions(AliHLTSystem*, const char*);

#ifdef __cplusplus
}
#endif

#endif //ALIHLT_EXTERNALINTERFACE_H 
