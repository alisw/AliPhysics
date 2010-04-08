// $Id$

#ifndef ALIHLT_C_COMPONENT_WARAPPERINTERFACE_H
#define ALIHLT_C_COMPONENT_WARAPPERINTERFACE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLT_C_Component_WrapperInterface.h
//  @author Matthias Richter, Timm Steinbeck
//  @date   
//  @brief  Old C interface to the AliRoot HLT component handler
//  @note   This interface is deprecated, the new interface is defined
//          in HLT/BASE/AliHLTExternalInterface

/** 
 * @defgroup alihlt_wrapper_interface_deprecated First version of the HLT wrapper
 * interface.
 * This is the first version of the wrapper interface and became deprecated in
 * July 2008. See @ref alihlt_wrapper_interface for current interface
 * description. The interface is fixed and can not be extended.
 * @ingroup alihlt_wrapper_interface
 */

#include <AliHLTDataTypes.h>
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

#ifdef __cplusplus
extern "C" {
#endif

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
int AliHLT_C_Component_InitSystem( AliHLTComponentEnvironment* environ );

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
int AliHLT_C_Component_DeinitSystem();

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
int AliHLT_C_Component_LoadLibrary( const char* libraryPath );

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
int AliHLT_C_Component_UnloadLibrary( const char* libraryPath );

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
int AliHLT_C_CreateComponent( const char* componentType, void* environ_param, int argc, const char** argv, AliHLTComponentHandle* handle );

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
void AliHLT_C_DestroyComponent( AliHLTComponentHandle );

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
int AliHLT_C_SetRunDescription(const AliHLTRunDesc* desc, const char* runType);

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
int AliHLT_C_ProcessEvent( AliHLTComponentHandle handle, const AliHLTComponent_EventData* evtData, const AliHLTComponent_BlockData* blocks, 
                           AliHLTComponent_TriggerData* trigData, AliHLTUInt8_t* outputPtr,
                           AliHLTUInt32_t* size, AliHLTUInt32_t* outputBlockCnt, 
                           AliHLTComponent_BlockData** outputBlocks,
                           AliHLTComponent_EventDoneData** edd );

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
int AliHLT_C_GetOutputDataType( AliHLTComponentHandle, AliHLTComponent_DataType* dataType );

/**
 * 
 * @ingroup alihlt_wrapper_interface_deprecated
 */
int AliHLT_C_GetOutputSize( AliHLTComponentHandle, unsigned long* constBase, double* inputMultiplier );

#ifdef __cplusplus
}
#endif

#endif //ALIHLT_C_COMPONENT_WARAPPERINTERFACE_H 
