// @(#) $Id$

#ifndef ALIHLT_C_COMPONENT_WARAPPERINTERFACE_H
#define ALIHLT_C_COMPONENT_WARAPPERINTERFACE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLT_C_Component_WrapperInterface.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Pure C interface to the AliRoot HLT component handler
    @note   Utilized by the HLT Online (PubSub) framework
*/

/** 
 * @defgroup alihlt_wrapper_interface The HLT wrapper interface
 * The wrapper interface is a pure C interface which allows to use the 
 * analysis components in external applications. The interface is utilized
 * to bind the analysis code to the PubSub framework. 
 *
 * \image html PubSub_WrapperComponent.png "Wrapper interface"
 *
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

#ifdef __cplusplus
}
#endif

#endif //ALIHLT_C_COMPONENT_WARAPPERINTERFACE_H 
