// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

//  @file   AliHLT_C_Component_WrapperInterface.cxx
//  @author Matthias Richter, Timm Steinbeck
//  @date   
//  @brief  Old C interface to the AliRoot HLT component handler
//  @note   This interface is deprecated, the new interface is defined
//          in HLT/BASE/AliHLTExternalInterface

#include "AliHLT_C_Component_WrapperInterface.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTComponent.h"
#include "AliHLTMisc.h"
#include <errno.h>

using namespace std;

static AliHLTComponentHandler *gComponentHandler_C = NULL;
static AliHLTRunDesc gRunDesc=kAliHLTVoidRunDesc;
static char* gRunType=NULL;

int AliHLT_C_Component_InitSystem( AliHLTComponentEnvironment* comenv )
{
  // init the HLT system
  if ( gComponentHandler_C )
    {
      return EINPROGRESS;
    }

  // July 2008  
  // Due to a bug in the SimpleComponentWrapper and AliRootWrapperSubscriber
  // the fStructSize member was never initialized and we can not use this
  // method of synchronizing different versions.
  // This interface is now deprecated, only kept for backward compatibility.
  // All function pointers are explicitely mapped to the new structure.

  AliHLTAnalysisEnvironment mappedEnv;
  memset(&mappedEnv, 0, sizeof(mappedEnv));
  mappedEnv.fStructSize=sizeof(mappedEnv);
  if (comenv) {
    mappedEnv.fParam               = comenv->fParam;
    mappedEnv.fAllocMemoryFunc     = comenv->fAllocMemoryFunc;
    mappedEnv.fGetEventDoneDataFunc= comenv->fGetEventDoneDataFunc;
    mappedEnv.fLoggingFunc         = comenv->fLoggingFunc;
  }

  gComponentHandler_C = new AliHLTComponentHandler(&mappedEnv);
  if ( !gComponentHandler_C )
    return EFAULT;
  gComponentHandler_C->InitAliLogTrap(gComponentHandler_C);
  gComponentHandler_C->AnnounceVersion();
  return 0;
}

int AliHLT_C_Component_DeinitSystem()
{
  // De-init the HLT system and clean-up internal memory
  if ( gComponentHandler_C )
    {
      delete gComponentHandler_C;
      gComponentHandler_C = NULL;
    }
  return 0;
}

int AliHLT_C_Component_LoadLibrary( const char* libraryPath )
{
  // load a component library
  if ( !gComponentHandler_C )
    return ENXIO;
  return gComponentHandler_C->LoadLibrary( libraryPath );
}

int AliHLT_C_Component_UnloadLibrary( const char* /*libraryPath*/ )
{
  // unload a component library
  if ( !gComponentHandler_C )
    return ENXIO;
  // Matthias 26.10.2007
  // Unloading of libraries has to be re-worked. It has been commented out here
  // since the libraries will be unloaded at the destruction of the component
  // handler instance anyway. So it has no effect to the operation in PubSub.
  // With the introduction of the dynamic component registration via module
  // agents we run into trouble when cleaning up the samples managed by the
  // component handler. Destruction of the sample objects is done AFTER
  // unloading of the library and thus the destructor is not present any 
  // more.
  //return gComponentHandler_C->UnloadLibrary( libraryPath );
  return 0;
}

int AliHLT_C_CreateComponent( const char* componentType, void* environParam, int argc, const char** argv, AliHLTComponentHandle* handle )
{
  // create a component
  if ( !gComponentHandler_C )
    return ENXIO;
  if ( !handle ) return EINVAL;
  AliHLTComponent* comp=NULL;
  const char* cdbPath = getenv("ALIHLT_HCDBDIR");
  if (!cdbPath) cdbPath = getenv("ALICE_ROOT");
  int ret = gComponentHandler_C->CreateComponent( componentType, comp);
  if (ret>=0 && comp) {
    AliHLTMisc::Instance().InitCDB(cdbPath, getenv("ALIHLT_HCDBSNAPSHOT"));
    AliHLTMisc::Instance().SetCDBRunNo(gRunDesc.fRunNo);
    comp->SetRunDescription(&gRunDesc, gRunType);
    const AliHLTAnalysisEnvironment* comenv=gComponentHandler_C->GetEnvironment();
    ret=comp->Init(comenv, environParam, argc, argv);
  }
  *handle = reinterpret_cast<AliHLTComponentHandle>( comp );

  return ret;
}

void AliHLT_C_DestroyComponent( AliHLTComponentHandle handle )
{
  // destroy a component
  if ( !handle )
    return;
  
  AliHLTComponent* pComp=reinterpret_cast<AliHLTComponent*>( handle );
  pComp->Deinit();
  delete pComp;
}

int AliHLT_C_SetRunDescription(const AliHLTRunDesc* desc, const char* runType)
{
  // set run description
  if (!desc) return -EINVAL;
  if (desc->fStructSize<sizeof(AliHLTUInt32_t)) return -EINVAL;
  if (!gComponentHandler_C) return ENXIO;

  memcpy(&gRunDesc, desc, desc->fStructSize<sizeof(gRunDesc)?desc->fStructSize:sizeof(gRunDesc));
  gRunDesc.fStructSize=sizeof(gRunDesc);
  if (gRunType) delete [] gRunType;
  gRunType=NULL;
  if (runType) {
    gRunType=new char[strlen(runType)+1];
    if (gRunType) strcpy(gRunType, runType);
  }
  return 0;
}

int AliHLT_C_ProcessEvent( AliHLTComponentHandle handle, const AliHLTComponentEventData* evtData, const AliHLTComponentBlockData* blocks, 
                           AliHLTComponentTriggerData* trigData, AliHLTUInt8_t* outputPtr,
                           AliHLTUInt32_t* size, AliHLTUInt32_t* outputBlockCnt, 
                           AliHLTComponentBlockData** outputBlocks,
                           AliHLTComponentEventDoneData** edd )
{
  // process one event
  if ( !handle )
    return ENXIO;
  AliHLTComponent* comp = reinterpret_cast<AliHLTComponent*>( handle );
  return comp->ProcessEvent( *evtData, blocks, *trigData, outputPtr, *size, *outputBlockCnt, *outputBlocks, *edd );
}

int AliHLT_C_GetOutputDataType( AliHLTComponentHandle handle, AliHLTComponentDataType* dataType )
{
  // get output data type of a component
  if ( !handle )
    return ENXIO;
  AliHLTComponent* comp = reinterpret_cast<AliHLTComponent*>( handle );
  *dataType = comp->GetOutputDataType();
  return 0;
}

int AliHLT_C_GetOutputSize( AliHLTComponentHandle handle, unsigned long* constBase, double* inputMultiplier )
{
  // get output data size of a component
  if ( !handle )
    return ENXIO;
  AliHLTComponent* comp = reinterpret_cast<AliHLTComponent*>( handle );
  comp->GetOutputDataSize( *constBase, *inputMultiplier );
  return 0;
}
