// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

//  @file   AliHLTExternalInterface.cxx
//  @author Matthias Richter, Timm Steinbeck
//  @date   
//  @brief  Pure C interface to the AliRoot HLT analysis framework
//          The C interface is used from external applications like the HLT
//          online framework

#include "AliHLTExternalInterface.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTComponent.h"
#include "AliHLTSystem.h"
#include "AliHLTMisc.h"
#include <cerrno>
#include <cstring>
/////////////////////////////////////////////////////////////////////////////////////
//
// AliHLT external interface functions
//

static AliHLTComponentHandler *gComponentHandler = NULL;
static AliHLTRunDesc gRunDesc=kAliHLTVoidRunDesc;
static char* gRunType=NULL;

int AliHLTAnalysisInitSystem( unsigned long version, AliHLTAnalysisEnvironment* externalEnv, unsigned long runNo, const char* runType )
{
  // init the HLT system
  if ( gComponentHandler ) {
    return EINPROGRESS;
  }

  if (!externalEnv) {
    return EINVAL;
  }

  if (!externalEnv->fStructSize) {
    return EPERM;
  }

  // map the environment structure to the internal structure to
  // support different versions of the interface
  AliHLTAnalysisEnvironment mappedEnv;
  memset(&mappedEnv, 0, sizeof(mappedEnv));
  memcpy(&mappedEnv, externalEnv, sizeof(mappedEnv)<externalEnv->fStructSize?sizeof(mappedEnv):externalEnv->fStructSize);
  mappedEnv.fStructSize=sizeof(mappedEnv);

  gComponentHandler = new AliHLTComponentHandler(&mappedEnv);
  if ( !gComponentHandler )
    return EFAULT;
  gComponentHandler->InitAliLogTrap(gComponentHandler);
  gComponentHandler->AnnounceVersion();

  if (version!=ALIHLT_DATA_TYPES_VERSION && externalEnv->fLoggingFunc) {
    const int messagelength=100;
    char message[messagelength+1];
    message[messagelength]=0;
    snprintf(message, messagelength, "interface definition does not match: internal %d - external %lu", ALIHLT_DATA_TYPES_VERSION, version);
    externalEnv->fLoggingFunc(externalEnv->fParam, kHLTLogWarning, "AliHLTAnalysisInitSystem", "interface version", message);
  }

  gRunDesc.fRunNo=runNo;
  if (runType && strlen(runType)>0) {
    int runtypelength=strlen(runType);
    gRunType=new char[runtypelength+1];
    if (gRunType) {
      gRunType[runtypelength]=0;
      strncpy(gRunType, runType, runtypelength);
    }
  }

  // the AliRoot dependent code is implemented by the
  // AliHLTMiscImplementation class in libHLTrec
  AliHLTMisc::Instance().InitCDB(getenv("ALIHLT_HCDBDIR"), getenv("ALIHLT_HCDBSNAPSHOT"));
  AliHLTMisc::Instance().SetCDBRunNo(gRunDesc.fRunNo);
  AliHLTMisc::Instance().InitMagneticField();

  return 0;
}

int AliHLTAnalysisDeinitSystem()
{
  // De-init the HLT system and clean-up internal memory
  if (gComponentHandler) delete gComponentHandler;
  gComponentHandler = NULL;

  if (gRunType) delete[] gRunType;
  gRunType=NULL;

  gRunDesc=kAliHLTVoidRunDesc;

  return 0;
}

int AliHLTAnalysisLoadLibrary( const char* libraryPath )
{
  // load a component library
  if ( !gComponentHandler )
    return ENXIO;
  return gComponentHandler->LoadLibrary( libraryPath );
}

int AliHLTAnalysisUnloadLibrary( const char* /*libraryPath*/ )
{
  // unload a component library
  if ( !gComponentHandler )
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
  //return gComponentHandler->UnloadLibrary( libraryPath );
  return 0;
}

int AliHLTAnalysisCreateComponent( const char* componentType, void* environParam, int argc, const char** argv, AliHLTComponentHandle* handle, const char* description )
{
  // create a component
  if ( !gComponentHandler ) return ENXIO;
  if (!handle) return EINVAL;

  AliHLTComponent* comp=NULL;

  int ret = gComponentHandler->CreateComponent( componentType, comp);
  if (ret>=0 && comp) {
    const AliHLTAnalysisEnvironment* comenv=gComponentHandler->GetEnvironment();
    comp->SetComponentEnvironment(comenv, environParam);
    if (comenv) {
      if (description) {
	comp->SetComponentDescription(description);
      }
      comp->SetRunDescription(&gRunDesc, gRunType);
    }
    ret=comp->Init(comenv, environParam, argc, argv);
  }
  *handle = reinterpret_cast<AliHLTComponentHandle>( comp );

  return ret;
}

int AliHLTAnalysisDestroyComponent( AliHLTComponentHandle handle )
{
  // destroy a component
  if ( !handle )
    return ENOENT;
  
  AliHLTComponent* pComp=reinterpret_cast<AliHLTComponent*>( handle );
  pComp->Deinit();
  delete pComp;
  return 0;
}

int AliHLTAnalysisProcessEvent( AliHLTComponentHandle handle, const AliHLTComponentEventData* evtData, const AliHLTComponentBlockData* blocks, 
                           AliHLTComponentTriggerData* trigData, AliHLTUInt8_t* outputPtr,
                           AliHLTUInt32_t* size, AliHLTUInt32_t* outputBlockCnt, 
                           AliHLTComponentBlockData** outputBlocks,
                           AliHLTComponentEventDoneData** edd )
{
  // process one event
  if ( !handle ) return EINVAL;
  AliHLTComponent* comp = reinterpret_cast<AliHLTComponent*>( handle );
  int ret=comp->ProcessEvent( *evtData, blocks, *trigData, outputPtr, *size, *outputBlockCnt, *outputBlocks, *edd );

  // internally, return values <0 are errors, >=0 are success with some
  // optional return vaue.
  // externally everthing !=0 is error. This is also according to the
  // common habit
  if (ret>0) ret=0;
  else if (ret<0) ret*=-1;
  return ret;
}

int AliHLTAnalysisGetOutputDataType( AliHLTComponentHandle handle, AliHLTComponentDataType* dataType )
{
  // get output data type of a component
  if ( !handle ) return EINVAL;
  AliHLTComponent* comp = reinterpret_cast<AliHLTComponent*>( handle );
  *dataType = comp->GetOutputDataType();
  return 0;
}

int AliHLTAnalysisGetOutputSize( AliHLTComponentHandle handle, unsigned long* constEventBase, unsigned long* constBlockBase, double* inputBlockMultiplier )
{
  // get output data size of a component
  if ( !handle ) return EINVAL;
  AliHLTComponent* comp = reinterpret_cast<AliHLTComponent*>( handle );
  // TODO: extend component interface
  if (constBlockBase) *constBlockBase=0;
  comp->GetOutputDataSize( *constEventBase, *inputBlockMultiplier );
  return 0;
}

struct AliHLTAnalysisInterfaceCall {
  const char* fSignature;
  void* fCall;
};

const AliHLTAnalysisInterfaceCall gAliHLTInterfaceCallSignatures[]={
  //int AliHLTAnalysisInitSystem( unsigned long version, AliHLTAnalysisEnvironment* externalEnv, unsigned long runNo, const char* runType )
  {"int AliHLTAnalysisInitSystem(unsigned long,AliHLTAnalysisEnvironment*,unsigned long,const char*)", (void*)AliHLTAnalysisInitSystem},

  //int AliHLTAnalysisDeinitSystem()
  {"int AliHLTAnalysisDeinitSystem()", (void*)AliHLTAnalysisDeinitSystem},

  //int AliHLTAnalysisLoadLibrary( const char* libraryPath )
  {"int AliHLTAnalysisLoadLibrary(const char*)", (void*)AliHLTAnalysisLoadLibrary},

  //int AliHLTAnalysisUnloadLibrary( const char* /*libraryPath*/ )
  {"int AliHLTAnalysisUnloadLibrary(const char*)", (void*)AliHLTAnalysisUnloadLibrary}, 

  //int AliHLTAnalysisCreateComponent( const char* componentType, void* environParam, int argc, const char** argv, AliHLTComponentHandle* handle, const char* description )
  {"int AliHLTAnalysisCreateComponent(const char*,void*,int,const char**,AliHLTComponentHandle*,const char*)", (void*)AliHLTAnalysisCreateComponent},

  //int AliHLTAnalysisDestroyComponent( AliHLTComponentHandle handle )
  {"int AliHLTAnalysisDestroyComponent(AliHLTComponentHandle)", (void*)AliHLTAnalysisDestroyComponent}, 

  //int AliHLTAnalysisProcessEvent( AliHLTComponentHandle handle, const AliHLTComponentEventData* evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData* trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t* size, AliHLTUInt32_t* outputBlockCnt, AliHLTComponentBlockData** outputBlocks, AliHLTComponentEventDoneData** edd )
  {"int AliHLTAnalysisProcessEvent(AliHLTComponentHandle,const AliHLTComponentEventData*,const AliHLTComponentBlockData*,AliHLTComponentTriggerData*,AliHLTUInt8_t*,AliHLTUInt32_t*,AliHLTUInt32_t*,AliHLTComponentBlockData**,AliHLTComponentEventDoneData**)", (void*)AliHLTAnalysisProcessEvent},

  //int AliHLTAnalysisGetOutputDataType( AliHLTComponentHandle handle, AliHLTComponentDataType* dataType )
  {"int AliHLTAnalysisGetOutputDataType(AliHLTComponentHandle,AliHLTComponentDataType*)", (void*)AliHLTAnalysisGetOutputDataType},

  //int AliHLTAnalysisGetOutputSize( AliHLTComponentHandle handle, unsigned long* constEventBase, unsigned long* constBlockBase, double* inputBlockMultiplier )
  {"int AliHLTAnalysisGetOutputSize(AliHLTComponentHandle,unsigned long*,unsigned long*,double*)", (void*)AliHLTAnalysisGetOutputSize},

  {NULL, NULL}
};

void* AliHLTAnalysisGetInterfaceCall(const char* signature)
{
  // get function pointer of an interface call
  // the function pointer is selected according to the requested function signature
  if (!signature) return NULL;
  for (int i=0; gAliHLTInterfaceCallSignatures[i].fSignature!=NULL; i++) {
    if (strcmp(gAliHLTInterfaceCallSignatures[i].fSignature, signature)==0) {
      return gAliHLTInterfaceCallSignatures[i].fCall;
    }
  }
  return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////
//
// AliHLTSystem interface functions
//

int AliHLTSystemSetOptions(AliHLTSystem* pInstance, const char* options)
{
  // deprecated helper function
  int iResult=0;
  if (pInstance) {
    AliHLTSystem* pSystem=reinterpret_cast<AliHLTSystem*>(pInstance);
    iResult=pSystem->ScanOptions(options);
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTSystemProcessHLTOUT(AliHLTSystem* pInstance, AliHLTOUT* pHLTOUT, AliESDEvent* esd)
{
  // deprecated helper function
  int iResult=0;
  if (pInstance) {
    AliHLTSystem* pSystem=reinterpret_cast<AliHLTSystem*>(pInstance);
    iResult=pSystem->ProcessHLTOUT(pHLTOUT, esd);
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}
