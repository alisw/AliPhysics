// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTComponentHandler.cxx
/// @author Matthias Richter, Timm Steinbeck
/// @date   
/// @brief  Implementation of HLT component handler.
///

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif
//#undef HAVE_DLFCN_H
#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#else
//#include <Riostream.h>
#include <TSystem.h>
#endif //HAVE_DLFCN_H
//#include "AliHLTStdIncludes.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTComponent.h"
#include "AliHLTDataTypes.h"
#include "AliHLTModuleAgent.h"
#include "TString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTComponentHandler)

AliHLTComponentHandler::AliHLTComponentHandler()
  :
  fComponentList(),
  fScheduleList(),
  fLibraryList(),
  fEnvironment(),
  fOwnedComponents(),
  fLibraryMode(kDynamic),
  fRunDesc(kAliHLTVoidRunDesc),
  fRunType(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  memset(&fEnvironment, 0, sizeof(AliHLTAnalysisEnvironment));
  fEnvironment.fStructSize=sizeof(AliHLTAnalysisEnvironment);
  AddStandardComponents();
}

AliHLTComponentHandler::AliHLTComponentHandler(AliHLTAnalysisEnvironment* pEnv)
  :
  AliHLTLogging(),
  fComponentList(),
  fScheduleList(),
  fLibraryList(),
  fEnvironment(),
  fOwnedComponents(),
  fLibraryMode(kDynamic),
  fRunDesc(kAliHLTVoidRunDesc),
  fRunType(NULL)
{
  // see header file for class documentation
  if (pEnv) {
    memcpy(&fEnvironment, pEnv, sizeof(AliHLTAnalysisEnvironment));
    if (pEnv->fLoggingFunc) {
      // the AliHLTLogging::Init method also sets the stream output
      // and notification handler to AliLog. This should only be done
      // if the logging environment contains a logging function
      // for redirection
      AliHLTLogging::Init(pEnv->fLoggingFunc);
    }
  }  else {
    memset(&fEnvironment, 0, sizeof(AliHLTAnalysisEnvironment));
    fEnvironment.fStructSize=sizeof(AliHLTAnalysisEnvironment);
  }
  //#ifndef __DEBUG
  //SetLocalLoggingLevel(kHLTLogError);
  //#else
  //SetLocalLoggingLevel(kHLTLogInfo);
  //#endif

  AddStandardComponents();
}

AliHLTComponentHandler::~AliHLTComponentHandler()
{
  // see header file for class documentation
  DeleteOwnedComponents();
  UnloadLibraries();
  if (fRunType) delete [] fRunType;
  fRunType=NULL;
}

AliHLTComponentHandler* AliHLTComponentHandler::fgpInstance=NULL;
int AliHLTComponentHandler::fgNofInstances=0;

AliHLTComponentHandler* AliHLTComponentHandler::CreateHandler()
{
  // see header file for class documentation
  if (!fgpInstance) fgpInstance=new AliHLTComponentHandler;
  fgNofInstances++;
  return fgpInstance;
}

int AliHLTComponentHandler::Destroy()
{
  // destroy/delete 'this', checks if 'this' pointer is the global instance,
  // reduce the instance counter and delete if there are no instances left
  // IMPORTANT: the object must be considered self-destroyed after the function
  int nofInstances=0;
  if (fgpInstance==this) {
    nofInstances=--fgNofInstances;
  }
  if (fgNofInstances==0) fgpInstance=NULL;
  if (nofInstances==0) delete this;
  return nofInstances;
}

int AliHLTComponentHandler::AnnounceVersion()
{
  // see header file for class documentation
  int iResult=0;
#ifdef PACKAGE_STRING
  extern void HLTbaseCompileInfo( const char*& date, const char*& time);
  const char* date="";
  const char* time="";
  HLTbaseCompileInfo(date, time);
  if (!date) date="unknown";
  if (!time) time="unknown";
  HLTImportant("%s build on %s (%s)", PACKAGE_STRING, date, time);
#else
  HLTImportant("ALICE High Level Trigger build on %s (%s) (embedded AliRoot build)", __DATE__, __TIME__);
#endif
  return iResult;
}

Int_t AliHLTComponentHandler::AddComponent(AliHLTComponent* pSample)
{
  // see header file for class documentation
  Int_t iResult=0;
  if (pSample==NULL) return -EINVAL;
  if ((iResult=RegisterComponent(pSample))>=0) {
    //HLTDebug("sample %s (%p) managed by handler", pSample->GetComponentID(), pSample);
    fOwnedComponents.push_back(pSample);
  }
  return iResult;
}

Int_t AliHLTComponentHandler::RegisterComponent(AliHLTComponent* pSample)
{
  // see header file for class documentation
  Int_t iResult=0;
  if (pSample) {
    if (FindComponent(pSample->GetComponentID())==NULL) {
      iResult=InsertComponent(pSample);
      if (iResult>=0) {
	HLTInfo("component %s registered", pSample->GetComponentID());
      }
    } else {
      // component already registered
      HLTDebug("component %s already registered, skipped", pSample->GetComponentID());
      iResult=-EEXIST;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponentHandler::DeregisterComponent( const char* componentID )
{
  // see header file for class documentation

  int iResult=0;
  if (componentID) {
    HLTWarning("not yet implemented, please notify the developers if you need this function");
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

Int_t AliHLTComponentHandler::ScheduleRegister(AliHLTComponent* pSample)
{
  // see header file for class documentation
  Int_t iResult=0;
  if (pSample) {
    fScheduleList.push_back(pSample);
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponentHandler::CreateComponent(const char* componentID, void* pEnvParam, int argc, const char** argv, AliHLTComponent*& component)
{
  // see header file for class documentation
  int iResult=CreateComponent(componentID, component);
  if (iResult>=0 && component) {
	HLTDebug("component \"%s\" created (%p)", componentID, component);
	if ((iResult=component->Init(&fEnvironment, pEnvParam, argc, argv))!=0) {
	  HLTError("Initialization of component \"%s\" failed with error %d", componentID, iResult);
	  delete component;
	  component=NULL;
	}
  }
  return iResult;
}

int AliHLTComponentHandler::CreateComponent(const char* componentID, AliHLTComponent*& component )
{
  // see header file for class documentation
  int iResult=0;
  if (componentID) {
    AliHLTComponent* pSample=FindComponent(componentID);
    if (pSample!=NULL) {
      component=pSample->Spawn();
      if (component) {
	HLTDebug("component \"%s\" created (%p)", componentID, component);
      } else {
	HLTError("can not spawn component \"%s\"", componentID);
	iResult=-ENOENT;
      }
    } else {
      HLTWarning("can not find component \"%s\"", componentID);
      iResult=-ENOENT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

Int_t AliHLTComponentHandler::FindComponentIndex(const char* componentID)
{
  // see header file for class documentation
  Int_t iResult=0;
  if (componentID) {
    AliHLTComponentPList::iterator element=fComponentList.begin();
    while (element!=fComponentList.end() && iResult>=0) {
      if (strcmp(componentID, (*element)->GetComponentID())==0) {
	break;
      }
      element++;
      iResult++;
    }
    if (element==fComponentList.end()) iResult=-ENOENT;
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

AliHLTComponent* AliHLTComponentHandler::FindComponent(const char* componentID)
{
  // see header file for class documentation
  AliHLTComponent* pSample=NULL;
  Int_t index=FindComponentIndex(componentID);
  if (index>=0) {
    pSample=(AliHLTComponent*)fComponentList.at(index);
  }
  return pSample;
}

Int_t AliHLTComponentHandler::InsertComponent(AliHLTComponent* pSample)
{
  // see header file for class documentation
  Int_t iResult=0;
  if (pSample!=NULL) {
    fComponentList.push_back(pSample);
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTComponentHandler::List() 
{
  // see header file for class documentation
  AliHLTComponentPList::iterator element=fComponentList.begin();
  int index=0;
  while (element!=fComponentList.end()) {
    HLTInfo("%d. %s", index++, (*element++)->GetComponentID());
  }
}

int AliHLTComponentHandler::HasOutputData( const char* componentID)
{
  // see header file for class documentation
  int iResult=0;
  AliHLTComponent* pSample=FindComponent(componentID);
  if (pSample) {
    AliHLTComponent::TComponentType ct=AliHLTComponent::kUnknown;
    ct=pSample->GetComponentType();
    iResult=(ct==AliHLTComponent::kSource || ct==AliHLTComponent::kProcessor);
  } else {
    iResult=-ENOENT;
  }
  return iResult;
}

void AliHLTComponentHandler::SetEnvironment(AliHLTAnalysisEnvironment* pEnv) 
{
  // see header file for class documentation
  if (pEnv) {
    memset(&fEnvironment, 0, sizeof(AliHLTAnalysisEnvironment));
    memcpy(&fEnvironment, pEnv, pEnv->fStructSize<sizeof(AliHLTAnalysisEnvironment)?pEnv->fStructSize:sizeof(AliHLTAnalysisEnvironment));
    fEnvironment.fStructSize=sizeof(AliHLTAnalysisEnvironment);
    if (fEnvironment.fLoggingFunc) {
      // the AliHLTLogging::Init method also sets the stream output
      // and notification handler to AliLog. This should only be done
      // if the logging environment contains a logging function
      // for redirection
      AliHLTLogging::Init(fEnvironment.fLoggingFunc);
    }
  }
}

AliHLTComponentHandler::TLibraryMode AliHLTComponentHandler::SetLibraryMode(TLibraryMode mode)
{
  // see header file for class documentation
  TLibraryMode old=fLibraryMode;
  fLibraryMode=mode;
  return old;
}

int AliHLTComponentHandler::LoadLibrary( const char* libraryPath, int bActivateAgents)
{
  // see header file for class documentation
  int iResult=0;
  if (libraryPath) {
    TString libName=libraryPath;
    int slash=-1;
    while ((slash=libName.Index("/"))>=0) {
      libName=libName.Remove(0, slash+1);
    }
    libName.ReplaceAll(".so","");

    // set the global component handler for static component registration
    AliHLTComponent::SetGlobalComponentHandler(this);

    AliHLTLibHandle hLib;
    AliHLTLibHandle* phSearch=FindLibrary(libraryPath);
    const char* loadtype="";
#ifdef HAVE_DLFCN_H
    // use interface to the dynamic linking loader

    // exeption does not help in Root context, the Root exeption
    // handler always catches the exeption before. Have to find out
    // how exeptions can be used in Root
    /*try*/ {
      hLib.fHandle=dlopen(libraryPath, RTLD_NOW);
      loadtype="dlopen";
    }
    /*
    catch (...) {
      // error message printed further down
      loadtype="dlopen exeption";
    }
    */
#else
    // use ROOT dynamic loader
    // check if the library was already loaded, as Load returns
    // 'failure' if the library was already loaded
    /*try*/ {
    if (phSearch) {
	int* pRootHandle=reinterpret_cast<int*>(phSearch->fHandle);
	(*pRootHandle)++;
	HLTDebug("instance %d of library %s loaded", (*pRootHandle), libraryPath);
	hLib.fHandle=pRootHandle;
    }
    
    if (hLib.fHandle==NULL && gSystem->Load(libraryPath)>=0) {
      int* pRootHandle=new int;
      if (pRootHandle) *pRootHandle=1;
      hLib.fHandle=pRootHandle;
      //HLTDebug("library %s loaded via gSystem", libraryPath);
    }
    loadtype="gSystem";
    }
    /*
    catch (...) {
      // error message printed further down
      loadtype="gSystem exeption";
    }
    */
#endif //HAVE_DLFCN_H
    if (hLib.fHandle!=NULL) {
      // create TString object to store library path and use pointer as handle 
      hLib.fName=new TString(libraryPath);
      hLib.fMode=fLibraryMode;
      fLibraryList.insert(fLibraryList.begin(), hLib);
      if (!phSearch) {
      typedef void (*CompileInfo)(const char*& date, const char*& time);
      CompileInfo fctInfo=(CompileInfo)FindSymbol(libraryPath, "CompileInfo");
      const char* date="";
      const char* time="";
      const char* buildOn="";
      if (fctInfo) {
	buildOn=" build on ";
	(*fctInfo)(date, time);
	if (!date) date="unknown";
	if (!time) time="unknown";
      }
      HLTImportant("using %s plugin%s%s %s (%s%s)", libraryPath, buildOn, date, time, hLib.fMode==kStatic?"persistent, ":"", loadtype);
      }

      // static registration of components when library is loaded
      iResult=RegisterScheduledComponents();

    } else {
      HLTError("can not load library %s (%s)", libraryPath, loadtype);
#ifdef HAVE_DLFCN_H
      HLTError("dlopen error: %s", dlerror());
#endif //HAVE_DLFCN_H
#ifdef __APPLE__
      iResult=-EFTYPE;
#else
      iResult=-ELIBACC;
#endif
    }
    AliHLTComponent::UnsetGlobalComponentHandler();
    
    if (iResult>=0) {
      // alternative dynamic registration by library agents
      // !!! has to be done after UnsetGlobalComponentHandler
      if (bActivateAgents) ActivateAgents(libName.Data());
    }

  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponentHandler::UnloadLibrary( const char* libraryPath )
{
  // see header file for class documentation
  int iResult=0;
  if (libraryPath) {
    vector<AliHLTLibHandle>::iterator element=fLibraryList.begin();
    while (element!=fLibraryList.end()) {
      TString* pName=reinterpret_cast<TString*>((*element).fName);
      if (pName->CompareTo(libraryPath)==0) {
	UnloadLibrary(*element);
	fLibraryList.erase(element);
	break;
      }
      element++;
  }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponentHandler::UnloadLibrary(AliHLTComponentHandler::AliHLTLibHandle &handle)
{
  // see header file for class documentation
  int iResult=0;
  fgAliLoggingFunc=NULL;
  TString* pName=reinterpret_cast<TString*>(handle.fName);
  if (handle.fMode!=kStatic) {
#ifdef HAVE_DLFCN_H
  // exeption does not help in Root context, the Root exeption
  // handler always catches the exeption before. Have to find out
  // how exeptions can be used in Root

  /*try*/ {
    dlclose(handle.fHandle);
  }
  /*
  catch (...) {
    HLTError("exeption caught during dlclose of library %s", pName!=NULL?pName->Data():"");
  }
  */
#else
  int* pCount=reinterpret_cast<int*>(handle.fHandle);
  if (--(*pCount)==0) {
    if (pName) {
      /** Matthias 26.04.2007
       * I spent about a week to investigate a bug which seems to be in ROOT.
       * Under certain circumstances, TSystem::Unload crashes. The crash occured
       * for the first time, when libAliHLTUtil was loaded from AliHLTSystem right
       * after the ComponentHandler was created. It does not occur when dlopen is
       * used. 
       * It has most likely to do with the garbage collection and automatic
       * cleanup in ROOT. The crash occurs when ROOT is terminated and before
       * an instance of AliHLTSystem was created.
       *   root [0] AliHLTSystem gHLT
       * It does not occur when the instance was created dynamically (but not even
       * deleted)
       *   root [0] AliHLTSystem* gHLT=new AliHLTSystem
       *
       * For that reason, the libraries are not unloaded here, even though there
       * will be memory leaks.
      gSystem->Unload(pName->Data());
       */
    }
    else {
      HLTError("missing library name, can not unload");
    }
    delete pCount;
  }
#endif //HAVE_DLFCN_H
  if (pName) {
    HLTDebug("unload library %s", pName->Data());
  } else {
    HLTWarning("missing name for unloaded library");
  }
  }
  handle.fName=NULL;
  handle.fHandle=NULL;
  if (pName) {
    delete pName;
  }
  pName=NULL;
  return iResult;
}

int AliHLTComponentHandler::UnloadLibraries()
{
  // see header file for class documentation
  int iResult=0;
  vector<AliHLTLibHandle>::iterator element=fLibraryList.begin();
  while (element!=fLibraryList.end()) {
    UnloadLibrary(*element);
    fLibraryList.erase(element);
    element=fLibraryList.begin();
  }
  return iResult;
}

AliHLTfctVoid AliHLTComponentHandler::FindSymbol(const char* library, const char* symbol)
{
  // see header file for class documentation
  AliHLTLibHandle* hLib=FindLibrary(library);
  if (hLib==NULL) return NULL;
  void (*pFunc)()=NULL;
#ifdef HAVE_DLFCN_H
  pFunc=(void (*)())dlsym(hLib->fHandle, symbol);
#else
  TString* name=reinterpret_cast<TString*>(hLib->fName);
  pFunc=(void (*)())gSystem->DynFindSymbol(name->Data(), symbol);
#endif
  return pFunc;
}

AliHLTComponentHandler::AliHLTLibHandle* AliHLTComponentHandler::FindLibrary(const char* library)
{
  // see header file for class documentation
  AliHLTLibHandle* hLib=NULL;
  vector<AliHLTLibHandle>::iterator element=fLibraryList.begin();
  while (element!=fLibraryList.end()) {
    TString* name=reinterpret_cast<TString*>((*element).fName);
    if (name->CompareTo(library)==0) {
      hLib=&(*element);
      break;
    }
    element++;
  }
  return hLib;
}

int AliHLTComponentHandler::AddStandardComponents()
{
  // see header file for class documentation
  int iResult=0;
  AliHLTComponent::SetGlobalComponentHandler(this);
  AliHLTComponent::UnsetGlobalComponentHandler();
  iResult=RegisterScheduledComponents();
  return iResult;
}

int AliHLTComponentHandler::RegisterScheduledComponents()
{
  // see header file for class documentation
  int iResult=0;
  AliHLTComponentPList::iterator element=fScheduleList.begin();
  int iLocalResult=0;
  while (element!=fScheduleList.end()) {
    iLocalResult=RegisterComponent(*element);
    if (iResult==0) iResult=iLocalResult;
    fScheduleList.erase(element);
    element=fScheduleList.begin();
  }
  return iResult;
}

int AliHLTComponentHandler::ActivateAgents(const char* library, const char* blackList)
{
  // see header file for class documentation
  vector<AliHLTModuleAgent*> agents;
  for (AliHLTModuleAgent* pAgent=AliHLTModuleAgent::GetFirstAgent(); 
       pAgent!=NULL;
       pAgent=AliHLTModuleAgent::GetNextAgent()) {

    // check if we found the agent for the specified library
    if (library) {
      TString check="libAliHLT"; check+=pAgent->GetModuleId();
      if (check.CompareTo(library)==0) {
	agents.clear();
	agents.push_back(pAgent);
	break;
      }
    }

    // check if the current agent is in the black list
    if (blackList) {
      const char* found=strstr(blackList, pAgent->GetModuleId());
      if (found) {
	found+=strlen(pAgent->GetModuleId());
	// skip this agent as it is in the blacklist
	if (*found==0 or *found==' ') continue;
      }
    }
    agents.push_back(pAgent);
  }

  for (vector<AliHLTModuleAgent*>::iterator element=agents.begin();
       element!=agents.end(); element++) {
    (*element)->ActivateComponentHandler(this);
  }

  return agents.size();
}

int AliHLTComponentHandler::DeleteOwnedComponents()
{
  // see header file for class documentation
  int iResult=0;
  AliHLTComponentPList::iterator element=fOwnedComponents.begin();
  while (element!=fOwnedComponents.end()) {
    //DeregisterComponent((*element)->GetComponentID());
    // exeption does not help in Root context, the Root exeption
    // handler always catches the exeption before. Have to find out
    // how exeptions can be used in Root
    /*try*/ {
      delete *element;
    }
    /*
    catch (...) {
      HLTError("delete managed sample %p", *element);
    }
    */
    fOwnedComponents.erase(element);
    element=fOwnedComponents.begin();
  }
  return iResult;
}

int AliHLTComponentHandler::SetRunDescription(const AliHLTRunDesc* desc, const char* runType)
{
  // see header file for class documentation
  if (!desc) return -EINVAL;
  if (desc->fStructSize!=sizeof(AliHLTRunDesc)) {
    HLTError("invalid size of RunDesc struct (%ul)", desc->fStructSize);
    return -EINVAL;
  }

  memcpy(&fRunDesc, desc, sizeof(AliHLTRunDesc));
  if (runType) {
    if (fRunType) delete [] fRunType;
    fRunType=new char[sizeof(runType)+1];
    if (fRunType) {
      strcpy(fRunType, runType);
    }
  }
  return 0;
}
