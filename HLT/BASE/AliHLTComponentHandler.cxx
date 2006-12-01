// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTComponentHandler.cxx
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Implementation of HLT component handler. */

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
#include "AliHLTStdIncludes.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTComponent.h"
#include "AliHLTDataTypes.h"
#include "AliHLTSystem.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTComponentHandler)

AliHLTComponentHandler::AliHLTComponentHandler()
  :
  fComponentList(),
  fScheduleList(),
  fLibraryList(),
  fEnvironment()
{
  memset(&fEnvironment, 0, sizeof(AliHLTComponentEnvironment));
}

AliHLTComponentHandler::~AliHLTComponentHandler()
{
  UnloadLibraries();
}

int AliHLTComponentHandler::AnnounceVersion()
{
  int iResult=0;
#ifdef PACKAGE_STRING
  void HLTbaseCompileInfo( char*& date, char*& time);
  char* date="";
  char* time="";
  HLTbaseCompileInfo(date, time);
  if (!date) date="unknown";
  if (!time) time="unknown";
  HLTInfo("%s build on %s (%s)", PACKAGE_STRING, date, time);
#else
  HLTInfo("ALICE High Level Trigger (embedded AliRoot build)");
#endif
  return iResult;
}

Int_t AliHLTComponentHandler::RegisterComponent(AliHLTComponent* pSample)
{
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
  int iResult=0;
  if (componentID) {
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

Int_t AliHLTComponentHandler::ScheduleRegister(AliHLTComponent* pSample)
{
  Int_t iResult=0;
  if (pSample) {
    fScheduleList.push_back(pSample);
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponentHandler::CreateComponent(const char* componentID, void* pEnv, int argc, const char** argv, AliHLTComponent*& component )
{
  int iResult=0;
  if (componentID) {
    AliHLTComponent* pSample=FindComponent(componentID);
    if (pSample!=NULL) {
      component=pSample->Spawn();
      if (component) {
	HLTDebug("component \"%s\" created (%p)", componentID, component);
	if ((iResult=component->Init(&fEnvironment, pEnv, argc, argv))!=0) {
	  HLTError("Initialization of component \"%s\" failed with error %d", componentID, iResult);
	  delete component;
	  component=NULL;
	}
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
  Int_t iResult=0;
  if (componentID) {
    vector<AliHLTComponent*>::iterator element=fComponentList.begin();
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
  AliHLTComponent* pSample=NULL;
  Int_t index=FindComponentIndex(componentID);
  if (index>=0) {
    pSample=(AliHLTComponent*)fComponentList.at(index);
  }
  return pSample;
}

Int_t AliHLTComponentHandler::InsertComponent(AliHLTComponent* pSample)
{
  Int_t iResult=0;
  if (pSample!=NULL) {
    fComponentList.push_back(pSample);
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTComponentHandler::List() {
  vector<AliHLTComponent*>::iterator element=fComponentList.begin();
  int index=0;
  while (element!=fComponentList.end()) {
    HLTInfo("%d. %s", index++, (*element++)->GetComponentID());
  }
}

void AliHLTComponentHandler::SetEnvironment(AliHLTComponentEnvironment* pEnv) {
  if (pEnv) {
    memcpy(&fEnvironment, pEnv, sizeof(AliHLTComponentEnvironment));
    AliHLTLogging::Init(fEnvironment.fLoggingFunc);
  }
}

int AliHLTComponentHandler::LoadLibrary( const char* libraryPath )
{
  int iResult=0;
  if (libraryPath) {
    AliHLTComponent::SetGlobalComponentHandler(this);
    AliHLTLibHandle hLib=NULL;
#ifdef HAVE_DLFCN_H
    // use interface to the dynamic linking loader
    hLib=dlopen(libraryPath, RTLD_NOW);
#else
    // use ROOT dynamic loader
    if (gSystem->Load(libraryPath)==0) {
      // create TString object to store library path and use pointer as handle 
      hLib=reinterpret_cast<AliHLTLibHandle>(new TString(libraryPath));
    }
#endif //HAVE_DLFCN_H
    if (hLib) {
      HLTInfo("library %s loaded", libraryPath);
      fLibraryList.push_back(hLib);
      vector<AliHLTComponent*>::iterator element=fScheduleList.begin();
      int iSize=fScheduleList.size();
      int iLocalResult=0;
      while (iSize-- > 0) {
	element=fScheduleList.begin();
	iLocalResult=RegisterComponent(*element);
	if (iResult==0) iResult=iLocalResult;
	fScheduleList.erase(element);
      }
    } else {
      HLTError("can not load library %s", libraryPath);
#ifdef HAVE_DLFCN_H
      HLTError("dlopen error: %s", dlerror());
#endif //HAVE_DLFCN_H
      iResult=-ELIBACC;
    }
    AliHLTComponent::UnsetGlobalComponentHandler();
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponentHandler::UnloadLibrary( const char* libraryPath )
{
  int iResult=0;
  if (libraryPath) {
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponentHandler::UnloadLibraries()
{
  int iResult=0;
  vector<AliHLTLibHandle>::iterator element=fLibraryList.begin();
  while (element!=fLibraryList.end()) {
#ifdef HAVE_DLFCN_H
    dlclose(*element);
#else
    TString* libraryPath=reinterpret_cast<TString*>(*element);
    gSystem->Unload(libraryPath->Data());
    delete libraryPath;
#endif //HAVE_DLFCN_H
    element++;
  }
  return iResult;
}
