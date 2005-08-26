// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          Artur Szostak <artursz@iafrica.com>                           *
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// global HLT module management                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__== 3
using namespace std;
#endif

#include <errno.h>
#include <string.h>
#include "AliL3StandardIncludes.h"
#include "AliHLTSystem.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTComponent.h"
#include "AliHLTConfiguration.h"

ClassImp(AliHLTSystem)

AliHLTSystem::AliHLTSystem()
{
  fpComponentHandler=new AliHLTComponentHandler();
  if (fpComponentHandler) {
    AliHLTComponentEnvironment env;
    memset(&env, 0, sizeof(AliHLTComponentEnvironment));
    env.fLoggingFunc=AliHLTLogging::Message;
    fpComponentHandler->SetEnvironment(&env);

    // init logging function in AliHLTLogging
    Init(AliHLTLogging::Message);
  } else {
    HLTFatal("can not create Component Handler");
  }
  fpConfigurationHandler=new AliHLTConfigurationHandler();
  if (fpConfigurationHandler) {
    AliHLTConfiguration::GlobalInit(fpConfigurationHandler);
  } else {
    HLTFatal("can not create Configuration Handler");
  }
}


AliHLTSystem::~AliHLTSystem()
{
    AliHLTConfiguration::GlobalDeinit();
    if (fpConfigurationHandler) {
      delete fpConfigurationHandler;
    }
    fpConfigurationHandler=NULL;

    if (fpComponentHandler) {
      delete fpComponentHandler;
    }
    fpComponentHandler=NULL;
}

int AliHLTSystem::AddConfiguration(AliHLTConfiguration* pConf)
{
  int iResult=0;
  return iResult;
}

// int AliHLTSystem::InsertConfiguration(AliHLTConfiguration* pConf, AliHLTConfiguration* pPrec)
// {
//   int iResult=0;
//   return iResult;
// }

int AliHLTSystem::DeleteConfiguration(AliHLTConfiguration* pConf)
{
  int iResult=0;
  return iResult;
}

int AliHLTSystem::BuildTaskList(AliHLTConfiguration* pConf)
{
  int iResult=0;
  if (pConf) {
    AliHLTTask* pTask=NULL;
    if ((pTask=FindTask(pConf->GetName()))!=NULL) {
      if (pTask->GetConf()!=pConf) {
	HLTError("configuration missmatch, there is already a task with configuration name \"%s\", but it is different. Most likely configuration %p is not registered properly", pConf->GetName(), pConf);
	iResult=-EEXIST;
	pTask=NULL;
      }
    } else if (pConf->SourcesResolved(1)!=1) {
	HLTError("configuration \"%s\" has unresolved sources, aborting ...", pConf->GetName());
	iResult=-ENOLINK;
    } else {
      pTask=new AliHLTTask(pConf, NULL);
      if (pTask==NULL) {
	iResult=-ENOMEM;
      }
    }
    if (pTask) {
      // check for ring dependencies
      if ((iResult=pConf->FollowDependency(pConf->GetName()))>0) {
	HLTError("detected ring dependency for configuration \"%s\"", pTask->GetName());
	pTask->PrintDependencyTree(pTask->GetName(), 1/*use the configuration list*/);
	HLTError("aborted ...");
	iResult=-ELOOP;
      }
      if (iResult>=0) {
	// check whether all dependencies are already in the task list
	// create the missing ones
	fTaskList.Add(pTask);
	AliHLTConfiguration* pDep=pConf->GetFirstSource();
	while (pDep!=NULL && iResult>=0) {
	  if (FindTask(pDep->GetName())==NULL) {
	    iResult=BuildTaskList(pDep);
	  }
	  pDep=pConf->GetNextSource();
	}
	fTaskList.Remove(pTask);

	// insert the task and set the cross-links
	if (iResult>=0) {
	  iResult=InsertTask(pTask);
	}
      } else {
	delete pTask;
	pTask=NULL;
      }
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTSystem::CleanTaskList()
{
  int iResult=0;
  TObjLink* lnk=NULL;
  while ((lnk=fTaskList.FirstLink())!=NULL) {
    fTaskList.Remove(lnk);
    delete (lnk->GetObject());
  }
  return iResult;
}

int AliHLTSystem::InsertTask(AliHLTTask* pTask)
{
  int iResult=0;
  TObjLink *lnk = NULL;
  if ((iResult=pTask->CheckDependencies())>0)
    lnk=fTaskList.FirstLink();
  while (lnk && iResult>0) {
    AliHLTTask* pCurr = (AliHLTTask*)lnk->GetObject();
    //HLTDebug("checking  \"%s\"", pCurr->GetName());
    iResult=pTask->Depends(pCurr);
    if (iResult>0) {
      iResult=pTask->SetDependency(pCurr);
      pCurr->SetTarget(pTask);
      HLTDebug("set dependency  \"%s\" for configuration \"%s\"", pCurr->GetName(), pTask->GetName());
    }
    if (pCurr->Depends(pTask)) {
      // ring dependency
      HLTError("ring dependency: can not resolve dependencies for configuration \"%s\"", pTask->GetName());
      iResult=-ELOOP;
    } else if ((iResult=pTask->CheckDependencies())>0) {
      lnk = lnk->Next();
    }
  }
  if (iResult==0) {
      if (lnk) {
	fTaskList.AddAfter(lnk, pTask);
      } else {
	fTaskList.AddFirst(pTask);
      }
      HLTDebug("task \"%s\" inserted", pTask->GetName());
  } else if (iResult>0) {
    HLTError("can not resolve dependencies for configuration \"%s\" (%d unresolved)", pTask->GetName(), iResult);
    iResult=-ENOLINK;
  }
  return iResult;
}

AliHLTTask* AliHLTSystem::FindTask(const char* id)
{
  AliHLTTask* pTask=NULL;
  if (id) {
    pTask=(AliHLTTask*)fTaskList.FindObject(id); 
  }
  return pTask;
}

void AliHLTSystem::PrintTaskList()
{
  HLTLogKeyword("task list");
  TObjLink *lnk = NULL;
  HLTMessage("Task List");
  lnk=fTaskList.FirstLink();
  while (lnk) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      HLTMessage("  %s - status:", obj->GetName());
      AliHLTTask* pTask=(AliHLTTask*)obj;
      pTask->PrintStatus();
    } else {
    }
    lnk = lnk->Next();
  }
}
