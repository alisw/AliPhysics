// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
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

/** @file   AliHLTSystem.cxx
    @author Matthias Richter
    @date   
    @brief  Implementation of HLT module management.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTStdIncludes.h"
#include "AliHLTSystem.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTComponent.h"
#include "AliHLTConfiguration.h"
#include "AliHLTConfigurationHandler.h"
#include "AliHLTTask.h"
#include "TString.h"

// #include <sstream>
// #include <iostream>
// #include "AliLog.h"

// ostringstream g_logstr;

// void LogNotification(AliLog::EType_t level, const char* message)
// {
//   cout << "notification handler" << endl;
//   cout << g_logstr.str() << endl;
//   g_logstr.clear();
// }

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSystem)

AliHLTSystem::AliHLTSystem()
  :
  fpComponentHandler(new AliHLTComponentHandler()),
  fpConfigurationHandler(new AliHLTConfigurationHandler()),
  fTaskList()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  if (fpComponentHandler) {
    AliHLTComponentEnvironment env;
    memset(&env, 0, sizeof(AliHLTComponentEnvironment));
    env.fAllocMemoryFunc=AliHLTSystem::AllocMemory;
    env.fLoggingFunc=AliHLTLogging::Message;
    fpComponentHandler->SetEnvironment(&env);

    // init logging function in AliHLTLogging
    Init(AliHLTLogging::Message);
  } else {
    HLTFatal("can not create Component Handler");
  }
  if (fpConfigurationHandler) {
    AliHLTConfiguration::GlobalInit(fpConfigurationHandler);
  } else {
    HLTFatal("can not create Configuration Handler");
  }
//   AliLog log;
//   log.SetLogNotification(LogNotification);
//   log.SetStreamOutput(&g_logstr);
//   AliInfo("this is a printf message");
//   AliInfoStream() << "this is a stream message" << endl;
}

AliHLTSystem::AliHLTSystem(const AliHLTSystem&)
  :
  AliHLTLogging(),
  fpComponentHandler(NULL),
  fpConfigurationHandler(NULL),
  fTaskList()
{
  // see header file for class documentation
  HLTFatal("copy constructor untested");
}

AliHLTSystem& AliHLTSystem::operator=(const AliHLTSystem&)
{ 
  // see header file for class documentation
  HLTFatal("assignment operator untested");
  return *this;
}

AliHLTSystem::~AliHLTSystem()
{
  // see header file for class documentation
  CleanTaskList();
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
  // see header file for class documentation
  int iResult=0;
  if (pConf) {
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTSystem::InsertConfiguration(AliHLTConfiguration* pConf, AliHLTConfiguration* pPrec)
{
  // see header file for class documentation
  int iResult=0;
  if (pConf) {
    if (pPrec) {
      // find the position
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTSystem::DeleteConfiguration(AliHLTConfiguration* pConf)
{
  // see header file for class documentation
  int iResult=0;
  if (pConf) {
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTSystem::BuildTaskList(AliHLTConfiguration* pConf)
{
  // see header file for class documentation
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
      pTask=new AliHLTTask(pConf);
      if (pTask==NULL) {
	iResult=-ENOMEM;
      }
    }
    if (pTask) {
      // check for circular dependencies
      if ((iResult=pConf->FollowDependency(pConf->GetName()))>0) {
	HLTError("detected circular dependency for configuration \"%s\"", pTask->GetName());
	pTask->PrintDependencyTree(pTask->GetName(), 1/*use the configuration list*/);
	HLTError("aborted ...");
	iResult=-ELOOP;
      }
      if (iResult>=0) {
	// check whether all dependencies are already in the task list
	// create the missing ones
	// this step is an iterative process which calls this function again for the missing
	// configurations, in order to avoid the currently processed task to be created
	// again it is added to the list temporarily and removed afterwards
	// This is of high importance to preserve the order of the tasks. Furthermore, the
	// InsertTask method has to be used in order to set all the cross links right 
	fTaskList.Add(pTask);
	AliHLTConfiguration* pDep=pConf->GetFirstSource();
	while (pDep!=NULL && iResult>=0) {
	  if (FindTask(pDep->GetName())==NULL) {
	    iResult=BuildTaskList(pDep);
	  }
	  pDep=pConf->GetNextSource();
	}
	// remove the temporarily added task
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
  // see header file for class documentation
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
  // see header file for class documentation
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
      // circular dependency
      HLTError("circular dependency: can not resolve dependencies for configuration \"%s\"", pTask->GetName());
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
  // see header file for class documentation
  AliHLTTask* pTask=NULL;
  if (id) {
    pTask=(AliHLTTask*)fTaskList.FindObject(id); 
  }
  return pTask;
}

void AliHLTSystem::PrintTaskList()
{
  // see header file for class documentation
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

int AliHLTSystem::Run(Int_t iNofEvents) 
{
  // see header file for class documentation
  int iResult=0;
  if ((iResult=InitTasks())>=0) {
    if ((iResult=StartTasks())>=0) {
      for (int i=0; i<iNofEvents && iResult>=0; i++) {
	iResult=ProcessTasks(i);
	if (iResult>=0) {
	  HLTInfo("Event %d successfully finished (%d)", i, iResult);
	  iResult=0;
	} else {
	  HLTError("Processing of event %d failed (%d)", i, iResult);
	  // TODO: define different running modes to either ignore errors in
	  // event processing or not
	  // currently ignored 
	  //iResult=0;
	}
      }
      StopTasks();
    } else {
      HLTError("can not start task list");
    }
    DeinitTasks();
  } else {
    HLTError("can not initialize task list");
  }
  return iResult;
}

int AliHLTSystem::InitTasks()
{
  // see header file for class documentation
  int iResult=0;
  TObjLink *lnk=fTaskList.FirstLink();
  while (lnk && iResult>=0) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      iResult=pTask->Init(NULL, fpComponentHandler);
    } else {
    }
    lnk = lnk->Next();
  }
  if (iResult<0) {
  }
  return iResult;
}

int AliHLTSystem::StartTasks()
{
  // see header file for class documentation
  int iResult=0;
  TObjLink *lnk=fTaskList.FirstLink();
  while (lnk && iResult>=0) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      iResult=pTask->StartRun();
    } else {
    }
    lnk = lnk->Next();
  }
  if (iResult<0) {
  }
  return iResult;
}

int AliHLTSystem::ProcessTasks(Int_t eventNo)
{
  // see header file for class documentation
  int iResult=0;
  HLTDebug("processing event no %d", eventNo);
  TObjLink *lnk=fTaskList.FirstLink();
  while (lnk && iResult>=0) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      iResult=pTask->ProcessTask(eventNo);
      HLTDebug("task %s finnished (%d)", pTask->GetName(), iResult);
    } else {
    }
    lnk = lnk->Next();
  }
  return iResult;
}

int AliHLTSystem::StopTasks()
{
  // see header file for class documentation
  int iResult=0;
  TObjLink *lnk=fTaskList.FirstLink();
  while (lnk && iResult>=0) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      iResult=pTask->EndRun();
    } else {
    }
    lnk = lnk->Next();
  }
  return iResult;
}

int AliHLTSystem::DeinitTasks()
{
  // see header file for class documentation
  int iResult=0;
  TObjLink *lnk=fTaskList.FirstLink();
  while (lnk && iResult>=0) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      iResult=pTask->Deinit();
    } else {
    }
    lnk = lnk->Next();
  }
  return iResult;
}

void* AliHLTSystem::AllocMemory( void* param, unsigned long size )
{
  // see header file for class documentation
  return (void*)new char[size];
}
