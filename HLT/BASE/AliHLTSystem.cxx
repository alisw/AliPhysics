// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTSystem.cxx
    @author Matthias Richter
    @date   
    @brief  Implementation of HLT module management.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include <cassert>
#include "AliHLTStdIncludes.h"
#include "AliHLTSystem.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTComponent.h"
#include "AliHLTConfiguration.h"
#include "AliHLTConfigurationHandler.h"
#include "AliHLTTask.h"
#include "AliHLTModuleAgent.h"
#include "AliHLTOfflineInterface.h"
#include <TObjArray.h>
#include <TObjString.h>
#include <TStopwatch.h>
#include <TROOT.h>
#include <TInterpreter.h>

/** HLT default component libraries */
const char* kHLTDefaultLibs[]= {
  "libAliHLTUtil.so", 
  "libAliHLTTPC.so", 
  //  "libAliHLTSample.so",
  "libAliHLTPHOS.so",
  //"libAliHLTMUON.so",
  "libAliHLTTRD.so",
  NULL
};

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSystem)

AliHLTSystem::AliHLTSystem()
  :
  fpComponentHandler(new AliHLTComponentHandler()),
  fpConfigurationHandler(new AliHLTConfigurationHandler()),
  fTaskList(),
  fState(0),
  fChains(),
  fStopwatches(new TObjArray),
  fEventCount(-1),
  fGoodEvents(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  if (fgNofInstances++>0)
    HLTWarning("multiple instances of AliHLTSystem, you should not use more than one at a time");

  SetGlobalLoggingLevel(kHLTLogDefault);
  if (fpComponentHandler) {
    AliHLTComponentEnvironment env;
    memset(&env, 0, sizeof(AliHLTComponentEnvironment));
    env.fAllocMemoryFunc=AliHLTSystem::AllocMemory;
    env.fLoggingFunc=NULL;
    fpComponentHandler->SetEnvironment(&env);
    InitAliLogFunc(fpComponentHandler);
    fpComponentHandler->AnnounceVersion();
  } else {
    HLTFatal("can not create Component Handler");
  }
  if (fpConfigurationHandler) {
    AliHLTConfiguration::GlobalInit(fpConfigurationHandler);
  } else {
    HLTFatal("can not create Configuration Handler");
  }
}

AliHLTSystem::~AliHLTSystem()
{
  // see header file for class documentation
  fgNofInstances--;
  CleanTaskList();
  AliHLTConfiguration::GlobalDeinit(fpConfigurationHandler);
  if (fpConfigurationHandler) {
    delete fpConfigurationHandler;
  }
  fpConfigurationHandler=NULL;
  
  if (fpComponentHandler) {
    delete fpComponentHandler;
  }
  fpComponentHandler=NULL;
}

int AliHLTSystem::fgNofInstances=0;

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

int AliHLTSystem::BuildTaskList(const char* id)
{
  // see header file for class documentation
  int iResult=0;
  if (id) {
    if (fpConfigurationHandler) {
      AliHLTConfiguration* pConf=fpConfigurationHandler->FindConfiguration(id);
      if (pConf) {
	iResult=BuildTaskList(pConf);
      } else {
	HLTError("unknown configuration \"%s\"", id);
	iResult=-EEXIST;
      }
    } else {
      iResult=-EFAULT;
    }
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
      }
      // task for this configuration exists, terminate
      pTask=NULL;
    } else if (pConf->SourcesResolved(1)!=1) {
	HLTError("configuration \"%s\" has unresolved sources, aborting ...", pConf->GetName());
	iResult=-ENOLINK;
    } else {
      pTask=new AliHLTTask(pConf);
      if (pTask==NULL) {
	iResult=-ENOMEM;
      }
    }
    static int iterationLevel=0;
    if (pTask && iResult>=0) {
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
	  HLTDebug("iteration %d: checking dependency %s (%p)", iterationLevel, pDep->GetName(), pDep);
	  if (FindTask(pDep->GetName())==NULL) {
	    HLTDebug("iteration %d: building task list for configuration %s (%p)", iterationLevel, pDep->GetName(), pDep);
	    iterationLevel++;
	    iResult=BuildTaskList(pDep);
	    iterationLevel--;
	  }
	  pDep=pConf->GetNextSource();
	}
	// remove the temporarily added task
	fTaskList.Remove(pTask);

	// insert the task and set the cross-links
	if (iResult>=0) {
	  HLTDebug("iteration %d: inserting task %s (%p)", iterationLevel, pTask->GetName(), pTask);
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
  while ((lnk=fTaskList.LastLink())!=NULL) {
    delete (lnk->GetObject());
    fTaskList.Remove(lnk);
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
      HLTDebug("task \"%s\" (%p) inserted (size %d)", pTask->GetName(), pTask, sizeof(AliHLTTask));
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
    pTask=dynamic_cast<AliHLTTask*>(fTaskList.FindObject(id)); 
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

int AliHLTSystem::Run(Int_t iNofEvents, int bStop)
{
  // see header file for class documentation
  int iResult=0;
  int iCount=0;
  SetStatusFlags(kRunning);
  if (fEventCount>=0 || (iResult=InitTasks())>=0) {
    if (fEventCount>=0 || (iResult=StartTasks())>=0) {
      if (fEventCount==0) {
	InitBenchmarking(fStopwatches);
      } else {
	//ResumeBenchmarking(fStopwatches);    
      }
      for (int i=fEventCount; i<fEventCount+iNofEvents && iResult>=0; i++) {
	if ((iResult=ProcessTasks(i))>=0) {
	  fGoodEvents++;
	  iCount++;
	} else {
	  // TODO: define different running modes to either ignore errors in
	  // event processing or not
	  // currently ignored 
	  iResult=0;
	}
      }
      fEventCount+=iNofEvents;
      if (bStop) StopTasks();
      //else PauseBenchmarking(fStopwatches);
    }
    if (bStop) DeinitTasks();
  }
  if (iResult>=0) {
    iResult=iCount;
  } else  if (iResult==-ENOKEY) {
    iResult=0; // do not propagate the error
  }
  ClearStatusFlags(kRunning);
  return iResult;
}

int AliHLTSystem::InitTasks()
{
  // see header file for class documentation
  int iResult=0;
  TObjLink *lnk=fTaskList.FirstLink();

  if (lnk==NULL) {
    HLTWarning("Task list is empty, aborting ...");
    return -ENOKEY;
  }
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
    HLTError("can not initialize task list, error %d", iResult);
  }

  return iResult;
}

int AliHLTSystem::InitBenchmarking(TObjArray* pStopwatches)
{
  // see header file for class documentation
  int iResult=0;
  if (pStopwatches==NULL) return 0;

  for (int i=0; i<(int)AliHLTComponent::kSWTypeCount; i++) {
    TStopwatch* pStopwatch= new TStopwatch;
    if (pStopwatch) {
      pStopwatch->Reset();
      pStopwatches->AddAt(pStopwatch, i);
    } else {
      iResult=-ENOMEM;
      break;
    }
  }

  TObjLink *lnk=fTaskList.FirstLink();
  while (lnk && iResult>=0) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      AliHLTComponent* pComp=NULL;
      if (iResult>=0 && (pComp=pTask->GetComponent())!=NULL) {
	switch (pComp->GetComponentType()) {
	case AliHLTComponent::kProcessor:
	  pComp->SetStopwatches(pStopwatches);
	  break;
	case AliHLTComponent::kSource:
	  {
	    // this switch determines whether the time consumption of the
	    // AliHLTComponent base methods should be counted to the input
	    // stopwatch or base stopwatch.
	    //int inputBase=(int)AliHLTComponent::kSWBase;
	    int inputBase=(int)AliHLTComponent::kSWInput;
	    pComp->SetStopwatch(pStopwatches->At(inputBase), AliHLTComponent::kSWBase);
	    pComp->SetStopwatch(pStopwatches->At((int)AliHLTComponent::kSWInput), AliHLTComponent::kSWDA);
	  }
	  break;
	case AliHLTComponent::kSink:
	  {
	    // this switch determines whether the time consumption of the
	    // AliHLTComponent base methods should be counted to the output
	    // stopwatch or base stopwatch.
	    //int outputBase=(int)AliHLTComponent::kSWBase;
	    int outputBase=(int)AliHLTComponent::kSWOutput;
	    pComp->SetStopwatch(pStopwatches->At(outputBase), AliHLTComponent::kSWBase);
	    pComp->SetStopwatch(pStopwatches->At((int)AliHLTComponent::kSWOutput), AliHLTComponent::kSWDA);
	  }
	  break;
	default:
	  HLTWarning("unknown component type %d", (int)pComp->GetComponentType());
	}
      }
    } else {
    }
    lnk = lnk->Next();
  }
  return iResult;
}

int AliHLTSystem::PrintBenchmarking(TObjArray* pStopwatches, int bClean)
{
  // see header file for class documentation
  int iInitialized=1;
  if (pStopwatches==NULL) return 0;

  for (int i=0; i<(int)AliHLTComponent::kSWTypeCount; i++) {
    if (!dynamic_cast<TStopwatch*>(pStopwatches->At(i))) {
      iInitialized=0;
      break;
    }
  }

  if (iInitialized!=0) {
    HLTInfo("HLT statistics:\n"
	    "    base:              R:%.3fs C:%.3fs\n"
	    "    input:             R:%.3fs C:%.3fs\n"
	    "    output:            R:%.3fs C:%.3fs\n"
	    "    event processing : R:%.3fs C:%.3fs"
	    , dynamic_cast<TStopwatch*>(pStopwatches->At(AliHLTComponent::kSWBase))->RealTime()
	    , dynamic_cast<TStopwatch*>(pStopwatches->At(AliHLTComponent::kSWBase))->CpuTime()
	    , dynamic_cast<TStopwatch*>(pStopwatches->At(AliHLTComponent::kSWInput))->RealTime()
	    , dynamic_cast<TStopwatch*>(pStopwatches->At(AliHLTComponent::kSWInput))->CpuTime()
	    , dynamic_cast<TStopwatch*>(pStopwatches->At(AliHLTComponent::kSWOutput))->RealTime()
	    , dynamic_cast<TStopwatch*>(pStopwatches->At(AliHLTComponent::kSWOutput))->CpuTime()
	    , dynamic_cast<TStopwatch*>(pStopwatches->At(AliHLTComponent::kSWDA))->RealTime()
	    , dynamic_cast<TStopwatch*>(pStopwatches->At(AliHLTComponent::kSWDA))->CpuTime()
	    );
  }

  if (bClean) {
    for (int i=0; i<(int)AliHLTComponent::kSWTypeCount; i++) {
      TObject* pObj=pStopwatches->RemoveAt(i);
      if (pObj) delete pObj;
    }
  }
  return 0;
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
    HLTError("can not start task list, error %d", iResult);
  } else {
    fEventCount=0;
    fGoodEvents=0;
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

  if (iResult>=0) {
    HLTInfo("Event %d successfully finished (%d)", eventNo, iResult);
    iResult=0;
  } else {
    HLTError("Processing of event %d failed (%d)", eventNo, iResult);
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
  PrintBenchmarking(fStopwatches, 1 /*clean*/);
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
  fEventCount=-1;
  fGoodEvents=-1;

  return iResult;
}

void* AliHLTSystem::AllocMemory( void* /*param*/, unsigned long size )
{
  // see header file for class documentation
  void* p=NULL;
  try {
    p=(void*)new char[size];
  }
  catch (...) {
    AliHLTLogging log;
    log.LoggingVarargs(kHLTLogError, "AliHLTSystem" , "AllocMemory" , __FILE__ , __LINE__ , "exeption during memory allocation" );
  }
  return p;
}

int AliHLTSystem::Reconstruct(int nofEvents, AliRunLoader* runLoader, 
			      AliRawReader* rawReader)
{
  // see header file for class documentation
  int iResult=0;
  if (runLoader || rawReader || nofEvents==0) {
    if (nofEvents>0) {HLTInfo("Run Loader %p, Raw Reader %p , %d event(s)", runLoader, rawReader, nofEvents);}
    if (CheckStatus(kReady)) {
      if (nofEvents==0) {
	// special case to close the reconstruction
	if (!CheckStatus(kError)) {
	StopTasks();
	DeinitTasks();
	}
      } else {
      if ((iResult=AliHLTOfflineInterface::SetParamsToComponents(runLoader, rawReader))>=0) {
	// the system always remains started after event processing, a specific
	// call with nofEvents==0 is needed to execute the stop sequence
	if ((iResult=Run(nofEvents, 0))<0) SetStatusFlags(kError);
      }
      }
    } else {
      HLTError("wrong state %#x, required flags %#x", GetStatusFlags(), kReady);
    }
  } else {
    HLTError("missing RunLoader (%p)/RawReader (%p) instance", runLoader, rawReader);
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTSystem::FillESD(int eventNo, AliRunLoader* runLoader, AliESDEvent* esd)
{
  // see header file for class documentation
  int iResult=0;
  if (runLoader || esd) {
    HLTInfo("Event %d: Run Loader %p, ESD %p", eventNo, runLoader, esd);
    iResult=AliHLTOfflineInterface::FillComponentESDs(eventNo, runLoader, esd);
  } else {
    HLTError("missing run loader/ESD instance(s)");
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTSystem::LoadComponentLibraries(const char* libraries)
{
  // see header file for class documentation
  int iResult=0;
  if (libraries) {
    if (fpComponentHandler) {
      TString libs(libraries);
      TObjArray* pTokens=libs.Tokenize(" ");
      if (pTokens) {
	int iEntries=pTokens->GetEntries();
	for (int i=0; i<iEntries && iResult>=0; i++) {
	  iResult=fpComponentHandler->LoadLibrary((((TObjString*)pTokens->At(i))->GetString()).Data());
	}
	delete pTokens;
      }
      if (iResult>=0) {
	SetStatusFlags(kLibrariesLoaded);
      } else {
	// lets see if we need this, probably not
	//fpComponentHandler->UnloadLibraries();
	ClearStatusFlags(kLibrariesLoaded);
      }
    } else {
      iResult=-EFAULT;
      HLTFatal("no component handler available");
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTSystem::Configure(AliRunLoader* runloader)
{
  // see header file for class documentation
  return Configure(NULL, runloader);
}

int AliHLTSystem::Configure(AliRawReader* rawReader, AliRunLoader* runloader)
{
  // see header file for class documentation
  int iResult=0;
  if (CheckStatus(kRunning)) {
    HLTError("HLT system in running state, can not configure");
    return -EBUSY;
  }
  ClearStatusFlags(kConfigurationLoaded|kTaskListCreated);
  if (CheckFilter(kHLTLogDebug))
    AliHLTModuleAgent::PrintStatus();
  if (CheckStatus(kConfigurationLoaded)==0) {
    iResult=LoadConfigurations(rawReader, runloader);
  } else {
    if (fChains.Length()==0) {
      HLTError("custom configuration(s) specified, but no configuration to run in local reconstruction, use \'localrec=<conf>\' option");
      iResult=-ENOENT;
    }
  }
  if (iResult>=0) {
    SetStatusFlags(kConfigurationLoaded);
    if (CheckFilter(kHLTLogDebug))
      fpConfigurationHandler->PrintConfigurations();
    iResult=BuildTaskListsFromReconstructionChains(rawReader, runloader);
    if (iResult>=0) {
      SetStatusFlags(kTaskListCreated);
    }
  }
  if (iResult<0) SetStatusFlags(kError);
  
  return iResult;
}

int AliHLTSystem::ScanOptions(const char* options)
{
  // see header file for class documentation
  int iResult=0;
  if (options) {
    TString libs("");
    TString alloptions(options);
    TObjArray* pTokens=alloptions.Tokenize(" ");
    if (pTokens) {
      int iEntries=pTokens->GetEntries();
      for (int i=0; i<iEntries; i++) {
	TString token=(((TObjString*)pTokens->At(i))->GetString());
	if (token.Contains("loglevel=")) {
	  TString param=token.ReplaceAll("loglevel=", "");
	  if (param.IsDigit()) {
	    SetGlobalLoggingLevel((AliHLTComponentLogSeverity)param.Atoi());
	  } else if (param.BeginsWith("0x") &&
		     param.Replace(0,2,"",0).IsHex()) {
	    int severity=0;
	    sscanf(param.Data(),"%x", &severity);
	    SetGlobalLoggingLevel((AliHLTComponentLogSeverity)severity);
	  } else {
	    HLTWarning("wrong parameter for option \'loglevel=\', (hex) number expected");
	  }
	} else if (token.Contains("alilog=off")) {
	  SwitchAliLog(0);
	} else if (token.Contains("config=")) {
	  TString param=token.ReplaceAll("config=", "");
	  Int_t error=0;
	  gROOT->Macro(param.Data(), &error);
	  if (error==0) {
	    SetStatusFlags(kConfigurationLoaded);
	  } else {
	    HLTError("can not execute macro \'%s\'", param.Data());
	    iResult=-EBADF;
	  }
	} else if (token.Contains("chains=")) {
	  TString param=token.ReplaceAll("chains=", "");
	  fChains=param.ReplaceAll(",", " ");
	} else if (token.BeginsWith("lib") && token.EndsWith(".so")) {
	  libs+=token;
	  libs+=" ";
	} else {
	  HLTWarning("unknown option \'%s\'", token.Data());
	}
      }
      delete pTokens;
    }

    if (iResult>=0) {
      if (libs.IsNull()) {
	const char** deflib=kHLTDefaultLibs;
	while (*deflib) {
	  libs+=*deflib++;
	  libs+=" ";
	}
      }
      if ((!CheckStatus(AliHLTSystem::kLibrariesLoaded)) &&
	  (LoadComponentLibraries(libs.Data())<0)) {
	HLTError("error while loading HLT libraries");
	iResult=-EFAULT;
      }
    }
  }
  return iResult;
}

int AliHLTSystem::Reset(int bForce)
{
  // see header file for class documentation
  int iResult=0;
  if (!bForce && CheckStatus(kRunning)) {
    HLTError("HLT system in running state, can not configure");
    return -EBUSY;
  }
  CleanTaskList();
  ClearStatusFlags(~kUninitialized);
  return iResult;
}

int AliHLTSystem::LoadConfigurations(AliRawReader* rawReader, AliRunLoader* runloader)
{
  // see header file for class documentation
  if (CheckStatus(kRunning)) {
    HLTError("HLT system in running state, can not configure");
    return -EBUSY;
  }
  int iResult=0;
  AliHLTModuleAgent* pAgent=AliHLTModuleAgent::GetFirstAgent();
  while (pAgent && iResult>=0) {
    const char* deplibs=pAgent->GetRequiredComponentLibraries();
    if (deplibs) {
      HLTDebug("load libraries \'%s\' for agent %s (%p)", deplibs, pAgent->GetName(), pAgent);
      iResult=LoadComponentLibraries(deplibs);
    }
    if (iResult>=0) {
      HLTDebug("load configurations for agent %s (%p)", pAgent->GetName(), pAgent);
      pAgent->CreateConfigurations(fpConfigurationHandler, rawReader, runloader);
      pAgent=AliHLTModuleAgent::GetNextAgent();
    }
  }
  return iResult;
}

int AliHLTSystem::BuildTaskListsFromReconstructionChains(AliRawReader* rawReader, AliRunLoader* runloader)
{
  // see header file for class documentation
  if (CheckStatus(kRunning)) {
    HLTError("HLT system in running state, can not configure");
    return -EBUSY;
  }
  if (!CheckStatus(kConfigurationLoaded)) {
    HLTWarning("configurations not yet loaded");
    return 0;
  }

  int iResult=0;
  int bHaveOutput=0;

  // query chains
  TString chains;
  if (fChains.Length()>0) {
    chains=fChains;
    HLTInfo("custom reconstruction chain: %s", chains.Data());
  } else {
    AliHLTModuleAgent* pAgent=AliHLTModuleAgent::GetFirstAgent();
    while ((pAgent || fChains.Length()>0) && iResult>=0) {
      const char* agentchains=pAgent->GetReconstructionChains(rawReader, runloader);
      if (agentchains) {
	if (!chains.IsNull()) chains+="";
	chains+=agentchains;
	HLTInfo("reconstruction chains for agent %s (%p): %s", pAgent->GetName(), pAgent, agentchains);
      }
      pAgent=AliHLTModuleAgent::GetNextAgent();
    }
  }

  // build task list for chains
  TObjArray* pTokens=chains.Tokenize(" ");
  if (pTokens) {
    int iEntries=pTokens->GetEntries();
    for (int i=0; i<iEntries && iResult>=0; i++) {
      const char* pCID=((TObjString*)pTokens->At(i))->GetString().Data();
      AliHLTConfiguration* pConf=fpConfigurationHandler->FindConfiguration(pCID);
      if (pConf) {
	iResult=BuildTaskList(pConf);
	if (runloader) {
	  assert(fpComponentHandler!=NULL);
	  TString cid=pConf->GetComponentID();
	  if (cid.CompareTo("HLTOUT")==0) {
	    // remove from the input of a global HLTOUT configuration
	    chains.ReplaceAll(pCID, "");
	  } else if (bHaveOutput==0) {
	    // check whether this configuration produces data output
	    if ((bHaveOutput=fpComponentHandler->HasOutputData(cid.Data()))<0) {
	      bHaveOutput=0;
	      chains.ReplaceAll(pCID, "");
	    }
	  }
	}
      } else {
	HLTWarning("can not find configuration %s", pCID);
      }
    }
    delete pTokens;
  }

  // build HLTOUT for simulation
  if (iResult>=0 && runloader) {
    if (bHaveOutput) {
      // there are components in the chain which produce data which need to be
      // piped to an HLTOUT
      if (fpComponentHandler->FindComponentIndex("HLTOUT")>=0 ||
	  LoadComponentLibraries("libHLTsim.so")>=0) {
	AliHLTConfiguration globalout("_globalout_", "HLTOUT", chains.Data(), NULL);
	iResult=BuildTaskList("_globalout_");
      } else {
	HLTError("can not load libHLTsim.so and HLTOUT component");
	iResult=-EFAULT;
      }
    }
  }

  if (iResult>=0) SetStatusFlags(kTaskListCreated);

  return iResult;
}

int AliHLTSystem::CheckStatus(int flag)
{
  // see header file for class documentation
  if (flag==kUninitialized && flag==fState) return 1;
  if ((fState&flag)==flag) return 1;
  return 0;
}

int AliHLTSystem::GetStatusFlags()
{
  // see header file for class documentation
  return fState;
}

int AliHLTSystem::SetStatusFlags(int flags)
{
  // see header file for class documentation
  fState|=flags;
  return fState;
}

int AliHLTSystem::ClearStatusFlags(int flags)
{
  // see header file for class documentation
  fState&=~flags;
  return fState;
}

void* AliHLTSystem::FindDynamicSymbol(const char* library, const char* symbol)
{
  // see header file for class documentation
  if (fpComponentHandler==NULL) return NULL;
  return fpComponentHandler->FindSymbol(library, symbol);
}
