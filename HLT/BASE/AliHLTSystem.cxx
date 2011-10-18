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
#include "AliHLTOnlineConfiguration.h"
#include "AliHLTTask.h"
#include "AliHLTModuleAgent.h"
#include "AliHLTOfflineInterface.h"
#include "AliHLTDataSource.h"
#include "AliHLTOUT.h"
#include "AliHLTOUTHandler.h"
#include "AliHLTOUTTask.h"
#include "AliHLTControlTask.h"
#include "AliHLTDataBuffer.h"
#include "AliHLTMisc.h"
#include <TObjArray.h>
#include <TObjString.h>
#include <TStopwatch.h>
#include <TList.h>
//#include <TSystem.h>
#include <TROOT.h>
//#include <TInterpreter.h>

/** HLT default component libraries */
const char* AliHLTSystem::fgkHLTDefaultLibs[]= {
  "libAliHLTUtil.so", 
  "libAliHLTRCU.so", 
  "libAliHLTTPC.so", 
  //  "libAliHLTSample.so",
  "libAliHLTCalo.so",
  "libAliHLTEMCAL.so",
  "libAliHLTPHOS.so",
  "libAliHLTMUON.so",
  "libAliHLTTRD.so",
  "libAliHLTITS.so",
  "libAliHLTVZERO.so",
  "libAliHLTZDC.so",
  "libAliHLTGlobal.so",
  "libAliHLTTrigger.so",
  NULL
};

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTSystem)

AliHLTSystem::AliHLTSystem(AliHLTComponentLogSeverity loglevel, const char* name,
			   AliHLTComponentHandler* pCompHandler,
			   AliHLTConfigurationHandler* pConfHandler
			   )
  : fpComponentHandler(pCompHandler==NULL?AliHLTComponentHandler::CreateHandler():pCompHandler)
  , fpConfigurationHandler(pConfHandler==NULL?AliHLTConfigurationHandler::CreateHandler():pConfHandler),
  fTaskList(),
  fState(0),
  fChains(),
  fStopwatches(new TObjArray),
  fEventCount(-1),
  fGoodEvents(-1),
  fpChainHandlers(NULL),
  fpEsdHandlers(NULL),
  fpProprietaryHandlers(NULL),
  fpHLTOUTTask(NULL),
  fpHLTOUT(NULL),
  fHLTOUTUse(0),
  fpControlTask(NULL),
  fName(name)
  , fECSParams()
  , fUseHLTOUTComponentTypeGlobal(true)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  if (fgNofInstances++>0) {
    // July 2008: multiple instances are now allowed
    // AliHLTSystem is used in multiple instances for the kChain HLTOUT handler
    //HLTWarning("multiple instances of AliHLTSystem, you should not use more than one at a time");
  }

  SetGlobalLoggingLevel(loglevel);
  SetFrameworkLog(loglevel);
  if (fpComponentHandler) {
    AliHLTAnalysisEnvironment env;
    memset(&env, 0, sizeof(AliHLTAnalysisEnvironment));
    env.fStructSize=sizeof(AliHLTAnalysisEnvironment);
    env.fAllocMemoryFunc=AliHLTSystem::AllocMemory;
    env.fGetEventDoneDataFunc=AliHLTSystem::AllocEventDoneData;
    env.fLoggingFunc=NULL;
    fpComponentHandler->SetEnvironment(&env);
    InitAliLogFunc(fpComponentHandler);
    if (fgNofInstances==1) {
    fpComponentHandler->AnnounceVersion();
    }
  } else {
    HLTFatal("can not create Component Handler");
  }
  if (fpConfigurationHandler==NULL) {
    HLTFatal("can not create Configuration Handler");
  }
}

AliHLTSystem::~AliHLTSystem()
{
  // see header file for class documentation
  fgNofInstances--;
  CleanupHLTOUTHandlers();
  CleanTaskList();
  if (fpConfigurationHandler) {
    fpConfigurationHandler->Destroy();
  }
  fpConfigurationHandler=NULL;
  
  if (fpComponentHandler) {
    fpComponentHandler->Destroy();
  }
  fpComponentHandler=NULL;
  delete fStopwatches;

  // note: fpHLTOUTTask and fpControlTask are deleted by
  // CleanTaskList
}

int AliHLTSystem::fgNofInstances=0;

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
	HLTError("configuration mismatch, there is already a task with configuration name \"%s\", but it is different. Most likely configuration %p is not registered properly", pConf->GetName(), pConf);
	iResult=-EEXIST;
      }
      // task for this configuration exists, terminate
      pTask=NULL;
    // check first if the configuration has all sources resolved, try to extract otherwise
    } else if (pConf->SourcesResolved()!=1 && pConf->ExtractSources()!=1) {
	HLTError("configuration \"%s\" has unresolved sources, aborting ...", pConf->GetName());
	iResult=-ENOLINK;
    } else {
      pTask=new AliHLTTask(pConf);
      if (pTask==NULL) {
	iResult=-ENOMEM;
      } else {
	pTask->SetLocalLoggingLevel(GetLocalLoggingLevel());
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
  fpHLTOUTTask=NULL;
  fpControlTask=NULL;
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
  if (fpControlTask==NULL) {
    fpControlTask=new AliHLTControlTask;
    if (!fpControlTask) return -ENOMEM;
    fTaskList.AddFirst(fpControlTask);
  }
  TObjLink *controlLnk=NULL;
  TObjLink *lnk = fTaskList.FirstLink();
  assert(!lnk || lnk->GetObject()==fpControlTask || fpControlTask==NULL);
  if (lnk && lnk->GetObject()==fpControlTask) {
    if (pTask->GetConf() && pTask->GetConf()->GetFirstSource()==NULL) {
      pTask->SetDependency(fpControlTask);
      fpControlTask->SetTarget(pTask);
    }
    controlLnk=lnk;
    lnk=lnk->Next();
  }
  if ((iResult=pTask->CheckDependencies())<=0)
    lnk=NULL;
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
      } else if (controlLnk) {
	fTaskList.AddAfter(controlLnk, pTask);
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

int AliHLTSystem::Run(Int_t iNofEvents, int bStop, AliHLTUInt64_t trgMask,
		      AliHLTUInt32_t timestamp, AliHLTUInt32_t eventtype,
		      AliHLTUInt32_t participatingDetectors)
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
	// Matthias Oct 11 2008 this is a bug
	// By resuming the stopwatches at this point, all continued counting, but the
	// starting and stopping is controlled by the AliHLTStopwatchGuard
	//ResumeBenchmarking(fStopwatches);    
      }
      for (int i=fEventCount; i<fEventCount+iNofEvents && iResult>=0; i++) {
	if (fpHLTOUTTask) {
	  if (iNofEvents>1 && i==fEventCount) {
	    HLTWarning("can not add more than one event to the HLTOUT, skipping all but last block");
	  }
	  // reset and prepare for new data
	  fpHLTOUTTask->Reset();
	}
	if (eventtype == 0) {
	  eventtype = gkAliEventTypeData;
	  participatingDetectors = 0x0;
	}
	if ((iResult=ProcessTasks(i, trgMask, timestamp, eventtype, participatingDetectors))>=0) {
	  fGoodEvents++;
	  iCount++;
	} else {
	  // TODO: define different running modes to either ignore errors in
	  // event processing or not
	  // currently ignored 
	  iResult=0;
	}
	AliHLTDataBuffer::SetGlobalEventCount(iCount);
      }
      fEventCount+=iNofEvents;
      if (bStop) StopTasks();
      else PauseBenchmarking(fStopwatches);
    }
    if (bStop) DeinitTasks();
  }
  if (iResult>=0) {
    iResult=iCount;
  } else  if (iResult==-126 /*ENOKEY*/) {
    iResult=0; // do not propagate the error
  }
  ClearStatusFlags(kRunning);
  AliHLTDataBuffer::PrintStatistics();
  return iResult;
}

int AliHLTSystem::InitTasks()
{
  // see header file for class documentation
  int iResult=0;
  TObjLink *lnk=fTaskList.FirstLink();

  if (lnk==NULL) {
    HLTInfo("Task list is empty, skipping HLT");
    return -126 /*ENOKEY*/;
  }
  while (lnk && iResult>=0) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      iResult=pTask->Init(NULL, fpComponentHandler);
//       ProcInfo_t ProcInfo;
//       gSystem->GetProcInfo(&ProcInfo);
//       HLTInfo("task %s initialized (%d), current memory usage %d %d", pTask->GetName(), iResult, ProcInfo.fMemResident, ProcInfo.fMemVirtual);
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

int AliHLTSystem::PauseBenchmarking(TObjArray* pStopwatches) const
{
  // see header file for class documentation
  if (pStopwatches==NULL) return 0;

  for (int i=0; i<(int)AliHLTComponent::kSWTypeCount; i++) {
    if (!pStopwatches->At(i)) continue;
    TStopwatch* pSw=dynamic_cast<TStopwatch*>(pStopwatches->At(i));
    if (pSw) pSw->Stop();
  }
  return 0;
}

int AliHLTSystem::ResumeBenchmarking(TObjArray* pStopwatches) const
{
  // see header file for class documentation
  if (pStopwatches==NULL) return 0;

  for (int i=0; i<(int)AliHLTComponent::kSWTypeCount; i++) {
    if (!pStopwatches->At(i)) continue;
    TStopwatch* pSw=dynamic_cast<TStopwatch*>(pStopwatches->At(i));
    if (pSw) pSw->Continue();
  }
  return 0;
}

int AliHLTSystem::PrintBenchmarking(TObjArray* pStopwatches, int bClean) const
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
    HLTImportant("HLT statistics:\n"
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
//       ProcInfo_t ProcInfo;
//       gSystem->GetProcInfo(&ProcInfo);
//       HLTInfo("task %s started (%d), current memory usage %d %d", pTask->GetName(), iResult, ProcInfo.fMemResident, ProcInfo.fMemVirtual);
    } else {
    }
    lnk = lnk->Next();
  }
  if (iResult<0) {
    HLTError("can not start task list, error %d", iResult);
  } else {
    SetStatusFlags(kStarted);
    fEventCount=0;
    fGoodEvents=0;
    if ((iResult=SendControlEvent(kAliHLTDataTypeSOR))<0) {
      HLTError("can not send SOR event: error %d", iResult);
    }
  }
  return iResult;
}

int AliHLTSystem::ProcessTasks(Int_t eventNo, AliHLTUInt64_t trgMask,
	  AliHLTUInt32_t timestamp, AliHLTUInt32_t eventtype,
	  AliHLTUInt32_t participatingDetectors)
{
  // see header file for class documentation
  int iResult=0;
  HLTDebug("processing event no %d", eventNo);
  TObjLink *lnk=fTaskList.FirstLink();
  while (lnk) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      if (iResult>=0) {
      iResult=pTask->ProcessTask(eventNo, eventtype, trgMask, timestamp, participatingDetectors);
//       ProcInfo_t ProcInfo;
//       gSystem->GetProcInfo(&ProcInfo);
//       HLTInfo("task %s processed (%d), current memory usage %d %d", pTask->GetName(), iResult, ProcInfo.fMemResident, ProcInfo.fMemVirtual);
      } else {
	pTask->SubscribeSourcesAndSkip();
      }
    } else {
    }
    lnk = lnk->Next();
  }

  if (iResult>=0) {
    HLTImportant("Event %d successfully finished (%d)", eventNo, iResult);
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
  if ((iResult=SendControlEvent(kAliHLTDataTypeEOR))<0) {
    HLTError("can not send EOR event");
  }

  // cleanup blocks from the last event. This is a bit awkward. All output
  // blocks from the chains need to be stored in the HLTOUT task. Though,
  // we do not know, whether HLTOUT is going to be processed or not.
  if (fpHLTOUTTask)
    fpHLTOUTTask->Reset();

  TObjLink *lnk=fTaskList.FirstLink();
  while (lnk) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      int locResult=pTask->EndRun();
      if (iResult>=0 && locResult<0) iResult=locResult;
//       ProcInfo_t ProcInfo;
//       gSystem->GetProcInfo(&ProcInfo);
//       HLTInfo("task %s stopped (%d), current memory usage %d %d", pTask->GetName(), iResult, ProcInfo.fMemResident, ProcInfo.fMemVirtual);
    } else {
    }
    lnk = lnk->Next();
  }
  PrintBenchmarking(fStopwatches, 1 /*clean*/);
  if (fEventCount!=fGoodEvents) {
    HLTError("%d out of %d event(s) failed", fEventCount-fGoodEvents, fEventCount);
  }
  ClearStatusFlags(kStarted);
  return iResult;
}

int AliHLTSystem::SendControlEvent(AliHLTComponentDataType dt)
{
  // see header file for class documentation
  int iResult=0;

  AliHLTComponentBlockDataList controlBlocks;
  AliHLTComponentBlockData bd;

  // run decriptor block of type kAliHLTDataTypeSOR/kAliHLTDataTypeEOR 
  AliHLTComponent::FillBlockData(bd);
  AliHLTRunDesc runDesc;
  memset(&runDesc, 0, sizeof(AliHLTRunDesc));
  runDesc.fStructSize=sizeof(AliHLTRunDesc);
  runDesc.fRunNo=AliHLTMisc::Instance().GetCDBRunNo();
  bd.fPtr=&runDesc;
  bd.fSize=sizeof(AliHLTRunDesc);
  bd.fDataType=dt;
  bd.fSpecification=kAliHLTVoidDataSpec;
  controlBlocks.push_back(bd);

  // ECS parameter of type kAliHLTDataTypeECSParam
  if (fECSParams.IsNull())
    fECSParams="CTP_TRIGGER_CLASS=00:DUMMY-TRIGGER-ALL:00-01-02-03-04-05-06-07-08-09-10-11-12-13-14-15-16-17";
  AliHLTComponent::FillBlockData(bd);
  bd.fPtr=(void*)fECSParams.Data();
  bd.fSize=fECSParams.Length()+1;
  bd.fDataType=kAliHLTDataTypeECSParam;
  bd.fSpecification=kAliHLTVoidDataSpec;
  controlBlocks.push_back(bd);  

  AliHLTControlTask::AliHLTControlEventGuard g(fpControlTask, controlBlocks);
  HLTDebug("sending event %s, run descriptor %p", AliHLTComponent::DataType2Text(dt).c_str(), &runDesc);
  TObjLink *lnk=fTaskList.FirstLink();
  while (lnk && iResult>=0) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      AliHLTUInt32_t eventType=gkAliEventTypeUnknown;
      if (dt==kAliHLTDataTypeSOR) eventType=gkAliEventTypeStartOfRun;
      else if (dt==kAliHLTDataTypeEOR) eventType=gkAliEventTypeEndOfRun;
      else HLTWarning("unknown control event %s", AliHLTComponent::DataType2Text(dt).c_str());
      iResult=pTask->ProcessTask(-1, eventType, 0, 0);
    } else {
    }
    lnk = lnk->Next();
  }

  // control events are not supposed to go into the HLTOUT
  if (fpHLTOUTTask)
    fpHLTOUTTask->Reset();

  HLTDebug("event %s done (%d)", AliHLTComponent::DataType2Text(dt).c_str(), iResult);
  return iResult;
}

int AliHLTSystem::DeinitTasks()
{
  // see header file for class documentation
  int iResult=0;
  TObjLink *lnk=fTaskList.LastLink();
  while (lnk) {
    TObject* obj=lnk->GetObject();
    if (obj) {
      AliHLTTask* pTask=(AliHLTTask*)obj;
      int localRes=pTask->Deinit();
      if (iResult>=0) iResult=localRes;
//       ProcInfo_t ProcInfo;
//       gSystem->GetProcInfo(&ProcInfo);
//       HLTInfo("task %s cleaned (%d), current memory usage %d %d", pTask->GetName(), iResult, ProcInfo.fMemResident, ProcInfo.fMemVirtual);
    } else {
    }
    lnk = lnk->Prev();
  }
  fEventCount=-1;
  fGoodEvents=-1;

  return iResult;
}

int AliHLTSystem::CleanupHLTOUTHandlers()
{
  // see header file for class documentation
  if (fpChainHandlers) {
    AliHLTOUT::AliHLTOUTHandlerListEntryVector* pHandlers=reinterpret_cast<AliHLTOUT::AliHLTOUTHandlerListEntryVector*>(fpChainHandlers);
    fpChainHandlers=NULL;
    if (pHandlers) {
      AliHLTOUT::InvalidateBlocks(*pHandlers);
      AliHLTOUT::RemoveEmptyDuplicateHandlers(*pHandlers);
    }
    assert(pHandlers->size()==0);
    delete pHandlers;
  }

  if (fpEsdHandlers) {
    AliHLTOUT::AliHLTOUTHandlerListEntryVector* pHandlers=reinterpret_cast<AliHLTOUT::AliHLTOUTHandlerListEntryVector*>(fpEsdHandlers);
    fpEsdHandlers=NULL;
    if (pHandlers) {
      AliHLTOUT::InvalidateBlocks(*pHandlers);
      AliHLTOUT::RemoveEmptyDuplicateHandlers(*pHandlers);
    }
    assert(pHandlers->size()==0);
    delete pHandlers;
  }

  if (fpProprietaryHandlers) {
    AliHLTOUT::AliHLTOUTHandlerListEntryVector* pHandlers=reinterpret_cast<AliHLTOUT::AliHLTOUTHandlerListEntryVector*>(fpProprietaryHandlers);
    fpProprietaryHandlers=NULL;
    if (pHandlers) {
      AliHLTOUT::InvalidateBlocks(*pHandlers);
      AliHLTOUT::RemoveEmptyDuplicateHandlers(*pHandlers);
    }
    assert(pHandlers->size()==0);
    delete pHandlers;
  }
  return 0;
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

int AliHLTSystem::AllocEventDoneData( void* /*param*/, AliHLTEventID_t /*eventID*/, unsigned long size, AliHLTComponentEventDoneData** edd )
{
  // see header file for class documentation
  unsigned long blocksize=sizeof(AliHLTComponentEventDoneData)+size;
  void* block=AllocMemory(NULL, blocksize);
  if (!block) return -ENOMEM;
  memset(block, 0, blocksize);
  *edd=reinterpret_cast<AliHLTComponentEventDoneData*>(block);
  (*edd)->fStructSize=sizeof(AliHLTComponentEventDoneData);
  (*edd)->fDataSize=size;
  (*edd)->fData=reinterpret_cast<AliHLTUInt8_t*>(block)+sizeof(AliHLTComponentEventDoneData);
  
  return 0;
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
	CleanupHLTOUTHandlers();
	}
      } else {
      if ((iResult=AliHLTOfflineInterface::SetParamsToComponents(runLoader, rawReader))>=0) {
	AliHLTUInt64_t trgMask=0x1;
	AliHLTUInt32_t timestamp=0;
	AliHLTUInt32_t eventtype=0;
	if (runLoader==NULL) {
	  // this is a quick workaround for the case of simulation
	  // the trigger framework is still under development, secondly, AliHLTSimulation
	  // does not yet add the emulated ECS parameters, so no CTP trigger is known in the HLT
	  // AliHLTTask will initialize one dummy CTP trigger class with bit 0, that's why the
	  // default trigger mask is 0x1
	  trgMask=AliHLTMisc::Instance().GetTriggerMask(rawReader);

	  // get the timestamp and type of the event from the raw reader
	  // this is currently only meaningfull for reconstruction (runloader==NULL)
	  timestamp=AliHLTMisc::Instance().GetTimeStamp(rawReader);
	  eventtype=AliHLTMisc::Instance().GetEventType(rawReader);
	}
	// the system always remains started after event processing, a specific
	// call with nofEvents==0 is needed to execute the stop sequence
	if ((iResult=Run(nofEvents, 0, trgMask, timestamp, eventtype))<0) SetStatusFlags(kError);
      }
      }

      // add the current HLTOUT task to the collection
      if (fpHLTOUTTask) {
	AliHLTOUT* pTask=dynamic_cast<AliHLTOUT*>(fpHLTOUTTask);
	if (pTask && (iResult=pTask->Init())>=0) {
	  if (pTask->GetNofDataBlocks()>0) {
	    AliHLTOUT* pHLTOUT=RequestHLTOUT();
	    if (pHLTOUT) {
	      pHLTOUT->AddSubCollection(pTask);
	      ReleaseHLTOUT(pHLTOUT);
	    } else {
	      HLTWarning("no HLTOUT instance available, output blocks of the chain are ignored");
	    }
	  }
	} else {
	  HLTWarning("can not initialize HLTOUT sub collection %s for reconstruction chain (%d), data blocks are lost", pTask?fpHLTOUTTask->GetName():"nil", iResult);
	  iResult=0;
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

int AliHLTSystem::ProcessHLTOUT(AliHLTOUT* pHLTOUT, AliESDEvent* esd)
{
  // see header file for class documentation
  int iResult=0;
  if (!pHLTOUT) return -EINVAL;
  HLTDebug("processing %d HLT data blocks", pHLTOUT->GetNofDataBlocks());

  //
  // process all kChain handlers first
  //
  if ((iResult=ProcessHLTOUTkChain(pHLTOUT))<0) {
    HLTWarning("Processing of kChain-type data blocks failed with error code %d", iResult);
    iResult=0;
  } 

  if (!fpEsdHandlers)
    fpEsdHandlers=new AliHLTOUT::AliHLTOUTHandlerListEntryVector;
  if (!fpProprietaryHandlers)
    fpProprietaryHandlers=new AliHLTOUT::AliHLTOUTHandlerListEntryVector;

  AliHLTOUT::AliHLTOUTHandlerListEntryVector* pEsdHandlers=reinterpret_cast<AliHLTOUT::AliHLTOUTHandlerListEntryVector*>(fpEsdHandlers);
  AliHLTOUT::AliHLTOUTHandlerListEntryVector* pProprietaryHandlers=reinterpret_cast<AliHLTOUT::AliHLTOUTHandlerListEntryVector*>(fpProprietaryHandlers);
  if (!pEsdHandlers || !pProprietaryHandlers) return -ENOMEM;

  // invalidate all blocks
  AliHLTOUT::InvalidateBlocks(*pEsdHandlers);
  AliHLTOUT::InvalidateBlocks(*pProprietaryHandlers);

  AliHLTComponentDataTypeList esdBlocks;

  for (iResult=pHLTOUT->SelectFirstDataBlock();
       iResult>=0;
       iResult=pHLTOUT->SelectNextDataBlock()) {
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    pHLTOUT->GetDataBlockDescription(dt, spec);
    AliHLTOUTHandler* pHandler=pHLTOUT->GetHandler();
    AliHLTModuleAgent::AliHLTOUTHandlerType handlerType=pHLTOUT->GetDataBlockHandlerType();

    // default handling for ESD data blocks does not require an explicite handler
    if (!pHandler && (dt==kAliHLTDataTypeESDObject || dt==kAliHLTDataTypeESDTree)) {
      handlerType=AliHLTModuleAgent::kEsd;
    }
    const char* pMsg="invalid";
    switch (handlerType) {
    case AliHLTModuleAgent::kEsd:
      {
	if (pHandler) {
	  // schedule for later processing
	  pHLTOUT->InsertHandler(*pEsdHandlers, pHLTOUT->GetDataBlockHandlerDesc());
	} else {
	  AliHLTComponentDataTypeList::iterator element=esdBlocks.begin();
	  for (; element!=esdBlocks.end(); element++) {
	    if (*element==dt) {
	      HLTWarning("multiple ESDs of identical data type %s, please add appropriate handler to merge ESDs", AliHLTComponent::DataType2Text(dt).c_str());
	      break;
	    }
	  }
	  if (element==esdBlocks.end()) esdBlocks.push_back(dt);

	  // write directly
	  const AliHLTUInt8_t* pBuffer=NULL;
	  AliHLTUInt32_t size=0;
	  if (pHLTOUT->GetDataBuffer(pBuffer, size)>=0) {
	    pHLTOUT->WriteESD(pBuffer, size, dt, esd);
	    pHLTOUT->ReleaseDataBuffer(pBuffer);
	  }
	  pHLTOUT->MarkDataBlockProcessed();
	}
      }
      break;
    case AliHLTModuleAgent::kRawReader:
      // handled in the AliRawReaderHLT
      break;
    case AliHLTModuleAgent::kRawStream:
      HLTWarning("HLTOUT handler type 'kRawStream' not yet implemented: agent %s, data type %s, specification %#x",
		 pHLTOUT->GetAgent()?pHLTOUT->GetAgent()->GetModuleId():"<invalid>",
		 AliHLTComponent::DataType2Text(dt).c_str(), spec);
      break;
    case AliHLTModuleAgent::kChain:
      HLTWarning("HLTOUT handler type 'kChain' has already been processed: agent %s, data type %s, specification %#x\n"
		 "New block of this type added by the chain? Skipping data block ...",
		 pHLTOUT->GetAgent()?pHLTOUT->GetAgent()->GetModuleId():"<invalid>",
		 AliHLTComponent::DataType2Text(dt).c_str(), spec);
      break;
    case AliHLTModuleAgent::kProprietary:
      HLTDebug("processing proprietary data: agent %s, data type %s, specification %#x",
		 pHLTOUT->GetAgent()?pHLTOUT->GetAgent()->GetModuleId():"<invalid>",
		 AliHLTComponent::DataType2Text(dt).c_str(), spec);
      if (pHandler) {
	AliHLTOUT::AliHLTOUTLockGuard g(pHLTOUT);
	int res=pHandler->ProcessData(pHLTOUT);
	if (res<0) {
	  HLTWarning("processing proprietary data failed (%d): agent %s, data type %s, specification %#x",
		     res, pHLTOUT->GetAgent()?pHLTOUT->GetAgent()->GetModuleId():"<invalid>",
		     AliHLTComponent::DataType2Text(dt).c_str(), spec);
	}
      }
      break;
    case AliHLTModuleAgent::kUnknownOutput:
      pMsg="unknown";
      // fall trough intended
    default:
      HLTWarning("%s handler type: agent %s, data type %s, specification %#x, ... skipping data block",
		 pMsg, pHLTOUT->GetAgent()?pHLTOUT->GetAgent()->GetModuleId():"<invalid>",
		 AliHLTComponent::DataType2Text(dt).c_str(), spec);
    }
  }
  // TODO: the return value of SelectFirst/NextDataBlock must be
  // changed in order to avoid this check
  if (iResult==-ENOENT) iResult=0;

  AliHLTOUT::AliHLTOUTHandlerListEntryVector::iterator handler;

  // process and write all esd data blocks
  for (handler=pEsdHandlers->begin(); handler!=pEsdHandlers->end() && iResult>=0; handler++) {
    AliHLTOUT::AliHLTOUTSelectionGuard g(pHLTOUT, &(*handler));	    
    AliHLTOUTHandler* pHandler=*handler;
    const AliHLTUInt8_t* pBuffer=NULL;
    AliHLTUInt32_t size=0;
    pHandler->ProcessData(pHLTOUT);
    if ((size=pHandler->GetProcessedData(pBuffer))>0) {
      AliHLTModuleAgent::AliHLTOUTHandlerDesc desc=*handler;
      AliHLTComponentDataType dt=desc;
      pHLTOUT->WriteESD(pBuffer, size, dt, esd);
      pHandler->ReleaseProcessedData(pBuffer, size);
    }
    pHLTOUT->MarkDataBlocksProcessed(&(*handler));
  }

  // process all kProprietary data blocks
  for (handler=pProprietaryHandlers->begin(); handler!=pProprietaryHandlers->end() && iResult>=0; handler++) {
    AliHLTOUT::AliHLTOUTSelectionGuard g(pHLTOUT, &(*handler));	    
    AliHLTOUTHandler* pHandler=*handler;
    const AliHLTUInt8_t* pBuffer=NULL;
    AliHLTUInt32_t size=0;
    pHandler->ProcessData(pHLTOUT);
    if ((size=pHandler->GetProcessedData(pBuffer))>0) {
      HLTWarning("data produced by kProprietary handler ignored");
      pHandler->ReleaseProcessedData(pBuffer, size);
    }
    pHLTOUT->MarkDataBlocksProcessed(&(*handler));
  }

  // remove all empty handlers form the list (handlers which did not get a block this time)
  AliHLTOUT::RemoveEmptyDuplicateHandlers(*pEsdHandlers);
  AliHLTOUT::RemoveEmptyDuplicateHandlers(*pProprietaryHandlers);

  return iResult;
}

int AliHLTSystem::ProcessHLTOUTkChain(AliHLTOUT* pHLTOUT)
{
  // see header file for class documentation
  int iResult=0;
  if (!pHLTOUT) return -EINVAL;

  if (!fpChainHandlers)
    fpChainHandlers=new AliHLTOUT::AliHLTOUTHandlerListEntryVector;

  AliHLTOUT::AliHLTOUTHandlerListEntryVector* pChainHandlers=reinterpret_cast<AliHLTOUT::AliHLTOUTHandlerListEntryVector*>(fpChainHandlers);
  if (!pChainHandlers) return -ENOMEM;

  // invalidate all blocks
  AliHLTOUT::InvalidateBlocks(*pChainHandlers);

  // fill the list
  pHLTOUT->FillHandlerList(*pChainHandlers, AliHLTModuleAgent::kChain);

  // process all defined chain handlers
  AliHLTOUT::AliHLTOUTHandlerListEntryVector::iterator chainHandler;
  for (chainHandler=pChainHandlers->begin(); chainHandler!=pChainHandlers->end() && iResult>=0; chainHandler++) {
    if (chainHandler->IsEmpty()) continue;
    AliHLTOUT::AliHLTOUTSelectionGuard g(pHLTOUT, &(*chainHandler));	    
    AliHLTOUTHandler* pHandler=*chainHandler;
    const AliHLTUInt8_t* pBuffer=NULL;
    AliHLTUInt32_t size=0;
    pHandler->ProcessData(pHLTOUT);
    if ((size=pHandler->GetProcessedData(pBuffer))>0) {
      AliHLTModuleAgent::AliHLTOUTHandlerDesc desc=*chainHandler;
      //AliHLTComponentDataType dt=desc;

      pHandler->ReleaseProcessedData(pBuffer, size);
    }
    pHLTOUT->MarkDataBlocksProcessed(&(*chainHandler));
  }

  // remove all empty handlers form the list (handlers which did not get a block this time)
  AliHLTOUT::RemoveEmptyDuplicateHandlers(*pChainHandlers);

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
	int iEntries=pTokens->GetEntriesFast();
	for (int i=0; i<iEntries && iResult>=0; i++) {
	  iResult=fpComponentHandler->LoadLibrary((((TObjString*)pTokens->At(i))->String()).Data());
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
  ClearStatusFlags(kTaskListCreated);
  if (CheckFilter(kHLTLogDebug))
    AliHLTModuleAgent::PrintStatus();
  if (CheckStatus(kConfigurationLoaded)==0) {
    iResult=LoadConfigurations(rawReader, runloader);
  } else {
    if (fChains.Length()==0) {
      HLTError("custom configuration(s) specified, but no configuration to run in local reconstruction, use \'chains=<chain,...>\' option");
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
    //AliHLTComponentHandler::TLibraryMode libMode=AliHLTComponentHandler::kDynamic;
    TString libs("");
    TString excludelibs("");
    TString alloptions(options);
    TObjArray* pTokens=alloptions.Tokenize(" ");
    if (pTokens) {
      int iEntries=pTokens->GetEntriesFast();
      for (int i=0; i<iEntries; i++) {
	TString token=(((TObjString*)pTokens->At(i))->String());
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
	} else if (token.Contains("frameworklog=")) {
	  TString param=token.ReplaceAll("frameworklog=", "");
	  if (param.IsDigit()) {
	    SetFrameworkLog((AliHLTComponentLogSeverity)param.Atoi());
	  } else if (param.BeginsWith("0x") &&
		     param.Replace(0,2,"",0).IsHex()) {
	    int severity=0;
	    sscanf(param.Data(),"%x", &severity);
	    SetFrameworkLog((AliHLTComponentLogSeverity)severity);
	  } else {
	    HLTWarning("wrong parameter for option \'loglevel=\', (hex) number expected");
	  }
	} else if (token.Contains("alilog=off")) {
	  SwitchAliLog(0);
	} else if (token.Contains("config=") || token.Contains("run-online-config")) {
	  if (!CheckStatus(kConfigurationLoaded)) {
	    Int_t error=0;
	    AliHLTOnlineConfiguration* pConf = NULL;
	    if (token.Contains("run-online-config")) {
	      AliCDBEntry* pEntry=AliHLTMisc::Instance().LoadOCDBEntry("HLT/Calib/OnlineConfig");
	      if (pEntry) {
		TObject* pObject=AliHLTMisc::Instance().ExtractObject(pEntry);
		if (pObject && pObject->IsA() == AliHLTOnlineConfiguration::Class())
		  pConf = (AliHLTOnlineConfiguration*)pObject;
	      }
	    }
	    if (token.Contains("config=")) {
	      TString param=token.ReplaceAll("config=", "");
	      if (token.EndsWith(".xml", TString::kIgnoreCase)) {
		Int_t filesize = 0;
		pConf = new AliHLTOnlineConfiguration;
		filesize = pConf->LoadConfiguration(param.Data());
		if (filesize <= 0) {
		  HLTError("cannot load config \'%s\'", param.Data());
		  iResult=-EBADF;
		}
	      } else {
		gROOT->Macro(param.Data(), &error);
		if (error==0) {
		  SetStatusFlags(kConfigurationLoaded);
		} else {
		  HLTError("cannot execute macro \'%s\'", param.Data());
		  iResult=-EBADF;
		}
	      }
	    }
	    if (pConf) {
		error = pConf->Parse();
		if (error==0) {
		  fChains = pConf->GetDefaultChains();
		  libs = pConf->GetComponentLibraries();
		  libs += " ";
		  SetStatusFlags(kConfigurationLoaded);
		} else {
		  HLTError("cannot parse online configuration");
		  iResult=-EBADF;
		}
	    }
	    delete pConf; pConf=NULL;
	  } else {
	    HLTWarning("HLT options has both a config file and run-online-config set");
	  }
	} else if (token.Contains("chains=")) {
	  TString param=token.ReplaceAll("chains=", "");
	  fChains=param.ReplaceAll(",", " ");
	  if (fChains.IsNull()) fChains=" "; // disable all chains
	} else if (token.Contains("libmode=")) {
	  TString param=token.ReplaceAll("libmode=", "");
	  param.ReplaceAll(",", " ");
	  if (fpComponentHandler) {
	    if (param.CompareTo("static")==0) {
	      fpComponentHandler->SetLibraryMode(AliHLTComponentHandler::kStatic);
	    } else if (param.CompareTo("dynamic")==0) {
	      fpComponentHandler->SetLibraryMode(AliHLTComponentHandler::kDynamic);
	    } else {
	      HLTWarning("wrong argument for option \'libmode=\', use \'static\' or \'dynamic\'");
	    }
	  }
	} else if (token.BeginsWith("ECS=")) {
	  fECSParams=token.ReplaceAll("ECS=", "");
	} else if (token.BeginsWith("hltout-mode=")) {
	  // The actual parameter for argument 'hltout-mode' is treated in AliSimulation.
	  // For AliHLTSystem the occurrence with parameter 'split' signals the use of the
	  // separated HLTOUTComponents for digit and raw data. All others indicate
	  // HLTOUTComponent type 'global' where the data generation is steered from global
	  // flags
	  fUseHLTOUTComponentTypeGlobal=token.CompareTo("hltout-mode=split")!=0;
	} else if (token.BeginsWith("lib") && token.EndsWith(".so")) {
	  libs+=token;
	  libs+=" ";
	} else if (token.BeginsWith("!lib") && token.EndsWith(".so")) {
	  excludelibs+=token;
	  excludelibs+=" ";
	} else {
	  HLTWarning("unknown option \'%s\'", token.Data());
	}
      }
      delete pTokens;
    }

    if (iResult>=0) {
      if (libs.IsNull()) {
	const char** deflib=fgkHLTDefaultLibs;
	for (;*deflib; deflib++) {
	  if (excludelibs.Contains(*deflib)) continue;
	  libs+=*deflib;
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
  AliHLTModuleAgent* pAgent=NULL;

  // first check for the required libraries and load those
  TString extralibs;
  for (pAgent=AliHLTModuleAgent::GetFirstAgent(); 
       pAgent && iResult>=0;
       pAgent=AliHLTModuleAgent::GetNextAgent()) {
    const char* deplibs=pAgent->GetRequiredComponentLibraries();
    if (deplibs) {
      HLTDebug("required libraries \'%s\' for agent %s (%p)", deplibs, pAgent->GetName(), pAgent);
      extralibs+=" ";
      extralibs+=deplibs;
    }
  }
  if (iResult>=0) {
    iResult=LoadComponentLibraries(extralibs.Data());
  }

  // in order to register the configurations in the correct sequence
  // all agents need to be ordered with respect to the required
  // libraries. Ordering relies on the naming convention
  // libAliHLT<Module>.so
  TList agents;
  for (pAgent=AliHLTModuleAgent::GetFirstAgent(); 
       pAgent && iResult>=0;
       pAgent=AliHLTModuleAgent::GetNextAgent()) {
    AliHLTModuleAgent* pPrevDep=NULL;
    TString dependencies=pAgent->GetRequiredComponentLibraries();
    TObjArray* pTokens=dependencies.Tokenize(" ");
    if (pTokens) {
      for (int n=0; n<pTokens->GetEntriesFast(); n++) {
	TString module=((TObjString*)pTokens->At(n))->String();
	HLTDebug("  checking %s", module.Data());
	module.ReplaceAll("libAliHLT", "");
	module.ReplaceAll(".so", "");
	
	for (AliHLTModuleAgent* pCurrent=dynamic_cast<AliHLTModuleAgent*>(pPrevDep==NULL?agents.First():agents.After(pPrevDep));
	     pCurrent!=NULL; pCurrent=dynamic_cast<AliHLTModuleAgent*>(agents.After(pCurrent))) {
	  HLTDebug("    checking %s == %s", module.Data(), pCurrent->GetModuleId());

	  if (module.CompareTo(pCurrent->GetModuleId())==0) {
	    pPrevDep=pCurrent;
	    break;
	  }
	}
      }
      delete pTokens;
    }

    if (pPrevDep) {
      // insert right after the last dependency
      agents.AddAfter(pPrevDep, pAgent);
      HLTDebug("insert %s after %s", pAgent->GetModuleId(), pPrevDep->GetModuleId());
    } else {
      // insert at the beginning
      agents.AddFirst(pAgent);
      HLTDebug("insert %s at beginning", pAgent->GetModuleId());
    }
  }

  // now we load the configurations
  if (agents.GetEntries()) {
    TIter next(&agents);
    while ((pAgent = dynamic_cast<AliHLTModuleAgent*>(next()))) {
      HLTDebug("load configurations for agent %s (%p)", pAgent->GetName(), pAgent);
      pAgent->CreateConfigurations(fpConfigurationHandler, rawReader, runloader);
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
    HLTImportant("custom reconstruction chain: %s", chains.Data());
  } else {
    for (AliHLTModuleAgent* pAgent=AliHLTModuleAgent::GetFirstAgent();
	 pAgent && iResult>=0;
	 pAgent=AliHLTModuleAgent::GetNextAgent()) {
      const char* agentchains=pAgent->GetReconstructionChains(rawReader, runloader);
      if (agentchains) {
	if (!chains.IsNull()) chains+=" ";
	chains+=agentchains;
	HLTInfo("reconstruction chains for agent %s (%p): %s", pAgent->GetName(), pAgent, agentchains);
      }
    }
  }

  // build task list for chains
  TObjArray* pTokens=chains.Tokenize(" ");
  if (pTokens) {
    int iEntries=pTokens->GetEntriesFast();
    for (int i=0; i<iEntries && iResult>=0; i++) {
      const char* pCID=((TObjString*)pTokens->At(i))->String().Data();
      AliHLTConfiguration* pConf=fpConfigurationHandler->FindConfiguration(pCID);
      if (pConf) {
	iResult=BuildTaskList(pConf);
	if (true) { // condition was deprecated but kept for sake of svn diff
	  // bHaveOutput variable has to be set for both running modes
	  // AliHLTSimulation and AliHLTReconstruction
	  assert(fpComponentHandler!=NULL);
	  TString cid=pConf->GetComponentID();
	  if (runloader!=NULL && cid.CompareTo("HLTOUT")==0) {
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

      // add the SchemaEvolutionComponent which analyzes the ROOT objects in
      // the output stream
      if (fpComponentHandler->FindComponentIndex("ROOTSchemaEvolutionComponent")>=0 ||
	  fpComponentHandler->LoadLibrary("libAliHLTUtil.so")>=0) {
	AliHLTConfiguration schemaevo("_schemaevolution_", "ROOTSchemaEvolutionComponent", 
				      chains.Data(), "-file=HLT.StreamerInfo.root");
	iResult=BuildTaskList("_schemaevolution_");
      } else {
	HLTWarning("can not load libAliHLTUtil.so and ROOTSchemaEvolutionComponent");
      }

      // add the HLTOUT component
      if (fpComponentHandler->FindComponentIndex("HLTOUT")>=0 ||
	  fpComponentHandler->LoadLibrary("libHLTsim.so")>=0) {
	// for the default HLTOUTComponent type 'global' the data generation is steered
	// by global flags from AliSimulation. This allows for emulation of the old
	// AliHLTSimulation behavior where only one chain is run on either digits or
	// simulated raw data and the HLT digits and raw files have been generated
	// depending on the configuration
	const char* HLTOUTComponentId="HLTOUT";
	if (!fUseHLTOUTComponentTypeGlobal) {
	  // choose the type of output depending  on the availability of
	  // the raw reader
	  if (rawReader) HLTOUTComponentId="HLTOUTraw";
	  else HLTOUTComponentId="HLTOUTdigits";
	}
	AliHLTConfiguration globalout("_globalout_", HLTOUTComponentId, chains.Data(), NULL);
	iResult=BuildTaskList("_globalout_");
      } else {
	HLTError("can not load libHLTsim.so and HLTOUT component");
	iResult=-EFAULT;
      }
    }
  }

  // build HLTOUT task for reconstruction
  // Matthias 08.07.2008 the rawReader is never set when running embedded into
  // AliReconstruction. The system is configured during AliHLTReconstructor::Init
  // where the RawReader is not available. It is available in the first invocation
  // of Reconstruct.
  // 
  // That means that policy is slightly changed:
  // - if the run loader is available -> AliSimulation
  // - no run loader available -> AliReconstruction
  if (iResult>=0 && !runloader) {
    if (bHaveOutput) {
      // there are components in the chain which produce data which need to be
      // piped to an HLTOUT sub-collection
      if (!fpHLTOUTTask) {
	iResult=AddHLTOUTTask(chains.Data());
      }
    }
  }

  if (iResult>=0) SetStatusFlags(kTaskListCreated);

  return iResult;
}

int AliHLTSystem::AddHLTOUTTask(const char* hltoutchains)
{
  // see header file for class documentation
  int iResult=0;
  if (!hltoutchains || hltoutchains[0]==0) return 0;

  // check chains for output
  TString chains=hltoutchains;
  TObjArray* pTokens=chains.Tokenize(" ");
  if (pTokens) {
    int iEntries=pTokens->GetEntriesFast();
    for (int i=0; i<iEntries && iResult>=0; i++) {
      const char* token=((TObjString*)pTokens->At(i))->String().Data();
      AliHLTConfiguration* pConf=fpConfigurationHandler->FindConfiguration(token);
      if (pConf) {
	TString cid=pConf->GetComponentID();
	if (fpComponentHandler->HasOutputData(cid.Data())) {
	  continue;
	}
      } else {
	HLTWarning("can not find configuration %s", token);
      }
      // remove from the list of hltout chains
      chains.ReplaceAll(token, "");
    }
    delete pTokens;
  }

  // do not create the HLTOUT task if none of the chains have output
  if (chains.IsNull()) return 0;

  // indicate the task to be available
  iResult=1;

  if (fpHLTOUTTask) {
    if (strcmp(chains.Data(), fpHLTOUTTask->GetSourceChains())==0) {
      HLTWarning("HLTOUT task already added for chains \"%s\" %p", chains.Data(), fpHLTOUTTask);
    } else {
      HLTError("HLTOUT task already added for chains \"%s\" %p, ignoring new chains  \"%s\"",
	       fpHLTOUTTask->GetSourceChains(), fpHLTOUTTask, chains.Data());
    }
    return iResult;
  }

  fpHLTOUTTask=new AliHLTOUTTask(chains);
  if (fpHLTOUTTask) {
    if (fpHLTOUTTask->GetConf() && 
	(fpHLTOUTTask->GetConf()->SourcesResolved()>0 ||
	 fpHLTOUTTask->GetConf()->ExtractSources()>0)) {
      iResult=InsertTask(fpHLTOUTTask);
    } else {
      HLTError("HLTOUT task (%s) sources not resolved", fpHLTOUTTask->GetName());
      iResult=-ENOENT;
    }

    if (iResult<0) {
      delete fpHLTOUTTask;
    }

  } else {
    iResult=-ENOMEM;
  }
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

AliHLTfctVoid AliHLTSystem::FindDynamicSymbol(const char* library, const char* symbol)
{
  // see header file for class documentation
  if (fpComponentHandler==NULL) return NULL;
  return fpComponentHandler->FindSymbol(library, symbol);
}

void AliHLTSystem::SetFrameworkLog(AliHLTComponentLogSeverity level) 
{
  // see header file for class documentation
  SetLocalLoggingLevel(level);
  if (fpComponentHandler) fpComponentHandler->SetLocalLoggingLevel(level);
  if (fpConfigurationHandler) fpConfigurationHandler->SetLocalLoggingLevel(level);
}

int AliHLTSystem::LoggingVarargs(AliHLTComponentLogSeverity severity, 
				 const char* originClass, const char* originFunc,
				 const char* file, int line, ... ) const
{
  // see header file for function documentation
  int iResult=0;

  va_list args;
  va_start(args, line);

  if (!fName.IsNull())
    AliHLTLogging::SetLogString(this, " (%p)", "%s_pfmt_: ", fName.Data());
  iResult=SendMessage(severity, originClass, originFunc, file, line, AliHLTLogging::BuildLogString(NULL, args, !fName.IsNull() /*append if non empty*/));
  va_end(args);

  return iResult;
}

int AliHLTSystem::InitHLTOUT(AliHLTOUT* instance)
{
  // Init the HLTOUT instance for the current event.
  // The instance can be used by other classes to get hold on the data
  // from HLTOUT.
  if (!instance) return -EINVAL;
  if (fpHLTOUT && fpHLTOUT!=instance) return -EBUSY;
  fpHLTOUT=instance;
  return 0;
}

int AliHLTSystem::InvalidateHLTOUT(AliHLTOUT** target)
{
  // Clear the HLTOUT instance.
  int iResult=0;
  if (fHLTOUTUse>0) {
    HLTWarning("HLTOUT instance still in use, potential problem due to invalid pointer ahead");
    fHLTOUTUse=0;
    iResult=-EBUSY;
  }
  if (target) *target=fpHLTOUT;
  fpHLTOUT=NULL;
  return iResult;
}

AliHLTOUT* AliHLTSystem::RequestHLTOUT()
{
  // Get the HLTOUT instance.
  // User method for processing classes. To be released after use.
  if (!fpHLTOUT) return NULL;
  fHLTOUTUse++;
  return fpHLTOUT;
}

int AliHLTSystem::ReleaseHLTOUT(const AliHLTOUT* instance)
{
  // Release the HLTOUT instance after use.
  if (!instance) return -EINVAL;
  if (instance!=fpHLTOUT) return -ENOENT;
  fHLTOUTUse--;
  return 0;
}
