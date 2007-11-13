// $Id$
// splitted from AliHLTConfiguration.cxx,v 1.25 2007/10/12 13:24:47
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

/** @file   AliHLTTask.cxx
    @author Matthias Richter
    @date   
    @brief  Implementation of HLT tasks.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

#include <cerrno>
#include <cassert>
#include <iostream>
#include <string>
#include "AliHLTTask.h"
#include "AliHLTConfiguration.h"
#include "AliHLTComponent.h"
#include "AliHLTComponentHandler.h"
#include "TList.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTask)

AliHLTTask::AliHLTTask()
  :
  fpConfiguration(NULL),
  fpComponent(NULL),
  fpDataBuffer(NULL),
  fListTargets(),
  fListDependencies(),
  fBlockDataArray()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTask::AliHLTTask(AliHLTConfiguration* pConf)
  :
  fpConfiguration(pConf),
  fpComponent(NULL),
  fpDataBuffer(NULL),
  fListTargets(),
  fListDependencies(),
  fBlockDataArray()
{
  // see header file for function documentation
}

AliHLTTask::~AliHLTTask()
{
  // see header file for function documentation
  TObjLink* lnk=fListDependencies.FirstLink();

  while (lnk!=NULL) {
    AliHLTTask* pTask=(AliHLTTask*)lnk->GetObject();
    pTask->UnsetTarget(this);
    lnk=lnk->Next();
  }
  lnk=fListTargets.FirstLink();

  while (lnk!=NULL) {
    AliHLTTask* pTask=(AliHLTTask*)lnk->GetObject();
    pTask->UnsetDependency(this);
    lnk=lnk->Next();
  }

  if (fpComponent) delete fpComponent;
  fpComponent=NULL;
}

int AliHLTTask::Init(AliHLTConfiguration* pConf, AliHLTComponentHandler* pCH)
{
  // see header file for function documentation
  int iResult=0;
  if (fpConfiguration!=NULL && pConf!=NULL && fpConfiguration!=pConf) {
    HLTWarning("overriding existing reference to configuration object %p (%s) by %p",
	       fpConfiguration, GetName(), pConf);
  }
  if (pConf!=NULL) fpConfiguration=pConf;
  if (fpConfiguration) {
    if (pCH) {
      int argc=0;
      const char** argv=NULL;
      if ((iResult=fpConfiguration->GetArguments(&argv))>=0) {
	argc=iResult; // just to make it clear
	// TODO: we have to think about the optional environment parameter,
	// currently just set to NULL. 
	iResult=pCH->CreateComponent(fpConfiguration->GetComponentID(), NULL, argc, argv, fpComponent);
	if (fpComponent || iResult<=0) {
	  //HLTDebug("component %s (%p) created", fpComponent->GetComponentID(), fpComponent); 
	} else {
	  //HLTError("can not find component \"%s\" (%d)", fpConfiguration->GetComponentID(), iResult);
	}
      } else {
	HLTError("can not get argument list for configuration %s (%s)", fpConfiguration->GetName(), fpConfiguration->GetComponentID());
	iResult=-EINVAL;
      }
    } else {
      HLTError("component handler instance needed for task initialization");
      iResult=-EINVAL;
    }
  } else {
    HLTError("configuration object instance needed for task initialization");
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::Deinit()
{
  // see header file for function documentation
  int iResult=0;
  AliHLTComponent* pComponent=GetComponent();
  fpComponent=NULL;
  if (pComponent) {
    //HLTDebug("delete component %s (%p)", pComponent->GetComponentID(), pComponent); 
    pComponent->Deinit();
    delete pComponent;
  } else {
    HLTWarning("task %s (%p) doesn't seem to be in initialized", GetName(), this);
  }
  return iResult;
}

const char *AliHLTTask::GetName() const
{
  // see header file for function documentation
  if (fpConfiguration)
    return fpConfiguration->GetName();
  return TObject::GetName();
}

AliHLTConfiguration* AliHLTTask::GetConf() const
{
  // see header file for function documentation
  return fpConfiguration;
}

AliHLTComponent* AliHLTTask::GetComponent() const
{
  // see header file for function documentation
  return fpComponent;
}

AliHLTTask* AliHLTTask::FindDependency(const char* id)
{
  // see header file for function documentation
  AliHLTTask* pTask=NULL;
  if (id) {
    pTask=(AliHLTTask*)fListDependencies.FindObject(id);
  }
  return pTask;
}

int AliHLTTask::FollowDependency(const char* id, TList* pTgtList)
{
  // see header file for function documentation
  int iResult=0;
  if (id) {
    AliHLTTask* pDep=NULL;
    if ((pDep=(AliHLTTask*)fListDependencies.FindObject(id))!=NULL) {
      if (pTgtList) pTgtList->Add(pDep);
      iResult++;
    } else {
      TObjLink* lnk=fListDependencies.FirstLink();
      while (lnk && iResult==0) {
	pDep=(AliHLTTask*)lnk->GetObject();
	if (pDep) {
	  if ((iResult=pDep->FollowDependency(id, pTgtList))>0) {
	    if (pTgtList) pTgtList->AddFirst(pDep);
	    iResult++;
	  }
	} else {
	  iResult=-EFAULT;
	}
	lnk=lnk->Next();
      }
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTTask::PrintDependencyTree(const char* id, int bFromConfiguration)
{
  // see header file for function documentation
  HLTLogKeyword("task dependencies");
  int iResult=0;
  TList tgtList;
  if (bFromConfiguration) {
    if (fpConfiguration)
      iResult=fpConfiguration->FollowDependency(id, &tgtList);
    else
      iResult=-EFAULT;
  } else
    iResult=FollowDependency(id, &tgtList);
  if (iResult>0) {
    HLTMessage("     task \"%s\": dependency level %d ", GetName(), iResult);
    TObjLink* lnk=tgtList.FirstLink();
    int i=iResult;
    char* pSpace = new char[iResult+1];
    if (pSpace) {
      memset(pSpace, 32, iResult);
      pSpace[i]=0;
      while (lnk) {
	TObject* obj=lnk->GetObject();
	HLTMessage("     %s^-- %s ", &pSpace[i--], obj->GetName());
	lnk=lnk->Next();
      }
      delete [] pSpace;
    } else {
      iResult=-ENOMEM;
    }
  }
}

int AliHLTTask::SetDependency(AliHLTTask* pDep)
{
  // see header file for function documentation
  int iResult=0;
  if (pDep) {
    if (FindDependency(pDep->GetName())==NULL) {
      fListDependencies.Add(pDep);
    } else {
      iResult=-EEXIST;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::UnsetDependency(AliHLTTask* pDep)
{
  // see header file for function documentation
  fListDependencies.Remove(pDep);
  if (fpConfiguration) {
    fpConfiguration->InvalidateSources();
  }
  return 0;
}

int AliHLTTask::CheckDependencies()
{
  // see header file for function documentation
  int iResult=0;
  AliHLTConfiguration* pSrc=fpConfiguration->GetFirstSource();
  while (pSrc) {
    if (FindDependency(pSrc->GetName())==NULL) {
      //HLTDebug("dependency \"%s\" unresolved", pSrc->GetName());
      iResult++;
    }
    pSrc=fpConfiguration->GetNextSource();
  }
  return iResult;
}


int AliHLTTask::Depends(AliHLTTask* pTask)
{
  // see header file for function documentation
  int iResult=0;
  if (pTask) {
    if (fpConfiguration) {
      iResult=fpConfiguration->GetSource(pTask->GetName())!=NULL;
      if (iResult>0) {
	//HLTDebug("task  \"%s\" depends on \"%s\"", GetName(), pTask->GetName());
      } else {
	//HLTDebug("task  \"%s\" independend of \"%s\"", GetName(), pTask->GetName());
      }
    } else {
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

AliHLTTask* AliHLTTask::FindTarget(const char* id)
{
  // see header file for function documentation
  AliHLTTask* pTask=NULL;
  if (id) {
    pTask=(AliHLTTask*)fListTargets.FindObject(id);
  }
  return pTask;
}

int AliHLTTask::SetTarget(AliHLTTask* pTgt)
{
  // see header file for function documentation
  int iResult=0;
  if (pTgt) {
    if (FindTarget(pTgt->GetName())==NULL) {
      fListTargets.Add(pTgt);
    } else {
      iResult=-EEXIST;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::UnsetTarget(AliHLTTask* pTarget)
{
  // see header file for function documentation
  fListTargets.Remove(pTarget);
  return 0;
}

int AliHLTTask::StartRun()
{
  // see header file for function documentation
  int iResult=0;
  int iNofInputDataBlocks=0;
  AliHLTComponent* pComponent=GetComponent();
  if (pComponent) {
    // determine the number of input data blocks provided from the source tasks
    TObjLink* lnk=fListDependencies.FirstLink();
    while (lnk && iResult>=0) {
      AliHLTTask* pSrcTask=(AliHLTTask*)lnk->GetObject();
      if (pSrcTask) {
	if ((iResult=pSrcTask->GetNofMatchingDataTypes(this))>0) {
	  iNofInputDataBlocks+=iResult;
	} else if (iResult==0) {
	  HLTWarning("source task %s (%p) does not provide any matching data type for task %s (%p)", pSrcTask->GetName(), pSrcTask, GetName(), this);
	} else {
	  HLTError("task %s (%p): error getting matching data types for source task %s (%p)", GetName(), this, pSrcTask->GetName(), pSrcTask);
	  iResult=-EFAULT;
	}
      }
      lnk=lnk->Next();
    }
    if (iResult>=0) {
      if (fBlockDataArray.size()>0) {
	HLTWarning("block data array for task %s (%p) was not cleaned", GetName(), this);
	fBlockDataArray.resize(0);
      }

      // component init
      // the initialization of the component is done by the ComponentHandler after creation
      // of the component.
      //iResult=Init( AliHLTComponentEnvironment* environ, void* environ_param, int argc, const char** argv );

      // allocate internal task variables for bookkeeping aso.
      // we allocate the BlockData array with at least one member
      if (iNofInputDataBlocks==0) iNofInputDataBlocks=1;
      AliHLTComponentBlockData init;
      memset(&init, 0, sizeof(AliHLTComponentBlockData));
      fBlockDataArray.resize(iNofInputDataBlocks, init);

      // allocate the data buffer, which controls the output buffer and subscriptions
      if (iResult>=0) {
	fpDataBuffer=new AliHLTDataBuffer;
	if (fpDataBuffer!=NULL) {
	  fpDataBuffer->SetLocalLoggingLevel(GetLocalLoggingLevel());
	  HLTDebug("created data buffer %p for task %s (%p)", fpDataBuffer, GetName(), this);
	  TObjLink* lnk=fListTargets.FirstLink();
	  while (lnk && iResult>=0) {
	    AliHLTTask* pTgtTask=(AliHLTTask*)lnk->GetObject();
	    if (pTgtTask) {
	      if ((iResult=fpDataBuffer->SetConsumer(pTgtTask->GetComponent()))>=0) {
	      }
	    } else {
	      break;
	      iResult=-EFAULT;
	    }
	    lnk=lnk->Next();
	  }
	} else {
	  HLTFatal("can not create data buffer object, memory allocation failed");
	  iResult=-ENOMEM;
	}
      }
    }
  } else {
    HLTError("task %s (%p) does not have a component", GetName(), this);
    iResult=-EFAULT;
  }
  return iResult;
}

int AliHLTTask::EndRun()
{
  // see header file for function documentation
  int iResult=0;
  if (fBlockDataArray.size()>0) {
    fBlockDataArray.resize(0);
  } else {
    HLTWarning("task %s (%p) doesn't seem to be in running mode", GetName(), this);
  }
  if (fpDataBuffer) {
    AliHLTDataBuffer* pBuffer=fpDataBuffer;
    fpDataBuffer=NULL;
    delete pBuffer;
  }
  return iResult;
}

int AliHLTTask::ProcessTask(Int_t eventNo)
{
  // see header file for function documentation
  int iResult=0;
  AliHLTComponent* pComponent=GetComponent();
  if (pComponent && fpDataBuffer) {
    HLTDebug("Processing task %s (%p) fpDataBuffer %p", GetName(), this, fpDataBuffer);
    fpDataBuffer->Reset();
    int iSourceDataBlock=0;
    int iInputDataVolume=0;

    AliHLTTask* pSrcTask=NULL;
    AliHLTTaskPList subscribedTaskList;
    TObjLink* lnk=fListDependencies.FirstLink();

    // subscribe to all source tasks
    while (lnk && iResult>=0) {
      pSrcTask=(AliHLTTask*)lnk->GetObject();
      if (pSrcTask) {
	int iMatchingDB=pSrcTask->GetNofMatchingDataBlocks(this);
	if (iMatchingDB>=0 && static_cast<unsigned int>(iMatchingDB)>fBlockDataArray.size()-iSourceDataBlock) {
	  AliHLTComponentBlockData init;
	  memset(&init, 0, sizeof(AliHLTComponentBlockData));
	  fBlockDataArray.resize(iSourceDataBlock+iMatchingDB, init);
	} else {
	  if (iMatchingDB<0) {
	    HLTError("task %s (%p): error getting no of matching data blocks from task %s (%p), error %d", GetName(), this, pSrcTask->GetName(), pSrcTask, iMatchingDB);
	    iResult=iMatchingDB;
	    break;
	  } else if (iMatchingDB==0) {
	    HLTDebug("source task %s (%p) does not provide any matching data type for task %s (%p)", pSrcTask->GetName(), pSrcTask, GetName(), this);
	  }
	}
	if ((iResult=pSrcTask->Subscribe(this, &fBlockDataArray[iSourceDataBlock],fBlockDataArray.size()-iSourceDataBlock))>=0) {
	  for (int i=0; i<iResult; i++) {
	    iInputDataVolume+=fBlockDataArray[i+iSourceDataBlock].fSize;
	    // put the source task as many times into the list as it provides data blocks
	    // makes the bookkeeping for the data release easier
	    subscribedTaskList.push_back(pSrcTask);
	  }
	  HLTDebug("Task %s (%p) successfully subscribed to %d data block(s) of task %s (%p)", GetName(), this, iResult, pSrcTask->GetName(), pSrcTask);
	  iSourceDataBlock+=iResult;	  
	  iResult=0;
	} else {
	  HLTError("Task %s (%p): subscription to task %s (%p) failed with error %d", GetName(), this, pSrcTask->GetName(), pSrcTask, iResult);
	  iResult=-EFAULT;
	}
      } else {
	HLTFatal("fatal internal error in ROOT list handling");
	iResult=-EFAULT;
      }
      lnk=lnk->Next();
    }

    // process the event
    int iNofTrial=0; // repeat processing if component returns -ENOSPC
    AliHLTUInt32_t size=0;
    if (iResult>=0) {
    do {
      long unsigned int iConstBase=0;
      double fInputMultiplier=0;
      if (pComponent->GetComponentType()!=AliHLTComponent::kSink)
	pComponent->GetOutputDataSize(iConstBase, fInputMultiplier);
      if (fInputMultiplier<0) {
	HLTWarning("ignoring negative input multiplier");
	fInputMultiplier=0;
      }
      long unsigned int iOutputDataSize=int(fInputMultiplier*iInputDataVolume) + iConstBase;
      //HLTDebug("task %s: reqired output size %d", GetName(), iOutputDataSize);
      if (iNofTrial>0) {
	// dont process again if the buffer size is the same
	if (size==iOutputDataSize) break;
	HLTInfo("processing task %s again with buffer size %d", GetName(), iOutputDataSize);
      }
      AliHLTUInt8_t* pTgtBuffer=NULL;
      if (iOutputDataSize>0) pTgtBuffer=fpDataBuffer->GetTargetBuffer(iOutputDataSize);
      //HLTDebug("provided raw buffer %p", pTgtBuffer);
      AliHLTComponentEventData evtData;
      AliHLTComponent::FillEventData(evtData);
      evtData.fEventID=(AliHLTEventID_t)eventNo;
      evtData.fBlockCnt=iSourceDataBlock;
      AliHLTComponentTriggerData trigData;
      size=iOutputDataSize;
      AliHLTUInt32_t outputBlockCnt=0;
      AliHLTComponentBlockData* outputBlocks=NULL;
      AliHLTComponentEventDoneData* edd;
      if (pTgtBuffer!=NULL || iOutputDataSize==0) {
	iResult=pComponent->ProcessEvent(evtData, &fBlockDataArray[0], trigData, pTgtBuffer, size, outputBlockCnt, outputBlocks, edd);
	HLTDebug("task %s: component %s ProcessEvent finnished (%d): size=%d blocks=%d", GetName(), pComponent->GetComponentID(), iResult, size, outputBlockCnt);
	if (iResult>=0 && outputBlocks) {
	  AliHLTComponentBlockDataList segments;
	  for (int oblock=0; oblock<outputBlockCnt; oblock++) {
	    int iblock=0;
	    for (; iblock<evtData.fBlockCnt; iblock++) {
	      if (fBlockDataArray[iblock].fPtr==outputBlocks[oblock].fPtr) {
		assert(subscribedTaskList[iblock]!=NULL);
		if (subscribedTaskList[iblock]==NULL) continue;
		HLTDebug("forward segment %d (source task %s %p) to data buffer %p", iblock, pSrcTask->GetName(), pSrcTask, fpDataBuffer);
		fpDataBuffer->Forward(subscribedTaskList[iblock], &fBlockDataArray[iblock]);
		subscribedTaskList[iblock]=NULL; // not to be released in the loop further down
		break;
	      }
	    }
	    if (iblock==evtData.fBlockCnt) segments.push_back(outputBlocks[oblock]);
	  }
	  if (pTgtBuffer && segments.size()>0) {
	    iResult=fpDataBuffer->SetSegments(pTgtBuffer, &segments[0], segments.size());
	  }
	  delete [] outputBlocks; outputBlocks=NULL; outputBlockCnt=0;
	} else {
	  fpDataBuffer->Reset();
	}
      } else {
	HLTError("task %s: no target buffer available", GetName());
	iResult=-EFAULT;
      }
    } while (iResult==-ENOSPC && iNofTrial++<1);
    }

    // now release all buffers which we have subscribed to
    iSourceDataBlock=0;
    AliHLTTaskPList::iterator element;
    while ((element=subscribedTaskList.begin())!=subscribedTaskList.end()) {
      pSrcTask=*element;
      if (pSrcTask) {
	int iTempRes=0;
	if ((iTempRes=pSrcTask->Release(&fBlockDataArray[iSourceDataBlock], this))>=0) {
	  HLTDebug("Task %s (%p) successfully released segment of task %s (%p)", GetName(), this, pSrcTask->GetName(), pSrcTask);
	} else {
	  HLTError("Task %s (%p): realease of task %s (%p) failed with error %d", GetName(), this, pSrcTask->GetName(), pSrcTask, iTempRes);
	}
      }
      subscribedTaskList.erase(element);
      iSourceDataBlock++;
    }
    if (subscribedTaskList.size()>0) {
      HLTError("task %s (%p): could not release all data buffers", GetName(), this);
    }
  } else {
    HLTError("task %s (%p): internal failure (not initialized component %p, data buffer %p)", GetName(), this, fpComponent, fpDataBuffer);
    iResult=-EFAULT;
  }
  return iResult;
}

int AliHLTTask::GetNofMatchingDataBlocks(const AliHLTTask* pConsumerTask) const
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumerTask) {
    if (fpDataBuffer) {
      iResult=fpDataBuffer->FindMatchingDataBlocks(pConsumerTask->GetComponent(), NULL);
    } else {
      HLTFatal("internal data buffer missing");
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::GetNofMatchingDataTypes(const AliHLTTask* pConsumerTask) const
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumerTask) {
    AliHLTComponent* pComponent=GetComponent();
    if (!pComponent) {
      // init ?
      HLTError("component not initialized");
      iResult=-EFAULT;
    }
    if (pComponent) {
      iResult=pComponent->FindMatchingDataTypes(pConsumerTask->GetComponent(), NULL);
    } else {
      HLTFatal("task initialization failed");
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::Subscribe(const AliHLTTask* pConsumerTask, AliHLTComponentBlockData* pBlockDesc, int iArraySize)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumerTask) {
    if (fpDataBuffer) {
      iResult=fpDataBuffer->Subscribe(pConsumerTask->GetComponent(), pBlockDesc, iArraySize);
    } else {
      HLTFatal("internal data buffer missing");
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTask::Release(AliHLTComponentBlockData* pBlockDesc, const AliHLTTask* pConsumerTask)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumerTask && pBlockDesc) {
    if (fpDataBuffer) {
      iResult=fpDataBuffer->Release(pBlockDesc, pConsumerTask->GetComponent(), this);
    } else {
      HLTFatal("internal data buffer missing");
      iResult=-EFAULT;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTTask::PrintStatus()
{
  // see header file for function documentation
  HLTLogKeyword("task properties");
  AliHLTComponent* pComponent=GetComponent();
  if (pComponent) {
    HLTMessage("     component: %s (%p)", pComponent->GetComponentID(), pComponent);
  } else {
    HLTMessage("     no component set!");
  }
  if (fpConfiguration) {
    AliHLTConfiguration* pSrc=fpConfiguration->GetFirstSource();
    while (pSrc) {
      const char* pQualifier="unresolved";
      if (FindDependency(pSrc->GetName()))
	pQualifier="resolved";
      HLTMessage("     source: %s (%s)", pSrc->GetName(), pQualifier);
      pSrc=fpConfiguration->GetNextSource();
    }
    TObjLink* lnk = fListTargets.FirstLink();
    while (lnk) {
      TObject *obj = lnk->GetObject();
      HLTMessage("     target: %s", obj->GetName());
      lnk = lnk->Next();
    }
  } else {
    HLTMessage("     task \"%s\" not initialized", GetName());
  }
}
