// $Id$
//**************************************************************************
//* This file is property of and copyright by the                          * 
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

/// @file   AliHLTTask.cxx
/// @author Matthias Richter
/// @date   
/// @brief  Implementation of HLT tasks.
///

#include <cerrno>
#include <cassert>
#include <iostream>
#include <string>
#include <ctime>
#include "AliHLTTask.h"
#include "AliHLTConfiguration.h"
#include "AliHLTConfigurationHandler.h"
#include "AliHLTComponent.h"
#include "AliHLTComponentHandler.h"
//#include "AliGRPManager.h"
//#include "AliGRPObject.h"
#include "TList.h"
#include "AliHLTErrorGuard.h"

using std::cout;

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
    HLTWarning("overriding existing reference to configuration object %p by %p",
	       fpConfiguration, pConf);
  }
  if (pConf!=NULL) fpConfiguration=pConf;
  iResult=CreateComponent(fpConfiguration, pCH, fpComponent);
  if (iResult>=0) {
    iResult=CustomInit(pCH);
  }
  return iResult;
}

int AliHLTTask::CreateComponent(AliHLTConfiguration* pConfiguration, AliHLTComponentHandler* pCH, AliHLTComponent*& pComponent) const
{
  // see header file for class documentation
  int iResult=0;
  if (!pConfiguration) return -EINVAL;

  const AliHLTConfiguration* pConf=AliHLTConfigurationHandler::FindSubstitution(*pConfiguration);
  if (!pConf) pConf=pConfiguration;
  if (pConf) {
    if (pCH) {
      int argc=0;
      const char** argv=NULL;
      if ((iResult=pConf->GetArguments(&argv))>=0) {
	argc=iResult; // just to make it clear
	// TODO: we have to think about the optional environment parameter,
	// currently just set to NULL.
	iResult=pCH->CreateComponent(pConf->GetComponentID(), pComponent);
	if (pComponent && iResult>=0) {
	  static bool grpInitialized = false;
	  //static const AliGRPObject *grpObj = NULL;
	  if (grpInitialized == false)
	  {
	    /*static AliGRPManager grp;
	    grp.ReadGRPEntry();
	    const AliGRPObject *grpObjTmp = grp.GetGRPData();
	    grpObj = grpObjTmp;*/
	    grpInitialized = true;
	  }
	  //if (grpObj) pComponent->SetTimeStamp(grpObj->GetTimeStart());
	
	  TString description;
	  description.Form("chainid=%s", GetName());
	  pComponent->SetComponentDescription(description.Data());
	  const AliHLTAnalysisEnvironment* pEnv=pCH->GetEnvironment();
	  if ((iResult=pComponent->Init(pEnv, NULL, argc, argv))>=0) {
	    //HLTDebug("component %s (%p) created", pComponent->GetComponentID(), pComponent); 
	  } else {
	    HLTError("Initialization of component \"%s\" failed with error %d", pComponent->GetComponentID(), iResult);
	  }
	} else {
	  //HLTError("can not find component \"%s\" (%d)", pConf->GetComponentID(), iResult);
	}
      } else {
	HLTError("can not get argument list for configuration %s (%s)", pConf->GetName(), pConf->GetComponentID());
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
  CustomCleanup();
  AliHLTComponent* pComponent=GetComponent();
  fpComponent=NULL;
  if (pComponent) {
    //HLTDebug("delete component %s (%p)", pComponent->GetComponentID(), pComponent); 
    static bool grpInitialized = false;
    //static const AliGRPObject *grpObj = NULL;
    if (grpInitialized == false)
    {
      /*static AliGRPManager grp;
      grp.ReadGRPEntry();
      const AliGRPObject *grpObjTmp = grp.GetGRPData();
      grpObj = grpObjTmp;*/
      grpInitialized = true;
    }
    //if (grpObj) pComponent->SetTimeStamp(grpObj->GetTimeEnd());

    pComponent->Deinit();
    delete pComponent;
  } else {
    HLTWarning("task doesn't seem to be in initialized");
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
    HLTMessage("     dependency level %d ", iResult);
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
    { // set scope for lnk as a local variable
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
    }
    if (iResult>=0) {
      if (fBlockDataArray.size()>0) {
	HLTWarning("block data array for task %s (%p) was not cleaned", GetName(), this);
	fBlockDataArray.clear();
      }

      // component init
      // the initialization of the component is done by the ComponentHandler after creation
      // of the component.
      //iResult=Init( AliHLTAnalysisEnvironment* environ, void* environ_param, int argc, const char** argv );

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
	      iResult=-EFAULT;
	      break;
	    }
	    lnk=lnk->Next();
	  }
	} else {
	  HLTFatal("can not create data buffer object, memory allocation failed");
	  iResult=-ENOMEM;
	}
      }
    }
    if (iResult>=0) {
      // send the SOR event
      
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
    fBlockDataArray.clear();
  }
  if (fpDataBuffer) {
    AliHLTDataBuffer* pBuffer=fpDataBuffer;
    fpDataBuffer=NULL;
    delete pBuffer;
  }
  return iResult;
}

int AliHLTTask::ProcessTask(Int_t eventNo, AliHLTUInt32_t eventType, AliHLTTriggerMask_t trgMask,
			    AliHLTUInt32_t timestamp, AliHLTUInt32_t participatingDetectors)
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

    // instances of SOR and EOR events to be kept
    int iSOR=-1;
    int iEOR=-1;
    // TODO 2009-09-30
    // generalize handling of the special blocks to be forwarded on SOR and EOR
    // just adding a new specific handling for the ECS parameter block as a quick
    // solution
    int iECS=-1;

    // subscribe to all source tasks
    fBlockDataArray.clear();
    while (lnk && iResult>=0) {
      pSrcTask=(AliHLTTask*)lnk->GetObject();
      if (pSrcTask) {
	int iMatchingDB=pSrcTask->GetNofMatchingDataBlocks(this);
	if (iMatchingDB<0) {
	  HLTError("task %s (%p): error getting no of matching data blocks from task %s (%p), error %d", GetName(), this, pSrcTask->GetName(), pSrcTask, iMatchingDB);
	  iResult=iMatchingDB;
	  break;
	} else if (iMatchingDB==0) {
	  HLTDebug("source task %s (%p) does not provide any matching data type for task %s (%p)", pSrcTask->GetName(), pSrcTask, GetName(), this);
	}
	if ((iResult=pSrcTask->Subscribe(this, fBlockDataArray))>=0) {
	  iSOR=iEOR=iECS=-1;
	  AliHLTComponentBlockDataList::iterator block=fBlockDataArray.begin();
	  for (int i=0; block!=fBlockDataArray.end(); i++) {
	    bool bRemove=0;
	    bRemove|=(*block).fDataType==kAliHLTDataTypeSOR && !(iSOR<0 && (iSOR=i)>=0);
	    bRemove|=(*block).fDataType==kAliHLTDataTypeEOR && !(iEOR<0 && (iEOR=i)>=0);
	    bRemove|=(*block).fDataType==kAliHLTDataTypeECSParam && !(iECS<0 && (iECS=i)>=0);
	    //HLTInfo("block %d, iSOR=%d iEOR=%d remove=%d", i, iSOR, iEOR, bRemove);
	    if (i<iSourceDataBlock) {
	      assert(!bRemove);
	    } else if (bRemove) {
	      HLTDebug("remove duplicated event %s (%d)", AliHLTComponent::DataType2Text((*block).fDataType).c_str(), i);
	      pSrcTask->Release(&(*block), this);
	      block=fBlockDataArray.erase(block);
	      continue;
	    } else {
	    iInputDataVolume+=(*block).fSize;
	    // put the source task as many times into the list as it provides data blocks
	    // makes the bookkeeping for the data release easier
	    subscribedTaskList.push_back(pSrcTask);
	    }
	    block++;
	  }
	  HLTDebug("Task %s (%p) successfully subscribed to %d data block(s) of task %s (%p)", GetName(), this, iResult, pSrcTask->GetName(), pSrcTask);
	  iSourceDataBlock=fBlockDataArray.size();	  
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
    AliHLTUInt32_t iLastOutputDataSize=0;
    if (iResult>=0) {
    do {
      long unsigned int iOutputDataSize=0;
      AliHLTConfiguration* pConf=GetConf();
      // check if there was a buffer size specified, query output size
      // estimator from component otherwize
      if (pConf && pConf->GetOutputBufferSize()>=0) {
	iOutputDataSize=pConf->GetOutputBufferSize();
      } else {
      long unsigned int iConstBase=0;
      double fInputMultiplier=0;
      if (pComponent->GetComponentType()!=AliHLTComponent::kSink) {
	pComponent->GetOutputDataSize(iConstBase, fInputMultiplier);
	// add a small margin to the buffer to allow optional component
	// statistics
	iConstBase+=100;
#if defined(__DEBUG) || defined(HLT_COMPONENT_STATISTICS)
	for (AliHLTComponentBlockDataList::iterator element=fBlockDataArray.begin();
	     element!=fBlockDataArray.end(); element++) {
	  if (element->fDataType==kAliHLTDataTypeComponentStatistics) {
	    iConstBase+=element->fSize;
	  }
	}
#endif
      }
      if (fInputMultiplier<0) {
	HLTWarning("ignoring negative input multiplier");
	fInputMultiplier=0;
      }
      iOutputDataSize=(long unsigned int)(fInputMultiplier*iInputDataVolume) + iConstBase;
      //HLTDebug("task %s: reqired output size %d", GetName(), iOutputDataSize);
      }
      if (iNofTrial == 0 && fpDataBuffer->GetMaxBufferSize() < iOutputDataSize) {
        //If the estimated buffer size exceeds the maximum buffer size of AliHLTRawBuffer, decrease the buffer size.
        //The estimation is often quite high, and GetMaxBufferSize should usually return a size that is sufficient.
        HLTImportant("Reducing estimated output buffer size of %lu to maximum output buffer size\n", iOutputDataSize);
        iOutputDataSize = fpDataBuffer->GetMaxBufferSize();
      }
      if (iNofTrial>0) {
	// dont process again if the buffer size is the same
	if (iLastOutputDataSize>=iOutputDataSize) break;
	HLTImportant("processing event %d again with buffer size %d", eventNo, iOutputDataSize);
      }
      AliHLTUInt8_t* pTgtBuffer=NULL;
      if (iOutputDataSize>0) pTgtBuffer=fpDataBuffer->GetTargetBuffer(iOutputDataSize);
      //HLTDebug("provided raw buffer %p", pTgtBuffer);
      AliHLTComponentEventData evtData;
      AliHLTComponent::FillEventData(evtData);
      if (eventNo>=0)
	evtData.fEventID=(AliHLTEventID_t)eventNo;
      if (timestamp < kMaxUInt) evtData.fEventCreation_s=timestamp;
      else
      evtData.fEventCreation_s=static_cast<AliHLTUInt32_t>(time(NULL));
      AliHLTComponentTriggerData trigData;
      AliHLTEventTriggerData evtTrigData;
      trigData.fStructSize=sizeof(trigData);
      trigData.fDataSize=sizeof(AliHLTEventTriggerData);
      memset(&evtTrigData, 0, trigData.fDataSize);
      // Setup the CDH in the trigger data, based on the event type, CTP trigger
      // mask and participating detectors.
      evtTrigData.fCommonHeaderWordCnt=gkAliHLTCommonHeaderCount;
      AliHLTUInt8_t l1msg = 0x0;
      switch (eventType)
      {
      case gkAliEventTypeData:        l1msg = 0x00; break;
      case gkAliEventTypeDataReplay:  l1msg = 0x00; break;
      case gkAliEventTypeStartOfRun:  l1msg = (0xE << 2) | 0x01; break;
      case gkAliEventTypeEndOfRun:    l1msg = (0xF << 2) | 0x01; break;
      case gkAliEventTypeCalibration: l1msg = (0x1 << 6) | 0x01; break;
      case gkAliEventTypeSoftware:    l1msg = 0x01; break;
      }
      evtTrigData.fCommonHeader[1] = (AliHLTUInt32_t(l1msg) << 14) | (0x03000000); //We have to set a CDH version, use 3 because we support 100 trigger classes below (0x03000000)
      evtTrigData.fCommonHeader[3] = ((l1msg & 0x1) == 0x1) ? (participatingDetectors & 0xFFFFFF) : 0x0;
      evtTrigData.fCommonHeader[5] = (trgMask & AliHLTTriggerMask_t(0xffffffff)).to_ulong();
      evtTrigData.fCommonHeader[6] = ((trgMask>>32) & AliHLTTriggerMask_t(0x3ffff)).to_ulong();
      evtTrigData.fCommonHeader[7] = ((trgMask>>50) & AliHLTTriggerMask_t(0xffffffff)).to_ulong();
      evtTrigData.fCommonHeader[8] = ((trgMask>>72) & AliHLTTriggerMask_t(0x3ffff)).to_ulong();

      trigData.fData=&evtTrigData;
      iLastOutputDataSize=iOutputDataSize;
      AliHLTUInt32_t size=iOutputDataSize;
      AliHLTUInt32_t outputBlockCnt=0;
      AliHLTComponentBlockData* outputBlocks=NULL;
      AliHLTComponentEventDoneData* edd=NULL;
      if (pTgtBuffer!=NULL || iOutputDataSize==0) {
	// add event type data block
	// the block is removed immediately after processing from the list
	AliHLTComponentBlockData eventTypeBlock;
	AliHLTComponent::FillBlockData(eventTypeBlock);
	// Note: no payload!
	eventTypeBlock.fDataType=kAliHLTDataTypeEvent;
	eventTypeBlock.fSpecification=eventType;
	fBlockDataArray.push_back(eventTypeBlock);

	if (CheckFilter(kHLTLogDebug)) Print("proc");

	AliHLTUInt32_t iblock=0;
	// check input and output buffers for consistency
	// to be enabled after fixing bug with DataBuffer and forwarded SOR/EOR
	//for (iblock=0; iblock<fBlockDataArray.size(); iblock++) {
	//  if ((AliHLTUInt8_t*)fBlockDataArray[iblock].fPtr >= pTgtBuffer+size) continue;
	//  if (pTgtBuffer >= (AliHLTUInt8_t*)fBlockDataArray[iblock].fPtr+fBlockDataArray[iblock].fSize) continue;
	//  HLTFatal("input and output buffer overlap for block descriptor %d (ptr %p size %d): output buffer %p %d",
	//	   iblock, fBlockDataArray[iblock].fPtr, fBlockDataArray[iblock].fSize,
	//	   pTgtBuffer, size);
	//}	

	// process
	evtData.fBlockCnt=fBlockDataArray.size();
	iResult=pComponent->ProcessEvent(evtData, &fBlockDataArray[0], trigData, pTgtBuffer, size, outputBlockCnt, outputBlocks, edd);
	HLTDebug("component %s ProcessEvent finnished (%d): size=%d blocks=%d", pComponent->GetComponentID(), iResult, size, outputBlockCnt);

	// EventDoneData is for the moment ignored in AliHLTSystem
	if (edd) {
	  HLTDebug("got EventDoneData size %d", edd->fDataSize);
	  delete [] reinterpret_cast<char*>(edd);
	  edd=NULL;
	}

	// remove event data block
	fBlockDataArray.pop_back();

	// check for forwarded blocks.
	// loop over all output blocks and check
	// 1. for duplicate blocks (pointing to same location in output buffer
	//    or to the same buffer)
	// 2. for blocks forwarded from the input.
	if (iResult>=0 && outputBlocks) {
	  if (fListTargets.First()!=NULL) {
	    AliHLTComponentBlockDataList segments;
	    for (AliHLTUInt32_t oblock=0; oblock<outputBlockCnt; oblock++) {
	      // consistency check for data reference
	      if (outputBlocks[oblock].fPtr!=NULL && outputBlocks[oblock].fPtr!=pTgtBuffer &&
		  outputBlocks[oblock].fOffset!=0) {
		HLTWarning("output block %s 0x%08x has inconsistent data reference ptr=%p offset=0x%08x: "
			   "for new blocks use offset only, forwarded blocks have fPtr set only",
			   AliHLTComponent::DataType2Text(outputBlocks[oblock].fDataType).c_str(),
			   outputBlocks[oblock].fSpecification,
			   outputBlocks[oblock].fPtr, outputBlocks[oblock].fOffset);
	      }

	      // check for duplicates in the output
	      // this check applies for forwarded data blocks where
	      // the ptr is not inside the target buffer
	      AliHLTUInt32_t checkblock=0;
	      for (; checkblock<oblock; checkblock++) {
		if (outputBlocks[oblock].fPtr!=NULL && outputBlocks[oblock].fPtr!=pTgtBuffer &&
		    outputBlocks[checkblock].fPtr==outputBlocks[oblock].fPtr) {
		  if (outputBlocks[checkblock].fSize!=outputBlocks[oblock].fSize ||
		      outputBlocks[checkblock].fDataType!=outputBlocks[oblock].fDataType) {
		    HLTWarning("output blocks %d (%s 0x%08x) and %d (%s 0x%08x) have identical data references ptr=%p "
			       "but differ in data type and/or size: %d vs. %d, ignoring block %d",
			       oblock,
			       AliHLTComponent::DataType2Text(outputBlocks[oblock].fDataType).c_str(),
			       outputBlocks[oblock].fSpecification,
			       checkblock,
			       AliHLTComponent::DataType2Text(outputBlocks[checkblock].fDataType).c_str(),
			       outputBlocks[checkblock].fSpecification,
			       outputBlocks[oblock].fPtr,
			       outputBlocks[oblock].fSize,
			       outputBlocks[checkblock].fSize,
			       checkblock);
		  }
		  // ignore from the second copy
		  break;
		}
	      }
	      if (checkblock<oblock) continue;

	      // search for the forwarded data blocks
	      // new data blocks are announced to the data buffer, forwarded data blocks
	      // to the publisher task. The publisher task of a forwarded data block is
	      // removed from the list in order to keep the buffer open. It will be releases
	      // when the subscribing task releases it
	      iblock=0;
	      for (; iblock<fBlockDataArray.size(); iblock++) {
		if (outputBlocks[oblock].fDataType==kAliHLTDataTypeEvent) {
		  // the event type data block is an artificial data block
		  // ignore if it was forwarded
		  break;
		}
		if (fBlockDataArray[iblock].fPtr==outputBlocks[oblock].fPtr) {
		  assert(subscribedTaskList[iblock]!=NULL);
		  if (subscribedTaskList[iblock]==NULL) {
		    ALIHLTERRORGUARD(1, "missing parent task for forwarded data block %s 0x%08x, original data block %s 0x%08x, subsequent errors are suppressed",
				     AliHLTComponent::DataType2Text(outputBlocks[oblock].fDataType).c_str(),
				     outputBlocks[oblock].fSpecification,
				     AliHLTComponent::DataType2Text(outputBlocks[iblock].fDataType).c_str(),
				     outputBlocks[iblock].fSpecification);
		    continue;
		  }
		  HLTDebug("forward segment %d (source task %s %p) to data buffer %p", iblock, pSrcTask->GetName(), pSrcTask, fpDataBuffer);
		  fpDataBuffer->Forward(subscribedTaskList[iblock], &fBlockDataArray[iblock]);
		  subscribedTaskList[iblock]=NULL; // not to be released in the loop further down
		  break;
		}
	      }
	      if (iblock==fBlockDataArray.size()) segments.push_back(outputBlocks[oblock]);
	    }
	    if (pTgtBuffer && segments.size()>0) {
	      iResult=fpDataBuffer->SetSegments(pTgtBuffer, &segments[0], segments.size());
	    }
	  } else {
	    // no forwarding, actually we dont even need to keep the data, this is a
	    // dead end (fListTargets empty)
	    //iResult=fpDataBuffer->SetSegments(pTgtBuffer, outputBlocks, outputBlockCnt);
	  }
	  delete [] outputBlocks; outputBlocks=NULL; outputBlockCnt=0;
	} else {
	  fpDataBuffer->Reset();
	}
	if (fListTargets.First()!=NULL) {
	  if (iSOR>=0 && subscribedTaskList[iSOR]!=NULL) {
	    HLTDebug("forward SOR event segment %d (source task %s %p) to data buffer %p", iSOR, pSrcTask->GetName(), pSrcTask, fpDataBuffer);
	    fpDataBuffer->Forward(subscribedTaskList[iSOR], &fBlockDataArray[iSOR]);
	    subscribedTaskList[iSOR]=NULL; // not to be released in the loop further down
	  }
	  if (iEOR>=0 && subscribedTaskList[iEOR]!=NULL) {
	    HLTDebug("forward EOR event (%s) segment %d (source task %s %p) to data buffer %p", AliHLTComponent::DataType2Text(fBlockDataArray[iEOR].fDataType).c_str(), iEOR, pSrcTask->GetName(), pSrcTask, fpDataBuffer);
	    fpDataBuffer->Forward(subscribedTaskList[iEOR], &fBlockDataArray[iEOR]);
	    subscribedTaskList[iEOR]=NULL; // not to be released in the loop further down
	  }
	  if (iECS>=0 && subscribedTaskList[iECS]!=NULL) {
	    HLTDebug("forward ECS event (%s) segment %d (source task %s %p) to data buffer %p", AliHLTComponent::DataType2Text(fBlockDataArray[iECS].fDataType).c_str(), iECS, pSrcTask->GetName(), pSrcTask, fpDataBuffer);
	    fpDataBuffer->Forward(subscribedTaskList[iECS], &fBlockDataArray[iECS]);
	    subscribedTaskList[iECS]=NULL; // not to be released in the loop further down
	  }
	}
      } else {
	HLTError("no target buffer available");
	iResult=-EFAULT;
      }
    } while (iResult==-ENOSPC && iNofTrial++<1);
    }

    fBlockDataArray.clear();
    if (CheckFilter(kHLTLogDebug)) Print("proc");

    // now release all buffers which we have subscribed to
    iSourceDataBlock=0;
    AliHLTTaskPList::iterator element;
    while ((element=subscribedTaskList.begin())!=subscribedTaskList.end()) {
      pSrcTask=*element;
      if (pSrcTask) {
	int iTempRes=0;
	if ((iTempRes=pSrcTask->Release(&fBlockDataArray[iSourceDataBlock], this))>=0) {
	  HLTDebug("successfully released segment of task %s (%p)", pSrcTask->GetName(), pSrcTask);
	} else {
	  HLTError("realease of task %s (%p) failed with error %d", pSrcTask->GetName(), pSrcTask, iTempRes);
	}
      }
      subscribedTaskList.erase(element);
      iSourceDataBlock++;
    }
    if (subscribedTaskList.size()>0) {
      HLTError("could not release all data buffers");
    }
  } else {
    HLTError("internal failure (not initialized component %p, data buffer %p)", fpComponent, fpDataBuffer);
    iResult=-EFAULT;
  }
  return iResult;
}

int AliHLTTask::SubscribeSourcesAndSkip()
{
  // function carries out the proper cleanup of the source components
  // by subscribing and releasing
  int iResult=0;  
  AliHLTTask* pSrcTask=NULL;
  AliHLTTaskPList subscribedTaskList;

  // cleanup the data buffer
  if (fpDataBuffer) fpDataBuffer->Reset();

  // subscribe to all source tasks
  fBlockDataArray.clear();
  for (TObjLink* lnk=fListDependencies.FirstLink(); lnk!=NULL; lnk=lnk->Next()) {
    if (!lnk->GetObject()) continue;
    pSrcTask=dynamic_cast<AliHLTTask*>(lnk->GetObject());
    if (!pSrcTask) continue;
    unsigned iPosition=fBlockDataArray.size();
    if ((iResult=pSrcTask->Subscribe(this, fBlockDataArray))>0) {
      for (unsigned i=iPosition; i<fBlockDataArray.size(); i++) {
	subscribedTaskList.push_back(pSrcTask);
      }
      HLTDebug("subscribed to %d blocks of task %s (%p)", iResult, pSrcTask->GetName(), pSrcTask, iResult);
    } else if (iResult<0) {
      HLTError("failed to subscribe to task %s (%p) with error %d", pSrcTask->GetName(), pSrcTask, iResult);
    }
  }

  unsigned iSourceDataBlock=0;
  AliHLTTaskPList::iterator element;
  while ((element=subscribedTaskList.begin())!=subscribedTaskList.end()) {
    assert(iSourceDataBlock<fBlockDataArray.size());
    pSrcTask=*element;
    if (pSrcTask && iSourceDataBlock<fBlockDataArray.size()) {
      if ((iResult=pSrcTask->Release(&fBlockDataArray[iSourceDataBlock], this))>=0) {
	HLTDebug("successfully released segment of task %s (%p)", pSrcTask->GetName(), pSrcTask);
      } else if (iSourceDataBlock>=fBlockDataArray.size()) {
	HLTError("mismatch between list of subscribed tasks and block list in task %s (%p), can not release task %s (%p)", GetName(), this, pSrcTask->GetName(), pSrcTask);
      } else {
	HLTError("realease of task %s (%p) failed with error %d", pSrcTask->GetName(), pSrcTask, iResult);
      }
    }
    subscribedTaskList.erase(element);
    iSourceDataBlock++;
  }
  if (iSourceDataBlock<fBlockDataArray.size()) {
    HLTWarning("not all subscriptions released for task %s (%p)", GetName(), this);
  }

  return 0;
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

int AliHLTTask::Subscribe(const AliHLTTask* pConsumerTask, AliHLTComponentBlockDataList& blockDescList)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumerTask) {
    if (fpDataBuffer) {
      iResult=fpDataBuffer->Subscribe(pConsumerTask->GetComponent(), blockDescList);
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
    HLTMessage("     task not initialized");
  }
}

void AliHLTTask::Print(const char* options) const
{
  // Overloaded from TObject
  if (strcmp(options, "proc")==0) {
    // print processing info
    HLTMessage("**********************************************");
    HLTMessage("******* AliHLTTask Processing info ***********");
    HLTMessage(" component: %p %s", fpComponent, (fpComponent?fpComponent->GetComponentID():""));
    HLTMessage(" data buffer: %p", fpDataBuffer);
    if (fpDataBuffer) fpDataBuffer->Print("");
    HLTMessage(" input block descriptors: %d", fBlockDataArray.size());
    for (unsigned i=0; i<fBlockDataArray.size(); i++) {
      HLTMessage("  %d: %s 0x%08x %p %d", i,
		 AliHLTComponent::DataType2Text(fBlockDataArray[i].fDataType).c_str(),
		 fBlockDataArray[i].fSpecification,
		 fBlockDataArray[i].fPtr,
		 fBlockDataArray[i].fSize
		 );
    }
    HLTMessage("**** end of AliHLTTask Processing info *******");
    HLTMessage("**********************************************");
    return;
  }

  cout << "AliHLTTask " << GetName() << " " << this
       << " component " << fpComponent << " "
       << (fpComponent?fpComponent->GetComponentID():"")
       << endl;
}


int AliHLTTask::CustomInit(AliHLTComponentHandler* /*pCH*/)
{
  // default implementation nothing to do
  return 0;
}

int AliHLTTask::CustomCleanup()
{
  // default implementation nothing to do
  return 0;
}

int AliHLTTask::LoggingVarargs(AliHLTComponentLogSeverity severity, 
				    const char* originClass, const char* originFunc,
				    const char* file, int line, ... ) const
{
  // see header file for function documentation
  int iResult=0;

  va_list args;
  va_start(args, line);

  AliHLTLogging::SetLogString(this, " (%p)", "%s_pfmt_: ", GetName());
  iResult=SendMessage(severity, originClass, originFunc, file, line, AliHLTLogging::BuildLogString(NULL, args, true /*append*/));
  va_end(args);

  return iResult;
}
