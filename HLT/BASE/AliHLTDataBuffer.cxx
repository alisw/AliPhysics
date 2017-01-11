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

//  @file   AliHLTDataBuffer.cxx
//  @author Matthias Richter
//  @date   
//  @brief  Handling of Data Buffers for HLT components.
//  @note   Only used in the AliRoot framework

#include "AliHLTDataBuffer.h"
#include "AliHLTConsumerDescriptor.h"
#include "AliHLTComponent.h"
#include "AliHLTTask.h"
#include <cerrno>
#include <cassert>
//#include <string>
//#include "AliHLTSystem.h"

using std::cout;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataBuffer)

AliHLTDataBuffer::AliHLTDataBuffer()
  :
  fSegments(),
  fConsumers(),
  fActiveConsumers(),
  fReleasedConsumers(),
  fpBuffer(NULL),
  fFlags(0),
  fForwardedSegmentSources(),
  fForwardedSegments()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  fSegments.empty();
  fConsumers.empty();
  fActiveConsumers.empty();
  fReleasedConsumers.empty();
  fgNofInstances++;
}

int AliHLTDataBuffer::fgNofInstances=0;
AliHLTDataBuffer::AliHLTRawBufferPList AliHLTDataBuffer::fgFreeBuffers;
AliHLTDataBuffer::AliHLTRawBufferPList AliHLTDataBuffer::fgActiveBuffers;
AliHLTUInt32_t AliHLTDataBuffer::fgMargin=1024;
AliHLTLogging AliHLTDataBuffer::fgLogging;
const Int_t AliHLTDataBuffer::fgkSafetyPatternSize=16;
const char AliHLTDataBuffer::fgkSafetyPattern[]={0x28, 0x63, 0x29, 0x4d, 0x52, 0x49, 0x43, 0x48, 0x54, 0x45, 0x52, 0x20, 0x32, 0x30, 0x30, 0x37};
AliHLTUInt32_t AliHLTDataBuffer::fgEventCount=0;

AliHLTDataBuffer::~AliHLTDataBuffer()
{
  // see header file for function documentation
  CleanupConsumerList();

  if (--fgNofInstances<=0) {
    DeleteRawBuffers();
  }
}

int AliHLTDataBuffer::SetConsumer(AliHLTComponent* pConsumer)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumer) {
    if (FindConsumer(pConsumer)) {
      HLTWarning("consumer %s (%p) already set to data buffer %p", pConsumer->GetComponentID(), pConsumer, this);
    }
    AliHLTConsumerDescriptor* pDesc=new AliHLTConsumerDescriptor(pConsumer);
    if (pDesc) {
      fConsumers.push_back(pDesc);
      HLTDebug("set consumer %s (%p) to data buffer %p", pConsumer->GetComponentID(), pConsumer, this);
    } else {
      HLTError("memory allocation failed");
      iResult=-ENOMEM;
    }
  } else {
    HLTError("invalid parameter: consumer component (nil)");
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::FindMatchingDataBlocks(const AliHLTComponent* pConsumer, AliHLTComponentDataTypeList* tgtList)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumer) {
    AliHLTDataSegmentList segments;
    if ((iResult=FindMatchingDataSegments(pConsumer, segments))>=0) {
      if (tgtList) {
	AliHLTDataSegmentList::iterator segment=segments.begin();
	while (segment!=segments.end()) {
	  tgtList->push_back((*segment).fDataType);
	  segment++;
	}
      }
      iResult=segments.size();
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::FindMatchingDataSegments(const AliHLTComponent* /*pConsumer*/,
					       AliHLTDataSegmentList& tgtList)
{
  // see header file for function documentation
  int iResult=0;

  // Matthias 26.09.2007 relax the restriction to matching data blocks
  // all blocks are passed to the consumer, which is the policy also in
  // PubSub
  tgtList.assign(fSegments.begin(), fSegments.end());

  // add all forwarded blocks
  tgtList.insert(tgtList.begin(), fForwardedSegments.begin(), fForwardedSegments.end());
  iResult=tgtList.size();
  return iResult;
}

int AliHLTDataBuffer::Subscribe(const AliHLTComponent* pConsumer, AliHLTComponentBlockDataList& blockDescList)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumer) {
    if (1/*fpBuffer*/) {
      AliHLTConsumerDescriptor* pDesc=FindConsumer(pConsumer, fConsumers);
      if (pDesc) {
	AliHLTDataSegmentList tgtList;
	// Matthias 26.07.2007 AliHLTSystem should behave the same way as PubSub
	// so it does not matter if there are matching data types or not, unless
	// we implement such a check in PubSub
	if ((iResult=FindMatchingDataSegments(pConsumer, tgtList))>=0) {
	  AliHLTDataSegmentList::iterator segment=tgtList.begin();
	  while (segment!=tgtList.end()) {
	    // fill the block data descriptor
	    AliHLTComponentBlockData bd;
	    AliHLTComponent::FillBlockData(bd);
	    // This models the behavior of PubSub.
	    // For incoming data blocks, fOffset must be ignored by the
	    // processing component. It is set for bookkeeping in the framework.
	    // fPtr always points to the beginning of the data.
	    bd.fOffset=0;
	    AliHLTUInt8_t* pTgt=*segment;
	    bd.fPtr=reinterpret_cast<void*>(pTgt);
	    bd.fSize=(*segment).fSegmentSize;
	    bd.fDataType=(*segment).fDataType;
	    bd.fSpecification=(*segment).fSpecification;
	    blockDescList.push_back(bd);
	    pDesc->SetActiveDataSegment(*segment);
	    HLTDebug("component %p (%s) subscribed to segment offset %d size %d data type %s %#x", 
		     pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID(), bd.fOffset,
		     bd.fSize, (AliHLTComponent::DataType2Text(bd.fDataType)).c_str(), 
		     bd.fSpecification);
	    segment++;
	  }
	  // move this consumer to the active list
	  if (tgtList.size()==0) {
	    ChangeConsumerState(pDesc, fConsumers, fReleasedConsumers);
	    HLTDebug("no input data for component %p (%s) available", pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID());
	  } else if (ChangeConsumerState(pDesc, fConsumers, fActiveConsumers)>=0) {
	    HLTDebug("component %p (%s) subscribed to data buffer %p", pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID(), this);
	  } else {
	    // TODO: cleanup the consumer descriptor correctly
	    segment=tgtList.begin();
	    while (segment!=tgtList.end()) {
	      blockDescList.pop_back();
	      segment++;
	    }
	    HLTError("can not activate consumer %p for data buffer %p", pConsumer, this);
	    iResult=-EACCES;
	  }
	} else {
	  HLTError("unresolved data segment(s) for component %p (%s)", pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID());
	  iResult=-EBADF;
	}
      } else {
	if (!FindConsumer(pConsumer)) {
	  HLTError("component %p is not a data consumer of data buffer %p", pConsumer, this);
	} else {
	  HLTError("component %p is a valid data consumer of data buffer %p, but did not release it's buffer subscription", pConsumer, this);
	}
	iResult=-ENOENT;
      }
    } else {
      // Matthias 26.07.2007 until now, data had to be present for successful subscription
      // in order to be consistent with the PubSub framework, this restiction has been
      // removed
      //HLTError("data buffer %p is empty", this);
      //iResult=-ENODATA;
    }
  } else {
    HLTError("invalid parameter");
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::Release(AliHLTComponentBlockData* pBlockDesc,
			      const AliHLTComponent* pConsumer,
			      const AliHLTTask* pOwnerTask)
{
  // see header file for function documentation
  int iResult=0;
  if (pBlockDesc && pConsumer) {
    AliHLTConsumerDescriptor* pDesc=FindConsumer(pConsumer, fActiveConsumers);
    if (pDesc) {
      if ((iResult=pDesc->CheckActiveDataSegment(AliHLTDataSegment(pBlockDesc->fPtr, pBlockDesc->fOffset, pBlockDesc->fSize)))!=1) {
	HLTWarning("data segment mismatch, component %p has not subscribed to a segment with offset %#x and size %d", pConsumer, pBlockDesc->fOffset, pBlockDesc->fSize);
	// TODO: appropriate error handling, but so far optional
	iResult=0;
      } else {
	pDesc->ReleaseActiveDataSegment(AliHLTDataSegment(pBlockDesc->fPtr, pBlockDesc->fOffset, pBlockDesc->fSize));
      }
      if (GetNofPendingConsumers()==0 && fForwardedSegments.size()>0) {
	// last consumer, release forwarded segments
	ReleaseForwardedBlock(pBlockDesc, pOwnerTask);
      }
      pBlockDesc->fOffset=0;
      pBlockDesc->fPtr=NULL;
      pBlockDesc->fSize=0;
      if (pDesc->GetNofActiveSegments()==0) {
	if ((iResult=ChangeConsumerState(pDesc, fActiveConsumers, fReleasedConsumers))>=0) {
	  if (GetNofActiveConsumers()==0 && GetNofPendingConsumers()==0) {
	    // this is the last consumer, reset the consumer list and release the raw buffer
	    ResetDataBuffer();
	  }
	} else {
	  HLTError("can not deactivate consumer %p for data buffer %p", pConsumer, this);
	  iResult=-EACCES;
	}
      }
    } else {
      HLTWarning("component %p has currently not subscribed to the data buffer %p", pConsumer, this);
      iResult=-ENOENT;
    }
  } else {
    HLTError("inavalid parameter: pBlockDesc=%p pConsumer=%p", pBlockDesc, pConsumer);
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::ReleaseForwardedBlock(AliHLTComponentBlockData* pBlockDesc,
					    const AliHLTTask* pOwnerTask)
{
  // see header file for function documentation
  int iResult=0;
  if (pBlockDesc && pOwnerTask) {
	assert(fForwardedSegments.size()==fForwardedSegmentSources.size());
	AliHLTDataSegmentList::iterator segment=fForwardedSegments.begin();
	AliHLTTaskPList::iterator src=fForwardedSegmentSources.begin();
	//HLTDebug("%p checking forwarded segments", this);
	for (; segment!=fForwardedSegments.end(); segment++, src++) {
	  //HLTDebug("segment ptr=%p offset=%d size=%d\n"
	  //   "block ptr=%p offset=%d size=%d", (*segment).fPtr, (*segment).fSegmentOffset, (*segment).fSegmentSize, pBlockDesc->fPtr, pBlockDesc->fOffset, pBlockDesc->fSize);
	  if ((*segment)==AliHLTDataSegment(pBlockDesc->fPtr, pBlockDesc->fOffset, pBlockDesc->fSize)) {
	    //HLTDebug("release segment of task %p", *src);
	    assert((*src)!=NULL);
	    if ((*src)!=NULL) {
	      if ((*src)->Release(pBlockDesc, pOwnerTask)>=0) {
		HLTDebug("task %s (%p) released forwarded segment %p size %d of task %s (%p)",
			 pOwnerTask->GetName(), pOwnerTask, (*segment).GetPtr(), (*segment).GetSize(),
			 (*src)->GetName(), *src);
	      } else {
		HLTError("task %s (%p) failed releasing forwarded segment %p size %d of task %s (%p)",
			 pOwnerTask->GetName(), pOwnerTask, (*segment).GetPtr(), (*segment).GetSize(),
			 (*src)->GetName(), *src);
	      }
	    }
	    fForwardedSegments.erase(segment);
	    fForwardedSegmentSources.erase(src);
	    break;
	  }
	}
  } else {
    HLTError("inavalid parameter: pBlockDesc=%p pOwnerTask=%p", pBlockDesc, pOwnerTask);
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::Forward(AliHLTTask* pSrcTask, AliHLTComponentBlockData* pBlockDesc)
{
  // see header file for function documentation
  if (pSrcTask==NULL || pBlockDesc==NULL) return -EINVAL;
  assert(fForwardedSegments.size()==fForwardedSegmentSources.size());
  if (fForwardedSegments.size()!=fForwardedSegmentSources.size()) return -EFAULT;
  fForwardedSegmentSources.push_back(pSrcTask);
  fForwardedSegments.push_back(AliHLTDataSegment(pBlockDesc->fPtr, pBlockDesc->fOffset, pBlockDesc->fSize, pBlockDesc->fDataType, pBlockDesc->fSpecification));
  return 0;
}

AliHLTUInt8_t* AliHLTDataBuffer::GetTargetBuffer(int iMinSize)
{
  // see header file for function documentation
  AliHLTUInt8_t* pTargetBuffer=NULL;
  if (fpBuffer!=NULL) {
    HLTWarning("data buffer not properly reset, possible memory leak\n");
  }
  fpBuffer=CreateRawBuffer(iMinSize);
  if (fpBuffer) {
    pTargetBuffer=*fpBuffer;
  } else {
    HLTError("can not create raw buffer of size %d", iMinSize);
  }
  return pTargetBuffer;
}

int AliHLTDataBuffer::SetSegments(AliHLTUInt8_t* pTgt, AliHLTComponentBlockData* arrayBlockData, int iSize)
{
  // see header file for function documentation
  int iResult=0;
  if (pTgt && arrayBlockData && iSize>=0) {
    if (fpBuffer) {
      if (*fpBuffer==pTgt) {
	AliHLTDataBuffer::AliHLTDataSegment segment;
	AliHLTUInt32_t maxSize=0;
	for (int i=0; i<iSize; i++) {
	  // This function has to model the behavior of PubSub
	  // For output blocks only the fOffset value is used, this must be the offset
	  // relative to the output pointer. fPtr must be either NULL or the output
	  // pointer. In either case it is 'ignored' and set to the beginning of the
	  // data buffer
	  if (arrayBlockData[i].fPtr==NULL ||
	      arrayBlockData[i].fPtr==*fpBuffer) {
	    arrayBlockData[i].fPtr=*fpBuffer;
	    if ((arrayBlockData[i].fOffset+arrayBlockData[i].fSize<=fpBuffer->GetUsedSize()) ||
		((arrayBlockData[i].fOffset==~(AliHLTUInt32_t)0) && arrayBlockData[i].fSize==0)) {
	      segment.fSegmentOffset=arrayBlockData[i].fOffset;
	      segment.fPtr=(AliHLTUInt8_t*)arrayBlockData[i].fPtr;
	      segment.fSegmentSize=arrayBlockData[i].fSize;
	      segment.fDataType=arrayBlockData[i].fDataType;
	      segment.fSpecification=arrayBlockData[i].fSpecification;
	      fSegments.push_back(segment);
	      HLTDebug("set segment %s with size %d at offset %d", AliHLTComponent::DataType2Text(segment.fDataType).data(), segment.fSegmentSize, segment.fSegmentOffset);

	      // find the actual size of the data
	      if ((arrayBlockData[i].fOffset!=~(AliHLTUInt32_t)0) &&
		  arrayBlockData[i].fOffset+arrayBlockData[i].fSize>maxSize) {
		maxSize=arrayBlockData[i].fOffset+arrayBlockData[i].fSize;
	      }
	    } else {
	      HLTError("block data specification %#d (%s) exceeds size of data buffer", i, AliHLTComponent::DataType2Text(arrayBlockData[i].fDataType).data());
	      HLTError("block offset=%d, block size=%d, buffer size=%d", arrayBlockData[i].fOffset, arrayBlockData[i].fSize, fpBuffer->GetUsedSize());
	      iResult=-E2BIG;
	    }
	  } else {
	    HLTError("invalid pointer (%p) in block data specification (buffer %p size %d)."
		     "please note: for output blocks only the fOffset value is valid and must "
		     "be relative to the output buffer", arrayBlockData[i].fPtr, fpBuffer->GetPointer(), fpBuffer->GetUsedSize());
	    iResult=-ERANGE;
	  }
	}
	// to be enabled if unit test is ready
	iResult=SetRawBufferDataSize(fpBuffer, maxSize);	
      } else {
	HLTError("this data buffer (%p) does not match the internal data buffer %p of raw buffer %p", pTgt, fpBuffer->GetPointer(), fpBuffer);
	iResult=-EINVAL;
      }
    } else {
      HLTFatal("internal data structur mismatch");
      iResult=-EFAULT;
    }
  } else {
    HLTError("invalid parameter: pTgtBuffer=%p arrayBlockData=%p", pTgt, arrayBlockData);
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::IsEmpty()
{
  // see header file for function documentation
  int iResult=(fpBuffer==NULL && fForwardedSegments.size()==0) || GetNofSegments()==0;
  return iResult;
}

int AliHLTDataBuffer::GetNofSegments() const
{
  // see header file for function documentation
  int iResult=fSegments.size() + fForwardedSegments.size();
  return iResult;
}

int AliHLTDataBuffer::GetNofConsumers() const
{
  // see header file for function documentation
  int iResult=fConsumers.size() + GetNofActiveConsumers() + fReleasedConsumers.size();
  return iResult;
}

int AliHLTDataBuffer::GetNofPendingConsumers() const
{
  // see header file for function documentation
  int iResult=fConsumers.size();
  return iResult;
}

int AliHLTDataBuffer::GetNofActiveConsumers() const
{
  // see header file for function documentation
  int iResult=fActiveConsumers.size();
  return iResult;
}

AliHLTDataBuffer::AliHLTRawBuffer* AliHLTDataBuffer::CreateRawBuffer(AliHLTUInt32_t size)
{
  // see header file for function documentation
  AliHLTRawBuffer* pRawBuffer=NULL;
  unsigned int reqSize=size+fgkSafetyPatternSize;
  pRawBuffer=AliHLTDataBuffer::AliHLTRawPage::GlobalAlloc(reqSize);
  if (pRawBuffer) {
    pRawBuffer->UseBuffer(size);
  }
  if (pRawBuffer!=NULL && fgkSafetyPatternSize>0) {
    //fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "writing safety pattern to %p offset %d", (*buffer)->GetPointer(), (*buffer)->GetUsedSize());
    pRawBuffer->WritePattern(fgkSafetyPattern, fgkSafetyPatternSize);
  }
  return pRawBuffer;
}

int AliHLTDataBuffer::SetRawBufferDataSize(AliHLTRawBuffer* pBuffer, AliHLTUInt32_t size) const
{
  // see header file for function documentation
  int iResult=0;
  if (!pBuffer) return -EINVAL;
  if (size>pBuffer->GetUsedSize()) {
    HLTError("indicated data size %d exceeds data buffer %p (%d)", size, pBuffer->GetPointer(), pBuffer->GetUsedSize());
    return -ENOSPC;
  }
  if (fgkSafetyPatternSize>0) {
    if (pBuffer->CheckPattern(fgkSafetyPattern, fgkSafetyPatternSize)) {
      HLTError("potential memory corruption: component has written beyond end of data buffer %p size %d", pBuffer->GetPointer(), pBuffer->GetUsedSize());
    }
  }
  // shrink the buffer and write new pattern at the end
  AliHLTDataBuffer::AliHLTRawPage* rawpage=AliHLTDataBuffer::AliHLTRawPage::FindPage(pBuffer);
  if (rawpage) {
    pBuffer->UseBuffer(size);
    if (rawpage->SetSize(pBuffer, size+fgkSafetyPatternSize)==0) {
      // nothing to do
    } else {
      fgLogging.Logging(kHLTLogError, "AliHLTDataBuffer::SetRawBufferDataSize", "data buffer handling", "failed to set size for raw buffer %p", pBuffer);
      iResult=-EFAULT;
    }
  } else {
    fgLogging.Logging(kHLTLogError, "AliHLTDataBuffer::SetRawBufferDataSize", "data buffer handling", "can not find raw page for buffer %p", pBuffer);
    iResult=-ENOENT;
  }
  if (fgkSafetyPatternSize>0) {
    pBuffer->WritePattern(fgkSafetyPattern, fgkSafetyPatternSize);
  }
  return iResult;
}

int AliHLTDataBuffer::ReleaseRawBuffer(AliHLTRawBuffer* pBuffer)
{
  // see header file for function documentation
  int iResult=0;
  if (pBuffer) {
    AliHLTDataBuffer::AliHLTRawPage* rawpage=AliHLTDataBuffer::AliHLTRawPage::FindPage(pBuffer);
    if (rawpage)
    {
      if (fgkSafetyPatternSize>0) {
	//fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "comparing safety pattern at %p offset %d", pBuffer->GetPointer(), reinterpret_cast<AliHLTUInt32_t>(pBuffer));
	if ((pBuffer)->CheckPattern(fgkSafetyPattern, fgkSafetyPatternSize)) {
	  fgLogging.Logging(kHLTLogFatal, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "component has written beyond end of data buffer %p size %d", pBuffer->GetPointer(), pBuffer->GetUsedSize());
	}
      }
      pBuffer->Reset();
      if (rawpage->Free(pBuffer)==0) {
      } else {
	fgLogging.Logging(kHLTLogError, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "failed to release raw buffer %p", pBuffer);
      }
    } else {
      fgLogging.Logging(kHLTLogError, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "can not find raw page for buffer %p", pBuffer);
      iResult=-ENOENT;
    }
  } else {
    fgLogging.Logging(kHLTLogError, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "invalid parameter");
    iResult=-EINVAL;
  }
  return iResult;
}


int AliHLTDataBuffer::DeleteRawBuffers() 
{
  // see header file for function documentation
  int iResult=0;
#ifdef ALIHLTSYSTEM_PROFILING
  int iTotalSize=0;
  int iCount=fgFreeBuffers.size()+fgActiveBuffers.size();
#endif //ALIHLTSYSTEM_PROFILING
  AliHLTRawBufferPList::iterator buffer;;
  while ((buffer=fgFreeBuffers.begin())!=fgFreeBuffers.end()) {
#ifdef ALIHLTSYSTEM_PROFILING
    iTotalSize+=(*buffer)->GetTotalSize();
#endif //ALIHLTSYSTEM_PROFILING
    delete *buffer;
    fgFreeBuffers.erase(buffer);
  }
  while ((buffer=fgActiveBuffers.begin())!=fgActiveBuffers.end()) {
#ifdef ALIHLTSYSTEM_PROFILING
    iTotalSize+=(*buffer)->GetTotalSize();
#endif //ALIHLTSYSTEM_PROFILING
    fgLogging.Logging(kHLTLogWarning, "AliHLTDataBuffer::DeleteRawBuffer", "data buffer handling", "request to delete active raw buffer container (raw buffer %p, size %d)", (*buffer)->GetPointer(), (*buffer)->GetTotalSize());
    delete *buffer;
    fgActiveBuffers.erase(buffer);
  }
#ifdef ALIHLTSYSTEM_PROFILING
  fgLogging.Logging(kHLTLogImportant, "AliHLTDataBuffer::DeleteRawBuffer", "data buffer handling", "Total memory allocation: %d byte in %d buffers", iTotalSize, iCount);
#endif //ALIHLTSYSTEM_PROFILING
  return iResult;
}

int AliHLTDataBuffer::PrintStatistics() 
{
  // see header file for function documentation
  int iResult=0;
  int nofPages=0;
  AliHLTUInt32_t totalSize=0;
  for (AliHLTDataBuffer::AliHLTRawPage* rawpage=AliHLTDataBuffer::AliHLTRawPage::NextPage(NULL);
       rawpage!=NULL; 
       rawpage=AliHLTDataBuffer::AliHLTRawPage::NextPage(rawpage)) {
    nofPages++;
    totalSize+=rawpage->Size();
    if (fgLogging.CheckFilter(kHLTLogDebug)) rawpage->Print("");
  }
  //if (rawpage) rawpage->Print("global");
  fgLogging.Logging(kHLTLogInfo, "AliHLTDataBuffer::PrintStatistics", "data buffer handling", "total number of memory pages: %d   total size %d", nofPages, totalSize);

  return iResult;
}

AliHLTConsumerDescriptor* AliHLTDataBuffer::FindConsumer(const AliHLTComponent* pConsumer, AliHLTConsumerDescriptorPList &list) const
{
  // see header file for function documentation
  AliHLTConsumerDescriptor* pDesc=NULL;
  AliHLTConsumerDescriptorPList::iterator desc=list.begin();
  while (desc!=list.end() && pDesc==NULL) {
    if ((pConsumer==NULL || (*desc)->GetComponent()==pConsumer)) {
      pDesc=*desc;
    }
    desc++;
  }
  return pDesc;
}

int AliHLTDataBuffer::ResetDataBuffer() 
{
  // see header file for function documentation
  int iResult=0;
  AliHLTRawBuffer* pBuffer=fpBuffer;
  fpBuffer=NULL;

  // cleanup forwarded segment lists
  assert(fForwardedSegments.size()==0);
  fForwardedSegments.clear();
  fForwardedSegmentSources.clear();

  // cleanup consumer states
  AliHLTConsumerDescriptorPList::iterator desc;
//   if (GetNofPendingConsumers()>0) {
//     desc=fConsumers.begin();
//     while (desc!=fConsumers.end()) {
//       AliHLTComponent* pComp=(*desc)->GetComponent();
//       HLTError("internal error: consumer %p (%s %p) did not get data from data buffer %p", *desc, pComp?pComp->GetComponentID():"", pComp, this);
//       desc++;
//     }
//   }
  desc=fReleasedConsumers.begin();
  while (desc!=fReleasedConsumers.end()) {
    AliHLTConsumerDescriptor* pDesc=*desc;
    fReleasedConsumers.erase(desc);
    desc=fReleasedConsumers.begin();
    fConsumers.push_back(pDesc);
  }
  desc=fActiveConsumers.begin();
  while (desc!=fActiveConsumers.end()) {
    AliHLTConsumerDescriptor* pDesc=*desc;
    HLTWarning("consumer %p (%s) was not released", pDesc, pDesc->GetComponent()?pDesc->GetComponent()->GetComponentID():"### invalid component ###");
    fActiveConsumers.erase(desc);
    desc=fActiveConsumers.begin();
    fConsumers.push_back(pDesc);
  }

  // cleanup segments
  AliHLTDataSegmentList::iterator segment=fSegments.begin();
  while (segment!=fSegments.end()) {
    fSegments.erase(segment);
    segment=fSegments.begin();
  }

  // cleanup raw buffer
  if (pBuffer) {
    ReleaseRawBuffer(pBuffer);
  }
  return iResult;
}

int AliHLTDataBuffer::Reset()
{
  // see header file for function documentation
  return ResetDataBuffer();
}

// this is the version which works on lists of components instead of consumer descriptors
// int AliHLTDataBuffer::ChangeConsumerState(AliHLTComponent* pConsumer, AliHLTComponentPList &srcList, AliHLTComponentPList &tgtList)
// {
//   int iResult=0;
//   if (pDesc) {
//     AliHLTComponentPList::iterator desc=srcList.begin();
//     while (desc!=srcList.end()) {
//       if ((*desc)==pConsumer) {
// 	srcList.erase(desc);
// 	tgtList.push_back(pConsumer);
// 	break;
//       }
//      desc++;
//     }
//     if (desc==srcList.end()) {
//       HLTError("can not find consumer component %p in list", pConsumer);
//       iResult=-ENOENT;
//     }
//   } else {
//     HLTError("invalid parameter");
//     iResult=-EINVAL;
//   }
//   return iResult;
// }

int AliHLTDataBuffer::ChangeConsumerState(AliHLTConsumerDescriptor* pDesc, AliHLTConsumerDescriptorPList &srcList, AliHLTConsumerDescriptorPList &tgtList)
{
  // see header file for function documentation
  int iResult=-ENOENT;
  if (pDesc) {
    AliHLTConsumerDescriptorPList::iterator desc=srcList.begin();
    while (desc!=srcList.end()) {
      if ((*desc)==pDesc) {
	srcList.erase(desc);
	tgtList.push_back(pDesc);
	iResult=0;
	break;
      }
      desc++;
    }
    if (iResult<0) {
      HLTError("can not find consumer descriptor %p in list", pDesc);
    }
  } else {
    HLTError("invalid parameter");
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::CleanupConsumerList() 
{
  // see header file for function documentation
  int iResult=0;
  ResetDataBuffer();
  AliHLTConsumerDescriptorPList::iterator desc=fConsumers.begin();
  while (desc!=fConsumers.end()) {
    delete *desc;
    fConsumers.erase(desc);
    desc=fConsumers.begin();
  }
  return iResult;
}

int AliHLTDataBuffer::FindConsumer(const AliHLTComponent* pConsumer, int bAllLists)
{
  // see header file for function documentation
  AliHLTConsumerDescriptorPList::iterator desc=fConsumers.begin();
  while (desc!=fConsumers.end()) {
    if ((*desc)->GetComponent()==pConsumer)
      return 1;
    desc++;
  }
  if (bAllLists==0) return 0;

  desc=fActiveConsumers.begin();
  while (desc!=fActiveConsumers.end()) {
    if ((*desc)->GetComponent()==pConsumer)
      return 1;
    desc++;
  }
  desc=fReleasedConsumers.begin();
  while (desc!=fReleasedConsumers.end()) {
    if ((*desc)->GetComponent()==pConsumer)
      return 1;
    desc++;
  }
  return 0;
}

AliHLTDataBuffer::AliHLTRawBuffer::AliHLTRawBuffer(AliHLTUInt32_t size)
  : fSize(0)
  , fTotalSize(size)
  , fExternalPtr(NULL)
  , fPtr(static_cast<AliHLTUInt8_t*>(malloc(size)))
  , fLastEventCount(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  if (fPtr==NULL) {
    fSize=0;
    fTotalSize=0;
  }
}

AliHLTDataBuffer::AliHLTRawBuffer::AliHLTRawBuffer(AliHLTUInt32_t size, AliHLTUInt8_t* buffer)
  : fSize(0)
  , fTotalSize(size)
  , fExternalPtr(buffer)
  , fPtr(fExternalPtr)
  , fLastEventCount(0)
{
  // see header file for class documentation
}

AliHLTDataBuffer::AliHLTRawBuffer::~AliHLTRawBuffer()
{
  // see header file for class documentation
  if (fExternalPtr==NULL && fPtr) {
    free(fPtr);
  }
  fPtr=NULL;
  fSize=0;
  fTotalSize=0;
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator==(void* ptr) const
{
  // see header file for function documentation
  return fPtr == static_cast<AliHLTUInt8_t*>(ptr);
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator<(void* ptr) const
{
  // see header file for function documentation
  int iResult=fPtr < static_cast<AliHLTUInt8_t*>(ptr);
  //printf("%p: %p <= %p (%d)\n", this, fPtr, ptr, iResult);
  return iResult;
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator<=(void* ptr) const
{
  // see header file for function documentation
  int iResult=fPtr <= static_cast<AliHLTUInt8_t*>(ptr);
  //printf("%p: %p <= %p (%d)\n", this, fPtr, ptr, iResult);
  return iResult;
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator>(void* ptr) const
{
  // see header file for function documentation
  int iResult=fPtr+fSize > static_cast<AliHLTUInt8_t*>(ptr);
  //printf("%p: %p + %d > %p (%d)\n", this, fPtr, fSize, ptr, iResult);
  return iResult;
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator-(void* ptr) const
{
  // see header file for function documentation
  return static_cast<int>(static_cast<AliHLTUInt8_t*>(ptr)-fPtr);
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator<(const AliHLTRawBuffer& op) const
{
  // see header file for function documentation
  return (fPtr+fSize < op.fPtr);
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator<=(const AliHLTRawBuffer& op) const
{
  // see header file for function documentation
  return (fPtr+fSize <= op.fPtr);
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator>(const AliHLTRawBuffer& op) const
{
  // see header file for function documentation
  return (fPtr >= op.fPtr+op.fSize);
}

AliHLTUInt8_t* AliHLTDataBuffer::AliHLTRawBuffer::UseBuffer(AliHLTUInt32_t size)
{
  // mark a portion of the buffer as used
  if (fTotalSize>=size) {
    fSize=size;
    fLastEventCount=AliHLTDataBuffer::fgEventCount;
    // only return pointer if there is a portion of the buffer used
    if (size>0) return fPtr;
  }
  return NULL;
}

AliHLTDataBuffer::AliHLTRawBuffer* AliHLTDataBuffer::AliHLTRawBuffer::Split(AliHLTUInt32_t size)
{
  // split a buffer at specified size
  // only possible for buffers with external memory
  if (fTotalSize>size && 
      (fSize==0 || fSize<=size) && // used size must fit into the first part
      fExternalPtr!=NULL) {
    AliHLTRawBuffer* part2=new AliHLTRawBuffer(fTotalSize-size, fPtr+size);
    if (part2) {
      fTotalSize=size;
    }
    return part2;
  }
  return NULL;
}

int AliHLTDataBuffer::AliHLTRawBuffer::CheckSize(AliHLTUInt32_t size) const
{
  // see header file for function documentation
  if (fTotalSize<size) return 0;
  unsigned adjust=0;
  if (fLastEventCount+1<AliHLTDataBuffer::fgEventCount) {
    adjust=AliHLTDataBuffer::fgEventCount-fLastEventCount;
  }
  return (adjust>2) || ((fTotalSize-size)<(fgMargin<<adjust));
}

int AliHLTDataBuffer::AliHLTRawBuffer::Reset()
{
  // see header file for function documentation
  fSize=0;
  return 0;
}

int AliHLTDataBuffer::AliHLTRawBuffer::WritePattern(const char* pattern, int size)
{
  // see header file for function documentation
  int iResult=0;
  if (pattern!=NULL && size>0) {
    if (fSize+size<=fTotalSize) {
      memcpy(((char*)fPtr)+fSize, pattern, size);
      iResult=size;
    } else {
      iResult=-ENOSPC;
    }
  }
  return iResult;
}

int AliHLTDataBuffer::AliHLTRawBuffer::CheckPattern(const char* pattern, int size) const
{
  // see header file for function documentation
  int iResult=0;
  if (pattern!=NULL && size>0) {
    if (fSize+size<=fTotalSize) {
      iResult=memcmp(((char*)fPtr)+fSize, pattern, size)!=0;
    } else {
      iResult=-ENOSPC;
    }
  }
  return iResult;
}

int AliHLTDataBuffer::AliHLTRawBuffer::Merge(const AliHLTDataBuffer::AliHLTRawBuffer& neighbor)
{
  // Merge buffer with neighboring buffer.
  // Only possible if the buffers are consecutive with out any gap.

  if (!fExternalPtr || !neighbor.fExternalPtr) return -EPERM;

  if (neighbor.fTotalSize==0 &&
      fPtr < neighbor.fPtr &&
      fPtr+fTotalSize > neighbor.fPtr) {
    // special case for a buffer of zero size embedded into this buffer
    // nothing to do
    return 0;
  }
  if (fTotalSize==0 &&
      neighbor.fPtr < fPtr &&
      neighbor.fPtr+neighbor.fTotalSize > fPtr) {
    // special case for this buffer of size zero embedded into another buffer
    fPtr=neighbor.fPtr;
    fExternalPtr=fPtr;
    fTotalSize+=neighbor.fTotalSize;
    fSize=0;
    return 0;
  }
  if (fPtr+fTotalSize == neighbor.fPtr) {
    fTotalSize+=neighbor.fTotalSize;
    fSize=0;
    return 0;
  }
  if (fPtr == neighbor.fPtr+neighbor.fTotalSize) {
    fPtr=neighbor.fPtr;
    fExternalPtr=fPtr;
    fTotalSize+=neighbor.fTotalSize;
    fSize=0;
    return 0;
  }
  return -EINVAL;
}

void AliHLTDataBuffer::AliHLTRawBuffer::Print(const char* option) const
{
  /// print buffer information
  if (strcmp(option, "min")!=0) {
    cout << "************* AliHLTRawBuffer status ***********" << endl;
  }
  printf("  %p: buffer %p%s size %d used %d\n", this, fPtr, fExternalPtr?" (external)":"", fTotalSize, fSize); fflush(stdout);
}

AliHLTDataBuffer::AliHLTRawPage::AliHLTRawPage(AliHLTUInt32_t pagesize)
  : fSize(pagesize)
  , fPtr(static_cast<AliHLTUInt8_t*>(malloc(pagesize)))
  , fFreeBuffers()
  , fUsedBuffers()
{
  // constructor
  if (fPtr) {
    fFreeBuffers.push_back(new AliHLTRawBuffer(fSize, fPtr));
  } else {
    fSize=0;
  }
}

AliHLTDataBuffer::AliHLTRawPage::~AliHLTRawPage()
{
  // destructor
  if (IsUsed()) {
    // do not free if the resources have not been completely freed
    HLTError("memory mismatch: not all allocated intances have been released");
  } else {
    if (IsFragmented()) {
      HLTWarning("page still fragmented");
    }
    AliHLTRawBufferPList::iterator element=fFreeBuffers.begin();
    while (element!=fFreeBuffers.end()) {
      if (*element) delete *element;
      element=fFreeBuffers.erase(element);
    }
    if (fPtr) {
      free(fPtr);
    }
    fPtr=NULL;
    fSize=0;
  }
}

AliHLTDataBuffer::AliHLTRawBuffer* AliHLTDataBuffer::AliHLTRawPage::Alloc(AliHLTUInt32_t size)
{
  /// alloc a buffer of specified size
  if (fFreeBuffers.size()==0) return NULL;
  
  for (AliHLTRawBufferPList::iterator iter=fFreeBuffers.begin();
       iter!=fFreeBuffers.end();
       iter++) {
    if ((*iter)->GetTotalSize()==size) {
      AliHLTRawBuffer* thisbuffer=*iter;
      fFreeBuffers.erase(iter);
      fUsedBuffers.push_back(thisbuffer);
      return thisbuffer;
    } else if ((*iter)->GetTotalSize()>size) {
      AliHLTRawBuffer* thisbuffer=*iter;
      AliHLTRawBuffer* newbuffer=thisbuffer->Split(size);
      if (newbuffer) {
	*iter=newbuffer;
	fUsedBuffers.push_back(thisbuffer);
	return thisbuffer;
      } else {
	HLTWarning("failed to alloc raw buffer: cannot split raw buffer %p of size %d (used %d) at size %d", *iter, (*iter)->GetTotalSize(), (*iter)->GetUsedSize(), size);
      }
    }
  }
  return NULL;
}

int AliHLTDataBuffer::AliHLTRawPage::Free(AliHLTRawBuffer* pBuffer)
{
  /// free a buffer and merge consecutive free buffers
  int iResult=0;
  for (AliHLTRawBufferPList::iterator iter=fUsedBuffers.begin();
       iter!=fUsedBuffers.end() && iResult>=0;
       iter++) {
    if ((*iter)==pBuffer) {
      fUsedBuffers.erase(iter);
      AliHLTRawBufferPList::iterator prev=fFreeBuffers.begin();
      for (; prev!=fFreeBuffers.end() && iResult>=0; prev++) {
	if ((*pBuffer)<(*(*prev)) ||
	    ((*prev)->GetTotalSize()==0 && pBuffer->GetPointer()<=(*prev)->GetPointer() && (*prev)->GetPointer()<=pBuffer->GetPointer()+pBuffer->GetTotalSize())) {
	  // check consecutive buffers
	  if ((*(*prev)) == (pBuffer->GetPointer()+pBuffer->GetTotalSize()) ||
	      ((*prev)->GetTotalSize()==0 && pBuffer->GetPointer()<=(*prev)->GetPointer() && (*prev)->GetPointer()<=pBuffer->GetPointer()+pBuffer->GetTotalSize())) {
	    // the buffer to be released has a consecutive free buffer -> merge them
	    if ((iResult=pBuffer->Merge(*(*prev)))>=0) {
	      delete *prev;
	      *prev=pBuffer;
	    } else {
	      HLTError("failed to merge consecutive/overlapping buffers %p and %p", pBuffer, (*prev));
	      pBuffer->Print("");
	      (*prev)->Print("");
	    }
	    break;
	  }
	  fFreeBuffers.insert(prev, pBuffer);
	  break;
	}
	if ((*pBuffer)>(*(*prev)) ||
	    (pBuffer->GetTotalSize()==0 && (*prev)->GetPointer()<=pBuffer->GetPointer() && pBuffer->GetPointer()<=(*prev)->GetPointer()+(*prev)->GetTotalSize())) {
	  // check consecutive buffers
	  if ((*pBuffer) == ((*prev)->GetPointer()+(*prev)->GetTotalSize())||
	      (pBuffer->GetTotalSize()==0 && (*prev)->GetPointer()<=pBuffer->GetPointer() && pBuffer->GetPointer()<=(*prev)->GetPointer()+(*prev)->GetTotalSize())) {
	    // the buffer to be released is consecutive to a free buffer -> merge them
	    if ((iResult=pBuffer->Merge(*(*prev)))>=0) {
	      AliHLTRawBufferPList::iterator succ=prev+1;
	      delete *prev;
	      *prev=pBuffer;
	      // check if the buffer and the following one are consecutive
	      if (succ!=fFreeBuffers.end() &&
		  (*(*succ)) == (pBuffer->GetPointer()+pBuffer->GetTotalSize())) {
		if ((iResult=pBuffer->Merge(*(*succ)))>=0) {
		  delete *succ;
		  fFreeBuffers.erase(succ);
		}
	      }
	    }
	    break;
	  }
	}
      }
      if (prev==fFreeBuffers.end()) {
	fFreeBuffers.push_back(pBuffer);
      }

      // merge consecutive free buffers
      prev=fFreeBuffers.begin();
      for (AliHLTRawBufferPList::iterator current=prev+1; current!=fFreeBuffers.end() && iResult>=0; ) {
	// check if the buffer is embedded into the previous one
	if ((*current)->GetTotalSize()==0 && (*prev)->GetPointer()<=(*current)->GetPointer() && (*current)->GetPointer()<(*prev)->GetPointer()+(*prev)->GetTotalSize())  {
	  if ((iResult=(*prev)->Merge(*(*current)))>=0) {
	    current=fFreeBuffers.erase(current);
	    continue;
	  } else {
	    HLTError("failed to merge embedded zero length buffer into preceeding buffer");
	    Print("");
	  }
	}
	// check if the buffer is consecutive to the previous one
	if ((*(*current)) == ((*prev)->GetPointer()+(*prev)->GetTotalSize())) {
	  if ((iResult=(*prev)->Merge(*(*current)))>=0) {
	    current=fFreeBuffers.erase(current);
	    continue;
	  } else {
	    HLTError("failed to merge consecutive free buffers");
	    Print("");
	  }
	}
	prev=current++;
      }

      // buffer was part of this page
      return 0;
    }
  }
  // buffer not found in this page
  return 1;
}

int AliHLTDataBuffer::AliHLTRawPage::SetSize(const AliHLTDataBuffer::AliHLTRawBuffer* pBuffer, AliHLTUInt32_t size)
{
  /// set the size of a raw buffer and release the remaining part
  int iResult=0;
  for (AliHLTRawBufferPList::iterator iter=fUsedBuffers.begin();
       iter!=fUsedBuffers.end() && iResult>=0;
       iter++) {
    if ((*iter)==pBuffer) {      // buffer was part of this page
      if ((*iter)->GetTotalSize()==size) return 0;
      if ((*iter)->GetTotalSize()<size) {
	HLTError("%d exceeds total size of buffer %p (%d used %d)\n", size, *iter, (*iter)->GetTotalSize(), (*iter)->GetUsedSize());
	return -ENOSPC;
      }
      AliHLTDataBuffer::AliHLTRawBuffer* freespace=(*iter)->Split(size);
      if (freespace) {
	fUsedBuffers.push_back(freespace);
	Free(freespace);
      } else {
	HLTWarning("failed to relase unused memory: cannot split raw buffer %p of size %d (used %d) at size %d", *iter, (*iter)->GetTotalSize(), (*iter)->GetUsedSize(), size);
      }
      return 0;
    }
  }
  // buffer not found in this page
  return 1;
}

bool AliHLTDataBuffer::AliHLTRawPage::HasBuffer(const AliHLTDataBuffer::AliHLTRawBuffer* pBuffer)
{
  /// check if the buffer is in this page
  for (AliHLTRawBufferPList::iterator iter=fUsedBuffers.begin();
       iter!=fUsedBuffers.end();
       iter++) {
    if ((*iter)==pBuffer) {      // buffer was part of this page
      return true;
    }
  }
  // buffer not found in this page
  return false;
}

AliHLTUInt32_t AliHLTDataBuffer::AliHLTRawPage::Capacity() const 
{
  /// get max available contiguous buffer
  AliHLTUInt32_t capacity=0;
  for (unsigned i=0; i<fFreeBuffers.size(); i++) {
    if (fFreeBuffers[i]->GetTotalSize()>capacity) 
      capacity=fFreeBuffers[i]->GetTotalSize();
  }
  return capacity;
}

void AliHLTDataBuffer::AliHLTRawPage::Print(const char* option)
{
  /// print page information
  if (strcmp(option, "global")==0) {
    cout << "number of global pages: " << fgGlobalPages.size() << endl;
    for (AliHLTRawPage* rawpage=NextPage(NULL);
	 rawpage!=NULL; 
	 rawpage=NextPage(rawpage)) {
      rawpage->Print("");
    }
    return;
  }
  cout << "************* AliHLTRawPage status ***********" << endl;
  cout << "  instance " << this << endl;
  printf("  buffer %p  size %d", fPtr, fSize);
  cout << "  used buffers: " << fUsedBuffers.size() << endl;
  AliHLTRawBufferPList::iterator iter=fUsedBuffers.begin();
  for (; iter!=fUsedBuffers.end(); iter++) {
    cout << "  "; (*iter)->Print("min");
  }
  cout << "  free buffers: " << fFreeBuffers.size() << endl;
  iter=fFreeBuffers.begin();
  for (; iter!=fFreeBuffers.end(); iter++) {
    cout << "  "; (*iter)->Print("min");
  }
}


vector<AliHLTDataBuffer::AliHLTRawPage*> AliHLTDataBuffer::AliHLTRawPage::fgGlobalPages;

AliHLTUInt32_t AliHLTDataBuffer::AliHLTRawPage::fgGlobalPageSize=30*1024*1024;

unsigned int AliHLTDataBuffer::GetMaxBufferSize()
{
	return(30 * AliHLTRawPage::GetGlobalPageSize()); //Use a max data size of 30 * GlobalPageSize, as this is the maximum AliHLTRawPage::GlobalAlloc will allocate
}

AliHLTDataBuffer::AliHLTRawBuffer* AliHLTDataBuffer::AliHLTRawPage::GlobalAlloc(AliHLTUInt32_t size, int verbosity)
{
  // alloc a buffer of specified size from the global pages
  AliHLTDataBuffer::AliHLTRawBuffer* rawbuffer=NULL;
  vector<AliHLTDataBuffer::AliHLTRawPage*>::iterator page=fgGlobalPages.begin();
  AliHLTLogging log;
  for (page=fgGlobalPages.begin();page!=fgGlobalPages.end(); page++) {
    if ((rawbuffer=(*page)->Alloc(size))!=NULL) {
      if (verbosity>1) {
	log.Logging(kHLTLogInfo, "AliHLTDataBuffer::AliHLTRawPage::GlobalAlloc", "data buffer handling", "allocated raw buffer %p from page %p\n", rawbuffer, *page);
	rawbuffer->Print("min");
      }
      break;
    }
  }
  if (!rawbuffer) {
    AliHLTUInt32_t rawPageSize=fgGlobalPageSize;
    if (rawPageSize<size) {
      if (rawPageSize*30+fgkSafetyPatternSize<size) {
	log.Logging(kHLTLogError, "AliHLTDataBuffer::AliHLTRawPage::GlobalAlloc", "data buffer handling", "refusing to allocate buffer of size %d", size);
	return NULL;
      }
      rawPageSize=size;
    }
    AliHLTDataBuffer::AliHLTRawPage* rawpage=new AliHLTDataBuffer::AliHLTRawPage(rawPageSize);
    if (!rawpage) {
      log.Logging(kHLTLogError, "AliHLTDataBuffer::AliHLTRawPage::GlobalAlloc", "data buffer handling", "can not create raw page");
      return NULL;
    }

    // check is there is at least one unused page which can be replaced by the newly created one
    for (page=fgGlobalPages.begin(); page!=fgGlobalPages.end(); page++) {
      if ((*page)->IsUsed()) continue;
      delete *page;
      fgGlobalPages.erase(page);
      break; // delete only one page to be replaced by the new page
    }
    fgGlobalPages.push_back(rawpage);
    if ((rawbuffer=rawpage->Alloc(size))!=NULL) {
      if (verbosity>1) {
	log.Logging(kHLTLogInfo, "AliHLTDataBuffer::AliHLTRawPage::GlobalAlloc", "data buffer handling", "allocated raw buffer %p from page %p\n", rawbuffer, rawpage);
	rawbuffer->Print("min");
      }
    }
  }

  return rawbuffer;
}

AliHLTDataBuffer::AliHLTRawPage* AliHLTDataBuffer::AliHLTRawPage::FindPage(AliHLTDataBuffer::AliHLTRawBuffer* buffer)
{
  // find buffer in the global pages
  vector<AliHLTDataBuffer::AliHLTRawPage*>::iterator page=fgGlobalPages.begin();
  for (; page!=fgGlobalPages.end(); page++) {
    if ((*page)->HasBuffer(buffer)) {
      return *page;
    }
  }

  return NULL;
}

int AliHLTDataBuffer::AliHLTRawPage::GlobalClean()
{
  // cleanup the global pages */
  vector<AliHLTDataBuffer::AliHLTRawPage*>::iterator page=fgGlobalPages.begin();
  while (page!=fgGlobalPages.end()) {
    if (!(*page)->IsUsed()) {
      delete *page;
      page=fgGlobalPages.erase(page);
      continue;
    }
    AliHLTLogging log;
    log.Logging(kHLTLogError, "AliHLTDataBuffer::AliHLTRawPage::GlobalClean", "data buffer handling", "HLT memory page still in use, skipping cleanup, potential memory leak");
    
    page++;
  }
  
  return 0;
}

AliHLTDataBuffer::AliHLTRawPage* AliHLTDataBuffer::AliHLTRawPage::NextPage(const AliHLTDataBuffer::AliHLTRawPage* prev)
{
  // get next global page
  vector<AliHLTDataBuffer::AliHLTRawPage*>::iterator page=fgGlobalPages.begin();
  for (; page!=fgGlobalPages.end(); page++) {
    if (prev==NULL) return *page;
    if (*page!=prev) continue;
    if (++page!=fgGlobalPages.end()) return *page;
    break;
  }
  return NULL;
}

void AliHLTDataBuffer::AliHLTDataSegment::Print(const char* /*option*/) const
{
  // print info for data segment
  cout << "AliHLTDataSegment " << this 
       << " " << AliHLTComponent::DataType2Text(fDataType)
       << " " << hex << fSpecification << dec
       << " Ptr " << (void*)fPtr
       << " offset " << fSegmentOffset
       << " size " << fSegmentSize
       << endl;
}

void AliHLTDataBuffer::AliHLTForwardedDataSegment::Print(const char* option) const
{
  // print info for data segment
  cout << "AliHLTForwardeDataSegment " << this << endl;
  cout << "    my    : "; AliHLTDataSegment::Print(option);
  cout << "    parent: "; fParentSegment.Print(option);
  cout << "    task  : "; 
  if (fParentTask) fParentTask->Print("");
  else cout << "nil" << endl;
}

void AliHLTDataBuffer::Print(const char* option) const
{
  // print info for data buffer
  unsigned i=0;
  cout << "AliHLTDataBuffer " << this << endl;
  cout << " raw buffer " << fpBuffer << endl;
  if (fpBuffer) {
    cout << " ";
    fpBuffer->Print(option);
  }

  cout << " total segments: " << GetNofSegments() << endl;
  cout << "   data segments: " << fSegments.size() << endl;
  for (i=0; i<fSegments.size(); i++) {
    cout << "     ";
    fSegments[i].Print(option);
  }

  cout << "   forwarded segments: " << fForwardedSegments.size() << endl;
  for (i=0; i<fForwardedSegments.size(); i++) {
    cout << "     ";
    fForwardedSegments[i].Print(option);
  }

  cout << " consumers: " << GetNofConsumers() << endl;
  for (i=0; i<fConsumers.size(); i++) {
    cout << "   ";
    fConsumers[i]->Print(option);
  }

  cout << " active consumers: " << GetNofActiveConsumers() << endl;
  for (i=0; i<fActiveConsumers.size(); i++) {
    cout << "   ";
    fActiveConsumers[i]->Print(option);
  }

  cout << " released consumers: " << fReleasedConsumers.size() << endl;
  for (i=0; i<fReleasedConsumers.size(); i++) {
    cout << "   ";
    fReleasedConsumers[i]->Print(option);
  }

}
