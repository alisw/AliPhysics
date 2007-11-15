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

/** @file   AliHLTDataBuffer.cxx
    @author Matthias Richter
    @date   
    @brief  Handling of Data Buffers for HLT components.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTDataBuffer.h"
#include "AliHLTConsumerDescriptor.h"
#include "AliHLTComponent.h"
#include "AliHLTTask.h"
#include <cerrno>
#include <cassert>
//#include <string>
//#include "AliHLTSystem.h"

typedef vector<AliHLTDataBuffer::AliHLTDataSegment> AliHLTDataSegmentList;
typedef vector<AliHLTDataBuffer::AliHLTRawBuffer*>  AliHLTRawBufferPList;

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
AliHLTRawBufferPList AliHLTDataBuffer::fgFreeBuffers;
AliHLTRawBufferPList AliHLTDataBuffer::fgActiveBuffers;
AliHLTUInt32_t AliHLTDataBuffer::fgMargin=1024;
AliHLTLogging AliHLTDataBuffer::fgLogging;
const Int_t AliHLTDataBuffer::fgkSafetyPatternSize=16;
const char AliHLTDataBuffer::fgkSafetyPattern[]={0x28, 0x63, 0x29, 0x4d, 0x52, 0x49, 0x43, 0x48, 0x54, 0x45, 0x52, 0x20, 0x32, 0x30, 0x30, 0x37};

AliHLTDataBuffer::~AliHLTDataBuffer()
{
  // see header file for function documentation
  if (--fgNofInstances<=0) {
    DeleteRawBuffers();
  }
  CleanupConsumerList();
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

int AliHLTDataBuffer::FindMatchingDataSegments(const AliHLTComponent* pConsumer, vector<AliHLTDataBuffer::AliHLTDataSegment>& tgtList)
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
  
  if (pConsumer) {
    AliHLTComponentDataTypeList dtlist;
    ((AliHLTComponent*)pConsumer)->GetInputDataTypes(dtlist);
    AliHLTDataSegmentList::iterator segment=fSegments.begin();
    while (segment!=fSegments.end()) {
      AliHLTComponentDataTypeList::iterator type=dtlist.begin();
      while (type!=dtlist.end()) {
	if ((*segment).fDataType==(*type) ||
	    (*type)==kAliHLTAnyDataType) {
	  tgtList.push_back(*segment);
	  iResult++;
	  break;
	}
	type++;
      }
      segment++;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::Subscribe(const AliHLTComponent* pConsumer, AliHLTComponentBlockData* arrayBlockDesc, int iArraySize)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumer && arrayBlockDesc) {
    if (1/*fpBuffer*/) {
      AliHLTConsumerDescriptor* pDesc=FindConsumer(pConsumer, fConsumers);
      if (pDesc) {
	AliHLTDataSegmentList tgtList;
	// Matthias 26.07.2007 AliHLTSystem should behave the same way as PubSub
	// so it does not matter if there are matching data types or not, unless
	// we implement such a check in PubSub
	if ((iResult=FindMatchingDataSegments(pConsumer, tgtList))>=0) {
	  int i =0;
	  AliHLTDataSegmentList::iterator segment=tgtList.begin();
	  while (segment!=tgtList.end() && i<iArraySize) {
	    // fill the block data descriptor
	    arrayBlockDesc[i].fStructSize=sizeof(AliHLTComponentBlockData);
	    // the shared memory key is not used in AliRoot
	    arrayBlockDesc[i].fShmKey.fStructSize=sizeof(AliHLTComponentShmData);
	    arrayBlockDesc[i].fShmKey.fShmType=gkAliHLTComponentInvalidShmType;
	    arrayBlockDesc[i].fShmKey.fShmID=gkAliHLTComponentInvalidShmID;
	    // This models the behavior of PubSub.
	    // For incoming data blocks, fOffset must be ignored by the
	    // processing component. It is set for bookkeeping in the framework.
	    // fPtr always points to the beginning of the data.
	    arrayBlockDesc[i].fOffset=0;
	    AliHLTUInt8_t* pTgt=*segment;
	    arrayBlockDesc[i].fPtr=reinterpret_cast<void*>(pTgt);
	    arrayBlockDesc[i].fSize=(*segment).fSegmentSize;
	    arrayBlockDesc[i].fDataType=(*segment).fDataType;
	    arrayBlockDesc[i].fSpecification=(*segment).fSpecification;
	    pDesc->SetActiveDataSegment(*segment);
	    HLTDebug("component %p (%s) subscribed to segment #%d offset %d size %d data type %s %#x", 
		     pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID(), i, arrayBlockDesc[i].fOffset,
		     arrayBlockDesc[i].fSize, (AliHLTComponent::DataType2Text(arrayBlockDesc[i].fDataType)).c_str(), 
		     arrayBlockDesc[i].fSpecification);
	    i++;
	    segment++;
	  }
	  // check whether there was enough space for the segments
	  if (i!=(int)tgtList.size()) {
	    HLTError("too little space in block descriptor array: required %d, available %d", tgtList.size(), iArraySize);
	    iResult=-ENOSPC;
	  } else {
	  // move this consumer to the active list
	  if (i==0) {
	    ChangeConsumerState(pDesc, fConsumers, fReleasedConsumers);
	    HLTDebug("no input data for component %p (%s) available", pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID());
	  } else if (ChangeConsumerState(pDesc, fConsumers, fActiveConsumers)>=0) {
	    HLTDebug("component %p (%s) subscribed to data buffer %p", pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID(), this);
	  } else {
	    // TODO: cleanup the consumer descriptor correctly
	    memset(arrayBlockDesc, 0, iArraySize*sizeof(AliHLTComponentBlockData));
	    HLTError("can not activate consumer %p for data buffer %p", pConsumer, this);
	    iResult=-EACCES;
	  }
	  }
	} else {
	  HLTError("unresolved data segment(s) for component %p (%s)", pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID());
	  iResult=-EBADF;
	}
      } else {
	HLTError("component %p is not a data consumer of data buffer %s", pConsumer, this);
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
	HLTWarning("data segment missmatch, component %p has not subscribed to a segment with offset %#x and size %d", pConsumer, pBlockDesc->fOffset, pBlockDesc->fSize);
	// TODO: appropriate error handling, but so far optional
	iResult=0;
      } else {
	pDesc->ReleaseActiveDataSegment(AliHLTDataSegment(pBlockDesc->fPtr, pBlockDesc->fOffset, pBlockDesc->fSize));
      }
      if (GetNofPendingConsumers()==0 && fForwardedSegments.size()>0) {
	// last consumer, release forwarded segments
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

int AliHLTDataBuffer::Forward(AliHLTTask* pSrcTask, AliHLTComponentBlockData* pBlockDesc)
{
  // see header file for function documentation
  if (pSrcTask==NULL || pBlockDesc==NULL) return -EINVAL;
  assert(fForwardedSegments.size()==fForwardedSegmentSources.size());
  if (fForwardedSegments.size()!=fForwardedSegmentSources.size()) return -EFAULT;
  fForwardedSegmentSources.push_back(pSrcTask);
  fForwardedSegments.push_back(AliHLTDataSegment(pBlockDesc->fPtr, pBlockDesc->fOffset, pBlockDesc->fSize));
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
    HLTError("can not create raw buffer");
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
	for (int i=0; i<iSize; i++) {
	  // This function has to model the behavior of PubSub
	  // For output blocks only the fOffset value is used, this must be the offset
	  // relative to the output pointer. fPtr must be either NULL or the output
	  // pointer. In either case it is 'ignored' and set to the beginning of the
	  // data buffer
	  if (arrayBlockData[i].fPtr==NULL ||
	      arrayBlockData[i].fPtr==*fpBuffer) {
	    arrayBlockData[i].fPtr=*fpBuffer;
	    if (arrayBlockData[i].fOffset+arrayBlockData[i].fSize<=fpBuffer->GetUsedSize()) {
	      segment.fSegmentOffset=arrayBlockData[i].fOffset;
	      segment.fPtr=(AliHLTUInt8_t*)arrayBlockData[i].fPtr;
	      segment.fSegmentSize=arrayBlockData[i].fSize;
	      segment.fDataType=arrayBlockData[i].fDataType;
	      segment.fSpecification=arrayBlockData[i].fSpecification;
	      fSegments.push_back(segment);
	      HLTDebug("set segment %s with size %d at offset %d", AliHLTComponent::DataType2Text(segment.fDataType).data(), segment.fSegmentSize, segment.fSegmentOffset);
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
      } else {
	HLTError("this data buffer (%p) does not match the internal data buffer %p of raw buffer %p", pTgt, fpBuffer->GetPointer(), fpBuffer);
	iResult=-EINVAL;
      }
    } else {
      HLTFatal("internal data structur missmatch");
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
  int iResult=fpBuffer==NULL || GetNofSegments()==0;
  return iResult;
}

int AliHLTDataBuffer::GetNofSegments()
{
  // see header file for function documentation
  int iResult=fSegments.size();
  return iResult;
}

int AliHLTDataBuffer::GetNofConsumers()
{
  // see header file for function documentation
  int iResult=fConsumers.size() + GetNofActiveConsumers() + fReleasedConsumers.size();
  return iResult;
}

int AliHLTDataBuffer::GetNofPendingConsumers()
{
  // see header file for function documentation
  int iResult=fConsumers.size();
  return iResult;
}

int AliHLTDataBuffer::GetNofActiveConsumers()
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
  AliHLTRawBufferPList::iterator buffer=fgFreeBuffers.begin();
  while (buffer!=fgFreeBuffers.end() && pRawBuffer==NULL) {
    if ((*buffer)->CheckSize(reqSize)) {
      // assign this element
      pRawBuffer=*buffer;
      pRawBuffer->UseBuffer(size);
      fgFreeBuffers.erase(buffer);
      fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "raw buffer container %p provided for request of %d bytes (total %d available in buffer %p)", pRawBuffer, size, pRawBuffer->GetTotalSize(), pRawBuffer->GetPointer());
      fgActiveBuffers.push_back(pRawBuffer);
      break;
    }
    buffer++;
  }
  if (pRawBuffer==NULL) {
    // no buffer found, create a new one
    pRawBuffer=new AliHLTRawBuffer(reqSize);
    if (pRawBuffer) {
      if (pRawBuffer->GetPointer()) {
	pRawBuffer->UseBuffer(size);
	fgActiveBuffers.push_back(pRawBuffer);
	fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "new raw buffer %p of size %d created (container %p)", pRawBuffer->GetPointer(), pRawBuffer->GetTotalSize(), pRawBuffer);
      } else {
	delete pRawBuffer;
	pRawBuffer=NULL;
	fgLogging.Logging(kHLTLogError, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "memory allocation failed");
      } 
    } else {
      fgLogging.Logging(kHLTLogError, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "memory allocation failed");
    }
  }
  if (pRawBuffer!=NULL && fgkSafetyPatternSize>0) {
    //fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "writing safety pattern to %p offset %d", (*buffer)->GetPointer(), (*buffer)->GetUsedSize());
    //int res=pRawBuffer->WritePattern(fgkSafetyPattern, fgkSafetyPatternSize);
    assert(res>=0);
  }
  return pRawBuffer;
}

int AliHLTDataBuffer::ReleaseRawBuffer(AliHLTRawBuffer* pBuffer)
{
  // see header file for function documentation
  int iResult=0;
  if (pBuffer) {
    AliHLTRawBufferPList::iterator buffer=fgActiveBuffers.begin();
    while (buffer!=fgActiveBuffers.end() && (*buffer)!=pBuffer) {
      buffer++;
    }
    if (buffer!=fgActiveBuffers.end()) {
      if (fgkSafetyPatternSize>0) {
	//fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "comparing safety pattern at %p offset %d", (*buffer)->GetPointer(), reinterpret_cast<AliHLTUInt32_t>(*buffer));
	if ((*buffer)->CheckPattern(fgkSafetyPattern, fgkSafetyPatternSize)) {
	  fgLogging.Logging(kHLTLogFatal, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "component has written beyond end of data buffer %p size %d", (*buffer)->GetPointer(), (*buffer)->GetUsedSize());
	}
      }
      (*buffer)->Reset();
      fgFreeBuffers.push_back(*buffer);
      fgActiveBuffers.erase(buffer);
    } else {
      fgLogging.Logging(kHLTLogWarning, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "can not find raw buffer container %p in the list of active containers", pBuffer);
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
//   int iTotalSize=0;
//   int iCount=fgFreeBuffers.size()+fgActiveBuffers.size();
  AliHLTRawBufferPList::iterator buffer;;
  while ((buffer=fgFreeBuffers.begin())!=fgFreeBuffers.end()) {
//     iTotalSize+=(*buffer)->GetTotalSize();
    delete *buffer;
    fgFreeBuffers.erase(buffer);
  }
  while ((buffer=fgActiveBuffers.begin())!=fgActiveBuffers.end()) {
//     iTotalSize+=(*buffer)->GetTotalSize();
    fgLogging.Logging(kHLTLogWarning, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "request to delete active raw buffer container (raw buffer %p, size %d)", (*buffer)->GetPointer(), (*buffer)->GetTotalSize());
    delete *buffer;
    fgActiveBuffers.erase(buffer);
  }
//   fgLogging.Logging(kHLTLogInfo, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "Total memory allocation: %d byte in %d buffers", iTotalSize, iCount);
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
    HLTWarning("consumer %p was not released", pDesc);
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

int AliHLTDataBuffer::FindConsumer(AliHLTComponent* pConsumer, int bAllLists)
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
  :
  fSize(0),
  fTotalSize(size),
  fPtr(static_cast<AliHLTUInt8_t*>(malloc(size)))
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

AliHLTDataBuffer::AliHLTRawBuffer::~AliHLTRawBuffer()
{
  if (fPtr) {
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

AliHLTUInt8_t* AliHLTDataBuffer::AliHLTRawBuffer::UseBuffer(AliHLTUInt32_t size)
{
  // see header file for function documentation
  if (size>0 && fTotalSize>=size) {
    fSize=size;
    return fPtr;
  }
  return NULL;
}

int AliHLTDataBuffer::AliHLTRawBuffer::CheckSize(AliHLTUInt32_t size) const
{
  // see header file for function documentation
  return fTotalSize>=size && ((fTotalSize-size)<fgMargin);
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
