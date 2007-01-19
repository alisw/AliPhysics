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

/** @file   AliHLTDataBuffer.cxx
    @author Matthias Richter
    @date   
    @brief  Handling of Data Buffers for HLT components.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTDataBuffer.h"
#include "AliHLTComponent.h"
#include <string>
#include "AliHLTSystem.h"

AliHLTConsumerDescriptor::AliHLTConsumerDescriptor()
  :
  fpConsumer(NULL),
  fSegments()
{
  fSegments.clear();
}

AliHLTConsumerDescriptor::AliHLTConsumerDescriptor(AliHLTComponent* pConsumer)
  :
  fpConsumer(pConsumer),
  fSegments()
{
  fSegments.clear();
}

AliHLTConsumerDescriptor::AliHLTConsumerDescriptor(const AliHLTConsumerDescriptor& desc)
  :
  TObject(),
  AliHLTLogging(),
  fpConsumer(desc.fpConsumer),
  fSegments()
{
  // we can simply transfer the pointer to th new object since there are no
  // release actions in the destructor
}

AliHLTConsumerDescriptor& AliHLTConsumerDescriptor::operator=(const AliHLTConsumerDescriptor& desc)
{ 
  // we can simply transfer the pointer to th new object since there are no
  // release actions in the destructor
  fpConsumer=desc.fpConsumer;
  return *this;
}

AliHLTConsumerDescriptor::~AliHLTConsumerDescriptor()
{
  if (fSegments.size()>0) {
    //HLTWarning("unreleased data segments found");
  }
}

int AliHLTConsumerDescriptor::SetActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size)
{
  int iResult=0;
  AliHLTDataSegment segment(offset, size);
  fSegments.push_back(segment);
  //HLTDebug("set active segment (%d:%d) for consumer %p", offset, size, this);
  return iResult;
}

int AliHLTConsumerDescriptor::CheckActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size)
{
  int iResult=0;
  if (fSegments.size()>0) {
    vector<AliHLTDataSegment>::iterator segment=fSegments.begin();
    while (segment!=fSegments.end()) {
      if ((iResult=((*segment).fSegmentOffset==offset && (*segment).fSegmentSize==size))>0) {
	break;
      }
      segment++;
    }
  } else {
    //HLTWarning("no data segment active for consumer %p", this);
    iResult=-ENODATA;
  }
  return iResult;
}

int AliHLTConsumerDescriptor::ReleaseActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size)
{
  int iResult=0;
  if (fSegments.size()>0) {
    vector<AliHLTDataSegment>::iterator segment=fSegments.begin();
    while (segment!=fSegments.end()) {
      if ((iResult=((*segment).fSegmentOffset==offset && (*segment).fSegmentSize==size))>0) {
	fSegments.erase(segment);
	break;
      }
      segment++;
    }
    if (iResult==0) {
      //HLTWarning("no data segment (%d:%d) active for consumer %p", offset, size, this);
      iResult=-ENOENT;
    }
  } else {
    //HLTWarning("no data segment active for consumer %p", this);
    iResult=-ENODATA;
  }
  return iResult;
}

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataBuffer)

AliHLTDataBuffer::AliHLTDataBuffer()
  :
  fSegments(),
  fConsumers(),
  fActiveConsumers(),
  fReleasedConsumers(),
  fpBuffer(NULL),
  fFlags(0)
{
  fSegments.empty();
  fConsumers.empty();
  fActiveConsumers.empty();
  fReleasedConsumers.empty();
  fNofInstances++;
}

AliHLTDataBuffer::AliHLTDataBuffer(const AliHLTDataBuffer&)
  :
  TObject(),
  AliHLTLogging(),
  fSegments(),
  fConsumers(),
  fActiveConsumers(),
  fReleasedConsumers(),
  fpBuffer(NULL),
  fFlags(0)
{
  HLTFatal("copy constructor untested");
}

AliHLTDataBuffer& AliHLTDataBuffer::operator=(const AliHLTDataBuffer&)
{ 
  HLTFatal("assignment operator untested");
  return *this;
}

int AliHLTDataBuffer::fNofInstances=0;
vector<AliHLTRawBuffer*> AliHLTDataBuffer::fFreeBuffers;
vector<AliHLTRawBuffer*> AliHLTDataBuffer::fActiveBuffers;
AliHLTUInt32_t AliHLTDataBuffer::fMargin=1024;
AliHLTLogging AliHLTDataBuffer::fgLogging;

AliHLTDataBuffer::~AliHLTDataBuffer()
{
  if (--fNofInstances<=0) {
    DeleteRawBuffers();
  }
  CleanupConsumerList();
}

int AliHLTDataBuffer::SetConsumer(AliHLTComponent* pConsumer)
{
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

int AliHLTDataBuffer::FindMatchingDataBlocks(const AliHLTComponent* pConsumer, vector<AliHLTComponentDataType>* tgtList)
{
  int iResult=0;
  if (pConsumer) {
    vector<AliHLTDataSegment> segments;
    if ((iResult=FindMatchingDataSegments(pConsumer, segments))>=0) {
      if (tgtList) {
	vector<AliHLTDataSegment>::iterator segment=segments.begin();
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

int AliHLTDataBuffer::FindMatchingDataSegments(const AliHLTComponent* pConsumer, vector<AliHLTDataSegment>& tgtList)
{
  int iResult=0;
  if (pConsumer) {
    vector<AliHLTComponentDataType> dtlist;
    ((AliHLTComponent*)pConsumer)->GetInputDataTypes(dtlist);
    vector<AliHLTDataSegment>::iterator segment=fSegments.begin();
    while (segment!=fSegments.end()) {
      vector<AliHLTComponentDataType>::iterator type=dtlist.begin();
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
  int iResult=0;
  if (pConsumer && arrayBlockDesc) {
    if (fpBuffer) {
      AliHLTConsumerDescriptor* pDesc=FindConsumer(pConsumer, fConsumers);
      if (pDesc) {
	vector<AliHLTDataSegment> tgtList;
	/* TODO: think about a good policy for this check
	 * is it enough that at least one segment is available, or have all to be available?
	 * or is it possible to have optional segments?
	 */
	if ((iResult=FindMatchingDataSegments(pConsumer, tgtList))>0) {
	  int i =0;
	  vector<AliHLTDataSegment>::iterator segment=tgtList.begin();
	  while (segment!=tgtList.end() && i<iArraySize) {
	    // fill the block data descriptor
	    arrayBlockDesc[i].fStructSize=sizeof(AliHLTComponentBlockData);
	    // the shared memory key is not used in AliRoot
	    arrayBlockDesc[i].fShmKey.fStructSize=sizeof(AliHLTComponentShmData);
	    arrayBlockDesc[i].fShmKey.fShmType=gkAliHLTComponentInvalidShmType;
	    arrayBlockDesc[i].fShmKey.fShmID=gkAliHLTComponentInvalidShmID;
	    arrayBlockDesc[i].fOffset=(*segment).fSegmentOffset;
	    arrayBlockDesc[i].fPtr=fpBuffer->fPtr;
	    arrayBlockDesc[i].fSize=(*segment).fSegmentSize;
	    arrayBlockDesc[i].fDataType=(*segment).fDataType;
	    arrayBlockDesc[i].fSpecification=(*segment).fSpecification;
	    pDesc->SetActiveDataSegment(arrayBlockDesc[i].fOffset, arrayBlockDesc[i].fSize);
	    HLTDebug("component %p (%s) subscribed to segment #%d offset %d", pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID(), i, arrayBlockDesc[i].fOffset);
	    i++;
	    segment++;
	  }
	  // move this consumer to the active list
	  if (ChangeConsumerState(pDesc, fConsumers, fActiveConsumers)>=0) {
	    HLTDebug("component %p (%s) subscribed to data buffer %p", pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID(), this);
	  } else {
	    // TODO: cleanup the consumer descriptor correctly
	    memset(arrayBlockDesc, 0, iArraySize*sizeof(AliHLTComponentBlockData));
	    HLTError("can not activate consumer %p for data buffer %p", pConsumer, this);
	    iResult=-EACCES;
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
      HLTError("data buffer %p is empty", this);
      iResult=-ENODATA;
    }
  } else {
    HLTError("invalid parameter");
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::Release(AliHLTComponentBlockData* pBlockDesc, const AliHLTComponent* pConsumer)
{
  int iResult=0;
  if (pBlockDesc && pConsumer) {
    AliHLTConsumerDescriptor* pDesc=FindConsumer(pConsumer, fActiveConsumers);
    if (pDesc) {
      if ((iResult=pDesc->CheckActiveDataSegment(pBlockDesc->fOffset, pBlockDesc->fSize))!=1) {
	HLTWarning("data segment missmatch, component %p has not subscribed to a segment with offset %#x and size %d", pConsumer, pBlockDesc->fOffset, pBlockDesc->fSize);
	// TODO: appropriate error handling, but so far optional
	iResult=0;
      } else {
	pDesc->ReleaseActiveDataSegment(pBlockDesc->fOffset, pBlockDesc->fSize);
	pBlockDesc->fOffset=0;
	pBlockDesc->fPtr=NULL;
	pBlockDesc->fSize=0;
      }
      if (pDesc->GetNofActiveSegments()==0) {
	if ((iResult=ChangeConsumerState(pDesc, fActiveConsumers, fReleasedConsumers))>=0) {
	  if (GetNofActiveConsumers()==0) {
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

AliHLTUInt8_t* AliHLTDataBuffer::GetTargetBuffer(int iMinSize)
{
  AliHLTUInt8_t* pTargetBuffer=NULL;
  fpBuffer=CreateRawBuffer(iMinSize);
  if (fpBuffer) {
    pTargetBuffer=(AliHLTUInt8_t*)fpBuffer->fPtr;
  } else {
    HLTError("can not create raw buffer");
  }
  return pTargetBuffer;
}

int AliHLTDataBuffer::SetSegments(AliHLTUInt8_t* pTgt, AliHLTComponentBlockData* arrayBlockData, int iSize)
{
  int iResult=0;
  if (pTgt && arrayBlockData && iSize>=0) {
    if (fpBuffer) {
      if (fpBuffer->fPtr==(void*)pTgt) {
	AliHLTDataSegment segment;
	memset(&segment, 0, sizeof(AliHLTDataSegment));
	for (int i=0; i<iSize; i++) {
	  if (arrayBlockData[i].fOffset+arrayBlockData[i].fSize<=fpBuffer->fSize) {
	    segment.fSegmentOffset=arrayBlockData[i].fOffset;
	    segment.fSegmentSize=arrayBlockData[i].fSize;
	    segment.fDataType=arrayBlockData[i].fDataType;
	    segment.fSpecification=arrayBlockData[i].fSpecification;
	    fSegments.push_back(segment);
	    HLTDebug("set segment %s with size %d at offset %d", AliHLTComponent::DataType2Text(segment.fDataType).data(), segment.fSegmentSize, segment.fSegmentOffset);
	  } else {
	    HLTError("block data specification %#d (%s) exceeds size of data buffer", i, AliHLTComponent::DataType2Text(arrayBlockData[i].fDataType).data());
	    HLTError("block offset=%d, block size=%d, buffer size=%d", arrayBlockData[i].fOffset, arrayBlockData[i].fSize, fpBuffer->fSize);
	  }
	}
      } else {
	HLTError("this data buffer (%p) does not match the internal data buffer %p of raw buffer %p", pTgt, fpBuffer->fPtr, fpBuffer);
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
  int iResult=fpBuffer==NULL || GetNofSegments()==0;
  return iResult;
}

int AliHLTDataBuffer::GetNofSegments()
{
  int iResult=fSegments.size();
  return iResult;
}

int AliHLTDataBuffer::GetNofConsumers()
{
  int iResult=fConsumers.size() + GetNofActiveConsumers() + fReleasedConsumers.size();
  return iResult;
}

int AliHLTDataBuffer::GetNofActiveConsumers()
{
  int iResult=fActiveConsumers.size();
  return iResult;
}

AliHLTRawBuffer* AliHLTDataBuffer::CreateRawBuffer(AliHLTUInt32_t size)
{
  AliHLTRawBuffer* pRawBuffer=NULL;
  vector<AliHLTRawBuffer*>::iterator buffer=fFreeBuffers.begin();
  while (buffer!=fFreeBuffers.end() && pRawBuffer==NULL) {
    if ((*buffer)->fTotalSize>=size && ((*buffer)->fTotalSize-size)<fMargin) {
      // assign this element
      pRawBuffer=*buffer;
      pRawBuffer->fSize=size;
      fFreeBuffers.erase(buffer);
      fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "raw buffer container %p provided for request of %d bytes (total %d available in buffer %p)", pRawBuffer, size, pRawBuffer->fTotalSize, pRawBuffer->fPtr);
      fActiveBuffers.push_back(pRawBuffer);
      break;
    }
    buffer++;
  }
  if (pRawBuffer==NULL) {
    // no buffer found, create a new one
    pRawBuffer=new AliHLTRawBuffer;
    if (pRawBuffer) {
      memset(pRawBuffer, 0, sizeof(AliHLTRawBuffer));
      pRawBuffer->fPtr=malloc(size);
      if (pRawBuffer->fPtr) {
	pRawBuffer->fSize=size;
	pRawBuffer->fTotalSize=size;
	fActiveBuffers.push_back(pRawBuffer);
	fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "new raw buffer %p of size %d created (container %p)", pRawBuffer->fPtr, pRawBuffer->fTotalSize, pRawBuffer);
      } else {
	delete pRawBuffer;
	pRawBuffer=NULL;
	fgLogging.Logging(kHLTLogError, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "memory allocation failed");
      } 
    } else {
      fgLogging.Logging(kHLTLogError, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "memory allocation failed");
    }
  }
  return pRawBuffer;
}

int AliHLTDataBuffer::ReleaseRawBuffer(AliHLTRawBuffer* pBuffer)
{
  int iResult=0;
  if (pBuffer) {
    vector<AliHLTRawBuffer*>::iterator buffer=fActiveBuffers.begin();
    while (buffer!=fActiveBuffers.end() && (*buffer)!=pBuffer) {
      buffer++;
    }
    if (buffer!=fActiveBuffers.end()) {
      (*buffer)->fSize=0;
      fFreeBuffers.push_back(*buffer);
      fActiveBuffers.erase(buffer);
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
  int iResult=0;
  vector<AliHLTRawBuffer*>::iterator buffer=fFreeBuffers.begin();
  while (buffer!=fFreeBuffers.end()) {
    free((*buffer)->fPtr);
    delete *buffer;
    fFreeBuffers.erase(buffer);
    buffer=fFreeBuffers.begin();
  }
  buffer=fActiveBuffers.begin();
  while (buffer!=fActiveBuffers.end()) {
    fgLogging.Logging(kHLTLogWarning, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "request to delete active raw buffer container (raw buffer %p, size %d)", (*buffer)->fPtr, (*buffer)->fTotalSize);
    free((*buffer)->fPtr);
    delete *buffer;
    fActiveBuffers.erase(buffer);
    buffer=fActiveBuffers.begin();
  }
  return iResult;
}

AliHLTConsumerDescriptor* AliHLTDataBuffer::FindConsumer(const AliHLTComponent* pConsumer, vector<AliHLTConsumerDescriptor*> &list)
{
  AliHLTConsumerDescriptor* pDesc=NULL;
  vector<AliHLTConsumerDescriptor*>::iterator desc=list.begin();
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
  int iResult=0;
  AliHLTRawBuffer* pBuffer=fpBuffer;
  fpBuffer=NULL;

  // cleanup consumer states
  vector<AliHLTConsumerDescriptor*>::iterator desc=fReleasedConsumers.begin();
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
  vector<AliHLTDataSegment>::iterator segment=fSegments.begin();
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
  return ResetDataBuffer();
}

// this is the version which works on lists of components instead of consumer descriptors
// int AliHLTDataBuffer::ChangeConsumerState(AliHLTComponent* pConsumer, vector<AliHLTComponent*> &srcList, vector<AliHLTComponent*> &tgtList)
// {
//   int iResult=0;
//   if (pDesc) {
//     vector<AliHLTComponent*>::iterator desc=srcList.begin();
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

int AliHLTDataBuffer::ChangeConsumerState(AliHLTConsumerDescriptor* pDesc, vector<AliHLTConsumerDescriptor*> &srcList, vector<AliHLTConsumerDescriptor*> &tgtList)
{
  int iResult=-ENOENT;
  if (pDesc) {
    vector<AliHLTConsumerDescriptor*>::iterator desc=srcList.begin();
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

int AliHLTDataBuffer::CleanupConsumerList() {
  int iResult=0;
  ResetDataBuffer();
  vector<AliHLTConsumerDescriptor*>::iterator desc=fConsumers.begin();
  while (desc!=fConsumers.end()) {
    delete *desc;
    fConsumers.erase(desc);
    desc=fConsumers.begin();
  }
  return iResult;
}

int AliHLTDataBuffer::FindConsumer(AliHLTComponent* pConsumer, int bAllLists) {
  vector<AliHLTConsumerDescriptor*>::iterator desc=fConsumers.begin();
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
