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
#include <cerrno>
//#include <string>
//#include "AliHLTSystem.h"

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
  // see header file for function documentation
  HLTFatal("copy constructor untested");
}

AliHLTDataBuffer& AliHLTDataBuffer::operator=(const AliHLTDataBuffer&)
{ 
  // see header file for function documentation
  HLTFatal("assignment operator untested");
  return *this;
}

int AliHLTDataBuffer::fgNofInstances=0;
vector<AliHLTDataBuffer::AliHLTRawBuffer*> AliHLTDataBuffer::fgFreeBuffers;
vector<AliHLTDataBuffer::AliHLTRawBuffer*> AliHLTDataBuffer::fgActiveBuffers;
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

int AliHLTDataBuffer::FindMatchingDataBlocks(const AliHLTComponent* pConsumer, vector<AliHLTComponentDataType>* tgtList)
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumer) {
    vector<AliHLTDataBuffer::AliHLTDataSegment> segments;
    if ((iResult=FindMatchingDataSegments(pConsumer, segments))>=0) {
      if (tgtList) {
	vector<AliHLTDataBuffer::AliHLTDataSegment>::iterator segment=segments.begin();
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
  iResult=tgtList.size();
  return iResult;
  
  if (pConsumer) {
    vector<AliHLTComponentDataType> dtlist;
    ((AliHLTComponent*)pConsumer)->GetInputDataTypes(dtlist);
    vector<AliHLTDataBuffer::AliHLTDataSegment>::iterator segment=fSegments.begin();
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
  // see header file for function documentation
  int iResult=0;
  if (pConsumer && arrayBlockDesc) {
    if (fpBuffer) {
      AliHLTConsumerDescriptor* pDesc=FindConsumer(pConsumer, fConsumers);
      if (pDesc) {
	vector<AliHLTDataBuffer::AliHLTDataSegment> tgtList;
	// Matthias 26.07.2007 AliHLTSystem should behave the same way as PubSub
	// so it does not matter if there are matching data types or not, unless
	// we implement such a check in PubSub
	if ((iResult=FindMatchingDataSegments(pConsumer, tgtList))>=0) {
	  int i =0;
	  vector<AliHLTDataBuffer::AliHLTDataSegment>::iterator segment=tgtList.begin();
	  while (segment!=tgtList.end() && i<iArraySize) {
	    // fill the block data descriptor
	    arrayBlockDesc[i].fStructSize=sizeof(AliHLTComponentBlockData);
	    // the shared memory key is not used in AliRoot
	    arrayBlockDesc[i].fShmKey.fStructSize=sizeof(AliHLTComponentShmData);
	    arrayBlockDesc[i].fShmKey.fShmType=gkAliHLTComponentInvalidShmType;
	    arrayBlockDesc[i].fShmKey.fShmID=gkAliHLTComponentInvalidShmID;
	    arrayBlockDesc[i].fOffset=(*segment).fSegmentOffset;
	    arrayBlockDesc[i].fPtr=*fpBuffer;
	    arrayBlockDesc[i].fSize=(*segment).fSegmentSize;
	    arrayBlockDesc[i].fDataType=(*segment).fDataType;
	    arrayBlockDesc[i].fSpecification=(*segment).fSpecification;
	    pDesc->SetActiveDataSegment(arrayBlockDesc[i].fOffset, arrayBlockDesc[i].fSize);
	    HLTDebug("component %p (%s) subscribed to segment #%d offset %d size %d data type %s %#x", 
		     pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID(), i, arrayBlockDesc[i].fOffset,
		     arrayBlockDesc[i].fSize, (AliHLTComponent::DataType2Text(arrayBlockDesc[i].fDataType)).c_str(), 
		     arrayBlockDesc[i].fSpecification);
	    i++;
	    segment++;
	  }
	  // check whether there was enough space for the segments
	  if (i!=tgtList.size()) {
	    HLTError("too little space in block descriptor array: required %d, available %d", tgtList.size(), iArraySize);
	    iResult=-ENOSPC;
	  } else {
	  // move this consumer to the active list
	  if (ChangeConsumerState(pDesc, fConsumers, fActiveConsumers)>=0) {
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

int AliHLTDataBuffer::Release(AliHLTComponentBlockData* pBlockDesc, const AliHLTComponent* pConsumer)
{
  // see header file for function documentation
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
	  // the pointer can be either NULL, than only the offset is considered, or a valid
	  // pointer in the range of the buffer
	  // The operator '>' includes the size of the buffer
	  if (arrayBlockData[i].fPtr==NULL ||
	      ((*fpBuffer)<=arrayBlockData[i].fPtr && (*fpBuffer)>arrayBlockData[i].fPtr)) {
	    int ptrOffset=0;
	    if (arrayBlockData[i].fPtr!=NULL) {
	      ptrOffset=(*fpBuffer)-arrayBlockData[i].fPtr;
	    }
	    if (arrayBlockData[i].fOffset+ptrOffset+arrayBlockData[i].fSize<=fpBuffer->fSize) {
	      segment.fSegmentOffset=arrayBlockData[i].fOffset+ptrOffset;
	      segment.fSegmentSize=arrayBlockData[i].fSize;
	      segment.fDataType=arrayBlockData[i].fDataType;
	      segment.fSpecification=arrayBlockData[i].fSpecification;
	      fSegments.push_back(segment);
	      HLTDebug("set segment %s with size %d at offset %d", AliHLTComponent::DataType2Text(segment.fDataType).data(), segment.fSegmentSize, segment.fSegmentOffset);
	    } else {
	      HLTError("block data specification %#d (%s) exceeds size of data buffer", i, AliHLTComponent::DataType2Text(arrayBlockData[i].fDataType).data());
	      HLTError("block offset=%d, block size=%d, buffer size=%d", arrayBlockData[i].fOffset, arrayBlockData[i].fSize, fpBuffer->fSize);
	      iResult=-E2BIG;
	    }
	  } else {
	    HLTError("invalid pointer (%p) in block data specification (buffer %p size %d)", arrayBlockData[i].fPtr, fpBuffer->fPtr, fpBuffer->fSize);
	    iResult=-ERANGE;
	  }
	}
      } else {
	HLTError("this data buffer (%p) does not match the internal data buffer %p of raw buffer %p", pTgt, fpBuffer->fPtr, fpBuffer);
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
  vector<AliHLTRawBuffer*>::iterator buffer=fgFreeBuffers.begin();
  while (buffer!=fgFreeBuffers.end() && pRawBuffer==NULL) {
    if ((*buffer)->fTotalSize>=reqSize && ((*buffer)->fTotalSize-reqSize)<fgMargin) {
      // assign this element
      pRawBuffer=*buffer;
      pRawBuffer->fSize=size;
      fgFreeBuffers.erase(buffer);
      fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::CreateRawBuffer", "data buffer handling", "raw buffer container %p provided for request of %d bytes (total %d available in buffer %p)", pRawBuffer, size, pRawBuffer->fTotalSize, pRawBuffer->fPtr);
      fgActiveBuffers.push_back(pRawBuffer);
      break;
    }
    buffer++;
  }
  if (pRawBuffer==NULL) {
    // no buffer found, create a new one
    pRawBuffer=new AliHLTRawBuffer;
    if (pRawBuffer) {
      pRawBuffer->fPtr=static_cast<AliHLTUInt8_t*>(malloc(reqSize));
      if (pRawBuffer->fPtr) {
	pRawBuffer->fSize=size;
	pRawBuffer->fTotalSize=reqSize;
	fgActiveBuffers.push_back(pRawBuffer);
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
  if (pRawBuffer!=NULL && fgkSafetyPatternSize>0) {
    //fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "writing safety pattern to %p offset %d", pRawBuffer->fPtr, pRawBuffer->fSize);
    memcpy(((char*)pRawBuffer->fPtr)+pRawBuffer->fSize, fgkSafetyPattern, fgkSafetyPatternSize);
  }
  return pRawBuffer;
}

int AliHLTDataBuffer::ReleaseRawBuffer(AliHLTRawBuffer* pBuffer)
{
  // see header file for function documentation
  int iResult=0;
  if (pBuffer) {
    vector<AliHLTRawBuffer*>::iterator buffer=fgActiveBuffers.begin();
    while (buffer!=fgActiveBuffers.end() && (*buffer)!=pBuffer) {
      buffer++;
    }
    if (buffer!=fgActiveBuffers.end()) {
      if (fgkSafetyPatternSize>0) {
	//fgLogging.Logging(kHLTLogDebug, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "comparing safety pattern at %p offset %d", (*buffer)->fPtr, (*buffer)->fSize);
	if (memcmp(((char*)(*buffer)->fPtr)+(*buffer)->fSize, fgkSafetyPattern, fgkSafetyPatternSize)!=0) {
	  fgLogging.Logging(kHLTLogFatal, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "component has written beyond end of data buffer %p size %d", (*buffer)->fPtr, (*buffer)->fSize);
	}
      }
      (*buffer)->fSize=0;
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
  vector<AliHLTRawBuffer*>::iterator buffer=fgFreeBuffers.begin();
  while (buffer!=fgFreeBuffers.end()) {
    free((*buffer)->fPtr);
    delete *buffer;
    fgFreeBuffers.erase(buffer);
    buffer=fgFreeBuffers.begin();
  }
  buffer=fgActiveBuffers.begin();
  while (buffer!=fgActiveBuffers.end()) {
    fgLogging.Logging(kHLTLogWarning, "AliHLTDataBuffer::ReleaseRawBuffer", "data buffer handling", "request to delete active raw buffer container (raw buffer %p, size %d)", (*buffer)->fPtr, (*buffer)->fTotalSize);
    free((*buffer)->fPtr);
    delete *buffer;
    fgActiveBuffers.erase(buffer);
    buffer=fgActiveBuffers.begin();
  }
  return iResult;
}

AliHLTConsumerDescriptor* AliHLTDataBuffer::FindConsumer(const AliHLTComponent* pConsumer, vector<AliHLTConsumerDescriptor*> &list) const
{
  // see header file for function documentation
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
  // see header file for function documentation
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
  vector<AliHLTDataBuffer::AliHLTDataSegment>::iterator segment=fSegments.begin();
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
  // see header file for function documentation
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

int AliHLTDataBuffer::CleanupConsumerList() 
{
  // see header file for function documentation
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

int AliHLTDataBuffer::FindConsumer(AliHLTComponent* pConsumer, int bAllLists)
{
  // see header file for function documentation
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

int AliHLTDataBuffer::AliHLTRawBuffer::operator==(void* ptr)
{
  return fPtr == static_cast<AliHLTUInt8_t*>(ptr);
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator<=(void* ptr)
{
  int iResult=fPtr <= static_cast<AliHLTUInt8_t*>(ptr);
  //printf("%p: %p <= %p (%d)\n", this, fPtr, ptr, iResult);
  return iResult;
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator>(void* ptr)
{
  int iResult=fPtr+fSize > static_cast<AliHLTUInt8_t*>(ptr);
  //printf("%p: %p + %d > %p (%d)\n", this, fPtr, fSize, ptr, iResult);
  return iResult;
}

int AliHLTDataBuffer::AliHLTRawBuffer::operator-(void* ptr)
{
  return static_cast<int>(static_cast<AliHLTUInt8_t*>(ptr)-fPtr);
}
