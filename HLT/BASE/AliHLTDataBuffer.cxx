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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// handling of HLT data buffers                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTDataBuffer.h"
#include <string>
#include "AliHLTSystem.h"

AliHLTConsumerDescriptor::AliHLTConsumerDescriptor()
{
  fpConsumer=NULL;
  memset(&fDataType, 0, sizeof(AliHLTComponent_DataType));
  fDataType.fStructSize=sizeof(AliHLTComponent_DataType);
  fpSegment=NULL;
}

AliHLTConsumerDescriptor::AliHLTConsumerDescriptor(AliHLTComponent* pConsumer, AliHLTComponent_DataType datatype)
{
  fpConsumer=pConsumer;
  fDataType=datatype;
  fpSegment=NULL;
}

AliHLTConsumerDescriptor::~AliHLTConsumerDescriptor()
{
}

int AliHLTConsumerDescriptor::SetActiveDataSegment(AliHLTDataSegment* pSegment)
{
  int iResult=0;
  return iResult;
}

int AliHLTConsumerDescriptor::CheckActiveDataSegment(AliHLTUInt32_t offset, AliHLTUInt32_t size)
{
  int iResult=0;
  return iResult;
}

int AliHLTConsumerDescriptor::ReleaseActiveDataSegment()
{
  int iResult=0;
  return iResult;
}

ClassImp(AliHLTDataBuffer)

int AliHLTDataBuffer::fNofInstances=0;
vector<AliHLTRawBuffer*> AliHLTDataBuffer::fFreeBuffers;
vector<AliHLTRawBuffer*> AliHLTDataBuffer::fActiveBuffers;
AliHLTUInt32_t AliHLTDataBuffer::fMargin=1024;


AliHLTDataBuffer::AliHLTDataBuffer()
{
  // TODO: do the right initialization 
  //fSegments.empty();
  //fConsumers;
  //fActiveConsumers;
  //fReleasedConsumers;
  fpBuffer=NULL;
  fFlags=0;
  fNofInstances++;
}

AliHLTDataBuffer::~AliHLTDataBuffer()
{
  if (--fNofInstances<=0) {
    DeleteRawBuffers();
  }
  CleanupConsumerList();
}

int AliHLTDataBuffer::SetConsumer(AliHLTComponent* pConsumer, AliHLTComponent_DataType datatype)
{
  int iResult=0;
  if (pConsumer) {
    AliHLTConsumerDescriptor* pDesc=new AliHLTConsumerDescriptor(pConsumer, datatype);
    if (pDesc) {
      fConsumers.push_back(pDesc);
    } else {
      HLTError("memory allocation failed");
      iResult=-ENOMEM;
    }
  } else {
    HLTError("invalid parameter");
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTDataBuffer::Subscribe(AliHLTComponent_DataType datatype, const AliHLTComponent* pConsumer, AliHLTComponent_BlockData* pBlockDesc)
{
  int iResult=0;
  if (pConsumer && pBlockDesc) {
    if (fpBuffer) {
      AliHLTConsumerDescriptor* pDesc=FindConsumer(pConsumer, datatype, fConsumers);
      if (pDesc) {
	AliHLTDataSegment* pSegment=FindDataSegment(datatype);
	if (pSegment) {
	  // move this consumer to the active list
	  if ((iResult=ChangeConsumerState(pDesc, fConsumers, fActiveConsumers))>=0) {
	    pDesc->SetActiveDataSegment(pSegment);
	    // fill the block data descriptor
	    pBlockDesc->fStructSize=sizeof(AliHLTComponent_BlockData);
	    // the shared memory key is not used in AliRoot
	    pBlockDesc->fShmKey.fStructSize=sizeof(AliHLTComponent_ShmData);
	    pBlockDesc->fShmKey.fShmType=0;
	    pBlockDesc->fShmKey.fShmID=0;
	    pBlockDesc->fOffset=pSegment->fSegmentOffset;
	    pBlockDesc->fPtr=fpBuffer->fPtr;
	    pBlockDesc->fSize=pSegment->fSegmentSize;
	    pBlockDesc->fDataType=pSegment->fDataType;
	    pBlockDesc->fSpecification=pSegment->fSpecification;
	    HLTDebug("component %p (%s) subscribed to data buffer %p (%s)", pConsumer, ((AliHLTComponent*)pConsumer)->GetComponentID(), this, datatype.fID);
	  } else {
	    HLTError("can not activate consumer %p for data buffer %p", pConsumer, this);
	    iResult=-EACCES;
	  }
	} else {
	  HLTError("unresolved data segment: %s::%s is not available", datatype.fID, datatype.fOrigin);
	  iResult=-EBADF;
	}
      } else {
	HLTWarning("can not find consumer %p in component list of data buffer %d", pConsumer, this);
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

int AliHLTDataBuffer::Release(AliHLTComponent_BlockData* pBlockDesc, const AliHLTComponent* pConsumer)
{
  int iResult=0;
  if (pBlockDesc && pConsumer) {
      AliHLTConsumerDescriptor* pDesc=FindConsumer(pConsumer, pBlockDesc->fDataType, fActiveConsumers);
      if (pDesc) {
	if ((iResult=pDesc->CheckActiveDataSegment(pBlockDesc->fOffset, pBlockDesc->fSize))!=1) {
	  HLTWarning("data segment missmatch, component %p has not subscribed to a segment with offset %#x and size %d", pBlockDesc->fOffset, pBlockDesc->fSize);
	  // TODO: appropriate error handling, but so far optional
	  iResult=0;
	}
	pDesc->ReleaseActiveDataSegment();
	pBlockDesc->fOffset=0;
	pBlockDesc->fPtr=NULL;
	pBlockDesc->fSize=0;
	if ((iResult=ChangeConsumerState(pDesc, fActiveConsumers, fReleasedConsumers))>=0) {
	  if (GetNofActiveConsumers()==0) {
	    // this is the last consumer, release the raw buffer
	    AliHLTRawBuffer* pBuffer=fpBuffer;
	    ResetDataBuffer();
	    ReleaseRawBuffer(pBuffer);
	  }
	} else {
	  HLTError("can not deactivate consumer %p for data buffer %p", pConsumer, this);
	  iResult=-EACCES;
	}
      } else {
	HLTWarning("component %p has currently not subscribed to the data buffer", pConsumer);
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
  pTargetBuffer=(AliHLTUInt8_t*)fpBuffer;
  return pTargetBuffer;
}

int AliHLTDataBuffer::SetSegments(AliHLTUInt8_t* pTgt, AliHLTComponent_BlockData* arrayBlockData, int iSize)
{
  int iResult=0;
  if (pTgt && arrayBlockData && iSize>=0) {
    AliHLTDataSegment segment;
    memset(&segment, 0, sizeof(AliHLTDataSegment));
    for (int i=0; i<iSize; i++) {
      if (arrayBlockData[i].fOffset+arrayBlockData[i].fSize<fpBuffer->fSize) {
      segment.fSegmentOffset=arrayBlockData[i].fOffset;
      segment.fSegmentSize=arrayBlockData[i].fSize;
      segment.fDataType=arrayBlockData[i].fDataType;
      segment.fSpecification=arrayBlockData[i].fSpecification;
      fSegments.push_back(segment);
      } else {
	HLTError("block data specification #%d (%s@%s) exceeds size of data buffer", i, arrayBlockData[i].fDataType.fOrigin, arrayBlockData[i].fDataType.fID);
      }
    }
  } else {
    HLTError("inavalid parameter: pTgtBuffer=%p arrayBlockData=%p", pTgt, arrayBlockData);
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
      //HLTDebug("raw buffer container %p provided for request of %d bytes (total %d available in buffer %p)", pRawBuffer, size, pRawBuffer->fTotalSize, pRawBuffer->fPtr);
      fActiveBuffers.push_back(pRawBuffer);
    } else {
      // check the next element
      buffer++;
    }
  }
  if (pRawBuffer==NULL) {
    pRawBuffer=new AliHLTRawBuffer;
    if (pRawBuffer) {
      memset(pRawBuffer, 0, sizeof(AliHLTRawBuffer));
      pRawBuffer->fPtr=malloc(size);
      if (pRawBuffer->fPtr) {
	pRawBuffer->fSize=size;
	pRawBuffer->fTotalSize=size;
	fActiveBuffers.push_back(pRawBuffer);
	//HLTDebug("new raw buffer %p of size %d created (container %p)", pRawBuffer->fPtr, pRawBuffer->fTotalSize, pRawBuffer);
      } else {
	delete pRawBuffer;
	pRawBuffer=NULL;
	//HLTError("memory allocation failed");
      } 
    } else {
      //HLTError("memory allocation failed");
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
      //HLTWarning("can not find raw buffer container %p in the list of active containers", pBuffer);
      iResult=-ENOENT;
    }
  } else {
    //HLTError("invalid parameter");
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
  while (buffer!=fFreeBuffers.end()) {
    //HLTWarning("request to delete active raw buffer container %d (raw buffer %p, size %d)", *buffer, *buffer->fPtr, *buffer->fTotalSize);
    free((*buffer)->fPtr);
    delete *buffer;
    fActiveBuffers.erase(buffer);
    buffer=fActiveBuffers.begin();
  }
  return iResult;
}

AliHLTConsumerDescriptor* AliHLTDataBuffer::FindConsumer(const AliHLTComponent* pConsumer, AliHLTComponent_DataType datatype, vector<AliHLTConsumerDescriptor*> &list)
{
  AliHLTConsumerDescriptor* pDesc=NULL;
  vector<AliHLTConsumerDescriptor*>::iterator desc=list.begin();
  while (desc!=list.end() && pDesc==NULL) {
    if ((pConsumer==NULL || (*desc)->GetComponent()==pConsumer) && (*desc)->GetDataType()==datatype) {
      pDesc=*desc;
    }
  }
  return pDesc;
}

int AliHLTDataBuffer::ResetDataBuffer() {
  int iResult=0;
  fpBuffer=NULL;
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
  return iResult;
}

int AliHLTDataBuffer::ChangeConsumerState(AliHLTConsumerDescriptor* pDesc, vector<AliHLTConsumerDescriptor*> &srcList, vector<AliHLTConsumerDescriptor*> &tgtList)
{
  int iResult=0;
  if (pDesc) {
    vector<AliHLTConsumerDescriptor*>::iterator desc=srcList.begin();
    while (desc!=srcList.end()) {
      if ((*desc)==pDesc) {
	srcList.erase(desc);
	tgtList.push_back(pDesc);
	break;
      }
    }
    if (desc==srcList.end()) {
      HLTError("can not find consumer descriptor %p in list", pDesc);
      iResult=-ENOENT;
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

AliHLTDataSegment* AliHLTDataBuffer::FindDataSegment(AliHLTComponent_DataType datatype)
{
  AliHLTDataSegment* pSegment=NULL;
  vector<AliHLTDataSegment>::iterator segment=fSegments.begin();
  while (segment!=fSegments.end() && pSegment==NULL) {
    if ((*segment).fDataType==datatype) {
      // TODO: check this use of the vector
      //pSegment=segment;
    }
  }
  return pSegment;
}
