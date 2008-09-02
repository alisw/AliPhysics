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

/** @file   AliHLTAltroChannelSelectorComponent.cxx
    @author Matthias Richter
    @date   
    @brief  A filter/selective readout component for Altro data.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include "AliHLTAltroChannelSelectorComponent.h"
#include "AliAltroDecoder.h"
#include "AliAltroData.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAltroChannelSelectorComponent)

AliHLTAltroChannelSelectorComponent::AliHLTAltroChannelSelectorComponent()
  :
  AliHLTProcessor(),
  fSkipCorrupted(false),
  fTalkative(false)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTAltroChannelSelectorComponent::~AliHLTAltroChannelSelectorComponent()
{
  // see header file for class documentation
}

const char* AliHLTAltroChannelSelectorComponent::GetComponentID()
{
  // see header file for class documentation
  return "AltroChannelSelector";
}

void AliHLTAltroChannelSelectorComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC);
  list.push_back(kAliHLTDataTypeHwAddr16);
}

AliHLTComponentDataType AliHLTAltroChannelSelectorComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC;
}

void AliHLTAltroChannelSelectorComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase=0;
  inputMultiplier=1.0;
}

AliHLTComponent* AliHLTAltroChannelSelectorComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTAltroChannelSelectorComponent;
}

int AliHLTAltroChannelSelectorComponent::DoInit(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  bool bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -skip-corrupted
    if (argument.CompareTo("-skip-corrupted")==0) {
      fSkipCorrupted=true;

    // -talkative
    } else if (argument.CompareTo("-talkative")==0) {
      fTalkative=true;
    } else {
      HLTError("unknown argument %s", argument.Data());
      iResult=-EINVAL;
    }
  }

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  return iResult;
}

int AliHLTAltroChannelSelectorComponent::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTAltroChannelSelectorComponent::DoEvent(const AliHLTComponentEventData& evtData,
						 const AliHLTComponentBlockData* blocks, 
						 AliHLTComponentTriggerData& /*trigData*/,
						 AliHLTUInt8_t* outputPtr, 
						 AliHLTUInt32_t& size,
						 AliHLTComponentBlockDataList& outputBlocks )
{
  // see header file for class documentation
  int iResult=0;
  const int cdhSize=32;

  // process the DLL input
  int blockno=0;
  const AliHLTComponentBlockData* pDesc=NULL;

  AliAltroDecoder* decoder=NULL;
  for (pDesc=GetFirstInputBlock(kAliHLTDataTypeDDLRaw); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    iResult=0;
    if (pDesc->fSize<=32) continue;

    // search for the active pad information
    AliHLTUInt16_t* pActiveHwAddressArray=NULL;
    int iArraySize=0;
    for (int i=0; i<(int)evtData.fBlockCnt; i++ ) {
      // search for selection data of hw address type
      // which matches the data specification of the block
      if (blocks[i].fDataType==kAliHLTDataTypeHwAddr16 && blocks[i].fSpecification==pDesc->fSpecification) {
	pActiveHwAddressArray=reinterpret_cast<AliHLTUInt16_t*>(blocks[i].fPtr);
	iArraySize=blocks[i].fSize/sizeof(AliHLTUInt16_t);
	break;
      }
    }
    if (pActiveHwAddressArray==NULL) {
      HLTWarning("no block of type %s for specification 0x%08x available, data block unchanged", 
		 DataType2Text(kAliHLTDataTypeHwAddr16).c_str(), 
		 pDesc->fSpecification);
      iResult=-EFAULT;
    }

    if (decoder) delete decoder;
    decoder=new AliAltroDecoder;
    if (decoder->SetMemory(reinterpret_cast<UChar_t*>(pDesc->fPtr), pDesc->fSize)<0) {
      HLTWarning("corrupted data block: initialization of decoder failed for block: %s specification %#x size %d",
		 DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification, pDesc->fSize);
      iResult=-EFAULT;
    } else {
      if (decoder->Decode()) {
	HLTDebug("init decoder %p size %d", pDesc->fPtr,pDesc->fSize);
      } else {
	HLTWarning("corrupted data block: decoding failed for raw data block: %s specification %#x size %d",
		   DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification, pDesc->fSize);
	iResult=-EFAULT;
      }
    }

    int rcuTrailerLength=decoder->GetRCUTrailerSize();
    if (rcuTrailerLength>pDesc->fSize-cdhSize) {
      HLTWarning("corrupted data block: RCU trailer length exceeds buffer size");
      iResult=-EFAULT;
    }

    if (iResult<0) {
      // forward the whole block
      outputBlocks.push_back(*pDesc);
      iResult=0;
      continue;
    }

    int iSelected=0;
    int iTotal=0;
    int iCorrupted=0;
    AliHLTUInt32_t iOutputSize=0;
    AliHLTUInt32_t iNofAltro40=0;
    AliHLTUInt32_t iCapacity=size;
    AliAltroData channel;
    while (decoder->NextChannel(&channel) && iResult>=0) {
      iTotal++;

      int hwAddress=channel.GetHadd();
      int active=0;
      for (active=0; active<iArraySize; active++) {
	if (pActiveHwAddressArray[active]==(AliHLTUInt16_t)hwAddress) {
	  break;
	}
      }
      if (active>=iArraySize) {
	HLTDebug("ALTRO block %#x (%d) discarded (inactive)", hwAddress, hwAddress);
	continue;
      }

      // no of 10 bit words is without the fill words to fill complete 40 bit words
      // in addition, align to complete 40 bit words (the '+3')
      // also, the 5 bytes of the Altro trailer must be added to get the full size
      int channelSize=((channel.GetDataSize()+3)/4)*5;
      channelSize+=5;
      HLTDebug("ALTRO block hwAddress 0x%08x (%d) selected (active), size %d", hwAddress, hwAddress, channelSize);
      if (iOutputSize==0) {
	// first add the RCU trailer
	AliHLTUInt8_t* pSrc=reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr);
	pSrc+=pDesc->fSize-rcuTrailerLength;
	if ((iResult=CopyBlockToEnd(outputPtr, iCapacity, iOutputSize, pSrc, rcuTrailerLength))>=0) {
	  assert(iResult==rcuTrailerLength);
	  iOutputSize+=rcuTrailerLength;
	} else {
	  HLTError("failed to write RCU trailer of length %d for block %d, too little space in output buffer?", rcuTrailerLength, blockno);
	  iResult=-ENOSPC;
	  break;
	}
      }

      if ((iResult=decoder->CopyBackward(outputPtr, iCapacity-iOutputSize))>=0) {
	if (channelSize == iResult) {
	  if (channelSize%5 == 0) {
	    iNofAltro40+=channelSize/5;
	    iOutputSize+=channelSize;
	  } else {
	    if (fTalkative) HLTWarning("corrupted ALTRO channel: incomplete 40 bit word");
	    iCorrupted++;
	    continue;
	  }
	} else {
	  if (fTalkative) HLTWarning("internal error: failed to copy full channel: %d out of %d bytes", iResult, channelSize);
	  iCorrupted++;
	  continue;
	}
      } else {
	if (fTalkative) HLTError("failed to write ALTRO channel of length %d for block %d", channelSize, blockno);
	// corrupted channel, but keep going
	iCorrupted++;
	iResult=0;
	continue;
      }
      iSelected++;
    }
    if (iResult>=0) {
      // write the common data header
      if ((iResult=CopyBlockToEnd(outputPtr, iCapacity, iOutputSize, pDesc->fPtr, cdhSize))>=0) {
	assert(iResult==cdhSize);
	iOutputSize+=cdhSize;

	// set new length of the data block
	AliHLTUInt32_t* pCdhSize=reinterpret_cast<AliHLTUInt32_t*>(outputPtr+iCapacity-iOutputSize);
	*pCdhSize=iOutputSize;

	// set number of Altro words
	AliHLTUInt32_t* pNofAltro40=reinterpret_cast<AliHLTUInt32_t*>(outputPtr+iCapacity-rcuTrailerLength);
	*pNofAltro40=iNofAltro40;

	// insert block descriptor
	AliHLTComponentBlockData bd;
	FillBlockData(bd);
	bd.fOffset=iCapacity-iOutputSize;
	bd.fSize=iOutputSize;
	bd.fDataType=pDesc->fDataType;
	bd.fSpecification=pDesc->fSpecification;
	outputBlocks.push_back(bd);
	iCapacity-=iOutputSize;
      } else {
	HLTError("failed to write CDH of length %d for block %d", cdhSize, blockno);
	break;
      }
    }
    HLTInfo("data block %d (0x%08x): selected %d out of %d ALTRO channel(s), %d corrupted channels skipped", blockno, pDesc->fSpecification, iSelected, iTotal, iCorrupted);
  }
  if (decoder) delete decoder;

  if (iResult<0) {
    outputBlocks.clear();
  }

  // !!! do not change the size since the output buffer is filled from the end !!!

  return iResult;
}

int AliHLTAltroChannelSelectorComponent::CopyBlockToEnd(AliHLTUInt8_t* pTgt, unsigned capacity, unsigned position, void* pSrc, unsigned size)
{
  int iResult=0;
  if (pTgt==NULL || pSrc==NULL) return -EINVAL;
  if (capacity-position<size) return -ENOSPC;
  
  memcpy(pTgt+(capacity-position-size), pSrc, size);
  iResult=size;
  
  return iResult;
}
