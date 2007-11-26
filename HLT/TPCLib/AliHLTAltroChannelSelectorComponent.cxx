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

/** @file   AliHLTAltroChannelSelectorComponent.cxx
    @author Matthias Richter
    @date   
    @brief  A filter/selective readout component for TPC Altro data. */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include "AliHLTAltroChannelSelectorComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDigitReaderRaw.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCPadArray.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAltroChannelSelectorComponent)

AliHLTAltroChannelSelectorComponent::AliHLTAltroChannelSelectorComponent()
  :
  AliHLTProcessor(),
  fRawreaderMode(0)
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
  //list.push_back(channel list);
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

    // -rawreadermode
    if (argument.CompareTo("-rawreadermode")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      int mode=AliHLTTPCDigitReaderRaw::DecodeMode(argv[i]);
      if (mode<0) {
	HLTError("invalid rawreadermode specifier '%s'", argv[i]);
	iResult=-EINVAL;
      } else {
	fRawreaderMode=static_cast<unsigned>(mode);
      }
    } else {
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

  // process the DLL input
  int blockno=0;
  for (; blockno<(int)evtData.fBlockCnt; blockno++ ) {
    if (blocks[blockno].fDataType != (kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC)) continue;

    // search for the active pad information
    AliHLTTPCPadArray::AliHLTTPCActivePads* pActivePadsArray=NULL;
    int iNofActivePads=0;
    for (int i=0; i<(int)evtData.fBlockCnt; i++ ) {
      if (blocks[i].fDataType == AliHLTTPCDefinitions::fgkActivePadsDataType &&
	  blocks[i].fSpecification==blocks[blockno].fSpecification) {
	pActivePadsArray=reinterpret_cast<AliHLTTPCPadArray::AliHLTTPCActivePads*>(blocks[i].fPtr);
	iNofActivePads=blocks[i].fSize/sizeof(AliHLTTPCPadArray::AliHLTTPCActivePads);
      }
    }
    if (pActivePadsArray==NULL) {
      HLTWarning("no block of type %s for specification 0x%08x available, data block unchanged", 
		 DataType2Text(AliHLTTPCDefinitions::fgkActivePadsDataType).c_str(), 
		 blocks[blockno].fSpecification);
      // forward the whole block
      outputBlocks.push_back(blocks[blockno]);
      continue;
    }

    int part=AliHLTTPCDefinitions::GetMinPatchNr(blocks[blockno]);
    assert(part==AliHLTTPCDefinitions::GetMaxPatchNr(blocks[blockno]));
    int slice=AliHLTTPCDefinitions::GetMinSliceNr(blocks[blockno]);
    assert(slice==AliHLTTPCDefinitions::GetMaxSliceNr(blocks[blockno]));
    int firstRow=AliHLTTPCTransform::GetFirstRow(part);
    int lastRow=AliHLTTPCTransform::GetLastRow(part);
    AliHLTTPCDigitReaderRaw reader(fRawreaderMode);
    reader.InitBlock(blocks[blockno].fPtr,blocks[blockno].fSize,firstRow,lastRow,part,slice);
    AliHLTUInt32_t iOutputSize=0;
    AliHLTUInt32_t iCapacity=size;
    while (reader.NextAltroBlock()) {
      int active=0;
      for (; active<iNofActivePads; active++) {
	if ((int)pActivePadsArray[active].fRow==reader.GetRow() &&
	    (int)pActivePadsArray[active].fPad==reader.GetPad()) {
	  break;
	}
      }
      if (active>=iNofActivePads) {
	HLTDebug("ALTRO block Row %d, Pad %d discarded (inactive)", reader.GetRow(), reader.GetPad());
	continue;
      }

      void* pChannel=NULL;
      AliHLTUInt16_t hwAddress=~(AliHLTUInt16_t)0;
      int channelSize=reader.GetAltroChannelRawData(pChannel, hwAddress);
      HLTDebug("ALTRO block Row/Pad %d/%d selected (active)", reader.GetRow(), reader.GetPad());
      if (channelSize>0 && pChannel!=NULL) {
	if (iOutputSize==0) {
	  // first add the RCU trailer
	  unsigned rcuTrailerLength=reader.GetRCUDataBlockLength();
	  AliHLTUInt8_t* pSrc=reinterpret_cast<AliHLTUInt8_t*>(blocks[blockno].fPtr);
	  pSrc+=blocks[blockno].fSize;
	  if ((iResult=CopyBlockToEnd(outputPtr, iCapacity, iOutputSize, pSrc, rcuTrailerLength))>=0) {
	    assert(iResult==rcuTrailerLength);
	    iOutputSize+=rcuTrailerLength;
	  } else {
	    HLTError("failed to writer RCU trailer of length %d for block %d", rcuTrailerLength, blockno);
	    break;
	  }
	}
      }
      if ((iResult=CopyBlockToEnd(outputPtr, iCapacity, iOutputSize, pChannel, channelSize))>=0) {
	assert(iResult==channelSize);
	iOutputSize+=channelSize;
      } else {
	HLTError("failed to writer ALTRO channel of length %d for block %d", channelSize, blockno);
	break;
      }
    }
    if (iResult>=0) {
      // write the common data header
      int cdhSize=reader.GetCommonDataHeaderSize();
      if ((iResult=CopyBlockToEnd(outputPtr, iCapacity, iOutputSize, blocks[blockno].fPtr, cdhSize))>=0) {
	assert(iResult==cdhSize);
	iOutputSize+=cdhSize;

	// set new length of the data block
	AliHLTUInt32_t* pCdhSize=reinterpret_cast<AliHLTUInt32_t*>(outputPtr+iCapacity-iOutputSize+1);
	*pCdhSize=iOutputSize;

	// insert block descriptor
	AliHLTComponentBlockData bd;
	FillBlockData(bd);
	bd.fOffset=iCapacity-iOutputSize;
	bd.fSize=iOutputSize;
	bd.fDataType=blocks[blockno].fDataType;
	bd.fSpecification=blocks[blockno].fSpecification;
	outputBlocks.push_back(bd);
	iCapacity-=iOutputSize;
      } else {
	HLTError("failed to write CDH of length %d for block %d", cdhSize, blockno);
	break;
      }
    }
  }

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
