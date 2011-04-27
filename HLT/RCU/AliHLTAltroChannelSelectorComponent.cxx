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

/// @file   AliHLTAltroChannelSelectorComponent.cxx
/// @author Matthias Richter
/// @date   
/// @brief  A filter/selective readout component for Altro data.
///

#include <cassert>
#include <memory>
#include "AliHLTAltroChannelSelectorComponent.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDAQ.h"
#include "AliRawReaderMemory.h"
#include "AliAltroRawStreamV3.h"
#include "TMath.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAltroChannelSelectorComponent)

AliHLTAltroChannelSelectorComponent::AliHLTAltroChannelSelectorComponent()
  :
  AliHLTProcessor(),
  fSkipCorrupted(true),
  fTalkative(false),
  fStartTimeBin(0),
  fEndTimeBin(0),
  fSignalThreshold(0),
  fRMSThreshold(0)
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
  char* cpErr=NULL;
  int i=0;
  for (; i<argc && iResult>=0; i++) {
    cpErr=NULL;
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -skip-corrupted, just for backward compatibility, not announced
    if (argument.CompareTo("-skip-corrupted")==0) {
      fSkipCorrupted=true;

    // -keep-corrupted
    } else if (argument.CompareTo("-keep-corrupted")==0) {
      fSkipCorrupted=false;

    // -talkative
    } else if (argument.CompareTo("-talkative")==0) {
      fTalkative=true;

    // -start-timebin
    } else if (argument.CompareTo("-start-timebin")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fStartTimeBin = strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

    // -end-timebin
    } else if (argument.CompareTo("-end-timebin")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fEndTimeBin = strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

    // -signal-threshold
    } else if (argument.CompareTo("-signal-threshold")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fSignalThreshold = strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

    // -rms-threshold
    } else if (argument.CompareTo("-rms-threshold")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fRMSThreshold = strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

    } else {
      HLTError("unknown argument %s", argument.Data());
      iResult=-EINVAL;
    }
  }

  if (cpErr && *cpErr) {
    HLTError("Cannot convert specifier '%s' for argument '%s'", argv[i], argument.Data());
    iResult=-EINVAL;
  } else if (bMissingParam) {
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

  if (!IsDataEvent()) {
    size=0;
    return 0;
  }

  // process the DLL input
  int blockno=0;
  const AliHLTComponentBlockData* pDesc=NULL;
  std::auto_ptr<AliRawReaderMemory> pRawReader(new AliRawReaderMemory);
  if (pRawReader.get()) return -ENOMEM;

  for (pDesc=GetFirstInputBlock(kAliHLTDataTypeDDLRaw); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    iResult=0;
    if (pDesc->fSize<=32) {
      continue;
    }

    // search for the active pad information
    AliHLTUInt16_t* pActiveHwAddressArray=NULL;
    int iArraySize=0;
    if (fSignalThreshold==0 && fRMSThreshold==0) {
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
      HLTWarning("no block of type %s for specification 0x%08x available, data block skipped", 
		 DataType2Text(kAliHLTDataTypeHwAddr16).c_str(), 
		 pDesc->fSpecification);
      break;
    }
    }

    static AliHLTErrorGuard required("AliHLTAltroChannelSelectorComponent", "DoEvent", "component commission required after major changes");
    (++required).Throw(1);

    pRawReader->Reset();
    int ddlid=AliHLTDAQ::DdlIDFromHLTBlockData(pDesc->fDataType.fOrigin, pDesc->fSpecification);
    if (ddlid<0) {
      HLTError("unable to extract DDL Id for data block %s 0x%08x", DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification);
      continue;
    }

    if (!pRawReader->AddBuffer((UChar_t*)pDesc->fPtr,pDesc->fSize, ddlid)) {
      ALIHLTERRORGUARD(1, "can not set up AltroDecoder for data block %s 0x%08x,"
		       " skipping data block and suppressing further messages",
		       DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification);
      continue;
    }

    std::auto_ptr<AliAltroRawStreamV3> altroRawStream(new AliAltroRawStreamV3(pRawReader.get()));

    if (!altroRawStream.get()) {
      iResult=-ENOMEM;
      break;
    }

    altroRawStream->Reset();
    if (!altroRawStream->NextDDL()) {
      ALIHLTERRORGUARD(1, "internal error, can not read data from AliRawReaderMemory");
      continue;
    }

    unsigned int rcuTrailerLength=0;
    if (iResult>=0 &&
	((rcuTrailerLength=altroRawStream->GetRCUTrailerSize())==0 ||
	 rcuTrailerLength>pDesc->fSize-cdhSize)) {
      if (rcuTrailerLength>0) {
	HLTWarning("corrupted data block: RCU trailer length exceeds buffer size");
      } else {
	HLTWarning("corrupted data block: RCU trailer of zero length");	
      }
      iResult=-EFAULT;
    }

    if (iResult<0) {
      // TODO: here the trigger has to come into play. It is up to
      // policy if a corrupted data block should be kept (original
      // DDL) or discarded. In any case, the block is not going to
      // be part of HLTOUT
      HLTWarning("skipping corrupted data block for event %lu, data block %s 0x%08x", evtData.fEventID,
		 DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification);
      iResult=0;
      continue;
    }

    int iSelected=0;
    int iTotal=0;
    int iCorrupted=0;
    AliHLTUInt32_t iOutputSize=0;
    AliHLTUInt32_t iNofAltro40=0;
    AliHLTUInt32_t iCapacity=size;

    // first add the RCU trailer
    AliHLTUInt8_t* pSrc=reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr);
    pSrc+=pDesc->fSize-rcuTrailerLength;
    if ((iResult=CopyBlockToEnd(outputPtr, iCapacity, iOutputSize, pSrc, rcuTrailerLength))>=0) {
      assert(iResult==(int)rcuTrailerLength);
      iOutputSize+=rcuTrailerLength;
    } else {
      HLTError("failed to write RCU trailer of length %d for block %d, too little space in output buffer?", rcuTrailerLength, blockno);
      iResult=-ENOSPC;
      break;
    }

    while (altroRawStream->NextChannel() && iResult>=0) {
      iTotal++;

      int hwAddress=altroRawStream->GetHWAddress();
      if (fSignalThreshold!=0) {
	// treshold by adc counts
	unsigned int sumSignals=0;
	unsigned int maxSignal=0;
	unsigned int nofSignals=0;
	while(altroRawStream->NextBunch()){
	  const UShort_t *bunchData=altroRawStream->GetSignals();
	  unsigned int time=altroRawStream->GetStartTimeBin();
	  for(Int_t i=0;i<altroRawStream->GetBunchLength();i++){
	    if(bunchData[i]>0){// disregarding 0 data.
	      if(time+i>=fStartTimeBin && time+i<=fEndTimeBin){
		sumSignals+=bunchData[i];
		if (maxSignal<bunchData[i]) maxSignal=bunchData[i];
		nofSignals++;
	      }
	    }
	  }
	}
	if (nofSignals==0 || maxSignal<=(sumSignals/nofSignals)+fSignalThreshold) {
	  continue;
	}

      } else if (fRMSThreshold!=0) {
	// treshold by adc counts
	unsigned int sumSignals=0;
	unsigned int maxSignal=0;
	unsigned int nofSignals=0;
	while(altroRawStream->NextBunch()){
	  const UShort_t *bunchData=altroRawStream->GetSignals();
	  unsigned int time=altroRawStream->GetStartTimeBin();
	  for(Int_t i=0;i<altroRawStream->GetBunchLength();i++){
	    if(bunchData[i]>0){// disregarding 0 data.
	      if(time+i>=fStartTimeBin && time+i<=fEndTimeBin){
		sumSignals+=bunchData[i]*bunchData[i];
		if (maxSignal<bunchData[i]) maxSignal=bunchData[i];
		nofSignals++;
	      }
	    }
	  }
	}
	if (nofSignals==0 || maxSignal<=TMath::Sqrt(sumSignals/nofSignals)*fRMSThreshold) {
	  continue;
	}
	
      } else {
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
      }

      // no of 10 bit words is without the fill words to fill complete 40 bit words
      // in addition, align to complete 40 bit words (the '+3')
      // also, the 5 bytes of the Altro trailer must be added to get the full size
      int channelSize=((altroRawStream->GetChannelPayloadSize()+2)/3)*4;
      if (channelSize==0) {
	if (fTalkative) HLTWarning("skipping zero length channel (hw address %d)", hwAddress);
	iCorrupted++;
	continue;
      }
      channelSize+=4;
      HLTDebug("ALTRO block hwAddress 0x%08x (%d) selected (active), size %d", hwAddress, hwAddress, channelSize);

      if (false) {
	if (channelSize == iResult) {
	  if (channelSize%4 == 0) {
	    iOutputSize+=channelSize;
	  } else {
	    if (fTalkative) HLTWarning("corrupted ALTRO channel: incomplete 40 bit word (channel hw address %d)", hwAddress);
	    iCorrupted++;
	    continue;
	  }
	} else {
	  if (fTalkative) HLTWarning("internal error: failed to copy full channel: %d out of %d bytes (hw address %d)", iResult, channelSize, hwAddress);
	  iCorrupted++;
	  continue;
	}
      } else {
	if (fTalkative) HLTError("failed to write ALTRO channel of length %d for block %d  (hw address %d)", channelSize, blockno, hwAddress);
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
    if (fTalkative) HLTImportant("data block %d (0x%08x): selected %d out of %d ALTRO channel(s), %d corrupted channels skipped", blockno, pDesc->fSpecification, iSelected, iTotal, iCorrupted);
  }

  if (iResult<0) {
    outputBlocks.clear();
  }

  // all data blocks need to be moved to the beginning of the
  // buffer because PubSub is not able to handle data blocks entirely
  // at the end of the buffer. The problem is that the component always
  // indicates to use the full size of the buffer
  if (outputBlocks.size()>0) {
    int offset=outputBlocks.back().fOffset;
    size-=offset;
    memmove(outputPtr, outputPtr+offset, size);
    for (AliHLTComponentBlockDataList::iterator block=outputBlocks.begin();
	 block!=outputBlocks.end();
	 block++) {
      block->fOffset-=offset;
    }
  }

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
