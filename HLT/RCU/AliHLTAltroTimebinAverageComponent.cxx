// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
//*                  Oystein Djuvsland
//*                  Matthias Richter                                      *
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

/** @file   AliHLTAltroTimebinAverageComponent.cxx
    @author Kalliopi Kanaki, Oystein Djuvsland, Matthias Richter
    @date   26.08.2008
    @brief  
*/

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTAltroTimebinAverageComponent.h"
#include "AliAltroDecoder.h"
#include "AliAltroData.h"
#include "AliAltroBunch.h"
#include "AliHLTAltroEncoder.h"
#include "AliRawDataHeader.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAltroTimebinAverageComponent)

AliHLTAltroTimebinAverageComponent::AliHLTAltroTimebinAverageComponent()
    :
    fStartTimeBin(0),
    fEndTimeBin(1024),
    fNTimeBins(1024)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTAltroTimebinAverageComponent::~AliHLTAltroTimebinAverageComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTAltroTimebinAverageComponent::GetComponentID()
{
  // see header file for class documentation
  return "AltroTimebinAverager";
}

void AliHLTAltroTimebinAverageComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC );
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginPHOS );
}

AliHLTComponentDataType AliHLTAltroTimebinAverageComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeDDLRaw;
}

int AliHLTAltroTimebinAverageComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeDDLRaw);
  return tgtList.size();
}

void AliHLTAltroTimebinAverageComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=0;
  inputMultiplier=1.0;
}

AliHLTComponent* AliHLTAltroTimebinAverageComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTAltroTimebinAverageComponent;
}
	
int AliHLTAltroTimebinAverageComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation

  Int_t i = 0;
  Char_t* cpErr;

  while ( i < argc ) {      

    // -- number of timebins
    if ( !strcmp( argv[i], "ntimebins" ) ) {
      fNTimeBins = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert ntimebins specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- first timebin
    if ( !strcmp( argv[i], "start-timebin" ) ) {
      fStartTimeBin = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert start-timebin specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- last timebin
    if ( !strcmp( argv[i], "end-timebin" ) ) {
      if(strtoul( argv[i+1], &cpErr ,0)<=1024){
	fEndTimeBin = strtoul( argv[i+1], &cpErr ,0);
      }
      if ( *cpErr ) {
	HLTError("Cannot convert end-timebin specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    HLTError("Unknown option '%s'", argv[i]);
    return -EINVAL;

  }

  return 0;
}

int AliHLTAltroTimebinAverageComponent::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTAltroTimebinAverageComponent::DoEvent( const AliHLTComponentEventData& evtData, 
						const AliHLTComponentBlockData* blocks, 
						AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
						AliHLTUInt32_t& size, 
						vector<AliHLTComponentBlockData>& outputBlocks )
{
  // see header file for class documentation
  int iResult=0;
  AliHLTUInt32_t capacity=size;
  size=0;
  AliHLTUInt32_t offset=0;

  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;

  AliAltroDecoder* decoder=NULL;
  AliHLTAltroEncoder* altroEncoder=NULL;

  for(ndx = 0; ndx < evtData.fBlockCnt; ndx++) {
    iter = blocks+ndx;
      
    if ( iter->fDataType != kAliHLTDataTypeDDLRaw) {
      continue;
    }
 
    if (decoder) delete decoder;
    decoder=new AliAltroDecoder;
    if (!decoder) {
      iResult=-ENOMEM;
      break;
    }

    if (altroEncoder) delete altroEncoder;

    int localResult=0;
    if ((localResult=decoder->SetMemory((UChar_t*)iter->fPtr,iter->fSize))<0) {
      HLTWarning("can not set up AltroDecoder for data block %s 0x%08x: error %d, skipping data block",
		 DataType2Text(iter->fDataType).c_str(), iter->fSpecification, localResult);
      continue;
    }
      
    if (!decoder->Decode()) {
      HLTWarning("can not decode data block %s 0x%08x: skipping data block",
		 DataType2Text(iter->fDataType).c_str(), iter->fSpecification);
      continue;
    }

    UChar_t *RCUTrailer=NULL;
    Int_t RCUTrailerSize=decoder->GetRCUTrailerSize();
    if (RCUTrailerSize<=0 || !decoder->GetRCUTrailerData(RCUTrailer) || RCUTrailer==NULL) {
      HLTWarning("can not find RCU trailer for data block %s 0x%08x: skipping data block",
		 DataType2Text(iter->fDataType).c_str(), iter->fSpecification);
      continue;
    }
      
    AliAltroData altrochannel;
    while (iResult>=0 && decoder->NextChannel(&altrochannel) && iResult>=0) {
      int hwadd=altrochannel.GetHadd();

      AliAltroBunch altrobunch;
      while (iResult>=0 && altrochannel.NextBunch(&altrobunch) && iResult>=0) {
	int bunchLength=altrobunch.GetBunchSize();
	int bunchEndTime=altrobunch.GetEndTimeBin();
	int time=bunchEndTime-bunchLength+1;
	const  UInt_t* bunchData=altrobunch.GetData();
	for (int bin=bunchLength && iResult>=0; bin>0; ) {
	  bin--;
	  if(bunchData[bin]>0){// disregarding 0 data.
	     
	    if(time+bin>=fStartTimeBin && time+bin<=fEndTimeBin){
	      if (!altroEncoder) {
		// set up the encoder
		altroEncoder=new AliHLTAltroEncoder;
		if (!altroEncoder) {
		  iResult=-ENOMEM;
		  break;
		}
		altroEncoder->SetBuffer(outputPtr+offset,capacity-offset);
		AliRawDataHeader cdh;
		altroEncoder->SetCDH((AliHLTUInt8_t*)iter->fPtr,sizeof(AliRawDataHeader));

		altroEncoder->SetRCUTrailer(RCUTrailer, RCUTrailerSize);
	      }
		
	      AliHLTUInt16_t signal=bunchData[bin];
	      if (bin-1>=0) signal+=bunchData[bin-1];
	      altroEncoder->AddSignal((time+bin)/2,signal/2);
	      bin--;
	    } // end if between start and end time bin
	  } // end if bunchData[i]>0
	} // for loop
      } //while loop over bunches
      if (altroEncoder) {
	altroEncoder->SetChannel(hwadd);
      }
    } // while loop over channels

    if (altroEncoder) {
     int sizeOfData=altroEncoder->SetLength();
     
     if (sizeOfData<0) {
       HLTError("data encoding failed");
       return sizeOfData;
     }
     if(sizeOfData>(int)capacity){
       HLTWarning("Buffer too small to add the altrodata: %d of %d byte(s) already used", sizeOfData, size);
       return -ENOSPC;
     }
   
     AliHLTComponentBlockData bd;
     FillBlockData( bd );
     bd.fOffset = offset;
     bd.fSize = sizeOfData;
     bd.fDataType = iter->fDataType;
     bd.fSpecification = iter->fSpecification;     
     outputBlocks.push_back( bd );
     
     offset+=bd.fSize;
    }

  } // while over data blocks

  if (iResult>=0) size=offset;
  return iResult;
}

