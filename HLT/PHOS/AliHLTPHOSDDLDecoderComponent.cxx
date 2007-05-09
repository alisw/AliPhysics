/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        * 
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                * 
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTPHOSDDLDecoderComponent.h"
#include <iostream>
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include "AliHLTPHOSRcuChannelDataStruct.h"
#include "AliHLTPHOSPulseGenerator.h"
#include "AliHLTPHOSDataCorruptor.h"


using namespace std;


AliHLTPHOSDDLDecoderComponent  gAliHLTPHOSDDLDecoderComponent;


AliHLTPHOSDDLDecoderComponent::AliHLTPHOSDDLDecoderComponent():AliHLTPHOSProcessor(), fPHOSRawStream(0), fRawMemoryReader(0), 
fOutPtr(0), fDataCorruptorPtr(0)
{
  //Default constructor
  fDataCorruptorPtr = new AliHLTPHOSDataCorruptor();
} 

AliHLTPHOSDDLDecoderComponent::~AliHLTPHOSDDLDecoderComponent()
{
  //Destructor
  delete  fRawMemoryReader;
  delete  fPHOSRawStream;
  delete  fDataCorruptorPtr;
}


AliHLTPHOSDDLDecoderComponent::AliHLTPHOSDDLDecoderComponent(const AliHLTPHOSDDLDecoderComponent & ) : AliHLTPHOSProcessor(), 
fPHOSRawStream(0),fRawMemoryReader(0), fOutPtr(0), fDataCorruptorPtr(0)
{

}


int 
AliHLTPHOSDDLDecoderComponent::Deinit()
{
  //Se html documentation of base class
  cout <<  "Deinit" << endl;
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRawAnalyzerComponen DoDeinit");

  if(fRawMemoryReader !=0)
    {
      delete fRawMemoryReader;
    }
    
  if(fPHOSRawStream != 0)
    {
      delete fPHOSRawStream;
    }
  return 0;

  return 0;
}


const char* 
AliHLTPHOSDDLDecoderComponent::GetComponentID()
{
  //Se html documentation of base class
  return "PhosDDLDecoder";
}



void
AliHLTPHOSDDLDecoderComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{
  //Se html documentation of base class
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSDDLDecoderComponent::GetOutputDataType()
{
  //See html documentation of base class
  return AliHLTPHOSDefinitions::fgkCellEnergyDataType;
}

void
AliHLTPHOSDDLDecoderComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //Se html documentation of base class
  constBase = 30;
  inputMultiplier = 1;
}


int
AliHLTPHOSDDLDecoderComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  //See html documentation of base class
  double testPulse[70];
  fDataCorruptorPtr->MakeCorruptedDataTest(testPulse,70);

  //  AliHLTUInt8_t tmpMod    = 0;
  //  AliHLTUInt8_t tmpZ      = 0;
  //  AliHLTUInt8_t tmpX      = 0;
  // AliHLTUInt8_t tmpGain   = 0;
  Int_t sampleCnt         = 0;
  Int_t processedChannels = 0;
  UInt_t offset           = 0; 
  UInt_t mysize           = 0;
  UInt_t tSize            = 0;

  Int_t tmpChannelCnt     = 0;
  //  Int_t tmpStartIndex     = 0;

  AliHLTUInt8_t* outBPtr;
  //unsigned long first;
  // unsigned long last;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      mysize = 0;
      tmpChannelCnt = 0;
      offset = tSize;

      if ( iter->fDataType != AliHLTPHOSDefinitions::fgkDDLPackedRawDataType )
	{
	  continue;
	}

      fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize );
      fOutPtr =  (AliHLTPHOSRcuChannelDataStruct*)outBPtr;
      mysize += sizeof(AliHLTPHOSRcuChannelDataStruct);

      fOutPtr->fRcuX = fRcuX;
      fOutPtr->fRcuZ = fRcuZ;
      fOutPtr->fModuleID = fModuleID;

      while(fPHOSRawStream->Next())
	{
	  
	  if (fPHOSRawStream->IsNewHWAddress())
	    {
	      sampleCnt = 0;
	      fOutPtr->fValidData[tmpChannelCnt].fZ = (AliHLTUInt8_t)fPHOSRawStream->GetColumn() - fRcuZOffset;
	      fOutPtr->fValidData[tmpChannelCnt].fX = (AliHLTUInt8_t)fPHOSRawStream->GetRow() - fRcuXOffset;
	      fOutPtr->fValidData[tmpChannelCnt].fGain = fPHOSRawStream->IsLowGain();
	      fOutPtr->fValidData[tmpChannelCnt].fNSamples = 0;
 	      tmpChannelCnt++;
	    }
	  fOutPtr->fValidData[tmpChannelCnt-1].fNSamples ++;
	  fOutPtr->fValidData[tmpChannelCnt-1].fChannelData[sampleCnt] = fPHOSRawStream->GetSignal();
	  sampleCnt ++; 
	}

      fOutPtr->fNValidChannels = tmpChannelCnt-1;;

      int tmpSampleCnt=0;
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;
      bd.fDataType = AliHLTPHOSDefinitions::fgkCellChannelDataDataType;
      bd.fSpecification = 0xeFFFFFFF;
      outputBlocks.push_back( bd);
      tSize += mysize;
      outBPtr += mysize;
    }

  
  if( tSize > size )
    {
      cout <<"kHLTLogFatal, HLT::AliHLTPHOSDDLDecoderComponent::DoEvent Too much data Data written over allowed buffer. Amount written:" << tSize << " allowed" << size << endl;
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSDDLDecoderComponent::DoEvent", "Too much data", "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",  tSize, size );
      return EMSGSIZE;
    }

  fPhosEventCount++; 
  
  if(fPrintInfo == kTRUE)
    {
      if(fPhosEventCount%fPrintInfoFrequncy == 0)
      	{
	  cout <<"Analyzing event " <<  fPhosEventCount  << "for Equippment " << fkEquippmentID << endl; 
	}  
    }
  size = tSize;
  return 0;
}//end DoEvent


int
AliHLTPHOSDDLDecoderComponent::DoInit( int argc, const char** argv )
{
  //See html documentation of base class 
  fPrintInfo = kFALSE;
  fRawMemoryReader = new AliRawReaderMemory();
  fPHOSRawStream = new  AliCaloRawStream(fRawMemoryReader,"PHOS");
  fPHOSRawStream->SetOldRCUFormat(kFALSE);
  int iResult=0;
  TString argument="";
  ScanArguments(argc, argv);

  if(fIsSetEquippmentID == kFALSE)

    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuHistogramProducerComponent::DoInt( int argc, const char** argv )", "Missing argument",
	       "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>");
      iResult = -2; 
    }

  return iResult;
}

AliHLTComponent*
AliHLTPHOSDDLDecoderComponent::Spawn()
{
  //See html documentation of base class 
  return new AliHLTPHOSDDLDecoderComponent;
}
