// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Per Thomas Hille, Oystein Djuvsland                   *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSRawAnalyzer.h"
#include "AliHLTPHOSRawAnalyzerComponentv2.h"
#include "AliHLTPHOSChannelDataHeaderStruct.h"
#include "AliHLTPHOSChannelDataStruct.h"
#include "AliHLTPHOSMapper.h"
#include "AliHLTPHOSSanityInspector.h"
#include  "AliAltroDecoder.h"
#include  "AliAltroData.h"   
#include  "AliAltroBunch.h"  


AliHLTPHOSRawAnalyzerComponentv2::AliHLTPHOSRawAnalyzerComponentv2():
  AliHLTPHOSRcuProcessor(), 
  fAnalyzerPtr(0), 
  fMapperPtr(0), 
  fSanityInspectorPtr(0),
  fDecoderPtr(0),  
  fAltroDataPtr(0),
  fAltroBunchPtr(0),
  fAlgorithm(0),
  fOffset(0),
  fBunchSizeCut(0),
  fMinPeakPosition(0),
  fMaxPeakPosition(100)
{
  //comment
  fMapperPtr = new AliHLTPHOSMapper();
  fAltroDataPtr = new AliAltroData();
  fAltroBunchPtr = new AliAltroBunch();
  fDecoderPtr = new AliAltroDecoder();
  fSanityInspectorPtr = new AliHLTPHOSSanityInspector();
}


AliHLTPHOSRawAnalyzerComponentv2::~AliHLTPHOSRawAnalyzerComponentv2()
{
  //comment
  Deinit();
}



int 
AliHLTPHOSRawAnalyzerComponentv2::Deinit()
{
  //comment
  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
  if(fMapperPtr)
    {
      delete  fMapperPtr;
      fMapperPtr = 0;
    }
  if(fAltroDataPtr)
    {
      delete fAltroDataPtr;
      fAltroDataPtr = 0;
    }
  if(fAltroBunchPtr)
    {
      delete fAltroBunchPtr;
      fAltroBunchPtr = 0;
    }
  if(fDecoderPtr)
    {
      delete fDecoderPtr;
      fDecoderPtr = 0;
    }
  return 0;
}

const char* 
AliHLTPHOSRawAnalyzerComponentv2::GetComponentID()
{
  //comment
  return "PhosRawAnalyzerv2";
}


void
AliHLTPHOSRawAnalyzerComponentv2::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //comment
  list.clear();
  list.push_back( AliHLTPHOSDefinitions::fgkDDLPackedRawDataType | kAliHLTDataOriginPHOS);
}

AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponentv2::GetOutputDataType()
{
  //comment
  return AliHLTPHOSDefinitions::fgkChannelDataType;
}

void
AliHLTPHOSRawAnalyzerComponentv2::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  //comment
  constBase = sizeof(AliHLTPHOSChannelDataHeaderStruct);
  inputMultiplier = 0.5;
}

int 
AliHLTPHOSRawAnalyzerComponentv2::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& /*trigData*/, 
					 AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  //comment
  Int_t blockSize          = 0;
  UInt_t totSize           = 0;

  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if ( iter->fDataType != AliHLTPHOSDefinitions::fgkDDLPackedRawDataType  )
	{
	  HLTDebug("Data block is not of type fgkDDLPackedRawDataType");
	  continue; 
	}

      blockSize = DoIt(iter, outputPtr, size, totSize); // Processing the block
      if(blockSize == -1) // If the processing returns -1 we are out of buffer and return an error msg.
	{
	  return -ENOBUFS;
	}

      totSize += blockSize; //Keeping track of the used size
      // HLTDebug("Output data size: %d - Input data size: %d", totSize, iter->fSize);

      //Filling the block data
      AliHLTComponentBlockData bdChannelData;
      FillBlockData( bdChannelData );
      bdChannelData.fOffset = 0; //FIXME
      bdChannelData.fSize = blockSize;
      bdChannelData.fDataType = AliHLTPHOSDefinitions::fgkChannelDataType;
      bdChannelData.fSpecification = iter->fSpecification;
      outputBlocks.push_back(bdChannelData);

      outputPtr += blockSize; //Updating position of the output buffer
    }

  fPhosEventCount++; 
  size = totSize; //telling the framework how much buffer space we have used.
  
  return 0;
}//end DoEvent

Int_t
AliHLTPHOSRawAnalyzerComponentv2::DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr, const AliHLTUInt32_t size, UInt_t& totSize)
{
 
  //comment
  Int_t crazyness          = 0;
  Int_t nSamples           = 0;
  Short_t channelCount     = 0;

  // Firs we want to write a header to the output
  AliHLTPHOSChannelDataHeaderStruct *channelDataHeaderPtr = reinterpret_cast<AliHLTPHOSChannelDataHeaderStruct*>(outputPtr); 

  // Then comes the channel data
  AliHLTPHOSChannelDataStruct *channelDataPtr = reinterpret_cast<AliHLTPHOSChannelDataStruct*>(outputPtr+sizeof(AliHLTPHOSChannelDataHeaderStruct)); 

  //Adding to the total size of data written
  totSize += sizeof(AliHLTPHOSChannelDataHeaderStruct);
 
  //Decoding data
  fDecoderPtr->SetMemory(reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize);
  fDecoderPtr->Decode();
      
  //Looping over data
  while(fDecoderPtr->NextChannel(fAltroDataPtr) == true )
    {          
      if(fAltroDataPtr->GetDataSize() != 0 )
	{
	  //Want to get the "first in time" "bunch"
	  while(fAltroDataPtr->NextBunch(fAltroBunchPtr) == true) {}
	  
	  //Skip strangely short bunches
	  if(fAltroBunchPtr->GetBunchSize() >  fBunchSizeCut)
	    {
	      totSize += sizeof(AliHLTPHOSChannelDataStruct);
	      if(totSize > size)
		{
		  HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", totSize, size);
		  return -1;
		}
	      
	      nSamples = fAltroBunchPtr->GetBunchSize();
	      
	      //crazyness = fSanityInspectorPtr->CheckInsanity(static_cast<const UInt_t*>(fAltroBunchPtr->GetData()), static_cast<const Int_t>(fAltroBunchPtr->GetBunchSize()));
	      
	      //Evalute signal amplitude and timing
	      fAnalyzerPtr->SetData(fAltroBunchPtr->GetData(), fAltroBunchPtr->GetBunchSize());   
	      fAnalyzerPtr->Evaluate(0, fAltroBunchPtr->GetBunchSize());  
	      
	      //Checking for sane timing...
	      if(fAnalyzerPtr->GetTiming() > fMinPeakPosition && fAnalyzerPtr->GetTiming() < fMaxPeakPosition)
		{
		  // Writing to the output buffer
		  channelDataPtr->fChannelID = fMapperPtr->GetChannelID(iter->fSpecification, fAltroDataPtr->GetHadd());
		  //channelDataPtr->fChannelID = fMapperPtr->GetChannelID(1, fAltroDataPtr->GetHadd());
		  channelDataPtr->fEnergy = static_cast<Float_t>(fAnalyzerPtr->GetEnergy()) - fOffset;
		  channelDataPtr->fTime = static_cast<Float_t>(fAnalyzerPtr->GetTiming());
		  channelDataPtr->fCrazyness = static_cast<Short_t>(crazyness);
		  channelCount++;
		  channelDataPtr++; // Updating position of the free output.
		}
	    }
	}
    }

  //Writing the header
  channelDataHeaderPtr->fNChannels  = channelCount;
  channelDataHeaderPtr->fAlgorithm  = fAlgorithm;
  channelDataHeaderPtr->fInfo       = 0;
  channelDataHeaderPtr->fHasRawData = false;

  HLTDebug("Number of channels: %d", channelCount);

  //returning the size used
  return sizeof(AliHLTPHOSChannelDataStruct)*channelCount + sizeof(AliHLTPHOSChannelDataHeaderStruct);
}

int
AliHLTPHOSRawAnalyzerComponentv2::DoInit( int argc, const char** argv )
{ 

  //See base class for documentation
  fPrintInfo = kFALSE;
  int iResult=0;
  fMapperPtr = new AliHLTPHOSMapper();

  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-offset", argv[i]))
	{
	  fOffset = atoi(argv[i+1]);
	}
      if(!strcmp("-bunchsizecut", argv[i]))
	{
	  fBunchSizeCut = atoi(argv[i+1]);
	}
      if(!strcmp("-minpeakposition", argv[i]))
	{
	  fMinPeakPosition = atoi(argv[i+1]);
	}
      if(!strcmp("-maxpeakposition", argv[i]))
	{
	  fMaxPeakPosition = atoi(argv[i+1]);
	}  
    }
 
  if(fMapperPtr->GetIsInitializedMapping() == false)
    {
      Logging(kHLTLogFatal, __FILE__ , IntToChar(  __LINE__ ) , "AliHLTPHOSMapper::Could not initialise mapping from file %s, aborting", fMapperPtr->GetFilePath());
      return -4;
    }

  return iResult;
}
