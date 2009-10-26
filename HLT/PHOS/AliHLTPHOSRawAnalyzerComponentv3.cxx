
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

#include "AliHLTPHOSRawAnalyzerComponentv3.h"
#include "AliHLTPHOSMapper.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSUtilities.h"

AliHLTPHOSRawAnalyzerComponentv3::AliHLTPHOSRawAnalyzerComponentv3():
  AliHLTCaloRawAnalyzerComponentv3()
{
  //comment
}


AliHLTPHOSRawAnalyzerComponentv3::~AliHLTPHOSRawAnalyzerComponentv3()
{
  //comment
}

int 
AliHLTPHOSRawAnalyzerComponentv3::Deinit()
{
  //comment
  return 0;
}

void
AliHLTPHOSRawAnalyzerComponentv3::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //comment
  list.clear();
  list.push_back( AliHLTPHOSDefinitions::fgkDDLPackedRawDataType | kAliHLTDataOriginPHOS );
}

AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponentv3::GetOutputDataType()
{
  //comment
  return AliHLTPHOSDefinitions::fgkChannelDataType;
}

void
AliHLTPHOSRawAnalyzerComponentv3::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  //comment
  constBase = sizeof(AliHLTPHOSChannelDataHeaderStruct);
  inputMultiplier = 0.5;
}

int 
AliHLTPHOSRawAnalyzerComponentv3::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& /*trigData*/, 
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
AliHLTPHOSRawAnalyzerComponentv3::DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr, const AliHLTUInt32_t size, UInt_t& totSize)
{
  //comment
  int tmpsize=  0;
  Int_t crazyness          = 0;
  Int_t nSamples           = 0;
  Short_t channelCount     = 0;

  // Firs we want to write a header to the output
  AliHLTPHOSChannelDataHeaderStruct *channelDataHeaderPtr = reinterpret_cast<AliHLTPHOSChannelDataHeaderStruct*>(outputPtr); 

  // Then comes the channel data
  AliHLTPHOSChannelDataStruct *channelDataPtr = reinterpret_cast<AliHLTPHOSChannelDataStruct*>(outputPtr+sizeof(AliHLTPHOSChannelDataHeaderStruct)); 
  //Adding to the total size of data written
  totSize += sizeof( AliHLTPHOSChannelDataHeaderStruct );

  fRawReaderMemoryPtr->SetMemory(         reinterpret_cast<UChar_t*>( iter->fPtr ),  static_cast<ULong_t>( iter->fSize )  );
  fRawReaderMemoryPtr->SetEquipmentID(    fMapperPtr->GetDDLFromSpec(  iter->fSpecification) + 1792  );
  fRawReaderMemoryPtr->Reset();
  fRawReaderMemoryPtr->NextEvent();
  
  if( fkDoPushRawData == true)
    {
      fRawDataWriter->NewEvent( );
    }
  if(fAltroRawStreamPtr != NULL)
    {
      delete fAltroRawStreamPtr;
      fAltroRawStreamPtr=NULL;
    }
  
  fAltroRawStreamPtr = new AliCaloRawStreamV3(fRawReaderMemoryPtr, TString("PHOS"));

  if(fAltroRawStreamPtr->NextDDL())
    {
      int cnt = 0;
      while( fAltroRawStreamPtr->NextChannel()  )
	{ 
	  if(  fAltroRawStreamPtr->GetHWAddress() < 128 || ( fAltroRawStreamPtr->GetHWAddress() ^ 0x800) < 128 ) 
	    {
	      continue; 
	    }
	  else
	    {
	      cnt ++;
	      UShort_t* firstBunchPtr = 0;
	      UShort_t chId = fMapperPtr->GetChannelID(iter->fSpecification, fAltroRawStreamPtr->GetHWAddress()); 
	    
	      if( fkDoPushRawData == true)
		{
		  fRawDataWriter->SetChannelId( chId );
		}
	      while( fAltroRawStreamPtr->NextBunch() == true )
		{
		  nSamples = fAltroRawStreamPtr->GetBunchLength();
		  
		  if( fkDoPushRawData == true)
		    {
		      //		      fRawDataWriter->WriteBunchData( fAltroRawStreamPtr->GetSignals(), nSamples,  fAltroRawStreamPtr->GetStartTimeBin()  );
		      fRawDataWriter->WriteBunchData( fAltroRawStreamPtr->GetSignals(), nSamples,  fAltroRawStreamPtr->GetEndTimeBin()  );
		    }
		  firstBunchPtr = const_cast< UShort_t* >(  fAltroRawStreamPtr->GetSignals()  );
		}
	      
	      totSize += sizeof( AliHLTPHOSChannelDataStruct );
	      if(totSize > size)
		{
		  HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", totSize, size);
		  return -1;
		}

	      fAnalyzerPtr->SetData( firstBunchPtr, nSamples);
	      fAnalyzerPtr->Evaluate(0, nSamples);  
	   
	      //	      if(fAnalyzerPtr->GetTiming() > fMinPeakPosition && fAnalyzerPtr->GetTiming() < fMaxPeakPosition)
	      {
		channelDataPtr->fChannelID =  chId;
		channelDataPtr->fEnergy = static_cast<Float_t>(fAnalyzerPtr->GetEnergy()) - fOffset;
		
		
		channelDataPtr->fTime = static_cast<Float_t>(fAnalyzerPtr->GetTiming());
		channelDataPtr->fCrazyness = static_cast<Short_t>(crazyness);
		channelCount++;
		channelDataPtr++; // Updating position of the free output.
	      }   
	    }
	  fRawDataWriter->NewChannel();
	}
    }
  
  //Writing the header
  channelDataHeaderPtr->fNChannels   =  channelCount;
  channelDataHeaderPtr->fAlgorithm   = fAlgorithm;
  channelDataHeaderPtr->fInfo        = 0;

  if( fkDoPushRawData == true)
    {
      tmpsize += fRawDataWriter->CopyBufferToSharedMemory( (UShort_t *)channelDataPtr, size, totSize);
    }

  // channelDataHeaderPtr->fHasRawData  = false;
  channelDataHeaderPtr->fHasRawData = fkDoPushRawData;

  HLTDebug("Number of channels: %d", channelCount);
  //returning the size used
  // delete fAltroRawStreamPtr;
  tmpsize += sizeof(AliHLTPHOSChannelDataStruct)*channelCount + sizeof(AliHLTPHOSChannelDataHeaderStruct); 

  //  return sizeof(AliHLTPHOSChannelDataStruct)*channelCount + sizeof(AliHLTPHOSChannelDataHeaderStruct);
  return  tmpsize;

}
 

int 
AliHLTPHOSRawAnalyzerComponentv3::WriteRawData(AliHLTPHOSChannelDataStruct*) //dtaPtr)
{
  return 0;
}

int
AliHLTPHOSRawAnalyzerComponentv3::InitMapping( const int spec)
{ 

  //See base class for documentation
  // fPrintInfo = kFALSE;

  if(fMapperPtr == 0)
    {
      fMapperPtr = new AliHLTPHOSMapper();
    }
 
  if(fMapperPtr->GetIsInitializedMapping() == false)
    {
      HLTError("%d:%d, ERROR, mapping not initialized ", __FILE__, __LINE__ );
      exit(-2);
    }

  return iResult;
}


