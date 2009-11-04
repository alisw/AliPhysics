
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
#include "AliHLTPHOSRawAnalyzerComponentv3.h"
#include "AliHLTPHOSChannelDataHeaderStruct.h"
#include "AliHLTPHOSChannelDataStruct.h"
#include "AliHLTPHOSMapper.h"
#include "AliHLTPHOSSanityInspector.h"

#include "AliAltroRawStreamV3.h"
#include "AliCaloRawStreamV3.h"
#include "AliRawReaderMemory.h"


#include "AliHLTPHOSUtilities.h"

AliHLTPHOSRawAnalyzerComponentv3::AliHLTPHOSRawAnalyzerComponentv3():
  AliHLTPHOSRcuProcessor(), 
  fAnalyzerPtr(0), 
  fMapperPtr(0), 
  fSanityInspectorPtr(0),
  fRawReaderMemoryPtr(0),
  fAltroRawStreamPtr(0),
  fAlgorithm(0),
  fOffset(0),
  fBunchSizeCut(0),
  fMinPeakPosition(0),
  fMaxPeakPosition(100),
  fDoPushRawData(false),
  fInspectSanity(false),
  fRawDataWriter(0)
{
  //comment
  fMapperPtr = new AliHLTPHOSMapper();

  fRawReaderMemoryPtr = new AliRawReaderMemory();

  fAltroRawStreamPtr = new AliAltroRawStreamV3(fRawReaderMemoryPtr);
  //  fAltroRawStreamPtr = new AliCaloRawStreamV3(fRawReaderMemoryPtr, TString("PHOS"));
  fSanityInspectorPtr = new AliHLTPHOSSanityInspector();

  if( fDoPushRawData == true  )
    {
      
      fRawDataWriter  = new AliHLTPHOSRawAnalyzerComponentv3::RawDataWriter();

    }

}


AliHLTPHOSRawAnalyzerComponentv3::~AliHLTPHOSRawAnalyzerComponentv3()
{
  //comment
  Deinit();
}

int 
AliHLTPHOSRawAnalyzerComponentv3::Deinit()
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
  if(fRawReaderMemoryPtr)
    {
      delete fRawReaderMemoryPtr;
      fRawReaderMemoryPtr = 0;
    }
  if(fAltroRawStreamPtr)
    {
      delete fAltroRawStreamPtr;
      fAltroRawStreamPtr = 0;
    }
  return 0;
}

const char* 
AliHLTPHOSRawAnalyzerComponentv3::GetComponentID()
{
  //comment
  return "PhosRawAnalyzerv3";
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

  fRawReaderMemoryPtr->SetMemory(reinterpret_cast<UChar_t*>(iter->fPtr), static_cast<ULong_t>(iter->fSize));
  fRawReaderMemoryPtr->SetEquipmentID(fMapperPtr->GetDDLFromSpec(iter->fSpecification) + 1792);
  fRawReaderMemoryPtr->Reset();
  fRawReaderMemoryPtr->NextEvent();
  
 //  if( fDoPushRawData == true)
//     {
//       fRawDataWriter->NewEvent( );
//     }
  if(fAltroRawStreamPtr != NULL)
    {
      delete fAltroRawStreamPtr;
      fAltroRawStreamPtr=NULL;
    }
  
  //  fAltroRawStreamPtr = new AliCaloRawStreamV3(fRawReaderMemoryPtr, TString("PHOS"));
  fAltroRawStreamPtr = new AliAltroRawStreamV3(fRawReaderMemoryPtr);

  if(fAltroRawStreamPtr->NextDDL())
    {
      int cnt = 0;
      while( fAltroRawStreamPtr->NextChannel()  )
	{ 
	  // Removing TRUs
	  if(  fAltroRawStreamPtr->GetHWAddress() < 128 || ( fAltroRawStreamPtr->GetHWAddress() ^ 0x800) < 128 ) 
	    {
	      continue; 
	    }
	  else
	    {
	      cnt ++;
	      UShort_t* firstBunchPtr = 0;
	      UShort_t chId = fMapperPtr->GetChannelID(iter->fSpecification, fAltroRawStreamPtr->GetHWAddress()); 
	    
// 	      if( fDoPushRawData == true)
// 		{
// 		  fRawDataWriter->SetChannelId( chId );
// 		}
	      while( fAltroRawStreamPtr->NextBunch() == true )
		{
		  nSamples = fAltroRawStreamPtr->GetBunchLength();
		  
// 		  if( fDoPushRawData == true)
// 		    {
// 		      fRawDataWriter->WriteBunchData( fAltroRawStreamPtr->GetSignals(), nSamples,  fAltroRawStreamPtr->GetEndTimeBin()  );
// 		    }
		  firstBunchPtr = const_cast< UShort_t* >(  fAltroRawStreamPtr->GetSignals()  );
		}
	      if(firstBunchPtr)
		{	      
		  totSize += sizeof( AliHLTPHOSChannelDataStruct );
		  if(totSize > size)
		    {
		      HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", totSize, size);
		      return -1;
		    }
// 		  if(fInspectSanity)
// 		    {
// 		      crazyness = fSanityInspectorPtr->CheckAndHealInsanity(firstBunchPtr, nSamples);
// 		    }

		  fAnalyzerPtr->SetData( firstBunchPtr, nSamples);
		  fAnalyzerPtr->Evaluate(0, nSamples);  
	   
		  channelDataPtr->fChannelID =  chId;
		  channelDataPtr->fEnergy = static_cast<Float_t>(fAnalyzerPtr->GetEnergy()) - fOffset;
		
		  channelDataPtr->fTime = static_cast<Float_t>(fAnalyzerPtr->GetTiming());
		  channelDataPtr->fCrazyness = static_cast<Short_t>(crazyness);
		  channelCount++;
		  channelDataPtr++; // Updating position of the free output.
		}
	    }
// 	  if(fDoPushRawData)
// 	    {
// 	      fRawDataWriter->NewChannel();
// 	    }
	}
    }
  
  //Writing the header
  channelDataHeaderPtr->fNChannels   =  channelCount;
  channelDataHeaderPtr->fAlgorithm   = fAlgorithm;
  channelDataHeaderPtr->fInfo        = 0;

//   if( fDoPushRawData == true)
//     {
//       tmpsize += fRawDataWriter->CopyBufferToSharedMemory( (UShort_t *)channelDataPtr, size, totSize);
//     }

  // channelDataHeaderPtr->fHasRawData  = false;
  channelDataHeaderPtr->fHasRawData = fDoPushRawData;

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
AliHLTPHOSRawAnalyzerComponentv3::DoInit( int argc, const char** argv )
{ 

  //See base class for documentation
  // fPrintInfo = kFALSE;
  int iResult=0;
  fMapperPtr = new AliHLTPHOSMapper();

  for(int i = 0; i < argc; i++)
    {
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
      if(!strcmp("-pushrawdata", argv[i]))
	{
	  fDoPushRawData = true;
	}
      if(!strcmp("-inspectsanity", argv[i]))
	{
	  fInspectSanity = true;
	}
    }
 
  if(fMapperPtr->GetIsInitializedMapping() == false)
    {
      Logging(kHLTLogFatal, __FILE__ , IntToChar(  __LINE__ ) , "AliHLTPHOSMapper::Could not initialise mapping from file %s, aborting", fMapperPtr->GetFilePath());
      return -4;
    }

  return iResult;
}




AliHLTPHOSRawAnalyzerComponentv3::RawDataWriter::RawDataWriter() :  //fIsFirstChannel(true),
								    fRawDataBuffer(0),
								    fCurrentChannelSize(0),
								    //    fIsFirstChannel(true),
								    fBufferIndex(0),
								    fBufferSize( NZROWSRCU*NXCOLUMNSRCU*ALTROMAXSAMPLES*NGAINS +1000 ),
								    fCurrentChannelIdPtr(0),
								    fCurrentChannelSizePtr(0),
								    fCurrentChannelDataPtr(0),
								    fTotalSize(0)
{
  fRawDataBuffer = new UShort_t[fBufferSize];
  Init();
}


   
void  
AliHLTPHOSRawAnalyzerComponentv3::RawDataWriter::Init()
{
  fCurrentChannelIdPtr = fRawDataBuffer;
  fCurrentChannelSizePtr = fRawDataBuffer +1;
  fCurrentChannelDataPtr = fRawDataBuffer +2;
  ResetBuffer();
}

 
void
AliHLTPHOSRawAnalyzerComponentv3::RawDataWriter::NewEvent()
{
  Init();
  fTotalSize = 0;
}


void
AliHLTPHOSRawAnalyzerComponentv3::RawDataWriter::NewChannel( )
{
  *fCurrentChannelSizePtr   = fCurrentChannelSize;
  fCurrentChannelIdPtr     += fCurrentChannelSize;
  fCurrentChannelSizePtr    += fCurrentChannelSize;
  fCurrentChannelDataPtr   += fCurrentChannelSize;
  fBufferIndex = 0;
  fCurrentChannelSize = 2;
  fTotalSize += 2;
}


void 
AliHLTPHOSRawAnalyzerComponentv3::RawDataWriter::WriteBunchData(const UShort_t *bunchdata,  const int length,   const UInt_t starttimebin )
{
  fCurrentChannelDataPtr[fBufferIndex] = starttimebin;
  fCurrentChannelSize ++;
  fBufferIndex++;
  fCurrentChannelDataPtr[fBufferIndex] = length;
  fCurrentChannelSize ++;
  fBufferIndex++;

  fTotalSize +=2;

  for(int i=0; i < length; i++)
    {
      fCurrentChannelDataPtr[ fBufferIndex + i ] =  bunchdata[i];
    }

  fCurrentChannelSize += length;
  fTotalSize += length;
  fBufferIndex += length;
}


void
AliHLTPHOSRawAnalyzerComponentv3::RawDataWriter::SetChannelId( const UShort_t channeldid  )
{
  *fCurrentChannelIdPtr =  channeldid;
}


void
AliHLTPHOSRawAnalyzerComponentv3::RawDataWriter::ResetBuffer()
{
  for(int i=0; i < fBufferSize ; i++ )
    {
      fRawDataBuffer[i] = 0;
    }
}


int
AliHLTPHOSRawAnalyzerComponentv3::RawDataWriter::CopyBufferToSharedMemory(UShort_t *memPtr, const int sizetotal, const int sizeused )
{
  int sizerequested =  (sizeof(int)*fTotalSize + sizeused);

  if(  sizerequested   > sizetotal  )
    {
      return 0;
  }
  else
    {
      for(int i=0; i < fTotalSize; i++)
	{
	  memPtr[i] = fRawDataBuffer[i]; 
  	}
      return fTotalSize;
   }
}
  
