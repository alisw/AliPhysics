
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

#include "AliHLTCaloRawAnalyzer.h"
#include "AliHLTCaloRawAnalyzerComponentv3.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTCaloChannelDataStruct.h"
#include "AliHLTCaloMapper.h"
#include "AliHLTCaloSanityInspector.h"
#include "AliAltroRawStreamV3.h"
#include "AliRawReaderMemory.h"
#include "AliAltroRawStreamV3.h"
#include "AliCaloRawStreamV3.h"

AliHLTCaloRawAnalyzerComponentv3::AliHLTCaloRawAnalyzerComponentv3():
  AliHLTProcessor(),
  fAnalyzerPtr(0),
  fMapperPtr(0),     
  fkDoPushRawData(true),
  fPhosEventCount(0),
  fSanityInspectorPtr(0),
  fRawReaderMemoryPtr(0),
  fAltroRawStreamPtr(0),
  fAlgorithm(0),  
  fOffset(0),
  fBunchSizeCut(0),
  fMinPeakPosition(0),
  fMaxPeakPosition(100),
  fRawDataWriter(0) 
{
  //comment

  //  fMapperPtr = new AliHLTCaloMapper();

  fRawReaderMemoryPtr = new AliRawReaderMemory();

  fAltroRawStreamPtr = new AliAltroRawStreamV3(fRawReaderMemoryPtr);

  fSanityInspectorPtr = new AliHLTCaloSanityInspector();
 
  if( fkDoPushRawData == true  )
    {
      fRawDataWriter  = new RawDataWriter(); 
    }
}


AliHLTCaloRawAnalyzerComponentv3::~AliHLTCaloRawAnalyzerComponentv3()
{
  //comment
  Deinit();
}



int 
AliHLTCaloRawAnalyzerComponentv3::Deinit()
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

/*
const char* 
AliHLTCaloRawAnalyzerComponentv3::GetComponentID()
{
  //comment
  return "CaloRawAnalyzerv3";
}
*/


 /*
void
AliHLTCaloRawAnalyzerComponentv3::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //comment
  list.clear();
  list.push_back( AliHLTCaloDefinitions::fgkDDLPackedRawDataType | kAliHLTDataOriginPHOS);
}
 */

AliHLTComponentDataType 
AliHLTCaloRawAnalyzerComponentv3::GetOutputDataType()
{
  //comment
  return AliHLTCaloDefinitions::fgkChannelDataType;
}

void
AliHLTCaloRawAnalyzerComponentv3::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  //comment
  constBase = sizeof(AliHLTCaloChannelDataHeaderStruct);
  inputMultiplier = 0.5;
}

int 
AliHLTCaloRawAnalyzerComponentv3::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& /*trigData*/, 
					 AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  Int_t blockSize          = 0;
  UInt_t totSize           = 0;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if( CheckInputDataType(  iter->fDataType ) == false )
	{
	  continue;
	}
      else
	{
	  InitMapping( iter->fSpecification); 
	  blockSize = DoIt(iter, outputPtr, size, totSize); // Processing the block
	  if(blockSize == -1) // If the processing returns -1 we are out of buffer and return an error msg.
	    {
	      return -ENOBUFS;
	    }
	  
	  totSize += blockSize; //Keeping track of the used size
	  AliHLTComponentBlockData bdChannelData;
	  FillBlockData( bdChannelData );
	  bdChannelData.fOffset = 0; //FIXME
	  bdChannelData.fSize = blockSize;
	  
	  //	  bdChannelData.fDataType = AliHLTPHOSDefinitions::fgkChannelDataType;
	  bdChannelData.fDataType = AliHLTCaloDefinitions::fgkChannelDataType;

	  bdChannelData.fSpecification = iter->fSpecification;
	  outputBlocks.push_back(bdChannelData);
	  outputPtr += blockSize; //Updating position of the output buffer
	}

      fPhosEventCount++; 
      size = totSize; //telling the framework how much buffer space we have used.
    }
  return 0;
  
}//end DoEvent



Int_t
AliHLTCaloRawAnalyzerComponentv3::DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr, const AliHLTUInt32_t size, UInt_t& totSize)
{
  //comment
  int tmpsize=  0;
  Int_t crazyness          = 0;
  Int_t nSamples           = 0;
  Short_t channelCount     = 0;


  //  Firs we want to write a header to the output
  //  cout << __FILE__ << __LINE__ << "The size of AliHLTCaloChannelDataHeaderStruct is " << sizeof(AliHLTCaloChannelDataHeaderStruct ) << endl;
  //  cout << __FILE__ << __LINE__ << "The size of AliHLTCaloChannelDataStruct is " << sizeof(AliHLTCaloChannelDataStruct ) << endl;


  AliHLTCaloChannelDataHeaderStruct *channelDataHeaderPtr = reinterpret_cast<AliHLTCaloChannelDataHeaderStruct*>(outputPtr); 
  AliHLTCaloChannelDataStruct *channelDataPtr = reinterpret_cast<AliHLTCaloChannelDataStruct*>(outputPtr+sizeof(AliHLTCaloChannelDataHeaderStruct)); 
  //Adding to the total size of data written
  totSize += sizeof( AliHLTCaloChannelDataHeaderStruct );
  fRawReaderMemoryPtr->SetMemory(         reinterpret_cast<UChar_t*>( iter->fPtr ),  static_cast<ULong_t>( iter->fSize )  );
  fRawReaderMemoryPtr->SetEquipmentID(    fMapperPtr->GetDDLFromSpec(  iter->fSpecification) + 4608  );
  fRawReaderMemoryPtr->Reset();
  fRawReaderMemoryPtr->NextEvent();
  //  cout << __FILE__  << ":" << __LINE__ << ":" <<__FUNCTION__ << "T1"  << endl; 
  if( fkDoPushRawData == true)
    {
     fRawDataWriter->NewEvent( );
    }
  if(fAltroRawStreamPtr != NULL)
    {
      delete fAltroRawStreamPtr;
      fAltroRawStreamPtr=NULL;
    }


  // fAltroRawStreamPtr = new AliCaloRawStreamV3(fRawReaderMemoryPtr, TString("PHOS"));

  fAltroRawStreamPtr = new AliCaloRawStreamV3(fRawReaderMemoryPtr, TString("EMCAL"));
  //  while( fAltroRawStreamPtr->NextDDL() )
  // cout << __FILE__  << ":" << __LINE__ << ":" <<__FUNCTION__ << "T2"  << endl; 
  fAltroRawStreamPtr->NextDDL();
  {
    int cnt = 0;
    while( fAltroRawStreamPtr->NextChannel()  )
      { 
	//	 cout << __FILE__  << ":" << __LINE__ << ":" <<__FUNCTION__ << "T3"  << endl; 
	if(  fAltroRawStreamPtr->GetHWAddress() < 128 || ( fAltroRawStreamPtr->GetHWAddress() ^ 0x800) < 128 ) 
	  {
	    continue; 
	  }
	else
	  {
	    ++ cnt;
	    UShort_t* firstBunchPtr = 0;

	    //   UShort_t chId = fMapperPtr->GetChannelID(iter->fSpecification, fAltroRawStreamPtr->GetHWAddress()); 
	    //   UInt_t chId = fMapperPtr->GetChannelID(iter->fSpecification, fAltroRawStreamPtr->GetHWAddress()); 
	    
	    int chId = fMapperPtr->GetChannelID(iter->fSpecification, fAltroRawStreamPtr->GetHWAddress()); 

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
	    //	    cout << __FILE__  << ":" << __LINE__ << ":" <<__FUNCTION__ << "T5"  << endl;  
	    totSize += sizeof( AliHLTCaloChannelDataStruct );
	    if(totSize > size)
	      {
		HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", totSize, size);
	
		return -1;
	      }

	    fAnalyzerPtr->SetData( firstBunchPtr, nSamples);
	    fAnalyzerPtr->Evaluate(0, nSamples);  
	    
	    //	      if(fAnalyzerPtr->GetTiming() > fMinPeakPosition && fAnalyzerPtr->GetTiming() < fMaxPeakPosition)
	    {
	      //	      cout << __FILE__  << ":" << __LINE__ << ":" <<__FUNCTION__ << "T7"  << endl; 
	      channelDataPtr->fChannelID =  chId;
	      channelDataPtr->fEnergy = static_cast<Float_t>(fAnalyzerPtr->GetEnergy()) - fOffset;
	      //	      cout << __FILE__  << ":" << __LINE__ << ":" <<__FUNCTION__ << "T8"  << endl;
	      
	      if( channelDataPtr->fEnergy > 70 )
		{
		  //		  cout << __FILE__ << __LINE__ << "The energy is of  channel  "<< chId << "  is "  <<  channelDataPtr->fEnergy << endl;
		}
		
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
  tmpsize += sizeof(AliHLTCaloChannelDataStruct)*channelCount + sizeof(AliHLTCaloChannelDataHeaderStruct); 

  //  return sizeof(AliHLTPHOSChannelDataStruct)*channelCount + sizeof(AliHLTPHOSChannelDataHeaderStruct);
  return  tmpsize;

}

int
AliHLTCaloRawAnalyzerComponentv3::DoInit( int argc, const char** argv )
{ 

  //See base class for documentation
  //  fPrintInfo = kFALSE;

  int iResult=0;
  
  //  fMapperPtr = new AliHLTCaloMapper();
  
  //InitMapping(); 

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
 
  /*
  if( fMapperPtr->GetIsInitializedMapping() == false)
    {
      Logging(kHLTLogFatal, __FILE__ , IntToChar(  __LINE__ ) , "AliHLTCaloMapper::Could not initialise mapping from file %s, aborting", fMapperPtr->GetFilePath());
      return -4;
    }
  */

  return iResult;
}




AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::RawDataWriter() :  //fIsFirstChannel(true),
								    fRawDataBuffer(0),
								    fCurrentChannelSize(0),
								    //    fIsFirstChannel(true),
								    fBufferIndex(0),
								    //	    fBufferSize( NZROWSRCU*NXCOLUMNSRCU*ALTROMAXSAMPLES*NGAINS +1000 ),
								    fBufferSize( 64*56*ALTROMAXSAMPLES*NGAINS +1000 ),
								    fCurrentChannelIdPtr(0),
								    fCurrentChannelSizePtr(0),
								    fCurrentChannelDataPtr(0),
								    fTotalSize(0)
{
  fRawDataBuffer = new UShort_t[fBufferSize];
  Init();
}


   
void  
AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::Init()
{
  fCurrentChannelIdPtr = fRawDataBuffer;
  fCurrentChannelSizePtr = fRawDataBuffer +1;
  fCurrentChannelDataPtr = fRawDataBuffer +2;
  ResetBuffer();
}

 
void
AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::NewEvent()
{
  Init();
  fTotalSize = 0;
}


void
AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::NewChannel( )
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
AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::WriteBunchData(const UShort_t *bunchdata,  const int length,   const UInt_t starttimebin )
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
AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::SetChannelId( const UShort_t channeldid  )
{
  *fCurrentChannelIdPtr =  channeldid;
}


void
AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::ResetBuffer()
{
  for(int i=0; i < fBufferSize ; i++ )
    {
      fRawDataBuffer[i] = 0;
    }
}


int
AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::CopyBufferToSharedMemory(UShort_t *memPtr, const int sizetotal, const int sizeused )
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
  
