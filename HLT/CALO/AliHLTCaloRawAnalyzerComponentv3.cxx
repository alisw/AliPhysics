
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

//#include "AliHLTCaloRawAnalyzer.h"
#include "AliCaloRawAnalyzer.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"

#include "AliHLTCaloRawAnalyzerComponentv3.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTCaloChannelDataStruct.h"
#include "AliHLTCaloMapper.h"
#include "AliHLTCaloSanityInspector.h"
#include "AliAltroRawStreamV3.h"
#include "AliRawReaderMemory.h"
#include "AliAltroRawStreamV3.h"
#include "AliCaloRawStreamV3.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTCaloCrazynessDefinitions.h"
#include "AliHLTCaloChannelRawDataStruct.h"
#include "AliHLTCaloCoordinate.h"
#include "AliLog.h"

#include <vector>
using namespace std;

//ClassImp(AliHLTCaloRawAnalyzerComponentv3);


#include  "TStopwatch.h"
TStopwatch  fgWatch; //CRAP PTH


AliHLTCaloRawAnalyzerComponentv3::AliHLTCaloRawAnalyzerComponentv3(TString det):AliHLTCaloProcessor(),
										AliHLTCaloConstantsHandler(det),
										fAnalyzerPtr(0),
										fMapperPtr(0),     
										fCurrentSpec(-1),
										fDebug(false),
										fSanityInspectorPtr(0),
										fRawReaderMemoryPtr(0),
										fAltroRawStreamPtr(0),
										fAlgorithm(0),  
										fOffset(0),
										fBunchSizeCut(0),
										fMinPeakPosition(0),
										fMaxPeakPosition(100),
										fDoPushBadRawData(false),
										fDoPushRawData(false),
										fRawDataWriter(0)
									
{
  //comment
  fRawReaderMemoryPtr = new AliRawReaderMemory();
  fAltroRawStreamPtr = new AliAltroRawStreamV3(fRawReaderMemoryPtr);
  fSanityInspectorPtr = new AliHLTCaloSanityInspector();
  
  if( fDoPushRawData == true  )
  {
    fRawDataWriter  = new RawDataWriter(fCaloConstants); 
  }
  
  fAltroRawStreamPtr = new AliCaloRawStreamV3(fRawReaderMemoryPtr, det);  
}


AliHLTCaloRawAnalyzerComponentv3::~AliHLTCaloRawAnalyzerComponentv3()
{
  delete fRawReaderMemoryPtr;
  delete fAltroRawStreamPtr;
  delete fRawDataWriter;
  delete fSanityInspectorPtr;
}



int
AliHLTCaloRawAnalyzerComponentv3::DoInit( int argc, const char** argv )
{ 
  //See base class for documentation
  //  fPrintInfo = kFALSE;
  int iResult=0;
  //  fMapperPtr = new AliHLTCaloMapper();
  
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
      if(!strcmp("-pushrawdata", argv[i]))
	{
	  fDoPushRawData = true;
	}
      if(!strcmp("-pushbaddata", argv[i]))
	{
	  fDoPushBadRawData = true;
	}
	if(fDoPushBadRawData && fDoPushRawData) 
	{
	   HLTWarning("fDoPushBadRawData and fDoPushRawData in conflict, using fDoPushRawData");
	   fDoPushBadRawData = false;
	}
	if(!strcmp("-suppressalilogwarnings", argv[i]))
	{
	    AliLog::SetGlobalLogLevel(AliLog::kError);  //PHOS sometimes produces bad data -> Fill up the HLT logs...
	}
    }
 
  return iResult;
}



int 
AliHLTCaloRawAnalyzerComponentv3::DoDeinit()
{
  //comment
  // if(fMapperPtr)
  //   {
  //     delete  fMapperPtr;
  //     fMapperPtr = 0;
  //   }

  // if(fRawReaderMemoryPtr)
  //   {
  //     delete fRawReaderMemoryPtr;
  //     fRawReaderMemoryPtr = 0;
  //   }

  if(fAltroRawStreamPtr)
    {
      delete fAltroRawStreamPtr;
      fAltroRawStreamPtr = 0;
    }

  return 0;
}



void 
AliHLTCaloRawAnalyzerComponentv3::PrintDebugInfo()
{
  //  static int evntcnt = 0;
  static double wlast = -1;
  static double wcurrent = 0;
  // evntcnt  ++;
  
  if( true == fDebug )
    {
      if( fCaloEventCount %1000 == 0  )
	{
	  cout << __FILE__ << __LINE__ << " : Processing event "  << fCaloEventCount << endl; 
	  wlast =  wcurrent;
	  wcurrent = fgWatch.RealTime();
	  ////cout << wlast << ":" << wcurrent << endl;
	  cout << __FILE__ << __LINE__ << "The event rate is " <<  
	    1000/( wcurrent  -  wlast ) << "  Hz" << endl; 	  fgWatch.Start(kFALSE); 
	  //     wlast =  fgWatch.RealTime(); 
	}
    }
}



bool 
AliHLTCaloRawAnalyzerComponentv3::CheckInputDataType(const AliHLTComponentDataType &datatype)
{
  vector <AliHLTComponentDataType> validTypes;
  GetInputDataTypes(validTypes);
  
  // if( true == fDebug)
  // 	{
  // 	  cout << __FILE__ << __LINE__ << "  :  fID  = " << datatype.fID <<  " : fOrigin = " 
  // 	       <<  datatype.fOrigin << endl;
  // 	}

  for(int i=0; i < validTypes.size(); i++ )
    {
      if ( datatype  ==  validTypes.at(i) )
	{
	  return true;
	}
    }
 
  HLTDebug("Invalid Datatype");
  return false;
}


int 
AliHLTCaloRawAnalyzerComponentv3::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					  AliHLTComponentTriggerData& /*trigData*/, 
					  AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  if(!IsDataEvent())
   {
     size = 0;
     return 0;
   }

  if( true == fDebug ) { PrintDebugInfo(); };
   
  
   Int_t blockSize          = -1;
   UInt_t totSize           = 0;
   const AliHLTComponentBlockData* iter = NULL; 
   unsigned long ndx;

   
  
   for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if(  ! CheckInputDataType(iter->fDataType) )
	//          CheckInputDataType	
	{
	  continue;
	}

      //    else
      //  {
      //     InitMapping( iter->fSpecification); 
   
      // blockSize = DoIt(iter, outputPtr, size, totSize); // Processing the block
      
      // if(blockSize == -1) // If the processing returns -1 we are out of buffer and return an error msg.
      // 	{
      // 	  return -ENOBUFS;
      // 	}

      if(iter->fSpecification != fCurrentSpec)
      {
	fCurrentSpec = iter->fSpecification;
	InitMapping(iter->fSpecification);
      }
  
      blockSize = DoIt(iter, outputPtr, size, totSize); // Processing the block
      
   
      totSize += blockSize; //Keeping track of the used size
      AliHLTComponentBlockData bdChannelData;
      FillBlockData( bdChannelData );
      bdChannelData.fOffset = 0; //FIXME
      bdChannelData.fSize = blockSize;
      bdChannelData.fDataType = GetOutputDataType();
      bdChannelData.fSpecification = iter->fSpecification;
      outputBlocks.push_back(bdChannelData);
      outputPtr += blockSize; //Updating position of the output buffer
      //      fCaloEventCount++; 
      //     size = totSize; //telling the framework how much buffer space we have used.
    }

   fCaloEventCount++; 
   size = totSize; //telling the framework how much buffer space we have used.
   
   return 0;
   
}//end DoEvent



Int_t
AliHLTCaloRawAnalyzerComponentv3::DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr, const AliHLTUInt32_t size, UInt_t& totSize)
{
  
  int tmpsize=  0;
  Int_t crazyness          = 0;
  Int_t nSamples           = 0;
  Short_t channelCount     = 0;

  AliHLTCaloChannelDataHeaderStruct *channelDataHeaderPtr = reinterpret_cast<AliHLTCaloChannelDataHeaderStruct*>(outputPtr); 
  AliHLTCaloChannelDataStruct *channelDataPtr = reinterpret_cast<AliHLTCaloChannelDataStruct*>(outputPtr+sizeof(AliHLTCaloChannelDataHeaderStruct)); 
  totSize += sizeof( AliHLTCaloChannelDataHeaderStruct );
  fRawReaderMemoryPtr->SetMemory(         reinterpret_cast<UChar_t*>( iter->fPtr ),  static_cast<ULong_t>( iter->fSize )  );
  fRawReaderMemoryPtr->SetEquipmentID(    fMapperPtr->GetDDLFromSpec(  iter->fSpecification) + fCaloConstants->GetDDLOFFSET() );
  fRawReaderMemoryPtr->Reset();
  fRawReaderMemoryPtr->NextEvent();

  if( fDoPushRawData == true)
    {
     fRawDataWriter->NewEvent( );
    }
  
  if(fAltroRawStreamPtr->NextDDL())
  {
    int cnt = 0;
    //    fOffset = ( fAltroRawStreamPtr->GetAltroCFG1() >> 10 ) & 0xf;
    fOffset = 0;
    while( fAltroRawStreamPtr->NextChannel()  )
      { 
	if(  fAltroRawStreamPtr->GetHWAddress() < 128 || ( fAltroRawStreamPtr->GetHWAddress() ^ 0x800) < 128 ) 
	  {
	    continue; 
	  }
	else
	  {
	    ++ cnt;
	    UShort_t* firstBunchPtr = 0;
	    int chId = fMapperPtr->GetChannelID(iter->fSpecification, fAltroRawStreamPtr->GetHWAddress()); 
	    
	    if( fDoPushRawData == true)
	      {
		fRawDataWriter->SetChannelId( chId );
	      }
	    
	    vector <AliCaloBunchInfo> bvctr;
	    while( fAltroRawStreamPtr->NextBunch() == true )
	      {
		bvctr.push_back( AliCaloBunchInfo( fAltroRawStreamPtr->GetStartTimeBin(), 
						   fAltroRawStreamPtr->GetBunchLength(),
						   fAltroRawStreamPtr->GetSignals() ) );	

		nSamples = fAltroRawStreamPtr->GetBunchLength();
		if( fDoPushRawData == true)
		  {
		    fRawDataWriter->WriteBunchData( fAltroRawStreamPtr->GetSignals(), 
						    nSamples,  fAltroRawStreamPtr->GetEndTimeBin()  );
		  }
		firstBunchPtr = const_cast< UShort_t* >(  fAltroRawStreamPtr->GetSignals()  );
	      }
	     
	    //return 1;
	  
	    totSize += sizeof( AliHLTCaloChannelDataStruct );
	    if(totSize > size)
	      {
		//HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", totSize, size);
		return -1;
	      }

	    //    fAnalyzerPtr->SetData( firstBunchPtr, nSamples);
	    AliCaloFitResults res = fAnalyzerPtr->Evaluate( bvctr,  fAltroRawStreamPtr->GetAltroCFG1(), fAltroRawStreamPtr->GetAltroCFG2() );  
 
	    HLTDebug("Channel energy: %f, max sig: %d, gain = %d, x = %d, z = %d", res.GetAmp(), res.GetMaxSig(), 
		     (chId >> 12)&0x1, chId&0x3f, (chId >> 6)&0x3f);
	    
	    //	      if(fAnalyzerPtr->GetTiming() > fMinPeakPosition && fAnalyzerPtr->GetTiming() < fMaxPeakPosition)
	    {
	      channelDataPtr->fChannelID =  chId;
	      channelDataPtr->fEnergy = static_cast<Float_t>( res.GetAmp()  ) - fOffset;
	      channelDataPtr->fTime = static_cast<Float_t>(  res.GetTof() );
	      channelDataPtr->fCrazyness = static_cast<Short_t>(crazyness);
	      channelCount++;
	      channelDataPtr++; // Updating position of the free output.
	    }   
	  }

	// This should be the functionality for pushing the raw data for the bad channels

	/*		  
            if ( bad )
	      {
                crazyness |= AliHLTCaloCrazynessDefinitions::fgkBadEstBit;
              }

	    channelDataPtr->fChannelID =  chId;
	    channelDataPtr->fEnergy = static_cast<Float_t>(fAnalyzerPtr->GetEnergy()) - fOffset;
	    AliHLTCaloCoordinate c;
	    fMapperPtr->ChannelId2Coordinate(chId, c); 
	    HLTDebug("Got channel, x: %d, z: %d, gain: %d, energy: %f", c.fX, c.fZ, c.fGain, channelDataPtr->fEnergy);
	    channelDataPtr->fTime = static_cast<Float_t>(fAnalyzerPtr->GetTiming());
	    channelDataPtr->fCrazyness = static_cast<Short_t>(crazyness);
	    channelCount++;
	    channelDataPtr++; // Updating position of the free output.
	    
	                     // Push the raw data if we have crazyness and feature to push raw data for bad channels
            if (crazyness && fDoPushBadRawData)

              {
		channelDataPtr->fCrazyness = crazyness | AliHLTCaloCrazynessDefinitions::fgkRawDataBit;
		UShort_t extraSize = sizeof ( AliHLTCaloChannelRawDataStruct ) + ( nSamples-1 ) *sizeof ( UShort_t );
		totSize += extraSize;
		
		if ( totSize > size )
		  {
		    HLTError ( "Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", totSize, size );
                    return -1;
                  }

		AliHLTCaloChannelRawDataStruct *tmpChannel = reinterpret_cast<AliHLTCaloChannelRawDataStruct*> ( channelDataPtr );
                memcpy ( & ( tmpChannel->fData ), firstBunchPtr, nSamples*sizeof ( UShort_t ) );
                channelDataPtr = reinterpret_cast<AliHLTCaloChannelDataStruct*> ( reinterpret_cast<UChar_t*> ( channelDataPtr ) + extraSize );
              }

	*/
      }
    if( fDoPushRawData == true)
      { 
	fRawDataWriter->NewChannel();
      }
  }

//Writing the header
  channelDataHeaderPtr->fNChannels   =  channelCount;
  channelDataHeaderPtr->fAlgorithm   = fAlgorithm;
  channelDataHeaderPtr->fInfo        = 0;
  
  if( fDoPushRawData == true)
    {
      tmpsize += fRawDataWriter->CopyBufferToSharedMemory( (UShort_t *)channelDataPtr, size, totSize);
    }

  channelDataHeaderPtr->fHasRawData = fDoPushRawData;
  HLTDebug("Number of channels: %d", channelCount);
  tmpsize += sizeof(AliHLTCaloChannelDataStruct)*channelCount + sizeof(AliHLTCaloChannelDataHeaderStruct); 
  return  tmpsize;

}



AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::RawDataWriter(AliHLTCaloConstants* cConst) : 
								    fRawDataBuffer(0),
								    fCurrentChannelSize(0),
								    fBufferIndex(0),
								    fBufferSize( 64*56*cConst->GetNGAINS()*cConst->GetALTROMAXSAMPLES() +1000 ),
								    fCurrentChannelIdPtr(0),
								    fCurrentChannelSizePtr(0),
								    fCurrentChannelDataPtr(0),
								    fTotalSize(0)
{
  fRawDataBuffer = new UShort_t[fBufferSize];
  Init();
}


AliHLTCaloRawAnalyzerComponentv3::RawDataWriter::~RawDataWriter()
{
  delete [] fRawDataBuffer;
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
  
