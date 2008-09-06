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

#include "AliHLTPHOSRawAnalyzer.h"
#include "AliHLTPHOSRawAnalyzerComponentv2.h"
#include "AliHLTPHOSChannelDataHeaderStruct.h"
#include "AliHLTPHOSChannelDataStruct.h"
#include "AliHLTPHOSMapper.h"
#include "AliHLTPHOSSanityInspector.h"
#include "AliHLTPHOSBaseline.h"
#include  "AliAltroDecoder.h"    // decoder for altro payload
#include  "AliAltroData.h"       // container for altro payload
#include  "AliAltroBunch.h"      // container for altro bunches


AliHLTPHOSRawAnalyzerComponentv2::AliHLTPHOSRawAnalyzerComponentv2():AliHLTPHOSRcuProcessor(), 
                                                                 fAnalyzerPtr(0), 
                                                                 fMapperPtr(0), 
                                                                 fSanityInspectorPtr(0),
                                                                 fDecoderPtr(0),  
                                                                 fAltroDataPtr(0),
                                                                 fAltroBunchPtr(0),
								 fNOKBlocks(0),
								 fAlgorithm(0)
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
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRawAnalyzerComponen Deinit");
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
  UInt_t blockSize         = 0;
  UInt_t totSize           = 0;

  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;
  bool droppedRaw = true;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if ( iter->fDataType != AliHLTPHOSDefinitions::fgkDDLPackedRawDataType  )
	{
	  HLTDebug("Data block is not of type fgkDDLPackedRawDataType");
	  continue; 
	}

      blockSize = DoIt(iter, outputPtr, size, totSize);
      if(blockSize == -1) 
	{
	  return -ENOBUFS;
	}

      totSize += blockSize;
      HLTDebug("Output data size: %d - Input data size: %d", totSize, iter->fSize);
      AliHLTComponentBlockData bdChannelData;
      FillBlockData( bdChannelData );
      //bdChannelData.fOffset = totSize-blockSize;
      bdChannelData.fOffset = 0; //CRAP
      bdChannelData.fSize = blockSize;
      bdChannelData.fDataType = AliHLTPHOSDefinitions::fgkChannelDataType;
      bdChannelData.fSpecification = iter->fSpecification;
      outputBlocks.push_back(bdChannelData);
    }

  fPhosEventCount++; 

  return 0;
}//end DoEvent

Int_t
AliHLTPHOSRawAnalyzerComponentv2::DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr, const AliHLTUInt32_t size, UInt_t& totSize)
{
 

  Int_t crazyness          = 0;
  Int_t nSamples           = 0;
  const int bunchsizecut   = 5;
 Short_t channelCount     = 0;
      
  AliHLTPHOSChannelDataHeaderStruct *channelDataHeaderPtr = reinterpret_cast<AliHLTPHOSChannelDataHeaderStruct*>(outputPtr); 

  AliHLTPHOSChannelDataStruct *channelDataPtr = reinterpret_cast<AliHLTPHOSChannelDataStruct*>(outputPtr+sizeof(AliHLTPHOSChannelDataHeaderStruct)); 
  Short_t channelSize = sizeof(AliHLTPHOSChannelDataStruct);
 
  totSize += sizeof(AliHLTPHOSChannelDataHeaderStruct);
 
  fDecoderPtr->SetMemory(reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize);
  fDecoderPtr->Decode();
      
  while(fDecoderPtr->NextChannel(fAltroDataPtr) == true )
    {          
      if(fAltroDataPtr->GetDataSize() != 0 )
	{
	  int bunchcount = 0;
	  
	  while(fAltroDataPtr->NextBunch(fAltroBunchPtr) == true);
	  
	  if(fAltroBunchPtr->GetBunchSize() > bunchsizecut && bunchcount == 0)
	    {
	      totSize += channelSize;
	      if(totSize > size)
		{
		  HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", totSize, size);
		  return -1;
		}
	      
	      nSamples = fAltroBunchPtr->GetBunchSize();
	      
	      fNOKBlocks ++;
	      
	      crazyness = fSanityInspectorPtr->CheckInsanity(static_cast<const UInt_t*>(fAltroBunchPtr->GetData()), static_cast<const Int_t>(fAltroBunchPtr->GetBunchSize()));
	      
	      fAnalyzerPtr->SetData(fAltroBunchPtr->GetData(), fAltroBunchPtr->GetBunchSize());   
	      fAnalyzerPtr->Evaluate(0, fAltroBunchPtr->GetBunchSize());  
	      
	      channelDataPtr->fChannelID = fMapperPtr->GetChannelID(iter->fSpecification, fAltroDataPtr->GetHadd());
	      channelDataPtr->fEnergy = static_cast<Float_t>(fAnalyzerPtr->GetEnergy());
	      channelDataPtr->fTime = static_cast<Float_t>(fAnalyzerPtr->GetTiming());
	      channelDataPtr->fCrazyness = static_cast<Short_t>(crazyness);
	      channelCount++;
	      channelDataPtr++;
	    }
	  bunchcount++;
	}
    }

  channelDataHeaderPtr->fNChannels  = channelCount;
  channelDataHeaderPtr->fAlgorithm  = fAlgorithm;
  channelDataHeaderPtr->fInfo       = 0;
  channelDataHeaderPtr->fHasRawData = false;

  HLTDebug("Number of channels: %d", channelCount);

  return channelSize*channelCount + sizeof(AliHLTPHOSChannelDataHeaderStruct);
}

int
AliHLTPHOSRawAnalyzerComponentv2::DoInit( int argc, const char** argv )
{ 

  //See base class for documentation
  fPrintInfo = kFALSE;
  int iResult=0;
  fMapperPtr = new AliHLTPHOSMapper();
 
  if(fMapperPtr->GetIsInitializedMapping() == false)
    {
      Logging(kHLTLogFatal, __FILE__ , IntToChar(  __LINE__ ) , "AliHLTPHOSMapper::Could not initial mapping from file %s, aborting", fMapperPtr->GetFilePath());
      return -4;
    }

  return iResult;
}
