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
#include "stdio.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include <cstdlib>
#include "AliHLTPHOSRcuChannelDataStruct.h"

 AliHLTPHOSDDLDecoderComponent  gAliHLTPHOSDDLDecoderComponent;

const AliHLTComponentDataType AliHLTPHOSDDLDecoderComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
int   AliHLTPHOSDDLDecoderComponent::fgEventCount = 0; 


AliHLTPHOSDDLDecoderComponent::AliHLTPHOSDDLDecoderComponent():AliHLTProcessor(), fEquippmentID(0), fRcuX(0), 
fRcuZ(0),fRcuZOffset(0), fRcuXOffset(0),  fModuleID(0), fSendChannelData(kFALSE), fPHOSRawStream(0), fRawMemoryReader(0), fOutPtr(0)
{

} 

AliHLTPHOSDDLDecoderComponent::~AliHLTPHOSDDLDecoderComponent()
{
  if(fRawMemoryReader != 0)
    {
      delete fRawMemoryReader;
    }
    if(fPHOSRawStream != 0)
    {
      delete fPHOSRawStream;
    }
}


AliHLTPHOSDDLDecoderComponent::AliHLTPHOSDDLDecoderComponent(const AliHLTPHOSDDLDecoderComponent & ) : AliHLTProcessor(), 
fEquippmentID(0), fRcuX(0), fRcuZ(0),fRcuZOffset(0), fRcuXOffset(0),  fModuleID(0),  fSendChannelData(kFALSE), fPHOSRawStream(0), fRawMemoryReader(0), fOutPtr(0)
{
}


int 
AliHLTPHOSDDLDecoderComponent::Deinit()
{
  cout <<  "Deinit" << endl;
  return 0;
}


int 
AliHLTPHOSDDLDecoderComponent::DoDeinit()
{
  cout << "DoDeinit" << endl;
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

}


const char* 
AliHLTPHOSDDLDecoderComponent::GetComponentID()
{
  return "PhosDDLDecoder";
}



void
AliHLTPHOSDDLDecoderComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSDDLDecoderComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkCellEnergyDataType;
}

void
AliHLTPHOSDDLDecoderComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  constBase = 30;
  inputMultiplier = 1;
}


  
int
AliHLTPHOSDDLDecoderComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  AliHLTUInt8_t tmpMod    = 0;
  AliHLTUInt8_t tmpZ      = 0;
  AliHLTUInt8_t tmpX      = 0;
  AliHLTUInt8_t tmpGain   = 0;
  Int_t sampleCnt         = 0;
  Int_t processedChannels = 0;
  UInt_t offset           = 0; 
  UInt_t mysize           = 0;
  UInt_t tSize            = 0;
  Int_t tmpChannelCnt     = 0;
  Int_t tmpStartIndex     = 0;
  AliHLTUInt8_t* outBPtr;
  unsigned long first;
  unsigned long last;
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

  fgEventCount++; 
  
  if(fPrintInfo == kTRUE)
    {
      if(fgEventCount%fPrintInfoFrequncy == 0)
      	{
	  cout <<"Analyzing event " <<  fgEventCount  << "for Equippment " << fEquippmentID << endl; 
	}  
    }
  size = tSize;
  return 0;
}//end DoEvent


int
AliHLTPHOSDDLDecoderComponent::DoInit( int argc, const char** argv )
{
  fSendChannelData = kFALSE;
  fPrintInfo = kFALSE;
  fRawMemoryReader = new AliRawReaderMemory();
  fPHOSRawStream = new  AliCaloRawStream(fRawMemoryReader,"PHOS");
  fPHOSRawStream->SetOldRCUFormat(kFALSE);
  int iResult=0;
  TString argument="";
  Bool_t isSetEquippmentID = kFALSE;

  for(int i=0; i<argc && iResult>=0; i++) 
    {
      argument=argv[i];
      
      if (argument.IsNull()) 
	{
	  continue;
	}
                         
    if (argument.CompareTo("-equipmentID") == 0) 
	{
	  cout << "AliHLTPHOSDDLDecoderComponent:DoInit  argument = -equipmentID   "  <<endl;  
	  if(i+1 <= argc)
	    {
	      fEquippmentID = atoi(argv[i+1]);
	      cout << "AliHLTPHOSDDLDecoderComponent:DoInit  setting equippment ID to  " << fEquippmentID <<endl;
	      fRawMemoryReader->SetEquipmentID(fEquippmentID); 
	      SetEquippmentID(fEquippmentID);
	      SetCoordinates(fEquippmentID);
	      isSetEquippmentID = kTRUE;
	    }
	  else
	    {
	       iResult= -1;
	       Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuHistogramProducerComponent::DoInt( int argc, const char** argv )", "Missing argument",
			"The argument -equippmentID expects a number");
	       return  iResult;   
	    }
	}
      
    
    if (argument.CompareTo("-datatype") == 0) 
      {
	if(i+1 <= argc)
	  {
	    //	    i++;
	    argument=argv[i+1];
	    if(argument.CompareTo("channeldata") == 0)
	      {
		cout << "AliHLTPHOSDDLDecoderComponent::DoIni  setting sendChannelData = kTRUE "<< endl; 
		fSendChannelData = kTRUE;
	      }
	  }	
      }
    

    if (argument.CompareTo("-printinfo") == 0) 
      {
	if(i+1 <= argc)
	  {
	    argument=argv[i+1];
	    fPrintInfoFrequncy = atoi(argv[i+1]);
	    fPrintInfo = kTRUE;
	    cout << "AliHLTPHOSDDLDecoderComponent::DoIni  setting printinfo = kTRUE, with update frequency every  "<< fPrintInfoFrequncy << "th event" <<endl; 
	  }
	else
	  {
	    cout << "WARNING: asking for event info, but no update frequency is specified, otipn is ignored" << endl;
	  }
	//     }	
      }
 
    }


  if(isSetEquippmentID == kFALSE)
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuHistogramProducerComponent::DoInt( int argc, const char** argv )", "Missing argument",
	       "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>");
      iResult = -2; 
    }


  return 0;
}



void 
AliHLTPHOSDDLDecoderComponent::SetEquippmentID(AliHLTUInt16_t id)
{
  fEquippmentID = id;
}


AliHLTUInt16_t
AliHLTPHOSDDLDecoderComponent::GetEquippmentID()
{
  return  fEquippmentID;
}


void 
AliHLTPHOSDDLDecoderComponent::SetCoordinates(AliHLTUInt16_t equippmentID)
{
  int rcuIndex =  (fEquippmentID - 1792)%N_RCUS_PER_MODULE;
  fModuleID = (fEquippmentID  -1792 -rcuIndex)/N_RCUS_PER_MODULE;
  
  if(rcuIndex == 0)
    {
      fRcuX = 0; 
      fRcuZ = 0;
    }

  if(rcuIndex == 1)
    {
      fRcuX = 0; 
      fRcuZ = 1;
    }
 
  if(rcuIndex == 2)
    {
      fRcuX = 1; 
      fRcuZ = 0;
    }

  if(rcuIndex == 3)
    {
      fRcuX = 1; 
      fRcuZ = 1;
    }



  fRcuZOffset =  N_ZROWS_RCU*fRcuZ;
  fRcuXOffset =  N_XCOLUMNS_RCU*fRcuX;

  cout <<"********InitInfo************"<< endl;
  cout <<"AliHLTPHOSDDLDecoderComponent::SetCoordinate"<< endl;
  cout <<"Equpippment ID =\t"<< fEquippmentID <<endl;
  cout <<"Module ID =\t"<<  (int)fModuleID<<endl;
  cout <<"RCUX =\t\t" << (int)fRcuX << endl;
  cout <<"RCUZ =\t\t" << (int)fRcuZ << endl;
  cout <<"RcuZOffset = \t" <<  (int)fRcuZOffset << endl;
  cout <<"RcuXOffset = \t" <<  (int)fRcuXOffset << endl << endl;
}

AliHLTComponent*
AliHLTPHOSDDLDecoderComponent::Spawn()
{
  return new AliHLTPHOSDDLDecoderComponent;
}
