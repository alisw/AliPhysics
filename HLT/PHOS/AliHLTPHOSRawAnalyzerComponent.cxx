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


/// @class AliHLTPHOSRawAnalyzerComponent
/// Base class of PHOS HLT online raw analysis component.
/// The class provides a common interface for the implementation of PHOS 
/// HLT raw data
/// processors components. The class is intended for processing of 
/// arrays of raw data samples to evaluate energy and timing.
/// The Energy will be given in entities of ADC leves ranging from 0 to
/// 1023. Timing will be given in entities of samples periods.
/// Drived clases  must implement the fucntions
/// - @ref GetComponentID
/// - @ref Spawn

#include "AliHLTPHOSRawAnalyzer.h"
#include "AliHLTPHOSRawAnalyzerComponent.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSRcuChannelDataStruct.h"
#include "AliHLTDDLDecoder.h"
#include "AliHLTAltroData.h"

#include "AliHLTPHOSMapper.h"

using namespace std;


//_________________________________________________________________________________________________
AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent():AliHLTPHOSProcessor(), fAnalyzerPtr(0), 
fSendChannelData(kFALSE),fPHOSRawStream(0), fRawMemoryReader(0), fOutPtr(0)
{
  fMapperPtr = new AliHLTPHOSMapper();
} 

//_________________________________________________________________________________________________
AliHLTPHOSRawAnalyzerComponent::~AliHLTPHOSRawAnalyzerComponent()
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

//_________________________________________________________________________________________________
AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & ) : AliHLTPHOSProcessor(), fAnalyzerPtr(0), 
fSendChannelData(kFALSE),fPHOSRawStream(0), fRawMemoryReader(0), fOutPtr(0)
{

}

//_________________________________________________________________________________________________
int 
AliHLTPHOSRawAnalyzerComponent::Deinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSRawAnalyzerComponen Deinit");

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

//_________________________________________________________________________________________________
const char* 
AliHLTPHOSRawAnalyzerComponent::GetComponentID()
{
  return "AliPhosTestRaw";
}


//_________________________________________________________________________________________________
void
AliHLTPHOSRawAnalyzerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

//_________________________________________________________________________________________________
AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkCellEnergyDataType;
}


//_________________________________________________________________________________________________
void
AliHLTPHOSRawAnalyzerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  constBase = 30;
  inputMultiplier = 1;
}

//_________________________________________________________________________________________________
int 
AliHLTPHOSRawAnalyzerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
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
      offset = tSize;

      if ( iter->fDataType != AliHLTPHOSDefinitions::fgkDDLPackedRawDataType )
	{
	  cout <<"WARNING: notAliHLTPHOSDefinitions::fgkDDLPackedRawDataTyp  "  << endl;
	  //	  continue;
	}

      fDecoderPtr->SetMemory(reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize);
      fDecoderPtr->Decode();
      fOutPtr =  (AliHLTPHOSRcuCellEnergyDataStruct*)outBPtr;

      //      int tmpCnt = 0;

      while( fDecoderPtr->NextChannel(fAltroDataPtr) == true )
	{

	  fAnalyzerPtr->SetData(fAltroDataPtr->fData);
	  fAnalyzerPtr->Evaluate(0, fAltroDataPtr->fDataSize -2);  

	  fOutPtr->fValidData[tmpChannelCnt].fGain = fMapperPtr->ALTRO_MAP[fAltroDataPtr->fHadd].gain;
	  fOutPtr->fValidData[tmpChannelCnt].fZ  = fMapperPtr->ALTRO_MAP[fAltroDataPtr->fHadd].row;
	  fOutPtr->fValidData[tmpChannelCnt].fX  = fMapperPtr->ALTRO_MAP[fAltroDataPtr->fHadd].col; 
	  fOutPtr->fValidData[tmpChannelCnt].fEnergy  = (float)fAnalyzerPtr->GetEnergy();
	  fOutPtr->fValidData[tmpChannelCnt].fTime    = (float)fAnalyzerPtr->GetTiming();
	  tmpChannelCnt ++;
	  //	  tmpCnt ++; 
	}

      mysize += sizeof(AliHLTPHOSRcuCellEnergyDataStruct);

      fOutPtr->fCnt =  tmpChannelCnt;
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;
 
      bd.fDataType = AliHLTPHOSDefinitions::fgkCellEnergyDataType;
      bd.fSpecification = 0xFFFFFFFF;
      outputBlocks.push_back( bd );
 
      tSize += mysize;
      outBPtr += mysize;
      
      if( tSize > size )
      	{
	  cout <<"kHLTLogFatal, HLT::AliHLTPHOSRawAnalyzerComponent::DoEvent Too much dataData written over allowed buffer. Amount written:"
	       << tSize << " allowed" << size << endl;
      	  Logging( kHLTLogFatal, "HLT::AliHLTPHOSRawAnalyzerComponent::DoEvent", "Too much data",
      		   "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
      		   , tSize, size );
      	  return EMSGSIZE;
      	}
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



//_________________________________________________________________________________________________
int
AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv )
{

cout <<"AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv ) "<< endl;

  fAltroDataPtr = new AliHLTAltroData();
  fDecoderPtr = new AliHLTDDLDecoder();
  fSendChannelData = kFALSE;
  fPrintInfo = kFALSE;
  Reset();
  fRawMemoryReader = new AliRawReaderMemory();
  fPHOSRawStream = new  AliCaloRawStream(fRawMemoryReader,"PHOS");

  //  fPHOSRawStream->SetOldRCUFormat(kFALSE);

  fPHOSRawStream->SetOldRCUFormat(kTRUE);

  int iResult=0;
  TString argument="";
  iResult = ScanArguments(argc, argv);


  if(fIsSetEquippmentID == kFALSE)
    {
      cout << "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>" << endl;
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuHistogramProducerComponent::DoInt( int argc, const char** argv )", "Missing argument",
	       "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>");
      iResult = -3; 
    }
  else
    {
      fRawMemoryReader->SetEquipmentID(fkEquippmentID);
    }

  return iResult;


return 0;

}

/*
//_________________________________________________________________________________________________
void
AliHLTPHOSRawAnalyzerComponent::DumpData(int gain) const
{
  for(int mod = 0; mod < N_MODULES; mod ++)
    {
      printf("\n ***********  MODULE %d ************\n", mod);
      for(int row = 0; row <  N_ROWS_MOD; row ++)
	{
	  for(int col = 0; col <  N_COLUMNS_MOD; col ++)
	    {
	      if( fMaxValues[mod][row][col][0] != 0)
		{ 
		  cout << fMaxValues[mod][row][col][gain] << "\t";
		}
	    }
	} 
    }
}
*/


//_________________________________________________________________________________________________

/*
void
AliHLTPHOSRawAnalyzerComponent::DumpChannelData(Double_t *data) const
{
  cout << endl;
  for(int i=0; i<  ALTRO_MAX_SAMPLES; i++)
    {
      if (data[i] != 0)
	{
	  cout <<i <<"\t";
	}
    }
  cout << endl;
  
  for(int i=0; i<  ALTRO_MAX_SAMPLES; i++)
    {
      if (data[i] != 0)
	{
	  cout <<data[i] <<"\t";
	}
    }
  
  cout << endl;
}
*/


//_________________________________________________________________________________________________
void
AliHLTPHOSRawAnalyzerComponent::Reset()
{
  for(int mod = 0; mod < N_MODULES; mod ++)
    {
      for(int row = 0; row < N_ROWS_MOD; row ++)
	{
	  for(int col = 0; col < N_COLUMNS_MOD; col ++)
	    {
	      for(int gain = 0; gain < N_GAINS; gain ++ )
		{
		  fMaxValues[mod][row][col][gain] = 0;
		}
	    }
	}
    }

  ResetDataPtr(0, ALTRO_MAX_SAMPLES);

} // end Reset




//_________________________________________________________________________________________________
void
AliHLTPHOSRawAnalyzerComponent::ResetDataPtr(int startindex, int sampleCnt)
{
  for(int i = startindex ; i< sampleCnt; i++)
    {
      fTmpChannelData[i] = 0;
    }
}

