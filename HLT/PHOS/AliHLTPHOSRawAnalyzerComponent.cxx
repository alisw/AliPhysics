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

#include "AliHLTPHOSRawAnalyzerComponent.h"
#include <iostream>
#include "stdio.h"

#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include <cstdlib>
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
//#include "AliHLTPHOSDataHeaderStruct.h"


const AliHLTComponentDataType AliHLTPHOSRawAnalyzerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
int   AliHLTPHOSRawAnalyzerComponent::fEventCount = 0; 


AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent():AliHLTProcessor(),fEquippmentID(0), fRcuX(0), 
fRcuZ(0),fRcuRowOffeset(0), fRcuColOffeset(0),  fModuleID(0), fPHOSRawStream(0), fRawMemoryReader(0), fOutPtr(0)
//AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent():AliHLTProcessor()
{

} 

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



AliHLTPHOSRawAnalyzerComponent::AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & ) : AliHLTProcessor(), 
fEquippmentID(0), fRcuX(0), fRcuZ(0),fRcuRowOffeset(0), fRcuColOffeset(0),  fModuleID(0), fPHOSRawStream(0), fRawMemoryReader(0), fOutPtr(0)
{
}


int 
AliHLTPHOSRawAnalyzerComponent::Deinit()
{
  return 0;
}

int 
AliHLTPHOSRawAnalyzerComponent::DoDeinit()
{
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
AliHLTPHOSRawAnalyzerComponent::GetComponentID()
{
  return "AliPhosTestRaw";
}

void
AliHLTPHOSRawAnalyzerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=inputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSRawAnalyzerComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::gkCellEnergyDataType;
}

void
AliHLTPHOSRawAnalyzerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  constBase = 30;
  inputMultiplier = 0.1;
}


int AliHLTPHOSRawAnalyzerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  Int_t tmpMod            = 0;
  Int_t tmpRow            = 0;
  Int_t tmpCol            = 0;
  Int_t tmpGain           = 0;
  Int_t sampleCnt         = 0;
  Int_t processedChannels = 0;
  UInt_t offset           = 0; 
  UInt_t mysize           = 0;
  UInt_t tSize            = 0;
  Int_t tmpChannelCnt     = 0;
  Int_t tmpStartIndex     = 0;
  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx;

  if((fEventCount % 100) == 0)
    {
      //      cout << "analyzing event: " << fEventCount << endl;
    }

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      mysize = 0;
      offset = tSize;

      if ( iter->fDataType != AliHLTPHOSDefinitions::gkDDLPackedRawDataType )
	{
	  //	  cout << "Warning: data type = is nOT gkDDLPackedRawDataType " << endl;
	  continue;
	}

      fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize );
      analyzerPtr->SetData(fTmpChannelData);
      fOutPtr =  (AliHLTPHOSRcuCellEnergyDataStruct*)outBPtr;
      mysize += sizeof(AliHLTPHOSRcuCellEnergyDataStruct);
      fOutPtr->fRcuX = fRcuX;
      fOutPtr->fRcuZ = fRcuZ;
      fOutPtr->fModuleID = fModuleID;
      tmpChannelCnt = 0;
 
      if(fEventCount%100 ==0)
	{
	  //	  cout <<"Analyzing event: " << fEventCount << endl; 
	}
 
      while(fPHOSRawStream->Next())
	{
	  if (fPHOSRawStream->IsNewHWAddress())
	    {
	      if(processedChannels > 0)
		{
		  analyzerPtr->SetData(fTmpChannelData);
		  analyzerPtr->Evaluate(0, sampleCnt);
		  sampleCnt = 0;
		  fOutPtr->fValidData[tmpChannelCnt].fGain = tmpGain;
		  fOutPtr->fValidData[tmpChannelCnt].fRow  = tmpRow;
		  fOutPtr->fValidData[tmpChannelCnt].fCol  = tmpCol; 
		  fOutPtr->fValidData[tmpChannelCnt].fEnergy  = analyzerPtr->GetEnergy();
		  fOutPtr->fValidData[tmpChannelCnt].fTime    = analyzerPtr->GetTiming();
		  tmpChannelCnt ++;
		  ResetDataPtr(tmpStartIndex, sampleCnt);
		  sampleCnt = 0;
		}

	      tmpMod  =  fPHOSRawStream->GetModule() ;
	      tmpRow  =  fPHOSRawStream->GetRow() - fRcuRowOffeset;
	      tmpCol  =  fPHOSRawStream->GetColumn() - fRcuColOffeset;
	      tmpGain =  fPHOSRawStream->IsLowGain(); 
	      processedChannels ++;
	    }
	  
	  if(sampleCnt == 0)
	    {
	      tmpStartIndex = fPHOSRawStream->GetTime();
	    }
	  
	  fTmpChannelData[fPHOSRawStream->GetTime()] =  fPHOSRawStream->GetSignal();
	  sampleCnt ++;

	}
   
      fOutPtr->fCnt =  tmpChannelCnt;
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;
      bd.fDataType = AliHLTPHOSDefinitions::gkCellEnergyDataType;
      bd.fSpecification = 0xFFFFFFFF;
      outputBlocks.push_back( bd );
      tSize += mysize;
      outBPtr += mysize;
      
      if( tSize > size )
	{
	  Logging( kHLTLogFatal, "HLT::AliHLTPHOSRawAnalyzerComponent::DoEvent", "Too much data",
		   "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
		   , tSize, size );
	  return EMSGSIZE;
	}
    }

  fEventCount++; 
  size = tSize;
  return 0;
}//end DoEvent



int
AliHLTPHOSRawAnalyzerComponent::DoInit( int argc, const char** argv )
{
  int equippmentID = atoi(argv[6]);
  Reset();
  fRawMemoryReader = new AliRawReaderMemory();
  fPHOSRawStream = new  AliCaloRawStream(fRawMemoryReader,"PHOS");
  fRawMemoryReader->SetEquipmentID(equippmentID); 
  SetEquippmentID(equippmentID);
  SetCoordinates(equippmentID);
  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}

void
AliHLTPHOSRawAnalyzerComponent::DumpData(int gain)
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


void
AliHLTPHOSRawAnalyzerComponent::DumpData()
{
  DumpData(0);
}

void
AliHLTPHOSRawAnalyzerComponent::DumpChannelData(Double_t *data)
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

  ResetDataPtr();

} // end Reset



void
AliHLTPHOSRawAnalyzerComponent::ResetDataPtr()
{
  for(int i = 0 ; i< ALTRO_MAX_SAMPLES; i++)
    {
      fTmpChannelData[i] = 0;
    }
}

void
AliHLTPHOSRawAnalyzerComponent::ResetDataPtr(int sampleCnt)
{
  for(int i = 0 ; i< sampleCnt; i++)
    {
      fTmpChannelData[i] = 0;
    }
}

void
AliHLTPHOSRawAnalyzerComponent::ResetDataPtr(int startindex, int sampleCnt)
{
  for(int i = startindex ; i< sampleCnt; i++)
    {
      fTmpChannelData[i] = 0;
    }
}


void 
AliHLTPHOSRawAnalyzerComponent::SetEquippmentID(AliHLTUInt32_t id)
{
  fEquippmentID = id;
}

int 
AliHLTPHOSRawAnalyzerComponent::GetEquippmentID()
{
  return  fEquippmentID;
}

void 
AliHLTPHOSRawAnalyzerComponent::SetCoordinates(AliHLTUInt32_t equippmentID)
{
  int rcuIndex =  (fEquippmentID - 1792)%4;
  fModuleID = (fEquippmentID  -1792 -rcuIndex)/5;

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

  fRcuRowOffeset = 32*fRcuX;
  fRcuColOffeset = 28*fRcuZ;

}
