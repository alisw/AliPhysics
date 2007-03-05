/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Boris Polichtchouk & Per Thomas Hille for the ALICE           *
 * offline/HLT Project. Contributors are mentioned in the code where      *
 * appropriate.                                                           *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSHistogramProducerComponent.h"
#include <iostream>
#include "stdio.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include <cstdlib>
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSModuleCellAccumulatedEnergyDataStruct.h"



const AliHLTComponentDataType  AliHLTPHOSHistogramProducerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType  AliHLTPHOSHistogramProducerComponent::outputDataType=kAliHLTVoidDataType;


AliHLTPHOSHistogramProducerComponent gAliHLTPHOSHistogramProducerComponent;

/*************************************************************************
* Class AliHLTPHOSHistogramProducerComponent accumulating histograms     *
* with amplitudes per PHOS channel                                       *
* It is intended to run at the HLT farm                                  *
* and it fills the histograms with amplitudes per channel.               * 
* Usage example see in PHOS/macros/Shuttle/AliPHOSCalibHistoProducer.C   *
**************************************************************************/
AliHLTPHOSHistogramProducerComponent:: AliHLTPHOSHistogramProducerComponent():AliHLTProcessor(), fEventCount(0),  fEquippmentID(0)
{
  Reset();
} 


AliHLTPHOSHistogramProducerComponent::~ AliHLTPHOSHistogramProducerComponent()
{

}


AliHLTPHOSHistogramProducerComponent::AliHLTPHOSHistogramProducerComponent(const  AliHLTPHOSHistogramProducerComponent & ) : AliHLTProcessor(), fEventCount(0),  fEquippmentID(0)
{

}


int 
AliHLTPHOSHistogramProducerComponent::Deinit()
{
  return 0;
}


int 
AliHLTPHOSHistogramProducerComponent::DoDeinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSHistogramProducer DoDeinit");
  return 0;
}


const char* 
AliHLTPHOSHistogramProducerComponent::GetComponentID()
{
  return "HistogramProducer";
}


void
 AliHLTPHOSHistogramProducerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=inputDataTypes;
  while (pType->fID!=0) 
    {
      list.push_back(*pType);
      pType++;
    }
}


AliHLTComponentDataType 
AliHLTPHOSHistogramProducerComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::gkCellEnergyDataType;
}


void
AliHLTPHOSHistogramProducerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  constBase = 30;
  inputMultiplier = 1;
}


int  AliHLTPHOSHistogramProducerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  unsigned long ndx;
  UInt_t offset           = 0; 
  UInt_t mysize           = 0;
  UInt_t tSize            = 0;

  const AliHLTComponentBlockData* iter = NULL;   
  AliHLTPHOSRcuCellEnergyDataStruct *cellDataPtr;
  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  fOutPtr =  (AliHLTPHOSModuleCellAccumulatedEnergyDataStruct*)outBPtr;

  int tmpCnt;
  AliHLTUInt8_t tmpModuleID;
  AliHLTUInt8_t tmpRcuX;
  AliHLTUInt8_t tmpRcuZ;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);
      tmpCnt =  cellDataPtr->fCnt;

      tmpModuleID = cellDataPtr->fModuleID;
      tmpRcuX     = cellDataPtr->fRcuX ;
      tmpRcuZ     = cellDataPtr->fRcuZ;

      fOutPtr->fModuleID = tmpModuleID;

      int tmpGain;
      int tmpZ;
      int tmpX;

      //      for(int i= 0; i< tmpCnt; i ++)
      for(int i= 0; i <= tmpCnt; i ++)
	{
	  tmpZ =  cellDataPtr->fValidData[i].fZ + N_ZROWS_RCU*tmpRcuZ;
	  tmpX =  cellDataPtr->fValidData[i].fX + N_XCOLUMNS_RCU*tmpRcuX;

	  if(cellDataPtr->fValidData[i].fGain == HIGH_GAIN)
	    {
	      fAccumulatedValues[tmpZ][tmpX][HIGH_GAIN] +=  cellDataPtr->fValidData[i].fEnergy;
	      fHits[tmpZ][tmpX][HIGH_GAIN] ++;
	    }
	  else if(cellDataPtr->fValidData[i].fGain == LOW_GAIN)
	    {
	      fAccumulatedValues[tmpZ][tmpX][LOW_GAIN] +=  cellDataPtr->fValidData[i].fEnergy;
	      fHits[tmpZ][tmpX][LOW_GAIN] ++;

	    }
	}
    }


  for(int z=0;  z < N_ZROWS_MOD; z ++ )
    {
      for(int x = 0; x < N_XCOLUMNS_MOD; x ++)
	{
	  for(int gain =0; gain < N_GAINS; gain ++)
	    {
	      fOutPtr->fAccumulatedEnergies[z][x][gain] = fAccumulatedValues[z][x][gain];
	      fOutPtr->fHits[z][x][gain] = fHits[z][x][gain];
	    }
	} 
    }


  //pushing data to shared output memory
  mysize += sizeof(AliHLTPHOSModuleCellAccumulatedEnergyDataStruct);
  AliHLTComponentBlockData bd;
  FillBlockData( bd );
  bd.fOffset = offset;
  bd.fSize = mysize;
  bd.fDataType = AliHLTPHOSDefinitions::gkCellAccumulatedEnergyDataType;
  bd.fSpecification = 0xFFFFFFFF;
  outputBlocks.push_back( bd );
  tSize += mysize;
  outBPtr += mysize;

  if( tSize > size )
    {
      Logging( kHLTLogFatal, "HLT::AliHLTHistogramProducerComponent::DoEvent", "Too much data",
	       "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
	       , tSize, size );
      return EMSGSIZE;
    }

  fEventCount++; 
  return 0;
}//end DoEvent


int
AliHLTPHOSHistogramProducerComponent::DoInit( int argc, const char** argv )
{
  Reset();
  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}


void
AliHLTPHOSHistogramProducerComponent::DumpData(int gain)
{

  if(gain < 0 || gain >  N_GAINS)
    {
      cout <<"AliHLTPHOSHistogramProducerComponent::DumpDat: Error, gain must be between " << 0 << "and" << N_GAINS << endl;
    }
  
  for(int mod = 0; mod < N_MODULES; mod ++)
    {
      if(gain == HIGH_GAIN)
	{
	  cout << endl <<" ***********  MODULE" << mod << "****HIGH_GAIN" <<"************" << endl;
	}
      else if(gain == LOW_GAIN)
	{
	  cout << endl <<" ***********  MODULE" << mod << "****LOW_GAIN" <<"************" << endl;
	}
      
      for(int row = 0; row < N_ROWS_MOD; row ++)
	{
	  for(int col = 0; col < N_COLUMNS_MOD; col ++)
	    {
	      if(fAccumulatedValues[row][col][0] != 0)
		{ 
		  cout << fAccumulatedValues[row][col][0] << "\t";
		}
	    }
	} 
    }
}



void
AliHLTPHOSHistogramProducerComponent::Reset()
{
  for(int mod = 0; mod < N_MODULES; mod ++)
    {
      for(int row = 0; row <  N_ROWS_MOD; row ++)
	{
	  for(int col = 0; col < N_COLUMNS_MOD; col ++)
	    {
	      for(int gain = 0; gain < N_GAINS; gain ++ )
		{ 
		  fAccumulatedValues[row][col][gain] = 0;
		  fHits[row][col][gain] = 0;	  
		}
	    }
	}
    }

  for(int i = 0 ; i<  ALTRO_MAX_SAMPLES; i++)
    {
      fTmpChannelData[i] = 0;
    }
} // end Reset

void
AliHLTPHOSHistogramProducerComponent::ResetDataPtr()
{
  for(int i = 0 ; i<  ALTRO_MAX_SAMPLES; i++)
    {
      fTmpChannelData[i] = 0;
    }
}


void 
AliHLTPHOSHistogramProducerComponent::SetEquippmentId(int id)
{
  fEquippmentID = id;
}

int 
AliHLTPHOSHistogramProducerComponent::GetEquippmentId()
{
  return  fEquippmentID;
}


AliHLTComponent*
AliHLTPHOSHistogramProducerComponent::Spawn()
{
  return new AliHLTPHOSHistogramProducerComponent;
}


