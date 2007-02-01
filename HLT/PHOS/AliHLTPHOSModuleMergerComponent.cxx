/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSModuleMergerComponent.h"
#include <iostream>
#include "stdio.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include <cstdlib>
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"


const AliHLTComponentDataType  AliHLTPHOSModuleMergerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType  AliHLTPHOSModuleMergerComponent::outputDataType=kAliHLTVoidDataType;


AliHLTPHOSModuleMergerComponent gAliHLTPHOSModuleMergerComponent;
AliHLTPHOSModuleMergerComponent:: AliHLTPHOSModuleMergerComponent():AliHLTProcessor(),  fEventCount(0),  fEquippmentID(0)
{

} 

AliHLTPHOSModuleMergerComponent::~ AliHLTPHOSModuleMergerComponent()
{

}

AliHLTPHOSModuleMergerComponent:: AliHLTPHOSModuleMergerComponent(const  AliHLTPHOSModuleMergerComponent & ) : AliHLTProcessor(),  fEventCount(0),  fEquippmentID(0)
{

}


int 
AliHLTPHOSModuleMergerComponent::Deinit()
{
  return 0;
}

int 
AliHLTPHOSModuleMergerComponent::DoDeinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSModuleMerger DoDeinit");
  return 0;

}

const char* 
AliHLTPHOSModuleMergerComponent::GetComponentID()
{
  return "ModuleMerger";
}

void
 AliHLTPHOSModuleMergerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=inputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSModuleMergerComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::gkCellEnergyDataType;
}

void
AliHLTPHOSModuleMergerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  constBase = 30;
  inputMultiplier = 1;
}


int  AliHLTPHOSModuleMergerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  unsigned long ndx;
  const AliHLTComponentBlockData* iter = NULL;   
  AliHLTPHOSRcuCellEnergyDataStruct *cellDataPtr;

  Reset();

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      int tmpModuleID = 0;
      int tmpRcuX = 0;
      int tmpRcuZ = 0;

      iter = blocks+ndx;
      AliHLTPHOSRcuCellEnergyDataStruct *cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);

      tmpModuleID = cellDataPtr->fModuleID;
      tmpRcuX     = cellDataPtr->fRcuX ;
      tmpRcuZ     = cellDataPtr->fRcuZ;

      for(int row = 0; row<32; row ++)
	{
	  for(int col = 0; col < 28; col ++)
	    {
	      for(int gain=0; gain <2; gain++)
		{
		  fMaxValues[tmpModuleID][row + 32*tmpRcuX][col + 28*tmpRcuZ][gain] =  cellDataPtr->fCellEnergies[row][col][gain];  
	 	}	  
	    }
	}

    }

  DumpData();
  fEventCount++; 
  return 0;
}//end DoEvent



int
AliHLTPHOSModuleMergerComponent::DoInit( int argc, const char** argv )
{
  Reset();

  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}

void
AliHLTPHOSModuleMergerComponent::DumpData()
{
  for(int mod = 0; mod <5; mod ++)
    {
      printf("\n ***********  MODULE %d ************\n", mod);
      for(int row = 0; row < 64; row ++)
	{
	  for(int col = 0; col < 56; col ++)
	    {
	      if( fMaxValues[mod][row][col][0] != 0)
		{ 
		  cout << fMaxValues[mod][row][col][0] << "\t";
		}
	    }
	} 
    }
}


void
AliHLTPHOSModuleMergerComponent::Reset()
{
  for(int mod = 0; mod <5; mod ++)
    {
      for(int row = 0; row < 64; row ++)
	{
	  for(int col = 0; col < 56; col ++)
	    {
	      for(int gain = 0; gain <2; gain ++ )
		{
		  fMaxValues[mod][row][col][gain] = 0;
		}
	    }
	}
    }

  for(int i = 0 ; i< 1008; i++)
    {
      fTmpChannelData[i] = 0;
    }
} // end Reset

void
AliHLTPHOSModuleMergerComponent::ResetDataPtr()
{
  for(int i = 0 ; i< 1008; i++)
    {
      fTmpChannelData[i] = 0;
    }
}


void 
AliHLTPHOSModuleMergerComponent::SetEquippmentId(int id)
{
  fEquippmentID = id;
}

int 
AliHLTPHOSModuleMergerComponent::GetEquippmentId()
{
  return  fEquippmentID;
}


AliHLTComponent*
AliHLTPHOSModuleMergerComponent::Spawn()
{
  return new AliHLTPHOSModuleMergerComponent;
}


