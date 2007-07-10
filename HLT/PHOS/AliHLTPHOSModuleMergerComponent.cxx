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
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include <cstdio>
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"


const AliHLTComponentDataType  AliHLTPHOSModuleMergerComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType  AliHLTPHOSModuleMergerComponent::fgkOutputDataType=kAliHLTVoidDataType;

AliHLTPHOSModuleMergerComponent gAliHLTPHOSModuleMergerComponent;

//_____________________________________________________________________________________________________
AliHLTPHOSModuleMergerComponent:: AliHLTPHOSModuleMergerComponent():AliHLTPHOSProcessor(),  fPhosEventCount(0),  fEquippmentID(0)
{

} 


//_____________________________________________________________________________________________________
AliHLTPHOSModuleMergerComponent::~ AliHLTPHOSModuleMergerComponent()
{

}


//_____________________________________________________________________________________________________
AliHLTPHOSModuleMergerComponent::AliHLTPHOSModuleMergerComponent(const  AliHLTPHOSModuleMergerComponent & ) : AliHLTPHOSProcessor(),  fPhosEventCount(0),  fEquippmentID(0)
{

}


//_____________________________________________________________________________________________________
int 
AliHLTPHOSModuleMergerComponent::Deinit()
{
  return 0;
}


//_____________________________________________________________________________________________________
int 
AliHLTPHOSModuleMergerComponent::DoDeinit()
{
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSModuleMerger DoDeinit");
  return 0;

}


//_____________________________________________________________________________________________________
const char* 
AliHLTPHOSModuleMergerComponent::GetComponentID()
{
  return "ModuleMerger";
}


//_____________________________________________________________________________________________________
void
AliHLTPHOSModuleMergerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) 
    {
      list.push_back(*pType);
      pType++;
    }
}


//_____________________________________________________________________________________________________
AliHLTComponentDataType 
AliHLTPHOSModuleMergerComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkCellEnergyDataType;
}


//_____________________________________________________________________________________________________
void
AliHLTPHOSModuleMergerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  ///
  constBase = 30;
  inputMultiplier = 1;
}


int  AliHLTPHOSModuleMergerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  //Merging of data from 4 RCUS to one module
  
  unsigned long ndx;
  const AliHLTComponentBlockData* iter = NULL;   
  AliHLTPHOSRcuCellEnergyDataStruct *cellDataPtr;

  Reset();

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      int tmpModuleID = 0;
      int tmpRcuX = 0;
      int tmpRcuZ = 0;
      int tmpCnt =  cellDataPtr->fCnt;
      iter = blocks+ndx;
      AliHLTPHOSRcuCellEnergyDataStruct *cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);
      tmpModuleID = cellDataPtr->fModuleID;
      tmpRcuX     = cellDataPtr->fRcuX ;
      tmpRcuZ     = cellDataPtr->fRcuZ;

      for(int i= 0; i< tmpCnt; tmpCnt ++)
	{
	  if(cellDataPtr->fValidData[i].fGain == HIGH_GAIN)
	    {
	      fMaxValues[tmpModuleID][ cellDataPtr->fValidData[i].fZ +  N_ZROWS_RCU*tmpRcuZ][ cellDataPtr->fValidData[i].fX + N_XCOLUMNS_RCU*tmpRcuX][HIGH_GAIN] =  cellDataPtr->fValidData[i].fEnergy;
	    }
	  else if(cellDataPtr->fValidData[i].fGain == LOW_GAIN)
	    {
	      fMaxValues[tmpModuleID][ cellDataPtr->fValidData[i].fZ +  N_ROWS_RCU*tmpRcuZ][ cellDataPtr->fValidData[i].fX +N_COLUMNS_RCU*tmpRcuX][LOW_GAIN] =  cellDataPtr->fValidData[i].fEnergy;
	    }
	}
      
    }

  DumpData(1);
  fPhosEventCount++; 
  return 0;
  
}//end DoEvent


//_____________________________________________________________________________________________________
int
AliHLTPHOSModuleMergerComponent::DoInit( int argc, const char** argv )
{
  //See base classs for documenation
  Reset();

  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}


//_____________________________________________________________________________________________________
void
AliHLTPHOSModuleMergerComponent::DumpData(int gain)
{
  if(gain < 0 || gain >  N_GAINS)
    {
      cout <<"AliHLTPHOSModuleMergerComponent::DumpDat: Error, gain must be between " << 0 << "and" << N_GAINS << endl;
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
	     	      if( fMaxValues[mod][row][col][0] != 0)
		{ 
		 		  cout << fMaxValues[mod][row][col][0] << "\t";
		}
	    }
	} 
    }
}



//_____________________________________________________________________________________________________
void
AliHLTPHOSModuleMergerComponent::Reset()
{
  for(int mod = 0; mod < N_MODULES; mod ++)
    {
      for(int row = 0; row <  N_ROWS_MOD; row ++)
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

  for(int i = 0 ; i<  ALTRO_MAX_SAMPLES; i++)
    {
      fTmpChannelData[i] = 0;
    }
} // end Reset


//_____________________________________________________________________________________________________
void
AliHLTPHOSModuleMergerComponent::ResetDataPtr()
{
  for(int i = 0 ; i<  ALTRO_MAX_SAMPLES; i++)
    {
      fTmpChannelData[i] = 0;
    }
}


//_____________________________________________________________________________________________________
void 
AliHLTPHOSModuleMergerComponent::SetEquippmentId(int id)
{
  fEquippmentID = id;
}


//_____________________________________________________________________________________________________
const int 
AliHLTPHOSModuleMergerComponent::GetEquippmentId() const
{
  return  fEquippmentID;
}


//_____________________________________________________________________________________________________
AliHLTComponent*
AliHLTPHOSModuleMergerComponent::Spawn()
{
  return new AliHLTPHOSModuleMergerComponent;
}


