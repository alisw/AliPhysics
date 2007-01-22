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

//#include "AliHLTPHOSModuleMergerComponent.h"

#include "AliHLTPHOSModuleMergerComponent.h"
#include <iostream>
#include "stdio.h"

#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include <cstdlib>
//#include "TH2.h"

//#include "AliHLTTPCDefinitions.h"

const AliHLTComponentDataType  AliHLTPHOSModuleMergerComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}}; //'zero' terminated array
const AliHLTComponentDataType  AliHLTPHOSModuleMergerComponent::outputDataType=kAliHLTVoidDataType;


// AliHLTPHOSModuleMergerComponent g AliHLTPHOSModuleMergerComponent;
//ClassImp( AliHLTPHOSModuleMergerComponent) 
AliHLTPHOSModuleMergerComponent:: AliHLTPHOSModuleMergerComponent():AliHLTProcessor(),  fEventCount(0),  fEquippmentId(0)
{

} 

AliHLTPHOSModuleMergerComponent::~ AliHLTPHOSModuleMergerComponent()
{

}

AliHLTPHOSModuleMergerComponent:: AliHLTPHOSModuleMergerComponent(const  AliHLTPHOSModuleMergerComponent & ) : AliHLTProcessor(),  fEventCount(0),  fEquippmentId(0)
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
  return outputDataType;
}

void
 AliHLTPHOSModuleMergerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  constBase = 30;
  inputMultiplier = 0.1;
}


int  AliHLTPHOSModuleMergerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{


  fEventCount++; 
  return 0;
}//end DoEvent



int
 AliHLTPHOSModuleMergerComponent::DoInit( int argc, const char** argv )
{
  cout << "DOINIT argc =" << argc << endl;
  cout << "DOINIT argv[0] =" << argv[0] << endl;
  cout << "DOINIT argv[1] =" << argv[1] << endl;
  cout << "DOINIT argv[2] =" << argv[2] << endl;
  cout << "DOINIT argv[3] =" << argv[3] << endl;
  cout << "DOINIT argv[4] =" << argv[4] << endl;
  cout << "DOINIT argv[5] =" << argv[5] << endl;
  cout << "DOINIT argv[6] =" << argv[6] << endl;
 
  int equippmentId = atoi(argv[6]);
  cout << "The equipment ID was set to " <<equippmentId << endl;
  

  Reset();
  cout << " AliHLTPHOSModuleMergerComponent::DoInit Creating new  AliRawReaderMemory()" << endl; 


  cout <<" AliHLTPHOSModuleMergerComponent::DoIni  DONE!" << endl;
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
  fEquippmentId = id;
}

int 
AliHLTPHOSModuleMergerComponent::GetEquippmentId()
{
  return  fEquippmentId;
}
