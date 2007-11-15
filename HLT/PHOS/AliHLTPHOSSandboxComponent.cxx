/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliHLTPHOSSandboxComponent.h"
#include "AliHLTPHOSProcessor.h"
#include "AliHLTPHOSChannelCounter.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSDefinitions.h"

const AliHLTComponentDataType AliHLTPHOSSandboxComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};

AliHLTPHOSSandboxComponent gAliHLTPHOSSandboxComponent;

AliHLTPHOSSandboxComponent::AliHLTPHOSSandboxComponent() :
  AliHLTPHOSProcessor(),
  fEvtCnt(0),
  fChannelCounterPtr(0)
{
  //comment
}

AliHLTPHOSSandboxComponent::~AliHLTPHOSSandboxComponent()
{
  //comment
}
/*
AliHLTPHOSSandboxComponent::AliHLTPHOSSandboxComponent(const AliHLTPHOSSandboxComponent &) :
  AliHLTDataSink(),
  fEvtCnt(0),
  fChannelCounterPtr(0)
{

  //Copy constructor, not implemented

}
*/

int 
AliHLTPHOSSandboxComponent::Deinit()
{
  //comment
  printf("AliHLTPHOSSandboxComponent::Deinit()\n");

  //fChannelCounterPtr->PrintOutOfSyncChannels(fEvtCnt);
  fChannelCounterPtr->FillHistograms(fEvtCnt);
  fChannelCounterPtr->WriteHistograms("/opt/HLT-public/rundir/channelcount.root");
  return 0;
}

Int_t
AliHLTPHOSSandboxComponent::DoDeinit()
{
  //comment
  //Deinitialize the component
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSSandboxComponent DoDeinit");

  return 0;
}



const char*
AliHLTPHOSSandboxComponent::GetComponentID()
{
  //comment
  return "PhosSandbox";
}

void
AliHLTPHOSSandboxComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{ 
 //Get datatypes for input
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType); 
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSSandboxComponent::GetOutputDataType()
{
  //comment
  return AliHLTPHOSDefinitions::fgkAliHLTSandboxDataType;
}


void 
AliHLTPHOSSandboxComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //comment
  constBase = 30;
  inputMultiplier = 1;
}

/*
int 
AliHLTPHOSSandboxComponent::DoProcessing(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
				    AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
				    vector<AliHLTComponentBlockData>& outputBlocks, AliHLTComponentEventDoneData *& edd)
*/
int 
AliHLTPHOSSandboxComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
				    AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/, //TODO: I think size should be set to zero when returning from this method if not data was written to the output buffer.
				    vector<AliHLTComponentBlockData>& /*outputBlocks*/)
{
   //Do event
  
  
  const AliHLTComponentBlockData* iter = 0; 
  unsigned long ndx; 

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      
      if(iter->fDataType != AliHLTPHOSDefinitions::fgkCellEnergyDataType)
	{
	//  cout << "Warning: data type is not fgkCellEnergyDataType " << endl;
	  continue;
	}
     // cout << "Data type is fgkCellEnergyDataType\n";
      fChannelCounterPtr->CountChannels(reinterpret_cast<AliHLTPHOSRcuCellEnergyDataStruct*>(iter->fPtr));
    }
  
  fEvtCnt++;

  if(fEvtCnt % 10 == 0)
   {
	printf("Event #: %d\n", fEvtCnt);   
	//fChannelCounterPtr->PrintOutOfSyncChannels(fEvtCnt);
    }
 
  if(fEvtCnt % 1000 == 0 && fEvtCnt != 0)
    {
//      fChannelCounterPtr->PrintOutOfSyncChannels(fEvtCnt);
      fChannelCounterPtr->FillHistograms(fEvtCnt);
      fChannelCounterPtr->WriteHistograms("/opt/HLT-public/rundir/channelcount.root");
    }

  //cout << "Doing event... \n";
  return 0;
}


int
AliHLTPHOSSandboxComponent::DoInit(int argc, const char** /*argv*/ )
{
  //Do initialization

  fChannelCounterPtr = new AliHLTPHOSChannelCounter();
  
  for(int i = 0; i < argc; i++)
    {
      /*
       */
    }

  return 0;
}

AliHLTComponent*
AliHLTPHOSSandboxComponent::Spawn()
{
  //comment
  return new AliHLTPHOSSandboxComponent();
}
