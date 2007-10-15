//insert copyright

#include "AliHLTPHOSSandboxComponent.h"
#include "AliHLTPHOSProcessor.h"
#include "AliHLTPHOSChannelCounter.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSDefinitions.h"

const AliHLTComponentDataType AliHLTPHOSSandboxComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};

AliHLTPHOSSandboxComponent gAliHLTPHOSSandboxComponent;

AliHLTPHOSSandboxComponent::AliHLTPHOSSandboxComponent() :
  AliHLTPHOSProcessor(),
  fEventCount(0),
  fChannelCounterPtr(0)
{

}

AliHLTPHOSSandboxComponent::~AliHLTPHOSSandboxComponent()
{
}
/*
AliHLTPHOSSandboxComponent::AliHLTPHOSSandboxComponent(const AliHLTPHOSSandboxComponent &) :
  AliHLTDataSink(),
  fEventCount(0),
  fChannelCounterPtr(0)
{

  //Copy constructor, not implemented

}
*/

int 
AliHLTPHOSSandboxComponent::Deinit()
{
  printf("AliHLTPHOSSandboxComponent::Deinit()\n");

  //fChannelCounterPtr->PrintOutOfSyncChannels(fEventCount);
  fChannelCounterPtr->FillHistograms(fEventCount);
  fChannelCounterPtr->WriteHistograms("/opt/HLT-public/rundir/channelcount.root");
  return 0;
}

Int_t
AliHLTPHOSSandboxComponent::DoDeinit()
{
  //Deinitialize the component
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSSandboxComponent DoDeinit");

  return 0;
}



const char*
AliHLTPHOSSandboxComponent::GetComponentID()
{
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
  return AliHLTPHOSDefinitions::fgkAliHLTSandboxDataType;
}


void 
AliHLTPHOSSandboxComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
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
				    AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
				    vector<AliHLTComponentBlockData>& outputBlocks)
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
  
  fEventCount++;

  if(fEventCount % 10 == 0)
   {
	printf("Event #: %d\n", fEventCount);   
	//fChannelCounterPtr->PrintOutOfSyncChannels(fEventCount);
    }
 
  if(fEventCount % 1000 == 0 && fEventCount != 0)
    {
//      fChannelCounterPtr->PrintOutOfSyncChannels(fEventCount);
      fChannelCounterPtr->FillHistograms(fEventCount);
      fChannelCounterPtr->WriteHistograms("/opt/HLT-public/rundir/channelcount.root");
    }

  //cout << "Doing event... \n";
  return 0;
}


int
AliHLTPHOSSandboxComponent::DoInit(int argc, const char** argv )
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
  return new AliHLTPHOSSandboxComponent();
}
