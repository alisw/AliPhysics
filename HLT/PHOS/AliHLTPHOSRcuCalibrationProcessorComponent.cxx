// $Id$


#include "AliHLTPHOSRcuCalibrationProcessorComponent.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSValidCellDataStruct.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSEmcCalibData.h"
#include "AliHLTPHOSRcuCalibrationProcessor.h"
#include "AliHLTPHOSSharedMemoryInterface.h"
#include "TObjArray.h"

using namespace std;

AliHLTPHOSRcuCalibrationProcessorComponent gAliHLTPHOSRcuCalibrationProcessorComponent;

AliHLTPHOSRcuCalibrationProcessorComponent::AliHLTPHOSRcuCalibrationProcessorComponent() :
  AliHLTCalibrationProcessor(),
  fCalibDataPtr(0),
  fRcuCalibProcessorPtr(0),
  fShmPtr(0)
{

}


AliHLTPHOSRcuCalibrationProcessorComponent::~AliHLTPHOSRcuCalibrationProcessorComponent() 
{

}


/*
AliHLTPHOSRcuCalibrationProcessorComponent::AliHLTPHOSRcuCalibrationProcessorComponent(const AliHLTPHOSRcuCalibrationProcessorComponent&) :
  AliHLTCalibrationProcessor(),
  fCalibDataPtr(0),
  fRcuCalibProcessorPtr(0),
  fShmPtr(0)
{
  HLTFatal("copy constructor untested");
}
*/

 /*
AliHLTPHOSRcuCalibrationProcessorComponent& AliHLTPHOSRcuCalibrationProcessorComponent::operator=(const AliHLTPHOSRcuCalibrationProcessorComponent&)
{
  HLTFatal("assignement operator untested");
}

const char* AliHLTPHOSRcuCalibrationProcessorComponent::GetComponentID()
{
  return "PhosRcuCalibrationProcessor";
}
*/

void AliHLTPHOSRcuCalibrationProcessorComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkCellEnergyDataType);
  /*
  list.push_back(AliHLTPHOSDefinitions::fgkCellEnergyHistogramDataType);
  list.push_back(AliHLTPHOSDefinitions::fgkCellAverageEnergyDataType);
  list.push_back(AliHLTPHOSDefinitions::fgkCellAccumulatedEnergyDataType);
  list.push_back(AliHLTPHOSDefinitions::fgkCellTimingHistogramDataType);
  list.push_back(AliHLTPHOSDefinitions::fgkCellTimingAverageDataType);
  */
}


AliHLTComponentDataType AliHLTPHOSRcuCalibrationProcessorComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkEmcCalibDataType;
}

                                     
void AliHLTPHOSRcuCalibrationProcessorComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  constBase = 0;
  inputMultiplier = 2;
}


AliHLTComponent* 
AliHLTPHOSRcuCalibrationProcessorComponent::Spawn()
{
  return new AliHLTPHOSRcuCalibrationProcessorComponent();
}


const char* 
AliHLTPHOSRcuCalibrationProcessorComponent::GetComponentID()
{
  return "PhosCalibrationProcessor";
}


Int_t 
AliHLTPHOSRcuCalibrationProcessorComponent::ScanArgument( Int_t argc, const char** argv)
{
  const char **c = argv;
  Int_t t= argc;
  c++;
  t++;

  return 0;
}

Int_t AliHLTPHOSRcuCalibrationProcessorComponent::InitCalibration()
{
  fCalibDataPtr = new TObjArray();
  fRcuCalibProcessorPtr = new AliHLTPHOSRcuCalibrationProcessor(2, 0, 0);
  fShmPtr = new AliHLTPHOSSharedMemoryInterface();
  return 0;
}

Int_t AliHLTPHOSRcuCalibrationProcessorComponent::DeinitCalibration()
{
  return 0;
}

Int_t AliHLTPHOSRcuCalibrationProcessorComponent::ProcessCalibration(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData)
{
  const  AliHLTComponentEventData eDta  = evtData;
  AliHLTComponentTriggerData  tDta =  trigData;

  UInt_t specification = 0;
  const AliHLTComponentBlockData* iter = 0;
  iter = GetFirstInputBlock( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC);
  AliHLTPHOSRcuCellEnergyDataStruct* cellDataPtr = 0;
  AliHLTPHOSValidCellDataStruct* currentChannel = 0;
  int totalSamples = 1;

  while(iter != 0)
    {
      specification = specification|iter->fSpecification;
      
      cellDataPtr = (AliHLTPHOSRcuCellEnergyDataStruct*)( iter->fPtr);
      fShmPtr->SetMemory(cellDataPtr);
      currentChannel = fShmPtr->NextChannel();
      
      Int_t* tmpDataPtr = 0;
      Int_t nSamples = 0;

       while(currentChannel != 0)
	{
	  tmpDataPtr = fShmPtr->GetRawData(nSamples);
	  fRcuCalibProcessorPtr->FillEnergy(currentChannel->fX, currentChannel->fZ, currentChannel->fGain, currentChannel->fEnergy);
	  fRcuCalibProcessorPtr->FillLiveChannels(tmpDataPtr, totalSamples, currentChannel->fX, currentChannel->fZ,currentChannel->fGain);
	  currentChannel = fShmPtr->NextChannel();
	}

      iter = GetNextInputBlock(); 
    }

  PushBack((TObject*) fCalibDataPtr,  AliHLTPHOSDefinitions::fgkEmcCalibDataType, specification);
 
  return 0; 
}
  
Int_t AliHLTPHOSRcuCalibrationProcessorComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation
 
  // ** PushBack data to FXS ...
  PushToFXS( (TObject*) fCalibDataPtr, "PHOS", "calibHistoHLT" ) ;
  
  return 0;
} 
