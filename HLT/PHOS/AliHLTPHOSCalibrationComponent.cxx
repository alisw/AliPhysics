// $Id$


#include "AliHLTPHOSCalibrationComponent.h"

// #include "AliHLTPHOSConstant.h"

#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSEmcCalibData.h"

using namespace std;

AliHLTPHOSCalibrationComponent gAliHLTPHOSCalibrationComponent;

AliHLTPHOSCalibrationComponent::AliHLTPHOSCalibrationComponent() :
  fEmcCalibData(0)
{
}

AliHLTPHOSCalibrationComponent::~AliHLTPHOSCalibrationComponent() 
{
}

AliHLTPHOSCalibrationComponent::AliHLTPHOSCalibrationComponent(const AliHLTPHOSCalibrationComponent&) :
  fEmcCalibData(0)
{
  HLTFatal("copy constructor untested");
}

AliHLTPHOSCalibrationComponent& AliHLTPHOSCalibrationComponent::operator=(const AliHLTPHOSCalibrationComponent&)
{
  HLTFatal("assignement operator untested");
}

const char* AliHLTPHOSCalibrationComponent::GetComponentID()
{
  return "PhosCalibrationComponent";
}

void AliHLTPHOSCalibrationComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkCellEnergyHistogramDataType);
  list.push_back(AliHLTPHOSDefinitions::fgkCellAverageEnergyDataType);
  list.push_back(AliHLTPHOSDefinitions::fgkCellAccumulatedEnergyDataType);
  list.push_back(AliHLTPHOSDefinitions::fgkCellTimingHistogramDataType);
  list.push_back(AliHLTPHOSDefinitions::fgkCellTimingAverageDataType);
}

AliHLTComponentDataType AliHLTPHOSCalibrationComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkEmcCalibDataType;
}
                                     
void AliHLTPHOSCalibrationComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  constBase = 0;
  inputMultiplier = 2;
}

AliHLTComponent* AliHLTPHOSCalibrationComponent::Spawn()
{
  return new AliHLTPHOSCalibrationComponent();
}


Int_t AliHLTPHOSCalibrationComponent::ScanArgument( Int_t argc, const char** argv)
{
  return 0;
}

Int_t AliHLTPHOSCalibrationComponent::InitCalibration()
{
  fEmcCalibData = new AliHLTPHOSEmcCalibData();

  return 0;
}

Int_t AliHLTPHOSCalibrationComponent::DeinitCalibration()
{
  return 0;
}

Int_t AliHLTPHOSCalibrationComponent::ProcessCalibration(const AliHLTComponent_EventData& evtData,
							 const AliHLTComponent_BlockData* blocks,
							 AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr,
							 AliHLTUInt32_t& size,
							 vector<AliHLTComponent_BlockData>& outputBlocks)
//Int_t AliHLTPHOSCalibrationComponent::ProcessCalibration(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData)
{

  const AliHLTComponentBlockData* iter = 0;

  for ( int ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks + ndx;
      if(iter->fDataType == AliHLTPHOSDefinitions::fgkCellEnergyHistogramDataType)
	{
	  for(Int_t mod = 0; mod < N_MODULES; mod++)
	    {
	      for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
		{
		  for(Int_t z = 0; z < N_ZROWS_MOD; z++)
		    {
		      for(Int_t gain = 0; gain < N_GAINS; gain++)
			{
			  fEmcCalibData->SetADCchannelEnergy(mod, x, z, gain, 1);
			  fEmcCalibData->SetADCpedestalEmcMeasured(mod, x, z, gain, 1);
			}
		    }
		}
	    }
	}
      if(iter->fDataType == AliHLTPHOSDefinitions::fgkCellAverageEnergyDataType)
	{

	}
      if(iter->fDataType == AliHLTPHOSDefinitions::fgkCellAccumulatedEnergyDataType)
	{

	}
      if(iter->fDataType == AliHLTPHOSDefinitions::fgkCellTimingHistogramDataType)
	{

	}
      if(iter->fDataType == AliHLTPHOSDefinitions::fgkCellTimingAverageDataType)
	{

	}
      else
	{

	}
    }

  PushBack((TObject*) fEmcCalibData,  AliHLTPHOSDefinitions::fgkEmcCalibDataType, 0);
 
  return 0; 
}
  
Int_t AliHLTPHOSCalibrationComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation
 
  // ** PushBack data to FXS ...
  PushToFXS( (TObject*) fEmcCalibData, "PHOS", "EmcCalibData" ) ;
  
  return 0;
} 
