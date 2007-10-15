//insert copyright

#include "AliHLTPHOSNoiseMapperComponent.h"
#include "AliHLTPHOSNoiseMapper.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSBaseline.h"
#include "AliHLTPHOSProcessor.h"
#include "TH2I.h"
#include "TFile.h" 
//#include <direct.h>
#include <sys/stat.h>
#include <sys/types.h>

const AliHLTComponentDataType AliHLTPHOSNoiseMapperComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};

AliHLTPHOSNoiseMapperComponent gAliHLTPHOSNoiseMapperComponent;

AliHLTPHOSNoiseMapperComponent::AliHLTPHOSNoiseMapperComponent() :
  AliHLTPHOSProcessor(),
  fNoiseMapperPtr(0),
  fWriteInterval(100),
  fFilename(0),
  fDirectory(0),
  fRunNb(0),
  fRateThreshold(10)	 
{
  
}

AliHLTPHOSNoiseMapperComponent::~AliHLTPHOSNoiseMapperComponent()
{
}

int 
AliHLTPHOSNoiseMapperComponent::Deinit()
{
 
  char filename [50];
  sprintf(filename, "%s/run%d_noisemap.root", fDirectory, fRunNb);
      
  FillHistograms();
      
  cout << "Writing file...";
  TFile *outfile = new TFile(filename,"recreate");
  fNoiseCountLowGainHistPtr->Write();
  fNoiseCountHighGainHistPtr->Write();
  fNoiseMapLowGainHistPtr->Write();     
  fNoiseMapHighGainHistPtr->Write();
  delete outfile;
  outfile = 0;
  cout << "Done!\n";
  
  if(fNoiseMapperPtr)
    {
      delete fNoiseMapperPtr;
      fNoiseMapperPtr = 0;
    }
  if(fFilename)
    {
      delete fFilename;
      fFilename = 0;
    }
}

const char*
AliHLTPHOSNoiseMapperComponent::GetComponentID()
{
  return "PhosNoiseMapper";
}

void

AliHLTPHOSNoiseMapperComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{ 
 //Get datatypes for input
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType); 
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSNoiseMapperComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkAliHLTNoiseMapDataType;
}


void 
AliHLTPHOSNoiseMapperComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  constBase = 30;
  inputMultiplier = 1;
}

int 
AliHLTPHOSNoiseMapperComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
					AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
					std::vector<AliHLTComponentBlockData>& outputBlocks)
{
   //Do event
     
  UInt_t tSize            = 0;
  UInt_t offset           = 0; 
  UInt_t mysize           = 0;
  Int_t nRecPoints        = 0;
  Int_t index             = 0;
  
  Int_t fileCount = 0;
  Int_t digitCount = 0;
  char filename [50];


  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = 0; 
  unsigned long ndx; 

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      
      if(iter->fDataType != AliHLTPHOSDefinitions::fgkCellEnergyDataType)
	{
	  //	  cout << "Warning: data type is not fgkCellEnergyDataType " << endl;
	  continue;

	}
      fNoiseMapperPtr->MapNoisyChannels(reinterpret_cast<AliHLTPHOSDigitContainerDataStruct*>(iter->fPtr));
    }
  
  fEventCount++;

  //PushBack(fDigitArrayPtr, kAliHLTAnyDataType, (AliHLTUInt32_t)0);
  
   if(fEventCount % 10 == 0)
    {
      cout << "Event #: " << fEventCount << endl;
    }

  if(fEventCount % fWriteInterval == 0)
    {
      char filename [50];
      sprintf(filename, "%s/run%d_noisemap.root", fDirectory, fRunNb);
      
      FillHistograms();
      
      cout << "Writing file...";
      TFile *outfile = new TFile(filename,"recreate");
      fNoiseCountLowGainHistPtr->Write();
      fNoiseCountLowGainHistPtr->Write();
      fNoiseMapLowGainHistPtr->Write(); 
      fNoiseMapHighGainHistPtr->Write();
      delete outfile;
      outfile = 0;
      cout << "Done!\n";	  
    }
  
  return 0;
}


int
AliHLTPHOSNoiseMapperComponent::DoInit(int argc, const char** argv )
{
  
  Bool_t pathSet = false;
  Bool_t nSamplesSet = false;
  
  fDirectory = new char[50];

  fNoiseMapperPtr = new AliHLTPHOSNoiseMapper();
  
  fNoiseCountLowGainHistPtr = new TH2I("noiseCountLowGain", "Noise count for low gain channels", 64, 0, 63, 56, 0, 55);
  fNoiseCountHighGainHistPtr = new TH2I("noiseCountHighGain", "Noise count for high gain channels", 64, 0, 63, 56, 0, 55);
  fNoiseMapLowGainHistPtr = new TH2I("noiseMapLowGain", "Noise map for low gain channels", 64, 0, 63, 56, 0, 55);
  fNoiseMapHighGainHistPtr = new TH2I("noiseMapHighGain", "Noise map for high gain channels", 64, 0, 63, 56, 0, 55);
  
  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-path", argv[i]))
	{
	  strcpy(fDirectory, argv[i+1]);
	  pathSet = true;
	}
      if(!strcmp("-noisethreshold", argv[i]))
	{
	  fNoiseMapperPtr->SetNoiseThreshold(atof(argv[i+1]));
	}
    }
	
  
  fWriteInterval = 100;

  return 0;
}

AliHLTComponent*
AliHLTPHOSNoiseMapperComponent::Spawn()
{
  return new AliHLTPHOSNoiseMapperComponent();
}

void 
AliHLTPHOSNoiseMapperComponent::FillHistograms()
{
  Int_t channelArray[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];
  fNoiseMapperPtr->GetChannelArray(channelArray);
  for(Int_t x = 0; x < N_XCOLUMNS_MOD; x++)
    {
      for(Int_t z = 0; z < N_ZROWS_MOD; z++)
	{
	  fNoiseCountLowGainHistPtr->SetBinContent(x, z, channelArray[x][z][0]);
	  fNoiseCountHighGainHistPtr->SetBinContent(x, z, channelArray[x][z][1]);
	  if((channelArray[x][z][0]/fEventCount) > fRateThreshold)
	  {
	    fNoiseMapLowGainHistPtr->SetBinContent(x, z, 10);
	  }
	  if((channelArray[x][z][1]/fEventCount) > fRateThreshold)
	  {
	    fNoiseMapHighGainHistPtr->SetBinContent(x, z, 10);
	  }
	}
    }
}
