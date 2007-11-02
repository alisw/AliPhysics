
//insert copyright

#include "AliHLTPHOSDigitMakerComponent.h"
#include "AliHLTPHOSDigitMaker.h"
#include "TTree.h"
#include "AliHLTPHOSProcessor.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"
#include "TClonesArray.h"
#include "TFile.h"
#include <sys/stat.h>
#include <sys/types.h>

const AliHLTComponentDataType AliHLTPHOSDigitMakerComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};

AliHLTPHOSDigitMakerComponent gAliHLTPHOSDigitMakerComponent;

AliHLTPHOSDigitMakerComponent::AliHLTPHOSDigitMakerComponent() :
  AliHLTPHOSProcessor(),
  fDigitMakerPtr(0),
  fEventCount(0)
{
 
}

AliHLTPHOSDigitMakerComponent::~AliHLTPHOSDigitMakerComponent()
{
}

int 
AliHLTPHOSDigitMakerComponent::Deinit()
{ 
  if(fDigitMakerPtr)
    {
      delete fDigitMakerPtr;
      fDigitMakerPtr = 0;
    }
  return 0;
}

const char*
AliHLTPHOSDigitMakerComponent::GetComponentID()
{
  return "PhosDigitMaker";
}

void

AliHLTPHOSDigitMakerComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{ 
 //Get datatypes for input
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType); 
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSDigitMakerComponent::GetOutputDataType()
{
  return AliHLTPHOSDefinitions::fgkAliHLTDigitDataType;
}


void 
AliHLTPHOSDigitMakerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  constBase = 30;
  inputMultiplier = 1;
}

int 
AliHLTPHOSDigitMakerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
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
  fDigitContainerPtr = (AliHLTPHOSDigitContainerDataStruct*)outBPtr;
  //fDigitMakerPtr->SetDigitContainerStruct(fDigitContainerPtr);
  fDigitMakerPtr->SetDigitContainerStruct((AliHLTPHOSDigitContainerDataStruct*)outBPtr);

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      
      if(iter->fDataType != AliHLTPHOSDefinitions::fgkCellEnergyDataType)
	{
	  //	  cout << "Warning: data type is not fgkCellEnergyDataType " << endl;
	  continue;

	}
      digitCount = fDigitMakerPtr->MakeDigits(reinterpret_cast<AliHLTPHOSRcuCellEnergyDataStruct*>(iter->fPtr));
    }
  fEventCount++;
  
  mysize = 0;
  offset = tSize;
      
  mysize += sizeof(AliHLTPHOSDigitContainerDataStruct);
  ((AliHLTPHOSDigitContainerDataStruct*)outBPtr)->fNDigits = digitCount;
  AliHLTComponentBlockData bd;
  FillBlockData( bd );
  bd.fOffset = offset;
  bd.fSize = mysize;
  bd.fDataType = AliHLTPHOSDefinitions::fgkAliHLTDigitDataType;
  bd.fSpecification = 0xFFFFFFFF;
  outputBlocks.push_back( bd );
       
  tSize += mysize;
  outBPtr += mysize;
      
  if( tSize > size )
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSDigitMakerComponent::DoEvent", "Too much data",
	       "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
	       , tSize, size );
      return EMSGSIZE;
    }
      
  fDigitMakerPtr->Reset();
  
  if(fEventCount % 10 == 0)
    {
      cout << "Event #: " << fEventCount << endl;
      cout << "  - Number of digits found: " << digitCount << endl;
    }
  
  return 0;
}


int
AliHLTPHOSDigitMakerComponent::DoInit(int argc, const char** argv )
{
  //Do initialization

  fDigitMakerPtr = new AliHLTPHOSDigitMaker();
  
  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-threshold", argv[i]))
	fDigitMakerPtr->SetDigitThreshold(atoi(argv[i+1]));
      if(!strcmp("-presamples", argv[i]))
	fDigitMakerPtr->SetNrPresamples(atoi(argv[i+1]));
    }
 
  //fDigitMakerPtr->SetDigitThreshold(2);

  return 0;
}

AliHLTComponent*
AliHLTPHOSDigitMakerComponent::Spawn()
{
  return new AliHLTPHOSDigitMakerComponent();
}
