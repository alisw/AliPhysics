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

#include "AliHLTPHOSHistogramProducerComponent.h"
#include "AliHLTPHOSProcessor.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "AliHLTPHOSHistogramProducer.h"
#include "AliHLTPHOSCaloClusterContainerStruct.h"

/** 
 * @file   AliHLTPHOSHistogramProducerComponent.cxx
 * @author Oystein Djuvsland
 * @date   
 * @brief  A digit maker component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


const AliHLTComponentDataType AliHLTPHOSHistogramProducerComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};

AliHLTPHOSHistogramProducerComponent gAliHLTPHOSHistogramProducerComponent;

AliHLTPHOSHistogramProducerComponent::AliHLTPHOSHistogramProducerComponent() :
  AliHLTPHOSProcessor(),
  fClusterEnergiesHistPtr(0),
  fMultiplicitiesHistPtr(0),
  fClusterNtuplePtr(0),
  fDoFillClusterEnergies(false),
  fDoFillMultiplicities(false),
  fDoFillNtuple(false),
  fHistogramProducerPtr(0)
{
  //see header file for documentation
}

AliHLTPHOSHistogramProducerComponent::~AliHLTPHOSHistogramProducerComponent()
{
  //see header file for documentation
}

int 
AliHLTPHOSHistogramProducerComponent::Deinit()
{ 
  //see header file for documentation
  return 0;
}

const char*
AliHLTPHOSHistogramProducerComponent::GetComponentID()
{
  //see header file for documentation
  return "PhosHistoProducer";
}


void
AliHLTPHOSHistogramProducerComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{ 
  //see header file for documentation
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkChannelDataType);

//   const AliHLTComponentDataType* pType=fgkInputDataTypes;
//   while (pType->fID!=0) {
//     list.push_back(*pType); 
//     pType++;
//   }
}

AliHLTComponentDataType 
AliHLTPHOSHistogramProducerComponent::GetOutputDataType()
{
  //see header file for documentation
  return AliHLTPHOSDefinitions::fgkPhosHistDataType;
}


void 
AliHLTPHOSHistogramProducerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //see header file for documentation
  constBase = 30;
  inputMultiplier = 5;
}

// int 
// AliHLTPHOSHistogramProducerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
// 					AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
// 					std::vector<AliHLTComponentBlockData>& outputBlocks)

int 
AliHLTPHOSHistogramProducerComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					      AliHLTComponentTriggerData& /*trigData*/)


{
  //see header file for documentation

  //  UInt_t specification = 0;

  //  AliHLTPHOSCaloClusterContainerStruct* tmpClusters = 0;

  const AliHLTComponentBlockData* block = GetFirstInputBlock(AliHLTPHOSDefinitions::fgkClusterDataType);
  
  while(block != 0)
    {
      fHistogramProducerPtr->Fill(reinterpret_cast<AliHLTPHOSCaloClusterContainerStruct*>(block->fPtr));
      block = GetNextInputBlock();
    }
  
  if(fDoFillClusterEnergies)
    {
      PushBack(fClusterEnergiesHistPtr, AliHLTPHOSDefinitions::fgkPhosHistDataType);
    }
  if(fDoFillMultiplicities)
    {
      PushBack(fMultiplicitiesHistPtr, AliHLTPHOSDefinitions::fgkPhosHistDataType);
    }
  if(fDoFillNtuple)
    {
      PushBack(fClusterNtuplePtr, AliHLTPHOSDefinitions::fgkPhosHistDataType);
    }
    
  return 0;
}


int
AliHLTPHOSHistogramProducerComponent::DoInit(int argc, const char** argv )
{
  //see header file for documentation

  fHistogramProducerPtr = new AliHLTPHOSHistogramProducer();
  
  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-dofillclusterenergies", argv[i]))
	{
	  fHistogramProducerPtr->SetFillClusterEnergies(true);
	  fDoFillClusterEnergies = true;
	}
      if(!strcmp("-dofillmultiplicities", argv[i]))
	{
	  fHistogramProducerPtr->SetFillMultiplicities(true);
	  fDoFillMultiplicities = true;
	}
      if(!strcmp("-dofillntuple", argv[i]))
	{
	  fHistogramProducerPtr->SetFillClusterNtuple(true);
	  fDoFillNtuple = true;
	}
      if(!strcmp("-maxntupleentries", argv[i]))
	{
	  fHistogramProducerPtr->SetMaxNtupleEntries(atoi(argv[i+1]));
	}
    }
 
  fHistogramProducerPtr->InitializeObjects();

  if(fDoFillClusterEnergies)
    {
      fClusterEnergiesHistPtr = fHistogramProducerPtr->GetClusterEnergiesHistogram();
    }
  if(fDoFillMultiplicities)
    {
      fMultiplicitiesHistPtr = fHistogramProducerPtr->GetMultiplicitiesHistogram();
    }
  if(fDoFillNtuple)
    {
      fClusterNtuplePtr = fHistogramProducerPtr->GetClusterNtuple();
    }

  //fDigitMakerPtr->SetDigitThreshold(2);

  return 0;
}

AliHLTComponent*
AliHLTPHOSHistogramProducerComponent::Spawn()
{
  //see header file for documentation
  return new AliHLTPHOSHistogramProducerComponent();
}

