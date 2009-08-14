// $Id$

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
#include "TFile.h"
#include "AliHLTPHOSHistogramProducer.h"
#include "AliHLTPHOSPhysicsHistogramProducer.h"
#include "AliHLTPHOSCaloClusterContainerStruct.h"
#include "TClonesArray.h"
/** 
 * @file   AliHLTPHOSHistogramProducerComponent.cxx
 * @author Oystein Djuvsland
 * @date   
 * @brief  A histogram producer component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


AliHLTPHOSHistogramProducerComponent gAliHLTPHOSHistogramProducerComponent;

AliHLTPHOSHistogramProducerComponent::AliHLTPHOSHistogramProducerComponent() :
  AliHLTPHOSProcessor(),
  fPhysicsHistogramProducerPtr(0),
  fPushModulo(1)
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
  if(fPhysicsHistogramProducerPtr != 0)
    {
      delete fPhysicsHistogramProducerPtr;
      fPhysicsHistogramProducerPtr = 0;
    }

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
  list.push_back(AliHLTPHOSDefinitions::fgkESDCaloClusterDataType);
  list.push_back(AliHLTPHOSDefinitions::fgkESDCaloCellsDataType);

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

  const AliHLTComponentBlockData* block = GetFirstInputBlock(AliHLTPHOSDefinitions::fgkESDCaloClusterDataType);
  
  if(block != 0)
    {
      HLTError("Number of clusters: %d", reinterpret_cast<TClonesArray*>(const_cast<AliHLTComponentBlockData*>(block))->GetEntries());
      fPhysicsHistogramProducerPtr->AnalyseClusters(reinterpret_cast<TClonesArray*>(const_cast<AliHLTComponentBlockData*>(block)));
    }
  
  //  block = GetFirstInputBlock(AliHLTPHOSDefinitions::fgkESDCaloCellsDataType);
			     

  if(block != 0)
    {
      
    }
  //  if(fEventCount%fPushModulo == 0)
  //  {

  HLTError("WRITING HISTOGRAMS!");
  TFile* file = TFile::Open("/tmp/hist.root", "RECREATE");
  (fPhysicsHistogramProducerPtr->GetHistograms())->Write();
  file->Close();
  PushBack(fPhysicsHistogramProducerPtr->GetHistograms(), AliHLTPHOSDefinitions::fgkPhysicsHistogramsDataType);
      // }
    
  return 0;
}


int
AliHLTPHOSHistogramProducerComponent::DoInit(int argc, const char** argv )
{
  //see header file for documentation

  fPhysicsHistogramProducerPtr = new AliHLTPHOSPhysicsHistogramProducer();
  
  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-pushmodulo", argv[i]))
	{
	
	  fPushModulo = atoi(argv[i+1]);
	}
    }

  return 0;
}

AliHLTComponent*
AliHLTPHOSHistogramProducerComponent::Spawn()
{
  //see header file for documentation
  return new AliHLTPHOSHistogramProducerComponent();
}

