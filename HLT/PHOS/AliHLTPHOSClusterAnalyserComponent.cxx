/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSClusterAnalyserComponent.h"
#include "AliHLTPHOSClusterAnalyser.h"
#include "AliHLTPHOSRecPointContainerStruct.h"
#include "AliHLTPHOSCaloClusterContainerStruct.h"

/** @file   AliHLTPHOSClusterAnalyserComponent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  A cluster analyser component for PHOS HLT
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

const AliHLTComponentDataType AliHLTPHOSClusterAnalyserComponent::fgkInputDataTypes[]=
  {
    kAliHLTVoidDataType,{0,"",""}
  };

AliHLTPHOSClusterAnalyserComponent gAliHLTPHOSClusterAnalyserComponent;


AliHLTPHOSClusterAnalyserComponent::AliHLTPHOSClusterAnalyserComponent(): AliHLTPHOSProcessor(), 
									  fClusterAnalyserPtr(0),
									  fDoDeconvolution(0),
									  fDoCalculateMoments(0)
{
  //See headerfile for documentation
}

AliHLTPHOSClusterAnalyserComponent::~AliHLTPHOSClusterAnalyserComponent()
{
  //See headerfile for documentation

  if (fClusterAnalyserPtr)
    {
      delete fClusterAnalyserPtr;
      fClusterAnalyserPtr = 0;
    }
}

int
AliHLTPHOSClusterAnalyserComponent::Deinit()
{
  //See headerfile for documentation

  if (fClusterAnalyserPtr)
    {
      delete fClusterAnalyserPtr;
      fClusterAnalyserPtr = 0;
    }
  return 0;
}

const Char_t*
AliHLTPHOSClusterAnalyserComponent::GetComponentID()
{
  //See headerfile for documentation

  return "PhosClusterAnalyser";
}

void
AliHLTPHOSClusterAnalyserComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //See headerfile for documentation
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkClusterDataType);
//   const AliHLTComponentDataType* pType=fgkInputDataTypes;
//   while (pType->fID!=0)
//     {
//       list.push_back(*pType);
//       pType++;
//     }
}

AliHLTComponentDataType
AliHLTPHOSClusterAnalyserComponent::GetOutputDataType()
{
  //See headerfile for documentation

  return AliHLTPHOSDefinitions::fgkClusterDataType;
}

void
AliHLTPHOSClusterAnalyserComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation

  constBase = sizeof(AliHLTPHOSCaloClusterContainerStruct);
  inputMultiplier = 1;
}

AliHLTComponent*
AliHLTPHOSClusterAnalyserComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSClusterAnalyserComponent();
}

int
AliHLTPHOSClusterAnalyserComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
                                        AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
                                        std::vector<AliHLTComponentBlockData>& outputBlocks)
{
  //See headerfile for documentation

  UInt_t tSize            = 0;
  UInt_t offset           = 0;
  UInt_t mysize           = 0;
  Int_t nClusters         = 0;

  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = 0;
  unsigned long ndx;

  UInt_t specification = 0;

  AliHLTPHOSRecPointContainerStruct* recPointContainerPtr = 0;
  fClusterAnalyserPtr->SetCaloClusterContainer((AliHLTPHOSCaloClusterContainerStruct*)outBPtr);
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx; 
      if (iter->fDataType != AliHLTPHOSDefinitions::fgkClusterDataType)
        {
          continue;
        }
      specification = specification|iter->fSpecification;
      recPointContainerPtr = reinterpret_cast<AliHLTPHOSRecPointContainerStruct*>(iter->fPtr);
      fClusterAnalyserPtr->SetRecPointContainer(recPointContainerPtr);
      if(fDoDeconvolution)
	{
	  fClusterAnalyserPtr->DeconvoluteClusters();
	}
      fClusterAnalyserPtr->CalculateCenterOfGravity();
      if(fDoCalculateMoments)
	{
	  fClusterAnalyserPtr->CalculateRecPointMoments();
	}
      nClusters = fClusterAnalyserPtr->CreateClusters();
    }

  mysize = 0; 
  offset = tSize;
  mysize += sizeof(AliHLTPHOSCaloClusterContainerStruct);
  
  AliHLTComponentBlockData bd;
  FillBlockData( bd );
  bd.fOffset = offset;
  bd.fSize = mysize;
  bd.fDataType = AliHLTPHOSDefinitions::fgkCaloClusterDataType;
  bd.fSpecification = specification;
  outputBlocks.push_back( bd );

  tSize += mysize;
  outBPtr += mysize;
 
  if ( tSize > size )
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSClusterAnalyserComponent::DoEvent", "Too much data","Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.", tSize, size );
      return EMSGSIZE;
    }

  fPhosEventCount++; 
  if(fPrintInfoModule == kTRUE)
    {
      if(fPhosEventCount%fPrintInfoFrequncyModule == 0)
      	{
	  Logging( kHLTLogInfo, __FILE__ , "Clusters analyzed",  "Has analyzed  %lu  for this event",  nClusters);
	}  
    }

  size = mysize;  
  return 0;

}

int
AliHLTPHOSClusterAnalyserComponent::DoInit(int argc, const char** argv )
{

  //See headerfile for documentation
  
  fClusterAnalyserPtr = new AliHLTPHOSClusterAnalyser();
  ScanArgumentsModule(argc, argv);
  for (int i = 0; i < argc; i++)
    {
      if(!strcmp("-dodeconvolution", argv[i]))
	{
	  fDoDeconvolution = true;
	}
      if(!strcmp("-doclusterfit", argv[i]))
	{
	  fClusterAnalyserPtr->SetDoClusterFit();
	  fDoCalculateMoments = true;
	}
      if(!strcmp("-haveCPV", argv[i]))
	{
	  fClusterAnalyserPtr->SetHaveCPVInfo();
	}
      if(!strcmp("-doPID", argv[i]))
	{
	  fClusterAnalyserPtr->SetDoPID();
	} 
      if(!strcmp("-havedistbadchannel", argv[i]))
	{
	  fClusterAnalyserPtr->SetHaveDistanceToBadChannel();
	}

    }

  return 0;
}
