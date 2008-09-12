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

#include <iostream>

#include "AliHLTPHOSClusterizerComponent.h"
#include "AliHLTPHOSClusterizer.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"



/** @file   AliHLTPHOSClusterizerComponent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  A clusterizer component for PHOS HLT
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

const AliHLTComponentDataType AliHLTPHOSClusterizerComponent::fgkInputDataTypes[]=
  {
    kAliHLTVoidDataType,{0,"",""}
  };

AliHLTPHOSClusterizerComponent gAliHLTPHOSClusterizerComponent;


AliHLTPHOSClusterizerComponent::AliHLTPHOSClusterizerComponent(): 
  AliHLTPHOSProcessor(), 
  fAllDigitsPtr(0),
  fClusterizerPtr(0),
  fRecPointStructArrayPtr(0),
  fDigitCount(0),
  fModuleClusterizationMode(false),
  fNoCrazyness(0)
{
  //See headerfile for documentation
}

AliHLTPHOSClusterizerComponent::~AliHLTPHOSClusterizerComponent()
{
  //See headerfile for documentation

  if(fClusterizerPtr)
    {
      delete fClusterizerPtr;
      fClusterizerPtr = 0;
    }
  if(fRecPointStructArrayPtr)
    {
      for (int i = 0; i < 1000; i++)
        {
          //	  fRecPointStructArrayPtr[i].Del();
        }
      delete fRecPointStructArrayPtr;
      fRecPointStructArrayPtr = 0;
    }
  if(fAllDigitsPtr)
    {
      delete fAllDigitsPtr;
      fAllDigitsPtr = 0;
    }
}


int
AliHLTPHOSClusterizerComponent::Deinit()
{
  //See headerfile for documentation

  if (fClusterizerPtr)
    {
      delete fClusterizerPtr;
      fClusterizerPtr = 0;
    }
  for (int i = 0; i < 1000; i++)
    {
      //    fRecPointStructArrayPtr[i].Del();
    }

  if (fRecPointStructArrayPtr)
    {
      for (int i = 0; i < 1000; i++)
        {
          //	  fRecPointStructArrayPtr[i].Del();
        }
      delete fRecPointStructArrayPtr;
      fRecPointStructArrayPtr = 0;
    }

  return 0;
}

const Char_t*
AliHLTPHOSClusterizerComponent::GetComponentID()
{
  //See headerfile for documentation

  return "PhosClusterizer";
}

void
AliHLTPHOSClusterizerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkDigitDataType);
}

AliHLTComponentDataType
AliHLTPHOSClusterizerComponent::GetOutputDataType()
{
  //See headerfile for documentation
  return AliHLTPHOSDefinitions::fgkClusterDataType;
}

void
AliHLTPHOSClusterizerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation
  constBase = sizeof(AliHLTPHOSRecPointContainerStruct);
  inputMultiplier = 1;
}

int
AliHLTPHOSClusterizerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
                                        AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
                                        std::vector<AliHLTComponentBlockData>& outputBlocks)
{
  //See headerfile for documentation
  
  UInt_t tSize            = 0;
  UInt_t offset           = 0;
  UInt_t mysize           = 0;
  Int_t nRecPoints        = 0;
  Int_t j = 0;

  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = 0;
  unsigned long ndx;

  UInt_t specification = 0;

  AliHLTPHOSDigitContainerDataStruct *digitContainerPtr = 0;
  fClusterizerPtr->SetRecPointContainer((AliHLTPHOSRecPointContainerStruct*)outBPtr);
  
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if (iter->fDataType == AliHLTPHOSDefinitions::fgkDigitDataType)
        {
	  specification = specification|iter->fSpecification;
	  digitContainerPtr = reinterpret_cast<AliHLTPHOSDigitContainerDataStruct*>(iter->fPtr);
	  if(fModuleClusterizationMode)
	    {
	      for (UInt_t i = 0; i < digitContainerPtr->fNDigits; i++)
		{
		  fAllDigitsPtr->fDigitDataStruct[j].fX = digitContainerPtr->fDigitDataStruct[i].fX;
		  fAllDigitsPtr->fDigitDataStruct[j].fZ = digitContainerPtr->fDigitDataStruct[i].fZ;
		  fAllDigitsPtr->fDigitDataStruct[j].fAmplitude = digitContainerPtr->fDigitDataStruct[i].fAmplitude;
		  fAllDigitsPtr->fDigitDataStruct[j].fEnergy = digitContainerPtr->fDigitDataStruct[i].fEnergy;
		  fAllDigitsPtr->fDigitDataStruct[j].fTime = digitContainerPtr->fDigitDataStruct[i].fTime;
		  fAllDigitsPtr->fDigitDataStruct[j].fCrazyness = digitContainerPtr->fDigitDataStruct[i].fCrazyness;
		  j++;
		}
	  
	    }
	}
    }
  fAllDigitsPtr->fNDigits = j;
  if(fModuleClusterizationMode && fAllDigitsPtr != 0)
    {
      fClusterizerPtr->SetDigitContainer(fAllDigitsPtr);
    }
  else if(digitContainerPtr != 0)
    {
      fClusterizerPtr->SetDigitContainer(digitContainerPtr);
    }      
  if(digitContainerPtr != 0)
    {
      nRecPoints = fClusterizerPtr->ClusterizeEvent();
      mysize = 0;
      offset = tSize;

      mysize += sizeof(AliHLTPHOSRecPointContainerStruct);
      
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;
      bd.fDataType = AliHLTPHOSDefinitions::fgkClusterDataType;
      bd.fSpecification = specification;
      outputBlocks.push_back( bd );
      tSize += mysize;
      outBPtr += mysize;
  
    }
  
  if ( tSize > size )
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSClusterizerComponent::DoEvent", "Too much data",
               "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
               , tSize, size );
      return EMSGSIZE;
    }

  size = mysize;
  
  return 0;
}

int
AliHLTPHOSClusterizerComponent::DoInit(int argc, const char** argv )
{
  //See headerfile for documentation

  fAllDigitsPtr = new AliHLTPHOSDigitContainerDataStruct();
  fClusterizerPtr = new AliHLTPHOSClusterizer();
  fClusterizerPtr->SetDigitContainer(fAllDigitsPtr);
  fNoCrazyness = false;
  //

  ScanArgumentsModule(argc, argv);
  for (int i = 0; i < argc; i++)
    {
      if(!strcmp("-digitthreshold", argv[i]))
	{
	  fClusterizerPtr->SetEmcMinEnergyThreshold(atof(argv[i+1]));
	}
      if(!strcmp("-recpointthreshold", argv[i]))
	{
	  fClusterizerPtr->SetEmcClusteringThreshold(atof(argv[i+1]));
	}
      if(!strcmp("-modulemode", argv[i]))
	{
	  fModuleClusterizationMode = true;
	}
    }

  return 0;
}

AliHLTComponent*
AliHLTPHOSClusterizerComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSClusterizerComponent();
}
