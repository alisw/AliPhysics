// $Id: AliHLTCaloClusterizerComponent.cxx 36709 2009-11-12 16:57:55Z odjuvsla $

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

#include "AliHLTCaloClusterizerComponent.h"
#include "AliHLTCaloClusterizer.h"
#include "AliHLTCaloRecPointDataStruct.h"
#include "AliHLTCaloRecPointHeaderStruct.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloDigitContainerDataStruct.h"
#include "AliHLTCaloDefinitions.h"

/** @file   AliHLTCaloClusterizerComponent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  A clusterizer component for PHOS HLT
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

AliHLTCaloClusterizerComponent gAliHLTCaloClusterizerComponent;

AliHLTCaloClusterizerComponent::AliHLTCaloClusterizerComponent(): 
  AliHLTPHOSProcessor(), 
  fAllDigitsPtr(0),
  fClusterizerPtr(0),
  fDigitCount(0),
  fNoCrazyness(0)
{
  //See headerfile for documentation
}

AliHLTCaloClusterizerComponent::~AliHLTCaloClusterizerComponent()
{
  //See headerfile for documentation

  if(fClusterizerPtr)
    {
      delete fClusterizerPtr;
      fClusterizerPtr = 0;
    }
  if(fAllDigitsPtr)
    {
      delete fAllDigitsPtr;
      fAllDigitsPtr = 0;
    }
}


int
AliHLTCaloClusterizerComponent::Deinit()
{
  //See headerfile for documentation

  if (fClusterizerPtr)
    {
      delete fClusterizerPtr;
      fClusterizerPtr = 0;
    }

  return 0;
}

void
AliHLTCaloClusterizerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //See headerfile for documentation
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkDigitDataType);
}

AliHLTComponentDataType
AliHLTCaloClusterizerComponent::GetOutputDataType()
{
  //See headerfile for documentation
  return AliHLTPHOSDefinitions::fgkRecPointDataType;
}

void
AliHLTCaloClusterizerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation
  constBase = sizeof(AliHLTCaloRecPointHeaderStruct) + sizeof(AliHLTPHOSRecPointDataStruct) + (sizeof(AliHLTCaloDigitDataStruct) << 7); //Reasonable estimate... ;
  inputMultiplier = 1.5;
}

int
AliHLTCaloClusterizerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
                                        AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
                                        std::vector<AliHLTComponentBlockData>& outputBlocks)
{
  //See headerfile for documentation

  if(blocks == 0) return 0;
  
  UInt_t offset           = 0;
  UInt_t mysize           = 0;
  Int_t nRecPoints        = 0;
  Int_t nDigits           = 0;

  UInt_t availableSize = size;
  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = 0;
  unsigned long ndx;
  
  UInt_t specification = 0;
  
  AliHLTCaloDigitDataStruct *digitDataPtr = 0;

  AliHLTCaloRecPointHeaderStruct* recPointHeaderPtr = reinterpret_cast<AliHLTCaloRecPointHeaderStruct*>(outBPtr);

  fClusterizerPtr->SetRecPointDataPtr(reinterpret_cast<AliHLTCaloRecPointDataStruct*>(outBPtr+sizeof(AliHLTCaloRecPointHeaderStruct)));

  // Adding together all the digits, should be put in standalone method  
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if (iter->fDataType == AliHLTCaloDefinitions::fgkDigitDataType)
	{
	  // Get the digit header

	  // Update the number of digits
	  nDigits += iter->fSize/sizeof(AliHLTCaloDigitDataStruct);;
	  // Get the specification
	  specification = specification|iter->fSpecification;

	  digitDataPtr = reinterpret_cast<AliHLTCaloDigitDataStruct*>(iter->fPtr);
	  for (Int_t i = 0; i < nDigits; i++)
	    {
	      fAllDigitsPtr->fDigitDataStruct[j].fX = digitDataPtr->fX;
	      fAllDigitsPtr->fDigitDataStruct[j].fZ = digitDataPtr->fZ;
	      fAllDigitsPtr->fDigitDataStruct[j].fEnergy = digitDataPtr->fEnergy;
	      //  HLTDebug("Digit energy: %f", digitDataPtr->fEnergy);
	      fAllDigitsPtr->fDigitDataStruct[j].fTime = digitDataPtr->fTime;
	      fAllDigitsPtr->fDigitDataStruct[j].fCrazyness = digitDataPtr->fCrazyness;
	      j++;
	      digitDataPtr++;
	    }

	}
    }
  
  fAllDigitsPtr->fNDigits = j;
  nRecPoints = fClusterizerPtr->ClusterizeEvent(size, mysize);

  if(nRecPoints == -1)
    {
      HLTError("Running out of buffer, exiting for safety.");
      return -ENOBUFS;
    }

  recPointHeaderPtr->fNRecPoints = nRecPoints;
  mysize += sizeof(AliHLTCaloRecPointHeaderStruct);
  
  HLTDebug("Number of clusters: %d", nRecPoints);

  AliHLTComponentBlockData bd;
  FillBlockData( bd );
  bd.fOffset = offset;
  bd.fSize = mysize;
  bd.fDataType = AliHLTPHOSDefinitions::fgkClusterDataType;
  bd.fSpecification = specification;
  outputBlocks.push_back( bd );
     
  size = mysize;
  
  return 0;
}

int 
AliHLTCaloClusterizerComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{  
  // see header file for class documentation

  const char* path="HLT/ConfigPHOS/ClusterizerComponent";

  if (cdbEntry) path = cdbEntry;

  return ConfigureFromCDBTObjString(cdbEntry);
}

int 
AliHLTCaloClusterizerComponent::ScanConfigurationArgument(int argc, const char **argv)
{
  //See header file for documentation

  if(argc <= 0) return 0;

  int i=0;

  TString argument=argv[i];

  if (argument.CompareTo("-digitthreshold") == 0)
    {
      if (++i >= argc) return -EPROTO;
      argument = argv[i];
      fClusterizerPtr->SetEmcMinEnergyThreshold(argument.Atof());
      return 1;
    }

  if (argument.CompareTo("-recpointthreshold") == 0)
    {
      if (++i >= argc) return -EPROTO;
      argument = argv[i];
      fClusterizerPtr->SetEmcClusteringThreshold(argument.Atof());
      return 1;
    }
  return 0;
}

int
AliHLTCaloClusterizerComponent::DoInit(int argc, const char** argv )
{
  //See headerfile for documentation

  fAllDigitsPtr = new AliHLTCaloDigitContainerDataStruct();
  fClusterizerPtr = new AliHLTCaloClusterizer();
  fClusterizerPtr->SetDigitContainer(fAllDigitsPtr);
  fNoCrazyness = false;
  //

  //  const char *path = "HLT/ConfigPHOS/ClusterizerComponent";

  //  ConfigureFromCDBTObjString(path);

  for (int i = 0; i < argc; i++)
    {
      ScanConfigurationArgument(i, argv);
    }

  return 0;
}

AliHLTComponent*
AliHLTCaloClusterizerComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTCaloClusterizerComponent();
}
