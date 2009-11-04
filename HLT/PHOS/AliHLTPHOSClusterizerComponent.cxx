// $Id$

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
#include "AliHLTPHOSRecPointHeaderStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
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

AliHLTPHOSClusterizerComponent gAliHLTPHOSClusterizerComponent;

AliHLTPHOSClusterizerComponent::AliHLTPHOSClusterizerComponent(): 
  AliHLTPHOSProcessor(), 
  fAllDigitsPtr(0),
  fClusterizerPtr(0),
  fDigitCount(0),
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
  //See headerfile for documentation
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkDigitDataType);
}

AliHLTComponentDataType
AliHLTPHOSClusterizerComponent::GetOutputDataType()
{
  //See headerfile for documentation
  return AliHLTPHOSDefinitions::fgkRecPointDataType;
}

void
AliHLTPHOSClusterizerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation
  constBase = sizeof(AliHLTPHOSRecPointHeaderStruct) + sizeof(AliHLTPHOSRecPointDataStruct) + (sizeof(AliHLTPHOSDigitDataStruct) << 7); //Reasonable estimate... ;
  inputMultiplier = 1.5;
}

int
AliHLTPHOSClusterizerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
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
  
  AliHLTPHOSDigitHeaderStruct *digitHeaderPtr = 0;
  AliHLTPHOSDigitHeaderStruct *outputDigitHeaderPtr = reinterpret_cast<AliHLTPHOSDigitHeaderStruct*>(outBPtr);

  AliHLTPHOSDigitDataStruct *firstDigitPtr = 0;
  AliHLTPHOSDigitDataStruct *lastDigitPtr = 0;
  
  // Adding together all the digits, should be put in standalone method  
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if (iter->fDataType == AliHLTPHOSDefinitions::fgkDigitDataType)
	{
	  // Update the number of digits
	  nDigits += digitHeaderPtr->fNDigits;

	  // Get the specification
	  specification = specification|iter->fSpecification;
	  
	  // Get the digit header
	  digitHeaderPtr = reinterpret_cast<AliHLTPHOSDigitHeaderStruct*>(iter->fPtr);
	  
	  // Check if we have the first buffer in the event
	  if(!firstDigitPtr)
	    {
	      if(availableSize < digitHeaderPtr->fNDigits*sizeof(AliHLTPHOSDigitDataStruct) + sizeof(AliHLTPHOSDigitHeaderStruct))
		{
		  HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", totSize, size);
		  return -1;
		}
	      // If so, lets copy the header and the corresponding digits to the output
	      memcpy(outBPtr, iter->fPtr, digitHeaderPtr->fNDigits*sizeof(AliHLTPHOSDigitDataStruct) + sizeof(AliHLTPHOSDigitHeaderStruct));

	      // Set the pointer to the first digit in the list
	      firstDigitPtr = reinterpret_cast<AliHLTPHOSDigitDataStruct*>(outBPtr + sizeof(AliHLTPHOSDigitHeaderStruct) + digitHeaderPtr->fFirstDigitOffset);

	      lastDigitPtr = reinterpret_cast<AliHLTPHOSDigitDataStruct*>(outBPtr + sizeof(AliHLTPHOSDigitHeaderStruct) + digitHeaderPtr->fLastDigitOffset);
	      
	      // Update the amount of the output buffer we have used
	      mysize += digitHeaderPtr->fNDigits*sizeof(AliHLTPHOSDigitDataStruct) + sizeof(AliHLTPHOSDigitHeaderStruct);
	    }
	  else
	    {
	      // Check if we have space for 
	      if(availableSize < digitHeaderPtr->fNDigits*sizeof(AliHLTPHOSDigitDataStruct))
		{
		  HLTError("Buffer overflow: Trying to write data of size: %d bytes. Output buffer available: %d bytes.", totSize, size);
		  return -1;
		}

	      // If we already have copied the first buffer to the output copy only the digits
	      memcpy(outBPtr, (iter->fPtr+sizeof(AliHLTPHOSDigitHeaderStruct)), digitHeaderPtr->fNDigits*sizeof(AliHLTPHOSDigitDataStruct));

	      // Check if the first digit in this buffer has a ID less than the first digit in the previous
	      if(firstDigitPtr->fID > reinterpret_cast<AliHLTPHOSDigitDataStruct*>(iter->fPtr + sizeof(AliHLTPHOSDigitDataStruct) + digitHeaderPtr->fFirstDigitOffset)->fID)
		{
		  // If that is the case we have to take care of the ordering
		  		  
		  // The last digit in the current buffer has to link to the first digit in the previous buffer
 		  lastDigitPtr->fMemOffsetNext = outBPtr + digitHeaderPtr->fFirstDigitOffset - reinterpret_cast<Int_t>(lastDigitPtr);

		  // Setting the pointer to the new first digit
		  firstDigitPtr = reinterpret_cast<AliHLTPHOSDigitDataStruct*>(outBPtr + digitHeaderPtr->fFirstDigitOffset);
		}
	      else
		{
		  // Previous last digit need to link to the current first digit
		  lastDigitPtr->fMemOffsetNext = outBPtr + digitHeaderPtr->fFirstDigitOffset;		  
		  
		  // We need to change the last digit pointer
		  lastDigitPtr = reinterpret_cast<AliHLTPHOSDigitDataStruct*>(outBPtr + digitHeaderPtr->fLastDigitOffset);
		}
	      // Update the amount of the output buffer we have used
	      mysize += digitHeaderPtr->fNDigits*sizeof(AliHLTPHOSDigitDataStruct);
	    }

	  outBPtr += mysize;
	}
    }
  
  // The digit header in the output needs to know about the position of the new first digit (the last digit is still the same)
  outputDigitHeaderPtr->fFirstDigitOffset = reinterpret_cast<Int_t>(firstDigitPtr) - reinterpret_cast<Int_t>(outputDigitHeaderPtr) + sizeof(AliHLTPHOSDigitHeaderStruct);
  
  // The digit header in the output needs to know about the position of the new last digit (the first digit is still the same)
  outputDigitHeaderPtr->fLastDigitOffset = reinterpret_cast<Int_t>(lastDigitPtr) - reinterpret_cast<Int_t>(outputDigitHeaderPtr) + sizeof(AliHLTPHOSDigitHeaderStruct);
		  

  HLTDebug("Number of digits: %d", nDigits);

  AliHLTPHOSRecPointHeaderStruct* recPointHeaderPtr = reinterpret_cast<AliHLTPHOSRecPointHeaderStruct*>(outBPtr);

  fClusterizerPtr->SetRecPointDataPtr(reinterpret_cast<AliHLTPHOSRecPointDataStruct*>(outBPtr+sizeof(AliHLTPHOSRecPointHeaderStruct)));

  nRecPoints = fClusterizerPtr->ClusterizeEvent(firstDigitPtr, size, mysize);

  recPointHeaderPtr->fNRecPoints = nRecPoints;

  mysize += sizeof(AliHLTPHOSRecPointHeaderStruct);
  
  HLTDebug("Number of clusters: %d", nRecPoints);

  AliHLTComponentBlockData clusterBd;
  FillBlockData( clusterBd );
  clusterBd.fOffset = offset;
  clusterBd.fSize = mysize;
  clusterBd.fDataType = AliHLTPHOSDefinitions::fgkRecPointDataType;
  clusterBd.fSpecification = specification;
  outputBlocks.push_back( clusterBd );

  if(false)
    {
      AliHLTComponentBlockData digitBd;
      FillBlockData(digitBd);
    }
       
  size = mysize;
  
  return 0;
}

int 
AliHLTPHOSClusterizerComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{  
  // see header file for class documentation

  const char* path="HLT/ConfigPHOS/ClusterizerComponent";

  if (cdbEntry) path = cdbEntry;

  return ConfigureFromCDBTObjString(cdbEntry);
}

int 
AliHLTPHOSClusterizerComponent::ScanConfigurationArgument(int argc, const char **argv)
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

  const char *path = "HLT/ConfigPHOS/ClusterizerComponent";

  //  ConfigureFromCDBTObjString(path);

  for (int i = 0; i < argc; i++)
    {
      ScanConfigurationArgument(i, argv);
    }

  return 0;
}

AliHLTComponent*
AliHLTPHOSClusterizerComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSClusterizerComponent();
}
