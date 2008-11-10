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

#include "AliHLTPHOSDigitMakerComponent.h"
#include "AliHLTPHOSDigitMaker.h"
#include "TTree.h"
#include "AliHLTPHOSProcessor.h"
#include "AliHLTPHOSRcuCellEnergyDataStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSChannelDataHeaderStruct.h"
#include "AliHLTPHOSChannelDataStruct.h"
#include "TClonesArray.h"
#include "TFile.h"
#include <sys/stat.h>
#include <sys/types.h>


/** 
 * @file   AliHLTPHOSDigitMakerComponent.cxx
 * @author Oystein Djuvsland
 * @date   
 * @brief  A digit maker component for PHOS HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


const AliHLTComponentDataType AliHLTPHOSDigitMakerComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};

AliHLTPHOSDigitMakerComponent gAliHLTPHOSDigitMakerComponent;

AliHLTPHOSDigitMakerComponent::AliHLTPHOSDigitMakerComponent() :
  AliHLTPHOSProcessor(),
  fDigitMakerPtr(0),
  fDigitContainerPtr(0)
  //  fEvtCnt(0)
{
  //see header file for documentation
}


AliHLTPHOSDigitMakerComponent::~AliHLTPHOSDigitMakerComponent()
{
  //see header file for documentation
}

int 
AliHLTPHOSDigitMakerComponent::Deinit()
{ 
  //see header file for documentation
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
  //see header file for documentation
  return "PhosDigitMaker";
}


void
AliHLTPHOSDigitMakerComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
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
AliHLTPHOSDigitMakerComponent::GetOutputDataType()
{
  //see header file for documentation
  return AliHLTPHOSDefinitions::fgkDigitDataType;
}


void 
AliHLTPHOSDigitMakerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //see header file for documentation
  constBase = 0;
  inputMultiplier = (float)sizeof(AliHLTPHOSDigitDataStruct)/sizeof(AliHLTPHOSChannelDataStruct) + 1;
}

int 
AliHLTPHOSDigitMakerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
					AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
					std::vector<AliHLTComponentBlockData>& outputBlocks)
{
  //see header file for documentation
  UInt_t tSize            = 0;
  UInt_t offset           = 0; 
  UInt_t mysize           = 0;
  Int_t digitCount        = 0;
  Int_t ret               = 0;

  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = 0; 
  unsigned long ndx; 

  UInt_t specification = 0;
  AliHLTPHOSChannelDataHeaderStruct* tmpChannelData = 0;
  
  fDigitMakerPtr->SetDigitDataPtr(reinterpret_cast<AliHLTPHOSDigitDataStruct*>(outputPtr));

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      
      if(iter->fDataType != AliHLTPHOSDefinitions::fgkChannelDataType)
	{
	  HLTDebug("Data block is not of type fgkChannelDataType");
	  continue;
	}

      specification |= iter->fSpecification;
      tmpChannelData = reinterpret_cast<AliHLTPHOSChannelDataHeaderStruct*>(iter->fPtr);
    
      ret = fDigitMakerPtr->MakeDigits(tmpChannelData, size-(digitCount*sizeof(AliHLTPHOSDigitDataStruct)));
      if(ret == -1) 
	{
	  HLTError("Trying to write over buffer size");
	  return -ENOBUFS;
	}
      digitCount += ret; 
    }
  
  mysize = digitCount*sizeof(AliHLTPHOSDigitDataStruct);

  HLTDebug("# of digits: %d, used memory size: %d, available size: %d", digitCount, mysize, size);
  
  AliHLTComponentBlockData bd;
  FillBlockData( bd );
  bd.fOffset = offset;
  bd.fSize = mysize;
  bd.fDataType = AliHLTPHOSDefinitions::fgkDigitDataType;
  bd.fSpecification = specification;
  outputBlocks.push_back(bd);
       
  if( tSize > size )
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSDigitMakerComponent::DoEvent", "Too much data", "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.", tSize, size );
      return EMSGSIZE;
    }

  fDigitMakerPtr->Reset();

  size = mysize; 

  return 0;
}


int
AliHLTPHOSDigitMakerComponent::DoInit(int argc, const char** argv )
{
  //see header file for documentation

  fDigitMakerPtr = new AliHLTPHOSDigitMaker();
  
  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-lowgainfactor", argv[i]))
	{
	  fDigitMakerPtr->SetGlobalLowGainFactor(atof(argv[i+1]));
	}
      if(!strcmp("-highgainfactor", argv[i]))
	{
	  fDigitMakerPtr->SetGlobalHighGainFactor(atof(argv[i+1]));
	}
      if(!strcmp("-reverseorder", argv[i]))
	{
	  fDigitMakerPtr->SetOrdered(false);
	}
    }
 
  //fDigitMakerPtr->SetDigitThreshold(2);

  return 0;
}

AliHLTComponent*
AliHLTPHOSDigitMakerComponent::Spawn()
{
  //see header file for documentation
  return new AliHLTPHOSDigitMakerComponent();
}
