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

#include "AliHLTPHOSESDMakerComponent.h"
#include "AliHLTPHOSESDMaker.h"
#include "AliHLTPHOSCaloClusterContainerStruct.h"
#include "AliESDEvent.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTDataTypes.h"
#include "AliESDCaloCluster.h"

/** @file   AliHLTPHOSESDMakerComponent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  An ESD maker component for PHOS HLT
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

const AliHLTComponentDataType AliHLTPHOSESDMakerComponent::fgkInputDataTypes[]=
  {
    kAliHLTVoidDataType,{0,"",""}
  };

AliHLTPHOSESDMakerComponent gAliHLTPHOSESDMakerComponent;


AliHLTPHOSESDMakerComponent::AliHLTPHOSESDMakerComponent(): 
  AliHLTPHOSProcessor(), 
  fESDMakerPtr(0),
  fESDEventPtr(0)
{
  //See headerfile for documentation
}

AliHLTPHOSESDMakerComponent::~AliHLTPHOSESDMakerComponent()
{
  //See headerfile for documentation

  if (fESDMakerPtr)
    {
      delete fESDMakerPtr;
      fESDMakerPtr = 0;
    }
  if (fESDEventPtr)
    {
      delete fESDEventPtr;
      fESDEventPtr = 0;
    }

}


int
AliHLTPHOSESDMakerComponent::Deinit()
{
  //See headerfile for documentation

  if (fESDMakerPtr)
    {
      delete fESDMakerPtr;
      fESDMakerPtr = 0;
    }
  if (fESDEventPtr)
    {
      delete fESDEventPtr;
      fESDEventPtr = 0;
    }

  return 0;
}

const Char_t*
AliHLTPHOSESDMakerComponent::GetComponentID()
{
  //See headerfile for documentation

  return "PhosEsdMaker";
}

void
AliHLTPHOSESDMakerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
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
AliHLTPHOSESDMakerComponent::GetOutputDataType()
{
  //See headerfile for documentation

  return kAliHLTDataTypeESDObject;
}

void
AliHLTPHOSESDMakerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation

  constBase = 30;
  inputMultiplier = 4;
}
 
Int_t 
AliHLTPHOSESDMakerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) 
{
  // see header file for class documentation
 
  
  AliHLTPHOSCaloClusterContainerStruct* caloClusterContainerPtr = 0;
  
  const AliHLTComponentBlockData* iter = 0;

  UInt_t specification = 0;

  fESDMakerPtr->ResetESD();

  fESDEventPtr = new AliESDEvent();
  fESDEventPtr->CreateStdContent();
  fESDMakerPtr->SetESDEvent(fESDEventPtr);
  
  for ( iter = GetFirstInputBlock(AliHLTPHOSDefinitions::fgkCaloClusterDataType); iter != 0; iter = GetNextInputBlock()) 
    {
      specification = specification|iter->fSpecification;
      caloClusterContainerPtr = reinterpret_cast<AliHLTPHOSCaloClusterContainerStruct*>(iter->fPtr);
      fESDMakerPtr->FillESDEvent(caloClusterContainerPtr);
      //      fESDMakerPtr->SetCaloClusterContainer(caloClusterContainerPtr);
    }
  
  PushBack(fESDEventPtr, kAliHLTDataTypeESDObject|kAliHLTDataOriginPHOS, specification);
 
  if(fESDEventPtr)
    {
      delete fESDEventPtr;
      fESDEventPtr = 0;
    }

  fPhosEventCount++; 
  if(fPrintInfoModule == kTRUE)
    {
      if(fPhosEventCount%fPrintInfoFrequncyModule == 0)
      	{
	  Logging(kHLTLogInfo, __FILE__ , "writing data" , "Made ESD from event %lu",  fPhosEventCount);
	}  
    }
  return 0;
}


int
AliHLTPHOSESDMakerComponent::DoInit(int argc, const char** argv )
{
  //See headerfile for documentation

  fESDMakerPtr = new AliHLTPHOSESDMaker();
  //fESDMakerPtr->SetESDEvent(fESDEventPtr);
  //
  ScanArgumentsModule(argc, argv);
  for (int i = 0; i < argc; i++)
    {
      
    }

  return 0;
}

AliHLTComponent*
AliHLTPHOSESDMakerComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSESDMakerComponent();
}
