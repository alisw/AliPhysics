
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

#include "AliHLTPHOSESDCaloClusterMakerComponent.h"
#include "AliHLTPHOSESDCaloClusterMaker.h"
#include "AliHLTPHOSCaloClusterContainerStruct.h"
#include "AliESDEvent.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTDataTypes.h"
#include "AliESDCaloCluster.h"

/** @file   AliHLTPHOSESDCaloClusterMakerComponent.cxx
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


AliHLTPHOSESDCaloClusterMakerComponent::AliHLTPHOSESDCaloClusterMakerComponent(): 
  AliHLTPHOSProcessor(), 
  fESDCaloClusterMakerPtr(0),
  fESDCaloClustersPtr(0)
{
  //See headerfile for documentation
}

AliHLTPHOSESDCaloClusterMakerComponent::~AliHLTPHOSESDCaloClusterMakerComponent()
{
  //See headerfile for documentation

  if (fESDCaloClusterMakerPtr)
    {
      delete fESDCaloClusterMakerPtr;
      fESDCaloClusterMakerPtr = 0;
    }
}


int
AliHLTPHOSESDCaloClusterMakerComponent::Deinit()
{
  //See headerfile for documentation

  if (fESDCaloClusterMakerPtr)
    {
      delete fESDCaloClusterMakerPtr;
      fESDCaloClusterMakerPtr = 0;
    }
  return 0;
}

const Char_t*
AliHLTPHOSESDCaloClusterMakerComponent::GetComponentID()
{
  //See headerfile for documentation

  return "PhosEsdCaloClusterMaker";
}

void
AliHLTPHOSESDCaloClusterMakerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //See headerfile for documentation
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkCaloClusterDataType);
//   const AliHLTComponentDataType* pType=fgkInputDataTypes;
//   while (pType->fID!=0)
//     {
//       list.push_back(*pType);
//       pType++;
//     }
}

AliHLTComponentDataType
AliHLTPHOSESDCaloClusterMakerComponent::GetOutputDataType()
{
  //See headerfile for documentation

  return AliHLTPHOSDefinitions::fgkCaloClusterDataType;
}

void
AliHLTPHOSESDCaloClusterMakerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation

  constBase = 30;
  inputMultiplier = 2;
}
 
Int_t 
AliHLTPHOSESDCaloClusterMakerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) 
{
  
  // see header file for class documentation
 
//   fESDCaloClustersPtr->Delete();
//   AliHLTPHOSCaloClusterContainerStruct* caloClusterContainerPtr = 0;
  
//   const AliHLTComponentBlockData* iter = 0;

//   UInt_t specification = 0;

//   for ( iter = GetFirstInputBlock(AliHLTPHOSDefinitions::fgkCaloClusterDataType); iter != 0; iter = GetNextInputBlock()) 
//     {
//       specification = specification|iter->fSpecification;
//       caloClusterContainerPtr = reinterpret_cast<AliHLTPHOSCaloClusterContainerStruct*>(iter->fPtr);
//       fESDCaloClusterMakerPtr->FillESDCaloClusters(fESDCaloClustersPtr, caloClusterContainerPtr);
//       //      fESDCaloClusterMakerPtr->SetCaloClusterContainer(caloClusterContainerPtr);
//     }
  
//   PushBack(fESDCaloClustersPtr, AliHLTPHOSDefinitions::fgkESDCaloClusterDataType|kAliHLTDataOriginPHOS, specification);

  return 0;
}


int
AliHLTPHOSESDCaloClusterMakerComponent::DoInit(int argc, const char** argv )
{
  //See headerfile for documentation

//   fESDCaloClusterMakerPtr = new AliHLTPHOSESDCaloClusterMaker();

//   fESDCaloClustersPtr = new TClonesArray("AliESDCaloCluster", 0);

//   ScanArgumentsModule(argc, argv);
//   for (int i = 0; i < argc; i++)
//     {
      
//     }

  return 0;
}

AliHLTComponent*
AliHLTPHOSESDCaloClusterMakerComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSESDCaloClusterMakerComponent();
}
