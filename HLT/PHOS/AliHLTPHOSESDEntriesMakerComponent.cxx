
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

#include "AliHLTPHOSESDEntriesMakerComponent.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliESDEvent.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTDataTypes.h"
#include "AliESDCaloCluster.h"
#include "AliHLTESDCaloClusterMaker.h"

/** @file   AliHLTPHOSESDEntriesMakerComponent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  An ESD entries maker component for PHOS HLT
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

AliHLTPHOSESDEntriesMakerComponent gAliHLTPHOSESDEntriesMakerComponent;

AliHLTPHOSESDEntriesMakerComponent::AliHLTPHOSESDEntriesMakerComponent(): 
  AliHLTPHOSProcessor(), 
  fESDCaloClusterMakerPtr(0),
  fESDCaloClustersPtr(0)
{
  //See headerfile for documentation
}

AliHLTPHOSESDEntriesMakerComponent::~AliHLTPHOSESDEntriesMakerComponent()
{
  //See headerfile for documentation
  if (fESDCaloClusterMakerPtr)
    {
      delete fESDCaloClusterMakerPtr;
      fESDCaloClusterMakerPtr = 0;
    }
}

int
AliHLTPHOSESDEntriesMakerComponent::Deinit()
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
AliHLTPHOSESDEntriesMakerComponent::GetComponentID()
{
  //See headerfile for documentation
  return "PhosEsdEntriesMaker";
}

void
AliHLTPHOSESDEntriesMakerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //See headerfile for documentation
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkCaloClusterDataType);
}

AliHLTComponentDataType
AliHLTPHOSESDEntriesMakerComponent::GetOutputDataType()
{
  //See headerfile for documentation
  return AliHLTPHOSDefinitions::fgkESDCaloClusterDataType;
}

void
AliHLTPHOSESDEntriesMakerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation
  constBase = 30;
  inputMultiplier = 2;
}
 
Int_t 
AliHLTPHOSESDEntriesMakerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) 
{
  
  // see header file for class documentation
 
  fESDCaloClustersPtr->Delete();
  //fESDCaloClustersPtr = new TClonesArray(AliESDCaloCluster::Class(), 10);
  AliHLTCaloClusterHeaderStruct* caloClusterHeaderPtr = 0;

  const AliHLTComponentBlockData* iter = 0;
  //  UInt_t nClusters = 0;

  UInt_t clusterSpecification = 0;
  //  UInt_t digitSpecification  = 0;
  //UInt_t nDigits = 0;

  for ( iter = GetFirstInputBlock(AliHLTPHOSDefinitions::fgkCaloClusterDataType); iter != 0; iter = GetNextInputBlock()) 
    {
      clusterSpecification = clusterSpecification|iter->fSpecification;
      caloClusterHeaderPtr = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(iter->fPtr);
      HLTDebug("%d HLT clusters", caloClusterHeaderPtr->fNClusters);
      //      nClusters = fESDCaloClusterMakerPtr->FillESDCaloClusters(fESDCaloClustersPtr, caloClusterHeaderPtr);
    }

//   for(iter = GetFirstInputBlock(AliHLTPHOSDefinitions::fgkDigitDataType); iter != 0; iter =GetNextInputBlock())
//     {
//       digitSpecification = digitSpecification|iter->fSpecification;
//       digitsPtr = reinterpret_cast<AliHLTPHOSDigitDataStruct*>(iter->fPtr);
//       nDigits = iter->fSize/sizeof(AliHLTPHOSDigitDataStruct);
//       fESDCaloCellMaker->FillESDCaloCells(fESDCaloCellsPtr, digitsPtr, nDigits);
//     }
  
  // for(int i = 0; i < nClusters; i++)
  //   {
  //     HLTDebug("Cluster energy: %f", ((AliESDCaloCluster*)(fESDCaloClustersPtr->At(i)))->E());
  //   }

  ///  HLTError("Number of ESD clusters: %d", nClusters);
  PushBack(fESDCaloClustersPtr, AliHLTPHOSDefinitions::fgkESDCaloClusterDataType|kAliHLTDataOriginPHOS);
  //  PushBack(fESDCaloCellsPtr, AliHLTPHOSDefinitions::fgkESDCaloCellsDataType|kAliHLTDataOriginPHOS, digitSpecification);
  
  //delete fESDCaloClustersPtr; 

  return 0;
}


int
AliHLTPHOSESDEntriesMakerComponent::DoInit(int argc, const char** argv )
{
  //See headerfile for documentation

  fESDCaloClusterMakerPtr = new AliHLTESDCaloClusterMaker();

  fESDCaloClustersPtr = new TClonesArray(AliESDCaloCluster::Class(), 10);

  ScanArgumentsModule(argc, argv);
  for (int i = 0; i < argc; i++)
    {
      
    }

  return 0;
}

AliHLTComponent*
AliHLTPHOSESDEntriesMakerComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSESDEntriesMakerComponent();
}
