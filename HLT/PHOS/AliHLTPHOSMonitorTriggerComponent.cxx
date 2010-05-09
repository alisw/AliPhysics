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
#include "AliHLTPHOSDefinitions.h" // PTH
#include "AliHLTPHOSMonitorTriggerComponent.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTDataTypes.h"

/** @file   AliHLTPHOSMonitorTriggerComponent.h
    @author Oystein Djuvsland
    @date   
    @brief  
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

const AliHLTComponentDataType AliHLTPHOSMonitorTriggerComponent::fgkInputDataTypes[]=
  {
    kAliHLTVoidDataType,{0,"",""}
  };

AliHLTPHOSMonitorTriggerComponent gAliHLTPHOSMonitorTriggerComponent;


AliHLTPHOSMonitorTriggerComponent::AliHLTPHOSMonitorTriggerComponent(): 
  AliHLTCaloProcessor(), 
  fCheckClusterEnergy(false),
  fCheckClusterMultiplicities(false),
  fClusterEnergyThreshold(1),
  fMultiplicityThreshold(5),
  fMultEnergyThreshold(0.5),
  fDigitMultiplicityThreshold(16),
  fMultDigitMultiplicityThreshold(9),
  fLowerCentrality(0.5),
  fUpperCentrality(0.9)
{
  //See headerfile for documentation
}

AliHLTPHOSMonitorTriggerComponent::~AliHLTPHOSMonitorTriggerComponent()
{
  //See headerfile for documentation

  
}


int
AliHLTPHOSMonitorTriggerComponent::Deinit()
{
  //See headerfile for documentation

  
  return 0;
}

const Char_t*
AliHLTPHOSMonitorTriggerComponent::GetComponentID()
{
  //See headerfile for documentation

  return "PhosMonitorTrigger";
}

void
AliHLTPHOSMonitorTriggerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkClusterDataType);
}

AliHLTComponentDataType
AliHLTPHOSMonitorTriggerComponent::GetOutputDataType()
{
  //See headerfile for documentation
  return AliHLTPHOSDefinitions::fgkSandboxDataType;
}

void
AliHLTPHOSMonitorTriggerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation
  constBase = sizeof(int);
  inputMultiplier = 0;
}

int
AliHLTPHOSMonitorTriggerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
					   AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/,
                                        std::vector<AliHLTComponentBlockData>& /*outputBlocks*/)
{
  //See headerfile for documentation

  Int_t specification = 0;
  Bool_t monitorflag = 0;
  UInt_t ndx = 0;
  
  const AliHLTComponentBlockData* iter = NULL; 

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      if ( iter->fDataType != AliHLTPHOSDefinitions::fgkDDLPackedRawDataType  )
	{
	  HLTDebug("Data block is not of type fgkDDLPackedRawDataType");
	  continue; 
	}
      specification |= iter->fSpecification;

      monitorflag += CheckClusters(reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(iter->fPtr));
	
    }

  if(monitorflag == true)
    {
      Int_t blockCount = 1;
      int ret = ReserveEventDoneData( 1+2+blockCount*4 ); 
      if(!ret)
	{
	  AliHLTUInt32_t eddWord;
	  
	  eddWord = 5; // kAliHLTEventDoneFlagMonitorEventCmd
	  PushEventDoneData( eddWord );
	  eddWord = 4; // kAliHLTEventDoneMonitorListCmd
	  PushEventDoneData( eddWord );
	  eddWord = blockCount;
	  PushEventDoneData( eddWord );
	  
	  //Int_t blockIndex = 0;
	  // low data type
	  eddWord = 0;
	 
	  for ( unsigned ii=0; ii<4; ii++ )
	    {
	      eddWord |= ((AliHLTUInt32_t)(AliHLTPHOSDefinitions::fgkDDLPackedRawDataType.fID[8-1-ii])) << (ii*8);
	    }
	  PushEventDoneData( eddWord );
	  
	  // high data type
	  eddWord = 0;
	  for ( unsigned ii=0; ii<4; ii++ )
	    {
	      eddWord |= ((AliHLTUInt32_t)(AliHLTPHOSDefinitions::fgkDDLPackedRawDataType.fID[8-5-ii])) << (ii*8);
	    }
	  PushEventDoneData( eddWord );
	  
	  // data spec
	  eddWord = 0;
	  for ( unsigned ii=0; ii<4; ii++ )
	    {
	      eddWord |= ((AliHLTUInt32_t)(AliHLTPHOSDefinitions::fgkDDLPackedRawDataType.fOrigin[4-1-ii])) << (ii*8);
	    }
	  PushEventDoneData( eddWord );
	  
	  eddWord = specification;
	  PushEventDoneData( eddWord );
	}
      else
	{
	  HLTWarning("Could not get Event Done Data blocks");
	}
    }

  return 0;
}

Bool_t
AliHLTPHOSMonitorTriggerComponent::CheckClusters(AliHLTCaloClusterHeaderStruct* clusterHeaderPtr)
{
  //See headerfile for documentation
  
  UInt_t nClusters = clusterHeaderPtr->fNClusters;
//  Float_t* ampFracPtr = 0;

  AliHLTCaloClusterDataStruct* clusterPtr = 0;

  clusterPtr = reinterpret_cast<AliHLTCaloClusterDataStruct*>(clusterHeaderPtr + sizeof(AliHLTCaloClusterHeaderStruct));
      
  for(UInt_t n = 0; n < nClusters; n++)
    {
      if(fCheckClusterEnergy == true && clusterPtr->fEnergy > fClusterEnergyThreshold && clusterPtr->fNCells > fDigitMultiplicityThreshold)
	{
/*	  ampFracPtr = &(clusterPtr->fCellsAmpFraction);
	  for(UInt_t i = 0; i < clusterPtr->fNCells; i++)
	    {
	      if(*ampFracPtr > fLowerCentrality && *ampFracPtr < fUpperCentrality) return true;
	      ampFracPtr += 6; // 6 = sizeof(Short_t) + sizeof(Float_t)
	    }*/
	}
      if(fCheckClusterMultiplicities == true && clusterPtr->fEnergy > fMultEnergyThreshold && clusterPtr->fNCells > fMultDigitMultiplicityThreshold)
	{
	  nClusters++;
	  if(nClusters > fMultiplicityThreshold)
	    {
	      return true;
	    }
	}
      clusterPtr += sizeof(AliHLTCaloClusterDataStruct) + 5*clusterPtr->fNCells; //5 = sizeof(Short_t) + sizeof(Float_t) - 1(pair already included)
    }

  return false;

}


int
AliHLTPHOSMonitorTriggerComponent::DoInit(int argc, const char** argv )
{
  //See headerfile for documentation

  for (int i = 0; i < argc; i++)
    {
      if(!strcmp("-checkenergy", argv[i]))
	{
	  fCheckClusterEnergy = true;
	}
      if(!strcmp("-checkmultiplicity", argv[i]))
	{
	  fCheckClusterMultiplicities = true;
	}
      if(!strcmp("-energythreshold", argv[i]))
	{
	  fClusterEnergyThreshold = atof(argv[i+1]);
	  fCheckClusterEnergy = true;
	}
      if(!strcmp("-multiplicityinclusterE", argv[i]))
	{
	  fDigitMultiplicityThreshold = atoi(argv[i+1]);
	}
      if(!strcmp("-multiplicitythreshold", argv[i]))
	{
	  fMultiplicityThreshold = atoi(argv[i+1]);
	  fCheckClusterMultiplicities = true;
	}
      if(!strcmp("-multiplicityinclusterMult", argv[i]))
	{
	  fMultDigitMultiplicityThreshold = atoi(argv[i+1]);
	}
      if(!strcmp("-centralitylimits", argv[i]))
	{
	  fLowerCentrality = atof(argv[i+1]);
	  fUpperCentrality = atof(argv[i+2]);
	} 
    }
  return 0;
}

AliHLTComponent*
AliHLTPHOSMonitorTriggerComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSMonitorTriggerComponent();
}
