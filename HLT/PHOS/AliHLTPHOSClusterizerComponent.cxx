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
#include "AliHLTCaloRecPointDataStruct.h"
#include "AliHLTCaloRecPointHeaderStruct.h"
#include "AliHLTPHOSGeometry.h"
#include "AliHLTCaloClusterAnalyser.h"




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
#include "AliHLTCaloDefinitions.h"
#include "AliHLTPHOSGeometry.h"
#include "AliHLTPHOSRecoParamHandler.h"

AliHLTPHOSClusterizerComponent gAliHLTPHOSClusterizerComponent;

AliHLTPHOSClusterizerComponent::AliHLTPHOSClusterizerComponent(): 
  AliHLTCaloClusterizerComponent("PHOS")
{
  //See headerfile for documentation

  fDataOrigin = const_cast<char*>(kAliHLTDataOriginPHOS);

  //AliHLTPHOSGeometry *geom = new AliHLTPHOSGeometry;
  
}

AliHLTPHOSClusterizerComponent::~AliHLTPHOSClusterizerComponent()
{
  //See headerfile for documentation
}

void
AliHLTPHOSClusterizerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //See headerfile for documentation
  list.clear();
  list.push_back(AliHLTCaloDefinitions::fgkDigitDataType|kAliHLTDataOriginPHOS);
}

AliHLTComponentDataType
AliHLTPHOSClusterizerComponent::GetOutputDataType()
{
  //See headerfile for documentation
  return kAliHLTDataTypeCaloCluster|kAliHLTDataOriginPHOS;
}

void
AliHLTPHOSClusterizerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation
  constBase = sizeof(AliHLTCaloRecPointHeaderStruct) + sizeof(AliHLTCaloRecPointDataStruct) + (sizeof(AliHLTCaloDigitDataStruct) << 7); //Reasonable estimate... ;
  inputMultiplier = 2.0;
}

const Char_t*
AliHLTPHOSClusterizerComponent::GetComponentID()
{
  //See headerfile for documentation
  return "PhosClusterizer";
}

AliHLTComponent*
AliHLTPHOSClusterizerComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSClusterizerComponent();
}

int AliHLTPHOSClusterizerComponent::DoInit(int argc, const char** argv)
{
   
   fRecoParamsPtr = new AliHLTPHOSRecoParamHandler(); 
    
    return AliHLTCaloClusterizerComponent::DoInit(argc, argv);
}

int AliHLTPHOSClusterizerComponent::DoDeinit()
{
    if(fRecoParamsPtr) 
    {
       delete fRecoParamsPtr;
       fRecoParamsPtr = 0;
    }
    return AliHLTCaloClusterizerComponent::DoDeinit();
}


Int_t AliHLTPHOSClusterizerComponent::InitialiseGeometry()
{
 
  fAnalyserPtr->SetGeometry(new AliHLTPHOSGeometry);
 
  return 0;
}
