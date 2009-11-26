// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors:                                                       *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTRDOfflineClusterizerComponent.cxx
    @author 
    @date   
    @brief  
*/

// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

#include "AliHLTTRDOfflineClusterizerComponent.h"
#include "AliCDBManager.h"

ClassImp(AliHLTTRDOfflineClusterizerComponent)
   
AliHLTTRDOfflineClusterizerComponent::AliHLTTRDOfflineClusterizerComponent():
  AliHLTTRDClusterizerComponent()
{
  // Default constructor
  fOffline = kTRUE;
}

AliHLTTRDOfflineClusterizerComponent::~AliHLTTRDOfflineClusterizerComponent()
{
  // Destructor
  // Work is Done in DoDeInit()
}

AliHLTComponent* AliHLTTRDOfflineClusterizerComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDOfflineClusterizerComponent;
};

const char* AliHLTTRDOfflineClusterizerComponent::GetComponentID()
{
  // Return the component ID const char *
  return "TRDOfflineClusterizer"; // The ID of this component
}

int AliHLTTRDOfflineClusterizerComponent::DoInit( int argc, const char** argv )
{
  SetOfflineParams();
  return AliHLTTRDClusterizerComponent::DoInit(argc, argv);
}

void AliHLTTRDOfflineClusterizerComponent::SetOfflineParams(){
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
    HLTFatal("You are resetting the Default Storage of the CDBManager!");
    HLTFatal("Let's hope that this program is NOT running on the HLT cluster!");
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  }
  if(AliCDBManager::Instance()->GetRun()<0){
    HLTFatal("You are resetting the CDB run number to 0!");
    HLTFatal("Let's hope that this program is NOT running on the HLT cluster!");
    AliCDBManager::Instance()->SetRun(0);
  }
}

int AliHLTTRDOfflineClusterizerComponent::DoDeinit()
{
  return AliHLTTRDClusterizerComponent::DoDeinit();
}

int AliHLTTRDOfflineClusterizerComponent::DoEvent(const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
						  AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
						  AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
{
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;
  return AliHLTTRDClusterizerComponent::DoEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks );
}

