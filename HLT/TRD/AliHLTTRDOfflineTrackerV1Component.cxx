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

/** @file   AliHLTTRDOfflineTrackerV1Component.cxx
    @author 
    @date   2009-08-31
    @brief  A TRDClusterizer processing component for the HLT. 
*/

#include "AliHLTTRDOfflineTrackerV1Component.h"
#include "AliCDBManager.h"

ClassImp(AliHLTTRDOfflineTrackerV1Component)
    
AliHLTTRDOfflineTrackerV1Component::AliHLTTRDOfflineTrackerV1Component():
  AliHLTTRDTrackerV1Component()
{
  // Default constructor

}

AliHLTTRDOfflineTrackerV1Component::~AliHLTTRDOfflineTrackerV1Component()
{
  // Destructor
  // Work is Done in DoDeInit()
}

AliHLTComponent* AliHLTTRDOfflineTrackerV1Component::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDOfflineTrackerV1Component;
};

int AliHLTTRDOfflineTrackerV1Component::DoInit( int argc, const char** argv )
{
  SetOfflineParams();
  return AliHLTTRDTrackerV1Component::DoInit(argc, argv);
}

const char* AliHLTTRDOfflineTrackerV1Component::GetComponentID()
{
  // Return the component ID const char *
  return "TRDOfflineTrackerV1"; // The ID of this component
}

void AliHLTTRDOfflineTrackerV1Component::SetOfflineParams(){
  HLTFatal("You have entered the OFFLINE configuration!");
  HLTFatal("This program shall NOT run on the HLT cluster like this!");
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
    HLTFatal("You are resetting the Default Storage of the CDBManager!");
    HLTFatal("Let's hope that this program is NOT running on the HLT cluster!");
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  }else{
    HLTError("DefaultStorage was already set!");
  }
  if(AliCDBManager::Instance()->GetRun()<0){
    HLTFatal("You are resetting the CDB run number to 0!");
    HLTFatal("Let's hope that this program is NOT running on the HLT cluster!");
    AliCDBManager::Instance()->SetRun(0);
  }else{
    HLTError("Run Number was already set!");
  }
}

int AliHLTTRDOfflineTrackerV1Component::DoDeinit()
{
  return AliHLTTRDTrackerV1Component::DoDeinit();
}

int AliHLTTRDOfflineTrackerV1Component::DoEvent(const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
						  AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
						  AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
{
  return AliHLTTRDTrackerV1Component::DoEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks );
}

