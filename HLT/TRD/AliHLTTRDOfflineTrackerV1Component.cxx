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
#include "AliTRDrecoParam.h"
#include "AliHLTTRDDefinitions.h"

ClassImp(AliHLTTRDOfflineTrackerV1Component)
    
AliHLTTRDOfflineTrackerV1Component::AliHLTTRDOfflineTrackerV1Component()
  :AliHLTTRDTrackerV1Component()
{
  // Default constructor
  fOffline=kTRUE;
}

AliHLTTRDOfflineTrackerV1Component::~AliHLTTRDOfflineTrackerV1Component()
{
  // Destructor
  // Work is Done in DoDeInit()
}

int AliHLTTRDOfflineTrackerV1Component::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // Get the output data types
  tgtList.clear();
  AliHLTTRDTrackerV1Component::GetOutputDataTypes(tgtList);
  tgtList.push_back(AliHLTTRDDefinitions::fgkTRDOffTracksDataType);
  return tgtList.size();
}

void AliHLTTRDOfflineTrackerV1Component::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = 1000000;
  inputMultiplier = 4*((double)fOutputPercentage);
}

AliHLTComponent* AliHLTTRDOfflineTrackerV1Component::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDOfflineTrackerV1Component;
};

int AliHLTTRDOfflineTrackerV1Component::DoInit( int argc, const char** argv )
{
  int iResult = 0;
  SetOfflineParams();
  iResult=AliHLTTRDTrackerV1Component::DoInit(argc, argv);
  fRecoParam->SetStreamLevel(AliTRDrecoParam::kTracker, 1); //in order to have the friends written
  return iResult;
}

const char* AliHLTTRDOfflineTrackerV1Component::GetComponentID()
{
  // Return the component ID const char *
  return "TRDOfflineTrackerV1"; // The ID of this component
}

void AliHLTTRDOfflineTrackerV1Component::SetOfflineParams(){
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

int AliHLTTRDOfflineTrackerV1Component::DoDeinit()
{
  return AliHLTTRDTrackerV1Component::DoDeinit();
}

int AliHLTTRDOfflineTrackerV1Component::DoEvent(const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
						  AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
						  AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
{
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;
  return AliHLTTRDTrackerV1Component::DoEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks );
}
