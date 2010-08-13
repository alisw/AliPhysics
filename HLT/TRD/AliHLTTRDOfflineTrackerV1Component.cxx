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
    @author Theodor Rascanu
    @date   
    @brief  Handles high level (object based) and low level (pod) input data (clusters). For debug purposes only
*/

#include "AliHLTTRDOfflineTrackerV1Component.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDUtils.h"
#include "AliTRDrecoParam.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDReconstructor.h"
#include "AliCDBManager.h"
#include "AliESDEvent.h"
#include "TClonesArray.h"
#include "TObjString.h"

#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"

ClassImp(AliHLTTRDOfflineTrackerV1Component)
    
AliHLTTRDOfflineTrackerV1Component::AliHLTTRDOfflineTrackerV1Component()
  :AliHLTTRDTrackerV1Component()
{
  // Default constructor
}

AliHLTTRDOfflineTrackerV1Component::~AliHLTTRDOfflineTrackerV1Component()
{
  // Destructor
  // Work is Done in DoDeInit()
}

void AliHLTTRDOfflineTrackerV1Component::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data
  list.clear();
  AliHLTTRDTrackerV1Component::GetInputDataTypes(list);
  list.push_back(AliHLTTRDDefinitions::fgkHiLvlClusterDataType);
}

AliHLTComponentDataType AliHLTTRDOfflineTrackerV1Component::GetOutputDataType()
{
  // Get the output data type
  return kAliHLTMultipleDataType;
}

int AliHLTTRDOfflineTrackerV1Component::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // Get the output data types
  tgtList.clear();
  tgtList.push_back(AliHLTTRDDefinitions::fgkHiLvlTracksDataType);
  return tgtList.size();
}

void AliHLTTRDOfflineTrackerV1Component::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  AliHLTTRDTrackerV1Component::GetOutputDataSize(constBase, inputMultiplier);
  constBase += 500;
}

AliHLTComponent* AliHLTTRDOfflineTrackerV1Component::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDOfflineTrackerV1Component;
};

const char* AliHLTTRDOfflineTrackerV1Component::GetComponentID()
{
  // Return the component ID const char *
  return "TRDOfflineTrackerV1"; // The ID of this component
}

int AliHLTTRDOfflineTrackerV1Component::SetParams()
{
  int iResult = AliHLTTRDTrackerV1Component::SetParams();
  fRecoParam->SetStreamLevel(AliTRDrecoParam::kTracker, 1); // in order to have the friends written
  return iResult;
}

int AliHLTTRDOfflineTrackerV1Component::DoEvent(const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
						  AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
						  AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
{
  if(!IsDataEvent())return 0;
  
  if(!GetFirstInputBlock(AliHLTTRDDefinitions::fgkHiLvlClusterDataType))
    return AliHLTTRDTrackerV1Component::DoEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks );

  for(const TObject *iter = GetFirstInputObject(AliHLTTRDDefinitions::fgkHiLvlClusterDataType); iter; iter = GetNextInputObject()) 
    {
      TClonesArray* clusterArray = dynamic_cast<TClonesArray*>(const_cast<TObject*>(iter));
      if(!clusterArray)continue;
      TObjString* strg = dynamic_cast<TObjString*>(const_cast<TObject*>(GetNextInputObject()));
      if(!strg)continue;
      
      fNtimeBins = strg->String().Atoi();
      fESD->Reset();
      AliTRDtrackerV1::SetNTimeBins(fNtimeBins);
      HLTDebug("TClonesArray of clusters: nbEntries = %i", clusterArray->GetEntriesFast());
      fTracker->LoadClusters(clusterArray);
      fTracker->Clusters2Tracks(fESD);
      Int_t nTracks = fESD->GetNumberOfTracks();
      HLTInfo("Number of tracks  == %d ==", nTracks);  
      TClonesArray* trdTracks = fTracker->GetListOfTracks();

      if(fEmulateHLTTracks && trdTracks){
	trdTracks = new TClonesArray(*trdTracks);
	AliHLTTRDUtils::EmulateHLTTracks(trdTracks);
      }

      AliHLTUInt32_t spec = GetSpecification(iter);
      if(trdTracks)
	PushBack(trdTracks, AliHLTTRDDefinitions::fgkHiLvlTracksDataType, spec);
      else{
	TClonesArray temp("AliTRDtrackV1");
	PushBack(&temp, AliHLTTRDDefinitions::fgkHiLvlTracksDataType, spec);
      }
      PushBack(strg, AliHLTTRDDefinitions::fgkHiLvlTracksDataType, spec);
      fTracker->UnloadClusters();
      AliTRDReconstructor::SetClusters(0x0);

      if(fEmulateHLTTracks && trdTracks){
	trdTracks->Delete();
	delete trdTracks;
      }
    }
  return 0;
}
