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
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDClusterizer.h"
#include "AliHLTTRDUtils.h"
#include "AliCDBManager.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TObjString.h"

ClassImp(AliHLTTRDOfflineClusterizerComponent)
   
AliHLTTRDOfflineClusterizerComponent::AliHLTTRDOfflineClusterizerComponent()
  :AliHLTTRDClusterizerComponent()
  ,fOffClusterizer(NULL)
{
  // Default constructor
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

void AliHLTTRDOfflineClusterizerComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data
  list.clear(); 
  AliHLTTRDClusterizerComponent::GetInputDataTypes(list);
  list.push_back(AliHLTTRDDefinitions::fgkDigitsDataType);
}

AliHLTComponentDataType AliHLTTRDOfflineClusterizerComponent::GetOutputDataType()
{
  // Get the output data type
  return kAliHLTMultipleDataType;
}

int AliHLTTRDOfflineClusterizerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // Get the output data types
  tgtList.clear();
  AliHLTTRDClusterizerComponent::GetOutputDataTypes(tgtList);
  tgtList.push_back(AliHLTTRDDefinitions::fgkHiLvlClusterDataType);
  return tgtList.size();
}

void AliHLTTRDOfflineClusterizerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  AliHLTTRDClusterizerComponent::GetOutputDataSize(constBase, inputMultiplier);
  constBase += 500;
  inputMultiplier *= 10;
}

int AliHLTTRDOfflineClusterizerComponent::DoInit( int argc, const char** argv )
{
  int iResult = 0;
  iResult=AliHLTTRDClusterizerComponent::DoInit(argc, argv);

  return iResult;
}

int AliHLTTRDOfflineClusterizerComponent::DoEvent(const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
						  AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
						  AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
{
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  if(!GetFirstInputBlock(AliHLTTRDDefinitions::fgkDigitsDataType))
    return AliHLTTRDClusterizerComponent::DoEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks );

  AliHLTUInt32_t offset = 0;

  for(const TObject *iter = GetFirstInputObject(AliHLTTRDDefinitions::fgkDigitsDataType); iter; iter = GetNextInputObject()) 
    {
      AliTRDclusterizer* clusterizer = new AliTRDclusterizer("TRDCclusterizer", "TRDCclusterizer");
      clusterizer->SetReconstructor(fReconstructor);
      clusterizer->SetUseLabels(kTRUE);

      TTree* digitsTree = dynamic_cast<TTree*>(const_cast<TObject*>(iter));
      clusterizer->ReadDigits(digitsTree);
      clusterizer->MakeClusters();
      TClonesArray* clusterArray = clusterizer->RecPoints();
      clusterizer->SetClustersOwner(kFALSE);
      
      if(fHighLevelOutput){
	if(fEmulateHLTClusters){
	  TClonesArray* temp = clusterArray;
	  clusterArray = new TClonesArray(*temp);
	  temp->Delete();
	  delete temp;
	  AliHLTTRDUtils::EmulateHLTClusters(clusterArray);
	}
	TObjString strg;
	strg.String() += clusterizer->GetNTimeBins();
	PushBack(clusterArray, AliHLTTRDDefinitions::fgkHiLvlClusterDataType, 0);
	PushBack(&strg, AliHLTTRDDefinitions::fgkHiLvlClusterDataType, 0);
      } else {
	Int_t nTimeBins = clusterizer->GetNTimeBins();
	Int_t addedSize = AliHLTTRDUtils::AddClustersToOutput(clusterArray, outputPtr+offset, nTimeBins);
	
	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset = offset;
	bd.fSize = addedSize;
	bd.fSpecification = 1;
	bd.fDataType = AliHLTTRDDefinitions::fgkClusterDataType;
	outputBlocks.push_back( bd );
	HLTDebug( "BD ptr 0x%x, offset %i, size %i, dataType %s, spec 0x%x ", bd.fPtr, bd.fOffset, bd.fSize, DataType2Text(bd.fDataType).c_str(), bd.fSpecification);
	offset += addedSize;
      }
      clusterArray->Delete();
      delete clusterArray;
      delete clusterizer;
    }

  return 0;

}
