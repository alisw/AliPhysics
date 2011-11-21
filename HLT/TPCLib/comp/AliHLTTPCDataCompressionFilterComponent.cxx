// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/// @file   AliHLTTPCDataCompressionFilterComponent.cxx
/// @author Matthias Richter
/// @date   2011-08-08
/// @brief  TPC component for data compression
///

#include "AliHLTTPCDataCompressionFilterComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTPluginBase.h"
#include "AliHLTSystem.h"
#include "AliHLTOUT.h"
#include "AliLog.h"

ClassImp(AliHLTTPCDataCompressionFilterComponent)

AliHLTTPCDataCompressionFilterComponent::AliHLTTPCDataCompressionFilterComponent()
  : AliHLTProcessor()
{
}

AliHLTTPCDataCompressionFilterComponent::~AliHLTTPCDataCompressionFilterComponent()
{
  /// destructor
}


const char* AliHLTTPCDataCompressionFilterComponent::GetComponentID()
{
  /// inherited from AliHLTComponent: id of the component
  return "TPCDataCompressorFilter";
}


void AliHLTTPCDataCompressionFilterComponent::GetInputDataTypes( AliHLTComponentDataTypeList& tgtList)
{
  /// inherited from AliHLTComponent: list of data types in the vector reference
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
  tgtList.push_back(AliHLTTPCDefinitions::RemainingClusterIdsDataType());
}

AliHLTComponentDataType AliHLTTPCDataCompressionFilterComponent::GetOutputDataType()
{
  /// inherited from AliHLTComponent: output data type of the component.
  return kAliHLTMultipleDataType;
}

int AliHLTTPCDataCompressionFilterComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  /// inherited from AliHLTComponent: multiple output data types of the component.
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
  tgtList.push_back(AliHLTTPCDefinitions::RemainingClusterIdsDataType());
  return tgtList.size();
}

void AliHLTTPCDataCompressionFilterComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  /// inherited from AliHLTComponent: output data size estimator
  constBase=0;
  inputMultiplier=1.;  // there should not be more data than input
}

AliHLTComponent* AliHLTTPCDataCompressionFilterComponent::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTTPCDataCompressionFilterComponent;
}

int AliHLTTPCDataCompressionFilterComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, 
						AliHLTComponentTriggerData& /*trigData*/)
{
  /// inherited from AliHLTProcessor: data processing
  if (!IsDataEvent()) return 0;

  if (GetFirstInputBlock(AliHLTTPCDefinitions::ClusterTracksCompressedDataType())!=NULL) {
    // This component is only used in conjunction with an emulation chain for compressed
    // partition cluster blocks. Blocks which are missing in HLTOUT but are existing in
    // the data stream from the parent, are forwarded. This scheme allows to automatically
    // create missing partitions in the compressed data from raw data (e.g. if an input
    // link of the HLT is broken and raw data recorded), and add it to the HLTOUT to
    // have a consistent data set. This requires individual cluster data blocks, track
    // model compression can not be used in the emulation because the clusters can not be
    // related to a particular partition.
    AliFatalClass("compressed track cluster data blocks can not be mixed. aborting");
    return -EBADF;
  }

  std::map<AliHLTUInt32_t, bool> hltoutmap;
  bool bHaveMap=false;

  for (const AliHLTComponentBlockData* pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
       pDesc!=NULL; pDesc=GetNextInputBlock()) {
    if (!bHaveMap) {
      InitMapFromHLTOUT(hltoutmap);
      bHaveMap=true;
    }
    if (hltoutmap.find(pDesc->fSpecification)!=hltoutmap.end()) {
      // block existing in HLTOUT
      continue;
    }
    Forward(pDesc);
    HLTInfo("inserting block 0x%08x", pDesc->fSpecification);
  }

  return 0;
}

int AliHLTTPCDataCompressionFilterComponent::InitMapFromHLTOUT(std::map<AliHLTUInt32_t, bool>& hltoutmap)
{
  // check the HLTOUT for availability of compressed data blocks
  AliHLTSystem* pSystem=AliHLTPluginBase::GetInstance();
  if (!pSystem) {
    // global system not initialized
    return -ENODEV;
  }
  AliHLTOUT* pHLTOUT=pSystem->RequestHLTOUT();
  if (!pHLTOUT) {
    // not HLTOUT, hence not clusters
    return 0;
  }

  for (bool bNextBlock=(pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::RemainingClustersCompressedDataType())>=0);
       bNextBlock; bNextBlock=(pHLTOUT->SelectNextDataBlock()>=0)) {
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    if (pHLTOUT->GetDataBlockDescription(dt, spec)<0)
      continue;

    hltoutmap[spec]=true;
  }

  for (bool bNextBlock=(pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::ClusterTracksCompressedDataType())>=0);
       bNextBlock; bNextBlock=(pHLTOUT->SelectNextDataBlock()>=0)) {
    // the first version of this component will not implement support for track model compression data blocks
    // to implement it
    // - decode data
    // - sort into index grid
    // - check if there is at least one cluster in a partition, that is a sufficient condition
    //   to decide whether a partition was included or not
    // The best would be to implement a class which supports the AliHLTTPCDataCompressionDecoder
    // interface and stores unpacked data in AliHLTTPCRawCluster format and fills the index
    // grid at the same time
    AliFatalClass("this functionality needs to be implemented");
  }
    
  return 0;
}

int AliHLTTPCDataCompressionFilterComponent::DoInit( int argc, const char** argv )
{
  /// inherited from AliHLTComponent: component initialisation and argument scan.
  int iResult=0;

  // component configuration
  //Stage 1: default initialization.
  //No default values until now.

  //Stage 2: OCDB. - disabled
  //TString cdbPath("HLT/ConfigTPC/");
  //cdbPath += GetComponentID();
  //
  //iResult = ConfigureFromCDBTObjString(cdbPath);
  //if (iResult < 0) 
  //  return iResult;

  //Stage 3: command line arguments.
  if (argc && (iResult = ConfigureFromArgumentString(argc, argv)) < 0)
    return iResult;

  return iResult;
}

int AliHLTTPCDataCompressionFilterComponent::DoDeinit()
{
  /// inherited from AliHLTComponent: component cleanup
  int iResult=0;

  return iResult;
}

int AliHLTTPCDataCompressionFilterComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  /// inherited from AliHLTComponent: argument scan
  int iResult=0;
  if (argc<1) return 0;
  int bMissingParam=0;
  int i=0;
  TString argument=argv[i];

  do {
    // -mode
    if (argument.CompareTo("-mode")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter=argv[i];
      if (parameter.IsDigit()) {

	return 2;
      } else {
	HLTError("invalid parameter for argument %s, expecting number instead of %s", argument.Data(), parameter.Data());
	return -EPROTO;
      }
    }

  } while (0); // using do-while only to have break available

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EPROTO;
  }

  return iResult;
}
