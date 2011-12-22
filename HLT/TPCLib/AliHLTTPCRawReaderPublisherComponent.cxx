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

/// @file   AliHLTTPCRawReaderPublisherComponent.cxx
/// @author Matthias Richter
/// @date   2011-08-08
/// @brief  
///

#include "AliHLTTPCRawReaderPublisherComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTPluginBase.h"
#include "AliHLTSystem.h"
#include "AliHLTOUT.h"
#include "AliLog.h"
#include <vector>

ClassImp(AliHLTTPCRawReaderPublisherComponent)

AliHLTTPCRawReaderPublisherComponent::AliHLTTPCRawReaderPublisherComponent()
  : AliHLTRawReaderPublisherComponent()
  , fArraySelected(NULL)
{
}

AliHLTTPCRawReaderPublisherComponent::~AliHLTTPCRawReaderPublisherComponent()
{
  /// destructor
}


const char* AliHLTTPCRawReaderPublisherComponent::GetComponentID()
{
  /// inherited from AliHLTComponent: id of the component
  return "TPCRawReaderPublisher";
}

AliHLTComponent* AliHLTTPCRawReaderPublisherComponent::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTTPCRawReaderPublisherComponent;
}

int AliHLTTPCRawReaderPublisherComponent::GetEvent(const AliHLTComponentEventData& evtData, 
						   AliHLTComponentTriggerData& trigData, 
						   AliHLTUInt8_t* outputPtr, 
						   AliHLTUInt32_t& size, 
						   vector<AliHLTComponentBlockData>& outputBlocks)
{
  /// inherited from AliHLTProcessor: data processing
  if (!IsDataEvent()) return 0;


  return AliHLTRawReaderPublisherComponent::GetEvent(evtData, trigData, outputPtr, size, outputBlocks);
}

int AliHLTTPCRawReaderPublisherComponent::InitMapFromHLTOUT(std::map<AliHLTUInt32_t, bool>& hltoutmap)
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

int AliHLTTPCRawReaderPublisherComponent::DoInit( int argc, const char** argv )
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

int AliHLTTPCRawReaderPublisherComponent::DoDeinit()
{
  /// inherited from AliHLTComponent: component cleanup
  int iResult=0;

  return iResult;
}

int AliHLTTPCRawReaderPublisherComponent::ScanConfigurationArgument(int argc, const char** argv)
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

int AliHLTTPCRawReaderPublisherComponent::GetSpecificationFromEquipmentId(int id, AliHLTUInt32_t &specification) const
{
  /// inherited from AliHLTRawReaderPublisherComponent: get specification

  // FIXME: add common functionality to AliHLTDAQ
  int partition;
  int slice;
  if (id < 840) {
    partition = id % 2;
    slice = (id - 768) / 2;
  } else {
    partition = (id % 4) + 2;
    slice = (id - 840) / 4;
  }
  specification=(slice<<24)|(slice<<16)|(partition<<8)|partition;

  return 0;
}

bool AliHLTTPCRawReaderPublisherComponent::IsSelected(int /*equipmentId*/) const
{
  /// inherited from AliHLTRawReaderPublisherComponent: check if a block is selected or not

  // TODO: implement logic
  return false;
}
