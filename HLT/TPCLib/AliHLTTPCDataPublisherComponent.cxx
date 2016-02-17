// $Id$
//**************************************************************************
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTPCDataPublisherComponent.cxx
/// @author Matthias Richter
/// @date   2011-08-08
/// @brief  
///

#include "AliHLTTPCDataPublisherComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTTPCClusterMCData.h"
#include "AliHLTTPCDataCompressionDecoder.h"
#include "AliHLTPluginBase.h"
#include "AliHLTSystem.h"
#include "AliHLTOUT.h"
#include "AliHLTDAQ.h"
#include "AliHLTTemplates.h"
#include "AliLog.h"
#include <vector>
#include <memory>
#include <algorithm>

ClassImp(AliHLTTPCDataPublisherComponent)

AliHLTTPCDataPublisherComponent::AliHLTTPCDataPublisherComponent()
  : AliHLTRawReaderPublisherComponent()
  , fMode(kPublisherModeDefault)
  , fArraySelected(NULL)
  , fClusters(NULL)
  , fpDecoder(NULL)
{
  /// constructor
}

AliHLTTPCDataPublisherComponent::~AliHLTTPCDataPublisherComponent()
{
  /// destructor
  if (fpDecoder) delete fpDecoder;
  fpDecoder=NULL;
}


const char* AliHLTTPCDataPublisherComponent::GetComponentID()
{
  /// inherited from AliHLTComponent: id of the component
  return "TPCDataPublisher";
}

AliHLTComponent* AliHLTTPCDataPublisherComponent::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTTPCDataPublisherComponent;
}

int AliHLTTPCDataPublisherComponent::GetEvent(const AliHLTComponentEventData& evtData, 
						   AliHLTComponentTriggerData& trigData, 
						   AliHLTUInt8_t* outputPtr, 
						   AliHLTUInt32_t& size, 
						   AliHLTComponentBlockDataList& outputBlocks)
{
  /// inherited from AliHLTProcessor: data processing
  if (!IsDataEvent()) return 0;

  int iResult=0;

  AliHLTComponentBlockDataList clusterBlocks;
  AliHLTUInt32_t offset=0;
  AliHLTUInt32_t capacity=size;
  size=0;
  if (fClusters) {
    fClusters->Clear();
    if (CheckMode(kPublishClustersAll)) {
      // set the target buffer only if the clusters should be published
      fClusters->SetTargetBuffer(outputPtr+offset, capacity-offset);
    } else if (CheckMode(kRegisterClusterBlocks)) {
      // data blocks are registered in the container, track model cluster blocks
      // are unpacked but not stored in order to find the included partitions
      //fClusters->
    }
    if (CheckMode(kPublishClustersAll) ||
	CheckMode(kRegisterClusterBlocks)) {
      if ((iResult=ReadClusterFromHLTOUT(fClusters))>=0) {
	if ((iResult=fClusters->GetState())>=0) {
	  if (fClusters->CopyBlockDescriptors(clusterBlocks)>0) {
	    for (AliHLTComponentBlockDataList::const_iterator bd=clusterBlocks.begin();
		 bd!=clusterBlocks.end(); bd++) {
	      if (offset<bd->fOffset+bd->fSize)
		offset=bd->fOffset+bd->fSize;
	    }
	  }
	} else if (iResult==-ENOSPC) {
	  offset=fClusters->GetBlockCount()*sizeof(AliHLTTPCRawClusterData)+
	    fClusters->GetClusterCount()*sizeof(AliHLTTPCRawCluster);
	  iResult=0; // keep going to also accumulate the size for raw data blocks
	}
      }
      if (iResult==-ENODATA) {
	// return indicates absence of compressed clusters in HLTOUT
	// but is not treated as an error further downstream
	iResult=0;
      }
    }
  }

  if (offset<=capacity) {
    size=capacity-offset;
    outputPtr+=offset;
  } else {
    // there is clearly not enough space, keep the full buffer to
    // publish the raw data blocks and determine the size of those
    // data will be overwritten
    size=capacity;
  }
  if (iResult>=0) {
    unsigned firstBlock=outputBlocks.size();
    iResult=AliHLTRawReaderPublisherComponent::GetEvent(evtData, trigData, outputPtr, size, outputBlocks);
    if (iResult==-ENOSPC) {
      // not enough space in the buffer, fMaxSize has been updated by base class
      fMaxSize+=offset;
    } else if (iResult>=0) {
      if (outputBlocks.size()>firstBlock && CheckMode(kPublishRawFiltered)) {
	AliInfo(Form("publishing %lu DDL(s) for emulation of compressed TPC clusters", outputBlocks.size()-firstBlock));
      }
      // correct for the shifted buffer which was provided to the
      // GetEvent method
      for (AliHLTComponentBlockDataList::iterator bd=outputBlocks.begin();
	   bd!=outputBlocks.end(); bd++) {
	if (firstBlock>0) {firstBlock--; continue;}
	bd->fOffset+=offset;
      }
      offset+=size;
    }
  }

  if (iResult>=0 && capacity<offset && fMaxSize<(int)offset) {
    // update the size requirement
    fMaxSize=offset;
    outputBlocks.clear();
    iResult=-ENOSPC;
  }

  if (iResult>=0) {
    size=offset;
    if ( clusterBlocks.size()>0 && CheckMode(kPublishClustersAll) ) { // write out cluster blocks with proper data type
      for (AliHLTComponentBlockDataList::iterator bd=clusterBlocks.begin(); bd!=clusterBlocks.end(); bd++) {
	bd->fDataType = (AliHLTTPCDefinitions::fgkRawClustersDataType  | kAliHLTDataOriginTPC );
      }
      outputBlocks.insert(outputBlocks.begin(), clusterBlocks.begin(), clusterBlocks.end());
    }
  }

  return iResult;
}

int AliHLTTPCDataPublisherComponent::ReadClusterFromHLTOUT(AliHLTTPCDataPublisherComponent::AliRawClusterContainer* pContainer)
{
  // check the HLTOUT for availability of compressed data blocks
  int iResult=0;
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

  if (!fpDecoder) {
    fpDecoder=new AliHLTTPCDataCompressionDecoder;
  }

  if (!fpDecoder) {
    AliError("failed to create decoder instance");
    return -ENODEV;
  }

  AliHLTTPCDataCompressionDecoder& decoder=*fpDecoder;
  decoder.Clear();
  decoder.SetVerbosity(GetVerbosity());  

  bool bHavePartitionRawData=false;
  bool bHavePartitionCompressedData=false;

  bool bNextBlock=false;
  // add cluster id and mc information data blocks
  for (bNextBlock=(pHLTOUT->SelectFirstDataBlock()>=0);
       bNextBlock; bNextBlock=(pHLTOUT->SelectNextDataBlock()>=0)) {
    AliHLTComponentBlockData desc;
    if ((iResult=pHLTOUT->GetDataBuffer(desc))<0) {
      continue;
    }
    if (desc.fDataType==AliHLTTPCDefinitions::DataCompressionDescriptorDataType()) {
      // compression header      
      if ((iResult=decoder.AddCompressionDescriptor(&desc))<0) {
	return iResult;
      }
      bHavePartitionCompressedData = true;
    }
    if (desc.fDataType==AliHLTTPCDefinitions::RawClustersDescriptorDataType()) {
      // raw clusters header      
      if ((iResult=decoder.AddRawClustersDescriptor(&desc))<0) {
	return iResult;
      }
      bHavePartitionRawData = true;
    }
    if (desc.fDataType==AliHLTTPCDefinitions::AliHLTDataTypeClusterMCInfo()) {
      // add mc information
      if ((iResult=decoder.AddClusterMCData(&desc))<0) {
	return iResult;
      }
    }
    if (desc.fDataType==AliHLTTPCDefinitions::RemainingClusterIdsDataType() ||
	desc.fDataType==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
      // add cluster ids
      if ((iResult=decoder.AddClusterIds(&desc))<0) {
	return iResult;
      }
    }
  }

  vector<bool> bHavePartitionData(216, false);

  // read data
  iResult=-ENODATA;
  int nExtractedClusters=0;
  for (bNextBlock=(pHLTOUT->SelectFirstDataBlock()>=0);
       bNextBlock; bNextBlock=(pHLTOUT->SelectNextDataBlock()>=0)) {
    decoder.SetPadShift(0.0);
    AliHLTComponentBlockData desc;
    if ((iResult=pHLTOUT->GetDataBuffer(desc))<0) {
      continue;
    }
    if (desc.fDataType==AliHLTTPCDefinitions::RawClustersDataType()) {
      // This is a special handling of data blocks produced with v5-01-Release
      // The pad shift by 0.5 was not included in the data but was applied in the
      // unpacking in this class. Changed in r51306, the next tag containing this
      // change in the online system is v5-01-Rev-07. There are only very few runs
      // of Sep 2011 with recorded clusters not containing the 0.5 shift
      // There was also a chenge in the data type of the compressed partition
      // cluster blocks which helps to identify the blocks which need the pad shift
      // here
      if (desc.fSize<sizeof(AliHLTTPCRawClusterData)) continue;
      const AliHLTTPCRawClusterData* clusterData = reinterpret_cast<const AliHLTTPCRawClusterData*>(desc.fPtr);
      if (!clusterData) continue;
      if (clusterData->fVersion==1) {
	// compressed clusters without the pad shift
	// no raw clusters (version==0) have ever been recorded
	decoder.SetPadShift(0.5);
      }
      AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(desc.fSpecification);
      AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(desc.fSpecification);
      if (slice!=AliHLTTPCDefinitions::GetMaxSliceNr(desc.fSpecification) ||
	  partition!=AliHLTTPCDefinitions::GetMaxPatchNr(desc.fSpecification)) {
	AliFatal(Form("inconsistent cluster data: can not handle blocks containing multiple partitions, "
		      "block specification 0x%08x", desc.fSpecification));
      }
      iResult=decoder.ReadClustersPartition(pContainer->BeginRemainingClusterBlock(0, desc.fSpecification),
					    reinterpret_cast<AliHLTUInt8_t*>(desc.fPtr),
					    desc.fSize,
					    desc.fSpecification);
      if (iResult>=0) nExtractedClusters+=iResult;
      else {
	AliFatal(Form("processing of cluster block 0x%08x failed with error code %d", desc.fSpecification, iResult));
      }
      unsigned index=slice*AliHLTTPCGeometry::GetNumberOfPatches()+partition;
      if (index>=bHavePartitionData.size()) bHavePartitionData.resize(index, false);
      if (bHavePartitionData[index]) {
	AliFatal(Form("inconsistent cluster data: multiple data blocks of identical specification indicate a failure "
		      "in the production of the data. Probably an HLT emulation chain is executed in the reconstruction "
		      "and produces data in addition to HLTOUT. Option 'ignore-hltout' is required in that case; "
		      "block specification 0x%08x", desc.fSpecification));
      }
      bHavePartitionData[index]=true;
      if (bHavePartitionCompressedData) {
	AliFatal(Form("inconsistent cluster data: both compressed and raw cluster blocks present in HLTOUT, indicates a failure "
		      "in the production of the data. Probably an HLT emulation chain is executed in the reconstruction "
		      "and produces data in addition to HLTOUT. Option 'ignore-hltout' is required in that case; "
		      "block specification 0x%08x", desc.fSpecification));
      }
      bHavePartitionRawData=true;
      continue;
    } else if (desc.fDataType==AliHLTTPCDefinitions::RemainingClustersCompressedDataType()) {
      AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(desc.fSpecification);
      AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(desc.fSpecification);
      if (slice!=AliHLTTPCDefinitions::GetMaxSliceNr(desc.fSpecification) ||
	  partition!=AliHLTTPCDefinitions::GetMaxPatchNr(desc.fSpecification)) {
	AliFatal(Form("inconsistent cluster data: can not handle blocks containing multiple partitions, "
		      "block specification 0x%08x", desc.fSpecification));
      }
      iResult=decoder.ReadClustersPartition(pContainer->BeginRemainingClusterBlock(0, desc.fSpecification),
					    reinterpret_cast<AliHLTUInt8_t*>(desc.fPtr),
					    desc.fSize,
					    desc.fSpecification);
      if (iResult>0) nExtractedClusters+=iResult;
      unsigned index=slice*AliHLTTPCGeometry::GetNumberOfPatches()+partition;
      if (index>=bHavePartitionData.size()) bHavePartitionData.resize(index, false);
      if (bHavePartitionData[index]) {
	AliFatal(Form("inconsistent cluster data: multiple data blocks of identical specification indicate a failure "
		      "in the production of the data. Probably an HLT emulation chain is executed in the reconstruction "
		      "and produces data in addition to HLTOUT. Option 'ignore-hltout' is required in that case; "
		      "block specification 0x%08x", desc.fSpecification));
      }
      bHavePartitionData[index]=true;
      bHavePartitionData[index]=true;
      if (bHavePartitionRawData) {
	AliFatal(Form("inconsistent cluster data: both compressed and raw cluster blocks present in HLTOUT, indicates a failure "
		      "in the production of the data. Probably an HLT emulation chain is executed in the reconstruction "
		      "and produces data in addition to HLTOUT. Option 'ignore-hltout' is required in that case; "
		      "block specification 0x%08x", desc.fSpecification));
      }
      bHavePartitionCompressedData=true;
      continue;
    } else if (desc.fDataType==AliHLTTPCDefinitions::ClusterTracksCompressedDataType()) {
      iResult=decoder.ReadTrackModelClustersCompressed(pContainer->BeginTrackModelClusterBlock(0),
							reinterpret_cast<AliHLTUInt8_t*>(desc.fPtr),
							desc.fSize,
							desc.fSpecification);
      continue;
    }
  }

  pSystem->ReleaseHLTOUT(pHLTOUT);
    
  if (iResult<0) return iResult;
  return nExtractedClusters;
}

int AliHLTTPCDataPublisherComponent::DoInit( int argc, const char** argv )
{
  /// inherited from AliHLTComponent: component initialisation and argument scan.
  int iResult=0;

  // component configuration
  //Stage 1: default initialization.
  const char* defaultArguments="-detector TPC -datatype 'DDL_RAW ' 'TPC ' -skipempty";
  if ((iResult = ConfigureFromArgumentString(1, &defaultArguments)) < 0)
    return iResult;

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
  if ((iResult=AliHLTRawReaderPublisherComponent::DoInit(0, NULL))<0)
    return iResult;

  auto_ptr<AliRawClusterContainer> container(new AliRawClusterContainer);
  if (!container.get()) return -ENOMEM;

  fClusters=container.release();

  return iResult;
}

int AliHLTTPCDataPublisherComponent::DoDeinit()
{
  /// inherited from AliHLTComponent: component cleanup
  int iResult=0;

  if (fpDecoder) delete fpDecoder;
  fpDecoder=NULL;

  return iResult;
}

int AliHLTTPCDataPublisherComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  /// inherited from AliHLTComponent: argument scan
  if (argc<1) return 0;
  int bMissingParam=0;
  int i=0;
  TString argument=argv[i];

  do {
    // -publish-raw
    if (argument.CompareTo("-publish-raw")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter=argv[i];
      if (parameter.CompareTo("all")==0) {
	fMode|=kPublishRawAll;
	return 2;
      } else if (parameter.CompareTo("filtered")==0) {
	fMode|=kPublishRawFiltered;
	fMode|=kRegisterClusterBlocks;
	fMode&=~kPublishRawAll;
	return 2;
      } else if (parameter.CompareTo("off")==0) {
	fMode&=~(kPublishRawAll|kPublishRawFiltered);
	return 2;
      } else {
	HLTError("invalid parameter for argument %s, expecting either 'all', 'filtered', or 'off' instead of %s", argument.Data(), parameter.Data());
	return -EPROTO;
      }
    }
    // -publish-clusters
    if (argument.CompareTo("-publish-clusters")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter=argv[i];
      if (parameter.CompareTo("all")==0) {
	fMode|=kPublishClustersAll;
	return 2;
      } else if (parameter.CompareTo("off")==0) {
	fMode&=~(kPublishClustersAll);
	return 2;
      } else {
	HLTError("invalid parameter for argument %s, expecting either 'all', or 'off' instead of %s", argument.Data(), parameter.Data());
	return -EPROTO;
      }
    }

  } while (0); // using do-while only to have break available

  return AliHLTRawReaderPublisherComponent::ScanConfigurationArgument(argc, argv);
}

int AliHLTTPCDataPublisherComponent::GetSpecificationFromEquipmentId(int id, AliHLTUInt32_t &specification) const
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

bool AliHLTTPCDataPublisherComponent::IsSelected(int equipmentId) const
{
  /// inherited from AliHLTRawReaderPublisherComponent: check if a block is selected or not
  /// check if a raw data block needs to be published. This is the case if
  /// there is no corresponding compressed data, i.e. function returns
  /// only false if the block can be found in the cluster container
  if (CheckMode(kPublishRawAll))
    return true;
  if (!CheckMode(kPublishRawFiltered))
    return false;

  if (!fClusters)
    return true;

  int offset=AliHLTDAQ::DdlIDOffset(3);
  int count=AliHLTDAQ::NumberOfDdls(3);
  if (offset<0 || count<0)
    return true;
  if (equipmentId<offset)
    return true;
  equipmentId-=offset;
  if (equipmentId>=count)
    return true;
  int slice=equipmentId<72?equipmentId/2:(equipmentId-72)/4;
  int partition=equipmentId<72?equipmentId%2:((equipmentId-72)%4)+2;
  AliHLTUInt32_t specification=AliHLTTPCDefinitions::EncodeDataSpecification(slice, slice, partition, partition);
  for (AliHLTComponentBlockDataList::const_iterator i=fClusters->GetBlockDescriptors().begin();
       i!=fClusters->GetBlockDescriptors().end(); i++) {
    if (i->fSpecification==specification)
      return false;
  }
  return true;
}

AliHLTTPCDataPublisherComponent::AliRawClusterContainer::AliRawClusterContainer()
  : AliHLTLogging()
  , fBlockCount(0)
  , fTotalClusterCount(0)
  , fBlockClusterCount(0)
  , fpBuffer(NULL)
  , fBufferSize(0)
  , fDescriptors()
  , fCurrentBlock(NULL)
  , fTrackModelClusters(NULL)
  , fTrackModelClusterMap()
  , fIterator()
  , fState(0)
{
  // constructor
}

AliHLTTPCDataPublisherComponent::AliRawClusterContainer::~AliRawClusterContainer()
{
  // destructor
}

int AliHLTTPCDataPublisherComponent::AliRawClusterContainer::SetTargetBuffer(AliHLTUInt8_t* pBuffer, int size)
{
  // set/reset the external target buffer
  Clear();
  fpBuffer=pBuffer;
  fBufferSize=pBuffer?size:0;
  return 0;  
}

int AliHLTTPCDataPublisherComponent::AliRawClusterContainer::Sort()
{
  // merge track model clusters into partition cluster blocks

  // TODO: implement merging
  // decoding of track model clusters needs to be done after all
  // partition blocks have been decoded. The track model clusters are
  // then at the end of the target buffer and have to be sorted into the
  // other blocks
  // 1) move track model cluster block by its own size back in buffer
  //    if not enough space, allocate temporary buffer and increase the
  //    size estimator for the next event
  // 2) fill the index grid
  // 3) make appropriate gaps between the partition cluster blocks
  // 4) copy clusters into the partitions and update descriptors
  return -ENOSYS;
}

int AliHLTTPCDataPublisherComponent::AliRawClusterContainer::CopyBlockDescriptors(AliHLTComponentBlockDataList& target) const
{
  // fill block descriptors of extracted partition cluster blocks to target list
  target.insert(target.begin(), fDescriptors.begin(), fDescriptors.end());
  return fDescriptors.size();
}

AliHLTTPCDataPublisherComponent::AliRawClusterContainer::iterator& AliHLTTPCDataPublisherComponent::AliRawClusterContainer::BeginPartitionClusterBlock(int count, AliHLTUInt32_t specification)
{
  /// iterator of partition clusters block of specification
  return ClusterIterator(count, AliHLTTPCDefinitions::RemainingClustersCompressedDataType(), specification, fCurrentBlock);
}

AliHLTTPCDataPublisherComponent::AliRawClusterContainer::iterator& AliHLTTPCDataPublisherComponent::AliRawClusterContainer::BeginTrackModelClusterBlock(int count)
{
  /// iterator of track model clusters
  return ClusterIterator(count, AliHLTTPCDefinitions::ClusterTracksCompressedDataType(), 0x23000500, fTrackModelClusters);
}

AliHLTTPCDataPublisherComponent::AliRawClusterContainer::iterator& AliHLTTPCDataPublisherComponent::AliRawClusterContainer::ClusterIterator(int /*count*/, AliHLTComponentDataType dt, AliHLTUInt32_t specification, AliHLTTPCRawClusterData* &pData)
{
  /// iterator of partition clusters block of specification
  fBlockCount++;
  fIterator.~iterator();
  fCurrentBlock=NULL;
  fTrackModelClusters=NULL;
  fTrackModelClusterMap.clear();
  fBlockClusterCount=0;
  AliHLTUInt32_t filled=0;
  for (AliHLTComponentBlockDataList::const_iterator desc=fDescriptors.begin();
       desc!=fDescriptors.end(); desc++) {
    filled+=desc->fSize;
    if (desc->fSpecification==specification &&
	desc->fDataType==dt) {
      HLTFatal("partition cluster block with data type %s and specification 0x%08x has been already processed",
	       AliHLTComponent::DataType2Text(dt).c_str(), specification);
      filled=fBufferSize;
    }
  }

  // insert an empty data block which is than updated later
  AliHLTComponentBlockData bd;
  AliHLTComponent::FillBlockData(bd);
  bd.fPtr=NULL;
  bd.fSize=0;
  bd.fOffset=filled;
  bd.fDataType=dt;
  bd.fSpecification=specification;
  fDescriptors.push_back(bd);

  // initialize only the header, during filling the cluster count of the header
  // and the block size will be incremented
  AliHLTUInt32_t blocksize=sizeof(AliHLTTPCRawClusterData);
  if (filled+blocksize>(unsigned)fBufferSize || fpBuffer==NULL) {
    new (&fIterator) iterator(this);
    return fIterator;
  }
  pData=reinterpret_cast<AliHLTTPCRawClusterData*>(fpBuffer+filled);
  pData->fVersion=0;
  pData->fCount=0;
  fDescriptors.back().fSize=blocksize;
  new (&fIterator) iterator(this);
  return fIterator;
}

AliHLTTPCRawCluster* AliHLTTPCDataPublisherComponent::AliRawClusterContainer::NextCluster(int slice, int partition)
{
  /// increment to next cluster
  fTotalClusterCount++;
  fBlockClusterCount++;
  if (!fCurrentBlock && !fTrackModelClusters)
    return NULL;
  if (fDescriptors.size()==0)
    return NULL;
  AliHLTTPCRawClusterData* data=fCurrentBlock?fCurrentBlock:fTrackModelClusters;
  if (int(fDescriptors.back().fOffset+fDescriptors.back().fSize+sizeof(AliHLTTPCRawCluster))>=fBufferSize) {
    fState=-ENOSPC;
    return NULL;
  }
  data->fCount++;
  fDescriptors.back().fSize+=sizeof(AliHLTTPCRawCluster);
  if (fTrackModelClusters)
    fTrackModelClusterMap.push_back(AliHLTTPCSpacePointData::GetID(slice, partition, fBlockClusterCount));
  return data->fClusters+(data->fCount-1);
}

void  AliHLTTPCDataPublisherComponent::AliRawClusterContainer::Clear(Option_t * /*option*/)
{
  /// internal cleanup
  fBlockCount=0;
  fTotalClusterCount=0;
  fBlockClusterCount=0;
  fpBuffer=NULL;
  fBufferSize=0;
  fCurrentBlock=NULL;
  fTrackModelClusters=NULL;
  fTrackModelClusterMap.clear();
  fDescriptors.clear();
  fState=0;
}

void AliHLTTPCDataPublisherComponent::AliRawClusterContainer::Print(Option_t */*option*/) const
{
  /// print info
}

AliHLTTPCDataPublisherComponent::AliRawClusterContainer::iterator& AliHLTTPCDataPublisherComponent::AliRawClusterContainer::iterator::Next(int slice, int partition)
{
  // increment iterator
  if (fContainer) {
    fCluster=fContainer->NextCluster(slice, partition);
    if (fCluster) memset(fCluster, 0, sizeof(AliHLTTPCRawCluster));
  } else {
    fCluster=NULL;
  }
  return *this;
}
