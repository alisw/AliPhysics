//**************************************************************************
//* This file is property of and copyright by the                          *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@scieq.net>         *
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

/// @file   AliHLTTPCDataCompressionUnpackerComponent.cxx
/// @author Matthias Richter
/// @date   2016-03-08
/// @brief  Unpacker component for compressed TPC cluster data
///

#include "AliHLTTPCDataCompressionUnpackerComponent.h"
#include "AliHLTTPCDataCompressionDecoder.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTErrorGuard.h"
#include <stdexcept>
#include <cerrno>

ClassImp(AliHLTTPCDataCompressionUnpackerComponent)

AliHLTTPCDataCompressionUnpackerComponent::AliHLTTPCDataCompressionUnpackerComponent()
: AliHLTProcessor()
  , fpDecoder(NULL)
  , fClusterWriter(NULL)
  , fRequiredSpace(216*sizeof(AliHLTTPCRawClusterData))
  , fInputMultiplier(5.)
{
  /// constructor
}

AliHLTTPCDataCompressionUnpackerComponent::~AliHLTTPCDataCompressionUnpackerComponent()
{
  /// destructor
}

const char* AliHLTTPCDataCompressionUnpackerComponent::GetComponentID()
{
  /// inherited from AliHLTComponent: id of the component
  return "TPCDataCompressorUnpacker";
}

void AliHLTTPCDataCompressionUnpackerComponent::GetInputDataTypes( AliHLTComponentDataTypeList& list)
{
  /// inherited from AliHLTComponent: list of data types in the vector reference
  list.push_back(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
  list.push_back(AliHLTTPCDefinitions::RemainingClusterIdsDataType());
  list.push_back(AliHLTTPCDefinitions::ClusterTracksCompressedDataType());
  list.push_back(AliHLTTPCDefinitions::ClusterIdTracksDataType());
}

AliHLTComponentDataType AliHLTTPCDataCompressionUnpackerComponent::GetOutputDataType()
{
  /// inherited from AliHLTComponent: output data type of the component.
  return kAliHLTMultipleDataType;
}

int AliHLTTPCDataCompressionUnpackerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
  /// inherited from AliHLTComponent: multiple output data types of the component.
  list.push_back(AliHLTTPCDefinitions::RawClustersDataType());
  return list.size();
}

void AliHLTTPCDataCompressionUnpackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  /// inherited from AliHLTComponent: output data size estimator
  constBase=fRequiredSpace;
  inputMultiplier=fInputMultiplier;
}

AliHLTComponent* AliHLTTPCDataCompressionUnpackerComponent::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTTPCDataCompressionUnpackerComponent;
}

int AliHLTTPCDataCompressionUnpackerComponent::DoEvent( const AliHLTComponentEventData& evtData,
                                                        const AliHLTComponentBlockData* blocks,
                                                        AliHLTComponentTriggerData& trigData,
                                                        AliHLTUInt8_t* outputPtr,
                                                        AliHLTUInt32_t& size,
                                                        AliHLTComponentBlockDataList& outputBlocks )
{
  /// inherited from AliHLTProcessor: data processing
  int iResult=0;
  AliHLTUInt32_t capacity=size;
  size=0;

  AliHLTUInt32_t eventType=gkAliEventTypeUnknown;
  if (!IsDataEvent(&eventType)) {
    return iResult;
  }

  const AliHLTComponentBlockData* pDesc=NULL;

  if (fClusterWriter && fpDecoder) {
    fClusterWriter->Init(outputPtr, capacity);
    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::RemainingClusterIdsDataType());
         pDesc!=NULL; pDesc=GetNextInputBlock()) {
      fClusterWriter->AddClusterIds(pDesc);
    }

    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::ClusterIdTracksDataType());
         pDesc!=NULL; pDesc=GetNextInputBlock()) {
      fClusterWriter->AddClusterIds(pDesc);
    }

    // read data
    AliHLTTPCDataCompressionDecoder& decoder=*fpDecoder;
    decoder.Clear();

    // first unpack track model clusters into temporary map and count clusters per
    // partition.
    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::ClusterTracksCompressedDataType());
         pDesc!=NULL && iResult>=0; pDesc=GetNextInputBlock()) {
      AliClusterWriter::iterator tmit=fClusterWriter->BeginTrackModelClusterBlock(0);
      iResult=decoder.ReadTrackModelClustersCompressed(tmit,
                                               reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr),
                                               pDesc->fSize,
                                               pDesc->fSpecification);
    }
    // count track model clusters in the individual partitions to allow correct
    // memory allocation of the partition cluster blocks
    int tmClusterCount=fClusterWriter->ProcessTrackModelClusterCount();
    HLTInfo("extracted %d cluster(s) from track model compressed data", tmClusterCount);

    // now unpack partition cluster blocks
    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
         pDesc!=NULL && iResult>=0; pDesc=GetNextInputBlock()) {
      if (pDesc->fSize<=sizeof(AliHLTTPCRawClusterData)) {
        ALIHLTERRORGUARD(5, "inconsistent size of cluster data block");
        continue;
      }
      AliHLTTPCRawClusterData* clusterData=reinterpret_cast<AliHLTTPCRawClusterData*>(pDesc->fPtr);
      AliClusterWriter::iterator pcit=fClusterWriter->BeginPartitionClusterBlock(clusterData->fCount, pDesc->fSpecification);
      iResult=decoder.ReadClustersPartition(pcit,
                                            reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr),
                                            pDesc->fSize,
                                            pDesc->fSpecification);
    }

    // merge track model clusters into partition cluster blocks and copy output block
    // descriptors to list
    if (iResult>=0) {
      iResult=fClusterWriter->Finish(outputBlocks);
      if (iResult==-ENOSPC) {
        fRequiredSpace=fClusterWriter->GetRequiredSpace();
        fInputMultiplier=1.;
      } else if (iResult>=0) {
        size=iResult;
      }
    }
    fClusterWriter->Clear();
  }

  return iResult;
}

int AliHLTTPCDataCompressionUnpackerComponent::DoInit(int argc, const char** argv )
{
  /// inherited from AliHLTComponent: component initialization
  auto_ptr<AliHLTTPCDataCompressionDecoder> decoder(new AliHLTTPCDataCompressionDecoder);
  if (!decoder.get()) {
    return -ENOMEM;
  }
  decoder->SetPadrowModeLocal();

  auto_ptr<AliClusterWriter> clw(new AliClusterWriter);
  if (!clw.get())
    return -ENOMEM;

  fpDecoder=decoder.release();
  fClusterWriter=clw.release();

  return 0;
}

int AliHLTTPCDataCompressionUnpackerComponent::DoDeinit()
{
  /// inherited from AliHLTComponent: component cleanup
  if (fpDecoder) {
    fpDecoder->Clear();
    delete fpDecoder;
    fpDecoder=NULL;
  }
  if (fClusterWriter) {
    fClusterWriter->Clear();
    delete fClusterWriter;
    fClusterWriter=NULL;
  }
  return 0;
}

int AliHLTTPCDataCompressionUnpackerComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  /// inherited from AliHLTComponent: argument scan
  return 0;
}

AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::AliClusterWriter()
  : AliHLTLogging()
  , fOutputBuffer(NULL)
  , fBufferSize(0)
  , fBufferFilled(0)
  , fRequiredSpace(0)
  , fTrackModelClusters()
  , fTrackModelClusterCounts()
  , fPartitionClusterTargets()
  , fPartitionClusterIds()
  , fTrackModelClusterIds()
{
}

AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::AliClusterWriter(AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size)
  : AliHLTLogging()
  , fOutputBuffer(pBuffer)
  , fBufferSize(size)
  , fBufferFilled(0)
  , fRequiredSpace(0)
  , fTrackModelClusters()
  , fTrackModelClusterCounts()
  , fPartitionClusterTargets()
  , fPartitionClusterIds()
  , fTrackModelClusterIds()
{
}

AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::~AliClusterWriter()
{
}

AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::iterator AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::BeginTrackModelClusterBlock(int /*count*/)
{
  /// iterator of track model clusters
  /// the iterator is placed BEFORE the position of the first cluster, this is a bit
  /// counter intuitive, but for the sake of the implementation of the decoder class
  /// AliHLTTPCDataCompressionDecoder. To provide the slice and partition parameters
  /// for every cluster it calls Next with these parameters before every cluster

  /// The count parameter can later be used to optimize the memory allocation, for the
  /// moment however it looks like the total cluster count of the track model block
  /// has not been set in AliHLTTPCDataCompressionComponent.
  AliClusterIdBlock* clusterIds=NULL;
  clusterIds=&fTrackModelClusterIds;

  AliHLTTPCRawClusterMapRef ref(fTrackModelClusters);
  return iterator(ref, -1, -1, clusterIds);
}

AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::iterator AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::BeginPartitionClusterBlock(int count, AliHLTUInt32_t specification)
{
  /// iterator of partition clusters block of specification
  /// see note above concerning initial position of the iterator
  AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetSingleSliceNr(specification);
  assert(slice>=0);
  AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetSinglePatchNr(specification);
  assert(partition>=0);

  if (fTrackModelClusterCounts.find(specification)!=fTrackModelClusterCounts.end()) {
    // add clusters from track model format
    count+=fTrackModelClusterCounts[specification];
  }

  if (fPartitionBlockDescriptors.find(specification)==fPartitionBlockDescriptors.end()) {
    // reserve a new block
    AliHLTComponentBlockData bd=ReservePartitionClusterBlock(count, specification);
    if (bd.fSize>0) {
      fPartitionBlockDescriptors[specification]=bd;
    }
  } else {
    AliHLTComponentBlockData bd=fPartitionBlockDescriptors[specification];
    assert(bd.fSize>=sizeof(AliHLTTPCRawClusterData));
    AliHLTTPCRawClusterData* clusterData=reinterpret_cast<AliHLTTPCRawClusterData*>(fOutputBuffer+bd.fOffset);
    assert(clusterData->fCount>=count);
    if (count>=0 && clusterData->fCount<(unsigned)count) {
      // not clear yet how to handle
      throw std::runtime_error("extension of pre-allocated cluster buffer not implemented");
    }
  }

  AliClusterIdBlock* clusterIds=NULL;
  if (fPartitionClusterIds.find(specification)!=fPartitionClusterIds.end())
    clusterIds=&fPartitionClusterIds[specification];

  if (fPartitionClusterTargets.find(specification)==fPartitionClusterTargets.end()) {
    throw std::runtime_error("allocation of target buffer failed");
  }
  AliHLTTPCRawClusterArrayRef ref(fPartitionClusterTargets[specification], count>0?count:0);
  return iterator(ref, slice, partition, clusterIds);
}

AliHLTComponentBlockData AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::ReservePartitionClusterBlock(int count, AliHLTUInt32_t specification)
{
  AliHLTComponentBlockData bd;
  memset (&bd, 0, sizeof(AliHLTComponentBlockData));
  bd.fStructSize=sizeof(AliHLTComponentBlockData);
  bd.fDataType=AliHLTTPCDefinitions::fgkRawClustersDataType;
  bd.fSpecification=specification;
  bd.fOffset=fBufferFilled;

  assert(count>0);
  if (count<=0) return bd;

  // check if a block of the same partition has already been reserved
  assert(fPartitionClusterTargets.find(specification)==fPartitionClusterTargets.end());
  if (fPartitionClusterTargets.find(specification)!=fPartitionClusterTargets.end()) {
    return bd;
  }

  int requiredSpace=sizeof(AliHLTTPCRawClusterData) + count*sizeof(AliHLTTPCRawCluster);
  fRequiredSpace+=requiredSpace;
  if (requiredSpace + fBufferFilled > fBufferSize) {
    // not enough space for cluster block, accumulate the required space and
    // return an empty block descriptor
    return bd;
  }

  bd.fSize=requiredSpace;

  // init the cluster data header and array
  AliHLTTPCRawClusterData* clusterData=reinterpret_cast<AliHLTTPCRawClusterData*>(fOutputBuffer+fBufferFilled);
  clusterData->fVersion=0;
  clusterData->fCount=count;
  AliHLTTPCRawCluster* firstCluster=reinterpret_cast<AliHLTTPCRawCluster*>(clusterData+1);
  memset(firstCluster, 0, count*sizeof(AliHLTTPCRawCluster));
  fPartitionClusterTargets[specification]=firstCluster;

  fBufferFilled+=requiredSpace;

  return bd;
}

int AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::ProcessTrackModelClusterCount()
{
  // count track model clusters in the individual partitions
  for (std::map<AliHLTUInt32_t, AliHLTTPCRawCluster>::const_iterator mapit=fTrackModelClusters.begin();
       mapit!=fTrackModelClusters.end(); ++mapit) {
    unsigned slice=AliHLTTPCSpacePointData::GetSlice(mapit->first);
    unsigned partition=AliHLTTPCSpacePointData::GetPatch(mapit->first);
    AliHLTUInt32_t specification=AliHLTTPCDefinitions::EncodeDataSpecification(slice, slice, partition, partition);
    if (fTrackModelClusterCounts.find(specification)==fTrackModelClusterCounts.end()) {
      fTrackModelClusterCounts[specification]=0;
    }
    fTrackModelClusterCounts[specification]+=1;
  }

  return fTrackModelClusters.size();
}

void AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::Init(AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size)
{
  /// set output buffer and init for writing
  fOutputBuffer=pBuffer;
  fBufferSize=size;
  fBufferFilled=0;
  fRequiredSpace=0;
}

void AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::Clear(Option_t * /*option*/)
{
  fOutputBuffer=NULL;
  fBufferSize=0;
  fBufferFilled=0;
  fRequiredSpace=0;

  fTrackModelClusters.clear();
  fTrackModelClusterCounts.clear();
  fPartitionClusterTargets.clear();
  fPartitionClusterIds.clear();
  fPartitionBlockDescriptors.clear();
  new (&fTrackModelClusterIds) AliClusterIdBlock;
}

int AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::AddClusterIds(const AliHLTComponentBlockData* pDesc)
{
  /// add cluster id block for partition or track model clusters
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
    fTrackModelClusterIds.fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fTrackModelClusterIds.fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  if (pDesc->fDataType==AliHLTTPCDefinitions::RemainingClusterIdsDataType()) {
    fPartitionClusterIds[pDesc->fSpecification].fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fPartitionClusterIds[pDesc->fSpecification].fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  return -ENODATA;
}

int AliHLTTPCDataCompressionUnpackerComponent::AliClusterWriter::Finish(AliHLTComponentBlockDataList& outputBlocks)
{
  /// finish unpacking of clusters: merge track model clusters, copy block descriptors
  std::map<AliHLTUInt32_t, int> writtenPartitionClusters;
  for (std::map<AliHLTUInt32_t, AliHLTTPCRawCluster>::const_iterator mapit=fTrackModelClusters.begin();
       mapit!=fTrackModelClusters.end(); ++mapit) {
    unsigned slice=AliHLTTPCSpacePointData::GetSlice(mapit->first);
    unsigned partition=AliHLTTPCSpacePointData::GetPatch(mapit->first);
    AliHLTUInt32_t specification=AliHLTTPCDefinitions::EncodeDataSpecification(slice, slice, partition, partition);
    if (writtenPartitionClusters.find(specification)==writtenPartitionClusters.end()) {
      // init map for counting of written clusters
      if (fTrackModelClusterIds.fIds!=NULL) {
	writtenPartitionClusters[specification]=-1;
      } else if (fPartitionBlockDescriptors.find(specification)!=fPartitionBlockDescriptors.end()) {
	// this is a bit awkward to determine the first location for the track
	// model clusters after all partition clusters have been written
	unsigned nTotalClusters=(fPartitionBlockDescriptors[specification].fSize-sizeof(AliHLTTPCRawClusterData))/sizeof(AliHLTTPCRawCluster);
	writtenPartitionClusters[specification]=nTotalClusters-fTrackModelClusterCounts[specification];
      } else {
	writtenPartitionClusters[specification]=0;
      }
    }
    if (fPartitionClusterTargets.find(specification)==fPartitionClusterTargets.end()) {
      // partition cluster block not yet allocated
      AliHLTComponentBlockData bd=ReservePartitionClusterBlock(fTrackModelClusterCounts[specification], specification);
      if (bd.fSize>0) fPartitionBlockDescriptors[specification]=bd;
      else {
	ALIHLTERRORGUARD(5, "allocation of partition cluster buffer for writing of track model clusters failed, specification 0x%08x", specification);
	continue;
      }
      // indicate that track model clusters are copied into a flat partition buffer
      // contiguously, ignoring the clusters ids
      writtenPartitionClusters[specification]=0;
    }
    unsigned clusterNo=AliHLTTPCSpacePointData::GetNumber(mapit->first);
    if (fTrackModelClusterIds.fIds!=NULL && (writtenPartitionClusters[specification]<0)) {
      clusterNo=AliHLTTPCSpacePointData::GetNumber(fTrackModelClusterIds.fIds[clusterNo]);
    } else {
      clusterNo=writtenPartitionClusters[specification]++;
    }
    if ((fPartitionClusterTargets[specification])[clusterNo].GetCharge()!=0) {
      throw std::runtime_error(" cluster position has already been filled");
    }
    (fPartitionClusterTargets[specification])[clusterNo]=mapit->second;
  }

  unsigned totalSize=0;
  AliHLTComponentBlockDataList ol;
  for( std::map<AliHLTUInt32_t, AliHLTComponentBlockData>::iterator it = fPartitionBlockDescriptors.begin();
       it != fPartitionBlockDescriptors.end();
       ++it ) {
    ol.push_back( it->second );
    totalSize+=it->second.fSize;
  }
  if (totalSize!=fBufferFilled) {
    throw std::runtime_error("inconsistent output size");
  }

  if (fRequiredSpace>totalSize) {
    return -ENOSPC;
  }

  outputBlocks.insert(outputBlocks.begin(), ol.begin(), ol.end());
  return totalSize;
}
