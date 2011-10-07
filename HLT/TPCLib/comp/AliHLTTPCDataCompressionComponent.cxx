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

/// @file   AliHLTTPCDataCompressionComponent.cxx
/// @author Matthias Richter
/// @date   2011-08-08
/// @brief  TPC component for data compression
///

#include "AliHLTTPCDataCompressionComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCTrackGeometry.h"
#include "AliHLTTPCSpacePointContainer.h"
#include "AliHLTTPCHWCFSpacePointContainer.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTComponentBenchmark.h"
#include "AliHLTDataDeflaterSimple.h"
#include "AliHLTDataDeflaterHuffman.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusterMCData.h"
#include "AliHLTTPCClusterTransformation.h"
#include "AliRawDataHeader.h"
#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "TH1F.h"
#include "TFile.h"
#include <memory>

ClassImp(AliHLTTPCDataCompressionComponent)

AliHLTTPCDataCompressionComponent::AliHLTTPCDataCompressionComponent()
  : AliHLTProcessor()
  , fMode(0)
  , fDeflaterMode(0)
  , fVerificationMode(0)
  , fMaxDeltaPad(AliHLTTPCDefinitions::GetMaxClusterDeltaPad())
  , fMaxDeltaTime(AliHLTTPCDefinitions::GetMaxClusterDeltaTime())
  , fRawInputClusters(NULL)
  , fInputClusters(NULL)
  , fTrackGrid(NULL)
  , fSpacePointGrid(NULL)
  , fpDataDeflater(NULL)
  , fHistoCompFactor(NULL)
  , fHistoResidualPad(NULL)
  , fHistoResidualTime(NULL)
  , fHistoClustersOnTracks(NULL)
  , fHistoClusterRatio(NULL)
  , fHistoTrackClusterRatio(NULL)
  , fHistogramFile()
  , fTrainingTableOutput()
  , fpBenchmark(NULL)
  , fpWrittenAssociatedClusterIds(NULL)
  , fDriftTimeFactorA(1.)
  , fDriftTimeOffsetA(0.)
  , fDriftTimeFactorC(1.)
  , fDriftTimeOffsetC(0.)
  , fVerbosity(0)
{
}

AliHLTTPCDataCompressionComponent::~AliHLTTPCDataCompressionComponent()
{
  /// destructor
  if (fpWrittenAssociatedClusterIds) delete fpWrittenAssociatedClusterIds;
}


const char* AliHLTTPCDataCompressionComponent::GetComponentID()
{
  /// inherited from AliHLTComponent: id of the component
  return "TPCDataCompressor";
}


void AliHLTTPCDataCompressionComponent::GetInputDataTypes( AliHLTComponentDataTypeList& tgtList)
{
  /// inherited from AliHLTComponent: list of data types in the vector reference
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkHWClustersDataType);
  tgtList.push_back(AliHLTTPCDefinitions::fgkClustersDataType);
  tgtList.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
}

AliHLTComponentDataType AliHLTTPCDataCompressionComponent::GetOutputDataType()
{
  /// inherited from AliHLTComponent: output data type of the component.
  return kAliHLTMultipleDataType;
}

int AliHLTTPCDataCompressionComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  /// inherited from AliHLTComponent: multiple output data types of the component.
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::RawClustersDataType());
  tgtList.push_back(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
  tgtList.push_back(AliHLTTPCDefinitions::RemainingClusterIdsDataType());
  tgtList.push_back(AliHLTTPCDefinitions::ClusterTracksCompressedDataType());
  tgtList.push_back(AliHLTTPCDefinitions::ClusterIdTracksDataType());
  return tgtList.size();
}

void AliHLTTPCDataCompressionComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  /// inherited from AliHLTComponent: output data size estimator
  constBase=0;
  inputMultiplier=1.;  // there should not be more data than input
  inputMultiplier+=.3; // slightly more data when using the old HWCF data with 20 Byte and raw clusters 22 Byte
  if (fpWrittenAssociatedClusterIds) inputMultiplier+=.3; // space for optional cluster id array
}

AliHLTComponent* AliHLTTPCDataCompressionComponent::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTTPCDataCompressionComponent;
}

int AliHLTTPCDataCompressionComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, 
						const AliHLTComponentBlockData* /*inputBlocks*/, 
						AliHLTComponentTriggerData& /*trigData*/,
						AliHLTUInt8_t* outputPtr,
						AliHLTUInt32_t& size,
						AliHLTComponentBlockDataList& outputBlocks )
{
  /// inherited from AliHLTProcessor: data processing
  int iResult=0;
  AliHLTUInt32_t capacity=size;
  size=0;

  if (!IsDataEvent()) return 0;

  if (!fRawInputClusters) {
    return -ENODEV;
  }

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->StartNewEvent();
    GetBenchmarkInstance()->Start(0);
  }

  // Process an event
  // Loop over all input blocks in the event
  bool bHaveMC=(GetFirstInputBlock(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo | kAliHLTDataOriginTPC))!=NULL;
  if ((bHaveMC || fVerificationMode>0) && fpWrittenAssociatedClusterIds==NULL) {
    fpWrittenAssociatedClusterIds=new vector<AliHLTUInt32_t>;
  }

  const AliHLTComponentBlockData* pDesc=NULL;

  AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
  AliHLTUInt32_t inputRawClusterSize=0;
  AliHLTUInt32_t outputDataSize=0;
  int allClusters=0;
  int associatedClusters=0;

  /// input track array
  vector<AliHLTGlobalBarrelTrack> inputTrackArray;

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->Start(2);
  }

  // transformed clusters
  if (fMode==10) { // FIXME: condition to be adjusted
    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType);
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      if (GetBenchmarkInstance()) {
	GetBenchmarkInstance()->AddInput(pDesc->fSize);
      }
      AliHLTUInt8_t slice = 0;
      AliHLTUInt8_t patch = 0;
      slice = AliHLTTPCDefinitions::GetMinSliceNr( pDesc->fSpecification );
      patch = AliHLTTPCDefinitions::GetMinPatchNr( pDesc->fSpecification );
      if ( minSlice==0xFF || slice<minSlice )	minSlice = slice;
      if ( maxSlice==0xFF || slice>maxSlice )	maxSlice = slice;
      if ( minPatch==0xFF || patch<minPatch )	minPatch = patch;
      if ( maxPatch==0xFF || patch>maxPatch )	maxPatch = patch;
      if (fInputClusters) {
	fInputClusters->AddInputBlock(pDesc);
      }
    }
    if (GetBenchmarkInstance()) {
      GetBenchmarkInstance()->Stop(2);
      GetBenchmarkInstance()->Start(3);
    }
  }

  // track data input
  if (fMode==2) {
    for (pDesc=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      if (GetBenchmarkInstance()) {
	GetBenchmarkInstance()->AddInput(pDesc->fSize);
      }
      AliHLTUInt8_t slice = 0;
      AliHLTUInt8_t patch = 0;
      slice = AliHLTTPCDefinitions::GetMinSliceNr( pDesc->fSpecification );
      patch = AliHLTTPCDefinitions::GetMinPatchNr( pDesc->fSpecification );
      if ( minSlice==0xFF || slice<minSlice )	minSlice = slice;
      if ( maxSlice==0xFF || slice>maxSlice )	maxSlice = slice;
      if ( minPatch==0xFF || patch<minPatch )	minPatch = patch;
      if ( maxPatch==0xFF || patch>maxPatch )	maxPatch = patch;
      const AliHLTTracksData* pTracks=reinterpret_cast<const AliHLTTracksData*>(pDesc->fPtr);
      if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(pTracks, pDesc->fSize, inputTrackArray))<0) {
	return iResult;
      }
    }
  }

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->Stop(3);
    GetBenchmarkInstance()->Start(4);
  }

  // processing
  for (vector<AliHLTGlobalBarrelTrack>::iterator track=inputTrackArray.begin();
       track!=inputTrackArray.end();
       track++) {
    int trackID=track->GetID();
    if (trackID<0) {
      // FIXME: error guard
      HLTError("invalid track ID");
      continue;
    }

    if (fVerbosity>0) {
      UInt_t nofPoints=track->GetNumberOfPoints();
      const UInt_t* points=track->GetPoints();
      for (unsigned i=0; i<nofPoints; i++) {
	int slice=AliHLTTPCSpacePointData::GetSlice(points[i]);
	int partition=AliHLTTPCSpacePointData::GetPatch(points[i]);
	int number=AliHLTTPCSpacePointData::GetNumber(points[i]);
	HLTInfo("track %d point %d id 0x%08x slice %d partition %d number %d", track->GetID(), i, points[i], slice, partition, number);
      }
    }

    AliHLTTPCTrackGeometry* trackpoints=new AliHLTTPCTrackGeometry;
    if (!trackpoints) continue;
    trackpoints->InitDriftTimeTransformation(fDriftTimeFactorA, fDriftTimeOffsetA, fDriftTimeFactorC, fDriftTimeOffsetC);
    trackpoints->SetTrackId(trackID);
    trackpoints->CalculateTrackPoints(*track);
    trackpoints->RegisterTrackPoints(fTrackGrid);
    track->SetTrackGeometry(trackpoints);
  }

  for (vector<AliHLTGlobalBarrelTrack>::const_iterator track=inputTrackArray.begin();
       track!=inputTrackArray.end();
       track++) {
    AliHLTTrackGeometry* trackpoints=track->GetTrackGeometry();
    if (!trackpoints) continue;
    trackpoints->FillTrackPoints(fTrackGrid);
  }
  if (fVerbosity>0) {
    fTrackGrid->Print();
  }

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->Stop(4);
    GetBenchmarkInstance()->Start(5);
  }

  // loop over raw cluster blocks, assign to tracks and write
  // unassigned clusters
  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkHWClustersDataType);
       pDesc!=NULL; pDesc=GetNextInputBlock()) {
    if (pDesc->fSize<=sizeof(AliRawDataHeader)) continue;
    if (GetBenchmarkInstance()) {
      GetBenchmarkInstance()->Start(1);
      GetBenchmarkInstance()->AddInput(pDesc->fSize);
    }
    AliHLTUInt8_t slice = 0;
    AliHLTUInt8_t patch = 0;
    slice = AliHLTTPCDefinitions::GetMinSliceNr( pDesc->fSpecification );
    patch = AliHLTTPCDefinitions::GetMinPatchNr( pDesc->fSpecification );
    if ( minSlice==0xFF || slice<minSlice )	minSlice = slice;
    if ( maxSlice==0xFF || slice>maxSlice )	maxSlice = slice;
    if ( minPatch==0xFF || patch<minPatch )	minPatch = patch;
    if ( maxPatch==0xFF || patch>maxPatch )	maxPatch = patch;
    inputRawClusterSize+=pDesc->fSize;

    // add the data and populate the index grid
    fRawInputClusters->AddInputBlock(pDesc);
    fRawInputClusters->PopulateAccessGrid(fSpacePointGrid, pDesc->fSpecification);
    if (fVerbosity>0 && fSpacePointGrid->GetNumberOfSpacePoints()>0) {
      HLTInfo("index grid slice %d partition %d", slice, patch);
      fSpacePointGrid->Print();
      for (AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid::iterator& cl=fSpacePointGrid->begin();
	   cl!=fSpacePointGrid->end(); cl++) {
	AliHLTUInt32_t id=cl.Data().fId;
	float row=fRawInputClusters->GetX(id);
	float pad=fRawInputClusters->GetY(id);
	float time=fRawInputClusters->GetZ(id);
	HLTInfo("    cluster id 0x%08x: row %f  pad %f  time %f", id, row, pad, time);
      }
    }
    if (GetBenchmarkInstance()) {
      GetBenchmarkInstance()->Stop(1);
      GetBenchmarkInstance()->Start(4);
    }

    // process the clusters per padrow and check the track grid
    // for tracks crossing that particular padrow
    if (GetBenchmarkInstance()) {
      GetBenchmarkInstance()->Stop(4);
      GetBenchmarkInstance()->Start(5);
    }
    allClusters+=fSpacePointGrid->GetNumberOfSpacePoints();
    iResult=ProcessTrackClusters(&inputTrackArray[0], inputTrackArray.size(), fTrackGrid, fSpacePointGrid, fRawInputClusters, slice, patch);
    int assignedInThisPartition=0;
    if (iResult>=0) {
      assignedInThisPartition=iResult;
      associatedClusters+=iResult;
    }
    iResult=ProcessRemainingClusters(&inputTrackArray[0], inputTrackArray.size(), fTrackGrid, fSpacePointGrid, fRawInputClusters, slice, patch);
    if (iResult>=0) {
      if (fSpacePointGrid->GetNumberOfSpacePoints()>0) {
	if (fVerbosity>0) HLTInfo("associated %d (%d) of %d clusters in slice %d partition %d", iResult+assignedInThisPartition, assignedInThisPartition, fSpacePointGrid->GetNumberOfSpacePoints(), slice, patch);
      }
      associatedClusters+=iResult;
    }

    // write all remaining clusters not yet assigned to tracks
    // the index grid is used to write sorted in padrow
    // FIXME: decoder index instead of data specification to be used
    // use an external access grid to reduce allocated memory
    // set to NULL after writing the clusters
    const char* writeoptions="";
    if (fpWrittenAssociatedClusterIds) {
      writeoptions="write-cluster-ids";
    }
    fRawInputClusters->SetSpacePointPropertyGrid(pDesc->fSpecification, fSpacePointGrid);
    iResult=fRawInputClusters->Write(outputPtr+size, capacity-size, outputBlocks, fpDataDeflater, writeoptions);
    fRawInputClusters->SetSpacePointPropertyGrid(pDesc->fSpecification, NULL);
    if (iResult>=0) {
      size+=iResult;
      outputDataSize+=iResult;
      // the size of the optional cluster id array must be subtracted
      if (fpWrittenAssociatedClusterIds && outputBlocks.size()>0 &&
	  outputBlocks.back().fDataType==AliHLTTPCDefinitions::AliHLTDataTypeClusterMCInfo()) {
	outputDataSize-=outputBlocks.back().fSize;
      }
      if (GetBenchmarkInstance()) GetBenchmarkInstance()->AddOutput(iResult);
    }
    if (GetBenchmarkInstance()) {
      GetBenchmarkInstance()->Stop(5);
    }

    fSpacePointGrid->Clear();
  }
  if (fHistoClusterRatio && allClusters>0) {
    if (fVerbosity>0) HLTInfo("associated %d of %d clusters to tracks", associatedClusters, allClusters);
    float ratio=associatedClusters; ratio/=allClusters;
    fHistoClusterRatio->Fill(ratio);
  }

  // output of track model clusters
  if (iResult>=0) do {
    AliHLTUInt32_t tracksBufferOffset=sizeof(AliHLTTPCTrackModelBlock);
    if (capacity-size<tracksBufferOffset) {
      iResult=-ENOSPC;
      break;
    }
    if (fpWrittenAssociatedClusterIds) fpWrittenAssociatedClusterIds->clear();
    AliHLTTPCTrackModelBlock* trackModelBlock=reinterpret_cast<AliHLTTPCTrackModelBlock*>(outputPtr+size);
    trackModelBlock->fVersion=1;
    trackModelBlock->fDeflaterMode=fpDataDeflater?fpDataDeflater->GetDeflaterVersion():0;
    trackModelBlock->fTrackCount=inputTrackArray.size();
    trackModelBlock->fClusterCount=0;
    trackModelBlock->fGlobalParameterCnt=5;
    tracksBufferOffset+=trackModelBlock->fGlobalParameterCnt*sizeof(trackModelBlock->fGlobalParameters);
    if (capacity-size<tracksBufferOffset) {
      iResult=-ENOSPC;
      break;
    }

    AliHLTUInt32_t parameterIndex=0;
    trackModelBlock->fGlobalParameters[parameterIndex++]=GetBz();
    trackModelBlock->fGlobalParameters[parameterIndex++]=fDriftTimeFactorA;
    trackModelBlock->fGlobalParameters[parameterIndex++]=fDriftTimeOffsetA;
    trackModelBlock->fGlobalParameters[parameterIndex++]=fDriftTimeFactorC;
    trackModelBlock->fGlobalParameters[parameterIndex++]=fDriftTimeOffsetC;
    if (parameterIndex!=trackModelBlock->fGlobalParameterCnt) {
      HLTError("internal error, size of parameter array has changed without providing all values");
      iResult=-EFAULT;
      break;
    }

    if (fMode==2) {
    iResult=WriteTrackClusters(inputTrackArray, fRawInputClusters, fpDataDeflater, outputPtr+size+tracksBufferOffset, capacity-size-tracksBufferOffset);
    if (iResult>=0) {
      AliHLTComponent_BlockData bd;
      FillBlockData(bd);
      bd.fOffset        = size;
      bd.fSize          = tracksBufferOffset+iResult;
      bd.fDataType      = AliHLTTPCDefinitions::ClusterTracksCompressedDataType();
      bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(minSlice, maxSlice, minPatch, maxPatch);
      outputBlocks.push_back(bd);
      size += bd.fSize;
      outputDataSize+=bd.fSize;
      HLTBenchmark("track data block of %d tracks: size %d", inputTrackArray.size(), bd.fSize);

      if (fpWrittenAssociatedClusterIds && fpWrittenAssociatedClusterIds->size()>0) {
	AliHLTComponent::FillBlockData(bd);
	bd.fOffset        = size;
	bd.fSize        = fpWrittenAssociatedClusterIds->size()*sizeof(vector<AliHLTUInt32_t>::value_type);
	if (capacity-size>bd.fSize) {
	  memcpy(outputPtr+bd.fOffset, &(*fpWrittenAssociatedClusterIds)[0], bd.fSize);
	  bd.fDataType    = AliHLTTPCDefinitions::ClusterIdTracksDataType();
	  bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(minSlice, maxSlice, minPatch, maxPatch);
	  outputBlocks.push_back(bd);    
	  size += bd.fSize;
	} else {
	  iResult=-ENOSPC;
	}
	
	fpWrittenAssociatedClusterIds->clear();
      }
    }
    }
  } while (0);

  fRawInputClusters->Clear();

  float compressionFactor=(float)inputRawClusterSize;
  if ((outputDataSize)>0) compressionFactor/=outputDataSize;
  else compressionFactor=0.;
  if (fHistoCompFactor) fHistoCompFactor->Fill(compressionFactor);

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->Stop(0);
    if (fDeflaterMode!=3) {
      HLTBenchmark("%s - compression factor %.2f", GetBenchmarkInstance()->GetStatistics(), compressionFactor);
    } else {
      HLTBenchmark("%s", GetBenchmarkInstance()->GetStatistics());
    }
  }

  if (fInputClusters) {
    fInputClusters->Clear();
  }
  if (fRawInputClusters) {
    fRawInputClusters->Clear();
  }
  if (fTrackGrid) {
    fTrackGrid->Clear();
  }

  // forward MC labels
  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo | kAliHLTDataOriginTPC);
       pDesc!=NULL; pDesc=GetNextInputBlock()) {
    outputBlocks.push_back(*pDesc);
  }

  return iResult;
}

int AliHLTTPCDataCompressionComponent::ProcessTrackClusters(AliHLTGlobalBarrelTrack* pTracks, unsigned nofTracks,
							    AliHLTTrackGeometry::AliHLTTrackGrid* pTrackIndex,
							    AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pClusterIndex,
							    AliHLTSpacePointContainer* pClusters,
							    int slice, int partition) const
{
  // process to assigned track clusters
  int assignedClusters=0;
  if (!pTracks || nofTracks==0) return 0;

  vector<AliHLTUInt32_t> processedTracks;
  for (AliHLTTrackGeometry::AliHLTTrackGrid::iterator& trackId=pTrackIndex->begin(slice, partition, -1);
       trackId!=pTrackIndex->end(); trackId++) {
    if (find(processedTracks.begin(), processedTracks.end(), trackId.Data())!=processedTracks.end()) {
      continue;
    }
    unsigned trackindex=0;
    for (; trackindex<nofTracks; trackindex++) {
      if ((unsigned)pTracks[trackindex].GetID()==trackId.Data()) break;
    }
    if (trackindex>=nofTracks) {
      HLTError("can not find track of id %d", trackId.Data());
      continue;
    }
    processedTracks.push_back(trackId.Data());
    AliHLTGlobalBarrelTrack& track=pTracks[trackindex];
    if (!track.GetTrackGeometry()) {
      HLTError("can not find track geometry for track %d", trackId.Data());
      continue;
    }
    AliHLTTPCTrackGeometry* pTrackPoints=dynamic_cast<AliHLTTPCTrackGeometry*>(track.GetTrackGeometry());
    if (!pTrackPoints) {
      HLTError("invalid track geometry type for track %d, expecting AliHLTTPCTrackGeometry", trackId.Data());
      continue;	
    }

    UInt_t nofTrackPoints=track.GetNumberOfPoints();
    const UInt_t* trackPoints=track.GetPoints();
    for (unsigned i=0; i<nofTrackPoints; i++) {
      const AliHLTUInt32_t& clusterId=trackPoints[i];
      if (AliHLTTPCSpacePointData::GetSlice(clusterId)!=(unsigned)slice ||
	  AliHLTTPCSpacePointData::GetPatch(clusterId)!=(unsigned)partition) {
	// not in the current partition;
	continue;
      }
	  
      AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid::iterator& cl=pClusterIndex->find(AliHLTSpacePointContainer::AliHLTSpacePointProperties(clusterId));
      if (cl==pClusterIndex->end()) {
	HLTError("can not find cluster no 0x%08x of track %d in index grid", clusterId, track.GetID());
	continue;
      }

      int clusterrow=(int)pClusters->GetX(clusterId);
      float clusterpad=pClusters->GetY(clusterId);
      float clustertime=pClusters->GetZ(clusterId);

      AliHLTUInt32_t pointId=AliHLTTPCSpacePointData::GetID(slice, partition, clusterrow);
      AliHLTTrackGeometry::AliHLTTrackPoint* point=pTrackPoints->GetRawTrackPoint(pointId);
      if (!point) {
	//HLTError("can not find track point slice %d partition %d padrow %d (0x%08x) of track %d", slice, partition, clusterrow, pointId, trackId.Data());
	continue;
      }
      float pad=point->GetU();
      float time=point->GetV();
      if (TMath::Abs(clusterpad-pad)<fMaxDeltaPad &&
	  TMath::Abs(clustertime-time)<fMaxDeltaTime) {
	// add this cluster to the track point and mark in the index grid
	cl.Data().fTrackId=track.GetID();
	point->AddAssociatedSpacePoint(clusterId, clusterpad-pad, clustertime-time);
	assignedClusters++;
      }
    }
  }
  return assignedClusters;
}

int AliHLTTPCDataCompressionComponent::ProcessRemainingClusters(AliHLTGlobalBarrelTrack* pTracks, unsigned nofTracks,
								AliHLTTrackGeometry::AliHLTTrackGrid* pTrackIndex,
								AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pClusterIndex,
								AliHLTSpacePointContainer* pClusters,
								int slice, int partition) const
{
  // assign remaining clusters to tracks
  int iResult=0;
  int associatedClusters=0;
  if (!pTracks || nofTracks==0) return 0;

  for (int padrow=0; padrow<AliHLTTPCTransform::GetNRows(partition); padrow++) {
    for (AliHLTTrackGeometry::AliHLTTrackGrid::iterator& trackId=pTrackIndex->begin(slice, partition, padrow);
	 trackId!=pTrackIndex->end(); trackId++) {
      unsigned i=0;
      for (; i<nofTracks; i++) {
	if ((unsigned)pTracks[i].GetID()==trackId.Data()) break;
      }
      if (i>=nofTracks) {
	HLTError("can not find track of id %d", trackId.Data());
	continue;
      }
      AliHLTGlobalBarrelTrack& track=pTracks[i];
      if (!track.GetTrackGeometry()) {
	HLTError("can not find track geometry for track %d", trackId.Data());
	continue;
      }
      AliHLTTPCTrackGeometry* pTrackPoints=dynamic_cast<AliHLTTPCTrackGeometry*>(track.GetTrackGeometry());
      if (!pTrackPoints) {
	HLTError("invalid track geometry type for track %d, expecting AliHLTTPCTrackGeometry", trackId.Data());
	continue;	
      }
      AliHLTUInt32_t pointId=AliHLTTPCSpacePointData::GetID(slice, partition, padrow);
      AliHLTTrackGeometry::AliHLTTrackPoint* point=pTrackPoints->GetRawTrackPoint(pointId);
      if (!point) {
	//HLTError("can not find track point slice %d partition %d padrow %d (0x%08x) of track %d", slice, partition, padrow, pointId, trackId.Data());
	continue;
      }
      float pad=point->GetU();
      float time=point->GetV();

      iResult=FindCellClusters(trackId.Data(), padrow, pad, time, pClusterIndex, pClusters, point);
      if (iResult>0) associatedClusters+=iResult;
      if (fVerbosity>0) {
	HLTInfo("trackpoint track %d slice %d partition %d padrow %d: %.3f \t%.3f - associated %d", track.GetID(), slice, partition, padrow, pad, time, iResult);
      }
    }
  }
  if (iResult<0) return iResult;
  return associatedClusters;
}

int AliHLTTPCDataCompressionComponent::FindCellClusters(int trackId, int padrow, float pad, float time,
							AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pClusterIndex,
							AliHLTSpacePointContainer* pClusters,
							AliHLTTrackGeometry::AliHLTTrackPoint* pTrackPoint) const
{
  // check index cell for entries and assign to track
  int count=0;
  // search a 4x4 matrix out of the 9x9 matrix around the cell addressed by
  // pad and time
  int rowindex=pClusterIndex->GetYIndex((float)padrow);
  int padindex=pClusterIndex->GetYIndex(pad);
  int timeindex=pClusterIndex->GetZIndex(time);
  if (pClusterIndex->GetCenterY(padindex)>pad) padindex--;
  if (pClusterIndex->GetCenterZ(timeindex)>pad) timeindex--;
  for (int padcount=0; padcount<2; padcount++, padindex++) {
    if (padindex<0) continue;
    if (padindex>=pClusterIndex->GetDimensionY()) break;
    for (int timecount=0; timecount<2; timecount++, timeindex++) {
      if (timeindex<0) continue;
      if (timeindex>=pClusterIndex->GetDimensionZ()) break;
      int cellindex=pClusterIndex->Index(rowindex, padindex, timeindex);
      pad=pClusterIndex->GetCenterY(cellindex);
      time=pClusterIndex->GetCenterZ(cellindex);
      for (AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid::iterator& cl=pClusterIndex->begin((float)padrow, pad, time);
	   cl!=pClusterIndex->end(); cl++) {
	if (cl.Data().fTrackId>=0) continue;
	float clusterpad=pClusters->GetY(cl.Data().fId);
	float clustertime=pClusters->GetZ(cl.Data().fId);
	if (TMath::Abs(clusterpad-pad)<fMaxDeltaPad &&
	    TMath::Abs(clustertime-time)<fMaxDeltaTime) {
	  // add this cluster to the track point and mark in the index grid
	  cl.Data().fTrackId=trackId;
	  pTrackPoint->AddAssociatedSpacePoint(cl.Data().fId, clusterpad-pad, clustertime-time);
	  count++;
	}
      }
    }
  }
  return count;
}

int AliHLTTPCDataCompressionComponent::WriteTrackClusters(const vector<AliHLTGlobalBarrelTrack>& tracks,
							  AliHLTSpacePointContainer* pSpacePoints,
							  AliHLTDataDeflater* pDeflater,
							  AliHLTUInt8_t* outputPtr,
							  AliHLTUInt32_t capacity) const
{
  // write the track data block including all associated clusters
  AliHLTUInt32_t size=0;
  for (vector<AliHLTGlobalBarrelTrack>::const_iterator track=tracks.begin();
       track!=tracks.end();
       track++) {
    if (!track->GetTrackGeometry()) {
      HLTError("can not find track geometry for track %d", track->GetID());
      return -EBADF;
    }
    AliHLTTPCTrackGeometry* pTrackPoints=dynamic_cast<AliHLTTPCTrackGeometry*>(track->GetTrackGeometry());
    if (!pTrackPoints) {
      HLTError("invalid track geometry type for track %d, expecting AliHLTTPCTrackGeometry", track->GetID());
      return -EBADF;
    }

    int result=pTrackPoints->Write(*track, pSpacePoints, pDeflater, outputPtr+size, capacity-size, fpWrittenAssociatedClusterIds);
    if (result<0) return result;
    size+=result;

    UInt_t nofTrackPoints=track->GetNumberOfPoints();
    const UInt_t* trackPoints=track->GetPoints();

    int assignedPoints=0;
    int assignedTrackPoints=0;
    const vector<AliHLTTrackGeometry::AliHLTTrackPoint>& rawPoints=pTrackPoints->GetRawPoints();
    for (vector<AliHLTTrackGeometry::AliHLTTrackPoint>::const_iterator point=rawPoints.begin();
	 point!=rawPoints.end(); point++) {
      const vector<AliHLTTrackGeometry::AliHLTTrackSpacepoint>& spacePoints=point->GetSpacepoints();
      for (vector<AliHLTTrackGeometry::AliHLTTrackSpacepoint>::const_iterator spacePoint=spacePoints.begin();
	   spacePoint!=spacePoints.end(); spacePoint++) {
	float dpad=spacePoint->GetResidual(0);
	float dtime=spacePoint->GetResidual(1);
	if (dpad>-1000 && dtime>-1000 && fHistoResidualPad && fHistoResidualTime) {
	  fHistoResidualPad->Fill(dpad);
	  fHistoResidualTime->Fill(dtime);
	}
	assignedPoints++;
	for (unsigned i=0; i<nofTrackPoints; i++) {
	  if (trackPoints[i]==spacePoint->fId) {
	    assignedTrackPoints++;
	    break;
	  }
	}
      }
    }
    if (fHistoClustersOnTracks) {
      fHistoClustersOnTracks->Fill(assignedPoints);
    }
    if (fHistoTrackClusterRatio && nofTrackPoints>0) {
      float ratio=assignedTrackPoints; ratio/=nofTrackPoints;
      fHistoTrackClusterRatio->Fill(ratio);
    }
  }
  return size;
}

int AliHLTTPCDataCompressionComponent::DoInit( int argc, const char** argv )
{
  /// inherited from AliHLTComponent: component initialisation and argument scan.
  int iResult=0;

  // component configuration
  //Stage 1: default initialization.
  //Default values.

  //Stage 2: OCDB.
  TString cdbPath("HLT/ConfigTPC/");
  cdbPath += GetComponentID();
  //
  iResult = ConfigureFromCDBTObjString(cdbPath);
  if (iResult < 0) 
    return iResult;

  //Stage 3: command line arguments.
  if (argc && (iResult = ConfigureFromArgumentString(argc, argv)) < 0)
    return iResult;

  std::auto_ptr<AliHLTComponentBenchmark> benchmark(new AliHLTComponentBenchmark);
  if (benchmark.get()) {
    benchmark->SetTimer(0,"total");
    benchmark->SetTimer(1,"rawclusterinput");
    benchmark->SetTimer(2,"clusterinput");
    benchmark->SetTimer(3,"trackinput");
    benchmark->SetTimer(4,"processing");
    benchmark->SetTimer(5,"output");
  } else {
    return -ENOMEM;
  }

  unsigned spacePointContainerMode=(fMode==2)?AliHLTTPCHWCFSpacePointContainer::kModeCreateMap:0;
  std::auto_ptr<AliHLTTPCHWCFSpacePointContainer> rawInputClusters(new AliHLTTPCHWCFSpacePointContainer(spacePointContainerMode));
  std::auto_ptr<AliHLTTPCSpacePointContainer> inputClusters(new AliHLTTPCSpacePointContainer);

  std::auto_ptr<TH1F> histoCompFactor(new TH1F("CompressionFactor",
					       "HLT TPC data compression factor",
					       100, 0., 10.));
  std::auto_ptr<TH1F> histoResidualPad(new TH1F("PadResidual",
						"HLT TPC pad residual",
						100, -fMaxDeltaPad, fMaxDeltaPad));
  std::auto_ptr<TH1F> histoResidualTime(new TH1F("TimeResidual",
						 "HLT TPC time residual",
						 100, -fMaxDeltaTime, fMaxDeltaTime));
  std::auto_ptr<TH1F> histoClustersOnTracks(new TH1F("ClustersOnTracks",
						 "Clusters in track model compression",
						 200, 0., 600));
  std::auto_ptr<TH1F> histoClusterRatio(new TH1F("ClusterRatio",
						 "Fraction of clusters in track model compression",
						 100, 0., 1.));
  std::auto_ptr<TH1F> histoTrackClusterRatio(new TH1F("UsedTrackClusters",
						 "Fraction of track clusters in track model compression",
						 100, 0., 1.));

  // track grid: 36 slices, each 6 partitions with max 33 rows
  fTrackGrid=new AliHLTTrackGeometry::AliHLTTrackGrid(36, 1, 6, 1, 33, 1, 20000);
  fSpacePointGrid=AliHLTTPCHWCFSpacePointContainer::AllocateIndexGrid();

  if (!rawInputClusters.get() ||
      !inputClusters.get() ||
      !fTrackGrid ||
      !fSpacePointGrid) {
    if (fTrackGrid) delete fTrackGrid; fTrackGrid=NULL;
    if (fSpacePointGrid) delete fSpacePointGrid; fSpacePointGrid=NULL;
    return -ENOMEM;
  }

  if (fDeflaterMode>0 && (iResult=InitDeflater(fDeflaterMode))<0)
    return iResult;

  fpBenchmark=benchmark.release();
  fRawInputClusters=rawInputClusters.release();
  fInputClusters=inputClusters.release();

  // initialize the histograms if stored at the end
  // condition might be extended
  if (!fHistogramFile.IsNull()) {
    fHistoCompFactor=histoCompFactor.release();
    fHistoResidualPad=histoResidualPad.release();
    fHistoResidualTime=histoResidualTime.release();
    fHistoClustersOnTracks=histoClustersOnTracks.release();
    fHistoClusterRatio=histoClusterRatio.release();
    fHistoTrackClusterRatio=histoTrackClusterRatio.release();
  }

  if (iResult>=0 && (iResult=InitDriftTimeTransformation())<0) return iResult;

  return iResult;
}

int AliHLTTPCDataCompressionComponent::InitDeflater(int mode)
{
  /// init the data deflater
  int iResult=0;
  if (mode==2 || mode==3) {
    // huffman deflater
    std::auto_ptr<AliHLTDataDeflaterHuffman> deflater(new AliHLTDataDeflaterHuffman(mode==3));
    if (!deflater.get()) return -ENOMEM;

    if (!deflater->IsTrainingMode()) {
      TString cdbPath("HLT/ConfigTPC/");
      cdbPath += GetComponentID();
      cdbPath += "HuffmanTables";
      TObject* pConf=LoadAndExtractOCDBObject(cdbPath);
      if (!pConf) return -ENOENT;
      if (dynamic_cast<TList*>(pConf)==NULL) {
	HLTError("huffman table configuration object of inconsistent type");
	return -EINVAL;
      }
      iResult=deflater->InitDecoders(dynamic_cast<TList*>(pConf));
      if (iResult<0) return iResult;
    }
    
    unsigned nofParameters=AliHLTTPCDefinitions::GetNumberOfClusterParameterDefinitions();
    unsigned p=0;
    for (; p<nofParameters; p++) {
      const AliHLTTPCDefinitions::AliClusterParameter& parameter=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[p];
      if (deflater->AddParameterDefinition(parameter.fName,
					   parameter.fBitLength)!=(int)parameter.fId) {
	// for performance reason the parameter id is simply used as index in the array of
	// definitions, the position must match the id
	HLTFatal("mismatch between parameter id and position in array for parameter %s, rearrange definitions!", parameter.fName);
	return -EFAULT;
      }
    }
    fpDataDeflater=deflater.release();
    return 0;
  }
  if (mode==1) {
    std::auto_ptr<AliHLTDataDeflaterSimple> deflater(new AliHLTDataDeflaterSimple);
    if (!deflater.get()) return -ENOMEM;

    unsigned nofParameters=AliHLTTPCDefinitions::GetNumberOfClusterParameterDefinitions();
    unsigned p=0;
    for (; p<nofParameters; p++) {
      const AliHLTTPCDefinitions::AliClusterParameter& parameter=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[p];
      if (deflater->AddParameterDefinition(parameter.fName,
					   parameter.fBitLength,
					   parameter.fOptional)!=(int)parameter.fId) {
	// for performance reason the parameter id is simply used as index in the array of
	// definitions, the position must match the id
	HLTFatal("mismatch between parameter id and position in array for parameter %s, rearrange definitions!", parameter.fName);
	return -EFAULT;
      }
    }
    fpDataDeflater=deflater.release();
    return 0;
  }
  HLTError("invalid deflater mode %d, allowed 1=simple 2=huffman", mode);
  return -EINVAL;
}

int AliHLTTPCDataCompressionComponent::DoDeinit()
{
  /// inherited from AliHLTComponent: component cleanup
  int iResult=0;
  if (fpBenchmark) delete fpBenchmark; fpBenchmark=NULL;
  if (fRawInputClusters) delete fRawInputClusters; fRawInputClusters=NULL;
  if (fInputClusters) delete fInputClusters; fInputClusters=NULL;
  if (!fHistogramFile.IsNull()) {
    TFile out(fHistogramFile, "RECREATE");
    if (!out.IsZombie()) {
      out.cd();
      if (fHistoCompFactor) fHistoCompFactor->Write();
      if (fHistoResidualPad) fHistoResidualPad->Write();
      if (fHistoResidualTime) fHistoResidualTime->Write();
      if (fHistoClusterRatio) fHistoClusterRatio->Write();
      if (fHistoClustersOnTracks) fHistoClustersOnTracks->Write();
      if (fHistoTrackClusterRatio) fHistoTrackClusterRatio->Write();
      out.Close();
    }
  }
  if (fHistoCompFactor) delete fHistoCompFactor;
  fHistoCompFactor=NULL;
  if (fHistoResidualPad) delete fHistoResidualPad;
  fHistoResidualPad=NULL;
  if (fHistoResidualTime) delete fHistoResidualTime;
  fHistoResidualTime=NULL;
  if (fHistoClustersOnTracks) delete fHistoClustersOnTracks;
  fHistoClustersOnTracks=NULL;
  if (fHistoClusterRatio) delete fHistoClusterRatio;
  fHistoClusterRatio=NULL;
  if (fHistoTrackClusterRatio) delete fHistoTrackClusterRatio;
  fHistoTrackClusterRatio=NULL;

  if (fpDataDeflater) {
    if (!fHistogramFile.IsNull()) {
      TString filename=fHistogramFile;
      filename.ReplaceAll(".root", "-deflater.root");
      fpDataDeflater->SaveAs(filename);
    }
    if (fDeflaterMode==3) {
      if (fTrainingTableOutput.IsNull()) {
	fTrainingTableOutput=GetComponentID();
	fTrainingTableOutput+="-huffman.root";
      }
      // TODO: currently, the code tables are also calculated in FindObject
      // check if a different function is more appropriate
      TObject* pConf=fpDataDeflater->FindObject("DeflaterConfiguration");
      if (pConf) {
	TString cdbEntryPath("HLT/ConfigTPC/");
	cdbEntryPath += GetComponentID();
	cdbEntryPath += "HuffmanTables";
	AliCDBPath cdbPath(cdbEntryPath);
	AliCDBId cdbId(cdbPath, AliCDBManager::Instance()->GetRun(), AliCDBRunRange::Infinity(), 0, 0);
	AliCDBMetaData* cdbMetaData=new AliCDBMetaData;
	cdbMetaData->SetResponsible("ALICE HLT Matthias.Richter@cern.ch");
	cdbMetaData->SetComment("Huffman encoder configuration");
	AliCDBEntry* entry=new AliCDBEntry(pConf, cdbId, cdbMetaData, kTRUE);

	entry->SaveAs(fTrainingTableOutput);
	// this is a small memory leak
	// seg fault in ROOT object handling if the two objects are deleted
	// investigate later
	//delete entry;
	//delete cdbMetaData;
      }
    }
    delete fpDataDeflater;
  }
  fpDataDeflater=NULL;


  if (fTrackGrid) delete fTrackGrid; fTrackGrid=NULL;
  if (fSpacePointGrid) delete fSpacePointGrid; fSpacePointGrid=NULL;

  return iResult;
}

int AliHLTTPCDataCompressionComponent::ScanConfigurationArgument(int argc, const char** argv)
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
	fMode=parameter.Atoi();
	return 2;
      } else {
	HLTError("invalid parameter for argument %s, expecting number instead of %s", argument.Data(), parameter.Data());
	return -EPROTO;
      }
    }

    // -deflater-mode
    if (argument.CompareTo("-deflater-mode")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter=argv[i];
      if (parameter.IsDigit()) {
	fDeflaterMode=parameter.Atoi();
	return 2;
      } else {
	HLTError("invalid parameter for argument %s, expecting number instead of %s", argument.Data(), parameter.Data());
	return -EPROTO;
      }
    }

    // -histogram-file
    if (argument.CompareTo("-histogram-file")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fHistogramFile=argv[i++];
      return 2;
    }
    // -save-histogram-table
    if (argument.CompareTo("-save-huffman-table")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fTrainingTableOutput=argv[i++];
      return 2;
    }
    // -cluster-verification
    if (argument.CompareTo("-cluster-verification")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter=argv[i];
      if (parameter.IsDigit()) {
	fVerificationMode=parameter.Atoi();
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

int AliHLTTPCDataCompressionComponent::InitDriftTimeTransformation()
{
  /// calculate correction factor and offset for a linear approximation of the
  /// drift time transformation, separately for A and C side
  int iResult=0;
  AliHLTTPCClusterTransformation transform;
  if ((iResult=transform.Init( GetBz(), GetTimeStamp()))<0) {
    HLTError("failed to init AliHLTTPCClusterTransformation: %d", iResult);
    return iResult;
  }

  if ((iResult=CalculateDriftTimeTransformation(transform, 0, 0, fDriftTimeFactorA, fDriftTimeOffsetA))<0) return iResult;
  if (fVerbosity>0) HLTInfo("drift time transformation A side: m=%f n=%f", fDriftTimeFactorA, fDriftTimeOffsetA);
  if ((iResult=CalculateDriftTimeTransformation(transform, 18, 0, fDriftTimeFactorC, fDriftTimeOffsetC))<0) return iResult;
  if (fVerbosity>0) HLTInfo("drift time transformation C side: m=%f n=%f", fDriftTimeFactorC, fDriftTimeOffsetC);

  return 0;
}

int AliHLTTPCDataCompressionComponent::CalculateDriftTimeTransformation(AliHLTTPCClusterTransformation& transform,
									int slice, int padrow,
									float& m, float& n) const
{
  /// calculate correction factor and offset for a linear approximation of the
  /// drift time transformation by just probing the range of timebins with
  /// AliHLTTPCClusterTransformation
  const int nofSteps=100;
  vector<float> zvalues;

  int nofTimebins=AliHLTTPCTransform::GetNTimeBins();
  int stepWidth=nofTimebins/nofSteps;
  int time=0;
  int count=0;
  float meanT=0.;
  float meanZ=0.;
  for (time=0; time<nofTimebins; time+=stepWidth, count++) {
    Float_t xyz[3];
    transform.Transform(slice, padrow, 0, time, xyz);
    zvalues.push_back(xyz[2]);
    meanT+=time;
    meanZ+=xyz[2];
  }
  meanT/=count;
  meanZ/=count;
  float sumTZ=.0;
  float sumT2=.0;
  time=0;
  for (vector<float>::const_iterator z=zvalues.begin();
       z!=zvalues.end(); z++, time+=stepWidth) {
    sumTZ+=(time-meanT)*((*z)-meanZ);
    sumT2+=(time-meanT)*(time-meanT);
  }
  m=sumTZ/sumT2;
  n=meanZ-m*meanT;

  return 0;
}
