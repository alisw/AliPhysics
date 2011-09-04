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
#include "AliHLTTPCClusterMCData.h"
#include "AliRawDataHeader.h"
#include "TH1F.h"
#include "TFile.h"
#include <memory>

ClassImp(AliHLTTPCDataCompressionComponent)

AliHLTTPCDataCompressionComponent::AliHLTTPCDataCompressionComponent()
  : AliHLTProcessor()
  , fMode(0)
  , fDeflaterMode(0)
  , fMaxDeltaPad(4)
  , fMaxDeltaTime(5)
  , fRawInputClusters(NULL)
  , fInputClusters(NULL)
  , fTrackGrid(NULL)
  , fSpacePointGrid(NULL)
  , fpDataDeflater(NULL)
  , fHistoCompFactor(NULL)
  , fHistogramFile()
  , fpBenchmark(NULL)
{
}

AliHLTTPCDataCompressionComponent::~AliHLTTPCDataCompressionComponent()
{
  /// destructor
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
  tgtList.push_back(AliHLTTPCDefinitions::fgkRawClustersDataType);
  return tgtList.size();
}

void AliHLTTPCDataCompressionComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  /// inherited from AliHLTComponent: output data size estimator
  constBase=0;
  inputMultiplier=1.3;
}

AliHLTComponent* AliHLTTPCDataCompressionComponent::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTTPCDataCompressionComponent;
}

int AliHLTTPCDataCompressionComponent::DoEvent( const AliHLTComponentEventData& evtData, 
						const AliHLTComponentBlockData* inputBlocks, 
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

  const AliHLTComponentBlockData* pDesc=NULL;

  AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
  AliHLTUInt32_t inputRawClusterSize=0;
  AliHLTUInt32_t outputDataSize=0;

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
  if (fMode==2) { // FIXME: condition to be adjusted
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
    AliHLTTrackGeometry* trackpoints=new AliHLTTPCTrackGeometry;
    if (!trackpoints) continue;
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
  for (vector<AliHLTGlobalBarrelTrack>::const_iterator track=inputTrackArray.begin();
       track!=inputTrackArray.end();
       track++) {
    AliHLTTrackGeometry* trackpoints=track->GetTrackGeometry();
    if (!trackpoints) continue;
    trackpoints->FillTrackPoints(fTrackGrid);
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
    // prototype not yet committed
    //ProcessTrackClusters(&inputTrackArray[0], inputTrackArray.size(), fTrackGrid, fSpacePointGrid, fRawInputClusters, slice, patch);

    // write all remaining clusters not yet assigned to tracks
    // the index grid is used to write sorted in padrow
    // FIXME: decoder index instead of data specification to be used
    // use an external access grid to reduce allocated memory
    // set to NULL after writing the clusters
    fRawInputClusters->SetSpacePointPropertyGrid(pDesc->fSpecification, fSpacePointGrid);
    iResult=fRawInputClusters->Write(outputPtr+size, capacity-size, outputBlocks, fpDataDeflater);
    fRawInputClusters->SetSpacePointPropertyGrid(pDesc->fSpecification, NULL);
    if (iResult>=0) {
      size+=iResult;
      outputDataSize+=iResult;
      if (GetBenchmarkInstance()) GetBenchmarkInstance()->AddOutput(iResult);
    }
    if (GetBenchmarkInstance()) {
      GetBenchmarkInstance()->Stop(5);
    }

    // forward MC labels
    if (bHaveMC) {
      // loop over input blocks and find MC block of current specification
      unsigned mcBlock=0;
      for (; mcBlock<evtData.fBlockCnt; mcBlock++) {
	if (inputBlocks[mcBlock].fDataType!=(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo | kAliHLTDataOriginTPC) ||
	    inputBlocks[mcBlock].fSpecification!=pDesc->fSpecification) {
	  continue;
	}
	if (size+inputBlocks[mcBlock].fSize>capacity) {
	  iResult=-ENOSPC;
	  break;
	}
	iResult=ForwardMCLabels(inputBlocks[mcBlock], fSpacePointGrid, outputPtr+size, capacity-size, size, outputBlocks);
	if (iResult>0) {
	  size+=iResult;
	}

	HLTDebug("forwarded MC data block of slice %d partition %d", slice, patch);
	break;
      }
      if (mcBlock==evtData.fBlockCnt) {
	HLTWarning("no mc data found for slice %d partition %d", slice, patch);
      }
    }

    fSpacePointGrid->Clear();
  }

  // output of track model clusters
  for (vector<AliHLTGlobalBarrelTrack>::const_iterator track=inputTrackArray.begin();
       track!=inputTrackArray.end();
       track++) {
    // prototype not yet committed
  }

  fRawInputClusters->Clear();

  float compressionFactor=(float)inputRawClusterSize;
  if ((outputDataSize)>0) compressionFactor/=outputDataSize;
  else compressionFactor=0.;
  if (fHistoCompFactor) fHistoCompFactor->Fill(compressionFactor);

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->Stop(0);
    HLTBenchmark("%s - compression factor %.2f", GetBenchmarkInstance()->GetStatistics(), compressionFactor);
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

  return iResult;
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

  std::auto_ptr<AliHLTTPCHWCFSpacePointContainer> rawInputClusters(new AliHLTTPCHWCFSpacePointContainer(2));
  std::auto_ptr<AliHLTTPCSpacePointContainer> inputClusters(new AliHLTTPCSpacePointContainer);
  std::auto_ptr<TH1F> histoCompFactor(new TH1F("factor", "HLT TPC data compression factor", 100, 0, 10));
  
  // track grid: 36 slices, each 6 partitions with max 33 rows
  fTrackGrid=new AliHLTTrackGeometry::AliHLTTrackGrid(36, 1, 6, 1, 33, 1, 20000);
  fSpacePointGrid=AliHLTTPCHWCFSpacePointContainer::AllocateIndexGrid();

  if (!rawInputClusters.get() ||
      !inputClusters.get() ||
      !histoCompFactor.get() ||
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
  fHistoCompFactor=histoCompFactor.release();

  return iResult;
}

int AliHLTTPCDataCompressionComponent::InitDeflater(int mode)
{
  /// init the data deflater
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
  } else if (mode==2) {
    // huffman deflater
    HLTError("huffman deflater to be implemented");
    return -ENOSYS; // change to 0 if implemented
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
  if (fHistoCompFactor) {
    if (!fHistogramFile.IsNull()) {
      TFile out(fHistogramFile, "RECREATE");
      if (!out.IsZombie()) {
	out.cd();
	fHistoCompFactor->Write();
	out.Close();
      }
    }
    delete fHistoCompFactor;
    fHistoCompFactor=NULL;
  }
  if (fpDataDeflater) delete fpDataDeflater; fpDataDeflater=NULL;
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
  } while (0); // using do-while only to have break available

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EPROTO;
  }

  return iResult;
}

int AliHLTTPCDataCompressionComponent::ForwardMCLabels(const AliHLTComponentBlockData& desc,
						       AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pIndex,
						       AliHLTUInt8_t* outputPtr,
						       AliHLTUInt32_t size, AliHLTUInt32_t offset,
						       vector<AliHLTComponentBlockData>& outputBlocks) const
{
  // forward the mc labels in the same order as the clusters
  // sorted according to the index grid
  if (!outputPtr || !pIndex) return -EINVAL;
  if (!desc.fPtr) return -ENODATA;
  if (size<desc.fSize) return -ENOSPC;

  int slice=AliHLTTPCDefinitions::GetMinSliceNr(desc.fSpecification);
  int part=AliHLTTPCDefinitions::GetMinPatchNr(desc.fSpecification);

  const AliHLTTPCClusterMCData* pInput = reinterpret_cast<const AliHLTTPCClusterMCData*>(desc.fPtr);
  Int_t nLabels = (Int_t) pInput->fCount;
  if (nLabels*sizeof(AliHLTTPCClusterMCLabel) + sizeof(AliHLTTPCClusterMCData) != desc.fSize) {
    HLTError("inconsistent cluster mc data block size, skipping block");
    return -EBADF;
  }
  const AliHLTTPCClusterMCLabel *pInputLabels = pInput->fLabels;

  AliHLTTPCClusterMCData* pOutput = reinterpret_cast<AliHLTTPCClusterMCData*>(outputPtr);
  AliHLTTPCClusterMCLabel *pOutputLabels = pOutput->fLabels;

  unsigned outIndex=0;
  for (AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid::iterator clusterID=pIndex->begin();
       clusterID!=pIndex->end(); clusterID++, outIndex++) {
      if ((unsigned)slice!=AliHLTTPCSpacePointData::GetSlice(clusterID.Data().fId) ||
      	  (unsigned)part!=AliHLTTPCSpacePointData::GetPatch(clusterID.Data().fId)) {
	HLTError("spacepoint index 0x%08x out of slice %d partition %d", clusterID.Data().fId, slice, part);
      }
      int index=AliHLTTPCSpacePointData::GetNumber(clusterID.Data().fId);
      pOutputLabels[outIndex]=pInputLabels[index];
  }
  if (outIndex==pInput->fCount) {
    pOutput->fCount=outIndex;
  } else {
    HLTError("failed to copy MC label data block 0x%08x: expecting %d, copied %d entries", desc.fSpecification, pInput->fCount, outIndex);
    return -EBADMSG;
  }

  AliHLTComponent_BlockData bd;
  AliHLTComponent::FillBlockData(bd);
  bd.fSize          = desc.fSize;
  bd.fOffset        = offset;
  bd.fSpecification = desc.fSpecification;
  bd.fDataType      = desc.fDataType;
  outputBlocks.push_back(bd);

  return bd.fSize;
}
