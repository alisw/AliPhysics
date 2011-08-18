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
#include "TH1F.h"
#include "TFile.h"
#include <memory>

AliHLTTPCDataCompressionComponent::AliHLTTPCDataCompressionComponent()
  : AliHLTProcessor()
  , fMode(0)
  , fDeflaterMode(0)
  , fRawInputClusters(NULL)
  , fInputClusters(NULL)
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

int AliHLTTPCDataCompressionComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, 
						const AliHLTComponentBlockData* /*blocks*/, 
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

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->StartNewEvent();
    GetBenchmarkInstance()->Start(0);
  }

  // Process an event
  // Loop over all input blocks in the event
  const AliHLTComponentBlockData* pDesc=NULL;

  AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
  AliHLTUInt32_t inputRawClusterSize=0;

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
  for (vector<AliHLTGlobalBarrelTrack>::const_iterator track=inputTrackArray.begin();
       track!=inputTrackArray.end();
       track++) {
    if (!fInputClusters) continue;
    int trackID=track->GetID();
    if (trackID<0) {
      // FIXME: error guard
      HLTError("invalid track ID");
      continue;
    }
    if ((iResult=fInputClusters->SetTrackID(trackID, track->GetPoints(), track->GetNumberOfPoints()))<0) {
      HLTError("failed to set cluster id for track %d: error %d", trackID, iResult);
      iResult=0;
      continue;
    }
    AliHLTTrackGeometry* trackpoints=new AliHLTTPCTrackGeometry;
    trackpoints->CalculateTrackPoints(const_cast<AliHLTGlobalBarrelTrack&>(*track));
    delete trackpoints;
  }

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->Stop(4);
    GetBenchmarkInstance()->Start(5);
  }

  // loop over raw cluster blocks, assign to tracks and write
  // unassigned clusters
  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkHWClustersDataType);
       pDesc!=NULL; pDesc=GetNextInputBlock()) {
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
    if (fRawInputClusters) {
      fRawInputClusters->AddInputBlock(pDesc);
    }
    if (GetBenchmarkInstance()) {
      GetBenchmarkInstance()->Stop(1);
      GetBenchmarkInstance()->Start(5);
    }
    inputRawClusterSize+=pDesc->fSize;
    iResult=fRawInputClusters->Write(outputPtr+size, capacity-size, outputBlocks, fpDataDeflater);
    if (iResult>=0) {
      size+=iResult;
      if (GetBenchmarkInstance()) GetBenchmarkInstance()->AddOutput(iResult);
    }
    if (GetBenchmarkInstance()) {
      GetBenchmarkInstance()->Stop(5);
    }
    fRawInputClusters->Clear();
  }

  float compressionFactor=(float)inputRawClusterSize;
  if ((size)>0) compressionFactor/=size;
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

  // forward MC labels
  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo | kAliHLTDataOriginTPC);
       pDesc!=NULL; pDesc=GetNextInputBlock()) {
    outputBlocks.push_back(*pDesc);
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

  std::auto_ptr<AliHLTTPCHWCFSpacePointContainer> rawInputClusters(new AliHLTTPCHWCFSpacePointContainer(1));
  std::auto_ptr<AliHLTTPCSpacePointContainer> inputClusters(new AliHLTTPCSpacePointContainer);
  std::auto_ptr<TH1F> histoCompFactor(new TH1F("factor", "HLT TPC data compression factor", 100, 0, 10));

  if (!rawInputClusters.get() || !inputClusters.get() || !histoCompFactor.get()) return -ENOMEM;

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

    const unsigned nofParameters=AliHLTTPCDefinitions::GetNumberOfClusterParameterDefinitions();
    for (unsigned p=0; p<nofParameters; p++) {
      const AliHLTTPCDefinitions::AliClusterParameter& parameter=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[p];
      if (deflater->AddParameterDefinition(parameter.fName,
					   parameter.fBitLength,
					   parameter.fOptional)!=(int)parameter.fId) {
	// for performance reason the parameter id is simply used as index in the array of
	// definitions, the position must match the id
	HLTFatal("mismatch between parameter id and position in array, rearrange definitions!");
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
  } while (0);

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EPROTO;
  }

  return iResult;
}
