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
#include "TString.h"

AliHLTTPCDataCompressionComponent::AliHLTTPCDataCompressionComponent()
  : AliHLTProcessor()
  , fMode(0)
  , fRawInputClusters(NULL)
  , fInputClusters(NULL)
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
  inputMultiplier=1.;
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
  AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
  const AliHLTComponentBlockData* pDesc=NULL;

  /// input track array
  vector<AliHLTGlobalBarrelTrack> inputTrackArray;

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->Start(1);
  }
  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkHWClustersDataType);
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
    if (fRawInputClusters) {
      fRawInputClusters->AddInputBlock(pDesc);
    }
  }
  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->Stop(1);
    GetBenchmarkInstance()->Start(2);
  }

  // transformed clusters
  if (fMode==1) { // FIXME: condition to be adjusted
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
  if (fMode==1) { // FIXME: condition to be adjusted
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

  // output
  if (fMode==0) {
    iResult=fRawInputClusters->Write(outputPtr, capacity-size, outputBlocks);
    if (iResult>=0) {
      size+=iResult;
      if (GetBenchmarkInstance()) GetBenchmarkInstance()->AddOutput(iResult);
    }
  }

  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->Stop(5);
    GetBenchmarkInstance()->Stop(0);
    HLTBenchmark(GetBenchmarkInstance()->GetStatistics());
  }

  if (fInputClusters) {
    fInputClusters->Clear();
  }
  if (fRawInputClusters) {
    fRawInputClusters->Clear();
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
  //iResult = ConfigureFromCDBTObjString(cdbPath);
  if (iResult < 0) 
    return iResult;

  //Stage 3: command line arguments.
  if (argc && (iResult = ConfigureFromArgumentString(argc, argv)) < 0)
    return iResult;

  fpBenchmark=new AliHLTComponentBenchmark;
  if (GetBenchmarkInstance()) {
    GetBenchmarkInstance()->SetTimer(0,"total");
    GetBenchmarkInstance()->SetTimer(1,"rawclusterinput");
    GetBenchmarkInstance()->SetTimer(2,"clusterinput");
    GetBenchmarkInstance()->SetTimer(3,"trackinput");
    GetBenchmarkInstance()->SetTimer(4,"processing");
    GetBenchmarkInstance()->SetTimer(5,"output");
  }

  fRawInputClusters=new AliHLTTPCHWCFSpacePointContainer;	  
  if (!fRawInputClusters) return -ENOMEM;

  fInputClusters=new AliHLTTPCSpacePointContainer;	  
  if (!fInputClusters) return -ENOMEM;

  return iResult;
}

int AliHLTTPCDataCompressionComponent::DoDeinit()
{
  /// inherited from AliHLTComponent: component cleanup
  int iResult=0;
  return iResult;
}

int AliHLTTPCDataCompressionComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  /// inherited from AliHLTComponent: argument scan
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  int i=0;
  for (; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -mode
    if (argument.CompareTo("-mode")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter=argv[i];
      if (parameter.IsDigit()) {
	fMode=parameter.Atoi();
      } else {
	HLTError("invalid parameter for argument %s, expecting number instead of %s", argument.Data(), parameter.Data());
	return -EPROTO;
      }
    }
  }

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EPROTO;
  }

  if (iResult>=0) return i;
  return iResult;
}
