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

/// @file   AliHLTTPCDataCompressionMonitorComponent.cxx
/// @author Matthias Richter
/// @date   2011-09-12
/// @brief  TPC component for monitoring of data compression
///

#include "AliHLTTPCDataCompressionMonitorComponent.h"
#include "AliHLTTPCDataCompressionComponent.h"
#include "AliHLTTPCDataCompressionDecoder.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCHWCFData.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackGeometry.h"
#include "AliHLTTPCHWCFSpacePointContainer.h"
#include "AliHLTDataInflaterSimple.h"
#include "AliHLTDataInflaterHuffman.h"
#include "AliRawDataHeader.h"
#include "AliTPCclusterMI.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TFile.h"
#include "TObjArray.h"
#include <memory>

ClassImp(AliHLTTPCDataCompressionMonitorComponent)

AliHLTTPCDataCompressionMonitorComponent::AliHLTTPCDataCompressionMonitorComponent()
  : AliHLTProcessor()
  , fpHWClusterDecoder(NULL)
  , fHistoHWCFDataSize(NULL)
  , fHistoHWCFReductionFactor(NULL)
  , fHistoNofClusters(NULL)
  , fHistogramFile("HLT.TPC-compression-statistics.root")
  , fMonitoringContainer(NULL)
  , fVerbosity(0)
  , fFlags(0)
{
}

AliHLTTPCDataCompressionMonitorComponent::~AliHLTTPCDataCompressionMonitorComponent()
{
  /// destructor
}


const char* AliHLTTPCDataCompressionMonitorComponent::GetComponentID()
{
  /// inherited from AliHLTComponent: id of the component
  return "TPCDataCompressorMonitor";
}


void AliHLTTPCDataCompressionMonitorComponent::GetInputDataTypes( AliHLTComponentDataTypeList& tgtList)
{
  /// inherited from AliHLTComponent: list of data types in the vector reference
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkHWClustersDataType);
  tgtList.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
  tgtList.push_back(AliHLTTPCDefinitions::fgkRawClustersDataType);
  tgtList.push_back(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
  tgtList.push_back(AliHLTTPCDefinitions::ClusterTracksCompressedDataType());  
}

AliHLTComponentDataType AliHLTTPCDataCompressionMonitorComponent::GetOutputDataType()
{
  /// inherited from AliHLTComponent: output data type of the component.
  return kAliHLTMultipleDataType;
}

int AliHLTTPCDataCompressionMonitorComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  /// inherited from AliHLTComponent: multiple output data types of the component.
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
  return tgtList.size();
}

void AliHLTTPCDataCompressionMonitorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  /// inherited from AliHLTComponent: output data size estimator
  constBase=0;
  inputMultiplier=1.0;
}

AliHLTComponent* AliHLTTPCDataCompressionMonitorComponent::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTTPCDataCompressionMonitorComponent;
}

int AliHLTTPCDataCompressionMonitorComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, 
						       const AliHLTComponentBlockData* /*inputBlocks*/, 
						       AliHLTComponentTriggerData& /*trigData*/,
						       AliHLTUInt8_t* /*outputPtr*/,
						       AliHLTUInt32_t& /*size*/,
						       AliHLTComponentBlockDataList& /*outputBlocks*/ )
{
  /// inherited from AliHLTProcessor: data processing
  int iResult=0;

  if (!IsDataEvent()) return 0;

  const AliHLTComponentBlockData* pDesc=NULL;
  unsigned rawDataSize=0;
  unsigned hwclustersDataSize=0;
  unsigned nofClusters=0;
  for (pDesc=GetFirstInputBlock(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC);
       pDesc!=NULL; pDesc=GetNextInputBlock()) {
    fFlags|=kHaveRawData;
    rawDataSize+=pDesc->fSize;
  }

  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkHWClustersDataType);
       pDesc!=NULL; pDesc=GetNextInputBlock()) {
    fFlags|=kHaveHWClusters;
    // FIXME: the decoding can now be handled via the data container
    if (pDesc->fSize<=sizeof(AliRawDataHeader)) continue;
    if (fpHWClusterDecoder) {
      hwclustersDataSize+=pDesc->fSize;
      AliHLTUInt8_t* pData=reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr);
      pData+=sizeof(AliRawDataHeader);
      if (fpHWClusterDecoder->Init(pData, pDesc->fSize-sizeof(AliRawDataHeader))<0 ||
	  (fpHWClusterDecoder->CheckVersion()<0 && (int)pDesc->fSize>fpHWClusterDecoder->GetRCUTrailerSize())) {
	HLTError("data block of type %s corrupted: can not decode format",
		 AliHLTComponent::DataType2Text(pDesc->fDataType).c_str());
      } else {
	nofClusters+=fpHWClusterDecoder->GetNumberOfClusters();
      }
    }
    if (fMonitoringContainer) {
      fMonitoringContainer->AddRawData(pDesc);
    }
  }

  if (fMonitoringContainer) {
    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::RemainingClusterIdsDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=fMonitoringContainer->AddClusterIds(pDesc);
    }

    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::ClusterIdTracksDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=fMonitoringContainer->AddClusterIds(pDesc);
    }

    // read data
    AliHLTTPCDataCompressionDecoder decoder;
    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=decoder.ReadRemainingClustersCompressed(fMonitoringContainer->BeginRemainingClusterBlock(0, pDesc->fSpecification),
					      reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr),
					      pDesc->fSize,
					      pDesc->fSpecification);
    }

    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::ClusterTracksCompressedDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=decoder.ReadTrackModelClustersCompressed(fMonitoringContainer->BeginTrackModelClusterBlock(0),
					       reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr),
					       pDesc->fSize,
					       pDesc->fSpecification);
    }

    fMonitoringContainer->Clear();
  }

  float ratio=0;
  if (hwclustersDataSize) {ratio=(float)rawDataSize; ratio/=hwclustersDataSize;}
  if (fHistoHWCFDataSize)  fHistoHWCFDataSize->Fill(rawDataSize/1024, hwclustersDataSize/1024);
  if (fHistoHWCFReductionFactor)  fHistoHWCFReductionFactor->Fill(rawDataSize/1024, ratio);
  if (fHistoNofClusters) fHistoNofClusters->Fill(rawDataSize/1024, nofClusters);
  HLTInfo("raw data %d, hwcf data %d, ratio %f, %d clusters", rawDataSize, hwclustersDataSize, ratio, nofClusters);

  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::DoInit( int argc, const char** argv )
{
  /// inherited from AliHLTComponent: component initialisation and argument scan.
  int iResult=0;

  // component configuration
  //Stage 1: default initialization.
  //Default values.
  fFlags=0;

  //Stage 2: OCDB.
  TString cdbPath("HLT/ConfigTPC/");
  cdbPath += GetComponentID();
  //
  // iResult = ConfigureFromCDBTObjString(cdbPath);
  // if (iResult < 0) 
  //   return iResult;

  //Stage 3: command line arguments.
  if (argc && (iResult = ConfigureFromArgumentString(argc, argv)) < 0)
    return iResult;

  std::auto_ptr<AliHLTTPCHWCFData> hwClusterDecoder(new AliHLTTPCHWCFData);
  std::auto_ptr<AliDataContainer> dataContainer(new AliDataContainer);

  std::auto_ptr<TH2I> histoHWCFDataSize(new TH2I("HWCFDataSize",
						 "HW ClusterFinder Size",
						 100, 0., 80000., 100, 0., 80000.));
  if (histoHWCFDataSize.get()) {
    TAxis* xaxis=histoHWCFDataSize->GetXaxis();
    if (xaxis) xaxis->SetTitle("raw event size [kB]");
    TAxis* yaxis=histoHWCFDataSize->GetYaxis();
    if (yaxis) yaxis->SetTitle("hwcf size");
  }

  std::auto_ptr<TH2I> histoHWCFReductionFactor(new TH2I("HWCFReductionFactor",
							"Data reduction HW ClusterFinder",
							100, 0., 80000., 100, 0., 10.));
  if (histoHWCFReductionFactor.get()) {
    TAxis* xaxis=histoHWCFReductionFactor->GetXaxis();
    if (xaxis) xaxis->SetTitle("raw event size [kB]");
    TAxis* yaxis=histoHWCFReductionFactor->GetYaxis();
    if (yaxis) yaxis->SetTitle("reduction factor");
  }

  std::auto_ptr<TH2I> histoNofClusters(new TH2I("NofClusters",
					       "Number of HLT TPC clusters",
					       100, 0., 80000., 100, 0., 3000000.));
  if (histoNofClusters.get()) {
    TAxis* xaxis=histoNofClusters->GetXaxis();
    if (xaxis) xaxis->SetTitle("event size [kB]");
    TAxis* yaxis=histoNofClusters->GetYaxis();
    if (yaxis) yaxis->SetTitle("count");
  }

  // initialize the histograms if stored at the end
  // condition might be extended
  if (!fHistogramFile.IsNull()) {
    fHistoHWCFDataSize=histoHWCFDataSize.release();
    fHistoHWCFReductionFactor=histoHWCFReductionFactor.release();
    fHistoNofClusters=histoNofClusters.release();
  }

  fpHWClusterDecoder=hwClusterDecoder.release();
  fMonitoringContainer=dataContainer.release();

  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::DoDeinit()
{
  /// inherited from AliHLTComponent: component cleanup
  int iResult=0;

  if (fpHWClusterDecoder) delete fpHWClusterDecoder;
  fpHWClusterDecoder=NULL;

  if (!fHistogramFile.IsNull()) {
    TFile out(fHistogramFile, "RECREATE");
    if (!out.IsZombie()) {
      out.cd();
      if (fHistoHWCFDataSize) fHistoHWCFDataSize->Write();
      if (fHistoHWCFReductionFactor) fHistoHWCFReductionFactor->Write();
      if (fHistoNofClusters) fHistoNofClusters->Write();
      if (fMonitoringContainer) {
	const TObject* o=fMonitoringContainer->FindObject("histograms");
	if (o) o->Write();
      }
      out.Close();
    }
  }
  if (fHistoHWCFDataSize) delete fHistoHWCFDataSize;
  fHistoHWCFDataSize=NULL;
  if (fHistoHWCFReductionFactor) delete fHistoHWCFReductionFactor;
  fHistoHWCFReductionFactor=NULL;
  if (fHistoNofClusters) delete fHistoNofClusters;
  fHistoNofClusters=NULL;
  if (fMonitoringContainer) {
    fMonitoringContainer->Clear();
    delete fMonitoringContainer;
  }
  fMonitoringContainer=NULL;


  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  /// inherited from AliHLTComponent: argument scan
  int iResult=0;
  if (argc<1) return 0;
  int bMissingParam=0;
  int i=0;
  TString argument=argv[i];

  do {
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

AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::AliDataContainer()
  : fHistograms(new TObjArray)
  , fHistogramPointers()
  , fRemainingClusterIds()
  , fTrackModelClusterIds()
  , fCurrentClusterIds(NULL)
  , fRawData(NULL)
  , fBegin()
  , fEnd()
{
  /// constructor
  if (fHistograms) {
    fHistograms->SetOwner(kTRUE);
    fHistogramPointers.resize(kNumberOfHistograms, NULL);
    for (const AliHistogramDefinition* definition=fgkHistogramDefinitions;
	 definition->fName!=NULL; definition++) {
      fHistogramPointers[definition->fId]=new TH1F(definition->fName,
						  definition->fTitle,
						  definition->fBins,
						  definition->fLowerBound,
						  definition->fUpperBound
						  );
      fHistograms->AddAt(fHistogramPointers[definition->fId], definition->fId);
    }
  }
}

const AliHLTTPCDataCompressionMonitorComponent::AliHistogramDefinition AliHLTTPCDataCompressionMonitorComponent::fgkHistogramDefinitions[] = {
  {kHistogramPadrow,        "padrow"   , "padrow"   ,  159,   0.,   158.},
  {kHistogramPad,           "pad"      , "pad"      ,  140,   0.,   139.},
  {kHistogramTime,          "time"     , "time"     , 1024,   0.,  1023.},
  {kHistogramSigmaY2,       "sigmaY2"  , "sigmaY2"  ,  100,   0.,     1.},
  {kHistogramSigmaZ2,       "sigmaZ2"  , "sigmaZ2"  ,  100,   0.,     1.},
  {kHistogramCharge,        "chareg"   , "charge"   , 1024,   0., 65535.},
  {kHistogramQMax,          "qmax"     , "qmax"     ,  128,   0.,  1023.},
  {kHistogramDeltaPadrow,   "d_padrow" , "d_padrow" , 1000,  -1.,     1.},
  {kHistogramDeltaPad,      "d_pad"    , "d_pad"    , 1000,  -1.,     1.},
  {kHistogramDeltaTime,     "d_time"   , "d_time"   , 1000,  -1.,     1.},
  {kHistogramDeltaSigmaY2,  "d_sigmaY2", "d_sigmaY2", 1000,  -1.,     1.},
  {kHistogramDeltaSigmaZ2,  "d_sigmaZ2", "d_sigmaZ2", 1000,  -1.,     1.},
  {kHistogramDeltaCharge,   "d_chareg" , "d_charge" , 1000,  -1.,     1.},
  {kHistogramDeltaQMax,     "d_qmax"   , "d_qmax"   , 1000,  -1.,     1.},
  {kNumberOfHistograms, NULL    ,  NULL    ,    0,  0.,     0.}
};

AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::~AliDataContainer()
{
  /// dectructor
  if (fRawData) delete fRawData;
  if (fHistograms) delete fHistograms;
}

AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::iterator& AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::BeginRemainingClusterBlock(int /*count*/, AliHLTUInt32_t specification)
{
  /// iterator of remaining clusters block of specification
  AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(specification);
  AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(specification);
  unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
  if (index<fRemainingClusterIds.size())
    fCurrentClusterIds=&fRemainingClusterIds[index];
  else
    fCurrentClusterIds=NULL;
  fBegin=iterator(this);
  return fBegin;
}

AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::iterator& AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::BeginTrackModelClusterBlock(int /*count*/)
{
  /// iterator of track model clusters
  if (fTrackModelClusterIds.fIds && fTrackModelClusterIds.fSize>0)
    fCurrentClusterIds=&fTrackModelClusterIds;
  else
    fCurrentClusterIds=NULL;
  fBegin=iterator(this);
  return fBegin;
}

int AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::AddRawData(const AliHLTComponentBlockData* pDesc)
{
    /// add raw data bloack
  if (pDesc->fDataType==AliHLTTPCDefinitions::HWClustersDataType()) {
    if (!fRawData) fRawData=new AliHLTTPCHWCFSpacePointContainer(AliHLTTPCHWCFSpacePointContainer::kModeCreateMap);
    if (!fRawData) return -ENOMEM;
    return fRawData->AddInputBlock(pDesc);
  }
  return -ENODATA;  
}

int AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::AddClusterIds(const AliHLTComponentBlockData* pDesc)
{
  /// add cluster id block for remaining or track model clusters
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
    fTrackModelClusterIds.fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fTrackModelClusterIds.fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  if (pDesc->fDataType==AliHLTTPCDefinitions::RemainingClusterIdsDataType()) {
    AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(pDesc->fSpecification);
    AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(pDesc->fSpecification);
    unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
    if (fRemainingClusterIds.size()<=index) {
      if ((int)fRemainingClusterIds.size()<AliHLTTPCTransform::GetNSlice()*AliHLTTPCTransform::GetNumberOfPatches()) {
	fRemainingClusterIds.resize(AliHLTTPCTransform::GetNSlice()*AliHLTTPCTransform::GetNumberOfPatches());
      } else {
	fRemainingClusterIds.resize(index+1);
      }
    }
    fRemainingClusterIds[index].fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fRemainingClusterIds[index].fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  return -ENODATA;
}

AliHLTUInt32_t AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::GetClusterId(int clusterNo) const
{
  /// get the cluster id from the current cluster id block (optional)
  if (!fCurrentClusterIds ||
      clusterNo<0 ||
      (int)fCurrentClusterIds->fSize<=clusterNo)
    return kAliHLTVoidDataSpec;
  return fCurrentClusterIds->fIds[clusterNo];
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillPadRow(int row, AliHLTUInt32_t clusterId)
{
  /// fill padrow histogram
  unsigned index=kHistogramPadrow;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(row);
  if (clusterId!=kAliHLTVoidDataSpec) {
    index=kHistogramDeltaPadrow;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(row-fRawData->GetX(clusterId));
      }
    }
  }
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillPad(float pad, AliHLTUInt32_t clusterId)
{
  /// fill pad histogram
  unsigned index=kHistogramPad;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(pad);
  if (clusterId!=kAliHLTVoidDataSpec) {
    index=kHistogramDeltaPad;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(pad-fRawData->GetY(clusterId));
      }
    }
  }
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillTime(float time, AliHLTUInt32_t clusterId)
{
  /// fill pad histogram
  unsigned index=kHistogramTime;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(time);
  if (clusterId!=kAliHLTVoidDataSpec) {
    index=kHistogramDeltaTime;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(time-fRawData->GetZ(clusterId));
      }
    }
  }
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillSigmaY2(float sigmaY2, AliHLTUInt32_t clusterId)
{
  /// fill sigmaY2 histogram
  unsigned index=kHistogramSigmaY2;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(sigmaY2);
  if (clusterId!=kAliHLTVoidDataSpec) {
    index=kHistogramDeltaSigmaY2;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(sigmaY2-fRawData->GetYWidth(clusterId));
      }
    }
  }
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillSigmaZ2(float sigmaZ2, AliHLTUInt32_t clusterId)
{
  /// fill sigmaZ2 histogram
  unsigned index=kHistogramSigmaZ2;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(sigmaZ2);
  if (clusterId!=kAliHLTVoidDataSpec) {
    index=kHistogramDeltaSigmaZ2;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(sigmaZ2-fRawData->GetZWidth(clusterId));
      }
    }
  }
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillCharge(unsigned charge, AliHLTUInt32_t clusterId)
{
  /// fill charge histogram
  unsigned index=kHistogramCharge;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(charge);
  if (clusterId!=kAliHLTVoidDataSpec) {
    index=kHistogramDeltaCharge;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(charge-fRawData->GetCharge(clusterId));
      }
    }
  }
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillQMax(unsigned qmax, AliHLTUInt32_t clusterId)
{
  /// fill qmax histogram
  unsigned index=kHistogramQMax;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(qmax);
  if (clusterId!=kAliHLTVoidDataSpec) {
    index=kHistogramDeltaQMax;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(qmax-fRawData->GetQMax(clusterId));
      }
    }
  }
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::Clear(Option_t * option)
{
  /// internal cleanup
  if (fRawData) fRawData->Clear(option);
}

TObject* AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FindObject(const char */*name*/) const
{
  /// get histogram object  
  return fHistograms;
}
