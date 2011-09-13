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
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCHWCFData.h"
#include "AliRawDataHeader.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TFile.h"
#include <memory>

ClassImp(AliHLTTPCDataCompressionMonitorComponent)

AliHLTTPCDataCompressionMonitorComponent::AliHLTTPCDataCompressionMonitorComponent()
  : AliHLTProcessor()
  , fpHWClusterDecoder(NULL)
  , fHistoHWCFDataSize(NULL)
  , fHistoHWCFReductionFactor(NULL)
  , fHistoNofClusters(NULL)
  , fHistogramFile("HLT.TPC-compression-statistics.root")
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
      out.Close();
    }
  }
  if (fHistoHWCFDataSize) delete fHistoHWCFDataSize;
  fHistoHWCFDataSize=NULL;
  if (fHistoHWCFReductionFactor) delete fHistoHWCFReductionFactor;
  fHistoHWCFReductionFactor=NULL;
  if (fHistoNofClusters) delete fHistoNofClusters;
  fHistoNofClusters=NULL;

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
