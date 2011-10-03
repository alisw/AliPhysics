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

    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=ReadRemainingClustersCompressed(fMonitoringContainer->BeginRemainingClusterBlock(0, pDesc->fSpecification),
					      reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr),
					      pDesc->fSize,
					      pDesc->fSpecification);
    }

    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::ClusterTracksCompressedDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=ReadTrackModelClustersCompressed(fMonitoringContainer->BeginTrackModelClusterBlock(0),
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

  if (fMonitoringContainer) {
    fMonitoringContainer->Clear();
    delete fMonitoringContainer;
  }
  fMonitoringContainer=NULL;

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

int AliHLTTPCDataCompressionMonitorComponent::ReadRemainingClustersCompressed(AliHLTTPCDataCompressionMonitorComponent::T& c, const AliHLTUInt8_t* pData, int dataSize, AliHLTUInt32_t specification)
{
  // read cluster data from AliHLTTPCClusterData
  int iResult=0;
  if (!pData  || dataSize<4) return -EINVAL;

  const AliHLTUInt8_t* pBuffer=pData;
  AliHLTUInt32_t size=dataSize;
  const AliHLTTPCRawClusterData* clusterData = reinterpret_cast<const AliHLTTPCRawClusterData*>(pBuffer);
  Int_t nCount = (Int_t) clusterData->fCount;

  AliHLTDataInflater* inflater=CreateInflater(clusterData->fVersion, 1);
  if (!inflater) return -ENODEV;

  if ((iResult=inflater->InitBitDataInput(reinterpret_cast<const AliHLTUInt8_t*>(clusterData->fClusters),
					  size-sizeof(AliHLTTPCRawClusterData)))<0) {
    return iResult;
  }

  iResult=ReadRemainingClustersCompressed(c, inflater, nCount, specification);

  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::ReadRemainingClustersCompressed(AliHLTTPCDataCompressionMonitorComponent::T& c, AliHLTDataInflater* pInflater, int nofClusters, AliHLTUInt32_t specification)
{
  // read cluster data

  int iResult=0;
  if (!pInflater) return -EINVAL;

  AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(specification);
  // the compressed format stores the difference of the local row number in
  // the partition to the row of the last cluster
  // add the first row in the partition to get global row number
  // offline uses row number in physical sector, inner sector consists of
  // partitions 0 and 1, outer sector of partition 2-5
  int rowOffset=AliHLTTPCTransform::GetFirstRow(partition);//-(partition<2?0:AliHLTTPCTransform::GetFirstRow(2));

  int parameterId=0;
  int outClusterCnt=0;
  AliHLTUInt64_t value=0;
  AliHLTUInt32_t length=0;
  AliHLTUInt32_t lastPadRow=0;
  while (outClusterCnt<nofClusters && pInflater->NextValue(value, length)) {
    const AliHLTTPCDefinitions::AliClusterParameter& parameter
      =AliHLTTPCDefinitions::fgkClusterParameterDefinitions[parameterId];

    if (parameter.fBitLength!=(int)length) {
      HLTError("decode error: expecting length %d for parameter %s, but got %d",
	       parameter.fBitLength, parameter.fName, length);
      break;
    }

    switch (parameterId) {
    case AliHLTTPCDefinitions::kPadRow:
      {c.SetPadRow(value+lastPadRow+rowOffset); lastPadRow+=value;break;}
    case AliHLTTPCDefinitions::kPad:
      {float pad=value; pad/=parameter.fScale; c.SetPad(pad); break;}
    case AliHLTTPCDefinitions::kTime:
      {float time=value; time/=parameter.fScale; c.SetTime(time); break;}
    case AliHLTTPCDefinitions::kSigmaY2:
      {float sigmaY2=value; sigmaY2/=parameter.fScale; c.SetSigmaY2(sigmaY2); break;}
    case AliHLTTPCDefinitions::kSigmaZ2:
      {float sigmaZ2=value; sigmaZ2/=parameter.fScale; c.SetSigmaZ2(sigmaZ2); break;}
    case AliHLTTPCDefinitions::kCharge:
      {c.SetCharge(value); break;}
    case AliHLTTPCDefinitions::kQMax:
      {c.SetQMax(value); break;}
    }
    if (parameterId>=AliHLTTPCDefinitions::kLast) {
      // switch to next cluster
      ++c;
      outClusterCnt++;
      parameterId=-1;
    }
    parameterId++;
  }
  pInflater->Pad8Bits();
  AliHLTUInt8_t bit=0;
  if (pInflater->InputBit(bit)) {
    HLTWarning("format error of compressed clusters, there is more data than expected");
  }
  pInflater->CloseBitDataInput();
  if (iResult>=0 && nofClusters!=outClusterCnt) {
    // is this a Fatal?
    HLTError("error reading compressed cluster format of block 0x%08x: expected %d, read only %d cluster(s)", specification, nofClusters, outClusterCnt);
    return -EPROTO;
  }
  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::ReadTrackModelClustersCompressed(AliHLTTPCDataCompressionMonitorComponent::T& c, const AliHLTUInt8_t* pData, int dataSize, AliHLTUInt32_t /*specification*/)
{
  // read cluster data from the track model data block
  int iResult=0;
  int dataOffset=sizeof(AliHLTTPCDataCompressionComponent::AliHLTTPCTrackModelBlock);
  if (!pData  || dataSize<dataOffset) return -EINVAL;

  const AliHLTTPCDataCompressionComponent::AliHLTTPCTrackModelBlock* trackModelBlock=reinterpret_cast<const AliHLTTPCDataCompressionComponent::AliHLTTPCTrackModelBlock*>(pData);
  if (trackModelBlock->fVersion!=1) {
    HLTError("unknown version %d", trackModelBlock->fVersion);
    return -EINVAL;
  }
  std::auto_ptr<AliHLTDataInflater> pInflater(CreateInflater(trackModelBlock->fDeflaterMode, 2));
  if (!pInflater.get()) {
    HLTError("failed to create the data inflater for mode %d", trackModelBlock->fDeflaterMode);
  }
  int nofTracks=trackModelBlock->fTrackCount;
  dataOffset+=trackModelBlock->fGlobalParameterCnt*sizeof(trackModelBlock->fGlobalParameters);
  if (dataSize<dataOffset) {
    HLTError("inconsistent data block, size %d, expecting at least %d to read AliHLTTPCTrackModelBlock with %d global parameters", dataSize, dataOffset, trackModelBlock->fGlobalParameterCnt);
    return -ENOSPC;
  }
  float bz=0.0;
  float driftTimeFactorA=0.;
  float driftTimeOffsetA=0.;
  float driftTimeFactorC=0.;
  float driftTimeOffsetC=0.;

  AliHLTUInt32_t parameterIndex=0;
  switch (trackModelBlock->fGlobalParameterCnt) {
  case 5:
    bz              =trackModelBlock->fGlobalParameters[parameterIndex++];
    driftTimeFactorA=trackModelBlock->fGlobalParameters[parameterIndex++];
    driftTimeOffsetA=trackModelBlock->fGlobalParameters[parameterIndex++];
    driftTimeFactorC=trackModelBlock->fGlobalParameters[parameterIndex++];
    driftTimeOffsetC=trackModelBlock->fGlobalParameters[parameterIndex++];
    break;
  default:
    HLTError("unknown version of global parameters %d", trackModelBlock->fGlobalParameterCnt);
    return -ENODATA;
  }

  if (parameterIndex!=trackModelBlock->fGlobalParameterCnt) {
    HLTError("internal error, size of parameter array has changed without providing all values");
    return -EFAULT;
  }

  for (int trackno=0; trackno<nofTracks; trackno++) {
    AliHLTTPCTrackGeometry trackpoints;
    trackpoints.InitDriftTimeTransformation(driftTimeFactorA, driftTimeOffsetA, driftTimeFactorC, driftTimeOffsetC);
    AliHLTUInt32_t  clusterBlockSize=0;
    if ((iResult=trackpoints.Read(pData+dataOffset, dataSize-dataOffset, bz, clusterBlockSize))<0) {
      return iResult;
    }
    dataOffset+=iResult;
    if (dataSize-dataOffset<(int)clusterBlockSize) {
      HLTError("to little data in buffer to read cluster block of size %d for track no %d", clusterBlockSize, trackno);
      return -ENODATA;
    }
    if ((iResult=pInflater->InitBitDataInput(pData+dataOffset, clusterBlockSize))<0) {
      return iResult;
    }
    if ((iResult=ReadTrackClustersCompressed(c, pInflater.get(), &trackpoints))<0) {
      HLTError("reading of associated clusters failed for track %d", trackno);
      return iResult;
    }
    pInflater->Pad8Bits();
    AliHLTUInt8_t bit=0;
    if (pInflater->InputBit(bit)) {
      HLTWarning("format error of compressed clusters, there is more data than expected");
    }
    pInflater->CloseBitDataInput();
    dataOffset+=clusterBlockSize;
  }

  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::ReadTrackClustersCompressed(AliHLTTPCDataCompressionMonitorComponent::T& c, AliHLTDataInflater* pInflater, AliHLTTPCTrackGeometry* pTrackPoints)
{
  // read cluster data

  int iResult=0;
  if (!pInflater || !pTrackPoints) return -EINVAL;

  const vector<AliHLTTrackGeometry::AliHLTTrackPoint>& rawTrackPoints=pTrackPoints->GetRawPoints();
  vector<AliHLTTrackGeometry::AliHLTTrackPoint>::const_iterator currentTrackPoint=rawTrackPoints.begin();

  bool bReadSuccess=true;
  AliHLTUInt32_t clusterCountBitLength=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kClusterCount].fBitLength;
  //unsigned long dataPosition=pInflater->GetCurrentByteInputPosition();
  for (unsigned row=0; row<159 && bReadSuccess; row++) {
    AliHLTUInt8_t haveClusters=0;
    // 1 bit for clusters on that padrow
    bReadSuccess=bReadSuccess && pInflater->InputBit(haveClusters);
    if (!haveClusters) continue;
    bool bEscape=false;
    do {
      if (currentTrackPoint==rawTrackPoints.end()) {
	if (bEscape || rawTrackPoints.begin()==rawTrackPoints.end()) break;
	currentTrackPoint=rawTrackPoints.begin();
	bEscape=true;
      }
      if (AliHLTTPCTransform::GetFirstRow(AliHLTTPCSpacePointData::GetPatch(currentTrackPoint->GetId())) +
	  AliHLTTPCSpacePointData::GetNumber(currentTrackPoint->GetId()) == row) {
	break;
      }
      currentTrackPoint++;
    } while (!bEscape);
    if (currentTrackPoint==rawTrackPoints.end()) {
      HLTError("decoding error, can not find track point on row %d", row);
      return -EFAULT;
    }
    AliHLTUInt8_t nofClusters=0;
    bReadSuccess=bReadSuccess && pInflater->InputBits(nofClusters, clusterCountBitLength);
    if (!bReadSuccess) break;

    static const AliHLTTPCDefinitions::AliClusterParameterId_t kParameterIdMapping[] = {
      AliHLTTPCDefinitions::kResidualPad,
      AliHLTTPCDefinitions::kResidualTime,
      AliHLTTPCDefinitions::kSigmaY2,
      AliHLTTPCDefinitions::kSigmaZ2,
      AliHLTTPCDefinitions::kCharge,
      AliHLTTPCDefinitions::kQMax,
    };

    int parameterId=0;
    int inClusterCnt=0;
    AliHLTUInt64_t value=0;
    AliHLTUInt32_t length=0;
    while (bReadSuccess && inClusterCnt<nofClusters && pInflater->NextValue(value, length)) {
      const AliHLTTPCDefinitions::AliClusterParameter& parameter
	=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[kParameterIdMapping[parameterId]];

      if (parameter.fBitLength!=(int)length) {
	HLTError("decode error: expecting length %d for parameter %s, but got %d",
		 parameter.fBitLength, parameter.fName, length);
	break;
      }

      static float deltapad=0.;
      static float deltatime=0.;
      bool lastParameter=false;
      switch (kParameterIdMapping[parameterId]) {
      case AliHLTTPCDefinitions::kResidualPad:
	{
	  AliHLTUInt8_t sign=0;
	  bReadSuccess=bReadSuccess && pInflater->InputBit(sign);
	  float pad=value*(sign?-1.:1.); pad/=parameter.fScale;
	  deltapad=pad;
	  pad+=currentTrackPoint->GetU();
	  c.SetPad(pad); 
	  break;
	}
      case AliHLTTPCDefinitions::kResidualTime:
	{
	  AliHLTUInt8_t sign=0;
	  bReadSuccess=bReadSuccess && pInflater->InputBit(sign);
	  float time=value*(sign?-1.:1.); time/=parameter.fScale;
	  deltatime=time;
	  time+=currentTrackPoint->GetV();
	  c.SetTime(time);
	  break;
	}
      case AliHLTTPCDefinitions::kSigmaY2:
	{float sigmaY2=value; sigmaY2/=parameter.fScale; c.SetSigmaY2(sigmaY2); break;}
      case AliHLTTPCDefinitions::kSigmaZ2:
	{float sigmaZ2=value; sigmaZ2/=parameter.fScale; c.SetSigmaZ2(sigmaZ2); break;}
      case AliHLTTPCDefinitions::kCharge:
	{c.SetCharge(value); break;}
      case AliHLTTPCDefinitions::kQMax:
	{c.SetQMax(value); lastParameter=true; break;}
      default:
	{
	  HLTError("parameter %d not expected", kParameterIdMapping[parameterId]);
	}
      }
      if (lastParameter) {
	// switch to next cluster
	c.SetPadRow(row);
	// cout << "  row "    << setfill(' ') << setw(3) << fixed << right                     << c.GetRow()
	//      << "  pad "    << setfill(' ') << setw(7) << fixed << right << setprecision (4) << c.GetPad()
	//      << "  dpad "   << setfill(' ') << setw(7) << fixed << right << setprecision (4) << deltapad
	//      << "  time "   << setfill(' ') << setw(7) << fixed << right << setprecision (4) << c.GetTimeBin()
	//      << "  dtime "  << setfill(' ') << setw(7) << fixed << right << setprecision (4) << deltatime
	//      << "  charge " << setfill(' ') << setw(5) << fixed << right << setprecision (0) << c.GetQ()
	//      << "  qmax "   << setfill(' ') << setw(4) << fixed << right << setprecision (0) << c.GetMax()
	//      << endl;
	++c;
	inClusterCnt++;
	parameterId=-1;
      }
      parameterId++;
    }
    if (iResult>=0 && nofClusters!=inClusterCnt) {
      // is this a Fatal?
      HLTError("error reading track model compressed cluster format of track: expected %d, read only %d cluster(s)", nofClusters, inClusterCnt);
      return -EPROTO;
    }
    currentTrackPoint++;
  }
  return iResult;
}

AliHLTDataInflater* AliHLTTPCDataCompressionMonitorComponent::CreateInflater(int deflater, int mode) const
{
  // create the inflater for the specified mode
  vector<AliHLTTPCDefinitions::AliClusterParameterId_t> parameterids;
  switch (mode) {
  case 1:
    parameterids.push_back(AliHLTTPCDefinitions::kPadRow );
    parameterids.push_back(AliHLTTPCDefinitions::kPad    );
    parameterids.push_back(AliHLTTPCDefinitions::kTime   );
    parameterids.push_back(AliHLTTPCDefinitions::kSigmaY2);
    parameterids.push_back(AliHLTTPCDefinitions::kSigmaZ2);
    parameterids.push_back(AliHLTTPCDefinitions::kCharge );
    parameterids.push_back(AliHLTTPCDefinitions::kQMax   );
    break;
  case 2:
    parameterids.push_back(AliHLTTPCDefinitions::kResidualPad );
    parameterids.push_back(AliHLTTPCDefinitions::kResidualTime);
    parameterids.push_back(AliHLTTPCDefinitions::kSigmaY2);
    parameterids.push_back(AliHLTTPCDefinitions::kSigmaZ2);
    parameterids.push_back(AliHLTTPCDefinitions::kCharge );
    parameterids.push_back(AliHLTTPCDefinitions::kQMax   );
    break;
  default:
    HLTError("invalid mode %d for inflater initialization", mode);
  }

  switch (deflater) {
  case 1:
    {
      std::auto_ptr<AliHLTDataInflaterSimple> inflatersimple(new AliHLTDataInflaterSimple);
      if (!inflatersimple.get()) return NULL;
      for (vector<AliHLTTPCDefinitions::AliClusterParameterId_t>::const_iterator id=parameterids.begin();
	   id!=parameterids.end(); id++) {
	const AliHLTTPCDefinitions::AliClusterParameter& parameter=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[*id];
	if (inflatersimple->AddParameterDefinition(parameter.fName,
						   parameter.fBitLength,
						   parameter.fOptional)<0) {
	  HLTError("error adding parameter definition %s to inflater", parameter.fName);
	  return NULL;
	}
      }
      return inflatersimple.release();
    }
    break;
  case 2:
    {
      std::auto_ptr<AliHLTDataInflaterHuffman> inflaterhuffman(new AliHLTDataInflaterHuffman);
      if (!inflaterhuffman.get()) return NULL;
      TString cdbPath("HLT/ConfigTPC/TPCDataCompressorHuffmanTables");
      TObject* pConf=LoadAndExtractOCDBObject(cdbPath);
      if (!pConf) {
	HLTError("can not load configuration object %s", cdbPath.Data());
	return NULL;
      }
      if (dynamic_cast<TList*>(pConf)==NULL) {
	HLTError("huffman table configuration object of inconsistent type");
	return NULL;
      }
      inflaterhuffman->InitDecoders(dynamic_cast<TList*>(pConf));
      for (vector<AliHLTTPCDefinitions::AliClusterParameterId_t>::const_iterator id=parameterids.begin();
	   id!=parameterids.end(); id++) {
	const AliHLTTPCDefinitions::AliClusterParameter& parameter=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[*id];
	if (inflaterhuffman->AddParameterDefinition(parameter.fName,
						    parameter.fBitLength)<0) {
	  HLTError("error adding parameter definition %s to inflater", parameter.fName);
	  return NULL;
	}
      }
      return inflaterhuffman.release();
    }
    break;
  default:
    HLTError("unknown inflater requested %d", deflater);
  }
  return NULL;
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
  fEnd=iterator();
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
  fEnd=iterator();
  return fBegin;
}

const AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::iterator& AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::End()
{
  /// get end iterator
  return fEnd;
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
