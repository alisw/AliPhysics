//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCDATACOMPRESSIONDECODER_H
#define ALIHLTTPCDATACOMPRESSIONDECODER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCDataCompressionDecoder.h
/// @author Matthias Richter
/// @date   2011-10-04
/// @brief  Generic decoder class for compressed TPC data, works on a container
///         class implementation which fills the actual target data struct

#include "AliHLTLogging.h"
#include "AliHLTMisc.h"
#include "AliHLTTPCDataCompressionComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackGeometry.h"
#include "AliHLTDataInflater.h"

/**
 * @class AliHLTTPCDataCompressionDecoder
 * Generic decoder class for compressed TPC data, works on a container
 * class implementation which fills the actual target data struct
 */
class AliHLTTPCDataCompressionDecoder : public AliHLTLogging {
 public:
  AliHLTTPCDataCompressionDecoder();
  ~AliHLTTPCDataCompressionDecoder();

  template<typename T>
  int ReadRemainingClustersCompressed(T& c, const AliHLTUInt8_t* pData, int dataSize, AliHLTUInt32_t specification);

  template<typename T>
  int ReadRemainingClustersCompressed(T& c, AliHLTDataInflater* pInflater, int nofClusters, AliHLTUInt32_t specification);

  template<typename T>
  int ReadTrackModelClustersCompressed(T& c, const AliHLTUInt8_t* pData, int dataSize, AliHLTUInt32_t specification);

  template<typename T>
  int ReadTrackClustersCompressed(T& c, AliHLTDataInflater* pInflater, AliHLTTPCTrackGeometry* pTrackPoints);

  AliHLTDataInflater* CreateInflater(int deflater, int mode) const;

 protected:
 private:
  ClassDef(AliHLTTPCDataCompressionDecoder, 0)
};

template<typename T>
int AliHLTTPCDataCompressionDecoder::ReadRemainingClustersCompressed(T& c, const AliHLTUInt8_t* pData, int dataSize, AliHLTUInt32_t specification)
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

template<typename T>
int AliHLTTPCDataCompressionDecoder::ReadRemainingClustersCompressed(T& c, AliHLTDataInflater* pInflater, int nofClusters, AliHLTUInt32_t specification)
{
  // read cluster data

  int iResult=0;
  if (!pInflater) return -EINVAL;

  AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(specification);
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
      c.Next(slice, partition);
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

template<typename T>
int AliHLTTPCDataCompressionDecoder::ReadTrackModelClustersCompressed(T& c, const AliHLTUInt8_t* pData, int dataSize, AliHLTUInt32_t /*specification*/)
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

template<typename T>
int AliHLTTPCDataCompressionDecoder::ReadTrackClustersCompressed(T& c, AliHLTDataInflater* pInflater, AliHLTTPCTrackGeometry* pTrackPoints)
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
    AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(currentTrackPoint->GetId());
    AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(currentTrackPoint->GetId());
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
	c.Next(slice, partition);
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
#endif
