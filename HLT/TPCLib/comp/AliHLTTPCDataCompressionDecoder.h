//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCDATACOMPRESSIONDECODER_H
#define ALIHLTTPCDATACOMPRESSIONDECODER_H
//* This file is property of and copyright by the                          * 
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
#include "AliHLTTPCGeometry.h"
#include "AliHLTTPCTrackGeometry.h"
#include "AliHLTDataInflater.h"
#include "AliHLTTPCHWClusterMerger.h"
#include "AliHLTTPCClusterFlagsData.h"
#include <vector>
#include <string>

/**
 * @class AliHLTTPCDataCompressionDecoder
 * Generic decoder class for compressed TPC data, works on a container
 * class implementation which fills the actual target data struct
 */
class AliHLTTPCDataCompressionDecoder : public AliHLTLogging {
 public:
  AliHLTTPCDataCompressionDecoder();
  ~AliHLTTPCDataCompressionDecoder();

  // extract the padrow number within partition
  void SetPadrowModeLocal() {fExtractGlobalPadrow=kFALSE;}
  // extract global padrow number including offset of first padrow in partition
  void SetPadrowModeGlobal() {fExtractGlobalPadrow=kTRUE;}

  template<typename T>
  int ReadRemainingClustersCompressed(T& c, const AliHLTUInt8_t* pData, int dataSize, AliHLTUInt32_t specification);

  template<typename T>
  int ReadRemainingClustersCompressed(T& c, AliHLTDataInflater* pInflater, int nofClusters, AliHLTUInt32_t specification, int formatVersion=0);

  template<typename T>
  int ReadTrackModelClustersCompressed(T& c, const AliHLTUInt8_t* pData, int dataSize, AliHLTUInt32_t specification);

  template<typename T>
  int ReadTrackClustersCompressed(T& c, AliHLTDataInflater* pInflater, AliHLTTPCTrackGeometry* pTrackPoints);

  template<typename T>
  int ReadClustersPartition(T& c, const AliHLTUInt8_t* pData, unsigned dataSize, AliHLTUInt32_t specification);

  AliHLTDataInflater* CreateInflater(int deflater, int mode) const;

  void SetPadShift(float padShift) {fPadShift=padShift;}
  float PadShift() const {return fPadShift;}
  void SetVerbosity(int verbosity) {fVerbosity=verbosity;}
  void SetHuffmanTableConfiguration(const char* name) {fHuffmanTableFile = name;}

  int InitPartitionClusterDecoding(AliHLTUInt32_t specification);
  int InitTrackModelClusterClusterDecoding();
  int AddCompressionDescriptor(const AliHLTComponentBlockData* pDesc);
  int AddRawClustersDescriptor(const AliHLTComponentBlockData* pDesc);
  int AddClusterMCData(const AliHLTComponentBlockData* pDesc);
  int AddClusterIds(const AliHLTComponentBlockData* pDesc);
  int AddClusterFlags(const AliHLTComponentBlockData* pDesc);
  AliHLTUInt32_t GetClusterId(int clusterNo) const;
  unsigned short GetNextClusterFlag();
  const AliHLTTPCClusterMCLabel* GetMCLabel(AliHLTUInt32_t clusterId) const;

  void Clear(const char* option="");

  struct AliClusterIdBlock {
    AliClusterIdBlock() : fIds(NULL), fSize(0) {}
    AliHLTUInt32_t* fIds; //!
    AliHLTUInt32_t  fSize; //!

    void Clear() {fIds=NULL; fSize=0;}
  };

 protected:
 private:
  AliHLTTPCDataCompressionDecoder(const AliHLTTPCDataCompressionDecoder&);
  AliHLTTPCDataCompressionDecoder& operator=(const AliHLTTPCDataCompressionDecoder&);

  float fPadShift; //! pad shift
  int fVerbosity; //! verbosity level
  Bool_t fUseClusterMerger; // flag to run the cluster merger
  Bool_t fExtractGlobalPadrow; //! extract global padrow including offset of first padrow in a partition
  AliHLTDataInflater* fpDataInflaterPartition; //! instance of inflater for partition clusters
  AliHLTDataInflater* fpDataInflaterTrack; //! instance of inflater for track clusters
  AliHLTTPCHWClusterMerger* fpClusterMerger; //! merger instance

  vector<AliClusterIdBlock> fPartitionClusterIds; //! clusters ids for clusters of individual partitions
  vector<AliHLTTPCClusterFlagsData*> fPartitionClusterFlags; //! same for cluster flags
  AliClusterIdBlock fTrackModelClusterIds; //! cluster ids for track model clusters
  AliClusterIdBlock* fCurrentClusterIds; //! id block currently active in the iteration
  AliHLTTPCClusterFlagsData* fCurrentClusterFlags; //! same for flags
  unsigned int fDecodeFlagsTmpFlag; //!
  unsigned int fDecodeFlagsTmpPos; //!
  unsigned int fDecodeFlagsTmpEntries; //!
  vector<const AliHLTTPCClusterMCData*> fClusterMCData; //! references to MC data blocks
  std::string fHuffmanTableFile; //! file for direct configuration of Huffman table

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

  int formatVersion=1;
  int deflaterMode=2;
  switch (clusterData->fVersion) {
  case 1: deflaterMode=1; formatVersion=0; break;
  case 2: deflaterMode=2; formatVersion=0; break;
  case 3: deflaterMode=1; formatVersion=1; break;
  case 4: deflaterMode=2; formatVersion=1; break;
  default:
    return -EBADF;
  }
  if (!fpDataInflaterPartition)
    fpDataInflaterPartition=CreateInflater(deflaterMode, 1);
  else
    fpDataInflaterPartition->Clear();
  if (!fpDataInflaterPartition) return -ENODEV;

  if ((iResult=fpDataInflaterPartition->InitBitDataInput(reinterpret_cast<const AliHLTUInt8_t*>(clusterData->fClusters),
						size-sizeof(AliHLTTPCRawClusterData)))<0) {
    return iResult;
  }

  iResult=ReadRemainingClustersCompressed(c, fpDataInflaterPartition, nCount, specification, formatVersion);
  if (iResult<0) {
    HLTError("cluster decoding of block 0x%08x failed, error %d", specification, iResult);
  }

  return iResult;
}

template<typename T>
int AliHLTTPCDataCompressionDecoder::ReadRemainingClustersCompressed(T& c, AliHLTDataInflater* pInflater, int nofClusters, AliHLTUInt32_t specification, int formatVersion)
{
  // read cluster data

  int iResult=0;
  if (!pInflater) return -EINVAL;

  if ((iResult= InitPartitionClusterDecoding(specification))<0)
    return iResult;

  int slice = AliHLTTPCDefinitions::GetMinSliceNr(specification);
  int partition = AliHLTTPCDefinitions::GetMinPatchNr(specification);
  if (slice<0 || slice>=AliHLTTPCGeometry::GetNSlice()) {
    HLTError("invalid slice %d decoded from specification 0x%08x", slice, specification);
    return -EINVAL;
  }
  if (partition<0 || partition>=AliHLTTPCGeometry::GetNumberOfPatches()) {
    HLTError("invalid partition %d decoded from specification 0x%08x", partition, specification);
    return -EINVAL;
  }
  // the compressed format stores the difference of the local row number in
  // the partition to the row of the last cluster
  // in padrow mode 'global', add the first row in the partition to get global row number
  int rowOffset=fExtractGlobalPadrow?AliHLTTPCGeometry::GetFirstRow(partition):0;

  int parameterId=pInflater->NextParameter();
  if (parameterId<0) return parameterId;
  int decodedClusterCnt=0;
  int outClusterCnt=0;
  AliHLTUInt64_t value=0;
  AliHLTUInt32_t length=0;
  AliHLTUInt32_t lastPadRow=0;
  AliHLTUInt64_t lastPad64=0;
  AliHLTUInt64_t lastTime64=0;
  AliHLTUInt8_t isSinglePad=0;
  AliHLTUInt8_t sign=0;
  bool bNextCluster=true;
  bool bReadSuccess=true;
  AliHLTTPCRawCluster rawCluster;

  if( fUseClusterMerger && !fpClusterMerger ){
    fpClusterMerger = new AliHLTTPCHWClusterMerger;
  }

  if (fpClusterMerger)
    fpClusterMerger->Clear();

  while (decodedClusterCnt<nofClusters && bReadSuccess && pInflater->NextValue(value, length)) {
    if (bNextCluster) {
      // switch to next cluster
      rawCluster.Clear();
      bNextCluster=false;
    }
    const AliHLTTPCDefinitions::AliClusterParameter& parameter
      =AliHLTTPCDefinitions::fgkClusterParameterDefinitions[parameterId];

    if (parameter.fBitLength!=(int)length) {
      HLTError("decode error: expecting length %d for parameter %s, but got %d",
	       parameter.fBitLength, parameter.fName, length);
      break;
    }

    switch (parameterId) {
    case AliHLTTPCDefinitions::kPadRow:
      {rawCluster.SetPadRow(value+lastPadRow); lastPadRow+=value;break;}
    case AliHLTTPCDefinitions::kPad: {
      if (formatVersion==1) {
	bReadSuccess=bReadSuccess && pInflater->InputBit(isSinglePad);
	if (isSinglePad==0) {
	  bReadSuccess=bReadSuccess && pInflater->InputBit(sign);
	  if (sign) {
	    value=lastPad64-value;
	  } else {
	    value+=lastPad64;
	  }
	  lastPad64=value;
	}
      }
      float pad=value;
      if (isSinglePad==0) pad/=parameter.fScale;
      else pad/=2; // for the sake of the 0.5 pad offset (see AliHLTTPCHWCFSpacePointContainer::WriteSorted for details)
      rawCluster.SetPad(pad+PadShift());
      break;
    }
    case AliHLTTPCDefinitions::kTime: {
      if (formatVersion==1) {
	bReadSuccess=bReadSuccess && pInflater->InputBit(sign);
	if (sign) {
	  value=lastTime64-value;
	} else {
	  value+=lastTime64;
	}
	lastTime64=value;
      }
      float time=value; time/=parameter.fScale;
      rawCluster.SetTime(time);
      break;
    }
    case AliHLTTPCDefinitions::kSigmaY2:
      {float sigmaY2=value; sigmaY2/=parameter.fScale; rawCluster.SetSigmaPad2(sigmaY2); break;}
    case AliHLTTPCDefinitions::kSigmaZ2:
      {float sigmaZ2=value; sigmaZ2/=parameter.fScale; rawCluster.SetSigmaTime2(sigmaZ2); break;}
    case AliHLTTPCDefinitions::kCharge:
      {rawCluster.SetCharge(value); break;}
    case AliHLTTPCDefinitions::kQMax:
      {rawCluster.SetQMax(value); break;}
    }
    if (parameterId>=AliHLTTPCDefinitions::kLast) {
      AliHLTUInt32_t id=GetClusterId(decodedClusterCnt);
      unsigned short flags = GetNextClusterFlag();
      rawCluster.SetFlags(flags);
      const AliHLTTPCClusterMCLabel* pMC=GetMCLabel(id);
      if (fUseClusterMerger && fpClusterMerger && fpClusterMerger->CheckCandidate(slice, partition, rawCluster)) {
	fpClusterMerger->AddCandidate(slice, partition, id, rawCluster, pMC);
      } else {
      c.Next(slice, partition);
      c.SetPadRow(rawCluster.GetPadRow()+rowOffset);
      c.SetPad(rawCluster.GetPad());
      c.SetTime(rawCluster.GetTime());
      c.SetSigmaY2(rawCluster.GetSigmaPad2());
      c.SetSigmaZ2(rawCluster.GetSigmaTime2());
      c.SetCharge(rawCluster.GetCharge());
      c.SetQMax(rawCluster.GetQMax());
      c.SetFlags(rawCluster.GetFlags());
      if (pMC) c.SetMC(pMC);
      outClusterCnt++;
      }
      bNextCluster=true;
      decodedClusterCnt++;
    }
    parameterId=pInflater->NextParameter();
    if (parameterId==AliHLTTPCDefinitions::kSigmaY2 && isSinglePad==1) {
      // skip sigmaY for single pad clusters in format version 1
      parameterId=pInflater->NextParameter();
      isSinglePad=0;
      rawCluster.SetSigmaPad2(0.);
    }
  }
  pInflater->Pad8Bits();
  AliHLTUInt8_t bit=0;
  if (pInflater->InputBit(bit)) {
    HLTWarning("format error of compressed clusters, there is more data than expected");
  }
  pInflater->CloseBitDataInput();
  int mergedClusterCnt=0;
  if (fUseClusterMerger && fpClusterMerger) {
    mergedClusterCnt=fpClusterMerger->Merge();
    int remainingCnt=0;
    if (mergedClusterCnt>=0) {
      for (AliHLTTPCHWClusterMerger::iterator i=fpClusterMerger->begin();
	   i!=fpClusterMerger->end(); i++) {
	c.Next((*i).GetSlice(), (*i).GetPartition());
	const AliHLTTPCRawCluster& mergedCluster=(*i).GetCluster();
	const AliHLTTPCClusterMCLabel& mc=(*i).MCLabel();
	c.SetPadRow(mergedCluster.GetPadRow()+rowOffset);
	c.SetPad(mergedCluster.GetPad());
	c.SetTime(mergedCluster.GetTime());
	c.SetSigmaY2(mergedCluster.GetSigmaPad2());
	c.SetSigmaZ2(mergedCluster.GetSigmaTime2());
	c.SetCharge(mergedCluster.GetCharge());
	c.SetQMax(mergedCluster.GetQMax());
	c.SetFlags(mergedCluster.GetFlags());
	c.SetMC(&mc);
	outClusterCnt++;
	remainingCnt++;
      }
    } else {
      iResult=mergedClusterCnt;
    }
    HLTDebug("copied %d cluster(s) from merger, %d merged, specification 0x%08x", remainingCnt, mergedClusterCnt, specification);
    fpClusterMerger->Clear();
  }
  if (iResult>=0 && nofClusters!=outClusterCnt+mergedClusterCnt) {
    // is this a Fatal?
    HLTError("error reading compressed cluster format of block 0x%08x: expected %d, read %d cluster(s), merged %d cluster(s)", specification, nofClusters, outClusterCnt, mergedClusterCnt);
    return -EPROTO;
  }
  if (iResult<0) return iResult;
  return outClusterCnt;
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
  if (!fpDataInflaterTrack)
    fpDataInflaterTrack=CreateInflater(trackModelBlock->fDeflaterMode, 2);
  else
    fpDataInflaterTrack->Clear();
  if (!fpDataInflaterTrack) {
    HLTError("failed to create the data inflater for mode %d", trackModelBlock->fDeflaterMode);
    return -ENODEV;
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
    if ((iResult=fpDataInflaterTrack->InitBitDataInput(pData+dataOffset, clusterBlockSize))<0) {
      return iResult;
    }
    if ((iResult=ReadTrackClustersCompressed(c, fpDataInflaterTrack, &trackpoints))<0) {
      HLTError("reading of associated clusters failed for track %d", trackno);
      return iResult;
    }
    fpDataInflaterTrack->Pad8Bits();
    AliHLTUInt8_t bit=0;
    if (fpDataInflaterTrack->InputBit(bit)) {
      HLTWarning("format error of compressed clusters, there is more data than expected");
    }
    fpDataInflaterTrack->CloseBitDataInput();
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

  if ((iResult= InitTrackModelClusterClusterDecoding())<0)
    return iResult;

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
      if (AliHLTTPCGeometry::GetFirstRow(AliHLTTPCSpacePointData::GetPatch(currentTrackPoint->GetId())) +
	  AliHLTTPCSpacePointData::GetNumber(currentTrackPoint->GetId()) == row) {
	break;
      }
      currentTrackPoint++;
    } while (!bEscape);
    if (currentTrackPoint==rawTrackPoints.end()) {
      HLTError("decoding error, can not find track point on row %d", row);
      return -EFAULT;
    }
    AliHLTUInt8_t slice = AliHLTTPCSpacePointData::GetSlice(currentTrackPoint->GetId());
    AliHLTUInt8_t partition = AliHLTTPCSpacePointData::GetPatch(currentTrackPoint->GetId());
    // subtract first row of partition if padrow mode 'local'
    unsigned firstRow=fExtractGlobalPadrow?0:AliHLTTPCGeometry::GetFirstRow(partition);
    AliHLTUInt8_t nofClusters=0;
    bReadSuccess=bReadSuccess && pInflater->InputBits(nofClusters, clusterCountBitLength);
    if (!bReadSuccess) break;
    HLTDebug("slice %02d partition %d row %03d: %d cluster(s)", slice, partition, row, nofClusters);

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
    bool bNextCluster=true;
    while (bReadSuccess && inClusterCnt<nofClusters && pInflater->NextValue(value, length)) {
      if (bNextCluster) {
	// switch to next cluster
	c.Next(slice, partition);
	c.SetPadRow(row-firstRow);
	bNextCluster=false;
      }
      const AliHLTTPCDefinitions::AliClusterParameter& parameter
	=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[kParameterIdMapping[parameterId]];

      if (parameter.fBitLength!=(int)length) {
	HLTError("decode error: expecting length %d for parameter %s, but got %d",
		 parameter.fBitLength, parameter.fName, length);
	break;
      }

      bool lastParameter=false;
      switch (kParameterIdMapping[parameterId]) {
      case AliHLTTPCDefinitions::kResidualPad:
	{
	  AliHLTUInt8_t sign=0;
	  bReadSuccess=bReadSuccess && pInflater->InputBit(sign);
	  AliHLTUInt64_t trackpad64=0;
	  double trackpad=currentTrackPoint->GetU();
	  trackpad*=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kResidualPad].fScale;
	  if (currentTrackPoint->GetU()>0.) trackpad64=(AliHLTUInt64_t)round(trackpad);
	  if (sign) {
	    if (trackpad64<value) {
	      bReadSuccess=false;
	      break;
	    }
	    value=trackpad64-value;
	  } else {
	    value+=trackpad64;
	  }
	  float pad=((float)value)/parameter.fScale;
	  c.SetPad(pad+PadShift()); 
	  break;
	}
      case AliHLTTPCDefinitions::kResidualTime:
	{
	  AliHLTUInt8_t sign=0;
	  bReadSuccess=bReadSuccess && pInflater->InputBit(sign);
	  AliHLTUInt64_t tracktime64=0;
	  double tracktime=currentTrackPoint->GetV();
	  tracktime*=AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kResidualTime].fScale;
	  if (currentTrackPoint->GetV()>0.) tracktime64=(AliHLTUInt64_t)round(tracktime);
	  if (sign) {
	    if (tracktime64<value) {
	      bReadSuccess=false;
	      break;
	    }
	    value=tracktime64-value;
	  } else {
	    value+=tracktime64;
	  }
	  float time=((float)value)/parameter.fScale;
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
	bNextCluster=true;
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

template<typename T>
int AliHLTTPCDataCompressionDecoder::ReadClustersPartition(T& c, const AliHLTUInt8_t* pData, unsigned dataSize, AliHLTUInt32_t specification)
{
  // read raw cluster data
  if (!pData) return -EINVAL;
  if (dataSize<sizeof(AliHLTTPCRawClusterData)) return -ENODATA;
  const AliHLTTPCRawClusterData* clusterData = reinterpret_cast<const AliHLTTPCRawClusterData*>(pData);
  Int_t nCount = (Int_t) clusterData->fCount;
  if (clusterData->fVersion!=0) {
    int iResult=ReadRemainingClustersCompressed(c, pData, dataSize, specification);
    if (iResult>=0 && fVerbosity>0) {
      HLTInfo("extracted %d cluster(s) from block 0x%08x", iResult, specification);
    }
    if (iResult<0) {
      HLTError("cluster decoding of block 0x%08x failed, error %d", specification, iResult);
    }
    return iResult;
  }
  if (nCount*sizeof(AliHLTTPCRawCluster) + sizeof(AliHLTTPCRawClusterData) != dataSize) return -EBADF;
  AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetSingleSliceNr(specification);
  AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetSinglePatchNr(specification);
  if (slice < 0 || partition < 0) {
    HLTError("can not decode cluster block of multiple partitions, specification 0x%08x", specification);
    return -1;
  }
  // in padrow mode 'global', add the first row in the partition to get global row number
  int rowOffset=fExtractGlobalPadrow?AliHLTTPCGeometry::GetFirstRow(partition):0;

  const AliHLTTPCRawCluster *clusters = clusterData->fClusters;
  for (int i=0; i<nCount; i++) {
    AliHLTUInt32_t id=GetClusterId(i);
    const AliHLTTPCClusterMCLabel* pMC=GetMCLabel(id);
    c.Next(slice, partition);
    c.SetPadRow(clusters[i].GetPadRow() + rowOffset);
    c.SetPad(clusters[i].GetPad()+PadShift());
    c.SetTime(clusters[i].GetTime());
    c.SetSigmaY2(clusters[i].GetSigmaPad2());
    c.SetSigmaZ2(clusters[i].GetSigmaTime2());
    c.SetCharge(clusters[i].GetCharge());
    c.SetQMax(clusters[i].GetQMax());
    if (pMC) c.SetMC(pMC);
  }
  return nCount;
}
#endif
