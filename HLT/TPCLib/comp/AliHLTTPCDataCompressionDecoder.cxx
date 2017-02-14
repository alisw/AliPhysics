// $Id$
//**************************************************************************
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTPCDataCompressionDecoder.cxx
/// @author Matthias Richter
/// @date   2011-10-04
/// @brief  Generic decoder class for compressed TPC data, works on a container
///         class implementation which fills the actual target data struct

#include "AliHLTTPCDataCompressionDecoder.h"
#include "AliHLTTPCDataCompressionDescriptor.h"
#include "AliHLTTPCRawClustersDescriptor.h"
#include "AliHLTDataInflaterSimple.h"
#include "AliHLTDataInflaterHuffman.h"
#include "AliCDBEntry.h"
#include "TList.h"
#include "TFile.h"
#include <memory>

ClassImp(AliHLTTPCDataCompressionDecoder)

AliHLTTPCDataCompressionDecoder::AliHLTTPCDataCompressionDecoder()
  : fPadShift(0.)
  , fVerbosity(0)
  , fUseClusterMerger(kTRUE)
  , fExtractGlobalPadrow(kTRUE)
  , fpDataInflaterPartition(NULL)
  , fpDataInflaterTrack(NULL)
  , fpClusterMerger(NULL)
  , fPartitionClusterIds()
  , fPartitionClusterFlags()
  , fTrackModelClusterIds()
  , fCurrentClusterIds(NULL)
  , fCurrentClusterFlags(NULL)
  , fDecodeFlagsTmpFlag(0)
  , fDecodeFlagsTmpPos(0)
  , fDecodeFlagsTmpEntries(0)
  , fClusterMCData()
  , fHuffmanTableFile()
{
  /// constructor
}

AliHLTTPCDataCompressionDecoder::~AliHLTTPCDataCompressionDecoder()
{
  ///destructor
  if (fpDataInflaterPartition) delete fpDataInflaterPartition;
  fpDataInflaterPartition=NULL;
  if (fpDataInflaterTrack) delete fpDataInflaterTrack;
  fpDataInflaterTrack=NULL;
  if (fpClusterMerger) delete fpClusterMerger;
  fpClusterMerger=NULL;
}

AliHLTDataInflater* AliHLTTPCDataCompressionDecoder::CreateInflater(int deflater, int mode) const
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
      TObject* pConf=NULL;
      if (fHuffmanTableFile.empty()) {
      TString cdbPath("HLT/ConfigTPC/TPCDataCompressorHuffmanTables");
      pConf=AliHLTMisc::Instance().ExtractObject(AliHLTMisc::Instance().LoadOCDBEntry(cdbPath));
      if (!pConf) {
	HLTError("can not load configuration object %s", cdbPath.Data());
	return NULL;
      }
      } else {
	// load huffman table directly from file
	TFile* tablefile = TFile::Open(fHuffmanTableFile.c_str());
	if (!tablefile || tablefile->IsZombie()) return NULL;
	TObject* obj = NULL;
	AliCDBEntry* cdbentry = NULL;
	tablefile->GetObject("AliCDBEntry", obj);
	if (obj == NULL || (cdbentry = dynamic_cast<AliCDBEntry*>(obj))==NULL) {
	  HLTError("can not read configuration object from file %s", fHuffmanTableFile.c_str());
	  return NULL;
	}
	HLTInfo("reading huffman table configuration object from file %s", fHuffmanTableFile.c_str());
	pConf = cdbentry->GetObject();
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

int AliHLTTPCDataCompressionDecoder::InitPartitionClusterDecoding(AliHLTUInt32_t specification)
{
  /// init the decoding of partition cluster block
  AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(specification);
  AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(specification);
  unsigned index=slice*AliHLTTPCGeometry::GetNumberOfPatches()+partition;
  if (index<fPartitionClusterIds.size())
    fCurrentClusterIds=&fPartitionClusterIds[index];
  else
    fCurrentClusterIds=NULL;

  if (index<fPartitionClusterFlags.size())
    fCurrentClusterFlags=fPartitionClusterFlags[index];
  else
    fCurrentClusterFlags=NULL;

  fDecodeFlagsTmpFlag=0;
  fDecodeFlagsTmpPos=0;
  fDecodeFlagsTmpEntries=0;

  return 0;
}

int AliHLTTPCDataCompressionDecoder::InitTrackModelClusterClusterDecoding()
{
  /// init the decoding of track model cluster block
  if (fTrackModelClusterIds.fIds && fTrackModelClusterIds.fSize>0)
    fCurrentClusterIds=&fTrackModelClusterIds;
  else
    fCurrentClusterIds=NULL;

  fCurrentClusterFlags=NULL;
  fDecodeFlagsTmpFlag=0;
  fDecodeFlagsTmpPos=0;
  fDecodeFlagsTmpEntries=0;
  return 0;
}

int AliHLTTPCDataCompressionDecoder::AddCompressionDescriptor(const AliHLTComponentBlockData* pDesc)
{
  /// read descriptor
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType!=AliHLTTPCDefinitions::DataCompressionDescriptorDataType()) return -ENODATA;
  const AliHLTTPCDataCompressionDescriptor* pHeader=reinterpret_cast<const AliHLTTPCDataCompressionDescriptor*>(pDesc->fPtr);
  if (! pHeader->CheckSize( pDesc->fSize ) ) return -EINVAL;
  if( pHeader->GetMergedClustersFlag() == 0 ){
    fUseClusterMerger = kTRUE;
  } else if( pHeader->GetMergedClustersFlag() == 1 ){
    fUseClusterMerger = kFALSE;
  } else return -EINVAL;
  return 0;
}

int AliHLTTPCDataCompressionDecoder::AddRawClustersDescriptor(const AliHLTComponentBlockData* pDesc)
{
  /// read descriptor
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType!=AliHLTTPCDefinitions::RawClustersDescriptorDataType()) return -ENODATA;
  const AliHLTTPCRawClustersDescriptor* pHeader=reinterpret_cast<const AliHLTTPCRawClustersDescriptor*>(pDesc->fPtr);
  if (! pHeader->CheckSize( pDesc->fSize ) ) return -EINVAL;
  if( pHeader->GetMergedClustersFlag() == 0 ){
    fUseClusterMerger = kTRUE;
  } else if( pHeader->GetMergedClustersFlag() == 1 ){
    fUseClusterMerger = kFALSE;
  } else return -EINVAL;
  return 0;
}

int AliHLTTPCDataCompressionDecoder::AddClusterMCData(const AliHLTComponentBlockData* pDesc)
{
  /// add cluster mc data block
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType==AliHLTTPCDefinitions::AliHLTDataTypeClusterMCInfo()) {
    AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(pDesc->fSpecification);
    AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(pDesc->fSpecification);
    unsigned index=slice*AliHLTTPCGeometry::GetNumberOfPatches()+partition;
    if (fClusterMCData.size()<=index) {
      if ((int)fClusterMCData.size()<AliHLTTPCGeometry::GetNSlice()*AliHLTTPCGeometry::GetNumberOfPatches()) {
	fClusterMCData.resize(AliHLTTPCGeometry::GetNSlice()*AliHLTTPCGeometry::GetNumberOfPatches(), NULL);
      } else {
	fClusterMCData.resize(index+1, NULL);
      }
    }
    if (pDesc->fSize<sizeof(AliHLTTPCClusterMCData)) return -EINVAL;
    const AliHLTTPCClusterMCData* pData=reinterpret_cast<const AliHLTTPCClusterMCData*>(pDesc->fPtr);
    unsigned nLabels = pData->fCount;
    if (nLabels*sizeof(AliHLTTPCClusterMCLabel) + sizeof(AliHLTTPCClusterMCData) != pDesc->fSize) {
      return -EINVAL;
    }
    fClusterMCData[index]=pData;
    return 0;
  }
  return -ENODATA;
}

int AliHLTTPCDataCompressionDecoder::AddClusterIds(const AliHLTComponentBlockData* pDesc)
{
  /// add cluster id block for partition or track model clusters
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
    fTrackModelClusterIds.fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fTrackModelClusterIds.fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  if (pDesc->fDataType==AliHLTTPCDefinitions::RemainingClusterIdsDataType()) {
    AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(pDesc->fSpecification);
    AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(pDesc->fSpecification);
    unsigned index=slice*AliHLTTPCGeometry::GetNumberOfPatches()+partition;
    if (fPartitionClusterIds.size()<=index) {
      if ((int)fPartitionClusterIds.size()<AliHLTTPCGeometry::GetNSlice()*AliHLTTPCGeometry::GetNumberOfPatches()) {
	fPartitionClusterIds.resize(AliHLTTPCGeometry::GetNSlice()*AliHLTTPCGeometry::GetNumberOfPatches());
      } else {
	fPartitionClusterIds.resize(index+1);
      }
    }
    fPartitionClusterIds[index].fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fPartitionClusterIds[index].fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  return -ENODATA;
}

int AliHLTTPCDataCompressionDecoder::AddClusterFlags(const AliHLTComponentBlockData* pDesc)
{
  /// add cluster flag block for partition or track model clusters
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType==AliHLTTPCDefinitions::ClustersFlagsDataType()) {
    AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(pDesc->fSpecification);
    AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(pDesc->fSpecification);
    unsigned index=slice*AliHLTTPCGeometry::GetNumberOfPatches()+partition;
    if (fPartitionClusterFlags.size()<=index) {
      if ((int)fPartitionClusterFlags.size()<AliHLTTPCGeometry::GetNSlice()*AliHLTTPCGeometry::GetNumberOfPatches()) {
	fPartitionClusterFlags.resize(AliHLTTPCGeometry::GetNSlice()*AliHLTTPCGeometry::GetNumberOfPatches());
      } else {
	fPartitionClusterFlags.resize(index+1);
      }
    }
    fPartitionClusterFlags[index]=reinterpret_cast<AliHLTTPCClusterFlagsData*>(pDesc->fPtr);
    return 0;
  }
  return -ENODATA;
}

AliHLTUInt32_t AliHLTTPCDataCompressionDecoder::GetClusterId(int clusterNo) const
{
  /// get the cluster id from the current cluster id block
  /// clusters ids correctly link the MC label from the separate MC data block
  /// to the cluster. The option is enabled by default in the simulation.
  if (!fCurrentClusterIds ||
      (int)fCurrentClusterIds->fSize<=clusterNo ||
      clusterNo<0)
    return kAliHLTVoidDataSpec;
  return fCurrentClusterIds->fIds[clusterNo];
}

unsigned short AliHLTTPCDataCompressionDecoder::GetNextClusterFlag()
{
  if (!fCurrentClusterFlags) return(0); //No data block with flags present
  if (fCurrentClusterFlags->fVersion == 0) return(0); //Shipping of cluster flags was disabled during reconstruction
  if (fCurrentClusterFlags->fVersion >= 2)
  {
    HLTError("Invalid version for compression of cluster flags. Version in the Data: %d, maximum version supported by this AliRoot Version: %d", fCurrentClusterFlags->fVersion, 1);
    return(0);
  }
  if (fCurrentClusterFlags->fVersion == 1)
  {
    if (fDecodeFlagsTmpPos == 0) fDecodeFlagsTmpFlag = ((unsigned int*) fCurrentClusterFlags->fData)[fDecodeFlagsTmpEntries++];
    unsigned short retVal = (fDecodeFlagsTmpFlag >> fDecodeFlagsTmpPos) & ((1 << fCurrentClusterFlags->fNumberOfFlags) - 1);
    fDecodeFlagsTmpPos += fCurrentClusterFlags->fNumberOfFlags;
    if (fDecodeFlagsTmpPos >= sizeof(fDecodeFlagsTmpFlag) * 8)
    {
      fDecodeFlagsTmpPos -= sizeof(fDecodeFlagsTmpFlag) * 8;
      if (fDecodeFlagsTmpPos)
      {
        fDecodeFlagsTmpFlag = ((unsigned int*) fCurrentClusterFlags->fData)[fDecodeFlagsTmpEntries++];
        retVal |= (fDecodeFlagsTmpFlag << (fCurrentClusterFlags->fNumberOfFlags - fDecodeFlagsTmpPos))& ((1 << fCurrentClusterFlags->fNumberOfFlags) - 1);
      }
    }
    return(retVal);
  }
  HLTFatal("Internal Error"); //We should never reach this code
  return(0);
}

const AliHLTTPCClusterMCLabel* AliHLTTPCDataCompressionDecoder::GetMCLabel(AliHLTUInt32_t clusterId) const
{
  /// get MC data for a cluster Id
  /// MC data is sent in a separate data block to keep the raw compressed
  /// format free from any overhead
  if (clusterId==kAliHLTVoidDataSpec) return NULL;

  unsigned slice=AliHLTTPCSpacePointData::GetSlice(clusterId);
  unsigned partition=AliHLTTPCSpacePointData::GetPatch(clusterId);
  unsigned number=AliHLTTPCSpacePointData::GetNumber(clusterId);
  if ((int)slice>=AliHLTTPCGeometry::GetNSlice() ||
      (int)partition>=AliHLTTPCGeometry::GetNumberOfPatches()) return NULL;
  unsigned index=slice*AliHLTTPCGeometry::GetNumberOfPatches()+partition;
  if (fClusterMCData.size()<=index ||
      fClusterMCData[index]==NULL ||
      fClusterMCData[index]->fCount<=number) return NULL;

  return &(fClusterMCData[index]->fLabels[number]);
}

void AliHLTTPCDataCompressionDecoder::Clear(const char* option)
{
  /// cleanup, tabula rase for next event
  if (fpDataInflaterPartition) fpDataInflaterPartition->Clear(option);
  if (fpDataInflaterTrack) fpDataInflaterTrack->Clear(option);
  if (fpClusterMerger) fpClusterMerger->Clear();
  fCurrentClusterIds=NULL;
  fCurrentClusterFlags=NULL;
  fDecodeFlagsTmpFlag=0;
  fDecodeFlagsTmpPos=0;
  fDecodeFlagsTmpEntries=0;
  fPartitionClusterIds.clear();
  fPartitionClusterFlags.clear();
  fTrackModelClusterIds.Clear();
  fClusterMCData.clear();
  fUseClusterMerger = kTRUE;
}
