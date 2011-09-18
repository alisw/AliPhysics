// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Definitions for the HLT TPC components                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTTPCDefinitions.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCDefinitions)

const AliHLTComponentDataType AliHLTTPCDefinitions::fgkDDLEncodedEntropyRawDataType = AliHLTComponentDataTypeInitializer("DDLENCEN", kAliHLTDataOriginTPC);
										      								  	      
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkPackedRawDataType = AliHLTComponentDataTypeInitializer("RAWPAKED", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkUnpackedRawDataType = AliHLTComponentDataTypeInitializer("RAWUNPAK", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkClustersDataType = AliHLTComponentDataTypeInitializer("CLUSTERS", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkRawClustersDataType = AliHLTComponentDataTypeInitializer("CLUSTRAW", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkHWClustersDataType = AliHLTComponentDataTypeInitializer("HWCLUST1", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkAlterClustersDataType = AliHLTComponentDataTypeInitializer("HWCL_ALT", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkVertexDataType = AliHLTComponentDataTypeInitializer("VERTEX  ", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkTrackSegmentsDataType = AliHLTComponentDataTypeInitializer("TRAKSEGS", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkTracksDataType = AliHLTComponentDataTypeInitializer("TRACKS  ", kAliHLTDataOriginTPC);

const AliHLTComponentDataType AliHLTTPCDefinitions::fgkClusterTracksModelDataType = AliHLTComponentDataTypeInitializer("CLSTRKMD", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkRemainingClustersModelDataType = AliHLTComponentDataTypeInitializer("REMCLSMD", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkClusterTracksCompressedDataType = AliHLTComponentDataTypeInitializer("CLSTRKCM", kAliHLTDataOriginTPC);
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkRemainingClustersCompressedDataType = AliHLTComponentDataTypeInitializer("REMCLSCM", kAliHLTDataOriginTPC);
										      								  	      
const AliHLTComponentDataType& AliHLTTPCDefinitions::DDLEncodedEntropyRawDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("DDLENCEN", kAliHLTDataOriginTPC);
  return dt;
}										      								  	      
const AliHLTComponentDataType& AliHLTTPCDefinitions::PackedRawDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("RAWPAKED", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::UnpackedRawDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("RAWUNPAK", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::ClustersDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CLUSTERS", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::RawClustersDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CLUSTRAW", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::HWClustersDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("HWCLUST1", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::AlterClustersDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("HWCL_ALT", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::VertexDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("VERTEX  ", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::TrackSegmentsDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("TRAKSEGS", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::TracksDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("TRACKS  ", kAliHLTDataOriginTPC);
  return dt;
}

const AliHLTComponentDataType& AliHLTTPCDefinitions::ClusterTracksModelDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CLSTRKMD", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::RemainingClustersModelDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("REMCLSMD", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::ClusterTracksCompressedDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CLSTRKCM", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::RemainingClustersCompressedDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("REMCLSCM", kAliHLTDataOriginTPC);
  return dt;
}

const AliHLTComponentDataType& AliHLTTPCDefinitions::CalibPedestalDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CAL_PED ", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::CalibPulserDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CAL_PULS", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::CalibCEDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CAL_CE  ", kAliHLTDataOriginTPC);
  return dt;
}

const AliHLTComponentDataType& AliHLTTPCDefinitions::OfflineCalibAlignDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CALALIGN", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType& AliHLTTPCDefinitions::OfflineCalibTracksDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CALTRACK", kAliHLTDataOriginTPC);
  return dt;
}

const AliHLTComponentDataType& AliHLTTPCDefinitions::OfflineCalibTracksGainDataType() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CALGAIN ", kAliHLTDataOriginTPC);
  return dt;
}

const AliHLTComponentDataType& AliHLTTPCDefinitions::AliHLTDataTypeClusterMCInfo() {
  static AliHLTComponentDataType dt = AliHLTComponentDataTypeInitializer("CLMCINFO", kAliHLTDataOriginTPC);
  return dt;
}
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo = AliHLTComponentDataTypeInitializer("CLMCINFO", kAliHLTDataOriginTPC);


const AliHLTTPCDefinitions::AliClusterParameter AliHLTTPCDefinitions::fgkClusterParameterDefinitions[]= {
  {AliHLTTPCDefinitions::kPadRow,  "padrow",   6,  1,   1}, // difference of rows, mostly 0 or 1
  {AliHLTTPCDefinitions::kPad,     "pad",     14, 12,  60}, // <100um for 6mm pads
  {AliHLTTPCDefinitions::kTime,    "time",    15, 13,  25}, // <100um for 2.5 mm timebin pitch
  {AliHLTTPCDefinitions::kSigmaY2, "sigmaY2",  8,  5,  25},
  {AliHLTTPCDefinitions::kSigmaZ2, "sigmaZ2",  8,  5,  10},
  {AliHLTTPCDefinitions::kCharge,  "charge",  16,  9,   1},
  {AliHLTTPCDefinitions::kQMax,    "qmax",    10,  6,   1},
  {AliHLTTPCDefinitions::kResidualPad, "respad",         9,  4, 60}, // <100um for 6mm pads, sign stored in separate bit
  {AliHLTTPCDefinitions::kResidualTime,"restime",        8,  4, 25}, // <100um for 2.5 mm timebin pitch, separate bit for sign
  {AliHLTTPCDefinitions::kClusterCount,"clustercount",   6,  3,  1}  // number of clusters on that row
};

unsigned AliHLTTPCDefinitions::GetNumberOfClusterParameterDefinitions()
{
  return sizeof(fgkClusterParameterDefinitions)/sizeof(AliClusterParameter);
}

// NOTE! those values are related to the number of bits in
// fgkClusterParameterDefinitions
const unsigned AliHLTTPCDefinitions::fgkMaxClusterDeltaPad=8;
const unsigned AliHLTTPCDefinitions::fgkMaxClusterDeltaTime=10;

AliHLTTPCDefinitions::AliHLTTPCDefinitions()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCDefinitions::~AliHLTTPCDefinitions()
{
  // see header file for class documentation
}

bool AliHLTTPCDefinitions::DDLIdToSlicePatch(AliHLTInt32_t ddlid, AliHLTUInt8_t& slice, AliHLTUInt8_t& patch)
{
	// Convert DDL ID to patch and slice numbers.
	
	if ((AliHLTUInt32_t(ddlid) >> 8) != 0x3) return false;  // Check that detector is TPC.
	AliHLTUInt32_t ddl = (AliHLTUInt32_t(ddlid) & 0xFF);
	if (ddl > 215) return false;
	if (ddl < 72)
	{
		slice = ddl / 2;
		patch = ddl % 2;
	}
	else
	{
		ddl -= 72;
		slice = ddl / 4;
		patch = ddl % 4 + 2;
	}
	return true;
}
    
