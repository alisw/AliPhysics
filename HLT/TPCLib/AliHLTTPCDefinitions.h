// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCDEFINITIONS_H
#define ALIHLTTPCDEFINITIONS_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliHLTDataTypes.h"
#include "Rtypes.h"

/**
 * @class AliHLTTPCDefinitions
 * Data type definitions for the libAliHLTTPC library.
 * 
 * @ingroup alihlt_tpc
 */
class AliHLTTPCDefinitions
{
public:
      AliHLTTPCDefinitions();
      virtual ~AliHLTTPCDefinitions();

	static AliHLTUInt8_t GetMinSliceNr( const AliHLTComponentBlockData& block )
		{
		return (AliHLTUInt8_t)( (block.fSpecification & 0x00FF0000) >> 16 );
		}
	static AliHLTUInt8_t GetMinSliceNr( ULong_t spec )
		{
		return (AliHLTUInt8_t)( (spec & 0x00FF0000) >> 16 );
		}
	static AliHLTUInt8_t GetMaxSliceNr( const AliHLTComponentBlockData& block )
		{
		return (AliHLTUInt8_t)( (block.fSpecification & 0xFF000000) >> 24 );
		}
	static AliHLTUInt8_t GetMaxSliceNr( ULong_t spec )
		{
		return (AliHLTUInt8_t)( (spec & 0xFF000000) >> 24 );
		}
	static AliHLTUInt8_t GetMinPatchNr( const AliHLTComponentBlockData& block )
		{
		return (AliHLTUInt8_t)( (block.fSpecification & 0x000000FF) );
		}
	static AliHLTUInt8_t GetMinPatchNr( ULong_t spec )
		{
		return (AliHLTUInt8_t)( (spec & 0x000000FF) );
		}
	static AliHLTUInt8_t GetMaxPatchNr( const AliHLTComponentBlockData& block )
		{
		return (AliHLTUInt8_t)( (block.fSpecification & 0x0000FF00) >> 8 );
		}
	static AliHLTUInt8_t GetMaxPatchNr( ULong_t spec )
		{
		return (AliHLTUInt8_t)( (spec & 0x0000FF00) >> 8 );
		}
	
	static AliHLTUInt32_t EncodeDataSpecification( AliHLTUInt8_t minSliceNr, 
						AliHLTUInt8_t maxSliceNr,
						AliHLTUInt8_t minPatchNr,
						AliHLTUInt8_t maxPatchNr )
		{
		return ((maxSliceNr & 0xFF) << 24) | ((minSliceNr & 0xFF) << 16) | ((maxPatchNr & 0xFF) << 8) | ((minPatchNr & 0xFF));
		}

  /** DDL packed RAW data */
  static const AliHLTComponentDataType fgkDDLPackedRawDataType;         // see above
  /** DDL entropy encoded data */
  static const AliHLTComponentDataType fgkDDLEncodedEntropyRawDataType; // see above
  /** packed RAW data */
  static const AliHLTComponentDataType fgkPackedRawDataType;            // see above
  /** unpacked RAW data */
  static const AliHLTComponentDataType fgkUnpackedRawDataType;          // see above
  /** cluster data */
  static const AliHLTComponentDataType fgkClustersDataType;             // see above
  /** track segments in local coordinates */
  static const AliHLTComponentDataType fgkTrackSegmentsDataType;        // see above
  /** tracks in global koordinates */
  static const AliHLTComponentDataType fgkTracksDataType;               // see above
  /** vertex data structure */
  static const AliHLTComponentDataType fgkVertexDataType;               // see above

  // Cluster & Tracks model data
  /** cluster tracks model data type */
  static const AliHLTComponentDataType fgkClusterTracksModelDataType;          // see above
  /** remaining clusters model data type */
  static const AliHLTComponentDataType fgkRemainingClustersModelDataType;      // see above
  /** cluster tracks compressed data type */
  static const AliHLTComponentDataType fgkClusterTracksCompressedDataType;     // see above
  /** remaining clusters compressed data type */
  static const AliHLTComponentDataType fgkRemainingClustersCompressedDataType; // see above

  // Calibration data
  /** pedestal calibration data */
  static const AliHLTComponentDataType fgkCalibPedestalDataType;   // see above
  /** signal calibration data */
  static const AliHLTComponentDataType fgkCalibPulserDataType;     // see above
  /** central electrode calibration data */
  static const AliHLTComponentDataType fgkCalibCEDataType;         // see above

  // offline calbration components

  /** alignment calibration data */
  static const AliHLTComponentDataType fgkOfflineCalibAlignDataType;         // see above
  /** track calibration data */
  static const AliHLTComponentDataType fgkOfflineCalibTracksDataType;        // see above
  /** gain calibration data */
  static const AliHLTComponentDataType fgkOfflineCalibTracksGainDataType;    // see above
  /** cluster monte carlo information */
  static const AliHLTComponentDataType fgkAliHLTDataTypeClusterMCInfo;    // see above

  ClassDef(AliHLTTPCDefinitions, 3)
};

#endif
