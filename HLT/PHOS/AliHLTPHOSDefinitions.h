// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTPHOSDEFINITIONS_H
#define ALIHLTPHOSDEFINITIONS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTTPCDefinitions
 */

#include "AliHLTDataTypes.h"
//#include "TObject.h"
#include "Rtypes.h"

class AliHLTPHOSDefinitions
    {
    public:

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

	static const AliHLTComponentDataType gkDDLPackedRawDataType;
	static const AliHLTComponentDataType gkPackedRawDataType;
	static const AliHLTComponentDataType gkUnpackedRawDataType;
	static const AliHLTComponentDataType gkClustersDataType;
	static const AliHLTComponentDataType gkTrackSegmentsDataType;
	static const AliHLTComponentDataType gkVertexDataType;

      //	ClassDef(AliHLTPHOSDefinitions, 0)

    };

#endif
