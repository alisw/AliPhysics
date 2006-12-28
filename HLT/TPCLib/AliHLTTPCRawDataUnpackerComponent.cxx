// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          for The ALICE Off-line Project.                               *
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
// a TPC cluster finder processing component for the HLT                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTTPCRawDataUnpackerComponent.h"
#include "AliTPCRawStream.h"
#include "AliRawDataHeader.h"
#include "AliRawReaderMemory.h"
#include "AliHLTTPCRawDataFormat.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCTransform.h"
#include <stdlib.h>
#include <errno.h>

// this is a global object used for automatic component registration, do not use this
AliHLTTPCRawDataUnpackerComponent gAliHLTTPCRawDataUnpackerComponent;

ClassImp(AliHLTTPCRawDataUnpackerComponent)

AliHLTTPCRawDataUnpackerComponent::AliHLTTPCRawDataUnpackerComponent()
    {
    fRawMemoryReader = NULL;
    fTPCRawStream = NULL;
    }

AliHLTTPCRawDataUnpackerComponent::~AliHLTTPCRawDataUnpackerComponent()
    {
    }

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCRawDataUnpackerComponent::GetComponentID()
    {
    return "TPCRawDataUnpacker";
    }

void AliHLTTPCRawDataUnpackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
    {
    list.clear();
    list.push_back( AliHLTTPCDefinitions::gkDDLPackedRawDataType );
    }

AliHLTComponentDataType AliHLTTPCRawDataUnpackerComponent::GetOutputDataType()
    {
    return AliHLTTPCDefinitions::gkUnpackedRawDataType;
    }

void AliHLTTPCRawDataUnpackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
    // XXX TODO: Find more realistic values.
    constBase = 0;
    inputMultiplier = 6.0;
    }

AliHLTComponent* AliHLTTPCRawDataUnpackerComponent::Spawn()
    {
    return new AliHLTTPCRawDataUnpackerComponent;
    }
	
int AliHLTTPCRawDataUnpackerComponent::DoInit( int argc, const char** argv )
    {
    if ( fRawMemoryReader || fTPCRawStream )
	return EINPROGRESS;

    int i = 0;
//     const char* tableFileBaseDir = NULL;
    while ( i < argc )
	{
// 	if ( !strcmp( argv[i], "table-dir" ) )
// 	    {
// 	    if ( i+1>=argc )
// 		{
// 		Logging(kHLTLogError, "HLT::TPCRawDataUnpacker::DoInit", "Missing Argument", "Missing -table-dir parameter");
// 		return ENOTSUP;
// 		}
// 	    tableFileBaseDir = argv[i+1];
// 	    i += 2;
// 	    continue;
// 	    }
	Logging(kHLTLogError, "HLT::TPCRawDataUnpacker::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
	return EINVAL;
	}

//     if ( !tableFileBaseDir )
// 	{
// 	Logging(kHLTLogError, "HLT::TPCRawDataUnpacker::DoInit", "Table file directory missing ", 
// 		"Table file base directory has to be set using the -table-dir component parameter." );
// 	return EINVAL;
// 	}

    fRawMemoryReader = new AliRawReaderMemory;
    //fTPCRawStream = new AliTPCRawStream( fRawMemoryReader, tableFileBaseDir );
    fTPCRawStream = new AliTPCRawStream( fRawMemoryReader );
    
    return 0;
    }

int AliHLTTPCRawDataUnpackerComponent::DoDeinit()
    {
    if ( fRawMemoryReader )
	delete fRawMemoryReader;
    fRawMemoryReader = NULL;
    if ( fTPCRawStream )
	delete fTPCRawStream;
    fTPCRawStream = NULL;
    return 0;
    }

int AliHLTTPCRawDataUnpackerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
    {
    const AliHLTComponentBlockData* iter = NULL;
    unsigned long ndx;

    AliHLTUInt8_t* outBPtr;
    AliHLTTPCUnpackedRawData* outPtr;
    AliHLTTPCDigitRowData* currentRow;
    AliHLTTPCDigitData* currentDigit;
    unsigned long long outputSize = 0;
    unsigned long blockOutputSize = 0;
    unsigned long rowSize = 0;
    Int_t slice, patch, rows[2];
    outBPtr = outputPtr;
    outPtr = (AliHLTTPCUnpackedRawData*)outputPtr;
    currentRow = outPtr->fDigits;
    currentDigit = currentRow->fDigitData;
    
    Logging( kHLTLogDebug, "HLT::TPCRawDataUnpackerSubscriber::DoEvent", "Event received", 
	     "Event 0x%08LX (%Lu) received with %lu blocks. Output data size: %lu",
	     evtData.fEventID, evtData.fEventID, evtData.fBlockCnt, size);
    for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
	{
	iter = blocks+ndx;
	char tmp1[14], tmp2[14];
	DataType2Text( iter->fDataType, tmp1 );
	DataType2Text( AliHLTTPCDefinitions::gkDDLPackedRawDataType, tmp2 );
	Logging( kHLTLogDebug, "HLT::TPCRawDataUnpackerSubscriber::DoEvent", "Event received", 
		 "Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
		 evtData.fEventID, evtData.fEventID, tmp1, tmp2 );
	if ( iter->fDataType != AliHLTTPCDefinitions::gkDDLPackedRawDataType )
	    {
	    continue;
	    }
	slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
	patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
	rows[0] = AliHLTTPCTransform::GetFirstRow( patch );
	rows[1] = AliHLTTPCTransform::GetLastRow( patch );
	blockOutputSize = 0;

	Logging( kHLTLogDebug, "HLT::TPCRawDataUnpackerSubscriber::DoEvent", "Input Raw Packed Data", 
		 "Input: Slice/Patch/RowMin/RowMax: %d/%d/%d/%d.",
		 slice, patch, rows[0], rows[1] );

	fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize );
	bool readValue = true;
	readValue = fTPCRawStream->Next();
	int row = -1, oldRow = -1;
	Int_t rowOffset = 0;
	if ( patch >= 2 ) // Outer sector, patches 2, 3, 4, 5
	    rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );

	while ( readValue )
	    {
	    row = fTPCRawStream->GetRow();
	    if ( row != oldRow )
		{
		if ( oldRow!=-1 && rowSize != sizeof(AliHLTTPCDigitRowData)+currentRow->fNDigit*sizeof(AliHLTTPCDigitData) )
		    {
		    Logging( kHLTLogFatal, "HLT::TPCRawDataUnpackerSubscriber::DoEvent", "Size inconsistency", 
			     "Size inconsistency for row %d data: %lu != %lu (%lu digits).", oldRow, rowSize, 
			     sizeof(AliHLTTPCDigitRowData)+currentRow->fNDigit*sizeof(AliHLTTPCDigitData), currentRow->fNDigit );
		    }
		rowSize = 0;
		if ( size < outputSize+sizeof(AliHLTTPCDigitRowData) )
		    {
		    Logging( kHLTLogFatal, "HLT::TPCRawDataUnpackerSubscriber::DoEvent", "Too much data", 
			     "Output data too big, output memory full. Aborting event 0x%08lX (%lu)" , 
			     evtData.fEventID, evtData.fEventID );
		    return 0;
		    }
		currentRow = (AliHLTTPCDigitRowData*)(outBPtr+outputSize);
		currentDigit = currentRow->fDigitData;
		currentRow->fRow = row+rowOffset;
		currentRow->fNDigit = 0;
		oldRow = row;
		outputSize += sizeof(AliHLTTPCDigitRowData);
		blockOutputSize += sizeof(AliHLTTPCDigitRowData);
		rowSize += sizeof(AliHLTTPCDigitRowData);
		}
	    if ( size < outputSize+sizeof(AliHLTTPCDigitData) )
		{
		Logging( kHLTLogFatal, "HLT::TPCRawDataUnpackerSubscriber::DoEvent", "Too much data", 
			 "Output data too big, output memory full. Aborting event 0x%08lX (%lu)" , 
			 evtData.fEventID, evtData.fEventID );
		return 0;
		}
	    currentDigit->fCharge = fTPCRawStream->GetSignal();
	    currentDigit->fPad = fTPCRawStream->GetPad();
	    currentDigit->fTime = fTPCRawStream->GetTime();
	    currentRow->fNDigit++;
	    currentDigit++;
	    outputSize += sizeof(AliHLTTPCDigitData);
	    blockOutputSize += sizeof(AliHLTTPCDigitData);
	    rowSize += sizeof(AliHLTTPCDigitData);
	    readValue = fTPCRawStream->Next();
	    }
	if ( oldRow!=-1 && rowSize != sizeof(AliHLTTPCDigitRowData)+currentRow->fNDigit*sizeof(AliHLTTPCDigitData) )
	    {
	    Logging( kHLTLogFatal, "HLT::TPCRawDataUnpackerSubscriber::DoEvent", "Size inconsistency", 
		     "Size inconsistency for row %d data: %lu != %lu (%lu digits).", oldRow, rowSize, 
		     sizeof(AliHLTTPCDigitRowData)+currentRow->fNDigit*sizeof(AliHLTTPCDigitData), currentRow->fNDigit );
	    }

	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset = outputSize-blockOutputSize;
	bd.fSize = blockOutputSize;
	bd.fSpecification = iter->fSpecification;
	Logging( kHLTLogDebug, "HLT::TPCRawDataUnpackerSubscriber::DoEvent", "Event received", 
		 "Event 0x%08LX (%Lu) output data block %lu of %lu bytes at offset %lu",
		 evtData.fEventID, evtData.fEventID, ndx, blockOutputSize, outputSize-blockOutputSize );
	//AliHLTSubEventDescriptor::FillBlockAttributes( bd.fAttributes );
	outputBlocks.push_back( bd );
	}

    
    size = outputSize;
    return 0;
    }

	
