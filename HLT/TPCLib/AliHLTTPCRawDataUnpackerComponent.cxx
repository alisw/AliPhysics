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
#include "AliL3DigitData.h"
#include "AliL3Transform.h"
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

void AliHLTTPCRawDataUnpackerComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
    {
    list.clear();
    list.push_back( AliHLTTPCDefinitions::gkPackedRawDataType );
    }

AliHLTComponent_DataType AliHLTTPCRawDataUnpackerComponent::GetOutputDataType()
    {
    return AliHLTTPCDefinitions::gkUnpackedRawDataType;
    }

void AliHLTTPCRawDataUnpackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
    // XXX TODO: Find more realistic values.
    constBase = 0;
    inputMultiplier = 3.0;
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
    const char* tableFileBaseDir = NULL;
    while ( i < argc )
	{
	if ( !strcmp( argv[i], "-table-dir" ) )
	    {
	    if ( i+1>=argc )
		{
		Logging(kHLTLogError, "HLT::TPCRawDataUnpacker::DoInit", "Missing Argument", "Missing -table-dir parameter");
		return ENOTSUP;
		}
	    tableFileBaseDir = argv[i+1];
	    i += 2;
	    continue;
	    }
	Logging(kHLTLogError, "HLT::TPCRawDataUnpacker::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
	return EINVAL;
	}

    if ( !tableFileBaseDir )
	{
	Logging(kHLTLogError, "HLT::TPCRawDataUnpacker::DoInit", "Table file directory missing ", 
		"Table file base directory has to be set using the -table-dir component parameter." );
	return EINVAL;
	}

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

int AliHLTTPCRawDataUnpackerComponent::DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
					      AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
    {
    const AliHLTComponent_BlockData* iter = NULL;
    unsigned long ndx;

    AliHLTUInt8_t* outBPtr;
    AliHLTTPCUnpackedRawData* outPtr;
    AliL3DigitRowData* currentRow;
    AliL3DigitData* currentDigit;
    unsigned long long outputSize = 0;
    unsigned long blockOutputSize = 0;
    unsigned long rowSize = 0;
    Int_t slice, patch, row[2];
    outBPtr = outputPtr;
    outPtr = (AliHLTTPCUnpackedRawData*)outputPtr;
    currentRow = outPtr->fDigits;
    currentDigit = currentRow->fDigitData;

    for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
	{
	iter = blocks+ndx;
	if ( iter->fDataType != AliHLTTPCDefinitions::gkDDLPackedRawDataType )
	    {
	    continue;
	    }
	slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
	patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
	row[0] = AliL3Transform::GetFirstRow( patch );
	row[1] = AliL3Transform::GetLastRow( patch );
	blockOutputSize = 0;

	fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize );
	bool readValue = true;
	readValue = fTPCRawStream->Next();
	int row = -1, oldRow = -1;
	Int_t rowOffset = 0;
	if ( patch >= 2 ) // Outer sector, patches 2, 3, 4, 5
	    rowOffset = AliL3Transform::GetFirstRow( 2 );

	while ( readValue )
	    {
	    row = fTPCRawStream->GetRow();
	    if ( row != oldRow )
		{
		if ( oldRow!=-1 && rowSize != sizeof(AliL3DigitRowData)+currentRow->fNDigit*sizeof(AliL3DigitData) )
		    {
		    Logging( kHLTLogFatal, "TPCRawDataUnpackerSubscriber::ProcessEvent", "Size inconsistency", 
			     "Size inconsistency for row %d data: %lu != %lu (%lu digits).", oldRow, rowSize, 
			     sizeof(AliL3DigitRowData)+currentRow->fNDigit*sizeof(AliL3DigitData), currentRow->fNDigit );
		    }
		rowSize = 0;
		if ( size < outputSize+sizeof(AliL3DigitRowData) )
		    {
		    Logging( kHLTLogFatal, "TPCRawDataUnpackerSubscriber::ProcessEvent", "Too much data", 
			     "Output data too big, output memory full. Aborting event 0x%08lX (%lu)" , 
			     evtData.fEventID, evtData.fEventID );
		    return 0;
		    }
		currentRow = (AliL3DigitRowData*)(outBPtr+outputSize);
		currentDigit = currentRow->fDigitData;
		currentRow->fRow = row+rowOffset;
		currentRow->fNDigit = 0;
		oldRow = row;
		outputSize += sizeof(AliL3DigitRowData);
		blockOutputSize += sizeof(AliL3DigitRowData);
		rowSize += sizeof(AliL3DigitRowData);
		}
	    if ( size < outputSize+sizeof(AliL3DigitData) )
		{
		Logging( kHLTLogFatal, "TPCRawDataUnpackerSubscriber::ProcessEvent", "Too much data", 
			 "Output data too big, output memory full. Aborting event 0x%08lX (%lu)" , 
			 evtData.fEventID, evtData.fEventID );
		return 0;
		}
	    currentDigit->fCharge = fTPCRawStream->GetSignal();
	    currentDigit->fPad = fTPCRawStream->GetPad();
	    currentDigit->fTime = fTPCRawStream->GetTime();
	    currentRow->fNDigit++;
	    currentDigit++;
	    outputSize += sizeof(AliL3DigitData);
	    blockOutputSize += sizeof(AliL3DigitData);
	    rowSize += sizeof(AliL3DigitData);
	    readValue = fTPCRawStream->Next();
	    }
	if ( oldRow!=-1 && rowSize != sizeof(AliL3DigitRowData)+currentRow->fNDigit*sizeof(AliL3DigitData) )
	    {
	    Logging( kHLTLogFatal, "TPCRawDataUnpackerSubscriber::ProcessEvent", "Size inconsistency", 
		     "Size inconsistency for row %d data: %lu != %lu (%lu digits).", oldRow, rowSize, 
		     sizeof(AliL3DigitRowData)+currentRow->fNDigit*sizeof(AliL3DigitData), currentRow->fNDigit );
	    }

	AliHLTComponent_BlockData bd;
	FillBlockData( bd );
	bd.fOffset = outputSize-blockOutputSize;
	bd.fSize = blockOutputSize;
	bd.fSpecification = iter->fSpecification;
	//AliHLTSubEventDescriptor::FillBlockAttributes( bd.fAttributes );
	outputBlocks.push_back( bd );
	}

    
    size = outputSize;
    return 0;
    }

	
