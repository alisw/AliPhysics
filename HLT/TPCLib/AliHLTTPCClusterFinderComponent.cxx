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

#include "AliHLTTPCClusterFinderComponent.h"
#include "AliL3ClustFinderNew.h"
#include "AliL3SpacePointData.h"
#include "AliHLTTPCRawDataFormat.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliL3Transform.h"
#include <stdlib.h>
#include <errno.h>

// this is a global object used for automatic component registration, do not use this
AliHLTTPCClusterFinderComponent gAliHLTTPCClusterFinderComponent;

ClassImp(AliHLTTPCClusterFinderComponent)

AliHLTTPCClusterFinderComponent::AliHLTTPCClusterFinderComponent()
    {
    fClusterFinder = NULL;
    fClusterDeconv = true;
    fXYClusterError = -1;
    fZClusterError = -1;
    }

AliHLTTPCClusterFinderComponent::~AliHLTTPCClusterFinderComponent()
    {
    }

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCClusterFinderComponent::GetComponentID()
    {
    return "TPCClusterFinder";
    }

void AliHLTTPCClusterFinderComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
    {
    list.clear();
    list.push_back( AliHLTTPCDefinitions::gkUnpackedRawDataType );
    }

AliHLTComponent_DataType AliHLTTPCClusterFinderComponent::GetOutputDataType()
    {
    return AliHLTTPCDefinitions::gkClustersDataType;
    }

void AliHLTTPCClusterFinderComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
    // XXX TODO: Find more realistic values.
    constBase = 0;
    inputMultiplier = 0.4;
    }

AliHLTComponent* AliHLTTPCClusterFinderComponent::Spawn()
    {
    return new AliHLTTPCClusterFinderComponent;
    }
	
int AliHLTTPCClusterFinderComponent::DoInit( int argc, const char** argv )
    {
    if ( fClusterFinder )
	return EINPROGRESS;
    fClusterFinder = new AliL3ClustFinderNew();
    fClusterDeconv = true;
    fXYClusterError = -1;
    fZClusterError = -1;
    int i = 0;
    while ( i < argc )
	{
	if ( !strcmp( argv[i], "-pp-run" ) )
	    {
	    fClusterDeconv = false;
	    i++;
	    continue;
	    }
	Logging(kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
	return EINVAL;
	}
    return 0;
    }

int AliHLTTPCClusterFinderComponent::DoDeinit()
    {
    if ( !fClusterFinder )
	return ECANCELED;
    if ( fClusterFinder )
	delete fClusterFinder;
    fClusterFinder = NULL;
    return 0;
    }

int AliHLTTPCClusterFinderComponent::DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
					      AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
    {
    const AliHLTComponent_BlockData* iter = NULL;
    unsigned long ndx;
    AliHLTTPCUnpackedRawData* inPtr;
    AliHLTTPCClusterData* outPtr;
    AliHLTUInt8_t* outBPtr;
    UInt_t offset, mysize, nSize, tSize = 0;
    outBPtr = outputPtr;
    outPtr = (AliHLTTPCClusterData*)outBPtr;
    Int_t slice, patch, row[2];
    unsigned long maxPoints, realPoints = 0;
    for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
	{
	iter = blocks+ndx;
	mysize = 0;
	offset = tSize;
	if ( iter->fDataType != AliHLTTPCDefinitions::gkUnpackedRawDataType )
	    continue;
	
	slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
	patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
	row[0] = AliL3Transform::GetFirstRow( patch );
	row[1] = AliL3Transform::GetLastRow( patch );
	
	Logging( kHLTLogDebug, "HLT::TPCClusterFinder::DoEvent", "Input Spacepoints", 
		 "Input: Number of spacepoints: %lu Slice/Patch/RowMin/RowMax: %d/%d/%d/%d.",
		 realPoints, slice, patch, row[0], row[1] );
	
	outPtr = (AliHLTTPCClusterData*)outBPtr;
	
	inPtr = (AliHLTTPCUnpackedRawData*)iter->fPtr;
	maxPoints = (size-tSize-sizeof(AliHLTTPCClusterData))/sizeof(AliL3SpacePointData);
	
	fClusterFinder->InitSlice( slice, patch, row[0], row[1], maxPoints );
	fClusterFinder->SetDeconv( fClusterDeconv );
	fClusterFinder->SetXYError( fXYClusterError );
	fClusterFinder->SetZError( fZClusterError );
	if ( (fXYClusterError>0) && (fZClusterError>0) )
	    fClusterFinder->SetCalcErr( false );
	fClusterFinder->SetOutputArray( outPtr->fSpacePoints );
	fClusterFinder->Read( maxPoints, inPtr->fDigits );
	fClusterFinder->ProcessDigits();
	realPoints = fClusterFinder->GetNumberOfClusters();
	
	Logging( kHLTLogDebug, "HLT::TPCClusterFinder::DoEvent", "Spacepoints", 
		 "Number of spacepoints found: %lu.", realPoints );
	
	outPtr->fSpacePointCnt = realPoints;
	nSize = sizeof(AliL3SpacePointData)*realPoints;
	mysize += nSize+sizeof(AliHLTTPCClusterData);
	
	Logging( kHLTLogDebug, "HLT::TPCClusterFinder::DoEvent", "Input Spacepoints", 
		 "Number of spacepoints: %lu Slice/Patch/RowMin/RowMax: %d/%d/%d/%d.",
		 realPoints, slice, patch, row[0], row[1] );
	
	
	AliHLTComponent_BlockData bd;
	FillBlockData( bd );
	bd.fOffset = offset;
	bd.fSize = mysize;
	bd.fSpecification = iter->fSpecification;
	//AliHLTSubEventDescriptor::FillBlockAttributes( bd.fAttributes );
	outputBlocks.push_back( bd );
	
	tSize += mysize;
	outBPtr += mysize;
	outPtr = (AliHLTTPCClusterData*)outBPtr;
	
	if ( tSize > size )
	    {
	    Logging( kHLTLogFatal, "HLT::TPCClusterFinder::DoEvent", "Too much data", 
		     "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
		     tSize, size );
	    return EMSGSIZE;
	    }
	}
    
    size = tSize;
    return 0;
    }

	
