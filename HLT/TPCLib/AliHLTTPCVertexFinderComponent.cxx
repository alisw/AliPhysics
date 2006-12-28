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
// a TPC vertex finder processing component for the HLT                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTTPCVertexFinderComponent.h"
#include "AliHLTTPCVertexFinder.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCVertexData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTransform.h"
#include <stdlib.h>
#include <errno.h>

// this is a global object used for automatic component registration, do not use this
AliHLTTPCVertexFinderComponent gAliHLTTPCVertexFinderComponent;

ClassImp(AliHLTTPCVertexFinderComponent)

AliHLTTPCVertexFinderComponent::AliHLTTPCVertexFinderComponent()
    {
    fVertexFinder = NULL;
    }

AliHLTTPCVertexFinderComponent::~AliHLTTPCVertexFinderComponent()
    {
    }

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCVertexFinderComponent::GetComponentID()
    {
    return "TPCVertexFinder";
    }

void AliHLTTPCVertexFinderComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
    {
    list.clear();
    list.push_back( AliHLTTPCDefinitions::gkClustersDataType );
    }

AliHLTComponentDataType AliHLTTPCVertexFinderComponent::GetOutputDataType()
    {
    return AliHLTTPCDefinitions::gkVertexDataType;
    }

void AliHLTTPCVertexFinderComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
    // XXX TODO: Find more realistic values.
    constBase = sizeof(AliHLTTPCVertexData);
    inputMultiplier = 0;
    }

AliHLTComponent* AliHLTTPCVertexFinderComponent::Spawn()
    {
    return new AliHLTTPCVertexFinderComponent;
    }
	
int AliHLTTPCVertexFinderComponent::DoInit( int argc, const char** argv )
    {
    if ( fVertexFinder )
	return EINPROGRESS;
    fVertexFinder = new AliHLTTPCVertexFinder();
    return 0;
    }

int AliHLTTPCVertexFinderComponent::DoDeinit()
    {
    if ( !fVertexFinder )
	return ECANCELED;
    if ( fVertexFinder )
	delete fVertexFinder;
    fVertexFinder = NULL;
    return 0;
    }

int AliHLTTPCVertexFinderComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
    {
    const AliHLTComponentBlockData* iter = NULL;
    unsigned long ndx;

    AliHLTTPCClusterData* inPtr;
    AliHLTTPCVertexData* outPtr;
    AliHLTUInt8_t* outBPtr;
    UInt_t offset, mysize, tSize = 0;
    outBPtr = outputPtr;
    Int_t slice, patch, row[2];
    AliHLTUInt32_t realPoints;

    for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
	{
	iter = blocks+ndx;
	mysize = 0;
	offset = tSize;
	if ( iter->fDataType != AliHLTTPCDefinitions::gkClustersDataType )
	    {
	    continue;
	    }
	
	inPtr = (AliHLTTPCClusterData*)(iter->fPtr);
	slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
	patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
	row[0] = AliHLTTPCTransform::GetFirstRow( patch );
	row[1] = AliHLTTPCTransform::GetLastRow( patch );
	realPoints = inPtr->fSpacePointCnt;

	Logging( kHLTLogDebug, "HLT::TPCVertexFinder::DoEvent", "Spacepoint count",
		 "realpoints: %lu.", realPoints );
	
	outPtr = (AliHLTTPCVertexData*)outBPtr;

	fVertexFinder->Reset();
	
        fVertexFinder->Read( realPoints, inPtr->fSpacePoints );
        fVertexFinder->Analyze();

        //publish Vertex
        fVertexFinder->Write( outPtr );


	mysize += sizeof(AliHLTTPCVertexData);
	
	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset = offset;
	bd.fSize = mysize;
	bd.fSpecification = iter->fSpecification;
	//AliHLTSubEventDescriptor::FillBlockAttributes( bd.fAttributes );
	outputBlocks.push_back( bd );

	tSize += mysize;
	outBPtr += mysize;

	if ( tSize > size )
	    {
	    Logging( kHLTLogFatal, "HLT::TPCVertexFinder::DoEvent", "Too much data",
		     "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
		     , tSize, size );
	    return EMSGSIZE;
	    }
	}
    
    size = tSize;
    return 0;
    }

	
