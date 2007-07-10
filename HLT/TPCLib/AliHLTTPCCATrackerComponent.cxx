

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  Ivan Kisel <kisel@kip.uni-heidelberg.de>              *
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
// a TPC tracker processing component for the HLT based on CA by Ivan Kisel  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCCATrackerComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCCATracker.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCVertexData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCDefinitions.h"
//#include "AliHLTTPC.h"
#include <stdlib.h>
#include <iostream>
#include <errno.h>

//static bool ask = true;
static bool ask = false;

// this is a global object used for automatic component registration, do not use this
AliHLTTPCCATrackerComponent gAliHLTTPCCATrackerComponent;

ClassImp(AliHLTTPCCATrackerComponent)

AliHLTTPCCATrackerComponent::AliHLTTPCCATrackerComponent()
  :
  fTracker(NULL),
  fVertex(NULL),
  fBField(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCCATrackerComponent::AliHLTTPCCATrackerComponent(const AliHLTTPCCATrackerComponent&)
  :
  fTracker(NULL),
  fVertex(NULL),
  fBField(0)
{
  // see header file for class documentation
  HLTFatal("copy constructor untested");
}

AliHLTTPCCATrackerComponent& AliHLTTPCCATrackerComponent::operator=(const AliHLTTPCCATrackerComponent&)
{
  // see header file for class documentation
  HLTFatal("assignment operator untested");
  return *this;
}

AliHLTTPCCATrackerComponent::~AliHLTTPCCATrackerComponent()
    {
  // see header file for class documentation
    }

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCCATrackerComponent::GetComponentID()
    {
  // see header file for class documentation
    return "TPCCATracker";
    }

void AliHLTTPCCATrackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
    {
  // see header file for class documentation
    list.clear();
    list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
    list.push_back( AliHLTTPCDefinitions::fgkVertexDataType );
    }

AliHLTComponentDataType AliHLTTPCCATrackerComponent::GetOutputDataType()
    {
  // see header file for class documentation
    return AliHLTTPCDefinitions::fgkTrackSegmentsDataType;
    }

void AliHLTTPCCATrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
  // see header file for class documentation
    // XXX TODO: Find more realistic values.
    constBase = 0;
    inputMultiplier = 0.2;
    }

AliHLTComponent* AliHLTTPCCATrackerComponent::Spawn()
    {
  // see header file for class documentation
    return new AliHLTTPCCATrackerComponent;
    }

int AliHLTTPCCATrackerComponent::DoInit( int argc, const char** argv )
    {
  // see header file for class documentation
    Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoInit", "DoInit", "DoInit()" );

    if ( fTracker || fVertex )
	return EINPROGRESS;

    fTracker = new AliHLTTPCCATracker();
    fTracker->CACreateHistos();

    fVertex = new AliHLTTPCVertex();


/* ---------------------------------------------------------------------------------
 * cmdline arguments not needed so far

    int i = 0;
    char* cpErr;

    while ( i < argc )
	{
	if ( !strcmp( argv[i], "bfield" ) )
	    {
	    if ( argc <= i+1 )
		{
		Logging( kHLTLogError, "HLT::TPCSliceTracker::DoInit", "Missing B-field", "Missing B-field specifier." );
		return ENOTSUP;
		}
	    fBField = strtod( argv[i+1], &cpErr );
	    if ( *cpErr )
		{
		Logging( kHLTLogError, "HLT::TPCSliceTracker::DoInit", "Missing multiplicity", "Cannot convert B-field specifier '%s'.", argv[i+1] );
		return EINVAL;
		}
	    i += 2;
	    continue;
	    }

	Logging(kHLTLogError, "HLT::TPCSliceTracker::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
	return EINVAL;
	}
--------------------------------------------------------------------------------- */

    return 0;
    }

int AliHLTTPCCATrackerComponent::DoDeinit()
    {
  // see header file for class documentation
    if ( fTracker )
	delete fTracker;
    fTracker = NULL;
    if ( fVertex )
	delete fVertex;
    fVertex = NULL;
    return 0;
    }

int AliHLTTPCCATrackerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
    {
  // see header file for class documentation
    Logging( kHLTLogDebug, "HLT::TPCCATracker::DoEvent", "DoEvent", "DoEvent()" );
    if ( evtData.fBlockCnt<=0 )
      {
	Logging( kHLTLogWarning, "HLT::TPCCATracker::DoEvent", "DoEvent", "no blocks in event" );
	return 0;
      }

    const AliHLTComponentBlockData* iter = NULL;
    unsigned long ndx;
    AliHLTTPCClusterData* inPtrSP;
    AliHLTTPCVertexData* inPtrV = NULL;
    const AliHLTComponentBlockData* vertexIter=NULL;
    AliHLTTPCTrackletData* outPtr;
    AliHLTUInt8_t* outBPtr;

    AliHLTUInt32_t vSize = 0;
    UInt_t offset=0, mysize, tSize = 0;
    outBPtr = outputPtr;

    outPtr = (AliHLTTPCTrackletData*)(outBPtr);
    fTracker->SetOutPtr( outPtr->fTracklets );

    // ------------------------------------------

    Int_t slice=-1, patch=-1, row[2];
    Int_t minPatch=INT_MAX, maxPatch = 0;
    offset = 0;
    std::vector<Int_t> slices;
    std::vector<Int_t>::iterator slIter, slEnd;
    std::vector<unsigned> sliceCnts;
    std::vector<unsigned>::iterator slCntIter;
    Int_t vertexSlice=-1;

    // Find min/max rows used in total and find and read out vertex if it is present
    // also determine correct slice number, if multiple slice numbers are present in event
    // (which should not happen in the first place) we use the one that occurs the most times
    row[0] = 0;
    row[1] = 0;
    bool found;

    for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
	{
	iter = blocks+ndx;


	slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
	found = false;
	slIter = slices.begin();
	slEnd = slices.end();
	slCntIter = sliceCnts.begin();
	while ( slIter != slEnd )
	    {
	    if ( *slIter == slice )
		{
		found = true;
		break;
		}
	    slIter++;
	    slCntIter++;
	    }
	if ( !found )
	    {
	    slices.insert( slices.end(), slice );
	    sliceCnts.insert( sliceCnts.end(), 1 );
	    }
	else
	    *slCntIter++;

	if ( iter->fDataType == AliHLTTPCDefinitions::fgkVertexDataType )
	    {
	    inPtrV = (AliHLTTPCVertexData*)(iter->fPtr);
	    vertexIter = iter;
	    vSize = iter->fSize;
	    fVertex->Read( inPtrV );
	    vertexSlice = slice;
	    }
	if ( iter->fDataType == AliHLTTPCDefinitions::fgkClustersDataType )
	    {
	    patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
	    if ( minPatch>patch )
		{
		minPatch = patch;
		row[0] = AliHLTTPCTransform::GetFirstRow( patch );
		}
	    if ( maxPatch<patch )
		{
		maxPatch = patch;
		row[1] = AliHLTTPCTransform::GetLastRow( patch );
		}
	    }
	}

    // Determine slice number to really use.
    if ( slices.size()>1 )
	{
	Logging( kHLTLogError, "HLT::TPCSliceTracker::DoEvent", "Multiple slices found in event",
		 "Multiple slice numbers found in event 0x%08lX (%lu). Determining maximum occuring slice number...",
		 evtData.fEventID, evtData.fEventID );
	unsigned maxCntSlice=0;
	slIter = slices.begin();
	slEnd = slices.end();
	slCntIter = sliceCnts.begin();
	while ( slIter != slEnd )
	    {
	    Logging( kHLTLogError, "HLT::TPCSliceTracker::DoEvent", "Multiple slices found in event",
		     "Slice %lu found %lu times.", *slIter, *slCntIter );
	    if ( maxCntSlice<*slCntIter )
		{
		maxCntSlice = *slCntIter;
		slice = *slIter;
		}
	    slIter++;
	    slCntIter++;
	    }
	Logging( kHLTLogError, "HLT::TPCSliceTracker::DoEvent", "Multiple slices found in event",
		 "Using slice %lu.", slice );
	}
    else if ( slices.size()>0 )
      {
	slice = *(slices.begin());
      }
    else
      {
	slice = -1;
      }
    
    if ( vertexSlice != slice )
	{
	// multiple vertex blocks in event and we used the wrong one...
	found = false;
	for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
	    {
	    iter = blocks+ndx;
	    if ( iter->fDataType == AliHLTTPCDefinitions::fgkVertexDataType && slice==AliHLTTPCDefinitions::GetMinSliceNr( *iter ) )
		{
		inPtrV = (AliHLTTPCVertexData*)(iter->fPtr);
		vertexIter = iter;
		vSize = iter->fSize;
		fVertex->Read( inPtrV );
		break;
		}
	    }
	}

    //    fTracker->InitSector( slice, row, fEta );
    //    fTracker->SetVertex(fVertex);

    mysize = 0;
    // read in all hits
    std::vector<unsigned long> patchIndices;
    std::vector<unsigned long>::iterator pIter, pEnd;
    for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
	{
	iter = blocks+ndx;

	if ( iter->fDataType == AliHLTTPCDefinitions::fgkClustersDataType && slice==AliHLTTPCDefinitions::GetMinSliceNr( *iter ) )
	    {
	    patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
	    pIter = patchIndices.begin();
	    pEnd = patchIndices.end();
	    while ( pIter!=pEnd && AliHLTTPCDefinitions::GetMinSliceNr( blocks[*pIter] ) < patch )
		pIter++;
	    patchIndices.insert( pIter, ndx );
	    }
	}

    fTracker->CAInitialize();

    pIter = patchIndices.begin();
    pEnd = patchIndices.end();
    while ( pIter!=pEnd )
	{
	ndx = *pIter;
	iter = blocks+ndx;

	patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
	inPtrSP = (AliHLTTPCClusterData*)(iter->fPtr);
	    
	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Reading hits",
		 "Reading hits for slice %d - patch %d", slice, patch );

	fTracker->CAReadPatchHits( patch, inPtrSP->fSpacePointCnt, inPtrSP->fSpacePoints );
	fTracker->CAFindPatchTracks( patch ); 

	pIter++;
	}

    fTracker->CAFindSliceTracks();

    //#ifdef XXX
    char symbol;
    if (ask){
      do{
	std::cin.get(symbol);
	if (symbol == 'r')
	  ask = false;
      } while (symbol != '\n');
    }
    //#endif //XXX


    UInt_t ntracks0=0;

    mysize = fTracker->GetOutputSize();
    ntracks0 = fTracker->GetOutputNTracks();
    outPtr->fTrackletCnt = fTracker->GetOutputNTracks();
    
    Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracks",
	     "Input: Number of tracks: %lu Slice/MinPatch/MaxPatch/RowMin/RowMax: %lu/%lu/%lu/%lu/%lu.", 
	     ntracks0, slice, minPatch, maxPatch, row[0], row[1] );

    tSize += mysize+sizeof(AliHLTTPCTrackletData);
    outBPtr += mysize+sizeof(AliHLTTPCTrackletData);
    
    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset = offset;
    bd.fSize = tSize;
    bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( slice, slice, minPatch, maxPatch );
    outputBlocks.push_back( bd );

#ifdef FORWARD_VERTEX_BLOCK
    if ( vertexIter )
	{
	// Copy the descriptor block for the vertex information.
	bd = *vertexIter;
	outputBlocks.push_back( bd );
	}
#endif // FORWARD_VERTEX_BLOCK

    size = tSize;

    fTracker->CAWriteHistos();

    return 0;
    }

	
