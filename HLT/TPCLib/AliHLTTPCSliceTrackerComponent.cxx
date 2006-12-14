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
// a TPC sector tracker processing component for the HLT                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCSliceTrackerComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCConfMapper.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCVertexData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
//#include "AliHLTTPC.h"
#include <stdlib.h>
#include <errno.h>

// this is a global object used for automatic component registration, do not use this
AliHLTTPCSliceTrackerComponent gAliHLTTPCSliceTrackerComponent;

ClassImp(AliHLTTPCSliceTrackerComponent)

AliHLTTPCSliceTrackerComponent::AliHLTTPCSliceTrackerComponent()
    {
    fTracker = NULL;
    fVertex = NULL;
    fEta[0] = 0.;
    fEta[1] = 1.1;
    fDoNonVertex = false;
    fMultiplicity = 4000;
    fBField = 0.4;
    fDoPP = false;
// BEGINN ############################################## MODIFIY JMT
    fnonvertextracking = kFALSE;   // enable NONVERTEX Tracking
    fmainvertextracking = kTRUE;   // enable MAINVERTEX Tracking
// END ################################################# MODIFIY JMT
    }

AliHLTTPCSliceTrackerComponent::~AliHLTTPCSliceTrackerComponent()
    {
    }

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCSliceTrackerComponent::GetComponentID()
    {
    return "TPCSliceTracker";
    }

void AliHLTTPCSliceTrackerComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
    {
    list.clear();
    list.push_back( AliHLTTPCDefinitions::gkClustersDataType );
    list.push_back( AliHLTTPCDefinitions::gkVertexDataType );
    }

AliHLTComponent_DataType AliHLTTPCSliceTrackerComponent::GetOutputDataType()
    {
    return AliHLTTPCDefinitions::gkTrackSegmentsDataType;
    }

void AliHLTTPCSliceTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
    // XXX TODO: Find more realistic values.
    constBase = 0;
    inputMultiplier = 0.2;
    }

AliHLTComponent* AliHLTTPCSliceTrackerComponent::Spawn()
    {
    return new AliHLTTPCSliceTrackerComponent;
    }

void AliHLTTPCSliceTrackerComponent::SetTrackerParam(Int_t phi_segments, Int_t eta_segments,
				   Int_t trackletlength, Int_t tracklength,
				   Int_t rowscopetracklet, Int_t rowscopetrack,
				   Double_t min_pt_fit, Double_t maxangle,
				   Double_t goodDist, Double_t hitChi2Cut,
				   Double_t goodHitChi2, Double_t trackChi2Cut,
				   Int_t maxdist, Double_t maxphi,Double_t maxeta, bool vertexConstraints )
    {
    //fTracker->SetClusterFinderParam( fXYClusterError, fZClusterError, kTRUE ); // ??
    //Set parameters input to the tracker
    //If no arguments are given, default parameters will be used
    
    fTracker->SetNSegments(phi_segments,eta_segments);
    fTracker->SetMaxDca(min_pt_fit);
    //   fTracker->MainVertexSettings(trackletlength,tracklength,rowscopetracklet,rowscopetrack);

// BEGINN ############################################## MODIFIY JMT
#if 1
    Logging( kHLTLogDebug, "HLT::TPCSliceTracker::SetTrackerParam", "Tracking", "==============================" );

    if ( fmainvertextracking == kTRUE && fnonvertextracking == kFALSE){
	fTracker->SetTrackCuts(hitChi2Cut,goodHitChi2,trackChi2Cut,maxdist,kTRUE);
	fTracker->SetTrackletCuts(maxangle,goodDist,kTRUE);

	fTracker->MainVertexSettings( trackletlength, tracklength, rowscopetracklet, rowscopetrack, maxphi, maxeta);	
	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::SetTrackerParam", "Tracking", "MAINVERTEXTRACKING" );
    }
    else if ( fmainvertextracking == kTRUE && fnonvertextracking == kTRUE){
	fTracker->SetTrackCuts(hitChi2Cut,goodHitChi2,trackChi2Cut,maxdist,kTRUE);
	fTracker->SetTrackCuts(hitChi2Cut,goodHitChi2,trackChi2Cut,maxdist,kFALSE);
	fTracker->SetTrackletCuts(maxangle,goodDist,kTRUE);
	fTracker->SetTrackletCuts(maxangle,goodDist,kFALSE);

	fTracker->MainVertexSettings( trackletlength, tracklength, rowscopetracklet, rowscopetrack, maxphi, maxeta);
	fTracker->NonVertexSettings( trackletlength, tracklength, rowscopetracklet, rowscopetrack);	
	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::SetTrackerParam", "Tracking", "MAINVERTEXTRACKING - NONVERTEXTRACKING" );
    }
    else if ( fmainvertextracking == kFALSE && fnonvertextracking == kTRUE){
	fTracker->SetTrackCuts(hitChi2Cut,goodHitChi2,trackChi2Cut,maxdist,kFALSE);
	fTracker->SetTrackletCuts(maxangle,goodDist,kFALSE);

	fTracker->NonVertexSettings( trackletlength, tracklength, rowscopetracklet, rowscopetrack);
	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::SetTrackerParam", "Tracking", "NONVERTEXTRACKING" );
    }
#else
    fTracker->SetTrackCuts(hitChi2Cut,goodHitChi2,trackChi2Cut,maxdist,vertexConstraints);
    fTracker->SetTrackletCuts(maxangle,goodDist,vertexConstraints);
    
    if( vertexConstraints )
	fTracker->MainVertexSettings( trackletlength, tracklength, rowscopetracklet, rowscopetrack, maxphi, maxeta);
    else
	fTracker->NonVertexSettings( trackletlength, tracklength, rowscopetracklet, rowscopetrack);
#endif
// END ################################################# MODIFIY JMT

    //fTracker->SetParamDone(true);
    /* Matthias 13.12.2006
     * the global variable AliHLTTPCS::fgDoVertexFit has never been used so far
     * and has always been kTRUE.
     * In order to remove the AliHLTTPC class (which is the old steering class for
     * HLT (TPC) tracking) from the compilation, this function can not be activated
     * again. We have to think about a more elegant way to specify the parameters
     * anyway. The following line was surely for some testing, but was never active
     * in a tested release.
     */
    //AliHLTTPC::SetVertexFit( kFALSE );
    
    fTracker->InitVolumes();
    }

void AliHLTTPCSliceTrackerComponent::SetTrackerParam( bool doPP, int multiplicity, double bField )
    {
    AliHLTTPCTransform::SetBField( bField );
    Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoInit", "BField", "Setting b field to %f T\n", bField );

    if ( doPP )
	{
	//tracker->SetClusterFinderParam(xyerror,zerror,kTRUE); // ??
	SetTrackerParam( 50, 100, 3, 10,
			 2, 2,
			 0, 0.1745, 5, 100,
			 5, 50, 50, 0.1, 0.1, kTRUE);
	}
    else
	{
	int mults[] = { 1000, 2000, 4000, 8000 };
	int multCount = 4;
	int closestMult = 0;
	int i;
	int multDist, tmpMultDist;
	if ( multiplicity>mults[closestMult] )
	    multDist = multiplicity-mults[closestMult];
	else
	    multDist = mults[closestMult]-multiplicity;
	for ( i = 1; i < multCount; i++ )
	    {
	    if ( multiplicity>mults[i] )
		tmpMultDist = multiplicity-mults[i];
	    else
		tmpMultDist = mults[i]-multiplicity;
	    if ( tmpMultDist < multDist )
		{
		closestMult = i;
		multDist = tmpMultDist;
		}
	    }
	
	double bfs[] = { 0.2, 0.4 };
	int bfCount = 2;
	int closestBf = 0;
	double bfDist, tmpBFDist;
	if ( bField>bfs[closestBf] )
	    bfDist = bField-bfs[closestBf];
	else
	    bfDist = bfs[closestBf]-bField;
	for ( i = 1; i < bfCount; i++ )
	    {
	    if ( bField>bfs[i] )
		tmpBFDist = bField-bfs[i];
	    else
		tmpBFDist = bfs[i]-bField;
	    if ( tmpBFDist < bfDist )
		{
		closestBf = i;
		bfDist = tmpBFDist;
		}
	    }

	switch ( closestMult )
	    {
	    case 0: // 1000
		switch ( closestBf )
		    {
		    case 0: // 0.2
			SetTrackerParam( 50, 100, 3, 10,
					2, 4,
					0, 0.1745, 5, 100,
					5, 50, 50, 0.1, 0.1, kTRUE );
			break;
		    case 1: // 0.4
			SetTrackerParam( 50, 100, 3, 10,
					 2, 4,
					 0, 0.1745, 5, 100,
					 5, 50, 50, 0.1, 0.1, kTRUE );
			break;
		    }
		break;
	    case 1: // 2000
		switch ( closestBf )
		    {
		    case 0: // 0.2
			SetTrackerParam( 50, 100, 3, 10,
					 2, 4,
					 0, 0.1745, 5, 30,
					 5, 20, 50, 0.1, 0.1, kTRUE );
			break;
		    case 1: // 0.4
			SetTrackerParam( 50, 100, 3, 10,
					 2, 5,
					 0, 0.1745, 5, 30,
					 5, 20, 50, 0.1, 0.1, kTRUE );
			break;
		    }
		break;
	    case 2: // 4000
		switch ( closestBf )
		    {
		    case 0: // 0.2
			SetTrackerParam( 50, 100, 3, 10,
					 2 , 10,
					 0, 0.1745, 5, 20,
					 5, 10 , 50, 0.1, 0.1, kTRUE );
			break;
		    case 1: // 0.4
			SetTrackerParam( 50, 100, 3, 10,
					 2, 10,
					 0, 0.1745, 5, 20,
					 5, 10, 50, 0.1, 0.1, kTRUE );
			break;
		    }
		break;
	    case 3: // 8000
		switch ( closestBf )
		    {
		    case 0: // 0.2
			SetTrackerParam( 50, 100, 3, 10,
					 3, 15,
					 0, 0.1745, 5, 10,
					 5, 5, 50, 0.1, 0.1, kTRUE );
			break;
		    case 1: // 0.4
			SetTrackerParam( 50, 100, 3, 10,
					 2, 15,
					 0, 0.1745, 5, 15,
					 5, 5, 50, 0.1, 0.1, kTRUE );
			break;
		    }
		break;
	    }
//	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoInit", "BField", "Setting b field to %f\n", bfs[closestBf] );
//	AliHLTTPCTransform::SetBField( bfs[closestBf] );
//	AliHLTTPCTransform::SetBField( bField );
//	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoInit", "BField", "Setting b field to %f\n", bField );
	}
    }


	
int AliHLTTPCSliceTrackerComponent::DoInit( int argc, const char** argv )
    {
    Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoInit", "DoInit", "DoInit()" );

    if ( fTracker || fVertex )
	return EINPROGRESS;
    fTracker = new AliHLTTPCConfMapper();
    fVertex = new AliHLTTPCVertex();
    fEta[0] = 0.;
    fEta[1] = 1.1;
    fDoNonVertex = false;
    fMultiplicity = 4000;
    fBField = 0.4;
    fDoPP = false;

    int i = 0;
    char* cpErr;
    while ( i < argc )
	{
	if ( !strcmp( argv[i], "pp-run" ) )
	    {
	    fDoPP = true;
	    i++;
	    continue;
	    }
	if ( !strcmp( argv[i], "multiplicity" ) )
	    {
	    if ( argc <= i+1 )
		{
		Logging( kHLTLogError, "HLT::TPCSliceTracker::DoInit", "Missing multiplicity", "Missing event multiplicity specifier." );
		return ENOTSUP;
		}
	    fMultiplicity = strtoul( argv[i+1], &cpErr, 0 );
	    if ( *cpErr )
		{
		Logging( kHLTLogError, "HLT::TPCSliceTracker::DoInit", "Missing multiplicity", "Cannot convert event multiplicity specifier '%s'.", argv[i+1] );
		return EINVAL;
		}
	    i += 2;
	    continue;
	    }
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

// BEGINN ############################################## MODIFIY JMT
	if ( !strcmp( argv[i], "nonvertextracking" ) ){
	    fnonvertextracking = kTRUE;
	    i++;
	    continue;	    
	}
	
	if ( !strcmp( argv[i], "mainvertextrackingoff" ) ){	
	    fmainvertextracking = kFALSE;
	    i++;
	    continue;	    
	}

	if ( !strcmp( argv[i], "etarange" ) ){	
	    if ( argc <= i+1 ){
		Logging( kHLTLogError, "HLT::TPCSliceTracker::DoInit", "Missing Eta range", "Missing Eta-range specifiers." );
		return ENOTSUP;
	    }
	    fEta[1] = strtod( argv[i+1], &cpErr );
	    if ( *cpErr ){
		Logging( kHLTLogError, "HLT::TPCSliceTracker::DoInit", "Missing Eta range", "Cannot convert Eta-range specifier '%s'.", argv[i+1] );
		return EINVAL;
	    }
   
	    i += 2;
	    continue;
	}
// END ################################################# MODIFIY JMT
	Logging(kHLTLogError, "HLT::TPCSliceTracker::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
	return EINVAL;
	}
// #### -B0-CHANGE-START == JMT
    if (fBField == 0.){
	// parameter for B=0 T 
	fDoPP = kTRUE;
	fnonvertextracking = kTRUE;
	fmainvertextracking = kFALSE;
    }
// #### -B0-CHANGE-END == JMT

    SetTrackerParam( fDoPP, fMultiplicity, fBField );
    return 0;
    }

int AliHLTTPCSliceTrackerComponent::DoDeinit()
    {
    if ( fTracker )
	delete fTracker;
    fTracker = NULL;
    if ( fVertex )
	delete fVertex;
    fVertex = NULL;
    return 0;
    }

int AliHLTTPCSliceTrackerComponent::DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
					      AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
    {
    Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "DoEvent", "DoEvent()" );
    if ( evtData.fBlockCnt<=0 )
      {
	Logging( kHLTLogWarning, "HLT::TPCSliceTracker::DoEvent", "DoEvent", "no blocks in event" );
	return 0;
      }
    const AliHLTComponent_BlockData* iter = NULL;
    unsigned long ndx;
    AliHLTTPCClusterData* inPtrSP;
    AliHLTTPCVertexData* inPtrV = NULL;
    const AliHLTComponent_BlockData* vertexIter=NULL;
    AliHLTTPCTrackletData* outPtr;
    AliHLTUInt8_t* outBPtr;
    AliHLTUInt32_t vSize = 0;
    UInt_t offset=0, mysize, tSize = 0;
    outBPtr = outputPtr;
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

	if ( iter->fDataType == AliHLTTPCDefinitions::gkVertexDataType )
	    {
	    inPtrV = (AliHLTTPCVertexData*)(iter->fPtr);
	    vertexIter = iter;
	    vSize = iter->fSize;
	    fVertex->Read( inPtrV );
	    vertexSlice = slice;
	    }
	if ( iter->fDataType == AliHLTTPCDefinitions::gkClustersDataType )
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
	    if ( iter->fDataType == AliHLTTPCDefinitions::gkVertexDataType && slice==AliHLTTPCDefinitions::GetMinSliceNr( *iter ) )
		{
		inPtrV = (AliHLTTPCVertexData*)(iter->fPtr);
		vertexIter = iter;
		vSize = iter->fSize;
		fVertex->Read( inPtrV );
		break;
		}
	    }
	}

    fTracker->InitSector( slice, row, fEta );
    fTracker->SetVertex(fVertex);
    mysize = 0;
    // read in all hits
    std::vector<unsigned long> patchIndices;
    std::vector<unsigned long>::iterator pIter, pEnd;
    for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
	{
	iter = blocks+ndx;

	if ( iter->fDataType == AliHLTTPCDefinitions::gkClustersDataType && slice==AliHLTTPCDefinitions::GetMinSliceNr( *iter ) )
	    {
	    patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
	    pIter = patchIndices.begin();
	    pEnd = patchIndices.end();
	    while ( pIter!=pEnd && AliHLTTPCDefinitions::GetMinSliceNr( blocks[*pIter] ) < patch )
		pIter++;
	    patchIndices.insert( pIter, ndx );
	    }
	}
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
	fTracker->ReadHits( inPtrSP->fSpacePointCnt, inPtrSP->fSpacePoints );
	pIter++;
	}

    outPtr = (AliHLTTPCTrackletData*)(outBPtr);
// BEGINN ############################################## MODIFIY JMT
#if 1
    fTracker->SetPointers();
    if ( fmainvertextracking == kTRUE && fnonvertextracking == kFALSE){	
	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracking", " ---MAINVERTEXTRACKING---");
	fTracker->MainVertexTracking_a();
	fTracker->MainVertexTracking_b();
	fTracker->FillTracks();
    }
    else if ( fmainvertextracking == kTRUE && fnonvertextracking == kTRUE){	
	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracking", " ---MAINVERTEXTRACKING---");
	fTracker->MainVertexTracking_a();
	fTracker->MainVertexTracking_b();
	fTracker->FillTracks();	
	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracking", " ---NONVERTEXTRACKING---");
	fTracker->NonVertexTracking();
    }
    else if ( fmainvertextracking == kFALSE && fnonvertextracking == kTRUE){
	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracking", " ---NONVERTEXTRACKING---");
	fTracker->NonVertexTracking();	
	fTracker->FillTracks();	
    }
#else
    fTracker->MainVertexTracking_a();
    fTracker->MainVertexTracking_b();
    fTracker->FillTracks();
    
    if ( fDoNonVertex )
	fTracker->NonVertexTracking();//Do a second pass for nonvertex tracks
#endif
// END ################################################# MODIFIY JMT
    // XXX Do track merging??
    
    UInt_t ntracks0=0;
    mysize = fTracker->GetTracks()->WriteTracks( ntracks0, outPtr->fTracklets );
    outPtr->fTrackletCnt = ntracks0;
    
    Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracks",
	     "Input: Number of tracks: %lu Slice/MinPatch/MaxPatch/RowMin/RowMax: %lu/%lu/%lu/%lu/%lu.", 
	     ntracks0, slice, minPatch, maxPatch, row[0], row[1] );

    tSize += mysize+sizeof(AliHLTTPCTrackletData);
    outBPtr += mysize+sizeof(AliHLTTPCTrackletData);
    
    AliHLTComponent_BlockData bd;
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
    return 0;
    }

	
