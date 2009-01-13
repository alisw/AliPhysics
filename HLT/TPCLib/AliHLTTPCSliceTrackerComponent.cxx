// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCSliceTrackerComponent.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  The TPC conformal mapping tracker component.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include <climits>
#include "AliHLTTPCSliceTrackerComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCConfMapper.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCVertexData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCInterMerger.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTCommonCDBEntries.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
//#include "AliHLTTPC.h"
//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCSliceTrackerComponent)

AliHLTTPCSliceTrackerComponent::AliHLTTPCSliceTrackerComponent()
  :
  fTracker(NULL),
  fVertex(NULL),
  fDoNonVertex(false),
  fDoPP(false),
  fDoPbPb(false),
  fMultiplicity(4000),
  fBField(0.4),
  fnonvertextracking(kFALSE),
  fmainvertextracking(kTRUE),
  fPhisegment(50),
  fEtasegment(100),
  fTrackletlength(3),
  fTracklength(60),
  fRowscopetracklet(6),
  fRowscopetrack(6),
  fMinPtFit(0),
  fMaxangle(0.1745),
  fGoodDist(5),
  fHitChi2Cut(100),
  fGoodHitChi2(5),
  fTrackChi2Cut(50),
  fMaxdist(50),
  fMaxphi(0.1),
  fMaxeta(0.1),
  fpInterMerger(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  fEta[0] = 0.;
  fEta[1] = 1.1;
}

AliHLTTPCSliceTrackerComponent::~AliHLTTPCSliceTrackerComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCSliceTrackerComponent::GetComponentID()
{
  // see header file for class documentation

  return "TPCSliceTracker";
}

void AliHLTTPCSliceTrackerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
  list.push_back( AliHLTTPCDefinitions::fgkVertexDataType );
}

AliHLTComponentDataType AliHLTTPCSliceTrackerComponent::GetOutputDataType()
{
  // see header file for class documentation
  return AliHLTTPCDefinitions::fgkTrackSegmentsDataType;
}

void AliHLTTPCSliceTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 0;
  inputMultiplier = 1;
}

AliHLTComponent* AliHLTTPCSliceTrackerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCSliceTrackerComponent;
}

void AliHLTTPCSliceTrackerComponent::SetTrackerParam(Int_t phiSegments, Int_t etaSegments,
				   Int_t trackletlength, Int_t tracklength,
				   Int_t rowscopetracklet, Int_t rowscopetrack,
				   Double_t minPtFit, Double_t maxangle,
				   Double_t goodDist, Double_t hitChi2Cut,
				   Double_t goodHitChi2, Double_t trackChi2Cut,
				   Int_t maxdist, Double_t maxphi,Double_t maxeta)
{
  // see header file for class documentation
    //fTracker->SetClusterFinderParam( fXYClusterError, fZClusterError, kTRUE ); // ??
    //Set parameters input to the tracker
    //If no arguments are given, default parameters will be used
    
    fTracker->SetNSegments(phiSegments,etaSegments);
    fTracker->SetMaxDca(minPtFit);
    //   fTracker->MainVertexSettings(trackletlength,tracklength,rowscopetracklet,rowscopetrack);

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


void AliHLTTPCSliceTrackerComponent::SetTrackerParam( Bool_t doPP, Bool_t doPbPb, Int_t multiplicity, 
						      Double_t bField, Int_t etasegment, Double_t hitchi2cut, 
						      Int_t rowscopetracklet, Int_t rowscopetrack, 
						      Int_t trackletlength, Int_t tracklength )
{
  // see header file for class documentation
  AliHLTTPCTransform::SetBField( bField );
  Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoInit", "BField", "Setting b field to %f T\n", bField );
  
  if ( doPP )
    {
      //tracker->SetClusterFinderParam(xyerror,zerror,kTRUE); // ??
      
      SetTrackerParam( 50, 100, 3, 60,
		       6, 6,
		       0, 0.1745, 5, 100,
		       5, 50, 50, 0.1, 0.1);
    }
  else if(doPbPb)
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
			       5, 50, 50, 0.1, 0.1);
	      break;
	    case 1: // 0.4
	      SetTrackerParam( 50, 100, 3, 10,
			       2, 4,
			       0, 0.1745, 5, 100,
			       5, 50, 50, 0.1, 0.1);
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
			       5, 20, 50, 0.1, 0.1);
	      break;
	    case 1: // 0.4
	      SetTrackerParam( 50, 100, 3, 10,
			       2, 5,
			       0, 0.1745, 5, 30,
			       5, 20, 50, 0.1, 0.1);
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
			       5, 10 , 50, 0.1, 0.1);
	      break;
	    case 1: // 0.4
	      SetTrackerParam( 50, 100, 3, 10,
			       2, 10,
			       0, 0.1745, 5, 20,
			       5, 10, 50, 0.1, 0.1);
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
			       5, 5, 50, 0.1, 0.1);
	      break;
	    case 1: // 0.4
	      SetTrackerParam( 50, 100, 3, 10,
			       2, 15,
			       0, 0.1745, 5, 15,
			       5, 5, 50, 0.1, 0.1);
	      break;
	    }
	  break;
	}
      //	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoInit", "BField", "Setting b field to %f\n", bfs[closestBf] );
      //	AliHLTTPCTransform::SetBField( bfs[closestBf] );
      //	AliHLTTPCTransform::SetBField( bField );
      //	Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoInit", "BField", "Setting b field to %f\n", bField );
    }
  else
    {
      SetTrackerParam( 50, etasegment, trackletlength, tracklength,
		       rowscopetracklet, rowscopetrack,
		       0, 0.1745, 5, hitchi2cut,
		       5, 50, 50, 0.1, 0.1);
    }
}
	
int AliHLTTPCSliceTrackerComponent::DoInit( int argc, const char** argv )
    {
  // see header file for class documentation
    Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoInit", "DoInit", "DoInit()" );

    if ( fTracker || fVertex )
	return EINPROGRESS;
    fTracker = new AliHLTTPCConfMapper();
    fVertex = new AliHLTTPCVertex();
    fEta[0] = 0.;
    fEta[1] = 1.1;
    fDoNonVertex = false;
    fMultiplicity = 4000;
    fBField = 0.5;
    fDoPP = false;
    fDoPbPb = false;
    fPhisegment=50;
    fEtasegment=100;
    fTrackletlength=3;
    fTracklength=60;
    fRowscopetracklet=6;
    fRowscopetrack=6;
    fMinPtFit=0;
    fMaxangle=0.1745;
    fGoodDist=5;
    fHitChi2Cut=100;
    fGoodHitChi2=5;
    fTrackChi2Cut=50;
    fMaxdist=50;
    fMaxphi=0.1;
    fMaxeta=0.1;
    int iResult=0;
        
    TString configuration="";
    TString argument="";
    for (int i=0; i<argc && iResult>=0; i++) {
      argument=argv[i];
      if (!configuration.IsNull()) configuration+=" ";
      configuration+=argument;
    }

    if (!configuration.IsNull()) {
      iResult=Configure(configuration.Data());
    } else {
      iResult=Reconfigure(NULL, NULL);
    }
    
    return iResult;
    }

int AliHLTTPCSliceTrackerComponent::DoDeinit()
{
  // see header file for class documentation
  if ( fTracker )
    delete fTracker;
  fTracker = NULL;
  if ( fVertex )
    delete fVertex;
  fVertex = NULL;
  if (fpInterMerger) {
    delete fpInterMerger;
  }
  fpInterMerger=NULL;
  return 0;
}

int AliHLTTPCSliceTrackerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, AliHLTComponentBlockDataList& outputBlocks )
{
  // see header file for class documentation
  int iResult=0;
  AliHLTUInt32_t capacity=size;
  size=0;

  if (!IsDataEvent()) return 0;

    if ( evtData.fBlockCnt<=0 )
      {
	Logging( kHLTLogWarning, "HLT::TPCSliceTracker::DoEvent", "DoEvent", "no blocks in event" );
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

	bool bIsClusterDataBlock=false;
	bool bIsVertexDataBlock=false;
	if(!(bIsClusterDataBlock=(iter->fDataType==AliHLTTPCDefinitions::fgkClustersDataType)) &&
	   !(bIsVertexDataBlock=(iter->fDataType==AliHLTTPCDefinitions::fgkVertexDataType))){
	  continue;
	}

	slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
	if (slice<0 || slice>=AliHLTTPCTransform::GetNSlice()) {
	  HLTError("invalid slice number %d extracted from specification 0x%08lx,  skipping block of type %s",
		   slice, iter->fSpecification, DataType2Text(iter->fDataType).c_str());
	  // just remember the error, if there are other valid blocks ignore the
	  // error, return code otherwise
	  iResult=-EBADF;
	  continue;
	}
	if (slice!=AliHLTTPCDefinitions::GetMaxSliceNr( *iter )) {
	  // the code was not written for/ never used with multiple slices
	  // in one data block/ specification
	  HLTWarning("specification 0x%08lx indicates multiple slices in data block %s: never used before, please audit the code",
		     iter->fSpecification, DataType2Text(iter->fDataType).c_str());
	}
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

	if (bIsVertexDataBlock)
	    {
	    inPtrV = (AliHLTTPCVertexData*)(iter->fPtr);
	    vertexIter = iter;
	    vSize = iter->fSize;
	    fVertex->Read( inPtrV );
	    vertexSlice = slice;
	    }
	if (bIsClusterDataBlock)
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
	// there is no sense in running the tracker without input, do not send an
	// empty output block
	return iResult;
      }
    
    iResult=0;

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

    fTracker->InitSector( slice, row, fEta );
    fTracker->SetVertex(fVertex);
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
	//fTracker->ReadHits( inPtrSP->fSpacePointCnt, inPtrSP->fSpacePoints );
	fTracker->ReadHitsChecked(inPtrSP->fSpacePointCnt, inPtrSP->fSpacePoints,iter->fSize );
	pIter++;
	}


    if ( fmainvertextracking == kTRUE && fnonvertextracking == kFALSE){	
      Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracking", " ---MAINVERTEXTRACKING---");
      fTracker->MainVertexTrackingA();
      fTracker->MainVertexTrackingB();
      fTracker->FillTracks();
    }
    else if ( fmainvertextracking == kTRUE && fnonvertextracking == kTRUE){	
      Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracking", " ---MAINVERTEXTRACKING---");
      fTracker->MainVertexTrackingA();
      fTracker->MainVertexTrackingB();
      fTracker->FillTracks();	
      Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracking", " ---NONVERTEXTRACKING---");
      fTracker->NonVertexTracking();
    }
    else if ( fmainvertextracking == kFALSE && fnonvertextracking == kTRUE){
      Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracking", " ---NONVERTEXTRACKING---");
      fTracker->NonVertexTracking();	
      fTracker->FillTracks();	
    }

    UInt_t ntracks0 =0;
    if(fpInterMerger){
      AliHLTTPCMemHandler memory;
      AliHLTTPCTrackSegmentData *trackdata0  = 
	(AliHLTTPCTrackSegmentData *) memory.Allocate(fTracker->GetTracks());
      memory.TrackArray2Memory(ntracks0,trackdata0,fTracker->GetTracks());
      fpInterMerger->Reset();
      fpInterMerger->Init(row,patch);
      fpInterMerger->FillTracks(ntracks0,trackdata0);
      fpInterMerger->Merge();
    } 
    ntracks0=0;
    AliHLTTPCTrackArray* pArray=fTracker->GetTracks();
    if (pArray->GetOutSize()+sizeof(AliHLTTPCTrackletData)<=capacity) {
    outPtr = (AliHLTTPCTrackletData*)(outBPtr);
    mysize = pArray->WriteTracks( ntracks0, outPtr->fTracklets );
    outPtr->fTrackletCnt = ntracks0;

    Logging( kHLTLogDebug, "HLT::TPCSliceTracker::DoEvent", "Tracks",
	     "Input: Number of tracks: %lu Slice/MinPatch/MaxPatch/RowMin/RowMax: %lu/%lu/%lu/%lu/%lu.", 
	     ntracks0, slice, minPatch, maxPatch, row[0], row[1] );

    fTracker->Reset();

    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset = offset;
    bd.fSize = mysize+sizeof(AliHLTTPCTrackletData);
    bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( slice, slice, minPatch, maxPatch );
    outputBlocks.push_back( bd );

    tSize += bd.fSize;
    outBPtr += bd.fSize;
    } else {
      iResult=-ENOSPC;
    }

#ifdef FORWARD_VERTEX_BLOCK
    if ( vertexIter )
	{
	// Copy the descriptor block for the vertex information.
	bd = *vertexIter;
	outputBlocks.push_back( bd );
	}
#endif // FORWARD_VERTEX_BLOCK

    size = tSize;
    return iResult;
    }

int AliHLTTPCSliceTrackerComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;
  Bool_t bDoMerger=kTRUE;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
      
      if (argument.CompareTo("-disable-merger")==0) {
	HLTInfo("Disabled Inter Merger");
	bDoMerger = kFALSE;
	continue;
      }
      else if (argument.CompareTo("-pp-run")==0) {
	HLTInfo("Using Trackparameters for pp-run");
	fDoPP = true;
	continue;
      }
      else if (argument.CompareTo("-PbPb-run")==0) {
	HLTInfo("Using Trackparameters for Pb-Pb-run");
	fDoPbPb = true;
	continue;
      }
      else if (argument.CompareTo("-nonvertextracking")==0) {
	HLTInfo("Doing Nonvertex Tracking");
	fnonvertextracking = kTRUE;
	continue;
      }     
      else if (argument.CompareTo("-mainvertextrackingoff")==0) {
	HLTInfo("Mainvertex Tracking off");
	fmainvertextracking = kFALSE;
	continue;
      }
      else if (argument.CompareTo("-multiplicity")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Multiplicity set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fMultiplicity=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      } 
      else if (argument.CompareTo("-solenoidBz")==0 || argument.CompareTo("-bfield")==0) {
	if(argument.CompareTo("-bfield")==0){
	  HLTWarning("-bfield is the old way. The field is set, but please use -solenoidBz.");
	}
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Magnetic Field set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fBField=((TObjString*)pTokens->At(i))->GetString().Atof();
	continue;
      } 
      else if (argument.CompareTo("-etarange")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Etarange set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fEta[1]=((TObjString*)pTokens->At(i))->GetString().Atof();
	continue;
      }
      else if (argument.CompareTo("-etasegment")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Number of Etasegment: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fEtasegment=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      }
      else if (argument.CompareTo("-chi2cut")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("chi2cut set to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fHitChi2Cut=((TObjString*)pTokens->At(i))->GetString().Atof();
	continue;
      }
      else if (argument.CompareTo("-rowscopetracklet")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Number of row to look for next cluster for tracklet: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fRowscopetracklet=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      }
      else if (argument.CompareTo("-rowscopetrack")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Number of row to look for next cluster for track: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fRowscopetrack=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      }
      else if (argument.CompareTo("-trackletlength")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Minimum number of clusters on a Tracklet: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fTrackletlength=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      }
      else if (argument.CompareTo("-tracklength")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Minimum number of clusters on a Track: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fTracklength=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      }
      else if (argument.CompareTo("-clusterZ")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Minimum number of clusters on a Track: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fTracker->SetClusterCutZ(((TObjString*)pTokens->At(i))->GetString().Atoi());
	continue;
      }
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  if (fBField == 0.){
    // parameter for B=0 T 
    fDoPP = kTRUE;
    fnonvertextracking = kTRUE;
    fmainvertextracking = kFALSE;
  }
 
  if (bDoMerger)
    fpInterMerger = new AliHLTTPCInterMerger();
  else
    fpInterMerger = NULL;

  SetTrackerParam(fDoPP,fDoPbPb,fMultiplicity,fBField,fEtasegment,fHitChi2Cut,fRowscopetracklet,fRowscopetrack,fTrackletlength,fTracklength);
  
  return iResult;
}

int AliHLTTPCSliceTrackerComponent::ReadPreprocessorValues(const char* modules)
{
  // see header file for class documentation
  
  int iResult = 0;
  TString str(modules);
  if(str.Contains("HLT") || str.Contains("TPC") || str.Contains("GRP")){
  
    const char* pathBField=kAliHLTCDBSolenoidBz;
    if (pathBField) {

      HLTInfo("reconfigure B-Field from entry %s, modules %s", pathBField,(modules!=NULL && modules[0]!=0)?modules:"<none>");
      //AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(pathBField/*,GetRunNo()*/);
      
      AliCDBPath path(pathBField);
      
      AliCDBStorage *stor = AliCDBManager::Instance()->GetDefaultStorage();
      Int_t version    = stor->GetLatestVersion(pathBField, GetRunNo());
      Int_t subVersion = stor->GetLatestSubVersion(pathBField, GetRunNo(), version);
      AliCDBEntry *pEntry = stor->Get(path,GetRunNo(), version, subVersion);
      
      HLTImportant("RunNo %d, Version %d, subversion %d", GetRunNo(), version, subVersion);
      
      if (pEntry) {
    	TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
    	if (pString) {
   	  HLTImportant("received configuration object string: \'%s\'", pString->GetString().Data());
   	  iResult=Configure(pString->GetString().Data());
    	} else {
   	  HLTError("configuration object \"%s\" has wrong type, required TObjString", pathBField);
    	}
      } else {
    	HLTError("cannot fetch object \"%s\" from CDB", pathBField);
      }
    }
  }  
  return iResult;
}

int AliHLTTPCSliceTrackerComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation

  int iResult=0;
  const char* path="HLT/ConfigTPC/SliceTrackerComponent";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("cannot fetch object \"%s\" from CDB", path);
    }
  }

  const char* pathBField=kAliHLTCDBSolenoidBz;

  if (pathBField) {
    HLTInfo("reconfigure B-Field from entry %s, chain id %s", path,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(pathBField/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("cannot fetch object \"%s\" from CDB", path);
    }
  }
  
  return iResult;

}

void AliHLTTPCSliceTrackerComponent::SetTrackerParam1()
{
  // see header file for class documentation
  SetTrackerParam( 10, 20, 5, 10, 2,2,
		   0, 1.31, 5, 100,
		   50, 100, 50, 0.1, 0.1);
}

	
