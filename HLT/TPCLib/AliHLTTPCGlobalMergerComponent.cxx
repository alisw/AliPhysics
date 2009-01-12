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

/** @file   AliHLTTPCGlobalMergerComponent.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  HLT TPC global merger component.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include <climits>
#include "AliHLTTPCGlobalMergerComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCGlobalMerger.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCVertexData.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrackletDataFormat.h"
//#include "AliHLTTPCSpacePointData.h"
//#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCDefinitions.h"
#include <cstdlib>
#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCGlobalMergerComponent);

AliHLTTPCGlobalMergerComponent::AliHLTTPCGlobalMergerComponent()
  :
  fGlobalMerger(NULL),
  fVertex(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCGlobalMergerComponent::~AliHLTTPCGlobalMergerComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCGlobalMergerComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCGlobalMerger";
}

void AliHLTTPCGlobalMergerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTPCDefinitions::fgkTrackSegmentsDataType );
  list.push_back( AliHLTTPCDefinitions::fgkVertexDataType );
}

AliHLTComponentDataType AliHLTTPCGlobalMergerComponent::GetOutputDataType()
{
  // see header file for class documentation
  return AliHLTTPCDefinitions::fgkTracksDataType;
}

void AliHLTTPCGlobalMergerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 0;
  inputMultiplier = 1.0;
}

AliHLTComponent* AliHLTTPCGlobalMergerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCGlobalMergerComponent;
}

void AliHLTTPCGlobalMergerComponent::SetMergerParameters(Double_t maxy,Double_t maxz,Double_t maxkappa,Double_t maxpsi,Double_t maxtgl)
{
  // see header file for class documentation
  fGlobalMerger->SetParameter( maxy, maxz, maxkappa, maxpsi, maxtgl );
}

int AliHLTTPCGlobalMergerComponent::DoInit( int /*argc*/, const char** /*argv*/ )
{
  // see header file for class documentation
  if ( fGlobalMerger || fVertex )
    return EINPROGRESS;
  fGlobalMerger = new AliHLTTPCGlobalMerger();
  fVertex = new AliHLTTPCVertex();
  SetMergerParameters();
  return 0;
}

int AliHLTTPCGlobalMergerComponent::DoDeinit()
{
  // see header file for class documentation
  if ( fGlobalMerger )
    delete fGlobalMerger;
  fGlobalMerger = NULL;
  if ( fVertex )
    delete fVertex;
  fVertex = NULL;
  return 0;
}

int AliHLTTPCGlobalMergerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, AliHLTComponentBlockDataList& outputBlocks )
{
  // see header file for class documentation
  int iResult=0;
  AliHLTUInt32_t capacity=size;
  size=0;

  if (!IsDataEvent()) return 0;

  const AliHLTComponentBlockData* iter = NULL;
  const AliHLTComponentBlockData* lastVertexBlock = NULL;
  unsigned long ndx;

  std::vector<SliceData> slices;
  std::vector<SliceData>::iterator sdIter, sdEnd;
  int minSlice = INT_MAX, maxSlice = 0;
  bool found;
  AliHLTTPCTrackletData* inPtr;
  AliHLTTPCTrackletData* outPtr;
  UInt_t tSize = 0;
  Int_t slice=0;

  // Create sorted (by slice number) list of data (tracks and vertex) for each slice present.
  // also note the min and max slice numbers
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      bool bIsTrackSegDataBlock=false;
      bool bIsVertexDataBlock=false;
      if(!(bIsTrackSegDataBlock=(iter->fDataType==AliHLTTPCDefinitions::fgkTrackSegmentsDataType)) &&
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
      found=false;
      sdIter = slices.begin();
      sdEnd = slices.end();
      while ( sdIter != sdEnd )
	{
	  if ( sdIter->fSlice > slice || sdIter->fSlice == slice )
	    break;
	  sdIter++;
	}
      if ( sdIter==sdEnd || sdIter->fSlice>slice )
	{
	  if ( sdIter == sdEnd )
	    maxSlice = slice;
	  if ( sdIter==slices.begin() )
	    minSlice = slice;
	  SliceData sd;
	  sd.fSlice = slice;
	  sd.fVertexBlock = NULL;
	  sd.fVertexBlockIndex = 0;
	  sd.fTrackletBlock = NULL;
	  sd.fTrackletBlockIndex = 0;
	  sdIter = slices.insert( sdIter, sd );
	}
      if ( sdIter->fSlice == slice )
	{
	  if (bIsTrackSegDataBlock)
	    {
	      if ( !sdIter->fTrackletBlock )
		{
		  sdIter->fTrackletBlock = iter;
		  sdIter->fTrackletBlockIndex = ndx;
		}
	      else
		{
		  Logging( kHLTLogError, "HLT::GlobalMerger::DoEvent", "Duplicate track data block",
			   "Duplicate track data block for slice %lu in event 0x%08lX (%lu) - previous block: %lu - new block: %lu.",
			   slice, evtData.fEventID, evtData.fEventID, sdIter->fTrackletBlockIndex, ndx );
		}
	    }
	  if (bIsVertexDataBlock)
	    {
	      lastVertexBlock = iter;
	      if ( !sdIter->fVertexBlock )
		{
		  sdIter->fVertexBlock = iter;
		  sdIter->fVertexBlockIndex = ndx;
		}
	      else
		{
		  Logging( kHLTLogError, "HLT::GlobalMerger::DoEvent", "Duplicate vertex data block",
			   "Duplicate vertex data block for slice %lu in event 0x%08lX (%lu) - previous block: %lu - new block: %lu.",
			   slice, evtData.fEventID, evtData.fEventID, sdIter->fVertexBlockIndex, ndx );
		}
	    }
	}
    }

  // skip the processing if there was no track segment data
  if (slices.size()==0) {
    return iResult;
  }
  iResult=0;

  //fGlobalMerger->Setup( minSlice, maxSlice );
  fGlobalMerger->Setup( 0, 35 );

  if ( !lastVertexBlock )
    {
      Logging( kHLTLogInfo, "HLT::GlobalMerger::DoEvent", "No vertex data block",
	       "No vertex data block found  for event  0x%08lX (%lu).", evtData.fEventID, evtData.fEventID );
      fVertex->SetZero();
    }
  else
    fVertex->Read( (AliHLTTPCVertexData*)( lastVertexBlock->fPtr ) );

  // Add all tracks into the merger
  sdIter = slices.begin();
  sdEnd = slices.end();
  int lastSlice = -1;
  while ( sdIter != sdEnd )
    {
      if ( sdIter->fVertexBlock )
	{
	  fVertex->Read( (AliHLTTPCVertexData*)( sdIter->fVertexBlock->fPtr ) );
	  fGlobalMerger->SetVertex( fVertex );
	}
      for ( int slNr=lastSlice+1; slNr<=sdIter->fSlice; slNr++ )
	fGlobalMerger->InitSlice( slNr );
      if ( sdIter->fTrackletBlock )
	{
	  inPtr = (AliHLTTPCTrackletData*)( sdIter->fTrackletBlock->fPtr );
	  if ( !inPtr )
	    {
	      Logging( kHLTLogError, "HLT::GlobalMerger::DoEvent", "No track data block",
		       "No track data block found  for event  0x%08lX (%lu).", evtData.fEventID, evtData.fEventID );
	    }
	  else
	    {
	      //fGlobalMerger->InitSlice( sdIter->fSlice );
	      fGlobalMerger->FillTracks( inPtr->fTrackletCnt, inPtr->fTracklets );
	    } 
	}
      lastSlice = sdIter->fSlice;
      sdIter++;
    }
  for ( int slNr=lastSlice+1; slNr<=35; slNr++ )
    fGlobalMerger->InitSlice( slNr );
    

  // Now we can really merge
  fGlobalMerger->Merge();
  fGlobalMerger->AddAllTracks();

  // check if there was enough space in the output buffer
  UInt_t ntracks0=fGlobalMerger->GetOutTracks()->GetNTracks();
  if (outputPtr==NULL || (sizeof(AliHLTTPCTrackletData)+ntracks0*sizeof(AliHLTTPCTrackSegmentData)>capacity)) {
    iResult=-ENOSPC;
  } else {
  outPtr = (AliHLTTPCTrackletData*)(outputPtr);

  tSize = fGlobalMerger->GetOutTracks()->WriteTracks( ntracks0, outPtr->fTracklets );
  outPtr->fTrackletCnt = ntracks0;

  tSize += sizeof(AliHLTTPCTrackletData);

  AliHLTComponentBlockData bd;
  FillBlockData( bd );
  bd.fOffset = 0;
  bd.fSize = tSize;
  bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( minSlice, maxSlice, 0, 5 );
  outputBlocks.push_back( bd );
  }

  size = tSize;
  return iResult;
}

	
