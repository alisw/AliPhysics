// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
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

/** @file   AliHLTTPCVertexFinderComponent.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  TPC vertex finder processing component
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// a TPC vertex finder processing component for the HLT                      //
//                                                                           //
// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCVertexFinderComponent.h"
#include "AliHLTTPCVertexFinder.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCVertexData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDefinitions.h"
#include <cstdlib>
#include <cerrno>

// this is a global object used for automatic component registration, do not use this
AliHLTTPCVertexFinderComponent gAliHLTTPCVertexFinderComponent;

ClassImp(AliHLTTPCVertexFinderComponent);

AliHLTTPCVertexFinderComponent::AliHLTTPCVertexFinderComponent()
  :
  fVertexFinder(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCVertexFinderComponent::AliHLTTPCVertexFinderComponent(const AliHLTTPCVertexFinderComponent&)
  :
  fVertexFinder(NULL)
{
  // see header file for class documentation
}

AliHLTTPCVertexFinderComponent& AliHLTTPCVertexFinderComponent::operator=(const AliHLTTPCVertexFinderComponent&)
{ 
  // see header file for class documentation
  return *this;
}

AliHLTTPCVertexFinderComponent::~AliHLTTPCVertexFinderComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCVertexFinderComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCVertexFinder";
}

void AliHLTTPCVertexFinderComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
}

AliHLTComponentDataType AliHLTTPCVertexFinderComponent::GetOutputDataType()
{
  // see header file for class documentation
  return AliHLTTPCDefinitions::fgkVertexDataType;
}

void AliHLTTPCVertexFinderComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = sizeof(AliHLTTPCVertexData);
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTPCVertexFinderComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCVertexFinderComponent;
}

int AliHLTTPCVertexFinderComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  if ( fVertexFinder )
    return EINPROGRESS;
  fVertexFinder = new AliHLTTPCVertexFinder();
  return 0;
}

int AliHLTTPCVertexFinderComponent::DoDeinit()
{
  // see header file for class documentation
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
  // see header file for class documentation
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
      if ( iter->fDataType != AliHLTTPCDefinitions::fgkClustersDataType )
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

	
