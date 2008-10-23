// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
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

/** @file   AliHLTTPCCompModelDeconverterComponent.cxx
    @author Timm Steinbeck
    @date   
    @brief  A copy processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompModelDeconverterComponent.h"
#include "AliHLTTPCDefinitions.h"
#include <stdlib.h>
#include <errno.h>

/**
 * An implementiation of a deconverter component that 
 * deconverts the tracks and clusters from the Vestbo-model
 * into the standard HLT cluster track format again 
 * in order to evaluate the loss of the model 
 * due to the Vestbo-compression 
 */

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCompModelDeconverterComponent)
    
AliHLTTPCCompModelDeconverterComponent::AliHLTTPCCompModelDeconverterComponent():
  fDeconverter(),
  fOutputTracks(kTRUE)
    {
      // see header file for class documentation
    }

AliHLTTPCCompModelDeconverterComponent::~AliHLTTPCCompModelDeconverterComponent()
    {
      // see header file for class documentation
    }

const char* AliHLTTPCCompModelDeconverterComponent::GetComponentID()
    {
      // see header file for class documentation
      return "TPCCompModelDeconverter"; // The ID of this component
    }

void AliHLTTPCCompModelDeconverterComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
    {
      // see header file for class documentation
      list.clear(); // We do not have any requirements for our input data type(s).
      list.push_back( AliHLTTPCDefinitions::fgkClusterTracksModelDataType );
      list.push_back( AliHLTTPCDefinitions::fgkRemainingClustersModelDataType );
    }

AliHLTComponent_DataType AliHLTTPCCompModelDeconverterComponent::GetOutputDataType()
    {
      // see header file for class documentation
      return AliHLTTPCDefinitions::fgkClustersDataType;
    }

void AliHLTTPCCompModelDeconverterComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
      // see header file for class documentation
      constBase = 8+216*4; // Track count + clusters count
      inputMultiplier = 4.;
    }

// Spawn function, return new instance of this class
AliHLTComponent* AliHLTTPCCompModelDeconverterComponent::Spawn()
    {
      // see header file for class documentation
      return new AliHLTTPCCompModelDeconverterComponent;
    };

int AliHLTTPCCompModelDeconverterComponent::DoInit( int argc, const char** argv )
    {
      // see header file for class documentation
      Int_t i = 0;
      //Char_t* cpErr;
      
      while ( i < argc )
	{      
	  if ( !strcmp( argv[i], "notracks" ) )
	    {
	      fOutputTracks = kFALSE;
	      ++i;
	      continue;
	    }
	  Logging(kHLTLogError, "HLT::TPCCompModelDeconverter::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
	  return EINVAL;
	}
      return 0;
    }

int AliHLTTPCCompModelDeconverterComponent::DoDeinit()
    {
      // see header file for class documentation
      return 0;
    }

int AliHLTTPCCompModelDeconverterComponent::DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
						   AliHLTComponent_TriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
    {
      // see header file for class documentation
      fDeconverter.Init();
      // Process an event
      // Loop over all input blocks in the event
      AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
      for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
	{
	  AliHLTUInt8_t slice, patch;
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkRemainingClustersModelDataType ||
	       blocks[n].fDataType == AliHLTTPCDefinitions::fgkClusterTracksModelDataType )
	    {
	      slice = AliHLTTPCDefinitions::GetMinSliceNr( blocks[n].fSpecification );
	      patch = AliHLTTPCDefinitions::GetMinPatchNr( blocks[n].fSpecification );
	      if ( minSlice==0xFF || slice<minSlice )
		minSlice = slice;
	      if ( maxSlice==0xFF || slice>maxSlice )
		maxSlice = slice;
	      if ( minPatch==0xFF || patch<minPatch )
		minPatch = patch;
	      if ( maxPatch==0xFF || patch>maxPatch )
		maxPatch = patch;
	      HLTDebug( "Slice: %u - Patch: %u", (unsigned)slice, (unsigned)patch );
	      slice = AliHLTTPCDefinitions::GetMaxSliceNr( blocks[n].fSpecification );
	      patch = AliHLTTPCDefinitions::GetMaxPatchNr( blocks[n].fSpecification );
	      if ( minSlice==0xFF || slice<minSlice )
		minSlice = slice;
	      if ( maxSlice==0xFF || slice>maxSlice )
		maxSlice = slice;
	      if ( minPatch==0xFF || patch<minPatch )
		minPatch = patch;
	      if ( maxPatch==0xFF || patch>maxPatch )
		maxPatch = patch;
	      HLTDebug( "Slice: %u - Patch: %u", (unsigned)slice, (unsigned)patch );
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClusterTracksModelDataType )
	    {
	      HLTDebug( "Tracks" );
	      fDeconverter.SetTrackClusterModelInputData( (AliHLTUInt8_t*)blocks[n].fPtr, blocks[n].fSize );
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkRemainingClustersModelDataType )
	    {
	      HLTDebug( "Clusters" );
	      fDeconverter.SetRemainingClustersModelInputData( (AliHLTUInt8_t*)blocks[n].fPtr, blocks[n].fSize );
	    }
	}
      
      HLTDebug( "min slice: %u - max slice: %u - min patch: %u - max patch: %u",
		(unsigned)minSlice, (unsigned)maxSlice, (unsigned)minPatch, (unsigned)maxPatch );
      
      UInt_t blockSize = size;
      UInt_t outputSize = 0;
      Int_t ret;
      if ( fOutputTracks )
	{
	  ret = fDeconverter.DeconvertTracks( outputPtr, blockSize );
	  if ( !ret )
	    {
	      if ( outputSize+blockSize > size )
		{
		  HLTError( "Output data too large. (%lu used instead of %u available)",
			    (unsigned long)blockSize, (unsigned long)size );
		  return ENOBUFS;
		}
	      
	      AliHLTComponent_BlockData ob;
	      // Let the structure be filled with the default values.
	      // This takes care of setting the shared memory and data type values to default values,
	      // so that they can be filled in by the calling code.
	      FillBlockData( ob );
	      // This block's start (offset) is after all other blocks written so far
	      ob.fOffset = outputSize;
	      // the size of this block's data.
	      ob.fSize = blockSize;
	      // The specification of the data is copied from the input block.
	      ob.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( minSlice, maxSlice, minPatch, maxPatch );
	      // The type of the data is copied from the input block.
	      ob.fDataType = AliHLTTPCDefinitions::fgkTracksDataType;
	      // Place this block into the list of output blocks
	      outputBlocks.push_back( ob );
	      outputSize += blockSize;
	    }
	  else
	    HLTError( "Error deconverting tracks: %s (%d)", strerror(ret), (int)ret );
	}
      
      for ( UInt_t slice=minSlice; slice<=maxSlice; slice++ )
	{
	  for ( UInt_t patch=minPatch; patch<=maxPatch; patch++ )
	    {
	      blockSize = size-outputSize;
	      ret = fDeconverter.DeconvertClusters( slice, patch, outputPtr+outputSize, blockSize );
	      if ( !ret )
		{
		  if ( outputSize+blockSize > size )
		    {
		      HLTError( "Output data too large. (%lu used instead of %u available)",
				(unsigned long)blockSize, (unsigned long)size );
		      return ENOBUFS;
		    }
		  
		  AliHLTComponent_BlockData ob;
		  // Let the structure be filled with the default values.
		  // This takes care of setting the shared memory and data type values to default values,
		  // so that they can be filled in by the calling code.
		  FillBlockData( ob );
		  // This block's start (offset) is after all other blocks written so far
		  ob.fOffset = outputSize;
		  // the size of this block's data.
		  ob.fSize = blockSize;
		  // The specification of the data is copied from the input block.
		  ob.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( slice, slice, patch, patch );
		  // The type of the data is copied from the input block.
		  ob.fDataType = AliHLTTPCDefinitions::fgkClustersDataType;
		  // Place this block into the list of output blocks
		  outputBlocks.push_back( ob );
		  outputSize += blockSize;
		}
	      else
		HLTError( "Error deconverting clusters for slice %u, patch %u: %s (%d)", 
			  (unsigned)slice, (unsigned)patch, strerror(ret), (int)ret );
	    }
	}
      
      // Finally we set the total size of output memory we consumed.
      size = outputSize;
      return 0;
    }
