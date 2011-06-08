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

/** @file   AliHLTTPCCompModelDeflaterComponent.cxx
    @author Timm Steinbeck
    @date   
    @brief  A copy processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompModelDeflaterComponent.h"
#include "AliHLTTPCDefinitions.h"
#include <stdlib.h>
#include <errno.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCompModelDeflaterComponent);
    
AliHLTTPCCompModelDeflaterComponent::AliHLTTPCCompModelDeflaterComponent():
  fModelDeflater(),
  fConverter(),
  fForwardIfUncompressed(true)
    {
      // see header file for class documentation
    }

AliHLTTPCCompModelDeflaterComponent::~AliHLTTPCCompModelDeflaterComponent()
    {
      // see header file for class documentation
    }

const char* AliHLTTPCCompModelDeflaterComponent::GetComponentID()
    {
      // see header file for class documentation
      return "TPCCompModelDeflater"; // The ID of this component
    }

void AliHLTTPCCompModelDeflaterComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
    {
      // see header file for class documentation
      list.clear(); // We do not have any requirements for our input data type(s).
      list.push_back( AliHLTTPCDefinitions::fgkClusterTracksModelDataType );
      list.push_back( AliHLTTPCDefinitions::fgkRemainingClustersModelDataType );
    }

AliHLTComponent_DataType AliHLTTPCCompModelDeflaterComponent::GetOutputDataType()
    {
      // see header file for class documentation
      return AliHLTTPCDefinitions::fgkClusterTracksCompressedDataType;
    }

void AliHLTTPCCompModelDeflaterComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
      // see header file for class documentation
      constBase = 1+1+216; // Format versions + cluster count per patch
      inputMultiplier = 4.;
      //#warning Adapt input Multiplier to something more realistic
    }

// Spawn function, return new instance of this class
AliHLTComponent* AliHLTTPCCompModelDeflaterComponent::Spawn()
    {
      // see header file for class documentation
      return new AliHLTTPCCompModelDeflaterComponent;
    };

int AliHLTTPCCompModelDeflaterComponent::DoInit( int argc, const char** argv )
    {
      // see header file for class documentation
      //char* cpErr;
      if ( argc )
	{
	  Logging( kHLTLogDebug, "HLT::TPCCompModelDeflater::DoInit", "Arguments", "argv[0] == %s", argv[0] );
	  Logging(kHLTLogError, "HLT::TPCCompModelDeflater::DoInit", "Unknown Option", "Unknown option '%s'", argv[0] );
	  return EINVAL;
	}
    return 0;
    }

int AliHLTTPCCompModelDeflaterComponent::DoDeinit()
    {
      // see header file for class documentation
      return 0;
    }

int AliHLTTPCCompModelDeflaterComponent::DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
						   AliHLTComponent_TriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
    {
      // see header file for class documentation
      // Process an event
      // Loop over all input blocks in the event
      AliHLTUInt8_t* writePtr = outputPtr;
      AliHLTUInt32_t outputSize = 0, blockSize;
      int ret;
      AliHLTComponent_BlockData ob;
      AliHLTUInt8_t minSlice=0xFF, maxSlice=0xFF, minPatch=0xFF, maxPatch=0xFF;
      fConverter.Init();
      unsigned long long totalNonModelDataSize=0;
      
      for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
	{
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClusterTracksModelDataType )
	    {
	      blockSize = size-outputSize;
	      ret = fModelDeflater.CompressTracks( (AliHLTUInt8_t*)blocks[n].fPtr, blocks[n].fSize, writePtr, blockSize );
	      if ( !ret && blockSize<blocks[n].fSize )
		{
		  // Let the structure be filled with the default values.
		  // This takes care of setting the shared memory and data type values to default values,
		  // so that they can be filled in by the calling code.
		  FillBlockData( ob );
		  // This block's start (offset) is after all other blocks written so far
		  ob.fOffset = outputSize;
		  // the size of this block's data.
		  ob.fSize = blockSize;
		  // The specification of the data is copied from the input block.
		  ob.fSpecification = blocks[n].fSpecification;
		  // The type of the data is copied from the input block.
		  ob.fDataType = AliHLTTPCDefinitions::fgkClusterTracksCompressedDataType;
		  // Place this block into the list of output blocks
		  outputBlocks.push_back( ob );
		  writePtr += blockSize;
		outputSize += blockSize;
		}
	      else if ( fForwardIfUncompressed )
		{
		  outputBlocks.push_back( blocks[n] );
		}
	      continue;
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkRemainingClustersModelDataType )
	    {
	    blockSize = size-outputSize;
	    ret = fModelDeflater.CompressRemainingClusters( (AliHLTUInt8_t*)blocks[n].fPtr, blocks[n].fSize, writePtr, blockSize );
	    HLTDebug( "ret: %d - blockSize: %u - blocks[%u].fSize: %u", ret, (unsigned)blockSize, (unsigned)n, (unsigned)blocks[n].fSize );
	    if ( !ret && blockSize<blocks[n].fSize )
	      {
		// Let the structure be filled with the default values.
		// This takes care of setting the shared memory and data type values to default values,
		// so that they can be filled in by the calling code.
		FillBlockData( ob );
		// This block's start (offset) is after all other blocks written so far
		ob.fOffset = outputSize;
		// the size of this block's data.
		ob.fSize = blockSize;
		// The specification of the data is copied from the input block.
		ob.fSpecification = blocks[n].fSpecification;
		// The type of the data is copied from the input block.
		ob.fDataType = AliHLTTPCDefinitions::fgkRemainingClustersCompressedDataType;
		// Place this block into the list of output blocks
		outputBlocks.push_back( ob );
		writePtr += blockSize;
		outputSize += blockSize;
	      }
	    else if ( fForwardIfUncompressed )
	      {
		outputBlocks.push_back( blocks[n] );
	      }
	    continue;
	    }
	  AliHLTUInt8_t slice = 0;
	  AliHLTUInt8_t patch = 0;
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClustersDataType ||
	       blocks[n].fDataType == AliHLTTPCDefinitions::fgkTracksDataType )
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
	      totalNonModelDataSize += blocks[n].fSize;
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClustersDataType )
	    {
	    fConverter.SetInputClusters( (AliHLTTPCClusterData*)blocks[n].fPtr, slice, patch );
	    }
	  if ( blocks[n].fDataType == (kAliHLTDataTypeTrack|kAliHLTDataOriginTPC) )
	    {
	      fConverter.SetInputTracks( (AliHLTTracksData*)blocks[n].fPtr, blocks[n].fSize );
	    }
	  
	}
      
      if ( totalNonModelDataSize>0 )
	{
	  fConverter.Convert();
	  
	  unsigned long trackSize = fConverter.GetOutputModelDataSize();
	  AliHLTUInt8_t* trackModelData = new AliHLTUInt8_t[ trackSize ];
	  if ( !trackModelData )
	    {
	      HLTError( "Out of memory trying to allocate %lu bytes of trackmodeldata", trackSize );
	      return ENOMEM;
	    }
	  
	  fConverter.OutputModelData( trackModelData, trackSize );
	  
	  unsigned long clusterSize = fConverter.GetRemainingClustersOutputDataSize();
	  AliHLTUInt8_t* remainingClustersModelData = new AliHLTUInt8_t[ clusterSize ];
	  if ( !remainingClustersModelData )
	    {
	      HLTError( "Out of memory trying to allocate %lu bytes of remaining cluster model data", clusterSize );
	      delete [] trackModelData;
	      return ENOMEM;
	    }
	  
	  fConverter.GetRemainingClusters( remainingClustersModelData, clusterSize );
	  
	  bool forwardUncompressed = false;
	  
	  blockSize = size-outputSize;
	  ret = fModelDeflater.CompressTracks( trackModelData, trackSize, writePtr, blockSize );
	  unsigned long long totalCompressedModelData = blockSize;
	  if ( !ret )
	    {
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
	      ob.fDataType = AliHLTTPCDefinitions::fgkClusterTracksCompressedDataType;
	      // Place this block into the list of output blocks
	      outputBlocks.push_back( ob );
	      writePtr += blockSize;
	      outputSize += blockSize;
	      
	      blockSize = size-outputSize;
	      ret = fModelDeflater.CompressRemainingClusters( remainingClustersModelData, clusterSize, writePtr, blockSize );
	      totalCompressedModelData += blockSize;
	      if ( !ret && totalCompressedModelData<totalNonModelDataSize )
		{
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
		  ob.fDataType = AliHLTTPCDefinitions::fgkRemainingClustersCompressedDataType;
		  // Place this block into the list of output blocks
		  outputBlocks.push_back( ob );
		  writePtr += blockSize;
		  outputSize += blockSize;
		}
	      else if ( fForwardIfUncompressed )
		{
		  outputSize -= (outputBlocks.end()-1)->fSize;
		  outputBlocks.erase( outputBlocks.end()-1 );
		  forwardUncompressed = true;
		}
	      
	      
	    }
	  else if ( fForwardIfUncompressed )
	    forwardUncompressed = true;
	  
	  delete [] trackModelData;
	  delete [] remainingClustersModelData;
	  
	  if ( forwardUncompressed )
	    {
	      for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
		{
		  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClustersDataType ||
		       blocks[n].fDataType == AliHLTTPCDefinitions::fgkTracksDataType )
		    {
		      outputBlocks.push_back( blocks[n] );
		    }
		  
		}
	    }
	}
      
      // Finally we set the total size of output memory we consumed.
      size = outputSize;
      return 0;
    }
