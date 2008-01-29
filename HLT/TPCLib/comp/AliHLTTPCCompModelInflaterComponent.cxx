// $Id$

/**************************************************************************
 * TPCCompModelInflaterright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
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

/** @file   AliHLTTPCCompModelInflaterComponent.cxx
    @author Timm Steinbeck
    @date   
    @brief  A copy processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompModelInflaterComponent.h"
#include "AliHLTTPCDefinitions.h"
#include <stdlib.h>
#include <errno.h>

// this is a global object used for automatic component registration, do not use this
AliHLTTPCCompModelInflaterComponent gAliHLTTPCCompClusterModelInflaterComponent;

ClassImp(AliHLTTPCCompModelInflaterComponent);
    
AliHLTTPCCompModelInflaterComponent::AliHLTTPCCompModelInflaterComponent()
  :
  fModelInflater()
    {
      // see header file for class documentation
    }

AliHLTTPCCompModelInflaterComponent::~AliHLTTPCCompModelInflaterComponent()
    {
      // see header file for class documentation
    }

const char* AliHLTTPCCompModelInflaterComponent::GetComponentID()
    {
      // see header file for class documentation
      return "TPCCompModelInflater"; // The ID of this component
    }

void AliHLTTPCCompModelInflaterComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
    {
      // see header file for class documentation
      list.clear(); // We do not have any requirements for our input data type(s).
      list.push_back( AliHLTTPCDefinitions::fgkClusterTracksCompressedDataType );
      list.push_back( AliHLTTPCDefinitions::fgkRemainingClustersCompressedDataType );
    }

AliHLTComponent_DataType AliHLTTPCCompModelInflaterComponent::GetOutputDataType()
    {
      // see header file for class documentation 
      return AliHLTTPCDefinitions::fgkClusterTracksModelDataType;
    }

void AliHLTTPCCompModelInflaterComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
      // see header file for class documentation
      constBase = 4+4; // Format versions
      inputMultiplier = 20.;
      //#warning Adapt input Multiplier to something more realistic
    }

// Spawn function, return new instance of this class
AliHLTComponent* AliHLTTPCCompModelInflaterComponent::Spawn()
    {
      // see header file for class documentation
      return new AliHLTTPCCompModelInflaterComponent;
    };

int AliHLTTPCCompModelInflaterComponent::DoInit( int argc, const char** argv )
    {
      // see header file for class documentation
      //char* cpErr;
      if ( argc )
	{
	  Logging( kHLTLogDebug, "HLT::TPCCompModelInflater::DoInit", "Arguments", "argv[0] == %s", argv[0] );
	  Logging(kHLTLogError, "HLT::TPCCompModelInflater::DoInit", "Unknown Option", "Unknown option '%s'", argv[0] );
	  return EINVAL;
	}
      return 0;
    }

int AliHLTTPCCompModelInflaterComponent::DoDeinit()
    {
      // see header file for class documentation
      return 0;
    }

int AliHLTTPCCompModelInflaterComponent::DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
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
      
      for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
	{
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkClusterTracksCompressedDataType )
	    {
	      blockSize = size-outputSize;
	      ret = fModelInflater.DecompressTracks( (AliHLTUInt8_t*)blocks[n].fPtr, blocks[n].fSize, writePtr, blockSize );
	      HLTDebug( "fModelInflater.DecompressTracks: ret: %d - blockSize: %u", ret, (unsigned)blockSize );
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
		  ob.fSpecification = blocks[n].fSpecification;
		  // The type of the data is copied from the input block.
		  ob.fDataType = AliHLTTPCDefinitions::fgkClusterTracksModelDataType;
		  // Place this block into the list of output blocks
		  outputBlocks.push_back( ob );
		  writePtr += blockSize;
		  outputSize += blockSize;
		}
	      continue;
	    }
	  if ( blocks[n].fDataType == AliHLTTPCDefinitions::fgkRemainingClustersCompressedDataType )
	    {
	      blockSize = size-outputSize;
	      ret = fModelInflater.DecompressRemainingClusters( (AliHLTUInt8_t*)blocks[n].fPtr, blocks[n].fSize, writePtr, blockSize );
	      HLTDebug( "ret: %d - blockSize: %u - blocks[%u].fSize: %u", ret, (unsigned)blockSize, (unsigned)n, (unsigned)blocks[n].fSize );
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
		  ob.fSpecification = blocks[n].fSpecification;
		  // The type of the data is copied from the input block.
		  ob.fDataType = AliHLTTPCDefinitions::fgkRemainingClustersModelDataType;
		  // Place this block into the list of output blocks
		  outputBlocks.push_back( ob );
		  writePtr += blockSize;
		  outputSize += blockSize;
		}
	      continue;
	    }
	}
      
      // Finally we set the total size of output memory we consumed.
      size = outputSize;
      return 0;
    }
