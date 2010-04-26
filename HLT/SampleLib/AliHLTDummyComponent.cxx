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

/** @file   AliHLTDummyComponent.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  A dummy processing component for the HLT. */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTSystem.h"
#include "AliHLTDummyComponent.h"
#include "AliHLTDefinitions.h"
#include <cstdlib>
#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDummyComponent)
    
AliHLTDummyComponent::AliHLTDummyComponent()
  :
    fOutputPercentage(100) // By default we copy to the output exactly what we got as input
    {
    // see header file for class documentation
    // or
    // refer to README to build package
    // or
    // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
    }

AliHLTDummyComponent::~AliHLTDummyComponent()
    {
    // see header file for class documentation
    }

const char* AliHLTDummyComponent::GetComponentID()
    {
    // see header file for class documentation
    return "Dummy"; // The ID of this component
    }

void AliHLTDummyComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
    {
    // see header file for class documentation
      list.push_back(kAliHLTAnyDataType); // We do not have any requirements for our input data type(s).

    }

AliHLTComponentDataType AliHLTDummyComponent::GetOutputDataType()
    {
    // see header file for class documentation
      /* in order to be backward compatible we have to keep the old code, at
       * least for a while. Remember to use the new const kAliHLTVoidDataType
       * if you are using a more recent AliRoot version (from Jan 07)
      return kAliHLTVoidDataType;
      */
      AliHLTComponentDataType dt = 
	{ sizeof(AliHLTComponentDataType),
	  {'\0','\0','\0','0','\0','\0','\0','\0'},
	  {'\0','\0','\0','\0'}};
      return dt;
    }

void AliHLTDummyComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
    {
    // see header file for class documentation
    constBase = 0;
    inputMultiplier = ((double)fOutputPercentage)/100.0;
    }



// Spawn function, return new instance of this class
AliHLTComponent* AliHLTDummyComponent::Spawn()
    {
    // see header file for class documentation
    return new AliHLTDummyComponent;
    }

int AliHLTDummyComponent::DoInit( int argc, const char** argv )
    {
    // perform initialization. We check whether our relative output size is specified in the arguments.
    fOutputPercentage = 100;
    int i = 0;
    char* cpErr;
    while ( i < argc )
	{
	HLTDebug("argv[%d] == %s", i, argv[i] );
	if ( !strcmp( argv[i], "output_percentage" ) ||
	     !strcmp( argv[i], "-output_percentage" ))
	    {
	    if ( i+1>=argc )
		{
		HLTError("Missing output_percentage parameter");
		return -EINVAL;
		}
	    HLTDebug("argv[%d+1] == %s", i, argv[i+1] );
	    fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	    if ( *cpErr )
		{
		HLTError("Cannot convert output_percentage parameter '%s'", argv[i+1] );
		return -EINVAL;
		}
	    HLTInfo("Output percentage set to %lu %%", fOutputPercentage );
	    i += 2;
	    continue;
	    }
	HLTError("Unknown option '%s'", argv[i] );
	return -EINVAL;
	}
    return 0;
    }

int AliHLTDummyComponent::DoDeinit()
    {
    // see header file for class documentation
    return 0;
    }

int AliHLTDummyComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
				      AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
    {
    // see header file for class documentation
    HLTInfo("Output percentage set to %lu %%", fOutputPercentage );
    // Process an event
    unsigned long totalSize = 0;
    // Loop over all input blocks in the event
    for ( unsigned long n = 0; n < evtData.fBlockCnt; n++ )
	{
	// Align the beginning of this  block to the required value.
	if ( totalSize % kAliHLTBlockAlignment )
	    totalSize += kAliHLTBlockAlignment-(totalSize % kAliHLTBlockAlignment);
	if ( totalSize > size )
	    break;
	// Determine the size we should use for the output for this block (the input block's size times the relative output size)
	unsigned long mySize = (blocks[n].fSize * fOutputPercentage) / 100;
	HLTInfo("HLT::Dummy::DoEvent", "mySize set (1)", "mySize == %lu B - blocks[%lu].fSize == %lu - fOutputPercentage == %lu", 
		 mySize, n, blocks[n].fSize, fOutputPercentage );
	// Check how much space we have left and adapt this output block's size accordingly.
	if ( totalSize + mySize > size )
	    mySize = size-totalSize;
	HLTInfo("HLT::Dummy::DoEvent", "mySize set (2)", "mySize == %lu B - totalSize == %lu - size == %lu", 
		 mySize, totalSize, size );
	if ( mySize<=0 )
	    continue; // No room left to write a further block.
	// Now copy the input block
	unsigned long copied = 0;
	// First copy all full multiples of the input block
	while ( copied+blocks[n].fSize <= mySize )
	    {
	    HLTInfo("Copying %lu B - Copied: %lu B - totalSize: %lu B", 
		     blocks[n].fSize, copied, totalSize );
	    memcpy( outputPtr+totalSize+copied, blocks[n].fPtr, blocks[n].fSize );
	    copied += blocks[n].fSize;
	    }
	// And the copy the remaining fragment of the block
	HLTInfo("Copying %lu B - Copied: %lu B - totalSize: %lu B", 
		 mySize-copied, copied, totalSize );
	memcpy( outputPtr+totalSize+copied, blocks[n].fPtr, mySize-copied );
	HLTInfo("Copied: %lu B - totalSize: %lu B", 
		 copied, totalSize );
	// Fill a block data structure for our output block.
	AliHLTComponentBlockData ob;
	// Let the structure be filled with the default values.
	// This takes care of setting the shared memory and data type values to default values,
	// so that they can be filled in by the calling code.
	FillBlockData( ob );
	// This block's start (offset) is after all other blocks written so far
	ob.fOffset = totalSize;
	// the size of this block's data.
	ob.fSize = mySize;
	// the data type of this block
	ob.fDataType=blocks[n].fDataType;
	// The specification of the data is copied from the input block.
	ob.fSpecification = blocks[n].fSpecification;
	// The data type is set automatically to the component's specified output data type.
	// Place this block into the list of output blocks
	outputBlocks.push_back( ob );
	// Increase the total amount of data written so far to our output memory
	totalSize += mySize;
	}
    // Finally we set the total size of output memory we consumed.
    size = totalSize;
    return 0;
    }
