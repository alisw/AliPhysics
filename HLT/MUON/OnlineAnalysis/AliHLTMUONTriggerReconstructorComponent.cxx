/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/**
 * @file   AliHLTMUONTriggerReconstructorComponent.cxx
 * @author Indranil Das <indra.das@saha.ac.in>, Artur Szostak <artursz@iafrica.com>
 * @date
 * @brief  Implementation of the trigger DDL reconstructor component.
 */

#include "AliHLTMUONTriggerReconstructorComponent.h"
#include "AliHLTMUONTriggerReconstructor.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockWriter.h"
#include <cstdlib>
#include <cerrno>
#include <cassert>
#include <fstream>

namespace
{
	// This is a global object used for automatic component registration,
	// do not use this for calculation.
	AliHLTMUONTriggerReconstructorComponent gAliHLTMUONTriggerReconstructorComponent;

} // end of namespace


ClassImp(AliHLTMUONTriggerReconstructorComponent)


AliHLTMUONTriggerReconstructorComponent::AliHLTMUONTriggerReconstructorComponent() :
	fTrigRec(NULL),
	fDDL(-1),
	fWarnForUnexpecedBlock(false),
	fSuppressPartialTrigs(false)
{
}


AliHLTMUONTriggerReconstructorComponent::~AliHLTMUONTriggerReconstructorComponent()
{
}


const char* AliHLTMUONTriggerReconstructorComponent::GetComponentID()
{
	return AliHLTMUONConstants::TriggerReconstructorId();
}


void AliHLTMUONTriggerReconstructorComponent::GetInputDataTypes( std::vector<AliHLTComponentDataType>& list)
{
	list.clear();
	list.push_back( AliHLTMUONConstants::TriggerDDLRawDataType() );
}


AliHLTComponentDataType AliHLTMUONTriggerReconstructorComponent::GetOutputDataType()
{
	return AliHLTMUONConstants::TriggerRecordsBlockDataType();
}


void AliHLTMUONTriggerReconstructorComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	constBase = sizeof(AliHLTMUONTriggerRecordsBlockWriter::HeaderType);
	inputMultiplier = 4;
}


AliHLTComponent* AliHLTMUONTriggerReconstructorComponent::Spawn()
{
	return new AliHLTMUONTriggerReconstructorComponent;
}


int AliHLTMUONTriggerReconstructorComponent::DoInit(int argc, const char** argv)
{
	// perform initialization.
	
	HLTInfo("Initialising dHLT trigger reconstructor component.");
	
	fWarnForUnexpecedBlock = false;
	fSuppressPartialTrigs = false;
	assert(fTrigRec == NULL);
	fTrigRec = new AliHLTMUONTriggerReconstructor();
	
	const char* lutFileName = NULL;
	
	for (int i = 0; i < argc; i++)
	{
		if ( !strcmp( argv[i], "-lut" ) )
		{
			if ( argc <= i+1 )
			{
				HLTError("LookupTable filename not specified." );
				return EINVAL; /* Invalid argument */ 
			}
			
			lutFileName = argv[i+1];
			
			i++;
			continue;
		}
		
		if ( !strcmp( argv[i], "-ddl" ) )
		{
			if ( argc <= i+1 )
			{
				HLTError("DDL number not specified." );
				return EINVAL;  /* Invalid argument */
			}
		
			char* cpErr = NULL;
			unsigned long num = strtoul(argv[i+1], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a DDL Number.", argv[i+1] );
				return EINVAL;
			}
			if (num < 21 or 22 < num)
			{
				HLTError("The DDL number must be in the range [21..22].");
				return EINVAL;
			}
			fDDL = num - 1; // Convert to DDL number in the range 0..21
			
			i++;
			continue;
		}
			
		if (not strcmp( argv[i], "-warn_on_unexpected_block" ))
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}
			
		if (not strcmp( argv[i], "-suppress_partial_triggers" ))
		{
			fSuppressPartialTrigs = true;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i] );
		return EINVAL;
			
	}//while loop
	
	if (fDDL == -1)
	{
		HLTWarning("DDL number not specified. Cannot check if incomming data is valid.");
	}
	
	if (lutFileName != NULL)
	{
		if (not ReadLookUpTable(lutFileName))
		{
			HLTError("Failed to read lut, lut cannot be read");
			return ENOENT ; /* No such file or directory */
		}
	}
	else
	{
		HLTWarning("The lookup table has not been specified. Output results will be invalid.");
	}
	
	return 0;
}


int AliHLTMUONTriggerReconstructorComponent::DoDeinit()
{
	HLTInfo("Deinitialising dHLT trigger reconstructor component.");

	if (fTrigRec != NULL)
	{
		delete fTrigRec;
		fTrigRec = NULL;
	}
	return 0;
}


int AliHLTMUONTriggerReconstructorComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		std::vector<AliHLTComponentBlockData>& outputBlocks
	)
{
	// Process an event
	unsigned long totalSize = 0; // Amount of memory currently consumed in bytes.

	HLTDebug("Processing event %llu with %u input data blocks.",
		evtData.fEventID, evtData.fBlockCnt
	);
	
	// Loop over all input blocks in the event and run the trigger DDL
	// reconstruction algorithm on the raw data.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
#ifdef __DEBUG
		char id[kAliHLTComponentDataTypefIDsize+1];
		for (int i = 0; i < kAliHLTComponentDataTypefIDsize; i++)
			id[i] = blocks[n].fDataType.fID[i];
		id[kAliHLTComponentDataTypefIDsize] = '\0';
		char origin[kAliHLTComponentDataTypefOriginSize+1];
		for (int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++)
			origin[i] = blocks[n].fDataType.fOrigin[i];
		origin[kAliHLTComponentDataTypefOriginSize] = '\0';
#endif // __DEBUG
		HLTDebug("Handling block: %u, with fDataType.fID = '%s',"
			  " fDataType.fID = '%s', fPtr = %p and fSize = %u bytes.",
			n, static_cast<char*>(id), static_cast<char*>(origin),
			blocks[n].fPtr, blocks[n].fSize
		);

		if (blocks[n].fDataType != AliHLTMUONConstants::TriggerDDLRawDataType())
		{
			// Log a message indicating that we got a data block that we
			// do not know how to handle.
			char id[kAliHLTComponentDataTypefIDsize+1];
			for (int i = 0; i < kAliHLTComponentDataTypefIDsize; i++)
				id[i] = blocks[n].fDataType.fID[i];
			id[kAliHLTComponentDataTypefIDsize] = '\0';
			char origin[kAliHLTComponentDataTypefOriginSize+1];
			for (int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++)
				origin[i] = blocks[n].fDataType.fOrigin[i];
			origin[kAliHLTComponentDataTypefOriginSize] = '\0';
			
			if (fWarnForUnexpecedBlock)
				HLTWarning("Received a data block of a type we can not handle: %s origin %s",
					static_cast<char*>(id), static_cast<char*>(origin)
				);
			else
				HLTDebug("Received a data block of a type we can not handle: %s origin %s",
					static_cast<char*>(id), static_cast<char*>(origin)
				);
			
			continue;
		}
		
		bool ddl[22];
		AliHLTMUONUtils::UnpackSpecBits(blocks[n].fSpecification, ddl);
		if (not ddl[fDDL])
		{
			HLTWarning("Received trigger DDL raw data from a DDL which we did not expect.");
		}
		
		// Create a new output data block and initialise the header.
		AliHLTMUONTriggerRecordsBlockWriter block(outputPtr+totalSize, size-totalSize);
		if (not block.InitCommonHeader())
		{
			HLTError("There is not enough space in the output buffer for the new data block.",
				 " We require at least %ufTrigRec->GetkDDLHeaderSize() bytes, but have %u bytes left.",
				sizeof(AliHLTMUONTriggerRecordsBlockWriter::HeaderType),
				block.BufferSize()
			);
			break;
		}

		AliHLTUInt32_t totalDDLSize = blocks[n].fSize / sizeof(AliHLTUInt32_t);
		AliHLTUInt32_t ddlRawDataSize = totalDDLSize - 8;
		AliHLTUInt32_t* buffer = reinterpret_cast<AliHLTUInt32_t*>(blocks[n].fPtr) + 8;
		AliHLTUInt32_t nofTrigRec = block.MaxNumberOfEntries();

		bool runOk = fTrigRec->Run(
				buffer, ddlRawDataSize,
				block.GetArray(), nofTrigRec,
				fSuppressPartialTrigs
			);
		if (not runOk)
		{
			HLTError("Error while processing of trigger DDL reconstruction algorithm.");
			size = totalSize; // Must tell the framework how much buffer space was used.
			return EIO;
		}
		
		// nofTrigRec should now contain the number of triggers actually found
		// and filled into the output data block, so we can set this number.
		assert( nofTrigRec <= block.MaxNumberOfEntries() );
		block.SetNumberOfEntries(nofTrigRec);
		
		HLTDebug("Number of trigger records found is %d", nofTrigRec);
		
		// Fill a block data structure for our output block.
		AliHLTComponentBlockData bd;
		FillBlockData(bd);
		bd.fPtr = outputPtr;
		// This block's start (offset) is after all other blocks written so far.
		bd.fOffset = totalSize;
		bd.fSize = block.BytesUsed();
		bd.fDataType = AliHLTMUONConstants::TriggerRecordsBlockDataType();
		bd.fSpecification = blocks[n].fSpecification;
		outputBlocks.push_back(bd);
		
		HLTDebug("Created a new output data block at fPtr = %p,"
			  " with fOffset = %u (0x%.X) and fSize = %u bytes.", 
			bd.fPtr, bd.fOffset, bd.fOffset, bd.fSize
		);
		
		// Increase the total amount of data written so far to our output memory.
		totalSize += block.BytesUsed();
	}
	
	// Finally we set the total size of output memory we consumed.
	size = totalSize;
	return 0;
}


bool AliHLTMUONTriggerReconstructorComponent::ReadLookUpTable(const char* lutpath)
{
	assert(fTrigRec != NULL);

	fstream file;
	file.open(lutpath, fstream::binary | fstream::in);
	if (not file)
	{
		HLTError("Could not open file: %s", lutpath);
		return false;
	}
	
	file.read(reinterpret_cast<char*>(fTrigRec->LookupTableBuffer()), fTrigRec->LookupTableSize());
	if (file.eof())
	{
		HLTError("The file %s was too short to contain a valid lookup table for this component.", lutpath);
		file.close();
		return false;
	}
	if (file.bad())
	{
		HLTError("Could not read from file: %s", lutpath);
		file.close();
		return false;
	}
	
	file.close();
	return true;
}
