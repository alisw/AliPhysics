// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
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

/// @file   AliHLTTPCDataCheckerComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   9 Aug 2010
/// @brief  Implementation of the AliHLTTPCDataCheckerComponent class.
///
/// The AliHLTTPCDataCheckerComponent is used to perform data sanity and integrity
/// checks on the TPC data. This component should be used for testing and debugging.

#include "AliHLTTPCDataCheckerComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliTPCRawStreamV3.h"
#include "AliRawReaderMemory.h"
#include <vector>
#include <cstring>
#include <cassert>

ClassImp(AliHLTTPCDataCheckerComponent)


AliHLTTPCDataCheckerComponent::AliHLTTPCDataCheckerComponent() :
	fRawStream(NULL),
	fRawReader(NULL),
	fForwardBadBlocks(false),
	fForwardGoodBlocks(false),
	fIgnoreType(false),
	fIgnoreOrigin(false),
	fIgnoreSpec(false),
	fHandleAllEvents(false)
{
	// Default constructor.
}


AliHLTTPCDataCheckerComponent::~AliHLTTPCDataCheckerComponent()
{
	// Default destructor.
	
	if (fRawStream != NULL) delete fRawStream;
	if (fRawReader != NULL) delete fRawReader;
}


const char* AliHLTTPCDataCheckerComponent::GetComponentID()
{
	// Returns the component ID.
	return "TPCDataChecker";
}


void AliHLTTPCDataCheckerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	// Returns the list of input data types that are handled.
	list.push_back(kAliHLTAnyDataType | kAliHLTDataOriginTPC);
}


AliHLTComponentDataType AliHLTTPCDataCheckerComponent::GetOutputDataType()
{
	// Returns kAliHLTMultipleDataType.
	return kAliHLTMultipleDataType;
}


int AliHLTTPCDataCheckerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
	// Returns the list of output data block types generated.
	list.push_back(kAliHLTAnyDataType);
	return int(list.size());
}


void AliHLTTPCDataCheckerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
	// Returns the buffer size requirements.
	constBase = 0;
	inputMultiplier = 0;
}


AliHLTComponent* AliHLTTPCDataCheckerComponent::Spawn()
{
	// Creates a new instance of the component.
	return new AliHLTTPCDataCheckerComponent;
}


Int_t AliHLTTPCDataCheckerComponent::DoInit(int argc, const char** argv)
{
	// Initialises the data checker component from the command line.
	
	HLTInfo("Starting TPC data checker component.");
	
	fForwardBadBlocks = true;
	fForwardGoodBlocks = true;
	fIgnoreType = false;
	fIgnoreOrigin = false;
	fIgnoreSpec = false;
	fHandleAllEvents = false;
	
	for (int i = 0; i < argc; ++i)
	{
		if (strcmp(argv[i], "-filter") == 0)
		{
			if (i+1 < argc)
			{
				if (strcmp(argv[i+1], "forwardbad") == 0)
				{
					fForwardBadBlocks = true;
					fForwardGoodBlocks = false;
					++i;
				}
				else if (strcmp(argv[i+1], "forwardgood") == 0)
				{
					fForwardBadBlocks = false;
					fForwardGoodBlocks = true;
					++i;
				}
			}
			fForwardBadBlocks = true;
			fForwardGoodBlocks = false;
			continue;
		}
		
		if (strcmp(argv[i], "-ignoretype") == 0)
		{
			fIgnoreType = true;
			continue;
		}
		
		if (strcmp(argv[i], "-ignoreorigin") == 0)
		{
			fIgnoreOrigin = true;
			continue;
		}
		
		if (strcmp(argv[i], "-ignorespec") == 0)
		{
			fIgnoreSpec = true;
			continue;
		}
		
		if (strcmp(argv[i], "-handle-all-events") == 0)
		{
			fHandleAllEvents = true;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	} // for loop
	
	if (fRawReader == NULL)
	{
		fRawReader = new AliRawReaderMemory();
		if (fRawReader == NULL)
		{
			HLTError("Could not allocate new AliRawReaderMemory object.");
			return -ENOMEM;
		}
	}
	if (fRawStream == NULL)
	{
		fRawStream = new AliTPCRawStreamV3(fRawReader);
		if (fRawStream == NULL)
		{
			HLTError("Could not allocate new AliTPCRawStreamV3 object.");
			return -ENOMEM;
		}
		fRawStream->SelectRawData("TPC");
	}
	
	return 0;
}


Int_t AliHLTTPCDataCheckerComponent::DoDeinit()
{
	// Cleans up the data checker component.
	HLTInfo("Stopping TPC data checker component.");
	if (fRawReader != NULL)
	{
		fRawReader->ClearBuffers();
		fRawReader->Reset();
	}
	if (fRawStream != NULL)
	{
		fRawStream->Reset();
	}
	return 0;
}


int AliHLTTPCDataCheckerComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks, 
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* /*outputPtr*/, 
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& /*outputBlocks*/
	)
{
	// Check all the input data blocks.
	
	if (fRawReader == NULL)
	{
		HLTError("The raw reader is not set.");
		size = 0;
		return -ENOENT;
	}
	if (fRawStream == NULL)
	{
		HLTError("The TPC raw stream is not set.");
		size = 0;
		return -ENOENT;
	}
	
	if (not IsDataEvent() and not fHandleAllEvents)
	{
		// Forward all data blocks if we are not supposed to handle this event.
		for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; ++n)
		{
			Forward(blocks + n);
		}
		size = 0;
		return 0;
	}
	
	AliHLTEventID_t event = evtData.fEventID;
	AliHLTComponentDataType anyPrivateType = AliHLTComponentDataTypeInitializer(
			kAliHLTAnyDataType, kAliHLTDataOriginPrivate
		);
	
	// Setup the markers indicating the bad blocks.
	std::vector<bool> badBlock(evtData.fBlockCnt, false);
	
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; ++n)
	{
		// Skip private data blocks from the framework.
		if (blocks[n].fDataType == anyPrivateType) continue;
		
		char ddltype[kAliHLTComponentDataTypefIDsize] = kAliHLTDDLRawDataTypeID;
		if (memcmp(&(blocks[n].fDataType.fID), &ddltype, sizeof(ddltype)) == 0)
		{
			badBlock[n] = not CheckRawDataBlock(event, n, blocks + n);
		}
		else if (not fIgnoreType)
		{
			HLTError("Received raw data block %d in event %lld that we do not know how to handle."
				" The data block has data type '%s' and specification 0x%8.8X.",
				n, event, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fSpecification
			);
		}
	}
	
	// Forward the different blocks.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; ++n)
	{
		if (badBlock[n] and fForwardBadBlocks)
		{
			//int result = Forward(blocks + n);
			int result = Forward(blocks + n);
			if (result != 0)
			{
				size = 0;
				return result;
			}
		}
		if (not badBlock[n] and fForwardGoodBlocks)
		{
			int result = Forward(blocks + n);
			if (result != 0)
			{
				size = 0;
				return result;
			}
		}
	}
	
	size = 0;
	return 0;
}


bool AliHLTTPCDataCheckerComponent::CheckRawDataBlock(
		AliHLTEventID_t event, AliHLTUInt32_t index,
		const AliHLTComponentBlockData* block
	)
{
	// Checks TPC raw DDL data blocks.
	
	assert(fRawReader != NULL);
	assert(fRawStream != NULL);
	
	// Check the origin field of the data block.
	if (not fIgnoreOrigin and
	    memcmp(&(block->fDataType.fOrigin), &kAliHLTDataOriginTPC, sizeof(kAliHLTDataOriginTPC)) != 0
	   )
	{
		char origin[kAliHLTComponentDataTypefOriginSize+1];
		char expectedOrigin[kAliHLTComponentDataTypefOriginSize+1];
		memcpy(&origin, &(block->fDataType.fOrigin), kAliHLTComponentDataTypefOriginSize);
		memcpy(&expectedOrigin, &(kAliHLTDataOriginTPC), kAliHLTComponentDataTypefOriginSize);
		origin[kAliHLTComponentDataTypefOriginSize] = '\0'; // remember the NULL character for the ANSI string.
		expectedOrigin[kAliHLTComponentDataTypefOriginSize] = '\0';
		HLTError("Received raw DDL data block %d in event %lld which has an origin '%s', but expected '%s'.",
			index, event, origin, expectedOrigin
		);
		return false;
	}
	
	// Decode and check the specification bits.
	AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(block->fSpecification);
	AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr(block->fSpecification);
	Int_t ddlid = AliHLTTPCDefinitions::SlicePatchToDDLId(slice, patch);
	if (slice != AliHLTTPCDefinitions::GetMaxSliceNr(block->fSpecification))
	{
		HLTError("Received raw DDL data block %d in event %lld which has a"
			" different minimum and maximum slice number (%d verse %d).",
			index, event, int(slice), int(AliHLTTPCDefinitions::GetMaxSliceNr(block->fSpecification))
		);
		return false;
	}
	if (patch != AliHLTTPCDefinitions::GetMaxPatchNr(block->fSpecification))
	{
		HLTError("Received raw DDL data block %d in event %lld which has a"
			" different minimum and maximum patch number (%d verse %d).",
			index, event, int(patch), int(AliHLTTPCDefinitions::GetMaxPatchNr(block->fSpecification))
		);
		return false;
	}
	if (ddlid == -1)
	{
		HLTError("Received raw DDL data block %d in event %lld which has an"
			" invalid specification 0x%8.8X. Cannot decode the DDL ID number.",
			index, event, block->fSpecification
		);
		return false;
	}
	
	// Now try decode the DDL data. Do it in a try catch block in case
	// the decoder segfaults. We want to know about this and log it.
	bool result = false;
	try
	{
		fRawReader->ClearBuffers();
		fRawReader->Reset();
		fRawStream->Reset();
		fRawReader->AddBuffer(reinterpret_cast<UChar_t*>(block->fPtr), block->fSize, ddlid);
		if (fRawStream->NextDDL())
		{
			while (fRawStream->NextChannel())
			{
				while (fRawStream->NextBunch())
				{
					Int_t bunchLength = fRawStream->GetBunchLength();
					const UShort_t* bunchData = fRawStream->GetSignals();
					for (Int_t i = 0; i < bunchLength; ++i)
					{
						// Check that the 10 bit signal is within range.
						if (bunchData[i] >= 1024)
						{
							HLTWarning("Data signal %d in sector %d row %d pad %d is a strange value %d,"
								" for data block %d (DDL ID = %d) in event %lld.",
								i, fRawStream->GetSector(), fRawStream->GetRow(),
								fRawStream->GetPad(), bunchData[i], index, ddlid, event
							);
						}
					}
				}
			}
			result = true;
		}
		else
		{
			HLTError("Cannot decode the raw DDL data (DDL ID = %d) from"
				" data block %d in event %lld.",
				ddlid, index, event
			);
		}
	}
	catch (...)
	{
		HLTFatal("Caught an exception when processing raw DDL data (DDL ID = %d)"
			" from data block %d in event %lld.",
			ddlid, index, event
		);
		throw;
	}
	return result;
}

