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

// $Id$

///
/// @file   AliHLTMUONTriggerReconstructorComponent.cxx
/// @author Indranil Das <indra.das@saha.ac.in>, Artur Szostak <artursz@iafrica.com>
/// @date   18 Sep 2007
/// @brief  Implementation of the trigger DDL reconstructor component.
///

#include "AliHLTMUONTriggerReconstructorComponent.h"
#include "AliHLTMUONTriggerReconstructor.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONCalculations.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliRawDataHeader.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMpLocalBoard.h"
#include "AliMpTriggerCrate.h"
#include <cstdlib>
#include <cerrno>
#include <cassert>
#include <fstream>


ClassImp(AliHLTMUONTriggerReconstructorComponent)


AliHLTMUONTriggerReconstructorComponent::AliHLTMUONTriggerReconstructorComponent() :
	AliHLTMUONProcessor(),
	fTrigRec(NULL),
	fDDL(-1),
	fWarnForUnexpecedBlock(false),
	fStopOnOverflow(false),
	fUseCrateId(true),
	fZmiddleSpecified(false),
	fBLSpecified(false),
	fLutInitialised(false)
{
	///
	/// Default constructor.
	///
}


AliHLTMUONTriggerReconstructorComponent::~AliHLTMUONTriggerReconstructorComponent()
{
	///
	/// Default destructor.
	///
	
	if (fTrigRec != NULL) delete fTrigRec;
}


const char* AliHLTMUONTriggerReconstructorComponent::GetComponentID()
{
	///
	/// Inherited from AliHLTComponent. Returns the component ID.
	///
	
	return AliHLTMUONConstants::TriggerReconstructorId();
}


void AliHLTMUONTriggerReconstructorComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	///
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	///
	
	list.clear();
	list.push_back( AliHLTMUONConstants::DDLRawDataType() );
}


AliHLTComponentDataType AliHLTMUONTriggerReconstructorComponent::GetOutputDataType()
{
	///
	/// Inherited from AliHLTComponent. Returns kAliHLTMultipleDataType.
	///
	
	return kAliHLTMultipleDataType;
}


int AliHLTMUONTriggerReconstructorComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTComponent. Returns the output data types.
	
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::TriggerRecordsBlockDataType() );
	list.push_back( AliHLTMUONConstants::TrigRecsDebugBlockDataType() );
	return list.size();
}


void AliHLTMUONTriggerReconstructorComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	///
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	///
	
	constBase = sizeof(AliHLTMUONTriggerRecordsBlockWriter::HeaderType) + 1024*1024;
	inputMultiplier = 4;
}


AliHLTComponent* AliHLTMUONTriggerReconstructorComponent::Spawn()
{
	///
	/// Inherited from AliHLTComponent. Creates a new object instance.
	///
	
	return new AliHLTMUONTriggerReconstructorComponent;
}


int AliHLTMUONTriggerReconstructorComponent::DoInit(int argc, const char** argv)
{
	///
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	///
	
	HLTInfo("Initialising dHLT trigger reconstructor component.");
	
	// Inherit the parents functionality.
	int result = AliHLTMUONProcessor::DoInit(argc, argv);
	if (result != 0) return result;
	
	// Make sure to cleanup fTrigRec if it is still there for some reason.
	if (fTrigRec != NULL)
	{
		delete fTrigRec;
		fTrigRec = NULL;
	}
	
	fDDL = -1;
	fWarnForUnexpecedBlock = false;
	fStopOnOverflow = false;
	fUseCrateId = true;
	fZmiddleSpecified = false;
	fBLSpecified = false;
	fLutInitialised = false;
	
	const char* lutFileName = NULL;
	bool useCDB = false;
	bool suppressPartialTrigs = true;
	bool tryRecover = false;
	bool useLocalId = true;
	bool makeDebugInfo = false;
	bool dontPrintWrongEventError = false;
	double zmiddle = 0;
	double bfieldintegral = 0;
	
	for (int i = 0; i < argc; i++)
	{
		// To keep the legacy behaviour we need to have the following check
		// for -cdbpath here, before ArgumentAlreadyHandled.
		if (strcmp(argv[i], "-cdbpath") == 0)
		{
			useCDB = true;
		}

		if (ArgumentAlreadyHandled(i, argv[i])) continue;

		if (strcmp( argv[i], "-lut" ) == 0)
		{
			if (lutFileName != NULL)
			{
				HLTWarning("LUT path was already specified."
					" Will replace previous value given by -lut."
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("The lookup table filename was not specified." );
				return -EINVAL;
			}
			
			lutFileName = argv[i+1];
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-ddl" ) == 0)
		{
			if (fDDL != -1)
			{
				HLTWarning("DDL number was already specified."
					" Will replace previous value given by -ddl or -ddlid."
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("DDL number not specified. It must be in the range [21..22]");
				return -EINVAL;
			}
		
			char* cpErr = NULL;
			unsigned long num = strtoul(argv[i+1], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a DDL Number.", argv[i+1]);
				return -EINVAL;
			}
			if (num < 21 or 22 < num)
			{
				HLTError("The DDL number must be in the range [21..22].");
				return -EINVAL;
			}
			fDDL = num - 1; // Convert to DDL number in the range 0..21
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-ddlid" ) == 0)
		{
			if (fDDL != -1)
			{
				HLTWarning("DDL number was already specified."
					" Will replace previous value given by -ddl or -ddlid."
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("DDL equipment ID number not specified. It must be in the range [2816..2817]");
				return -EINVAL;
			}
		
			char* cpErr = NULL;
			unsigned long num = strtoul(argv[i+1], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a DDL equipment ID Number.", argv[i+1]);
				return -EINVAL;
			}
			fDDL = AliHLTMUONUtils::EquipIdToDDLNumber(num); // Convert to DDL number in the range 0..21
			if (fDDL < 20 or 21 < fDDL)
			{
				HLTError("The DDL equipment ID number must be in the range [2816..2817].");
				return -EINVAL;
			}
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-cdb" ) == 0)
		{
			useCDB = true;
			continue;
		}
	
		if (strcmp( argv[i], "-zmiddle" ) == 0)
		{
			if (fZmiddleSpecified)
			{
				HLTWarning("The Z coordinate for the middle of the dipole was already specified."
					" Will replace previous value given by -zmiddle."
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("The Z coordinate for the middle of the dipole was not specified.");
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			zmiddle = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid floating point number.",
					argv[i+1]
				);
				return -EINVAL;
			}
			
			fZmiddleSpecified = true;
			i++;
			continue;
		}
	
		if (strcmp( argv[i], "-bfieldintegral" ) == 0)
		{
			if (fBLSpecified)
			{
				HLTWarning("The magnetic field integral was already specified."
					" Will replace previous value given by -bfieldintegral."
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("The magnetic field integral was not specified." );
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			bfieldintegral = strtod(argv[i+1], &cpErr);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid floating point number.",
					argv[i+1]
				);
				return -EINVAL;
			}
			
			fBLSpecified = true;
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-warn_on_unexpected_block" ) == 0)
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}
		
		if (strcmp( argv[i], "-suppress_partial_triggers" ) == 0)
		{
			suppressPartialTrigs = true;
			continue;
		}
		
		if (strcmp( argv[i], "-generate_partial_triggers" ) == 0)
		{
			suppressPartialTrigs = false;
			continue;
		}
		
		if (strcmp( argv[i], "-stop_on_buffer_overflow" ) == 0)
		{
			fStopOnOverflow = true;
			continue;
		}
		
		if (strcmp( argv[i], "-tryrecover" ) == 0)
		{
			tryRecover = true;
			continue;
		}
		
		if (strcmp( argv[i], "-dont_use_crateid" ) == 0)
		{
			fUseCrateId = false;
			continue;
		}
		
		if (strcmp( argv[i], "-dont_use_localid" ) == 0)
		{
			useLocalId = false;
			continue;
		}
		
		if (strcmp( argv[i], "-makedebuginfo" ) == 0)
		{
			makeDebugInfo = true;
			continue;
		}
		
		if (strcmp( argv[i], "-dontprintwrongeventerror" ) == 0)
		{
			dontPrintWrongEventError = true;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
		
	} // for loop
	
	try
	{
		fTrigRec = new AliHLTMUONTriggerReconstructor();
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the trigger reconstructor component.");
		return -ENOMEM;
	}
	
	if (fZmiddleSpecified and useCDB)
	{
		HLTWarning("The -cdb or -cdbpath parameter was specified, which indicates that"
			" this component should read from the CDB, but then the -zmiddle argument"
			" was also used. Will override the value from CDB with the command"
			" line parameter given."
		);
	}
	if (fBLSpecified and useCDB)
	{
		HLTWarning("The -cdb or -cdbpath parameter was specified, which indicates that"
			" this component should read from the CDB, but then the -bfieldintegral"
			" argument was also used. Will override the value from CDB with the"
			" command line parameter given."
		);
	}
	
	if (lutFileName != NULL and useCDB == true)
	{
		HLTWarning("The -cdb or -cdbpath parameter was specified, which indicates that"
			" this component should read from the CDB, but then the -lut argument"
			" was also used. Will ignore the -lut option and load from CDB anyway."
		);
	}
	
	if (lutFileName == NULL) useCDB = true;
	
	if (fDDL == -1 and not DelaySetup())
	{
		HLTWarning("DDL number not specified. Cannot check if incomming data is valid.");
	}
	
	if (useCDB)
	{
		if (not DelaySetup())
		{
			HLTInfo("Loading lookup table information from CDB for DDL %d (ID = %d).",
				fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
			);
			result = ReadLutFromCDB();
			if (result != 0)
			{
				// Error messages already generated in ReadLutFromCDB.
				delete fTrigRec; // Make sure to delete fTrigRec to avoid partial initialisation.
				fTrigRec = NULL;
				return result;
			}
			fLutInitialised = true;
		}
	}
	else
	{
		HLTInfo("Loading lookup table information from file %s.", lutFileName);
		result = ReadLookUpTable(lutFileName);
		if (result != 0)
		{
			// Error messages already generated in ReadLookUpTable.
			delete fTrigRec; // Make sure to delete fTrigRec to avoid partial initialisation.
			fTrigRec = NULL;
			return result;
		}
		fLutInitialised = true;
	}
	
	if (fZmiddleSpecified) AliHLTMUONCalculations::Zf(zmiddle);
	if (fBLSpecified) AliHLTMUONCalculations::QBL(bfieldintegral);
	
	if (not DelaySetup())
	{
		if (not fZmiddleSpecified or not fBLSpecified)
		{
			HLTInfo("Loading configuration parameters from CDB for DDL %d (ID = %d).",
				fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
			);
			
			result = ReadConfigFromCDB(not fZmiddleSpecified, not fBLSpecified);
			if (result != 0)
			{
				// Error messages already generated in ReadConfigFromCDB.
				delete fTrigRec; // Make sure to delete fTrigRec to avoid partial initialisation.
				fTrigRec = NULL;
				return result;
			}
		}
		else
		{
			// Print the debug messages here since ReadConfigFromCDB does not get called,
			// in-which the debug messages would have been printed.
			HLTDebug("Using the following configuration parameters:");
			HLTDebug("  Middle of dipole Z coordinate = %f cm", AliHLTMUONCalculations::Zf());
			HLTDebug("        Magnetic field integral = %f T.m", AliHLTMUONCalculations::QBL());
		}
	}
	
	fTrigRec->SuppressPartialTriggers(suppressPartialTrigs);
	fTrigRec->TryRecover(tryRecover);
	fTrigRec->UseCrateId(fUseCrateId);
	fTrigRec->UseLocalId(useLocalId);
	fTrigRec->StoreDebugInfo(makeDebugInfo);
	fTrigRec->DontPrintWrongEventError(dontPrintWrongEventError);
	
	return 0;
}


int AliHLTMUONTriggerReconstructorComponent::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	
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
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks
	)
{
	///
	/// Inherited from AliHLTProcessor. Processes the new event data.
	///
	
	// Initialise the LUT and configuration parameters from CDB if we were
	// requested to initialise only when the first event was received.
	if (DelaySetup())
	{
		// Use the specification given by the first data block if we
		// have not been given a DDL number on the command line.
		if (fDDL == -1)
		{
			bool blockFound = false;
			for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt and not blockFound; n++)
			{
				if (blocks[n].fDataType != AliHLTMUONConstants::DDLRawDataType()) continue;
				blockFound = true;

				fDDL = AliHLTMUONUtils::SpecToDDLNumber(blocks[n].fSpecification);
			
				if (fDDL == -1)
				{
					HLTError("Received a data block with a specification (0x%8.8X)"
						" indicating multiple DDL data sources, but we must only"
						" receive raw DDL data from one trigger station DDL.",
						blocks[n].fSpecification
					);
					return -EPROTO;
				}
			}

			if (not blockFound)
			{
				HLTError("The initialisation from CDB of the component has"
					" been delayed to the first received event. However,"
					" no raw DDL data blocks have been found in the first event."
				);
				return -ENOENT;
			}
		}
		
		// Check that the LUT was not already loaded in DoInit.
		if (not fLutInitialised)
		{
			HLTInfo("Loading lookup table information from CDB for DDL %d (ID = %d).",
				fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
			);
			int result = ReadLutFromCDB();
			if (result != 0) return result;
			fLutInitialised = true;
		}
		
		// Load the configuration paramters from CDB if they have not been given
		// on the command line.
		if (not fZmiddleSpecified or not fBLSpecified)
		{
			HLTInfo("Loading configuration parameters from CDB for DDL %d (ID = %d).",
				fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
			);
			int result = ReadConfigFromCDB(not fZmiddleSpecified, not fBLSpecified);
			if (result != 0) return result;
		}
		
		DoneDelayedSetup();
	}
	
	// Process an event
	unsigned long totalSize = 0; // Amount of memory currently consumed in bytes.

	HLTDebug("Processing event %llu with %u input data blocks.",
		evtData.fEventID, evtData.fBlockCnt
	);
	
	// Loop over all input blocks in the event and run the trigger DDL
	// reconstruction algorithm on the raw data.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
			n, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fPtr, blocks[n].fSize
		);

		if (blocks[n].fDataType != AliHLTMUONConstants::DDLRawDataType()
		    or not AliHLTMUONUtils::IsTriggerDDL(blocks[n].fSpecification)
		   )
		{
			// Log a message indicating that we got a data block that we
			// do not know how to handle.
			if (fWarnForUnexpecedBlock)
				HLTWarning("Received a data block of a type we cannot handle: '%s', spec: 0x%X",
					DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fSpecification
				);
#ifdef __DEBUG
			else
				HLTDebug("Received a data block of a type we cannot handle: '%s', spec: 0x%X",
					DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fSpecification
				);
#endif
			
			continue;
		}
		
		AliHLTInt32_t receivedDDL = AliHLTMUONUtils::SpecToDDLNumber(blocks[n].fSpecification);
		if (fDDL != -1)
		{
			if (receivedDDL != fDDL)
			{
				HLTWarning("Received raw data from DDL %d (ID = %d),"
					" but expect data only from DDL %d (ID = %d).",
					receivedDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(receivedDDL),
					fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
				);
			}
		}
		
		// Create a new output data block and initialise the header.
		AliHLTMUONTriggerRecordsBlockWriter block(outputPtr+totalSize, size-totalSize);
		if (not block.InitCommonHeader())
		{
			HLTError("There is not enough space in the output buffer for the new data block."
				 " We require at least %u bytes, but have %u bytes left.",
				sizeof(AliHLTMUONTriggerRecordsBlockWriter::HeaderType),
				block.BufferSize()
			);
			break;
		}

		AliHLTUInt32_t totalDDLSize = blocks[n].fSize;
		if (totalDDLSize < sizeof(AliRawDataHeader))
		{
			HLTError("Raw data block %d is %d bytes in size and is too short to"
				 " possibly contain valid DDL raw data. We expect it to have"
				 " at least %d bytes for the commond data header.",
				n, totalDDLSize, sizeof(AliRawDataHeader)
			);
			continue;
		}
		AliRawDataHeader* header = reinterpret_cast<AliRawDataHeader*>(blocks[n].fPtr);
		AliHLTUInt32_t payloadSize = totalDDLSize - sizeof(AliRawDataHeader);
		AliHLTUInt8_t* buffer = reinterpret_cast<AliHLTUInt8_t*>(header + 1);
		AliHLTUInt32_t nofTrigRec = block.MaxNumberOfEntries();
		
		// Decode if this is a scalar event or not.
		bool scalarEvent = ((header->GetL1TriggerMessage() & 0x1) == 0x1);
		
		// Remember: the following does NOT change the mapping!
		// It is just to generate unique trigger record IDs.
		fTrigRec->SetDDL(receivedDDL);
		
		bool runOk = fTrigRec->Run(
				buffer, payloadSize, scalarEvent,
				block.GetArray(), nofTrigRec
			);
		if (not runOk)
		{
			HLTError("Error while processing the trigger DDL reconstruction algorithm.");
			if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
			if (not fTrigRec->OverflowedOutputBuffer()
			    or (fTrigRec->OverflowedOutputBuffer() and fStopOnOverflow)
			   )
			{
				size = totalSize; // Must tell the framework how much buffer space was used.
				return -EIO;
			}
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
		
		if (fTrigRec->StoreDebugInfo())
		{
			// Create a new output data block and initialise the header.
			AliHLTMUONTrigRecsDebugBlockWriter infoblock(outputPtr+totalSize, size-totalSize);
			if (not infoblock.InitCommonHeader())
			{
				HLTError("There is not enough space in the output buffer for the new debug"
					" data block. We require at least %u bytes, but have %u bytes left.",
					sizeof(AliHLTMUONTrigRecsDebugBlockWriter::HeaderType),
					infoblock.BufferSize()
				);
				break;
			}
			
			infoblock.SetNumberOfEntries(fTrigRec->InfoBufferCount());
			for (AliHLTUInt32_t i = 0; i < fTrigRec->InfoBufferCount(); ++i)
			{
				infoblock[i] = fTrigRec->InfoBuffer()[i];
			}
			fTrigRec->ZeroInfoBuffer();
			
			// Fill the block data structure for our output block.
			AliHLTComponentBlockData bd2;
			FillBlockData(bd2);
			bd2.fPtr = outputPtr;
			// This block's start (offset) is after all other blocks written so far.
			bd2.fOffset = totalSize;
			bd2.fSize = infoblock.BytesUsed();
			bd2.fDataType = AliHLTMUONConstants::TrigRecsDebugBlockDataType();
			bd2.fSpecification = blocks[n].fSpecification;
			outputBlocks.push_back(bd2);
			
			HLTDebug("Created a new output data block for debug information at fPtr = %p,"
				" with fOffset = %u (0x%.X) and fSize = %u bytes.",
				bd2.fPtr, bd2.fOffset, bd2.fOffset, bd2.fSize
			);
			
			// Increase the total amount of data written so far to our output memory.
			totalSize += infoblock.BytesUsed();
		}
	}
	
	// Finally we set the total size of output memory we consumed.
	size = totalSize;
	return 0;
}


int AliHLTMUONTriggerReconstructorComponent::Reconfigure(
		const char* cdbEntry, const char* componentId
	)
{
	/// Inherited from AliHLTComponent. This method will reload CDB configuration
	/// entries for this component from the CDB.
	/// \param cdbEntry If this is NULL then it is assumed that all CDB entries should
	///      be reloaded. Otherwise a particular value for 'cdbEntry' will trigger
	///      reloading of the LUT if the path contains 'MUON/', but the other
	///      configuration parameters will be loaded if 'cdbEntry' contains
	///      "HLT/ConfigMUON/TriggerReconstructor".
	/// \param componentId  The name of the component in the current chain.
	
	TString path = cdbEntry;
	bool startsWithMUON = path.Index("MUON/", 5, 0, TString::kExact) == 0;
	bool givenConfigPath = (path == AliHLTMUONConstants::TriggerReconstructorCDBPath());
	
	if (cdbEntry == NULL or startsWithMUON or givenConfigPath)
	{
		HLTInfo("Reading new configuration entries from CDB for component '%s'.", componentId);
	}
	
	if (cdbEntry == NULL or startsWithMUON)
	{
		int result = ReadLutFromCDB();
		if (result != 0) return result;
	}
	
	if (cdbEntry == NULL or givenConfigPath)
	{
		int result = ReadConfigFromCDB();
		if (result != 0) return result;
	}
	
	return 0;
}


int AliHLTMUONTriggerReconstructorComponent::ReadPreprocessorValues(const char* modules)
{
	/// Inherited from AliHLTComponent. 
	/// Updates the configuration of this component if either HLT or MUON have
	/// been specified in the 'modules' list.

	TString mods = modules;
	if (mods.Contains("ALL") or (mods.Contains("HLT") and mods.Contains("MUON")))
	{
		return Reconfigure(NULL, GetComponentID());
	}
	if (mods.Contains("HLT"))
	{
		return Reconfigure(AliHLTMUONConstants::TriggerReconstructorCDBPath(), GetComponentID());
	}
	if (mods.Contains("MUON"))
	{
		return Reconfigure("MUON/*", GetComponentID());
	}
	return 0;
}


int AliHLTMUONTriggerReconstructorComponent::ReadLookUpTable(const char* lutpath)
{
	///
	/// Read in the lookup table from file.
	///
	
	assert(fTrigRec != NULL);

	fstream file;
	file.open(lutpath, fstream::binary | fstream::in);
	if (not file)
	{
		HLTError("Could not open file: %s", lutpath);
		return -ENOENT;
	}
	
	file.read(reinterpret_cast<char*>(fTrigRec->LookupTableBuffer()), fTrigRec->LookupTableSize());
	if (file.eof())
	{
		HLTError("The file %s was too short to contain a valid lookup table for this component.", lutpath);
		file.close();
		return -EIO;
	}
	if (file.fail())
	{
		HLTError("Could not read from file: %s", lutpath);
		file.close();
		return -EIO;
	}
	
	file.close();
	return 0;
}


int AliHLTMUONTriggerReconstructorComponent::ReadConfigFromCDB(
		bool setZmiddle, bool setBL
	)
{
	/// Reads this component's configuration parameters from the CDB.
	/// These include the middle of the dipole Z coordinate (zmiddle) and the
	/// integrated magnetic field of the dipole.
	/// \param setZmiddle Indicates if the zmiddle parameter should be set
	///       (default true).
	/// \param setBL Indicates if the integrated magnetic field parameter should
	///       be set (default true).
	/// \return 0 if no errors occured and negative error code compatible with
	///       the HLT framework on errors.

	const char* pathToEntry = AliHLTMUONConstants::TriggerReconstructorCDBPath();
	
	TMap* map = NULL;
	int result = FetchTMapFromCDB(pathToEntry, map);
	if (result != 0) return result;
	
	Double_t value = 0;
	if (setZmiddle)
	{
		result = GetFloatFromTMap(map, "zmiddle", value, pathToEntry, "dipole middle Z coordinate");
		if (result != 0) return result;
		AliHLTMUONCalculations::Zf(value);
	}
	
	if (setBL)
	{
		Double_t bfieldintegral;
		result = FetchFieldIntegral(bfieldintegral);
		if (result == 0)
		{
			AliHLTMUONCalculations::QBL(bfieldintegral);
		}
		else
		{
			HLTWarning("Failed to load the magnetic field integral from GRP information.");
			result = GetFloatFromTMap(map, "bfieldintegral", value, pathToEntry, "integrated magnetic field");
			if (result != 0) return result;
			HLTWarning(Form("Using deprecated magnetic field integral value of %f T.m.", value));
			AliHLTMUONCalculations::QBL(value);
		}
	}
	
	HLTDebug("Using the following configuration parameters:");
	HLTDebug("  Middle of dipole Z coordinate = %f cm", AliHLTMUONCalculations::Zf());
	HLTDebug("        Magnetic field integral = %f T.m", AliHLTMUONCalculations::QBL());
	
	return 0;
}


int AliHLTMUONTriggerReconstructorComponent::ReadLutFromCDB()
{
	/// Loads the lookup table containing channel and geometrical position
	/// information about trigger strips from CDB.
	///
	/// \note To override the default CDB path and / or run number the
	/// SetCDBPathAndRunNo(cdbPath, run) method should be called before this
	/// method.
	///
	/// \return 0 on success and non zero codes for errors.

	if (fDDL == -1)
	{
		HLTError("No DDL number specified for which to load LUT data from CDB.");
		return -EINVAL;
	}

	int result = FetchMappingStores();
	// Error message already generated in FetchMappingStores.
	if (result != 0) return result;
	AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
	
	AliMpSegmentation* segmentation = AliMpSegmentation::Instance();
	if (segmentation == NULL)
	{
		HLTError("Could not find segmentation mapping (AliMpSegmentation) instance.");
		return -EIO;
	}
	
	// Only load geometry if not already loaded.
	if (AliGeomManager::GetGeometry() == NULL)
	{
		AliGeomManager::LoadGeometry();
	}
	AliMUONGeometryTransformer transformer;
	if (not transformer.LoadGeometryData())
	{
		HLTError("Could not load geometry into transformer.");
		return -ENOENT;
	}
	
	AliHLTMUONTriggerRecoLookupTable* lookupTable = fTrigRec->LookupTableBuffer();
	
	for (Int_t i = 0; i < 16; i++)
	for (Int_t j = 0; j < 16; j++)
	for (Int_t k = 0; k < 4; k++)
	for (Int_t n = 0; n < 2; n++)
	for (Int_t m = 0; m < 16; m++)
	{
		lookupTable->fRow[i][j][k][n][m].fIdFlags = 0x0;
		lookupTable->fRow[i][j][k][n][m].fX = 0;
		lookupTable->fRow[i][j][k][n][m].fY = 0;
		lookupTable->fRow[i][j][k][n][m].fZ = 0;
	}
	
	AliMpDEIterator detElemIter;
	for (Int_t iReg = 0; iReg < 8; iReg++)
	{
		AliMpTriggerCrate* crate = ddlStore->GetTriggerCrate(fDDL, iReg);
		if (crate == NULL)
		{
			HLTError("Could not get crate mapping for regional header = %d"
				" and DDL %d (ID = %d).",
				iReg, fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
			);
			continue;
		}
		// Depending on the value of fUseCrateId, use either the crate ID as would
		// be found in the regional header structures or the sequencial index number
		// of the structure.
		UInt_t crateId = (fUseCrateId ? crate->GetId() : iReg);
		if (crateId >= 16)
		{
			HLTError("The crate ID number (%d) for regional header = %d and"
				" DDL %d (ID = %d) is too big. It should be in the range [0..15]",
				crateId, iReg, fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
			);
			continue;
		}
		
		for (Int_t iLocBoard = 0; iLocBoard < 16; iLocBoard++)
		{
			Int_t boardId = crate->GetLocalBoardId(iLocBoard);
			if (boardId == 0) continue;
			
			AliMpLocalBoard* localBoard = ddlStore->GetLocalBoard(boardId);
			if (localBoard == NULL)
			{
				HLTError("Could not get local board: %d.", boardId);
				continue;
			}

			// skip copy cards
			if (! localBoard->IsNotified()) continue;
		
			for (Int_t iChamber = 0; iChamber < 4; iChamber++)
			{
				Int_t detElemId = ddlStore->GetDEfromLocalBoard(boardId, iChamber);
				
				AliHLTUInt32_t idflags = AliHLTMUONUtils::PackRecHitFlags(iChamber+10, detElemId);
				
				const AliMUONGeometryDetElement* detElemTransform = transformer.GetDetElement(detElemId);
				if (detElemTransform == NULL)
				{
					HLTError("Got NULL pointer for geometry transformer"
						" for detection element ID = %d.",
						detElemId
					);
					continue;
				}
				
				for (Int_t iCathode = 0; iCathode <= 1; iCathode++)
				{
					const AliMpVSegmentation* seg = segmentation->GetMpSegmentation(
							detElemId, AliMp::GetCathodType(iCathode)
						);
					
					for (Int_t bitxy = 0; bitxy < 16; bitxy++)
					{
						Int_t offset = 0;
						if (iCathode && localBoard->GetSwitch(6)) offset = -8;
						
// just use this switch for simplicity for two different but shortly after each other applied changes
#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
						AliMpPad pad = seg->PadByLocation(boardId, bitxy+offset, kFALSE);
#else // old AliMpPad functionality < r 31742
						AliMpPad pad = seg->PadByLocation(AliMpIntPair(boardId, bitxy+offset), kFALSE);
#endif //HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					
						if (! pad.IsValid())
						{
							// There is no pad associated with the given local board and bit pattern.
							continue;
						}
						
						// Get the global coodinates of the pad.
#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
						Float_t lx = pad.GetPositionX();
						Float_t ly = pad.GetPositionY();
#else // old AliMpPad functionality < r 31769
						Float_t lx = pad.Position().X();
						Float_t ly = pad.Position().Y();
#endif //HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
						Float_t gx, gy, gz;
						detElemTransform->Local2Global(lx, ly, 0, gx, gy, gz);
						
						// Fill the LUT
						lookupTable->fRow[crateId][iLocBoard][iChamber][iCathode][bitxy].fIdFlags = idflags;
						lookupTable->fRow[crateId][iLocBoard][iChamber][iCathode][bitxy].fX = gx;
						lookupTable->fRow[crateId][iLocBoard][iChamber][iCathode][bitxy].fY = gy;
						lookupTable->fRow[crateId][iLocBoard][iChamber][iCathode][bitxy].fZ = gz;
					}
				}
			}
		}
	}
	
	return 0;
}


bool AliHLTMUONTriggerReconstructorComponent::GenerateLookupTable(
		AliHLTInt32_t ddl, const char* filename,
		const char* cdbPath, Int_t run, bool useCrateId
	)
{
	/// Generates a binary file containing the lookup table (LUT) from the
	/// CDB, which can be used for the trigger reconstructor component later.
	/// @param ddl  Must be the DDL for which to generate the DDL,
	///             in the range [20..21].
	/// @param filename  The name of the LUT file to generate.
	/// @param cdbPath  The CDB path to use.
	/// @param run  The run number to use for the CDB.
	/// @param useCrateId  Indicates that the crate ID should be used rather
	///             than a sequencial number.
	/// @return  True if the generation of the LUT file succeeded.
	
	AliHLTMUONTriggerReconstructorComponent comp;
	try
	{
		// Perform minimal initialisation to be able to fill the LUT buffer.
		comp.fTrigRec = new AliHLTMUONTriggerReconstructor();
	}
	catch (const std::bad_alloc&)
	{
		std::cerr << "ERROR: Could not allocate more memory for the trigger reconstructor component." << std::endl;
		return false;
	}
	
	if (ddl < 20 or 21 < ddl)
	{
		std::cerr << "ERROR: the DDL number must be in the range [20..21]." << std::endl;
		return false;
	}
	
	comp.fDDL = ddl;
	comp.fUseCrateId = useCrateId;
	if (comp.SetCDBPathAndRunNo(cdbPath, run) != 0) return false;
	if (comp.ReadLutFromCDB() != 0) return false;
	
	std::fstream file(filename, std::ios::out);
	if (not file)
	{
		std::cerr << "ERROR: could not open file: " << filename << std::endl;
		return false;
	}
	
	file.write(
			reinterpret_cast<char*>(comp.fTrigRec->LookupTableBuffer()),
			comp.fTrigRec->LookupTableSize()
		);
	if (not file)
	{
		std::cerr << "ERROR: There was a problem writing to the file: " << filename << std::endl;
		return false;
	}
	file.close();
	
	return true;
}
