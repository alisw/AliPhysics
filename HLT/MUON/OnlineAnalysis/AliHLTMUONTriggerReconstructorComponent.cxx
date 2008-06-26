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
	fUseCrateId(true)
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


void AliHLTMUONTriggerReconstructorComponent::GetInputDataTypes( std::vector<AliHLTComponentDataType>& list)
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
	/// Inherited from AliHLTComponent. Returns the output data type.
	///
	
	return AliHLTMUONConstants::TriggerRecordsBlockDataType();
}


void AliHLTMUONTriggerReconstructorComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	///
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	///
	
	constBase = sizeof(AliHLTMUONTriggerRecordsBlockWriter::HeaderType);
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
	
	// perform initialization.
	
	HLTInfo("Initialising dHLT trigger reconstructor component.");
	
	// Make sure to cleanup fTrigRec if it is still there for some reason.
	if (fTrigRec != NULL)
	{
		delete fTrigRec;
		fTrigRec = NULL;
	}
	
	try
	{
		fTrigRec = new AliHLTMUONTriggerReconstructor();
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the trigger reconstructor component.");
		return -ENOMEM;
	}
	
	fDDL = -1;
	fWarnForUnexpecedBlock = false;
	fStopOnOverflow = false;
	fUseCrateId = true;
	
	const char* lutFileName = NULL;
	const char* cdbPath = NULL;
	Int_t run = -1;
	bool useCDB = false;
	bool suppressPartialTrigs = true;
	bool tryRecover = false;
	
	for (int i = 0; i < argc; i++)
	{
		if (strcmp( argv[i], "-lut" ) == 0)
		{
			if ( argc <= i+1 )
			{
				HLTError("LookupTable filename not specified." );
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
				return -EINVAL;
			}
			
			lutFileName = argv[i+1];
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-ddl" ) == 0)
		{
			if ( argc <= i+1 )
			{
				HLTError("DDL number not specified. It must be in the range [21..22]" );
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
				return -EINVAL;
			}
		
			char* cpErr = NULL;
			unsigned long num = strtoul(argv[i+1], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a DDL Number.", argv[i+1]);
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
				return -EINVAL;
			}
			if (num < 21 or 22 < num)
			{
				HLTError("The DDL number must be in the range [21..22].");
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
				return -EINVAL;
			}
			fDDL = num - 1; // Convert to DDL number in the range 0..21
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-ddlid" ) == 0)
		{
			if ( argc <= i+1 )
			{
				HLTError("DDL equipment ID number not specified. It must be in the range [2816..2817]" );
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
				return -EINVAL;
			}
		
			char* cpErr = NULL;
			unsigned long num = strtoul(argv[i+1], &cpErr, 0);
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a DDL equipment ID Number.", argv[i+1]);
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
				return -EINVAL;
			}
			fDDL = AliHLTMUONUtils::EquipIdToDDLNumber(num); // Convert to DDL number in the range 0..21
			if (fDDL < 20 or 21 < fDDL)
			{
				HLTError("The DDL equipment ID number must be in the range [2816..2817].");
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
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
		
		if (strcmp( argv[i], "-cdbpath" ) == 0)
		{
			if ( argc <= i+1 )
			{
				HLTError("The CDB path was not specified." );
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
				return -EINVAL;
			}
			cdbPath = argv[i+1];
			useCDB = true;
			i++;
			continue;
		}
	
		if (strcmp( argv[i], "-run" ) == 0)
		{
			if ( argc <= i+1 )
			{
				HLTError("The RUN number was not specified." );
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			run = Int_t( strtoul(argv[i+1], &cpErr, 0) );
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to a valid run number."
					" Expected an integer value.", argv[i+1]
				);
				// Make sure to delete fTrigRec to avoid partial initialisation.
				delete fTrigRec;
				fTrigRec = NULL;
				return -EINVAL;
			}
			
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
		
		HLTError("Unknown option '%s'.", argv[i] );
		// Make sure to delete fTrigRec to avoid partial initialisation.
		delete fTrigRec;
		fTrigRec = NULL;
		return -EINVAL;
		
	} // for loop
	
	if (lutFileName == NULL) useCDB = true;
	
	if (fDDL == -1)
	{
		HLTWarning("DDL number not specified. Cannot check if incomming data is valid.");
	}
	
	int result = 0;
	if (useCDB)
	{
		HLTInfo("Loading lookup table information from CDB for DDL %d (ID = %d).",
			fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
		);
		if (fDDL == -1)
			HLTWarning("DDL number not specified. The lookup table loaded from CDB will be empty!");
		result = ReadCDB(cdbPath, run);
	}
	else
	{
		HLTInfo("Loading lookup table information from file %s.", lutFileName);
		result = ReadLookUpTable(lutFileName);
	}
	if (result != 0)
	{
		// Error messages already generated in ReadCDB or ReadLookUpTable.
		
		// Make sure to delete fTrigRec to avoid partial initialisation.
		delete fTrigRec;
		fTrigRec = NULL;
		return result;
	}
	
	fTrigRec->SuppressPartialTriggers(suppressPartialTrigs);
	fTrigRec->TryRecover(tryRecover);
	fTrigRec->UseCrateId(fUseCrateId);
	
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
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		std::vector<AliHLTComponentBlockData>& outputBlocks
	)
{
	///
	/// Inherited from AliHLTProcessor. Processes the new event data.
	///
	
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
			else
				HLTDebug("Received a data block of a type we cannot handle: '%s', spec: 0x%X",
					DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fSpecification
				);
			
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
				 " We require at least %ufTrigRec->GetkDDLHeaderSize() bytes, but have %u bytes left.",
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
	}
	
	// Finally we set the total size of output memory we consumed.
	size = totalSize;
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


int AliHLTMUONTriggerReconstructorComponent::ReadCDB(const char* cdbPath, Int_t run)
{
	/// Loads the lookup table containing channel and geometrical position
	/// information about trigger strips from CDB.
	/// \param cdbPath  This specifies the CDB path to use to load from.
	///                 Can be set to NULL meaning the default storage is used.
	/// \param run  Specifies the run number to use. If set to -1 then the
	///             default / current run number set for the CDB is used.
	/// \return 0 on success and non zero codes for errors.

	if (fDDL == -1)
	{
		HLTError("No DDL number specified for which to load LUT data from CDB.");
		return -EINVAL;
	}

	int result = FetchMappingStores(cdbPath, run);
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
						
						AliMpPad pad = seg->PadByLocation(AliMpIntPair(boardId, bitxy+offset), kFALSE);
					
						if (! pad.IsValid())
						{
							// There is no pad associated with the given local board and bit pattern.
							continue;
						}
						
						// Get the global coodinates of the pad.
						Float_t lx = pad.Position().X();
						Float_t ly = pad.Position().Y();
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
		const char* cdbPath, Int_t run
		//TODO add option fCrateId
	)
{
	/// Generates a binary file containing the lookup table (LUT) from the
	/// CDB, which can be used for the trigger reconstructor component later.
	/// @param ddl  Must be the DDL for which to generate the DDL,
	///             in the range [20..21].
	/// @param filename  The name of the LUT file to generate.
	/// @param cdbPath  The CDB path to use.
	/// @param run  The run number to use for the CDB.
	/// @return  True if the generation of the LUT file succeeded.
	
	AliHLTMUONTriggerReconstructorComponent comp;
	
	if (ddl < 20 or 21 < ddl)
	{
		std::cerr << "ERROR: the DDL number must be in the range [20..21]." << std::endl;
		return false;
	}
	
	char ddlNum[32];
	char runNum[32];
	sprintf(ddlNum, "%d", ddl+1);
	sprintf(runNum, "%d", run);
	const char* argv[7] = {"-ddl", ddlNum, "-cdbpath", cdbPath, "-run", runNum, NULL};
	int result = comp.DoInit(6, argv);
	if (result != 0)
	{
		// Error message already generated in DoInit.
		return false;
	}
	
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
	
	comp.DoDeinit();
	
	return true;
}
