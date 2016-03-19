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
///  @file   AliHLTMUONHitReconstructorComponent.cxx
///  @author Indranil Das <indra.das@saha.ac.in> | <indra.ehep@gmail.com>, Artur Szostak <artursz@iafrica.com>
///  @date   28 May 2007
///  @brief  Implementation of the hit Reconstruction processing component.
///
///  The HitRec Component is designed to deal the rawdata inputfiles to findout the 
///  the reconstructed hits. The output is send to the output block for further 
///  processing.
///
///  Author : Indranil Das ( indra.das@saha.ac.in || indra.ehep@gmail.com )
///

#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONHitReconstructorComponent.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliHLTMUONHitReconstructor.h"
#include "AliHLTLogging.h"
#include "AliHLTSystem.h"
#include "AliHLTDefinitions.h"
#include "AliHLTCDHWrapper.h"
#include <cstdlib>
#include <cerrno>
#include <cassert>
#include <fstream>

#include "TMap.h"

//STEER
#include "AliCDBManager.h"
#include "AliGeomManager.h"

//MUON
#include "AliMUONGeometryTransformer.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"

//MUON/mapping
#include "AliMpCDB.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"

ClassImp(AliHLTMUONHitReconstructorComponent)


AliHLTMUONHitReconstructorComponent::AliHLTMUONHitReconstructorComponent() :
	AliHLTMUONProcessor(),
	fHitRec(NULL),
	fDDL(-1),
	fLutSize(0),
	fLut(NULL),
	fIdToEntry(),
	fMaxEntryPerBusPatch(),
	fWarnForUnexpecedBlock(false),
	fWarnIfPadSkipped(false)
{
	///
	/// Default constructor.
	///
}


AliHLTMUONHitReconstructorComponent::~AliHLTMUONHitReconstructorComponent()
{
	///
	/// Default destructor.
	///
	
	if (fHitRec != NULL)
	{
		delete fHitRec;
	}
	if (fLut != NULL)
	{
		delete [] fLut;
	}
}

const char* AliHLTMUONHitReconstructorComponent::GetComponentID()
{
	///
	/// Inherited from AliHLTComponent. Returns the component ID.
	///
	
	return AliHLTMUONConstants::HitReconstructorId();
}


void AliHLTMUONHitReconstructorComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	///
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	///
	
	list.clear();
	list.push_back( AliHLTMUONConstants::DDLRawDataType() );
}


AliHLTComponentDataType AliHLTMUONHitReconstructorComponent::GetOutputDataType()
{
	///
	/// Inherited from AliHLTComponent. Returns the output data type.
	///
	
	return kAliHLTMultipleDataType;
}


int AliHLTMUONHitReconstructorComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTComponent. Returns the output data types.
	
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::RecHitsBlockDataType() );
	list.push_back( AliHLTMUONConstants::ClusterBlockDataType() );
	list.push_back( AliHLTMUONConstants::ChannelBlockDataType() );
	return list.size();
}


void AliHLTMUONHitReconstructorComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
	///
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	///
	
	constBase = sizeof(AliHLTMUONRecHitsBlockWriter::HeaderType) + 1024*1024;
	inputMultiplier = 1;
}


AliHLTComponent* AliHLTMUONHitReconstructorComponent::Spawn()
{
	///
	/// Inherited from AliHLTComponent. Creates a new object instance.
	///
	
	return new AliHLTMUONHitReconstructorComponent;
}


int AliHLTMUONHitReconstructorComponent::DoInit(int argc, const char** argv)
{
	///
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	///

	HLTInfo("Initialising dHLT hit reconstruction component.");

	// Inherit the parents functionality.
	int result = AliHLTMUONProcessor::DoInit(argc, argv);
	if (result != 0) return result;
	
	// Must make sure that fHitRec and fLut is deleted if it is still
	// allocated for whatever reason.
	FreeMemory();
	
	// Initialise fields with default values then parse the command line.
	fDDL = -1;
	fIdToEntry.clear();
	fMaxEntryPerBusPatch.clear();
	fWarnForUnexpecedBlock = false;
	fWarnIfPadSkipped = false;
	const char* lutFileName = NULL;
	bool useCDB = false;
	typedef AliHLTMUONHitReconstructor HR;
	HR::ERecoveryMode recoveryMode = HR::kDontTryRecover;
	AliHLTInt32_t dccut = -1;
	bool skipParityErrors = false;
	bool dontPrintParityErrors = false;
	bool makeClusters = false;
	bool makeChannels = false;
	
	for (int i = 0; i < argc; i++)
	{
		// To keep the legacy behaviour we need to have the following check
		// for -cdbpath here, before ArgumentAlreadyHandled.
		if (strcmp(argv[i], "-cdbpath") == 0)
		{
			useCDB = true;
		}

		if (ArgumentAlreadyHandled(i, argv[i])) continue;
		
		if (strcmp( argv[i], "-ddl" ) == 0)
		{
			if (fDDL != -1)
			{
				HLTWarning("DDL number was already specified."
					" Will replace previous value given by -ddl or -ddlid."
				);
			}
			
			if (argc <= i+1)
			{
				HLTError("The DDL number was not specified. Must be in the range [1..20].");
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			unsigned long num = strtoul( argv[i+1], &cpErr, 0 );
			if (cpErr == NULL or *cpErr != '\0')
			{
				HLTError("Cannot convert '%s' to DDL a number.", argv[i+1] );
				return -EINVAL;
			}
			if (num < 1 or 20 < num)
			{
				HLTError("The DDL number must be in the range [1..20].");
				return -EINVAL;
			}
			fDDL = num - 1;  // convert to range [0..19]
			
			i++;
			continue;
		} // -ddl argument
		
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
				HLTError("DDL equipment ID number not specified. It must be in the range [2560..2579]" );
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
			if (fDDL < 0 or 19 < fDDL)
			{
				HLTError("The DDL equipment ID number must be in the range [2560..2579].");
				return -EINVAL;
			}
			
			i++;
			continue;
		}
		
		if (strcmp( argv[i], "-lut" ) == 0)
		{
			if (lutFileName != NULL)
			{
				HLTWarning("LUT path was already specified."
					" Will replace previous value given by -lut."
				);
			}
			
			if (argc <= i+1)
			{
				HLTError("The lookup table filename was not specified.");
				return -EINVAL;
			}
			lutFileName = argv[i+1];
			i++;
			continue;
		} // -lut argument
		
		if (strcmp( argv[i], "-cdb" ) == 0)
		{
			useCDB = true;
			continue;
		} // -cdb argument
		
		if (strcmp( argv[i], "-dccut" ) == 0)
		{
			if (dccut != -1)
			{
				HLTWarning("DC cut parameter was already specified."
					" Will replace previous value given by -dccut."
				);
			}
			
			if ( argc <= i+1 )
			{
				HLTError("No DC cut value was specified. It should be a positive integer value." );
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			dccut = AliHLTInt32_t( strtol(argv[i+1], &cpErr, 0) );
			if (cpErr == NULL or *cpErr != '\0' or dccut < 0)
			{
				HLTError("Cannot convert '%s' to a valid DC cut value."
					" Expected a positive integer value.", argv[i+1]
				);
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
		
		if (strcmp( argv[i], "-tryrecover" ) == 0)
		{
			if (argc > i+1)
			{
				// There might be an optional parameter so check
				// if it is a recognised one. If not then assume it
				// is the next argument, so no error message.
				if (strcmp(argv[i+1], "full") == 0)
				{
					recoveryMode = HR::kRecoverFull;
					i++;
				}
				else if (strcmp(argv[i+1], "skip") == 0)
				{
					recoveryMode = HR::kRecoverJustSkip;
					i++;
				}
				else if (strcmp(argv[i+1], "parityerrors") == 0)
				{
					recoveryMode = HR::kRecoverFromParityErrorsOnly;
					i++;
				}
				else
				{
					recoveryMode = HR::kRecoverFull;
				}
			}
			else
			{
				recoveryMode = HR::kRecoverFull;
			}
			continue;
		}
		
		if (strcmp( argv[i], "-skipparityerrors" ) == 0)
		{
			skipParityErrors = true;
			continue;
		}
		
		if (strcmp( argv[i], "-dontprintparityerrors" ) == 0)
		{
			dontPrintParityErrors = true;
			continue;
		}
		
		if (strcmp( argv[i], "-useidealgain" ) == 0)
    {
      HLTWarning("DEPRECATED option -useidealgain");
			continue;
		}
		
		if (strcmp( argv[i], "-makeclusters" ) == 0)
		{
			makeClusters = true;
			continue;
		}
		
		if (strcmp( argv[i], "-makechannels" ) == 0)
		{
			makeChannels = true;
			continue;
		}
		
		if (strcmp( argv[i], "-warnifpadskipped" ) == 0)
		{
			fWarnIfPadSkipped = true;
			continue;
		}
	
		HLTError("Unknown option '%s'", argv[i]);
		return -EINVAL;
	
	} // for loop
	
	try
	{
		fHitRec = new AliHLTMUONHitReconstructor();
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the hit reconstructor component.");
		return -ENOMEM;
	}
	
	if (dccut != -1 and useCDB)
	{
		HLTWarning("The -cdb or -cdbpath parameter was specified, which indicates that"
			" this component should read from the CDB, but then the -dccut argument"
			" was also used. Will override the value from CDB with the command"
			" line DC cut parameter given."
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
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return result;
			}
			fHitRec->SetLookUpTable(fLut, &fIdToEntry, &fMaxEntryPerBusPatch);
		}
	}
	else
	{
		HLTInfo("Loading lookup table information from file %s.", lutFileName);
		result = ReadLookUpTable(lutFileName);
		if (result != 0)
		{
			// Error messages already generated in ReadLookUpTable.
			FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
			return result;
		}
		fHitRec->SetLookUpTable(fLut, &fIdToEntry, &fMaxEntryPerBusPatch);
	}
	
	if (dccut == -1)
	{
		if (not DelaySetup())
		{
			HLTInfo("Loading DC cut parameters from CDB for DDL %d (ID = %d).",
				fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
			);
			result = ReadDCCutFromCDB();
			if (result != 0)
			{
				// Error messages already generated in ReadDCCutFromCDB.
				FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
				return result;
			}
		}
		else
		{
			// Print the debug messages here since ReadDCCutFromCDB does not get called,
			// in-which the debug messages would have been printed.
			HLTDebug("Using DC cut parameter of %d ADC channels.", fHitRec->GetDCCut());
		}
	}
	else
	{
		fHitRec->SetDCCut(dccut);
		HLTDebug("Using DC cut parameter of %d ADC channels.", fHitRec->GetDCCut());
	}
	
	fHitRec->TryRecover(recoveryMode);
	fHitRec->SkipParityErrors(skipParityErrors);
	fHitRec->DontPrintParityErrors(dontPrintParityErrors);
	fHitRec->GenerateClusterInfo(makeClusters);
	fHitRec->GenerateChannelInfo(makeChannels);
	fHitRec->DDLNumber(fDDL);
	//The DDL number has to be set before the following InitDetElemInDDLArray() method
	fHitRec->InitDetElemInDDLArray();
	HLTDebug("dHLT hit reconstruction component is initialized.");
	return 0;
}


int AliHLTMUONHitReconstructorComponent::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	fHitRec->DeInitDetElemInDDLArray();
	HLTInfo("Deinitialising dHLT hit reconstruction component.");
	FreeMemory();
	return 0;
}


int AliHLTMUONHitReconstructorComponent::Reconfigure(
		const char* cdbEntry, const char* componentId
	)
{
	/// Inherited from AliHLTComponent. This method will reload CDB configuration
	/// entries for this component from the CDB.
	/// \param cdbEntry If this is NULL then it is assumed that all CDB entries should
	///      be reloaded. Otherwise a particular value for 'cdbEntry' will trigger
	///      reloading of the LUT if the path contains 'MUON/' and reloading of the DC
	///      cut parameter if 'cdbEntry' equals "HLT/ConfigMUON/HitReconstructor".
	/// \param componentId  The name of the component in the current chain.
	
	TString path = cdbEntry;
	bool startsWithMUON = path.Index("MUON/", 5, 0, TString::kExact) == 0;
	bool givenConfigPath = (path == AliHLTMUONConstants::HitReconstructorCDBPath());
	
	if (cdbEntry == NULL or startsWithMUON or givenConfigPath)
	{
		HLTInfo("Reading new configuration entries from CDB for component '%s'.", componentId);
	}
		
	if (cdbEntry == NULL or startsWithMUON)
	{
		// First clear the current LUT data and then load in the new values.
		if (fLut != NULL)
		{
			delete [] fLut;
			fLut = NULL;
			fLutSize = 0;
		}
		
		fIdToEntry.clear();
		fMaxEntryPerBusPatch.clear();	
		int result = ReadLutFromCDB();
		if (result != 0) return result;
		fHitRec->SetLookUpTable(fLut, &fIdToEntry, &fMaxEntryPerBusPatch);
	}
	
	if (cdbEntry == NULL or not startsWithMUON)
	{
		int result = ReadDCCutFromCDB();
		if (result != 0) return result;
	}
	
	return 0;
}


int AliHLTMUONHitReconstructorComponent::ReadPreprocessorValues(const char* modules)
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
		return Reconfigure(AliHLTMUONConstants::HitReconstructorCDBPath(), GetComponentID());
	}
	if (mods.Contains("MUON"))
	{
		return Reconfigure("MUON/*", GetComponentID());
	}
	return 0;
}


int AliHLTMUONHitReconstructorComponent::DoEvent(
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

	// Initialise the LUT and DC cut parameter from CDB if we were requested
	// to initialise only when the first event was received.
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
						" receive raw DDL data from one tracking station DDL.",
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
		if (fLut == NULL)
		{
			HLTInfo("Loading lookup table information from CDB for DDL %d (ID = %d).",
				fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
			);
			int result = ReadLutFromCDB();
			if (result != 0) return result;
			
			fHitRec->SetLookUpTable(fLut, &fIdToEntry, &fMaxEntryPerBusPatch);
		}
		
		// Check that the DC cut was not already loaded in DoInit.
		if (fHitRec->GetDCCut() == -1)
		{
			HLTInfo("Loading DC cut parameters from CDB for DDL %d (ID = %d).",
				fDDL+1, AliHLTMUONUtils::DDLNumberToEquipId(fDDL)
			);
			int result = ReadDCCutFromCDB();
			if (result != 0) return result;
		}
		
		DoneDelayedSetup();
	}
	
	if (fLut == NULL)
	{
		HLTFatal("Lookup table not loaded! Cannot continue processing data.");
		return -ENOENT;
	}
	
	// Process an event
	unsigned long totalSize = 0; // Amount of memory currently consumed in bytes.

	HLTDebug("Processing event %llu with %u input data blocks.",
		evtData.fEventID, evtData.fBlockCnt
	);

	// Loop over all input blocks in the event
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
			n, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fPtr, blocks[n].fSize
		);

		if (blocks[n].fDataType != AliHLTMUONConstants::DDLRawDataType()
		    or not AliHLTMUONUtils::IsTrackerDDL(blocks[n].fSpecification)
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
		
		if (fDDL != -1)
		{
			AliHLTInt32_t receivedDDL = AliHLTMUONUtils::SpecToDDLNumber(blocks[n].fSpecification);
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
		AliHLTMUONRecHitsBlockWriter block(outputPtr+totalSize, size-totalSize);
		if (not block.InitCommonHeader())
		{
			HLTError("There is not enough space in the output buffer for the new data block."
				 " We require at least %u bytes, but have %u bytes left.",
				sizeof(AliHLTMUONRecHitsBlockWriter::HeaderType),
				block.BufferSize()
			);
			break;
		}
		
	    AliHLTCDHWrapper cdh(blocks[n].fPtr);

		AliHLTUInt32_t headerSize = cdh.GetHeaderSize();

		AliHLTUInt32_t totalDDLSize = blocks[n].fSize / sizeof(AliHLTUInt32_t);
		AliHLTUInt32_t ddlRawDataSize = totalDDLSize - headerSize / sizeof(AliHLTUInt32_t);
		AliHLTUInt32_t* buffer = reinterpret_cast<AliHLTUInt32_t*>(blocks[n].fPtr)
    + headerSize/sizeof(AliHLTUInt32_t);
		AliHLTUInt32_t nofHit = block.MaxNumberOfEntries();

#ifdef DEBUG
		HLTDebug("=========== Dumping DDL payload buffer ==========");
		for (AliHLTUInt32_t j = 0; j < totalDDLSize; j++)
			HLTDebug("buffer[%d] : %x",j,buffer[j]);
		HLTDebug("================== End of dump =================");
#endif // DEBUG

		if (not fHitRec->Run(buffer, ddlRawDataSize, block.GetArray(), nofHit))
		{
			HLTError("Error while processing the hit reconstruction algorithm.");
			if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
			size = totalSize; // Must tell the framework how much buffer space was used.
			return -EIO;
		}
		
		// nofHit should now contain the number of reconstructed hits actually found
		// and filled into the output data block, so we can set this number.
		assert( nofHit <= block.MaxNumberOfEntries() );
		block.SetNumberOfEntries(nofHit);
		
		HLTDebug("Number of reconstructed hits found is %d", nofHit);

		// Fill a block data structure for our output block.
		AliHLTComponentBlockData bd;
		FillBlockData(bd);
		bd.fPtr = outputPtr;
		// This block's start (offset) is after all other blocks written so far.
		bd.fOffset = totalSize;
		bd.fSize = block.BytesUsed();
		bd.fDataType = AliHLTMUONConstants::RecHitsBlockDataType();
		bd.fSpecification = blocks[n].fSpecification;
		outputBlocks.push_back(bd);

		// Increase the total amount of data written so far to our output memory
		totalSize += block.BytesUsed();
		
		if (fHitRec->GenerateClusterInfo())
		{
			// Create a new output clusters data block and initialise the header.
			AliHLTMUONClustersBlockWriter clustblock(outputPtr+totalSize, size-totalSize);
			if (not clustblock.InitCommonHeader())
			{
				HLTError("There is not enough space in the output buffer for the new clusters data block."
					" We require at least %u bytes, but have %u bytes left.",
					sizeof(AliHLTMUONClustersBlockWriter::HeaderType),
					clustblock.BufferSize()
				);
				break;
			}
			
			AliHLTUInt32_t nofClusters = clustblock.MaxNumberOfEntries();
			bool filledOk = fHitRec->FillClusterData(clustblock.GetArray(), nofClusters);
			// nofClusters should now contain the number of clusters filled.
			assert( nofClusters <= clustblock.MaxNumberOfEntries() );
			clustblock.SetNumberOfEntries(nofClusters);
			
			// Fill a block data structure for our output block.
			AliHLTComponentBlockData bdc;
			FillBlockData(bdc);
			bdc.fPtr = outputPtr;
			// This block's start (offset) is after all other blocks written so far.
			bdc.fOffset = totalSize;
			bdc.fSize = clustblock.BytesUsed();
			bdc.fDataType = AliHLTMUONConstants::ClusterBlockDataType();
			bdc.fSpecification = blocks[n].fSpecification;
			outputBlocks.push_back(bdc);
	
			// Increase the total amount of data written so far to our output memory
			totalSize += clustblock.BytesUsed();
			
			if (not filledOk)
			{
				HLTError("We have overflowed the output buffer space for the new clusters data block.");
				break;
			}
		}
		
		if (fHitRec->GenerateChannelInfo())
		{
			// Create a new output channels data block and initialise the header.
			AliHLTMUONChannelsBlockWriter channelblock(outputPtr+totalSize, size-totalSize);
			if (not channelblock.InitCommonHeader())
			{
				HLTError("There is not enough space in the output buffer for the new channels data block."
					" We require at least %u bytes, but have %u bytes left.",
					sizeof(AliHLTMUONChannelsBlockWriter::HeaderType),
					channelblock.BufferSize()
				);
				break;
			}
			
			AliHLTUInt32_t nofChannels = channelblock.MaxNumberOfEntries();
			bool filledOk = fHitRec->FillChannelData(channelblock.GetArray(), nofChannels);
			// nofChannels should now contain the number of channels filled.
			assert( nofChannels <= channelblock.MaxNumberOfEntries() );
			channelblock.SetNumberOfEntries(nofChannels);
			
			// Fill a block data structure for our output block.
			AliHLTComponentBlockData bdc;
			FillBlockData(bdc);
			bdc.fPtr = outputPtr;
			// This block's start (offset) is after all other blocks written so far.
			bdc.fOffset = totalSize;
			bdc.fSize = channelblock.BytesUsed();
			bdc.fDataType = AliHLTMUONConstants::ChannelBlockDataType();
			bdc.fSpecification = blocks[n].fSpecification;
			outputBlocks.push_back(bdc);
	
			// Increase the total amount of data written so far to our output memory
			totalSize += channelblock.BytesUsed();
			
			if (not filledOk)
			{
				HLTError("We have overflowed the output buffer space for the new channels data block.");
				break;
			}
		}
	}
	// Finally we set the total size of output memory we consumed.
	size = totalSize;

	return 0;
}


void AliHLTMUONHitReconstructorComponent::FreeMemory()
{
	/// Deletes any allocated objects if they are allocated else nothing is
	/// done for objects not yet allocated.
	/// This is used as a helper method to make sure the corresponding pointers
	/// are NULL and we get back to a well defined state.

	if (fHitRec != NULL)
	{
		delete fHitRec;
		fHitRec = NULL;
	}
	if (fLut != NULL)
	{
		delete [] fLut;
		fLut = NULL;
		fLutSize = 0;
	}
	
	fIdToEntry.clear();
	fMaxEntryPerBusPatch.clear();
}


int AliHLTMUONHitReconstructorComponent::ReadLookUpTable(const char* lutFileName)
{
	/// Read in the lookup table from a text file.
	/// Note that this method could leave fLut allocated which is cleaned up
	/// by DoInit with a call to FreeMemory().
	
	assert( fLut == NULL );
	assert( fLutSize == 0 );
	assert( fIdToEntry.empty() );
	assert( fMaxEntryPerBusPatch.empty() );
	
	std::ifstream file(lutFileName);
	if (not file.good())
	{
		HLTError("Could not open the LUT file %s", lutFileName);
		return -ENOENT;
	}
	
	// First count the number of lines of text in the LUT file before decoding.
	// This is not the most optimal. It would be better to read and decode at the
	// same time but we are not allowed to use STL and ROOT containers are too
	// heavy for this task. At least this is done only at the start of run.
	std::string str;
	AliHLTUInt32_t lineCount = 0;
	while (std::getline(file, str)) lineCount++;
	if (not file.eof())
	{
		HLTError("There was a problem reading the LUT file %s", lutFileName);
		return -EIO;
	}
	if (lineCount == 0)
	{
		HLTWarning("The LUT file %s was empty.", lutFileName);
	}
	
	// Add one extra LUT line for the first element which is used as a sentinel value.
	lineCount++;
	
	try
	{
		fLut = new AliHLTMUONHitRecoLutRow[lineCount];
		fLutSize = lineCount;
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the lookuptable.");
		return -ENOMEM;
	}
	
	// Initialise the sentinel value.
	fLut[0].fDetElemId = 0;
	fLut[0].fIX = 0;
	fLut[0].fIY = 0;
	fLut[0].fRealX = 0.0;
	fLut[0].fRealY = 0.0;
	fLut[0].fRealZ = 0.0;
	fLut[0].fHalfPadSize = 0.0;
	fLut[0].fPlane = -1;
	fLut[0].fPed = -1;
	fLut[0].fSigma = -1;
	fLut[0].fA0 = -1;
	fLut[0].fA1 = -1;
	fLut[0].fThres = -1;
	fLut[0].fSat = -1;
	
	// Clear the eof flag and start reading from the beginning of the file again.
	file.clear();
	file.seekg(0, std::ios::beg);
	if (not file.good())
	{
		HLTError("There was a problem seeking in the LUT file %s", lutFileName);
		return -EIO;
	}
	
	AliHLTInt32_t idManuChannel,buspatchId;
	for (AliHLTUInt32_t i = 1; i < fLutSize; i++)
	{
		if (std::getline(file, str).fail())
		{
			HLTError("There was a problem reading line %d of LUT file %s", i, lutFileName);
			return -EIO;
		}
		
		int result = sscanf(
			str.c_str(), "%14d\t%14d\t%14d\t%14d\t%23e\t%23e\t%23e\t%23e\t%14d\t%23e\t%23e\t%23e\t%23e\t%14d\t%14d",
			&idManuChannel, &fLut[i].fDetElemId, &fLut[i].fIX,
			&fLut[i].fIY, &fLut[i].fRealX,
			&fLut[i].fRealY, &fLut[i].fRealZ,
			&fLut[i].fHalfPadSize, &fLut[i].fPlane,
			&fLut[i].fPed, &fLut[i].fSigma, &fLut[i].fA0,
			&fLut[i].fA1, &fLut[i].fThres, &fLut[i].fSat
		);
		
		if (result != 15)
		{
			HLTError("Line %d in LUT file %s does not contain 15 elements.", i, lutFileName);
			return -EIO;
		}
		buspatchId = (idManuChannel>>17) & 0x7FF;
		fIdToEntry[idManuChannel] = i;
		fMaxEntryPerBusPatch[buspatchId] = fMaxEntryPerBusPatch[buspatchId] + 1;  
		
	}

	MaxEntryPerBusPatch::iterator it;
	for(it=fMaxEntryPerBusPatch.begin(); it!=fMaxEntryPerBusPatch.end(); it++){
		HLTDebug("fMaxEntryPerBusPatch[%d] : %d",it->first,it->second);
		fMaxEntryPerBusPatch[it->first] = AliHLTInt32_t(0.05*(it->second));///< for 10% occupancy 
		HLTDebug("fMaxEntryPerBusPatch[%d] : %d",it->first,it->second);
	}
	
	return 0;
}


int AliHLTMUONHitReconstructorComponent::ReadLutFromCDB()
{
	/// Reads LUT from CDB.
	/// To override the default CDB path and / or run number the
	/// SetCDBPathAndRunNo(cdbPath, run) method should be called before this
	/// method.

	assert( fLut == NULL );
	assert( fLutSize == 0 );
	assert( fIdToEntry.empty() );
	assert( fMaxEntryPerBusPatch.empty() );
	
	if (fDDL == -1)
	{
		HLTError("No DDL number specified for which to load LUT data from CDB.");
		return -EINVAL;
	}
	
	std::vector<AliHLTMUONHitRecoLutRow> lutList;
	AliHLTMUONHitRecoLutRow lut;
	AliHLTUInt32_t iEntry = 0;
	
	int result = FetchMappingStores();
	// Error message already generated in FetchMappingStores.
	if (result != 0) return result;
	AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
	
	AliMpSegmentation* mpSegFactory = AliMpSegmentation::Instance();
	if (mpSegFactory == NULL)
	{
		HLTError("Could not find segmentation mapping (AliMpSegmentation) instance.");
		return -EIO;
	}
	
	// Only load geometry if not already loaded.
	if (AliGeomManager::GetGeometry() == NULL)
	{
		AliGeomManager::LoadGeometry();
	}
	AliMUONGeometryTransformer chamberGeometryTransformer;
	if (not chamberGeometryTransformer.LoadGeometryData())
	{
		HLTError("Failed to load geomerty data.");
		return -ENOENT;
	}

	AliMUONCalibrationData calibData(AliCDBManager::Instance()->GetRun());
	
	bool skippedPads = false;
	Int_t chamberId;
	
	for(Int_t iCh = 0; iCh < 10; iCh++)
	{
		chamberId = iCh;
		
		AliMpDEIterator it;
		for ( it.First(chamberId); ! it.IsDone(); it.Next() )
		{
			Int_t detElemId = it.CurrentDEId();
			int iDDL = ddlStore->GetDetElement(detElemId)->GetDdlId();
			if (iDDL != fDDL) continue;
		
			for (Int_t iCath = 0 ; iCath <= 1 ; iCath++)
			{
				AliMp::CathodType cath;
				
				if(iCath == 0)
					cath = AliMp::kCath0;
				else
					cath = AliMp::kCath1;
				
				const AliMpVSegmentation* seg = mpSegFactory->GetMpSegmentation(detElemId, cath);
				AliMp::PlaneType plane = seg->PlaneType();
				Int_t maxIX = seg->MaxPadIndexX();
				Int_t maxIY = seg->MaxPadIndexY();
				
				Int_t idManuChannel, manuId, channelId, buspatchId;
				AliHLTFloat32_t padSizeX, padSizeY;
				AliHLTFloat32_t halfPadSize;
				AliHLTFloat32_t padSizeXY;
				Double_t realX, realY, realZ;
				Double_t localX, localY, localZ;
				Float_t calibA0Coeff,calibA1Coeff,pedestal,sigma;
				Float_t thresold,saturation;
				
				// Pad Info of a slat to print in lookuptable
				for (Int_t iX = 0; iX<= maxIX ; iX++)
				for (Int_t iY = 0; iY<= maxIY ; iY++)
				{
#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					if (not seg->HasPadByIndices(iX,iY)) continue;
#else // old AliMpPad functionality < r 31742
					if (not seg->HasPad(AliMpIntPair(iX,iY))) continue;
#endif

#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					AliMpPad pad = seg->PadByIndices(iX,iY, kFALSE);
#else // old AliMpPad functionality < r 31742
					AliMpPad pad = seg->PadByIndices(AliMpIntPair(iX,iY), kFALSE);
#endif
					
					// Getting Manu id
#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					manuId = pad.GetManuId();
#else // old AliMpPad functionality < r 31742
					manuId = pad.GetLocation().GetFirst();
#endif
					manuId &= 0x7FF; // 11 bits 
					if (calibData.Pedestals(detElemId, manuId) == NULL) continue;
			
					buspatchId = ddlStore->GetBusPatchId(detElemId,manuId);
					
					// Getting channel id
#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					channelId =  pad.GetManuChannel();
#else // old AliMpPad functionality < r 31742
					channelId =  pad.GetLocation().GetSecond();
#endif
					channelId &= 0x3F; // 6 bits
					
					idManuChannel = buspatchId << 11;
					idManuChannel = (idManuChannel | manuId) << 6;
					idManuChannel |= channelId;
					
#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					localX = pad.GetPositionX();
					localY = pad.GetPositionY();
#else // old AliMpPad functionality < r 31769
					localX = pad.Position().X();
					localY = pad.Position().Y();
#endif //HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					localZ = 0.0;
					
					chamberGeometryTransformer.Local2Global(
						detElemId,localX,localY,localZ,
						realX,realY,realZ
					);
					
#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					padSizeX = AliHLTFloat32_t( pad.GetDimensionX() );
					padSizeY = AliHLTFloat32_t( pad.GetDimensionY() );
#else // old AliMpPad functionality < r 31769
					padSizeX = AliHLTFloat32_t( pad.Dimensions().X() );
					padSizeY = AliHLTFloat32_t( pad.Dimensions().Y() );
#endif //HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					
          calibA0Coeff = 1.0;
          calibA1Coeff = 0.0;
          thresold = 0;
          saturation = 0;
					
					pedestal = (calibData.Pedestals(detElemId, manuId))->ValueAsFloat(channelId, 0);
					sigma = (calibData.Pedestals(detElemId, manuId))->ValueAsFloat(channelId, 1);
					
					// Check if any of the values fetched from the calibration data are
					// invalid. If they are then skip this pad.
					if (calibA0Coeff == AliMUONVCalibParam::InvalidFloatValue() or
					    calibA1Coeff == AliMUONVCalibParam::InvalidFloatValue() or
					    thresold == AliMUONVCalibParam::InvalidFloatValue() or
					    saturation == AliMUONVCalibParam::InvalidFloatValue() or
					    pedestal == AliMUONVCalibParam::InvalidFloatValue() or
					    sigma == AliMUONVCalibParam::InvalidFloatValue()
					   )
					{
						if (fWarnIfPadSkipped)
						{
							HLTWarning("Skipping pad on detection element %d, MANU %d, channel %d, since"
								" the calibration data contains invalid values in that channel.",
								detElemId, manuId, channelId
							);
						}
						skippedPads = true;
						continue;
					}
					

					if (plane == 0){
						halfPadSize = padSizeX;
						padSizeXY = padSizeY;
					}else{
						halfPadSize = padSizeY;
						padSizeXY = padSizeX;
					}

					fIdToEntry[idManuChannel] = iEntry+1;
					fMaxEntryPerBusPatch[buspatchId] = fMaxEntryPerBusPatch[buspatchId] + 1;  
					
					lut.fDetElemId = detElemId;
					lut.fIX = iX;
					lut.fIY = iY;
					lut.fRealX = realX;
					lut.fRealY = realY;
					lut.fRealZ = realZ;
					lut.fHalfPadSize = halfPadSize;
					lut.fPadSizeXY = padSizeXY;
					lut.fPlane = plane;
					lut.fPed = pedestal;
					lut.fSigma = sigma;
					lut.fA0 = calibA0Coeff;
					lut.fA1 = calibA1Coeff;
					lut.fThres = Int_t(thresold);
					lut.fSat = Int_t(saturation);
					
					HLTDebug("lut : detele : %d, id : %d, manu : %d, channel : %d, iX : %d, iY: %d, (X,Y) : (%f, %f), padsize : %f, plane : %d, ped : %f, sigma : %f",
						lut.fDetElemId,idManuChannel,manuId,channelId,lut.fIX,lut.fIY,lut.fRealX,lut.fRealY,lut.fHalfPadSize,lut.fPlane,lut.fPed,lut.fSigma
					);
					
					lutList.push_back(lut);
					iEntry++;
				} // iX, iY loop
			} // iCath loop
		} // detElemId loop
	} // ichamber loop
	
	if (skippedPads and not fWarnIfPadSkipped)
	{
		HLTWarning("Skipped pads since they contained invalid calibration values."
			" Use the -warnifpadskipped argument to generate detailed information"
			" about which pads were skipped."
		);
	}

	try
	{
		// Use iEntry+1 since we add one extra LUT line for the first element
		// which is used as a sentinel value.
		fLut = new AliHLTMUONHitRecoLutRow[iEntry+1];
		fLutSize = iEntry+1;
		HLTDebug("Address of new LUT buffer at fLut = %p", fLut);
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the lookuptable.");
		//lutList.clear();  not necessary, implicitly done during stack cleanup.
		return -ENOMEM;
	}
	
	// Initialise the sentinel value.
	fLut[0].fDetElemId = 0;
	fLut[0].fIX = 0;
	fLut[0].fIY = 0;
	fLut[0].fRealX = 0.0;
	fLut[0].fRealY = 0.0;
	fLut[0].fRealZ = 0.0;
	fLut[0].fHalfPadSize = 0.0;
	fLut[0].fPadSizeXY = 0.0;
	fLut[0].fPlane = -1;
	fLut[0].fPed = -1;
	fLut[0].fSigma = -1;
	fLut[0].fA0 = -1;
	fLut[0].fA1 = -1;
	fLut[0].fThres = -1;
	fLut[0].fSat = -1;
	
	for (AliHLTUInt32_t i = 0; i < iEntry; i++)
		fLut[i+1] = lutList[i];
	lutList.clear();

	MaxEntryPerBusPatch::iterator it;
	for(it=fMaxEntryPerBusPatch.begin(); it!=fMaxEntryPerBusPatch.end(); it++){
		HLTDebug("fMaxEntryPerBusPatch[%d] : %d",it->first,it->second);
		fMaxEntryPerBusPatch[it->first] = AliHLTInt32_t(0.05*(it->second));///< for 10% occupancy
		HLTDebug("fMaxEntryPerBusPatch[%d] : %d",it->first,it->second);
	}

	return 0;
}


int AliHLTMUONHitReconstructorComponent::ReadDCCutFromCDB()
{
	/// Reads the DC cut parameter from the CDB.

	const char* pathToEntry = AliHLTMUONConstants::HitReconstructorCDBPath();
	
	TMap* map = NULL;
	int result = FetchTMapFromCDB(pathToEntry, map);
	if (result != 0) return result;
	
	Int_t value = 0;
	result = GetPositiveIntFromTMap(map, "dccut", value, pathToEntry, "DC cut");
	if (result != 0) return result;
	
	assert(fHitRec != NULL);
	fHitRec->SetDCCut(value);
	
	HLTDebug("Using DC cut parameter of %d ADC channels.", fHitRec->GetDCCut());
	
	return 0;
}


bool AliHLTMUONHitReconstructorComponent::GenerateLookupTable(
		AliHLTInt32_t ddl, const char* filename,
		const char* cdbPath, Int_t run
	)
{
	/// Generates a ASCII text file containing the lookup table (LUT) from
	/// the CDB, which can be used for the hit reconstructor component later.
	/// @param ddl  Must be the DDL for which to generate the DDL,
	///             in the range [0..19].
	/// @param filename  The name of the LUT file to generate.
	/// @param cdbPath  The CDB path to use.
	/// @param run  The run number to use for the CDB.
	/// @return  True if the generation of the LUT file succeeded.
	
	AliHLTMUONHitReconstructorComponent comp;
	
	if (ddl < 0 or 19 < ddl)
	{
		std::cerr << "ERROR: the DDL number must be in the range [0..19]." << std::endl;
		return false;
	}
	
	comp.fDDL = ddl;
	if (comp.SetCDBPathAndRunNo(cdbPath, run) != 0) return false;
	if (comp.ReadLutFromCDB() != 0) return false;
	
	char str[1024*4];
	std::fstream file(filename, std::ios::out);
	if (not file)
	{
		std::cerr << "ERROR: could not open file: " << filename << std::endl;
		return false;
	}
	
	assert( comp.fLut != NULL );
	
	for (IdManuChannelToEntry::iterator id = comp.fIdToEntry.begin();
	     id != comp.fIdToEntry.end();
	     id++
	    )
	{
		AliHLTInt32_t idManuChannel = id->first;
		AliHLTInt32_t row = id->second;
		
		assert( AliHLTUInt32_t(row) < comp.fLutSize );
		
		sprintf(str, "%d\t%d\t%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\t%d\t%.15e\t%.15e\t%.15e\t%.15e\t%d\t%d",
			idManuChannel, comp.fLut[row].fDetElemId, comp.fLut[row].fIX,
			comp.fLut[row].fIY, comp.fLut[row].fRealX,
			comp.fLut[row].fRealY, comp.fLut[row].fRealZ,
			comp.fLut[row].fHalfPadSize, comp.fLut[row].fPlane,
			comp.fLut[row].fPed, comp.fLut[row].fSigma, comp.fLut[row].fA0,
			comp.fLut[row].fA1, comp.fLut[row].fThres, comp.fLut[row].fSat
		);
		
		file << str << endl;
		if (file.fail())
		{
			std::cerr << "ERROR: There was an I/O error when writing to the file: "
				<< filename << std::endl;
			return false;
		}
	}
	
	return true;
}
