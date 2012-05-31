/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
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
/// @file   AliHLTMUONRecHitsSource.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 May 2007
/// @brief  Implementation of the AliHLTMUONRecHitsSource component.
///

#include "AliHLTMUONRecHitsSource.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONDataInterface.h"
#include "AliMUONHit.h"
#include "AliMUONVCluster.h"
#include "AliMUONConstants.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVHitStore.h"
#include "mapping/AliMpCDB.h"
#include "mapping/AliMpDEManager.h"
#include "mapping/AliMpDetElement.h"
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <cassert>
#include <new>

ClassImp(AliHLTMUONRecHitsSource);


AliHLTMUONRecHitsSource::AliHLTMUONRecHitsSource() :
	AliHLTOfflineDataSource(),
	fMCDataInterface(NULL),
	fDataInterface(NULL),
	fSelection(kWholePlane),
	fServeChamber(),
	fCurrentEventIndex(0)
{
	///
	/// Default constructor.
	///

	for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++)
		fServeChamber[i] = false;
}


AliHLTMUONRecHitsSource::~AliHLTMUONRecHitsSource()
{
	///
	/// Default destructor.
	///
	
	if (fMCDataInterface != NULL) delete fMCDataInterface;
	if (fDataInterface != NULL) delete fDataInterface;
}


int AliHLTMUONRecHitsSource::DoInit(int argc, const char** argv)
{
	///
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	///
	
	HLTInfo("Initialising dHLT reconstructed hit source component.");

	if (fMCDataInterface != NULL)
	{
		delete fMCDataInterface;
		fMCDataInterface = NULL;
	}
	if (fDataInterface != NULL)
	{
		delete fDataInterface;
		fDataInterface = NULL;
	}
	
	// Parse the command line arguments:
	bool simdata = false;
	bool recdata = false;
	bool chamberWasSet = false;
	fCurrentEventIndex = 0;
	bool firstEventSet = false;
	bool eventNumLitSet = false;
	
	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-simdata") == 0)
		{
			simdata = true;
		}
		else if (strcmp(argv[i], "-recdata") == 0)
		{
			recdata = true;
		}
		else if (strcmp(argv[i], "-plane") == 0)
		{
			i++;
			if (i >= argc)
			{
				Logging(kHLTLogError,
					"AliHLTMUONRecHitsSource::DoInit",
					"Missing parameter",
					"Expected one of 'left', 'right' or 'all' after '-plane'."
				);
				return -EINVAL;
			}
			if (strcmp(argv[i], "left") == 0)
				fSelection = kLeftPlane;
			else if (strcmp(argv[i], "right") == 0)
				fSelection = kRightPlane;
			else if (strcmp(argv[i], "all") == 0)
				fSelection = kWholePlane;
			else
			{
				Logging(kHLTLogError,
					"AliHLTMUONRecHitsSource::DoInit",
					"Invalid parameter",
					"The parameter '%s' is invalid and must be one of 'left',"
					  " 'right' or 'all'.",
					argv[i]
				);
				return -EINVAL;
			}
		}
		else if (strcmp(argv[i], "-chamber") == 0)
		{
			i++;
			if (i >= argc)
			{
				Logging(kHLTLogError,
					"AliHLTMUONRecHitsSource::DoInit",
					"Missing parameter",
					"Expected a chamber number, range eg. '1-10' or list eg."
					  " '1,2,3' after '-chamber'."
				);
				return -EINVAL;
			}
			int result = ParseChamberString(argv[i]);
			if (result != 0) return result;
			chamberWasSet = true;
		}
		else if (strcmp(argv[i], "-firstevent") == 0)
		{
			if (eventNumLitSet)
			{
				HLTWarning("The -firstevent flag is overridden by a"
					" previous use of -event_number_literal."
				);
			}
			i++;
			if (i >= argc)
			{
				HLTError("Expected a positive number after -firstevent.");
				return -EINVAL;
			}
			char* end = NULL;
			long num = strtol(argv[i], &end, 0);
			if ((end != NULL and *end != '\0') or num < 0) // Check if the conversion is OK.
			{
				HLTError(
					"Expected a positive number after -firstevent"
					" but got: %s", argv[i]
				);
				return -EINVAL;
			}
			fCurrentEventIndex = Int_t(num);
			firstEventSet = true;
		}
		else if (strcmp(argv[i], "-event_number_literal") == 0)
		{
			if (firstEventSet)
			{
				HLTWarning("The -event_number_literal option will"
					" override -firstevent."
				);
			}
			fCurrentEventIndex = -1;
			eventNumLitSet = true;
		}
		else
		{
			Logging(kHLTLogError,
				"AliHLTMUONRecHitsSource::DoInit",
				"Unknown argument",
				"The argument '%s' is invalid.",
				argv[i]
			);
			return -EINVAL;
		}
	}

	// Check the parameters we have parsed.
	if (simdata and recdata)
	{
		Logging(kHLTLogError,
			"AliHLTMUONRecHitsSource::DoInit",
			"Invalid arguments",
			"Cannot have both -simdata and -recdata set."
		);
		return -EINVAL;
	}
	
	if (not simdata and not recdata)
	{
		Logging(kHLTLogError,
			"AliHLTMUONRecHitsSource::DoInit",
			"Missing arguments",
			"Must have either -simdata or -recdata specified."
		);
		return -EINVAL;
	}
	
	if (not chamberWasSet)
	{
		Logging(kHLTLogInfo,
			"AliHLTMUONRecHitsSource::DoInit",
			"Setting Parameters",
			"No chambers were selected so we will publish for all chambers."
		);
		for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++)
			fServeChamber[i] = true;
	}
	
	// Must load the mapping data for AliMpDetElement::GetDdlId()
	// to return useful information later on.
	AliMpCDB::LoadDDLStore();
		
	// Now we can initialise the data interface objects and loaders.
	if (simdata)
	{
		Logging(kHLTLogDebug,
			"AliHLTMUONRecHitsSource::DoInit",
			"Data interface",
			"Loading simulated GEANT hits with AliMUONMCDataInterface."
		);

		try
		{
			fMCDataInterface = new AliMUONMCDataInterface("galice.root");
		}
		catch (const std::bad_alloc&)
		{
			Logging(kHLTLogError,
				"AliHLTMUONRecHitsSource::DoInit",
				"Out of memory",
				"Not enough memory to allocate AliMUONMCDataInterface."
			);
			return -ENOMEM;
		}
	}
	else if (recdata)
	{
		Logging(kHLTLogDebug,
			"AliHLTMUONRecHitsSource::DoInit",
			"Data interface",
			"Loading reconstructed clusters with AliMUONDataInterface."
		);
		
		try
		{
			fDataInterface = new AliMUONDataInterface("galice.root");
		}
		catch (const std::bad_alloc&)
		{
			Logging(kHLTLogError,
				"AliHLTMUONRecHitsSource::DoInit",
				"Out of memory",
				"Not enough memory to allocate AliMUONDataInterface."
			);
			return -ENOMEM;
		}
	}
	
	// Check that the fCurrentEventIndex number falls within the correct range.
	UInt_t maxevent = 0;
	if (fMCDataInterface != NULL)
		maxevent = UInt_t(fMCDataInterface->NumberOfEvents());
	else if (fDataInterface != NULL)
		maxevent = UInt_t(fDataInterface->NumberOfEvents());
	if (fCurrentEventIndex != -1 and UInt_t(fCurrentEventIndex) >= maxevent and maxevent != 0)
	{
		fCurrentEventIndex = 0;
		HLTWarning(Form("The selected first event number (%d) was larger than"
			" the available number of events (%d). Resetting the event"
			" counter to zero.", fCurrentEventIndex, maxevent
		));
	}
	
	return 0;
}


int AliHLTMUONRecHitsSource::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	
	HLTInfo("Deinitialising dHLT reconstructed hit source component.");
	
	if (fMCDataInterface != NULL)
	{
		delete fMCDataInterface;
		fMCDataInterface = NULL;
	}
	if (fDataInterface != NULL)
	{
		delete fDataInterface;
		fDataInterface = NULL;
	}
	return 0;
}


const char* AliHLTMUONRecHitsSource::GetComponentID()
{
	///
	/// Inherited from AliHLTComponent. Returns the component ID.
	///
	
	return AliHLTMUONConstants::RecHitsSourceId();
}


AliHLTComponentDataType AliHLTMUONRecHitsSource::GetOutputDataType()
{
	///
	/// Inherited from AliHLTComponent. Returns the output data type.
	///
	
	return AliHLTMUONConstants::RecHitsBlockDataType();
}


void AliHLTMUONRecHitsSource::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	///
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	///
	
	constBase = sizeof(AliHLTMUONRecHitsBlockStruct)
		+ 256*16*sizeof(AliHLTMUONRecHitStruct);
	inputMultiplier = 0;
}


AliHLTComponent* AliHLTMUONRecHitsSource::Spawn()
{
	///
	/// Inherited from AliHLTComponent. Creates a new object instance.
	///
	
	return new AliHLTMUONRecHitsSource();
}


int AliHLTMUONRecHitsSource::GetEvent(
		const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks
	)
{
	///
	/// Inherited from AliHLTOfflineDataSource. Creates new event data blocks.
	///
	
	assert( fMCDataInterface != NULL or fDataInterface != NULL );
	
	if (not IsDataEvent()) return 0;  // ignore non data events.

	// Check the size of the event descriptor structure.
	if (evtData.fStructSize < sizeof(AliHLTComponentEventData))
	{
		Logging(kHLTLogError,
			"AliHLTMUONRecHitsSource::GetEvent",
			"Invalid event descriptor",
			"The event descriptor (AliHLTComponentEventData) size is"
			  " smaller than expected. It claims to be %d bytes, but"
			  " we expect it to be %d bytes.",
			evtData.fStructSize,
			sizeof(AliHLTComponentEventData)
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -EINVAL;
	}
	
	// Use the fEventID as the event number to load if fCurrentEventIndex == -1,
	// check it and load that event with the runloader.
	// If fCurrentEventIndex is a positive number then use it instead and
	// increment it.
	UInt_t eventnumber = UInt_t(evtData.fEventID);
	UInt_t maxevent = 0;
	if (fMCDataInterface != NULL)
		maxevent = UInt_t(fMCDataInterface->NumberOfEvents());
	else if (fDataInterface != NULL)
		maxevent = UInt_t(fDataInterface->NumberOfEvents());
	if (fCurrentEventIndex != -1)
	{
		eventnumber = UInt_t(fCurrentEventIndex);
		fCurrentEventIndex++;
		if (UInt_t(fCurrentEventIndex) >= maxevent)
			fCurrentEventIndex = 0;
	}
	if ( eventnumber >= maxevent )
	{
		Logging(kHLTLogError,
			"AliHLTMUONRecHitsSource::GetEvent",
			"Bad event ID",
			"The event number (%d) is larger than the available number"
			  " of events on file (%d).",
			eventnumber,
			maxevent
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -EINVAL;
	}
	
	// Create and initialise a new data block.
	AliHLTMUONRecHitsBlockWriter block(outputPtr, size);
	if (not block.InitCommonHeader())
	{
		Logging(kHLTLogError,
			"AliHLTMUONRecHitsSource::GetEvent",
			"Buffer too small",
			"There is not enough buffer space to create a new data block."
			  " We require at least %d bytes but the buffer is only %d bytes.",
			sizeof(AliHLTMUONRecHitsBlockWriter::HeaderType),
			block.BufferSize()
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -ENOBUFS;
	}
	
	// Initialise the DDL list containing the DDLs which contributed to the
	// data block. These are required to create the specification word later.
	bool ddlList[22];
	for (Int_t i = 0; i < 22; i++)
		ddlList[i] = false;
	
	if (fMCDataInterface != NULL)
	{
		Logging(kHLTLogDebug,
			"AliHLTMUONRecHitsSource::GetEvent",
			"Filling hits",
			"Filling data block with GEANT hits for event %d.",
			eventnumber
		);
		
		// Loop over all tracks, extract the hits and write them to the
		// data block.
		Int_t ntracks = fMCDataInterface->NumberOfTracks(eventnumber);
		for (Int_t i = 0; i < ntracks; ++i)
		{
			AliMUONVHitStore* hitStore = fMCDataInterface->HitStore(eventnumber,i);
			AliMUONHit* hit;
			TIter next(hitStore->CreateIterator());
			while ( ( hit = static_cast<AliMUONHit*>(next()) ) )
			{
				// Select only hits on selected chambers.
				Int_t chamber = hit->Chamber() - 1;
				if (chamber > AliMUONConstants::NTrackingCh()) continue;
				if (not fServeChamber[chamber]) continue;
				
				// Only select hits from the given part of the plane
				if (fSelection == kLeftPlane and not (hit->Xref() < 0)) continue;
				if (fSelection == kRightPlane and not (hit->Xref() >= 0)) continue;
				
				AliHLTMUONRecHitStruct* rechit = block.AddEntry();
				if (rechit == NULL)
				{
					Logging(kHLTLogError,
						"AliHLTMUONRecHitsSource::GetEvent",
						"Buffer overflow",
						"There is not enough buffer space to add more hits."
						  " We overflowed the buffer which is only %d bytes.",
						block.BufferSize()
					);
					size = 0; // Important to tell framework that nothing was generated.
					return -ENOBUFS;
				}
				
				rechit->fX = hit->Xref();
				rechit->fY = hit->Yref();
				rechit->fZ = hit->Zref();
				
				// Workout which DDL this hit will be readout of.
				AliMpDetElement* de = AliMpDEManager::GetDetElement(hit->DetElemId());
				if (de != NULL and (0 <= de->GetDdlId() and de->GetDdlId() < 22))
					ddlList[de->GetDdlId()] = true;
				else
					Logging(kHLTLogError,
						"AliHLTMUONRecHitsSource::GetEvent",
						"No DDL ID",
						"Could not find the DDL ID from which readout would take place."
					);
			}
		}
	}
	else if (fDataInterface != NULL)
	{
		Logging(kHLTLogDebug,
			"AliHLTMUONRecHitsSource::GetEvent",
			"Filling hits",
			"Filling data block with reconstructed raw clusters for event %d.",
			eventnumber
		);
		
		AliMUONVClusterStore* clusterStore = fDataInterface->ClusterStore(eventnumber);
    
		// Loop over selected chambers and extract the raw clusters.
		for (Int_t chamber = 0; chamber < AliMUONConstants::NTrackingCh(); chamber++)
		{
			// Select only hits on selected chambers.
			if (not fServeChamber[chamber]) continue;
			
			TIter next(clusterStore->CreateChamberIterator(chamber,chamber));
			AliMUONVCluster* cluster;
			while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) )
			{
				// Only select hits from the given part of the plane
				if (fSelection == kLeftPlane and not (cluster->GetX() < 0)) continue;
				if (fSelection == kRightPlane and not (cluster->GetX() >= 0)) continue;
			
				AliHLTMUONRecHitStruct* rechit = block.AddEntry();
				if (rechit == NULL)
				{
					Logging(kHLTLogError,
						"AliHLTMUONRecHitsSource::GetEvent",
						"Buffer overflow",
						"There is not enough buffer space to add more hits."
						  " We overflowed the buffer which is only %d bytes.",
						block.BufferSize()
					);
					size = 0; // Important to tell framework that nothing was generated.
					return -ENOBUFS;
				}
				
				rechit->fX = cluster->GetX();
				rechit->fY = cluster->GetY();
				rechit->fZ = cluster->GetZ();
				
				// Workout which DDL this hit will be readout of.
				AliMpDetElement* de = AliMpDEManager::GetDetElement(cluster->GetDetElemId());
				if (de != NULL and (0 <= de->GetDdlId() and de->GetDdlId() < 22))
					ddlList[de->GetDdlId()] = true;
				else
					Logging(kHLTLogError,
						"AliHLTMUONRecHitsSource::GetEvent",
						"No DDL ID",
						"Could not find the DDL ID from which readout would take place."
					);
			}
		}
	}
	else
	{
		Logging(kHLTLogError,
			"AliHLTMUONRecHitsSource::GetEvent",
			"Missing data interface",
			"Neither AliMUONDataInterface nor AliMUONMCDataInterface were created."
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -EFAULT;
	}
	
	AliHLTComponentBlockData bd;
	FillBlockData(bd);
	bd.fPtr = outputPtr;
	bd.fOffset = 0;
	bd.fSize = block.BytesUsed();
	bd.fDataType = AliHLTMUONConstants::RecHitsBlockDataType();
	bd.fSpecification = AliHLTMUONUtils::PackSpecBits(ddlList);
	outputBlocks.push_back(bd);
	size = block.BytesUsed();

	return 0;
}


int AliHLTMUONRecHitsSource::ParseChamberString(const char* str)
{
	///
	/// Parses a string with the following format:
	///   <number>|<number>-<number>[,<number>|<number>-<number>]...
	/// For example: 1  1,2,3  1-2   1,2-4,5  etc...
	/// Flags in the fServeChamber will be set to 'true' for all appropriate
	/// values parsed.
	/// @param str  The string to parse.
	/// @return  Zero on success and EINVAL if there is a parse error.
	///
	
	char* end = const_cast<char*>(str);
	long lastChamber = -1;
	do
	{
		// Parse the next number.
		char* current = end;
		long chamber = strtol(current, &end, 0);
		
		// Check for parse errors of the number.
		if (current == end)
		{
			Logging(kHLTLogError,
				"AliHLTMUONRecHitsSource::GetEvent",
				"Parse error",
				"Expected a number in the range [1..%d] but got '%s'.",
				AliMUONConstants::NTrackingCh(), current
			);
			return -EINVAL;
		}
		if (chamber < 1 or AliMUONConstants::NTrackingCh() < chamber)
		{
			Logging(kHLTLogError,
				"AliHLTMUONRecHitsSource::GetEvent",
				"Parse error",
				"Got the chamber number %d which is outside the valid range of [1..%d].",
				chamber, AliMUONConstants::NTrackingCh()
			);
			return -EINVAL;
		}
		
		// Skip any whitespace after the number
		while (*end != '\0' and (*end == ' ' or *end == '\t' or *end == '\r' or *end == '\n')) end++;
		
		// Check if we are dealing with a list or range, or if we are at
		// the end of the string.
		if (*end == '-')
		{
			lastChamber = chamber;
			end++;
			continue;
		}
		else if (*end == ',')
		{
			assert( 1 <= chamber and chamber <= 10 );
			fServeChamber[chamber-1] = true;
			end++;
		}
		else if (*end == '\0')
		{
			assert( 1 <= chamber and chamber <= 10 );
			fServeChamber[chamber-1] = true;
		}
		else
		{
			Logging(kHLTLogError,
				"AliHLTMUONRecHitsSource::GetEvent",
				"Parse error",
				"Could not understand parameter list '%s'. Expected '-', ','"
				  " or end of line but got '%c' at character %d.",
				str, *end, (int)(end - str) +1
			);
			return -EINVAL;
		}
		
		// Set the range of chambers to publish for.
		if (lastChamber > 0)
		{
			Int_t min, max;
			if (lastChamber < chamber)
			{
				min = lastChamber;
				max = chamber;
			}
			else
			{
				min = chamber;
				max = lastChamber;
			}
			assert( min >= 1 );
			assert( max <= 10 );
			for (Int_t i = min; i <= max; i++)
				fServeChamber[i-1] = true;
		}
		lastChamber = -1;
	}
	while (*end != '\0');
	return 0;
}
