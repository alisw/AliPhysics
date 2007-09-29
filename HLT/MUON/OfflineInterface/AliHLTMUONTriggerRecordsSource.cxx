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

/* $Id$ */

/**
 * @file   AliHLTMUONTriggerRecordsSource.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of the AliHLTMUONTriggerRecordsSource component.
 */

#include "AliHLTMUONTriggerRecordsSource.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliHLTMUONCalculations.h"
#include "AliMUONMCDataInterface.h"
#include "AliMUONDataInterface.h"
#include "AliMUONHit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONConstants.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVHitStore.h"
#include "mapping/AliMpCDB.h"
#include "mapping/AliMpDDLStore.h"
#include "mapping/AliMpLocalBoard.h"
#include "mapping/AliMpTriggerCrate.h"
#include "mapping/AliMpDEManager.h"
#include "mapping/AliMpDetElement.h"
#include "AliLog.h"
#include "TClonesArray.h"
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <cassert>
#include <new>

namespace
{
	// The global object used for automatic component registration.
	// Note DO NOT use this component for calculation!
	AliHLTMUONTriggerRecordsSource gAliHLTMUONTriggerRecordsSource;
	
	//TODO: The following method should be in MUON/mapping
	Int_t FindDDLOfDetElement(Int_t detElemId)
	{
		// Find what the DDL ID number is for a detector element from
		// trigger chambers 11 to 14. We first have to find the local
		// board associated with the detector element and then we can
		// associate that local board to the trigger crate which has
		// the DDL number specified.
		AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
		if (ddlStore == NULL) return -1;
		Int_t ddl = -1, boardIndex = 1;
		do
		{
			AliMpLocalBoard* board = ddlStore->GetLocalBoard(boardIndex++);
			if (board == NULL) break;
			if (board->HasDEId(detElemId))
			{
				AliMpTriggerCrate* crate = ddlStore->GetTriggerCrate(board->GetCrate());
				if (crate == NULL) continue;
				ddl = crate->GetDdlId();
				break;
			}
		}
		while (ddl == -1);
		return ddl;
	}

}


ClassImp(AliHLTMUONTriggerRecordsSource);


AliHLTMUONTriggerRecordsSource::AliHLTMUONTriggerRecordsSource() :
	AliHLTOfflineDataSource(),
	fMCDataInterface(NULL),
	fDataInterface(NULL),
	fBuildFromHits(false),
	fSelection(kWholePlane),
	fCurrentEvent(0)
{
}


AliHLTMUONTriggerRecordsSource::~AliHLTMUONTriggerRecordsSource()
{
	assert( fMCDataInterface == NULL );
	assert( fDataInterface == NULL );
}


int AliHLTMUONTriggerRecordsSource::DoInit(int argc, const char** argv)
{
	assert( fMCDataInterface == NULL );
	assert( fDataInterface == NULL );
	
	// Parse the command line arguments:
	bool hitdata = false;
	bool simdata = false;
	bool recdata = false;
	fCurrentEvent = 0;
	bool firstEventSet = false;
	bool eventNumLitSet = false;
	
	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-hitdata") == 0)
		{
			hitdata = true;
		}
		else if (strcmp(argv[i], "-simdata") == 0)
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
					"AliHLTMUONTriggerRecordsSource::DoInit",
					"Missing parameter",
					"Expected one of 'left', 'right' or 'all' after '-plane'."
				);
				return EINVAL;
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
					"AliHLTMUONTriggerRecordsSource::DoInit",
					"Invalid parameter",
					"The parameter '%s' is invalid and must be one of 'left',"
					  " 'right' or 'all'.",
					argv[i]
				);
				return EINVAL;
			}
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
				return EINVAL;
			}
			char* end = "";
			long num = strtol(argv[i], &end, 0);
			if (*end != '\0' or num < 0) // Check if the conversion is OK.
			{
				HLTError(Form(
					"Expected a positive number after -firstevent"
					" but got: %s", argv[i]
				));
				return EINVAL;
			}
			fCurrentEvent = Int_t(num);
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
			fCurrentEvent = -1;
			eventNumLitSet = true;
		}
		else
		{
			Logging(kHLTLogError,
				"AliHLTMUONTriggerRecordsSource::DoInit",
				"Unknown argument",
				"The argument '%s' is invalid.",
				argv[i]
			);
			return EINVAL;
		}
	}

	// Check that one and only one of the the -hitdata, -simdata or
	// -recdata parameters was specified on the command line.
	if ((not hitdata and not simdata and not recdata) or
	    (not hitdata and simdata and recdata) or
	    (hitdata and not simdata and recdata) or
	    (hitdata and simdata and not recdata) or
	    (hitdata and simdata and recdata)
	   )
	{
		Logging(kHLTLogError,
			"AliHLTMUONTriggerRecordsSource::DoInit",
			"Missing arguments",
			"Must have one and only one of -hitdata, -simdata or -recdata specified."
		);
		return EINVAL;
	}
	
	// Must load the mapping data for AliMpTriggerCrate::GetDdlId()  //TODO AliMpTriggerCrate => AliMpDetElement
	// to return useful information later on.
	AliMpCDB::LoadDDLStore();
	
	// Now we can initialise the data interface objects and loaders.
	fBuildFromHits = hitdata;
	if (hitdata or simdata)
	{
		const char* message = fBuildFromHits ?
			"Loading simulated GEANT hits with AliMUONMCDataInterface."
			: "Loading simulated local trigger objects with AliMUONMCDataInterface.";
				
		Logging(kHLTLogDebug, "AliHLTMUONTriggerRecordsSource::DoInit",
			"Data interface", message
		);
		
		try
		{
			fMCDataInterface = new AliMUONMCDataInterface("galice.root");
		}
		catch (const std::bad_alloc&)
		{
			Logging(kHLTLogError,
				"AliHLTMUONTriggerRecordsSource::DoInit",
				"Out of memory",
				"Not enough memory to allocate AliMUONMCDataInterface."
			);
			return ENOMEM;
		}
	}
	else if (recdata)
	{
		Logging(kHLTLogDebug,
			"AliHLTMUONTriggerRecordsSource::DoInit",
			"Data interface",
			"Loading reconstructed local trigger objects with AliMUONDataInterface."
		);
		
		try
		{
			fDataInterface = new AliMUONDataInterface("galice.root");
		}
		catch (const std::bad_alloc&)
		{
			Logging(kHLTLogError,
				"AliHLTMUONTriggerRecordsSource::DoInit",
				"Out of memory",
				"Not enough memory to allocate AliMUONDataInterface."
			);
			return ENOMEM;
		}
	}
	
	// Check that the fCurrentEvent number falls within the correct range.
	UInt_t maxevent = 0;
	if (fMCDataInterface != NULL)
		maxevent = UInt_t(fMCDataInterface->NumberOfEvents());
	else if (fDataInterface != NULL)
		maxevent = UInt_t(fDataInterface->NumberOfEvents());
	if (fCurrentEvent != -1 and UInt_t(fCurrentEvent) >= maxevent and maxevent != 0)
	{
		fCurrentEvent = 0;
		HLTWarning(Form("The selected first event number (%d) was larger than"
			" the available number of events (%d). Resetting the event"
			" counter to zero.", fCurrentEvent, maxevent
		));
	}
	
	return 0;
}


int AliHLTMUONTriggerRecordsSource::DoDeinit()
{
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


const char* AliHLTMUONTriggerRecordsSource::GetComponentID()
{
	return AliHLTMUONConstants::TriggerRecordsSourceId();
}


AliHLTComponentDataType AliHLTMUONTriggerRecordsSource::GetOutputDataType()
{
	return AliHLTMUONConstants::TriggerRecordsBlockDataType();
}


void AliHLTMUONTriggerRecordsSource::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	constBase = sizeof(AliHLTMUONTriggerRecordsBlockStruct) +
		sizeof(AliHLTMUONTriggerRecordStruct) * AliMUONConstants::NTriggerCircuit();
	inputMultiplier = 0;
}


AliHLTComponent* AliHLTMUONTriggerRecordsSource::Spawn()
{
	return new AliHLTMUONTriggerRecordsSource();
}


int AliHLTMUONTriggerRecordsSource::GetEvent(
		const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& /*trigData*/,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& outputBlocks
	)
{
	assert( fMCDataInterface != NULL or fDataInterface != NULL );

	AliHLTInt32_t trigRecId = 0;

	// Check the size of the event descriptor structure.
	if (evtData.fStructSize < sizeof(AliHLTComponentEventData))
	{
		Logging(kHLTLogError,
			"AliHLTMUONTriggerRecordsSource::GetEvent",
			"Invalid event descriptor",
			"The event descriptor (AliHLTComponentEventData) size is"
			  " smaller than expected. It claims to be %d bytes, but"
			  " we expect it to be %d bytes.",
			evtData.fStructSize,
			sizeof(AliHLTComponentEventData)
		);
		size = 0; // Important to tell framework that nothing was generated.
		return EINVAL;
	}
	
	// Use the fEventID as the event number to load if fCurrentEvent == -1,
	// check it and load that event with the runloader.
	// If fCurrentEvent is a positive number then us it instead and
	// increment it.
	UInt_t eventnumber = UInt_t(evtData.fEventID);
	UInt_t maxevent = fMCDataInterface != NULL ?
		UInt_t(fMCDataInterface->NumberOfEvents())
		: UInt_t(fDataInterface->NumberOfEvents());
	if (fCurrentEvent != -1)
	{
		eventnumber = UInt_t(fCurrentEvent);
		fCurrentEvent++;
		if (UInt_t(fCurrentEvent) >= maxevent)
			fCurrentEvent = 0;
	}
	if ( eventnumber >= maxevent )
	{
		Logging(kHLTLogError,
			"AliHLTMUONTriggerRecordsSource::GetEvent",
			"Bad event ID",
			"The event number (%d) is larger than the available number"
			  " of events on file (%d).",
			eventnumber,
			maxevent
		);
		size = 0; // Important to tell framework that nothing was generated.
		return EINVAL;
	}
	
	// Create and initialise a new data block.
	AliHLTMUONTriggerRecordsBlockWriter block(outputPtr, size);
	if (not block.InitCommonHeader())
	{
		Logging(kHLTLogError,
			"AliHLTMUONTriggerRecordsSource::GetEvent",
			"Buffer too small",
			"There is not enough buffer space to create a new data block."
			  " We require at least %d bytes but the buffer is only %d bytes.",
			sizeof(AliHLTMUONTriggerRecordsBlockWriter::HeaderType),
			block.BufferSize()
		);
		size = 0; // Important to tell framework that nothing was generated.
		return ENOBUFS;
	}
	
	// Initialise the DDL list containing the DDLs which contributed to the
	// data block. These are required to create the specification word later.
	bool ddlList[22];
	for (Int_t i = 0; i < 22; i++)
		ddlList[i] = false;
	
	if (fMCDataInterface != NULL and fBuildFromHits)
	{
		Logging(kHLTLogDebug,
			"AliHLTMUONTriggerRecordsSource::GetEvent",
			"Filling triggers",
			"Filling data block with trigger records from GEANT hits for event %d.",
			eventnumber
		);
		
		// Loop over all tracks, extract the hits from chambers 11 to 14 and
		// create trigger records from them to write to the data block.
		Int_t ntracks = fMCDataInterface->NumberOfTracks(eventnumber);
		for (Int_t i = 0; i < ntracks; ++i)
		{
			AliMUONHit* hit11 = NULL;
			AliMUONHit* hit12 = NULL;
			AliMUONHit* hit13 = NULL;
			AliMUONHit* hit14 = NULL;
			Int_t ddl11 = -1;
			Int_t ddl12 = -1;
			Int_t ddl13 = -1;
			Int_t ddl14 = -1;
			
			AliMUONVHitStore* hitStore = fMCDataInterface->HitStore(eventnumber,i);
			AliMUONHit* hit;
			TIter next(hitStore->CreateIterator());
			while ( ( hit = static_cast<AliMUONHit*>(next()) ) )
			{
				// Select only hits on trigger chambers.
				if (hit->Chamber() <= AliMUONConstants::NTrackingCh()) continue;
				
				// Only select hits from the given part of the plane
				if (fSelection == kLeftPlane and not (hit->Xref() < 0)) continue;
				if (fSelection == kRightPlane and not (hit->Xref() >= 0)) continue;
				
				// Workout which DDL this hit should be readout of.
				Int_t ddl = FindDDLOfDetElement(hit->DetElemId());
				if (not (0 <= ddl and ddl < 22))
				{
					ddl = -1;
					Logging(kHLTLogError,
						"AliHLTMUONTriggerRecordsSource::GetEvent",
						"No DDL ID",
						"Could not find the DDL ID from which readout would take place."
					);
				}
				
				switch (hit->Chamber())
				{
				case 11: hit11 = hit; ddl11 = ddl; break;
				case 12: hit12 = hit; ddl12 = ddl; break;
				case 13: hit13 = hit; ddl13 = ddl; break;
				case 14: hit14 = hit; ddl14 = ddl; break;
				default: break;
				}
			}
			
			// Check that there are at least 3 of 4 hits on the trigger chambers.
			Int_t hitCount = 0;
			if (hit11 != NULL) hitCount++;
			if (hit12 != NULL) hitCount++;
			if (hit13 != NULL) hitCount++;
			if (hit14 != NULL) hitCount++;
			if (hitCount < 3) continue;
				
			AliHLTMUONTriggerRecordStruct* trigRec = block.AddEntry();
			if (trigRec == NULL)
			{
				Logging(kHLTLogError,
					"AliHLTMUONTriggerRecordsSource::GetEvent",
					"Buffer overflow",
					"There is not enough buffer space to add more trigger records."
					  " We overflowed the buffer which is only %d bytes.",
					block.BufferSize()
				);
				size = 0; // Important to tell framework that nothing was generated.
				return ENOBUFS;
			}
			
			// Fill the new trigger record with the hit information.
			bool hitset[4] = {false, false, false, false};
			AliHLTFloat32_t x1 = 0, y1 = 0, y2 = 0, z1 = 0, z2 = 0;
			if (hit11 != NULL)
			{
				trigRec->fHit[0].fX = hit11->Xref();
				trigRec->fHit[0].fY = hit11->Yref();
				trigRec->fHit[0].fZ = hit11->Zref();
				hitset[0] = true;
				x1 = hit11->Xref();
				y1 = hit11->Yref();
				z1 = hit11->Zref();
			}
			if (hit12 != NULL)
			{
				trigRec->fHit[1].fX = hit12->Xref();
				trigRec->fHit[1].fY = hit12->Yref();
				trigRec->fHit[1].fZ = hit12->Zref();
				hitset[1] = true;
				x1 = hit12->Xref();
				y1 = hit12->Yref();
				z1 = hit12->Zref();
			}
			if (hit13 != NULL)
			{
				trigRec->fHit[2].fX = hit13->Xref();
				trigRec->fHit[2].fY = hit13->Yref();
				trigRec->fHit[2].fZ = hit13->Zref();
				hitset[2] = true;
				y2 = hit13->Yref();
				z2 = hit13->Zref();
			}
			if (hit14 != NULL)
			{
				trigRec->fHit[3].fX = hit14->Xref();
				trigRec->fHit[3].fY = hit14->Yref();
				trigRec->fHit[3].fZ = hit14->Zref();
				hitset[3] = true;
				y2 = hit14->Yref();
				z2 = hit14->Zref();
			}
			
			bool calculated = AliHLTMUONCalculations::ComputeMomentum(x1, y1, y2, z1, z2);
			if (not calculated)
				Logging(kHLTLogDebug,
					"AliHLTMUONTriggerRecordsSource::GetEvent",
					"Calculation failure",
					"Something went wrong when calculating the momentum from"
					  " x1 = %f, y1 = %f, y2 = %f, z1 = %f, z2 = %f.",
					x1, y1, y2, z1, z2
				);
			
			trigRec->fId = trigRecId++;
			trigRec->fFlags = AliHLTMUONUtils::PackTriggerRecordFlags(
					AliHLTMUONCalculations::Sign(), hitset
				);
			trigRec->fPx = AliHLTMUONCalculations::Px();
			trigRec->fPy = AliHLTMUONCalculations::Py();
			trigRec->fPz = AliHLTMUONCalculations::Pz();
			
			// Mark the DDLs over which this trigger record would be readout.
			if (ddl11 != -1) ddlList[ddl11] = true;
			if (ddl12 != -1) ddlList[ddl12] = true;
			if (ddl13 != -1) ddlList[ddl13] = true;
			if (ddl14 != -1) ddlList[ddl14] = true;
		}
	}
	else if (fMCDataInterface != NULL and not fBuildFromHits)
	{
		Logging(kHLTLogDebug,
			"AliHLTMUONTriggerRecordsSource::GetEvent",
			"Filling triggers",
			"Filling data block with simulated local triggers for event %d.",
			eventnumber
		);
		
		AliFatal("Sorry, -simdata option not yet implemented!");
		// TODO
	}
	else if (fDataInterface != NULL)
	{
		Logging(kHLTLogDebug,
			"AliHLTMUONTriggerRecordsSource::GetEvent",
			"Filling triggers",
			"Filling data block with reconstructed local triggers for event %d.",
			eventnumber
		);
		// TODO
		AliFatal("Sorry, -recdata option not yet implemented!");
	}
	else
	{
		Logging(kHLTLogError,
			"AliHLTMUONTriggerRecordsSource::GetEvent",
			"Missing data interface",
			"Neither AliMUONDataInterface nor AliMUONMCDataInterface were created."
		);
		size = 0; // Important to tell framework that nothing was generated.
		return EFAULT;
	}
	
	AliHLTComponentBlockData bd;
	FillBlockData(bd);
	bd.fPtr = outputPtr;
	bd.fOffset = 0;
	bd.fSize = block.BytesUsed();
	bd.fDataType = AliHLTMUONConstants::TriggerRecordsBlockDataType();
	bd.fSpecification = AliHLTMUONUtils::PackSpecBits(ddlList);
	outputBlocks.push_back(bd);
	size = block.BytesUsed();

	return 0;
}
