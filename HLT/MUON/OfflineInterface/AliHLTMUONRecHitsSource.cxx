/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
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
 * @file   AliHLTMUONRecHitsSource.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of the AliHLTMUONRecHitsSource component.
 */

#include "AliHLTMUONRecHitsSource.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliMUONSimData.h"
#include "AliMUONRecData.h"
#include "AliMUONHit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONConstants.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "TClonesArray.h"
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <cassert>
#include <new>

#include <iostream>
using namespace std;


//namespace
//{
	// The global object used for automatic component registration.
	// Note DO NOT use this component for calculation!
	AliHLTMUONRecHitsSource gAliHLTMUONRecHitsSource;
//}


ClassImp(AliHLTMUONRecHitsSource);


AliHLTMUONRecHitsSource::AliHLTMUONRecHitsSource() :
	AliHLTOfflineDataSource(),
	fSimData(NULL), fRecData(NULL),
	fRunLoader(NULL), fLoader(NULL)
{
}

AliHLTMUONRecHitsSource::~AliHLTMUONRecHitsSource()
{
}


int AliHLTMUONRecHitsSource::DoInit(int argc, const char** argv)
{
	// Parse the command line arguments:
	bool simdata = false;
	bool recdata = false;
	int i = 0;
	while (i < argc)
	{
		if (strcmp(argv[i], "-simdata") != 0)
			simdata = true;
		else if (strcmp(argv[i], "-recdata") != 0)
			recdata = true;
		else
		{
			Logging(kHLTLogError,
				"AliHLTMUONRecHitsSource::DoInit",
				"Unknown argument",
				"The argument '%s' is invalid.",
				argv[i]
			);
			return EINVAL;
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
		return EINVAL;
	}
	
	if (not simdata and not recdata)
	{
		Logging(kHLTLogError,
			"AliHLTMUONRecHitsSource::DoInit",
			"Missing arguments",
			"Must have either -simdata or -recdata specified."
		);
		return EINVAL;
	}
	
	// Now we can initialise the data interface objects and loaders.
	if (simdata)
	{
		try
		{
			fSimData = new AliMUONSimData("galice.root");
		}
		catch (const std::bad_alloc&)
		{
			Logging(kHLTLogError,
				"AliHLTMUONRecHitsSource::DoInit",
				"Out of memory",
				"Not enough memory to allocate AliMUONSimData."
			);
			return ENOMEM;
		}
		fLoader = fSimData->GetLoader();
		fLoader->LoadHits("READ");
	}
	else if (recdata)
	{
		try
		{
			fRecData = new AliMUONRecData("galice.root");
		}
		catch (const std::bad_alloc&)
		{
			Logging(kHLTLogError,
				"AliHLTMUONRecHitsSource::DoInit",
				"Out of memory",
				"Not enough memory to allocate AliMUONRecData."
			);
			return ENOMEM;
		}
		fLoader = fRecData->GetLoader();
		fLoader->LoadRecPoints("READ");
	}
	
	fRunLoader = AliRunLoader::GetRunLoader();
	
	return 0;
}


int AliHLTMUONRecHitsSource::DoDeinit()
{
	if (fSimData != NULL)
	{
		fLoader->UnloadHits();
		delete fSimData;
		fSimData = NULL;
	}
	if (fRecData != NULL)
	{
		fLoader->UnloadRecPoints();
		delete fRecData;
		fRecData = NULL;
	}
	fRunLoader = NULL;
	fLoader = NULL;
	return 0;
}


const char* AliHLTMUONRecHitsSource::GetComponentID()
{
	return AliHLTMUONConstants::RecHitsSourceId();
}

AliHLTComponentDataType AliHLTMUONRecHitsSource::GetOutputDataType()
{
	return AliHLTMUONConstants::RecHitsBlockDataType();
}

void AliHLTMUONRecHitsSource::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	constBase = sizeof(AliHLTMUONRecHitsBlockStruct);
	inputMultiplier = 1;
}

AliHLTComponent* AliHLTMUONRecHitsSource::Spawn()
{
	return new AliHLTMUONRecHitsSource();
}


int AliHLTMUONRecHitsSource::GetEvent(
		const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		vector<AliHLTComponentBlockData>& outputBlocks
	)
{
	assert( fSimData != NULL or fRecData != NULL );
	assert( fRunLoader != NULL );
	assert( fLoader != NULL );

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
		return EINVAL;
	}
	
	// Use the fEventID as the event number to load, check it and load that
	// event with the runloader.
	Int_t eventnumber = Int_t(evtData.fEventID);
	if (eventnumber >= fRunLoader->GetNumberOfEvents())
	{
		Logging(kHLTLogError,
			"AliHLTMUONRecHitsSource::GetEvent",
			"Bad event ID",
			"The event number (%d) is larger than the available number"
			  " of events on file (%d).",
			eventnumber,
			fRunLoader->GetNumberOfEvents()
		);
		size = 0; // Important to tell framework that nothing was generated.
		return EINVAL;
	}
	fRunLoader->GetEvent(eventnumber);
	
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
		return ENOBUFS;
	}
	
	if (fSimData != NULL)
	{
		// Loop over all tracks, extract the hits and write them to the
		// data block.
		fSimData->SetTreeAddress("H");
		for (Int_t i = 0; i < fSimData->GetNtracks(); i++)
		{
			fSimData->GetTrack(i);
			assert( fSimData->Hits() != NULL );
			Int_t nhits = fSimData->Hits()->GetEntriesFast();
			for (Int_t j = 0; j < nhits; j++)
			{
				AliMUONHit* hit = static_cast<AliMUONHit*>(
						fSimData->Hits()->At(j)
					);
				hit->Print();
				
				// Select only hits on a certain chamber.
				if (hit->Chamber() != 7) continue;
				
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
					fSimData->ResetHits();
					size = 0; // Important to tell framework that nothing was generated.
					return ENOBUFS;
				}
				
				rechit->fX = hit->Xref();
				rechit->fY = hit->Yref();
				rechit->fZ = hit->Zref();
			}
			fSimData->ResetHits();
		}
	}
	else if (fRecData != NULL)
	{
		//Int_t nchambers = AliMUONConstants::NTrackingCh();
		fRecData->SetTreeAddress("RC,TC"); 
		fRecData->GetRawClusters();
		
		// Select a specific chamber.
		Int_t chamber = 7;
		char branchname[32];
		sprintf(branchname, "MUONRawClusters%d", chamber);

		TClonesArray* clusterarray = fRecData->RawClusters(chamber);
		Int_t nrecpoints = clusterarray->GetEntriesFast();
		for (Int_t i = 0; i < nrecpoints; i++)
		{
			AliMUONRawCluster* cluster = static_cast<AliMUONRawCluster*>(clusterarray->At(i));
			cluster->GetX();
			
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
				fRecData->ResetRawClusters();
				size = 0; // Important to tell framework that nothing was generated.
				return ENOBUFS;
			}
			
			rechit->fX = cluster->GetX();
			rechit->fY = cluster->GetY();
			rechit->fZ = cluster->GetZ();
		}
		
		fRecData->ResetRawClusters();
	}
	else
	{
		Logging(kHLTLogError,
			"AliHLTMUONRecHitsSource::GetEvent",
			"Missing data interface",
			"Neither AliMUONSimData or AliMUONRecData were created."
		);
		size = 0; // Important to tell framework that nothing was generated.
		return EFAULT;
	}

	return 0;
}
