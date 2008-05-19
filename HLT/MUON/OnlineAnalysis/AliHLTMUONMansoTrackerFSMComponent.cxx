/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
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
///  @file   AliHLTMUONMansoTrackerFSMComponent.cxx
///  @author Artur Szostak <artursz@iafrica.com>,
///          Indranil Das <indra.das@saha.ac.in>
///  @date   
///  @brief  Implementation of AliHLTMUONMansoTrackerFSMComponent class.
///

#include "AliHLTMUONMansoTrackerFSMComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONMansoTrackerFSM.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONDataBlockWriter.h"
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <new>

ClassImp(AliHLTMUONMansoTrackerFSMComponent);


AliHLTMUONMansoTrackerFSMComponent::AliHLTMUONMansoTrackerFSMComponent() :
	AliHLTMUONProcessor(),
	AliHLTMUONMansoTrackerFSMCallback(),
	fTracker(NULL),
	fTrackCount(0),
	fBlock(NULL),
	fRecHitBlockArraySize(0),
	fWarnForUnexpecedBlock(false)
{
	///
	/// Default constructor.
	///
	
	for (Int_t i = 0; i < 4; i++)
	{
		fRecHitBlockCount[i] = 0;
		fRecHitBlock[i] = NULL;
	}
}


AliHLTMUONMansoTrackerFSMComponent::~AliHLTMUONMansoTrackerFSMComponent()
{
	///
	/// Default destructor.
	///
	
	// Should never have the following 2 pointers non-NULL since DoDeinit
	// should have been called before, but handle this case anyway.
	if (fTracker != NULL) delete fTracker;
	
	// Remember that only fRecHitBlock[0] stores the pointer to the allocated
	// memory. The other pointers are just reletive to this.
	if (fRecHitBlock[0] != NULL) delete [] fRecHitBlock[0];
}


const char* AliHLTMUONMansoTrackerFSMComponent::GetComponentID()
{
	///
	/// Inherited from AliHLTComponent. Returns the component ID.
	///
	
	return AliHLTMUONConstants::MansoTrackerFSMId();
}


void AliHLTMUONMansoTrackerFSMComponent::GetInputDataTypes(
		vector<AliHLTComponentDataType>& list
	)
{
	///
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	///
	
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::TriggerRecordsBlockDataType() );
	list.push_back( AliHLTMUONConstants::RecHitsBlockDataType() );
}


AliHLTComponentDataType AliHLTMUONMansoTrackerFSMComponent::GetOutputDataType()
{
	///
	/// Inherited from AliHLTComponent. Returns the output data type.
	///
	
	return AliHLTMUONConstants::MansoTracksBlockDataType();
}


void AliHLTMUONMansoTrackerFSMComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	///
	/// Inherited from AliHLTComponent. Returns an estimate of the expected output data size.
	///
	
	constBase = sizeof(AliHLTMUONMansoTracksBlockStruct);
	inputMultiplier = 1;
}


AliHLTComponent* AliHLTMUONMansoTrackerFSMComponent::Spawn()
{
	///
	/// Inherited from AliHLTComponent. Creates a new object instance.
	///
	
	return new AliHLTMUONMansoTrackerFSMComponent;
}


int AliHLTMUONMansoTrackerFSMComponent::DoInit(int argc, const char** argv)
{
	///
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	///
	
	HLTInfo("Initialising dHLT manso tracker FSM component.");
	
	// Just in case for whatever reason we still have some of the internal
	// object allocated previously still hanging around delete them now.
	FreeMemory();
	
	try
	{
		fTracker = new AliHLTMUONMansoTrackerFSM();
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the tracker component.");
		return -ENOMEM;
	}
	fTracker->SetCallback(this);
	
	fWarnForUnexpecedBlock = false;
	
	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-warn_on_unexpected_block") == 0)
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}

		HLTError("Unknown option '%s'.", argv[i]);
		FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
		return -EINVAL;
	}
	
	const int initArraySize = 10;
	// Allocate some initial memory for the reconstructed hit arrays.
	try
	{
		fRecHitBlock[0] = new AliRecHitBlockInfo[initArraySize*4];
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for the reconstructed hit arrays.");
		FreeMemory(); // Make sure we cleanup to avoid partial initialisation.
		return -ENOMEM;
	}
	// Only set the arrays' size once we have successfully allocated the memory for the arrays.
	fRecHitBlockArraySize = initArraySize;
	// Now we need to set the pointers fRecHitBlock[i] {i>0} relative to fRecHitBlock[0].
	for (Int_t i = 1; i < 4; i++)
	{
		fRecHitBlock[i] = fRecHitBlock[i-1] + fRecHitBlockArraySize;
	}
	// And reset the number of records actually stored in the arrays.
	for (Int_t i = 0; i < 4; i++)
	{
		fRecHitBlockCount[i] = 0;
	}
	
	return 0;
}


int AliHLTMUONMansoTrackerFSMComponent::DoDeinit()
{
	///
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	///
	
	HLTInfo("Deinitialising dHLT manso tracker FSM component.");
	FreeMemory();
	return 0;
}


int AliHLTMUONMansoTrackerFSMComponent::DoEvent(
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
	
	Reset();
	AliHLTUInt32_t specification = 0;  // Contains the output data block spec bits.
	
	// Resize the rec hit arrays if we possibly will need more space.
	// To guarantee that they will not overflow we need to make sure each
	// array is at least as big as the number of input data blocks.
	if (fRecHitBlockArraySize < evtData.fBlockCnt)
	{
		// Release the old memory block and allocate more memory.
		if (fRecHitBlock[0] != NULL)
		{
			delete [] fRecHitBlock[0];
		}
		
		// Reset the number of records actually stored in the arrays.
		for (Int_t i = 0; i < 4; i++)
		{
			fRecHitBlockCount[i] = 0;
		}
		
		try
		{
			fRecHitBlock[0] = new AliRecHitBlockInfo[evtData.fBlockCnt*4];
		}
		catch (const std::bad_alloc&)
		{
			HLTError("Could not allocate more memory for the reconstructed hit arrays.");
			// Ok so now we need to clear all the pointers because we actually
			// deleted the memory.
			fRecHitBlockArraySize = 0;
			for (Int_t i = 0; i < 4; i++)
			{
				fRecHitBlock[i] = NULL;
			}
			return -ENOMEM;
		}
		// Only set the arrays' size once we have successfully allocated the memory for the arrays.
		fRecHitBlockArraySize = evtData.fBlockCnt;
		// Now we need to set the pointers fRecHitBlock[i] {i>0} relative to fRecHitBlock[0].
		for (Int_t i = 1; i < 4; i++)
		{
			fRecHitBlock[i] = fRecHitBlock[i-1] + fRecHitBlockArraySize;
		}
	}
	
	AliHLTMUONMansoTracksBlockWriter block(outputPtr, size);
	fBlock = &block;
	
	if (not block.InitCommonHeader())
	{
		Logging(kHLTLogError,
			"AliHLTMUONMansoTrackerFSMComponent::DoEvent",
			"Buffer overflow",
			"The buffer is only %d bytes in size. We need a minimum of %d bytes.",
			size, sizeof(AliHLTMUONMansoTracksBlockWriter::HeaderType)
		);
		size = 0; // Important to tell framework that nothing was generated.
		return -ENOBUFS;
	}

	// Loop over all input blocks in the event and add the ones that contain
	// reconstructed hits into the hit buffers. The blocks containing trigger
	// records are ignored for now and will be processed later.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
			n, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fPtr, blocks[n].fSize
		);
		
		if (blocks[n].fDataType == AliHLTMUONConstants::RecHitsBlockDataType())
		{
			specification |= blocks[n].fSpecification;
			
			AliHLTMUONRecHitsBlockReader inblock(blocks[n].fPtr, blocks[n].fSize);
			if (not BlockStructureOk(inblock)) continue;
			
			if (inblock.Nentries() != 0)
				AddRecHits(blocks[n].fSpecification, inblock.GetArray(), inblock.Nentries());
			else
			{
				Logging(kHLTLogDebug,
					"AliHLTMUONMansoTrackerFSMComponent::DoEvent",
					"Block empty",
					"Received a reconstructed hits data block which contains no entries."
				);
			}
		}
		else if (blocks[n].fDataType != AliHLTMUONConstants::TriggerRecordsBlockDataType())
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
		}
	}
  
	// Again loop over all input blocks in the event, but this time look for
	// the trigger record blocks and process these.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		if (blocks[n].fDataType != AliHLTMUONConstants::TriggerRecordsBlockDataType())
			continue;
		
		AliHLTMUONTriggerRecordsBlockReader inblock(blocks[n].fPtr, blocks[n].fSize);
		if (not BlockStructureOk(inblock)) continue;
		
		DebugTrace("Processing a trigger block with "
			<< inblock.Nentries() << " entries."
		);
		
		specification |= blocks[n].fSpecification;
		
		for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
		{
			fTracker->FindTrack(inblock[i]);
			
			// Reset the tracker so that we do not double count tracks.
			fTracker->Reset();
		}
	}
	
	AliHLTComponentBlockData bd;
	FillBlockData(bd);
	bd.fPtr = outputPtr;
	bd.fOffset = 0;
	bd.fSize = block.BytesUsed();
	bd.fDataType = AliHLTMUONConstants::MansoTracksBlockDataType();
	bd.fSpecification = specification;
	outputBlocks.push_back(bd);
	size = block.BytesUsed();

	return 0;
}


void AliHLTMUONMansoTrackerFSMComponent::Reset()
{
	///
	/// Reset the track count and reconstructed hit data block arrays.
	///
	
	DebugTrace("Resetting AliHLTMUONMansoTrackerFSMComponent.");

	//fTracker->Reset();  // Not necessary here because it is done after every FindTrack call.
	fTrackCount = 0;
	fBlock = NULL;  // Do not delete. Already done implicitly at the end of DoEvent.
	for (int i = 0; i < 4; i++)
	{
		fRecHitBlockCount[i] = 0;
	}
}


void AliHLTMUONMansoTrackerFSMComponent::FreeMemory()
{
	/// Deletes any objects and arrays allocated by this component and releases
	/// the memory used. This is called as a helper routine by the init and deinit
	/// methods. If some or all of the object pointers are already NULL then
	/// nothing is done for those. This method guarantees that all the relevant
	/// pointers will be NULL after returning from this method.

	if (fTracker != NULL)
	{
		delete fTracker;
		fTracker = NULL;
	}
	
	// Remember that only fRecHitBlock[0] stores the pointer to the allocated memory.
	// The other pointers are just reletive to this.
	if (fRecHitBlock[0] != NULL)
		delete [] fRecHitBlock[0];
	
	fRecHitBlockArraySize = 0;
	for (Int_t i = 0; i < 4; i++)
	{
		fRecHitBlockCount[i] = 0;
		fRecHitBlock[i] = NULL;
	}
}


void AliHLTMUONMansoTrackerFSMComponent::AddRecHits(
		AliHLTUInt32_t specification,
		const AliHLTMUONRecHitStruct* recHits,
		AliHLTUInt32_t count
	)
{
	///
	/// Adds a new reconstructed hit data block to the internal list of blocks
	/// for the tracker to process.
	/// These lists will later be used when the tracker requests them through
	/// the callback method 'RequestClusters'.
	///
	
	DebugTrace("AliHLTMUONMansoTrackerFSMComponent::AddRecHits called with specification = 0x"
		 << std::hex << specification << std::dec << " and count = "
		 << count << " rec hits."
	);
	
	AliHLTUInt8_t chamberMap[20] = {
		1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10
	};
	
	// Identify the chamber the rec hits came from using the specifications field.
	bool gotDataFromDDL[22];
	AliHLTMUONUtils::UnpackSpecBits(specification, gotDataFromDDL);
		
	AliHLTInt8_t chamber = -1;
	for (int i = 0; i < 20; i++)
	{
		if (not gotDataFromDDL[i]) continue;
		if (7 <= chamberMap[i] and chamberMap[i] <= 10)
		{
			if (chamber != -1 and chamber != chamberMap[i])
			{
				Logging(kHLTLogError,
					"AliHLTMUONMansoTrackerFSMComponent::AddRecHits",
					"Invalid block",
					"Received a data block with data from multiple chambers."
					  " This component cannot handle such a case."
				);
				return;
			}
			else
				chamber = chamberMap[i];
		}
		else
		{
			Logging(kHLTLogError,
				"AliHLTMUONMansoTrackerFSMComponent::AddRecHits",
				"Invalid chamber",
				"Received a data block with data from chamber %d"
				  " which is outside the expected range: [7..10].",
				chamberMap[i]
			);
			return;
		}
	}
	
	// Make sure we got one chamber number.
	if (chamber < 7 or 10 < chamber)
	{
		Logging(kHLTLogError,
			"AliHLTMUONMansoTrackerFSMComponent::AddRecHits",
			"Invalid block",
			"Received a reconstructed hit data block with a null specification."
			 " Cannot know which chamber the data comes from."
		);
		return;
	}
	
	DebugTrace("Added " << count << " reconstructed hits from chamber "
		<< (int)chamber	<< " to the internal arrays."
	);
	
	assert( fRecHitBlockCount[chamber-7] < fRecHitBlockArraySize );
	AliRecHitBlockInfo info(count, recHits);
	fRecHitBlock[chamber-7][fRecHitBlockCount[chamber-7]] = info;
	fRecHitBlockCount[chamber-7]++;
}


void AliHLTMUONMansoTrackerFSMComponent::RequestClusters(
		AliHLTMUONMansoTrackerFSM* tracker,
		AliHLTFloat32_t left, AliHLTFloat32_t right,
		AliHLTFloat32_t bottom, AliHLTFloat32_t top,
		AliHLTMUONChamberName chamber, const void* tag
	)
{
	///
	/// Inherited from AliHLTMUONMansoTrackerFSMCallback.
	/// This is the call back method used by the tracker algorithm to request
	/// clusters on a certain chamber.
	///

	DebugTrace("AliHLTMUONMansoTracker::RequestClusters(chamber = " << chamber << ")");
	void* ctag = const_cast<void*>(tag);
	int chNo = -1;
	AliHLTUInt32_t recHitsCount = 0;
	AliRecHitBlockInfo* recHitsBlock = NULL;
	switch (chamber)
	{
	case kChamber7:
		recHitsCount = fRecHitBlockCount[0];
		recHitsBlock = fRecHitBlock[0];
		chNo = 7;
		break;

	case kChamber8:
		recHitsCount = fRecHitBlockCount[1];
		recHitsBlock = fRecHitBlock[1];
		chNo = 8;
		break;

	case kChamber9:
		recHitsCount = fRecHitBlockCount[2];
		recHitsBlock = fRecHitBlock[2];
		chNo = 9;
		break;

	case kChamber10:
		recHitsCount = fRecHitBlockCount[3];
		recHitsBlock = fRecHitBlock[3];
		chNo = 10;
		break;

	default: return;
	}
	
	DebugTrace("Returning requested hits for chamber " << chNo << ":");
	for (AliHLTUInt32_t i = 0; i < recHitsCount; i++)
	for (AliHLTUInt32_t j = 0; j < recHitsBlock[i].Count(); j++)
	{
		const AliHLTMUONRecHitStruct* hit = &(recHitsBlock[i].Data()[j]);
		if (left < hit->fX and hit->fX < right and bottom < hit->fY and hit->fY < top)
			tracker->ReturnClusters(ctag, hit, 1);
	}
	DebugTrace("Done returning hits from chamber " << chNo << ".");
	tracker->EndOfClusters(ctag);
}


void AliHLTMUONMansoTrackerFSMComponent::EndOfClusterRequests(
		AliHLTMUONMansoTrackerFSM* /*tracker*/
	)
{
	///
	/// Inherited from AliHLTMUONMansoTrackerFSMCallback.
	/// Nothing special to do here.
	///
	
	DebugTrace("End of cluster requests.");
}


void AliHLTMUONMansoTrackerFSMComponent::FoundTrack(AliHLTMUONMansoTrackerFSM* tracker)
{
	///
	/// Inherited from AliHLTMUONMansoTrackerFSMCallback.
	/// This is the call back method used by the tracker algorithm to declare
	/// that a new track has been found.
	///
	
	DebugTrace("AliHLTMUONMansoTrackerFSMComponent::FoundTrack()");
	
	AliHLTMUONMansoTracksBlockWriter* block =
		reinterpret_cast<AliHLTMUONMansoTracksBlockWriter*>(fBlock);
	
	AliHLTMUONMansoTrackStruct* track = block->AddEntry();
	if (track == NULL)
	{
		Logging(kHLTLogError,
			"AliHLTMUONMansoTrackerFSMComponent::FoundTrack",
			"Buffer overflow",
			"We have overflowed the output buffer for Manso track data."
			  " The output buffer size is only %d bytes.",
			block->BufferSize()
		);
		return;
	}
 
	fTrackCount++;
	tracker->FillTrackData(*track);
	DebugTrace("\tTrack data = " << *track);
}


void AliHLTMUONMansoTrackerFSMComponent::NoTrackFound(AliHLTMUONMansoTrackerFSM* /*tracker*/)
{
	///
	/// Inherited from AliHLTMUONMansoTrackerFSMCallback.
	/// Nothing special to do here.
	///
	
	DebugTrace("No track found.");
}

