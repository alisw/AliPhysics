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

/**
 *  @file   AliHLTMUONMansoTrackerFSMComponent.cxx
 *  @author Artur Szostak <artursz@iafrica.com>,
 *          Indranil Das <indra.das@saha.ac.in>
 *  @date   
 *  @brief  Implementation of AliHLTMUONMansoTrackerFSMComponent class.
 */

#include "AliHLTMUONMansoTrackerFSMComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONMansoTrackerFSM.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONDataBlockWriter.h"
#include <cstdlib>
#include <cstring>
#include <cerrno>

namespace
{
	// The global object used for automatic component registration.
	// Note DO NOT use this component for calculation!
	AliHLTMUONMansoTrackerFSMComponent gAliHLTMUONMansoTrackerFSMComponent;

} // end of namespace


ClassImp(AliHLTMUONMansoTrackerFSMComponent);


AliHLTMUONMansoTrackerFSMComponent::AliHLTMUONMansoTrackerFSMComponent() :
	AliHLTProcessor(),
	AliHLTMUONMansoTrackerFSMCallback(),
	fTracker(NULL),
	fTrackCount(0),
	fBlock(NULL),
	fWarnForUnexpecedBlock(false)
{
}


AliHLTMUONMansoTrackerFSMComponent::~AliHLTMUONMansoTrackerFSMComponent()
{
	assert( fTracker == NULL );
}


const char* AliHLTMUONMansoTrackerFSMComponent::GetComponentID()
{
	return AliHLTMUONConstants::MansoTrackerFSMId();
}


void AliHLTMUONMansoTrackerFSMComponent::GetInputDataTypes(
		vector<AliHLTComponentDataType>& list
	)
{
	assert( list.empty() );
	list.push_back( AliHLTMUONConstants::TriggerRecordsBlockDataType() );
	list.push_back( AliHLTMUONConstants::RecHitsBlockDataType() );
}


AliHLTComponentDataType AliHLTMUONMansoTrackerFSMComponent::GetOutputDataType()
{
	return AliHLTMUONConstants::MansoTracksBlockDataType();
}


void AliHLTMUONMansoTrackerFSMComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	constBase = sizeof(AliHLTMUONMansoTracksBlockStruct);
	inputMultiplier = 1;
}


AliHLTComponent* AliHLTMUONMansoTrackerFSMComponent::Spawn()
{
	return new AliHLTMUONMansoTrackerFSMComponent;
}


int AliHLTMUONMansoTrackerFSMComponent::DoInit(int argc, const char** argv)
{
	fTracker = new AliHLTMUONMansoTrackerFSM();
	fTracker->SetCallback(this);
	
	fWarnForUnexpecedBlock = false;
	
	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-warn_on_unexpected_block") == 0)
			fWarnForUnexpecedBlock = true;
	}
	
	return 0;
}


int AliHLTMUONMansoTrackerFSMComponent::DoDeinit()
{
	if (fTracker != NULL)
	{
		delete fTracker;
		fTracker = NULL;
	}
	return 0;
}


int AliHLTMUONMansoTrackerFSMComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks, 
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr, 
		AliHLTUInt32_t& size,
		std::vector<AliHLTComponentBlockData>& outputBlocks
	)
{
	Reset();
	AliHLTUInt32_t specification = 0;  // Contains the output data block spec bits.
	
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
		return ENOBUFS;
	}

	// Loop over all input blocks in the event and add the ones that contain
	// reconstructed hits into the hit buffers. The blocks containing trigger
	// records are ignored for now and will be processed later.
	for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
	{
		if (blocks[n].fDataType == AliHLTMUONConstants::RecHitsBlockDataType())
		{
			specification |= blocks[n].fSpecification;
			
			AliHLTMUONRecHitsBlockReader inblock(blocks[n].fPtr, blocks[n].fSize);
			if (not inblock.BufferSizeOk())
			{
				size_t headerSize = sizeof(AliHLTMUONRecHitsBlockReader::HeaderType);
				if (blocks[n].fSize < headerSize)
				{
					HLTError("Received a reconstructed hits data block with a size of %d bytes,"
						" which is smaller than the minimum valid header size of %d bytes."
						" The block must be corrupt.",
						blocks[n].fSize, headerSize
					);
					continue;
				}
				
				size_t expectedWidth = sizeof(AliHLTMUONRecHitsBlockReader::ElementType);
				if (inblock.CommonBlockHeader().fRecordWidth != expectedWidth)
				{
					HLTError("Received a reconstructed hits data block with a record"
						" width of %d bytes, but the expected value is %d bytes."
						" The block might be corrupt.",
						blocks[n].fSize, headerSize
					);
					continue;
				}
				
				HLTError("Received a reconstructed hits data block with a size of %d bytes,"
					" but the block header claims the block should be %d bytes."
					" The block might be corrupt.",
					blocks[n].fSize, inblock.BytesUsed()
				);
				continue;
			}
			
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
			char id[kAliHLTComponentDataTypefIDsize+1];
			for (int i = 0; i < kAliHLTComponentDataTypefIDsize; i++)
				id[i] = blocks[n].fDataType.fID[i];
			id[kAliHLTComponentDataTypefIDsize] = '\0';
			char origin[kAliHLTComponentDataTypefOriginSize+1];
			for (int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++)
				origin[i] = blocks[n].fDataType.fOrigin[i];
			origin[kAliHLTComponentDataTypefOriginSize] = '\0';
			
			if (fWarnForUnexpecedBlock)
				HLTWarning("Received a data block of a type we cannot handle: %s origin: %s",
					static_cast<char*>(id), static_cast<char*>(origin)
				);
			else
				HLTDebug("Received a data block of a type we cannot handle: %s origin: %s",
					static_cast<char*>(id), static_cast<char*>(origin)
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
		if (not inblock.BufferSizeOk())
		{
			size_t headerSize = sizeof(AliHLTMUONTriggerRecordsBlockReader::HeaderType);
			if (blocks[n].fSize < headerSize)
			{
				HLTError("Received a trigger records data block with a size of %d bytes,"
					" which is smaller than the minimum valid header size of %d bytes."
					" The block must be corrupt.",
					blocks[n].fSize, headerSize
				);
				continue;
			}
			
			size_t expectedWidth = sizeof(AliHLTMUONTriggerRecordsBlockReader::ElementType);
			if (inblock.CommonBlockHeader().fRecordWidth != expectedWidth)
			{
				HLTError("Received a trigger records data block with a record"
					" width of %d bytes, but the expected value is %d bytes."
					" The block might be corrupt.",
					blocks[n].fSize, headerSize
				);
				continue;
			}
			
			HLTError("Received a trigger records data block with a size of %d bytes,"
				" but the block header claims the block should be %d bytes."
				" The block might be corrupt.",
				blocks[n].fSize, inblock.BytesUsed()
			);
			continue;
		}
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
	DebugTrace("Resetting AliHLTMUONMansoTrackerFSMComponent.");

	//fTracker->Reset();  // Not necessary here because it is done after every FindTrack call.
	fTrackCount = 0;
	fBlock = NULL;  // Do not delete. Already done implicitly at the end of DoEvent.
	for (int i = 0; i < 4; i++)
	{
		fRecHitBlock[i].erase(fRecHitBlock[i].begin(), fRecHitBlock[i].end());
	}
}


void AliHLTMUONMansoTrackerFSMComponent::AddRecHits(
		AliHLTUInt32_t specification,	
		const AliHLTMUONRecHitStruct* recHits,
		AliHLTUInt32_t count
	)
{
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
	
	RecHitBlockInfo info;
	info.fCount = count;
	info.fData = recHits;
	fRecHitBlock[chamber-7].push_back(info);
}


void AliHLTMUONMansoTrackerFSMComponent::RequestClusters(
		AliHLTMUONMansoTrackerFSM* tracker,
		AliHLTFloat32_t /*left*/, AliHLTFloat32_t /*right*/,
		AliHLTFloat32_t /*bottom*/, AliHLTFloat32_t /*top*/,
		AliHLTMUONChamberName chamber, const void* tag
	)
{
	DebugTrace("AliHLTMUONMansoTracker::RequestClusters(chamber = " << chamber << ")");
	void* ctag = const_cast<void*>(tag);
	int chNo = -1;
	std::vector<RecHitBlockInfo>* recHitsBlock = NULL;
	switch (chamber)
	{
	case kChamber7:
		recHitsBlock = &fRecHitBlock[0];
		chNo = 7;
		break;

	case kChamber8:
		recHitsBlock = &fRecHitBlock[1];
		chNo = 8;
		break;

	case kChamber9:
		recHitsBlock = &fRecHitBlock[2];
		chNo = 9;
		break;

	case kChamber10:
		recHitsBlock = &fRecHitBlock[3];
		chNo = 10;
		break;

	default: return;
	}
	
	DebugTrace("Returning requested hits for chamber " << chNo << ":");
	for (AliHLTUInt32_t i = 0; i < recHitsBlock->size(); i++)
	{
		tracker->ReturnClusters(
				ctag,
				(*recHitsBlock)[i].fData,
				(*recHitsBlock)[i].fCount
			);
	}
	DebugTrace("Done returning hits from chamber " << chNo << ".");
	tracker->EndOfClusters(ctag);
}


void AliHLTMUONMansoTrackerFSMComponent::EndOfClusterRequests(
		AliHLTMUONMansoTrackerFSM* tracker
	)
{
	DebugTrace("End of cluster requests.");
}


void AliHLTMUONMansoTrackerFSMComponent::FoundTrack(AliHLTMUONMansoTrackerFSM* tracker)
{
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


void AliHLTMUONMansoTrackerFSMComponent::NoTrackFound(AliHLTMUONMansoTrackerFSM* tracker)
{
	DebugTrace("No track found.");
}
