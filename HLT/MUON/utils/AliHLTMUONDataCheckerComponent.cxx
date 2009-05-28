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

/* $Id: AliHLTMUONDataCheckerComponent.cxx 26179 2008-05-29 22:27:27Z aszostak $ */

///
/// @file   AliHLTMUONDataCheckerComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   27 May 2008
/// @brief  Implementation of the dHLT data integrity checker component.
///
/// This component is used to check the data integrity of dHLT raw internal data
/// blocks. If there are any problems found then an appropriate error message is
/// logged.
///

#include "AliHLTMUONDataCheckerComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTLogging.h"
#include "AliHLTSystem.h"
#include "AliHLTDefinitions.h"
#include "AliRawDataHeader.h"
#include "AliMUONConstants.h"
#include "AliMUONTrackerDDLDecoder.h"
#include "AliMUONTrackerDDLDecoderEventHandler.h"
#include "AliMUONTriggerDDLDecoder.h"
#include "AliMUONTriggerDDLDecoderEventHandler.h"
#include "AliMpDDLStore.h"
#include "AliMpDEStore.h"
#include "AliMpDEManager.h"
#include "AliMpBusPatch.h"
#include "AliMpDetElement.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpPad.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <cassert>


namespace
{
	/**
	 * Routine to check if at least one corresponding DDL has been marked
	 * for a particular chamber.
	 */
	bool ChamberMarkedInDDLList(AliHLTInt32_t chamber, bool ddl[22])
	{
		if (chamber < 0 or chamber > 21) return false;
		switch (chamber)
		{
		case 0:  return ddl[0] or ddl[1];
		case 1:  return ddl[2] or ddl[3];
		case 2:  return ddl[4] or ddl[5];
		case 3:  return ddl[6] or ddl[7];
		case 4:  return ddl[8] or ddl[9] or ddl[10] or ddl[11];
		case 5:  return ddl[8] or ddl[9] or ddl[10] or ddl[11];
		case 6:  return ddl[12] or ddl[13];
		case 7:  return ddl[14] or ddl[15];
		case 8:  return ddl[16] or ddl[17];
		case 9:  return ddl[18] or ddl[19];
		case 10: return ddl[20] or ddl[21];
		case 11: return ddl[20] or ddl[21];
		case 12: return ddl[20] or ddl[21];
		case 13: return ddl[20] or ddl[21];
		default: return false;
		}
	}

} // end of namespace


ClassImp(AliHLTMUONDataCheckerComponent)


AliHLTMUONDataCheckerComponent::AliHLTMUONDataCheckerComponent() :
	AliHLTMUONProcessor(),
	fIgnoreType(false),
	fIgnoreSpec(false),
	fDontForward(false),
	fFilterBadBlocks(false),
	fNoGlobalChecks(false),
	fWarnForUnexpecedBlock(false),
	fReturnError(false)
{
	/// Default constructor.
}


AliHLTMUONDataCheckerComponent::~AliHLTMUONDataCheckerComponent()
{
	/// Default destructor.
}

const char* AliHLTMUONDataCheckerComponent::GetComponentID()
{
	/// Inherited from AliHLTComponent. Returns the component ID.
	
	return AliHLTMUONConstants::DataCheckerComponentId();
}


void AliHLTMUONDataCheckerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	/// Inherited from AliHLTProcessor. Returns the list of expected input data types.
	/// At the moment this list is "any data type" since it is not known before
	/// hand what kind of input blocks we will get.
	
	assert( list.empty() );
	list.push_back( kAliHLTAnyDataType );
}


AliHLTComponentDataType AliHLTMUONDataCheckerComponent::GetOutputDataType()
{
	/// Inherited from AliHLTComponent. Returns the output data type of
	/// "any data type" with MUON origin.
	
	return kAliHLTAnyDataType | kAliHLTDataOriginMUON;
}


void AliHLTMUONDataCheckerComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	/// Inherited from AliHLTComponent.
	/// Returns an estimate of the expected output data size.
	
	// Both of these are zero because we will only ever pass on input data blocks
	// and never generate data in this component.
	constBase = 0;
	inputMultiplier = 0;
}


AliHLTComponent* AliHLTMUONDataCheckerComponent::Spawn()
{
	/// Inherited from AliHLTComponent. Creates a new object instance.
	
	return new AliHLTMUONDataCheckerComponent;
}


int AliHLTMUONDataCheckerComponent::DoInit(int argc, const char** argv)
{
	/// Inherited from AliHLTComponent.
	/// Parses the command line parameters and initialises the component.
	
	HLTInfo("Initialising dHLT data checker component.");
	
	// Inherit the parents functionality.
	int result = AliHLTMUONProcessor::DoInit(argc, argv);
	if (result != 0) return result;

	// Initialise flags with default values.
	fIgnoreType = false;
	fIgnoreSpec = false;
	fDontForward = false;
	fFilterBadBlocks = false;
	fNoGlobalChecks = false;
	fWarnForUnexpecedBlock = false;
	fReturnError = false;

	for (int i = 0; i < argc; i++)
	{
		if (ArgumentAlreadyHandled(i, argv[i])) continue;

		if (strcmp(argv[i], "-ignoretype") == 0)
		{
			fIgnoreType = true;
			HLTInfo("Ignoring data type of data blocks as given by framework.");
			continue;
		}
		if (strcmp(argv[i], "-ignorespec") == 0)
		{
			fIgnoreSpec = true;
			HLTInfo("Ignoring data specification of data blocks as given by framework.");
			continue;
		}
		if (strcmp(argv[i], "-dontforward") == 0)
		{
			fDontForward = true;
			HLTInfo("Not forwarding input data blocks.");
			continue;
		}
		if (strcmp(argv[i], "-filter") == 0)
		{
			fFilterBadBlocks = true;
			HLTInfo("Passing only bad blocks to output.");
			continue;
		}
		if (strcmp(argv[i], "-no_global_check") == 0)
		{
			fNoGlobalChecks = true;
			HLTInfo("Only per block data consistancy checks will be applied,"
				" but no global checks will be made."
			);
			continue;
		}
		if (strcmp(argv[i], "-warn_on_unexpected_block") == 0)
		{
			fWarnForUnexpecedBlock = true;
			continue;
		}
		if (strcmp(argv[i], "-return_error") == 0)
		{
			fReturnError = true;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	}

	return 0;
}


int AliHLTMUONDataCheckerComponent::DoDeinit()
{
	/// Inherited from AliHLTComponent. Performs a cleanup of the component.
	
	HLTInfo("Deinitialising dHLT data checker component.");
	return 0;
}


int AliHLTMUONDataCheckerComponent::DoEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* /*outputPtr*/,
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks
	)
{
	/// Inherited from AliHLTProcessor. Processes the new event data.
	/// Here we go through the list of input data blocks and apply extensive
	/// data integrity checking on the data found.
	
	if (not IsDataEvent()) return 0;
	
	HLTDebug("Processing event %llu with %u input data blocks.",
		evtData.fEventID, evtData.fBlockCnt
	);
	
	// Allocate an array of flags indicating if the data block is OK or not,
	// also arrays to store specific.
	bool dataProblems = false;
	bool* blockOk = NULL;
	typedef const AliHLTComponentBlockData* PAliHLTComponentBlockData;
	PAliHLTComponentBlockData* trigRecBlocks = NULL;
	PAliHLTComponentBlockData* trigRecDebugBlocks = NULL;
	PAliHLTComponentBlockData* hitBlocks = NULL;
	PAliHLTComponentBlockData* clusterBlocks = NULL;
	PAliHLTComponentBlockData* channelBlocks = NULL;
	PAliHLTComponentBlockData* mansoTrackBlocks = NULL;
	PAliHLTComponentBlockData* mansoCandidateBlocks = NULL;
	PAliHLTComponentBlockData* singleDecisionBlocks = NULL;
	PAliHLTComponentBlockData* pairDecisionBlocks = NULL;
	AliHLTUInt32_t trigRecBlocksCount = 0;
	AliHLTUInt32_t trigRecDebugBlocksCount = 0;
	AliHLTUInt32_t hitBlocksCount = 0;
	AliHLTUInt32_t clusterBlocksCount = 0;
	AliHLTUInt32_t channelBlocksCount = 0;
	AliHLTUInt32_t mansoTrackBlocksCount = 0;
	AliHLTUInt32_t mansoCandidateBlocksCount = 0;
	AliHLTUInt32_t singleDecisionBlocksCount = 0;
	AliHLTUInt32_t pairDecisionBlocksCount = 0;
	try
	{
		blockOk = new bool[evtData.fBlockCnt];
		trigRecBlocks = new PAliHLTComponentBlockData[evtData.fBlockCnt];
		trigRecDebugBlocks = new PAliHLTComponentBlockData[evtData.fBlockCnt];
		hitBlocks = new PAliHLTComponentBlockData[evtData.fBlockCnt];
		clusterBlocks = new PAliHLTComponentBlockData[evtData.fBlockCnt];
		channelBlocks = new PAliHLTComponentBlockData[evtData.fBlockCnt];
		mansoTrackBlocks = new PAliHLTComponentBlockData[evtData.fBlockCnt];
		mansoCandidateBlocks = new PAliHLTComponentBlockData[evtData.fBlockCnt];
		singleDecisionBlocks = new PAliHLTComponentBlockData[evtData.fBlockCnt];
		pairDecisionBlocks = new PAliHLTComponentBlockData[evtData.fBlockCnt];
	}
	catch (const std::bad_alloc&)
	{
		HLTError("Could not allocate more memory for internal arrays.");
		// Make sure to clean up if partially allocated memory.
		if (blockOk != NULL) delete [] blockOk;
		if (trigRecBlocks != NULL) delete [] trigRecBlocks;
		if (trigRecDebugBlocks != NULL) delete [] trigRecDebugBlocks;
		if (hitBlocks != NULL) delete [] hitBlocks;
		if (clusterBlocks != NULL) delete [] clusterBlocks;
		if (channelBlocks != NULL) delete [] channelBlocks;
		if (mansoTrackBlocks != NULL) delete [] mansoTrackBlocks;
		if (mansoCandidateBlocks != NULL) delete [] mansoCandidateBlocks;
		if (singleDecisionBlocks != NULL) delete [] singleDecisionBlocks;
		if (pairDecisionBlocks != NULL) delete [] pairDecisionBlocks;
		return -ENOMEM;
	}
	
	AliHLTComponentDataType anyPrivateType = AliHLTComponentDataTypeInitializer(
			kAliHLTAnyDataType, kAliHLTDataOriginPrivate
		);
	
	try
	{
		// Clear all the flags indicating if the blocks are ok.
		for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
		{
			blockOk[n] = false;
		}
	
		for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
		{
			HLTDebug("Handling block: %u, with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
				n, DataType2Text(blocks[n].fDataType).c_str(), blocks[n].fPtr, blocks[n].fSize
			);
			
			AliHLTMUONDataBlockType blockType = kUnknownDataBlock;
			
			if (fIgnoreType)
			{
				// Decode the block type if we must ignore the block type
				// as given by the HLT framework.
				if (blocks[n].fSize >= sizeof(AliHLTMUONDataBlockHeader))
				{
					const AliHLTMUONDataBlockHeader* header =
						reinterpret_cast<const AliHLTMUONDataBlockHeader*>(blocks[n].fPtr);
					blockType = AliHLTMUONDataBlockType(header->fType);
				}
			}
			else
			{
				if (blocks[n].fDataType == anyPrivateType)
				{
					// Completely ignore any private HLT internal block types.
					blockOk[n] = true;
					continue;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::DDLRawDataType())
				{
					blockOk[n] = CheckRawDataBlock(blocks[n], n);
					continue;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::TriggerRecordsBlockDataType())
				{
					blockType = kTriggerRecordsDataBlock;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::TrigRecsDebugBlockDataType())
				{
					blockType = kTrigRecsDebugDataBlock;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::RecHitsBlockDataType())
				{
					blockType = kRecHitsDataBlock;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::ClusterBlockDataType())
				{
					blockType = kClustersDataBlock;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::ChannelBlockDataType())
				{
					blockType = kChannelsDataBlock;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::MansoTracksBlockDataType())
				{
					blockType = kMansoTracksDataBlock;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::MansoCandidatesBlockDataType())
				{
					blockType = kMansoCandidatesDataBlock;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::SinglesDecisionBlockDataType())
				{
					blockType = kSinglesDecisionDataBlock;
				}
				else if (blocks[n].fDataType == AliHLTMUONConstants::PairsDecisionBlockDataType())
				{
					blockType = kPairsDecisionDataBlock;
				}
				else
				{
					// Log a message indicating that we got a data block that we
					// do not know how to handle.
					if (fWarnForUnexpecedBlock)
						HLTWarning("Received a data block of a type we cannot"
							" handle: '%s', spec: 0x%8.8X",
							DataType2Text(blocks[n].fDataType).c_str(),
							blocks[n].fSpecification
						);
#ifdef __DEBUG
					else
						HLTDebug("Received a data block of a type we cannot"
							" handle: '%s', spec: 0x%8.8X",
							DataType2Text(blocks[n].fDataType).c_str(),
							blocks[n].fSpecification
						);
#endif
				}
			}
			
			switch (blockType)
			{
			case kTriggerRecordsDataBlock:
				blockOk[n] = CheckTriggerRecordsBlock(blocks[n], n);
				trigRecBlocks[trigRecBlocksCount++] = &blocks[n];
				break;
			case kTrigRecsDebugDataBlock:
				blockOk[n] = CheckTrigRecsDebugBlock(blocks[n], n);
				trigRecDebugBlocks[trigRecDebugBlocksCount++] = &blocks[n];
				break;
			case kRecHitsDataBlock:
				blockOk[n] = CheckRecHitsBlock(blocks[n], n);
				hitBlocks[hitBlocksCount++] = &blocks[n];
				break;
			case kClustersDataBlock:
				blockOk[n] = CheckClustersBlock(blocks[n], n);
				clusterBlocks[clusterBlocksCount++] = &blocks[n];
				break;
			case kChannelsDataBlock:
				blockOk[n] = CheckChannelsBlock(blocks[n], n);
				channelBlocks[channelBlocksCount++] = &blocks[n];
				break;
			case kMansoTracksDataBlock:
				blockOk[n] = CheckMansoTracksBlock(blocks[n], n);
				mansoTrackBlocks[mansoTrackBlocksCount++] = &blocks[n];
				break;
			case kMansoCandidatesDataBlock:
				blockOk[n] = CheckMansoCandidatesBlock(blocks[n], n);
				mansoCandidateBlocks[mansoCandidateBlocksCount++] = &blocks[n];
				break;
			case kSinglesDecisionDataBlock:
				blockOk[n] = CheckSinglesDecisionBlock(blocks[n], n);
				singleDecisionBlocks[singleDecisionBlocksCount++] = &blocks[n];
				break;
			case kPairsDecisionDataBlock:
				blockOk[n] = CheckPairsDecisionBlock(blocks[n], n);
				pairDecisionBlocks[pairDecisionBlocksCount++] = &blocks[n];
				break;
			default:
				HLTDebug("Received a data block for which we could not decode the data type."
					" fDataType = '%s', fSpecification = 0x%8.8X, fSize = %u bytes.",
					DataType2Text(blocks[n].fDataType).c_str(),
					blocks[n].fSpecification,
					blocks[n].fSize
				);
				break;
			}
		}
		
		// Apply the global data consistancy checks if not suppressed by the user.
		if (not fNoGlobalChecks)
		{
			MakeGlobalChecks(
				blocks, blockOk, evtData.fBlockCnt,
				trigRecBlocks, trigRecBlocksCount,
				trigRecDebugBlocks, trigRecDebugBlocksCount,
				hitBlocks, hitBlocksCount,
				clusterBlocks, clusterBlocksCount,
				channelBlocks, channelBlocksCount,
				mansoTrackBlocks, mansoTrackBlocksCount,
				mansoCandidateBlocks, mansoCandidateBlocksCount,
				singleDecisionBlocks, singleDecisionBlocksCount,
				pairDecisionBlocks, pairDecisionBlocksCount
			);
		}
		
		// Forward the input data blocks if we have not been asked to drop them.
		// Also remember to filter for bad blocks if so specified.
		if (not fDontForward)
		{
			for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
			{
				if (fFilterBadBlocks and blockOk[n]) continue;
				outputBlocks.push_back(blocks[n]);
			}
		}
		
		// Set dataProblems flag is there was at least one block with problems.
		for (AliHLTUInt32_t n = 0; n < evtData.fBlockCnt; n++)
		{
			if (not blockOk[n]) dataProblems = true;
		}
	}
	finally
	(
		// make sure to cleanup memory
		delete [] blockOk;
		delete [] trigRecBlocks;
		delete [] trigRecDebugBlocks;
		delete [] hitBlocks;
		delete [] clusterBlocks;
		delete [] channelBlocks;
		delete [] mansoTrackBlocks;
		delete [] mansoCandidateBlocks;
		delete [] singleDecisionBlocks;
		delete [] pairDecisionBlocks;
	)
	
	// Finally we set the total size of output memory we consumed, which is
	// zero since we just copied the input descriptors to output if anything.
	size = 0;

	if (dataProblems and DumpDataOnError()) DumpEvent(evtData, trigData);
	
	if (fReturnError)
	{
		// If we were requested to return errors if there were integrity
		// problems then check if any data blocks had problems and return
		// an error code.
		if (dataProblems) return -EFAULT;
	}
	return 0;
}


bool AliHLTMUONDataCheckerComponent::IsSpecificationValid(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber,
		const char* name
	) const
{
	/// Checks if the specification bits are valid.
	/// \param  block The block whose specification should be checked.
	/// \param  blockNumber The block index number being checked.
	/// \param  name The name of the type of block being checked.
	/// \returns true if the specification is valid and false otherwise.

	if (AliHLTMUONUtils::IsSpecValid(block.fSpecification))
		return true;
	
	HLTError("Problem found with data block %d, fDataType = '%s',"
		 " fPtr = %p and fSize = %u bytes."
		 " Assuming this is a %s data block."
		 " Problem: The specification does not contain a valid pattern,"
		 " received 0x%8.8X for the specification.",
		blockNumber,
		DataType2Text(block.fDataType).c_str(),
		block.fPtr,
		block.fSize,
		name,
		block.fSpecification
	);
	return false;
}


bool AliHLTMUONDataCheckerComponent::IsFromTrackerOnly(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber,
		const char* name
	) const
{
	/// Checks if the specification bits are valid and indicate the block
	/// contains data or information only from the tracker DDLs.
	/// \param  block The block whose specification should be checked.
	/// \param  blockNumber The block index number being checked.
	/// \param  name The name of the type of block being checked.
	/// \returns true if the specification indicates data is only from tracker.
	
	bool result = IsSpecificationValid(block, blockNumber, name);
	
	if (AliHLTMUONUtils::ContainsDataFromTracker(block.fSpecification) and
	    not AliHLTMUONUtils::ContainsDataFromTrigger(block.fSpecification)
	   )
	{
		return result;
	}
	
	HLTError("Problem found with data block %d, fDataType = '%s',"
		 " fPtr = %p and fSize = %u bytes."
		 " Assuming this is a %s data block."
		 " Problem: The data block does not contain data only from the"
		 " tracker DDLs as expected."
		 " Received 0x%8.8X for the specification.",
		blockNumber,
		DataType2Text(block.fDataType).c_str(),
		block.fPtr,
		block.fSize,
		name,
		block.fSpecification
	);
	return false;
}


bool AliHLTMUONDataCheckerComponent::IsFromTriggerOnly(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber,
		const char* name
	) const
{
	/// Checks if the specification bits are valid and indicate the block
	/// contains data or information only from the trigger DDLs.
	/// \param  block The block whose specification should be checked.
	/// \param  blockNumber The block index number being checked.
	/// \param  name The name of the type of block being checked.
	/// \returns true if the specification indicates data is only from trigger.
	
	bool result = IsSpecificationValid(block, blockNumber, name);
	
	if (AliHLTMUONUtils::ContainsDataFromTrigger(block.fSpecification) and
	    not AliHLTMUONUtils::ContainsDataFromTracker(block.fSpecification)
	   )
	{
		return result;
	}
	
	HLTError("Problem found with data block %d, fDataType = '%s',"
		 " fPtr = %p and fSize = %u bytes."
		 " Assuming this is a %s data block."
		 " Problem: The data block does not contain data only from the"
		 " trigger DDLs as expected."
		 " Received 0x%8.8X for the specification.",
		blockNumber,
		DataType2Text(block.fDataType).c_str(),
		block.fPtr,
		block.fSize,
		name,
		block.fSpecification
	);
	return false;
}


bool AliHLTMUONDataCheckerComponent::IsMomentumVectorOk(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber,
		const char* name,
		AliHLTUInt32_t entryNumber,
		AliHLTFloat32_t px,
		AliHLTFloat32_t py,
		AliHLTFloat32_t pz
	) const
{
	/// Checks if the momentum vector is reasonable.
	/// \param  block The block from which the momentum vector data comes from.
	/// \param  blockNumber The block index number.
	/// \param  name The name of the type of block.
	/// \param  entryNumber The entry index number of the structure holding
	///      the momentum vector data.
	/// \param  px The X coordinate of the momentum vector (GeV/c).
	/// \param  py The Y coordinate of the momentum vector (GeV/c).
	/// \param  pz The Z coordinate of the momentum vector (GeV/c).
	/// \returns true if the momentum vector is valid and false otherwise.
	
	// If the momentum vector is nil then ignore it.
	if (px == 0 and py == 0 and pz == 0) return true;
	
	bool result = true;
	
	// If the momentum vector is sane then we should not have a particle with
	// more energy than 14 TeV and momentum should be in the negative direction.
	double momentum = sqrt(px*px + py*py + pz*pz);
	if (momentum > 14e3)
	{
		// Just warn since this is a data sanity problem rather
		// than a data integrity problem.
		HLTWarning("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The momentum vector"
			" p = {%f, %f, %f}, |p| = %f looks too big.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			px, py, pz,
			momentum
		);
		result = false;
	}
	
	if (pz > 0.)
	{
		// Just warn since this is a data sanity problem rather
		// than a data integrity problem.
		HLTWarning("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The momentum vector"
			" p = {%f, %f, %f} points away from the dimuon"
			" spectrometer (p_z > 0).",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			px, py, pz
		);
		result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::AreMomentumCalcParamsOk(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber,
		const char* name,
		AliHLTUInt32_t entryNumber,
		AliHLTFloat32_t zmiddle,
		AliHLTFloat32_t bl
	) const
{
	/// Checks if the parameters for the momentum calculation are reasonable.
	/// \param  block The block from which the parameter data comes from.
	/// \param  blockNumber The block index number.
	/// \param  name The name of the type of block.
	/// \param  entryNumber The entry index number of the structure holding
	///      the parameter data data.
	/// \param  zmiddle The z-coordinate of the middle of the magnetic field (cm).
	/// \param  bl The integrated magnetic field (T.m).
	/// \returns true if the parameters are valid and false otherwise.
	
	bool result = true;
	
	// Check that the value of the fZmiddle value is somewhere
	// within the tracking / dipole magnetic field area.
	if (zmiddle < AliMUONConstants::AbsZEnd() or
		zmiddle < AliMUONConstants::MuonFilterZBeg()
		)
	{
		// Just warn since this is a data sanity problem rather
		// than a data integrity problem.
		HLTWarning("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The Z coordinate %f cm"
			" used as the middle of the magnetic field in the momentum"
			" calculation is outside the dimuon spectrometers dipole"
			" magnetic field volume.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			zmiddle
		);
		result = false;
	}
	
	// Also check that the value of the 'bl' value is within a
	// reasonable range: |bl| < Lmax * Bmax, where
	// Lmax = max length from vertex to end of spectrometer, and
	// Bmax = max magnetic field of dipole, taken as 1 tesla.
	// Approximating Lmax * Bmax as 20 T.m
	if (fabs(bl) > 20.)
	{
		// Just warn since this is a data sanity problem rather
		// than a data integrity problem.
		HLTWarning("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The integrated magnetic"
			" field value %f T.m used in the momentum calculation"
			" has an unreasonably large absolute value.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			bl
		);
		result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::IsHitCoordinateOk(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber,
		const char* name,
		AliHLTUInt32_t entryNumber,
		const AliHLTMUONRecHitStruct& hit,
		AliHLTInt32_t minChamber,
		AliHLTInt32_t maxChamber,
		AliHLTInt32_t expectedChamber,
		bool ddl[22]
	) const
{
	/// Checks if the hit coordinate is compatible with a the location of a
	/// dimuon spectrometer chamber. Also, if expectedChamber is not -1, then
	/// the hit coordinate is checked if to comes from that chamber.
	/// We also check if the fFlags containing the chamber number and detector
	/// element ID are correct.
	/// \param  block The block from which the hit data comes from.
	/// \param  blockNumber The block index number.
	/// \param  name The name of the type of block.
	/// \param  entryNumber The entry index number of the hit.
	/// \param  hit The hit data being checked.
	/// \param  minChamber The minimum valid chamber number to check for.
	/// \param  maxChamber The maximum valid chamber number to check for.
	/// \param  expectedChamber If not -1 then this is the chamber number to
	///      check against.
	/// \param  ddl  The array decoded by AliHLTMUONUtils::UnpackSpecBits.
	/// \returns true if the hit is valid and false otherwise.
	
	assert( 0 <= minChamber and minChamber < 14 );
	assert( 0 <= maxChamber and maxChamber < 14 );
	
	bool result = true;
	
	AliHLTUInt8_t chNum = 0xFF;
	AliHLTUInt16_t detElemId = 0xFFFF;
	AliHLTMUONUtils::UnpackRecHitFlags(hit.fFlags, chNum, detElemId);
	
	Int_t chamber = AliMUONConstants::ChamberNumber(hit.fZ, false); // false = do not warn.
	if (chamber < minChamber or maxChamber < chamber)
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The hit {x = %f, y = %f,"
			" z = %f} cm has a z-coordinate that does not correspond"
			" to the nominal position of any chambers in the range"
			" [%d..%d].",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			hit.fX, hit.fY, hit.fZ,
			minChamber+1,
			maxChamber+1
		);
		return false;
	}
	
	if (chNum != chamber)
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The hit {x = %f, y = %f,"
			" z = %f} cm has a chamber number %d that does not correspond"
			" to the expected chamber %d given by the z-coordinate.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			hit.fX, hit.fY, hit.fZ,
			chNum+1,
			chamber+1
		);
		result = false;
		if (minChamber <= Int_t(chNum) and Int_t(chNum) <= maxChamber)
		{
			// Rather use the explicit value in the data if it
			// is in range.
			chamber = chNum;
		}
	}
	
	if (expectedChamber != -1 and chamber != expectedChamber)
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The hit {x = %f, y = %f,"
			" z = %f} cm has a position that corresponds to chamber %d,"
			" but expected it to be on chamber %d.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			hit.fX, hit.fY, hit.fZ,
			chamber+1,
			expectedChamber+1
		);
		result = false;
	}
	
	AliHLTFloat32_t rmin = AliMUONConstants::Rmin(chamber / 2);
	AliHLTFloat32_t rmax = AliMUONConstants::Rmax(chamber / 2);
	AliHLTFloat32_t radius = sqrt(hit.fX*hit.fX + hit.fY*hit.fY);
	if (radius < rmin or rmax < radius)
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The hit {x = %f, y = %f,"
			" z = %f} cm has a position in the X-Y plane that does not"
			" correspond to the nominal position of chamber %d.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			hit.fX, hit.fY, hit.fZ,
			chamber+1
		);
		result = false;
	}
	
	if (not fIgnoreSpec and not ChamberMarkedInDDLList(chamber, ddl))
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The hit {x = %f, y = %f,"
			" z = %f} cm has a position that corresponds to chamber %d"
			" but the data block specification 0x%8.8X does have a"
			" corresponding DDL bit set.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			hit.fX, hit.fY, hit.fZ,
			chamber+1,
			block.fSpecification
		);
		result = false;
	}
	
	// Check that the detector element ID is valid and it corresponds to
	// the chamber number.
	if (FetchMappingStores() == 0)  // are stores loaded?
	{
		Bool_t warn = kFALSE;
		AliMpDEStore* store = AliMpDEStore::Instance(warn);
		AliMpDetElement* de = store->GetDetElement(Int_t(detElemId), warn);
		if (de == NULL)
		{
			HLTError("Problem found with data block %d, fDataType = '%s',"
				" fPtr = %p and fSize = %u bytes."
				" Assuming this is a %s data block."
				" Problem with entry %d in block: The hit {x = %f, y = %f,"
				" z = %f} cm has a detector element ID %d,"
				" which is not valid.",
				blockNumber,
				DataType2Text(block.fDataType).c_str(),
				block.fPtr,
				block.fSize,
				name,
				entryNumber,
				hit.fX, hit.fY, hit.fZ,
				detElemId
			);
			result = false;
		}
			
		// Check that the chamber number from the detector element number
		// has the expected value.
		Int_t ch = AliMpDEManager::GetChamberId(Int_t(detElemId), warn);
		if (ch != chamber)
		{
			HLTError("Problem found with data block %d, fDataType = '%s',"
				" fPtr = %p and fSize = %u bytes."
				" Assuming this is a %s data block."
				" Problem with entry %d in block: The hit {x = %f, y = %f,"
				" z = %f} cm has a detector element ID %d,"
				" which does not correspond to the chamber %d.",
				blockNumber,
				DataType2Text(block.fDataType).c_str(),
				block.fPtr,
				block.fSize,
				name,
				entryNumber,
				hit.fX, hit.fY, hit.fZ,
				detElemId,
				chamber+1
			);
			result = false;
		}
	}
	else
	{
		HLTWarning("Cannot check a hit's detector element ID information"
			" without being able to load the mapping from CDB."
		);
		result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::IsMansoTrackOk(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber,
		const char* name,
		AliHLTUInt32_t entryNumber,
		const AliHLTMUONMansoTrackStruct& track,
		bool ddl[22]
	) const
{
	/// Checks if the Manso track structure is Ok.
	/// \param  block The block from which the track data comes from.
	/// \param  blockNumber The block index number.
	/// \param  name The name of the type of block.
	/// \param  entryNumber The entry index number of the structure in the
	///      block being checked.
	/// \param  track The Manso track data being checked.
	/// \param  ddl  The array decoded by AliHLTMUONUtils::UnpackSpecBits.
	/// \returns true if the hit is valid and false otherwise.
	
	bool result = true;
	
	// Chi^2 should not be greater than the worst fit possible, estimated
	// as the diameter of largest chamber times the number of points
	// findable in a track. Max points is 10 tracker chambers times
	// 2 cathodes + 4 trigger chambers.
	if (track.fChi2 > AliMUONConstants::Dmax(6)*AliMUONConstants::Dmax(6)*(10*2+4))
	{
		// Just a warning since this is not technically an
		// integrity problem.
		HLTWarning("Problem found with data block %d, fDataType = '%s',"
			" fPtr = %p and fSize = %u bytes."
			" Assuming this is a %s data block."
			" Problem with entry %d in block: The Manso track has"
			" the chi squared value of %f that unreasonably big.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			name,
			entryNumber,
			track.fChi2
		);
		result = false;
	}
	
	// Check if the momentum vector is reasonable.
	bool momOk = IsMomentumVectorOk(
			block, blockNumber, name, entryNumber,
			track.fPx, track.fPy, track.fPz
		);
	if (not momOk) result = false;
	
	AliHLTMUONParticleSign sign;
	bool hitset[4];
	AliHLTMUONUtils::UnpackMansoTrackFlags(track.fFlags, sign, hitset);
	
	// Min and max allowed chamber numbers for hits:
	Int_t minCh = 0;
	Int_t maxCh = AliMUONConstants::NTrackingCh() - 1;
	
	// Check that this hit coordinates are OK.
	for (AliHLTUInt32_t i = 0; i < 4; i++)
	{
		if (not hitset[i]) continue; // ignore hits that are not initialised.
		bool hitOk = IsHitCoordinateOk(
				block, blockNumber, name, entryNumber, track.fHit[i],
				minCh, maxCh, i+6, ddl
			);
		if (not hitOk) result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckDetElemIds(
		const AliHLTComponentBlockData& infoBlock,
		AliHLTUInt32_t infoBlockNumber,
		AliHLTUInt32_t infoEntryNumber,
		const AliHLTMUONTrigRecInfoStruct& info,
		const AliHLTComponentBlockData& trBlock,
		AliHLTUInt32_t trBlockNumber,
		AliHLTUInt32_t trEntryNumber,
		const AliHLTMUONTriggerRecordStruct& tr
	) const
{
	/// Checks if the detector element IDs are the same in the debug
	/// information structure and the trigger record structure.
	/// \param  infoBlock The debug information block from which the 'info'
	///      data comes from.
	/// \param  infoBlockNumber The debug information block index number.
	/// \param  infoEntryNumber The entry index number of the 'info'
	///      structure in the debug information data block.
	/// \param  info  The debug information structure being checked.
	/// \param  trBlock The trigger record block from which the 'tr' data
	///      comes from.
	/// \param  trBlockNumber The trigger record block index number.
	/// \param  trEntryNumber The entry index number of the 'tr' structure
	///      in the trigger record data block.
	/// \param  tr  The trigger record structure being checked.
	/// \returns true if the detector element IDs are the same and false
	///      otherwise.
	
	bool result = true;
	
	for (int i = 0; i < 4; i++)
	{
		AliHLTUInt8_t chamber = 0xFF;
		AliHLTUInt16_t detElemId = 0xFFFF;
		AliHLTMUONUtils::UnpackRecHitFlags(tr.fHit[i].fFlags, chamber, detElemId);
		if (info.fDetElemId[i] == detElemId) continue;
		
		HLTError("Problem found with trigger record debug information %d"
			" in data block %d (fDataType = '%s', fPtr = %p, fSize"
			" = %u bytes) and trigger record %d in data block %d"
			" (fDataType = '%s', fPtr = %p, fSize = %u bytes):"
			" The detection element ID %d for chamber %d in the debug"
			" information, is not the same as %d"
			" found in the trigger record.",
			infoEntryNumber,
			infoBlockNumber,
			DataType2Text(infoBlock.fDataType).c_str(),
			infoBlock.fPtr,
			infoBlock.fSize,
			trEntryNumber,
			trBlockNumber,
			DataType2Text(trBlock.fDataType).c_str(),
			trBlock.fPtr,
			trBlock.fSize,
			info.fDetElemId[i],
			i+11,
			detElemId
		);
		result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckDetElemIds(
		const AliHLTComponentBlockData& clusterBlock,
		AliHLTUInt32_t clusterBlockNumber,
		AliHLTUInt32_t clusterEntryNumber,
		const AliHLTMUONClusterStruct& cluster,
		const AliHLTComponentBlockData& hitBlock,
		AliHLTUInt32_t hitBlockNumber,
		AliHLTUInt32_t hitEntryNumber,
		const AliHLTMUONRecHitStruct& hit
	) const
{
	/// Checks if the detector element IDs are the same in the cluster
	/// structure and the reconstructed hit structure.
	/// \param  clusterBlock The cluster block from which the 'cluster' data
	///      comes from.
	/// \param  clusterBlockNumber The cluster block index number.
	/// \param  clusterEntryNumber The entry index number of the 'cluster'
	///      structure in the cluster data block.
	/// \param  cluster  The cluster structure being checked.
	/// \param  hitBlock The reconstructed hit block from which the 'hit'
	///      data comes from.
	/// \param  hitBlockNumber The reconstructed hit block index number.
	/// \param  hitEntryNumber The entry index number of the 'hit' structure
	///      in the reconstructed hit data block.
	/// \param  hit  The trigger record structure being checked.
	/// \returns true if the detector element IDs are the same and false
	///      otherwise.
	
	bool result = true;
	
	AliHLTUInt8_t chamber = 0xFF;
	AliHLTUInt16_t detElemId = 0xFFFF;
	AliHLTMUONUtils::UnpackRecHitFlags(hit.fFlags, chamber, detElemId);
	if (cluster.fDetElemId != detElemId)
	{
		HLTError("Problem found with cluster %d in data block %d"
			" (fDataType = '%s', fPtr = %p, fSize = %u bytes)"
			" and reconstructed hit %d in data block %d"
			" (fDataType = '%s', fPtr = %p, fSize = %u bytes):"
			" The detection element ID %d in the cluster, is not"
			" the same as %d found in the reconstructed hit.",
			clusterEntryNumber,
			clusterBlockNumber,
			DataType2Text(clusterBlock.fDataType).c_str(),
			clusterBlock.fPtr,
			clusterBlock.fSize,
			hitEntryNumber,
			hitBlockNumber,
			DataType2Text(hitBlock.fDataType).c_str(),
			hitBlock.fPtr,
			hitBlock.fSize,
			cluster.fDetElemId,
			detElemId
		);
		result = false;
	}
	
	return result;
}


namespace
{
	/**
	 * Class for logging errors found in raw DDL data.
	 */
	class AliHLTMUONDecoderHandler : public AliHLTLogging
	{
	public:
	
		/// Default constructor
		AliHLTMUONDecoderHandler() :
			AliHLTLogging(),
			fBufferStart(NULL),
			fDescriptor(NULL),
			fBlockNumber(0)
		{
		}
		
		/// Default destructor.
		virtual ~AliHLTMUONDecoderHandler() {}
		
		/// Sets the DDL raw data block descriptor.
		void SetDescriptor(const AliHLTComponentBlockData* b) { fDescriptor = b; }
		
		/// Sets the block number of the raw data block descriptor.
		void SetBlockNumber(AliHLTUInt32_t n) { fBlockNumber = n; }
		
		/// Logs an error message describing the problem with the DDL raw data.
		template <typename ErrorCode, class DecoderHandler>
		void LogError(ErrorCode code, const void* location, DecoderHandler& handler);
	
	protected:
		// Do not allow copying of this class.
		/// Not implemented
		AliHLTMUONDecoderHandler(const AliHLTMUONDecoderHandler& rhs); // copy constructor
		/// Not implemented
		AliHLTMUONDecoderHandler& operator = (const AliHLTMUONDecoderHandler& rhs); // assignment operator
		
		const void* fBufferStart; ///< Pointer to the start of the current DDL payload buffer.
		const AliHLTComponentBlockData* fDescriptor; ///< Descriptor for the DDL raw data block corresponding to the buffer.
		AliHLTUInt32_t fBlockNumber;  ///< The number / index of the block descriptor.
	};
	
	
	template <typename ErrorCode, class DecoderHandler>
	void AliHLTMUONDecoderHandler::LogError(ErrorCode code, const void* location, DecoderHandler& handler)
	{
		/// Logs a HLT error message describing the problem with the raw DDL data.
		/// \param code  The error code describing the problem.
		/// \param location  A pointer to the location in the raw data buffer
		///      where the problem was found.
		/// \param handler  The decoder handler object.
		
		long bytepos = long(location) - long(fBufferStart) + sizeof(AliRawDataHeader);
		
		// create data type string.
		char dataType[kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2];
		memset( dataType, 0, kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2 );
		strncat( dataType, fDescriptor->fDataType.fOrigin, kAliHLTComponentDataTypefOriginSize );
		strcat( dataType, ":" );
		strncat( dataType, fDescriptor->fDataType.fID, kAliHLTComponentDataTypefIDsize );
		
		HLTError("Problem found with data block %d, fDataType = '%s',"
			 " fPtr = %p and fSize = %u bytes."
			 " Assuming this is a DDL raw data block."
			 " Problem: %s (Error code: %d, at byte %d)",
			fBlockNumber,
			&dataType[0],
			fDescriptor->fPtr,
			fDescriptor->fSize,
			handler.ErrorCodeToMessage(code),
			code,
			bytepos
		);
	};
	
	
	/**
	 * Class for logging decoding errors when checking tracker raw DDL data.
	 * Used in the AliHLTMUONDataCheckerComponent::CheckRawDataBlock method.
	 */
	class AliHLTMUONTrackerDecoderHandler :
		public AliMUONTrackerDDLDecoderEventHandler, public AliHLTMUONDecoderHandler
	{
	public:
		AliHLTMUONTrackerDecoderHandler() :
			AliMUONTrackerDDLDecoderEventHandler(),
			AliHLTMUONDecoderHandler(),
			fMaxDigits(0),
			fDigitCount(0),
			fDigits(NULL),
			fCurrentBusPatch(0),
			fDataProblems(false)
		{}
		
		virtual ~AliHLTMUONTrackerDecoderHandler()
		{
			if (fDigits != NULL) delete [] fDigits;
		}
		
		/// Structure to store raw data words found in the raw data.
		struct AliDigit
		{
			UInt_t fBusPatchId;  ///< Bus patch ID for the data word.
			UInt_t fDataWord;   ///< Raw data word found in the DDL payload.
		};
		
		/// Returns the number of digits found.
		UInt_t DigitCount() const { return fDigitCount; }
		
		/// Returns the array of digits found.
		const AliDigit* Digits() const { return fDigits; }
		
		/// Returns true if there were problems with the data.
		bool DataProblems() const { return fDataProblems; }
		
		// Methods inherited from AliMUONTrackerDDLDecoderEventHandler:
		
		/// Called for each new buffer.
		void OnNewBuffer(const void* buffer, UInt_t bufferSize);
		
		/// Called for each new DSP header.
		void OnNewDSP(const AliMUONDSPHeaderStruct* header, const void* /*data*/);
		
		/// Called for each new bus patch. Just marks the current bus patch ID.
		void OnNewBusPatch(const AliMUONBusPatchHeaderStruct* header, const void* /*data*/)
		{
			fCurrentBusPatch = header->fBusPatchId;
		}
		
		/// Called for each new data word found.
		void OnData(UInt_t data, bool /*parityError*/);
		
		/// Logs an error message if there was a decoding problem with the DDL payload.
		void OnError(ErrorCode code, const void* location)
		{
			fDataProblems = true;
			LogError(code, location, *this);
		}
	
	private:
	
		// Do not allow copying of this object.
		/// Not implemented.
		AliHLTMUONTrackerDecoderHandler(const AliHLTMUONTrackerDecoderHandler& obj);
		/// Not implemented.
		AliHLTMUONTrackerDecoderHandler& operator = (const AliHLTMUONTrackerDecoderHandler& obj);
		
		UInt_t fMaxDigits;  ///< Maximum number of digits that can be stored in fDigits.
		UInt_t fDigitCount;  ///< The number of digits currently stored in fDigits.
		AliDigit* fDigits;  ///< The array of digits found in the DDL data.
		UInt_t fCurrentBusPatch;  ///< The current bus patch ID being processed.
		bool fDataProblems;  ///< flag indicating there were problems with the data.
	};
	
	
	void AliHLTMUONTrackerDecoderHandler::OnNewBuffer(const void* buffer, UInt_t bufferSize)
	{
		/// Called for a new buffer. It will reset internal counters and
		/// resize the digits array if necessary.
		
		fDataProblems = false;
		fDigitCount = 0;
		fBufferStart = buffer;
		
		// Resize the fDigits array to be able to store
		// all the digits in the data buffer.
		UInt_t maxSize = bufferSize / sizeof(UInt_t) + 1;
		if (maxSize > fMaxDigits)
		{
			if (fDigits != NULL)
			{
				delete [] fDigits;
				fDigits = NULL;
				fMaxDigits = 0;
			}
			try
			{
				fDigits = new AliDigit[maxSize];
				fMaxDigits = maxSize;
			}
			catch (const std::bad_alloc&)
			{
				HLTError("Could not allocate enough buffer space for internal arrays.");
				return;
			}
		}
	}
	
	
	void AliHLTMUONTrackerDecoderHandler::OnNewDSP(
			const AliMUONDSPHeaderStruct* header, const void* /*data*/
		)
	{
		if (header->fPaddingWord != 0 and header->fPaddingWord != 1)
		{
			// create data type string.
			char dataType[kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2];
			memset( dataType, 0, kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2 );
			strncat( dataType, fDescriptor->fDataType.fOrigin, kAliHLTComponentDataTypefOriginSize );
			strcat( dataType, ":" );
			strncat( dataType, fDescriptor->fDataType.fID, kAliHLTComponentDataTypefIDsize );
				
			HLTError("Problem found with data block %d, fDataType = '%s',"
				" fPtr = %p and fSize = %u bytes."
				" Assuming this is a tracker DDL raw data block."
				" Problem: Found padding word marker 0x%8.8X in DSP"
				" header with DSP ID %d which has an invalid value.",
				fBlockNumber,
				&dataType[0],
				fDescriptor->fPtr,
				fDescriptor->fSize,
				header->fPaddingWord,
				header->fDSPId
			);
			fDataProblems = true;
			return;
		}
	}
	
	
	void AliHLTMUONTrackerDecoderHandler::OnData(UInt_t data, bool /*parityError*/)
	{
		/// Called for each new data word found. This method will add
		/// these to the list of digits and check if they are not duplicated.
		
		assert( fDigits != NULL );
		
		if ((data & 0x60000000) != 0)
		{
			// create data type string.
			char dataType[kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2];
			memset( dataType, 0, kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2 );
			strncat( dataType, fDescriptor->fDataType.fOrigin, kAliHLTComponentDataTypefOriginSize );
			strcat( dataType, ":" );
			strncat( dataType, fDescriptor->fDataType.fID, kAliHLTComponentDataTypefIDsize );
				
			HLTError("Problem found with data block %d, fDataType = '%s',"
				" fPtr = %p and fSize = %u bytes."
				" Assuming this is a tracker DDL raw data block."
				" Problem: Found a data word 0x%8.8X for bus patch %d"
				" whose bits 29 or 30 are not zero.",
				fBlockNumber,
				&dataType[0],
				fDescriptor->fPtr,
				fDescriptor->fSize,
				data,
				fCurrentBusPatch
			);
			fDataProblems = true;
			return;
		}
		
		// Check if the data word + bus patch have been duplicated.
		for (UInt_t i = 0; i < fDigitCount; i++)
		{
			if (fDigits[i].fDataWord == data and fDigits[i].fBusPatchId == fCurrentBusPatch)
			{
				// create data type string.
				char dataType[kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2];
				memset( dataType, 0, kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2 );
				strncat( dataType, fDescriptor->fDataType.fOrigin, kAliHLTComponentDataTypefOriginSize );
				strcat( dataType, ":" );
				strncat( dataType, fDescriptor->fDataType.fID, kAliHLTComponentDataTypefIDsize );
				
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a tracker DDL raw data block."
					" Problem: Found a duplicate data word 0x%8.8X for bus patch %d.",
					fBlockNumber,
					&dataType[0],
					fDescriptor->fPtr,
					fDescriptor->fSize,
					data,
					fCurrentBusPatch
				);
				fDataProblems = true;
				return;
			}
		}
		
		// Add the data word + bus patch to the list of decoded digits.
		if (fDigitCount < fMaxDigits)
		{
			fDigits[fDigitCount].fBusPatchId = fCurrentBusPatch;
			fDigits[fDigitCount].fDataWord = data;
			fDigitCount++;
		}
	}
	
	/**
	 * Class for logging decoding errors when checking trigger raw DDL data.
	 * Used in the AliHLTMUONDataCheckerComponent::CheckRawDataBlock method.
	 */
	class AliHLTMUONTriggerDecoderHandler :
		public AliMUONTriggerDDLDecoderEventHandler, public AliHLTMUONDecoderHandler
	{
	public:
		// Methods inherited from AliMUONTriggerDDLDecoderEventHandler:
		
		/// Called for each new buffer.
		void OnNewBuffer(const void* buffer, UInt_t /*bufferSize*/)
		{
			fBufferStart = buffer;
		}
		
		/// Logs an error message if there was a decoding problem with the DDL payload.
		void OnError(ErrorCode code, const void* location)
		{
			LogError(code, location, *this);
		}
	};
	
} // end of namespace


bool AliHLTMUONDataCheckerComponent::CheckRawDataBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a raw data block.
	
	bool result = true;

	if (fIgnoreSpec)
	{
		HLTWarning("Not able to check DDL raw data if -ignorespec is specified.");
		return false;
	}
	
	bool ddl[22];
	AliHLTMUONUtils::UnpackSpecBits(block.fSpecification, ddl);
	
	// Check that only one DDL was marked in the specification.
	int ddlIndex = -1;
	for (int i = 0; i < 22; i++)
	{
		if (not ddl[i]) continue;
		
		if (ddlIndex == -1)
		{
			ddlIndex = i;
			continue;
		}
		
		HLTError("Problem found with data block %d, fDataType = '%s',"
			 " fPtr = %p and fSize = %u bytes."
			 " Assuming this is a DDL raw data block."
			 " Problem: The specification indicates multiple"
			 " DDL sources, DDL %d and %d.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			ddlIndex,
			i
		);
		result = false;
	}
	
	// Check the DDL common data header.
	AliHLTUInt32_t totalDDLSize = block.fSize;
	if (totalDDLSize < sizeof(AliRawDataHeader))
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			 " fPtr = %p and fSize = %u bytes."
			 " Assuming this is a DDL raw data block."
			 " Problem: The size of the data block is too short to contain"
			 " a valid common DDL data header. Size of buffer is only %d"
			 " bytes, but expected at least %d bytes.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			totalDDLSize,
			sizeof(AliRawDataHeader)
		);
		return false;
	}
	
	const AliRawDataHeader* header =
		reinterpret_cast<const AliRawDataHeader*>(block.fPtr);
	
	if (header->GetVersion() != 2)
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			 " fPtr = %p and fSize = %u bytes."
			 " Assuming this is a DDL raw data block."
			 " Problem: The common DDL data header indicates an"
			 " incorrect version number. Expected 2 but got %d.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			int( header->GetVersion() )
		);
		result = false;
	}
	
	if (header->fSize != 0xFFFFFFFF and header->fSize != block.fSize)
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			 " fPtr = %p and fSize = %u bytes."
			 " Assuming this is a DDL raw data block."
			 " Problem: The common DDL data header indicates an"
			 " incorrect DDL buffer size. Expected %d bytes but"
			 " size reported in header is %d bytes.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			block.fSize,
			header->fSize
		);
		result = false;
	}
	
	if (header->fSize != 0xFFFFFFFF and header->fSize != block.fSize)
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			 " fPtr = %p and fSize = %u bytes."
			 " Assuming this is a DDL raw data block."
			 " Problem: The common DDL data header indicates an"
			 " incorrect DDL buffer size. Expected %d bytes but"
			 " size reported in header is %d bytes.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			block.fSize,
			header->fSize
		);
		result = false;
	}
	
	// Check that the bits that should be zero in the CDH are infact zero.
	if ((header->fWord2 & 0x00C03000) != 0 or
	    (header->fEventID2 & 0xFF000000) != 0 or
	    (header->fStatusMiniEventID & 0xF0000000) != 0 or
	    (header->fROILowTriggerClassHigh & 0x0FFC0000) != 0
	   )
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			 " fPtr = %p and fSize = %u bytes."
			 " Assuming this is a DDL raw data block."
			 " Problem: The common DDL data header has non-zero"
			 " bits that are reserved and must be set to zero.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize
		);
		result = false;
	}
	
	AliHLTUInt32_t payloadSize = block.fSize - sizeof(AliRawDataHeader);
	const AliHLTUInt8_t* payload =
		reinterpret_cast<const AliHLTUInt8_t*>(header + 1);
	
	if (AliHLTMUONUtils::IsTriggerDDL(block.fSpecification))
	{
		bool scalarEvent = ((header->GetL1TriggerMessage() & 0x1) == 0x1);
		AliMUONTriggerDDLDecoder<AliHLTMUONTriggerDecoderHandler> decoder;
		decoder.ExitOnError(false);
		decoder.TryRecover(false);
		decoder.AutoDetectScalars(false);
		decoder.GetHandler().SetDescriptor(&block);
		decoder.GetHandler().SetBlockNumber(blockNumber);
		result = decoder.Decode(payload, payloadSize, scalarEvent);
	}
	else if (AliHLTMUONUtils::IsTrackerDDL(block.fSpecification))
	{
		AliMUONTrackerDDLDecoder<AliHLTMUONTrackerDecoderHandler> decoder;
		decoder.ExitOnError(false);
		decoder.TryRecover(false);
		decoder.SendDataOnParityError(true);
		decoder.AutoDetectTrailer(true);
		decoder.CheckForTrailer(true);
		decoder.GetHandler().SetDescriptor(&block);
		decoder.GetHandler().SetBlockNumber(blockNumber);
		result = decoder.Decode(payload, payloadSize);
		if (decoder.GetHandler().DataProblems()) result = false;
		
		if (FetchMappingStores() == 0)  // are stores loaded?
		{
			Bool_t warn = kFALSE;
			AliMpDDLStore* ddlStore = AliMpDDLStore::Instance(warn);
			
			// Check that the bus patch, manu ID and channel addresses are valid
			// for each raw data word.
			for (UInt_t i = 0; i < decoder.GetHandler().DigitCount(); i++)
			{
				UInt_t busPatchId = decoder.GetHandler().Digits()[i].fBusPatchId;
				UInt_t dataWord = decoder.GetHandler().Digits()[i].fDataWord;
				
				UShort_t manuId; UChar_t channelId; UShort_t adc;
				AliMUONTrackerDDLDecoderEventHandler::UnpackADC(
						dataWord, manuId, channelId, adc
					);
				
				// Check if the bus patch is valid.
				AliMpBusPatch* busPatch = ddlStore->GetBusPatch(busPatchId, warn);
				if (busPatch == NULL)
				{
					HLTError("Problem found with data block %d, fDataType = '%s',"
						 " fPtr = %p and fSize = %u bytes."
						 " Assuming this is a tracker DDL raw data block."
						 " Problem: Found a bus patch identifier %d that"
						 " is not valid.",
						blockNumber,
						DataType2Text(block.fDataType).c_str(),
						block.fPtr,
						block.fSize,
						busPatchId
					);
					result = false;
					continue;
				}
				
				// We can check that the bus patch is for the DDL
				// which is also indicated by the specification bits.
				if (not fIgnoreSpec and busPatch->GetDdlId() != ddlIndex)
				{
					HLTError("Problem found with data block %d, fDataType = '%s',"
						 " fPtr = %p and fSize = %u bytes."
						 " Assuming this is a tracker DDL raw data block."
						 " Problem: Found a bus patch identifier %d for"
						 " DDL %d, but the data block specification 0x%8.8X"
						 " indicates a different DDL of %d.",
						blockNumber,
						DataType2Text(block.fDataType).c_str(),
						block.fPtr,
						block.fSize,
						busPatchId,
						busPatch->GetDdlId(),
						block.fSpecification,
						ddlIndex
					);
					result = false;
					continue;
				}
				
				// Check if the MANU ID is valid.
				if (not busPatch->HasManu(manuId))
				{
					HLTError("Problem found with data block %d, fDataType = '%s',"
						 " fPtr = %p and fSize = %u bytes."
						 " Assuming this is a tracker DDL raw data block."
						 " Problem: Found a MANU identifier %d on bus patch %d"
						 " that is not valid.",
						blockNumber,
						DataType2Text(block.fDataType).c_str(),
						block.fPtr,
						block.fSize,
						manuId,
						busPatchId
					);
					result = false;
					continue;
				}
				
				// Now try to fetch the detector element to check the MANU channel.
				AliMpDetElement* de = ddlStore->GetDetElement(busPatch->GetDEId(), warn);
				if (de == NULL)
				{
					HLTError("Problem found with data block %d, fDataType = '%s',"
						 " fPtr = %p and fSize = %u bytes."
						 " Assuming this is a tracker DDL raw data block."
						 " Problem: Found a bus patch identifier %d that"
						 " does not correspond to a detector element.",
						blockNumber,
						DataType2Text(block.fDataType).c_str(),
						block.fPtr,
						block.fSize,
						busPatchId
					);
					result = false;
					continue;
				}
				
				if (not de->IsConnectedChannel(manuId, channelId))
				{
					// Just a warning because this is marked not
					// to be an error in the AliMUONDigitMaker.
					HLTWarning("Problem found with data block %d, fDataType = '%s',"
						 " fPtr = %p and fSize = %u bytes."
						 " Assuming this is a tracker DDL raw data block."
						 " Problem: Found a channel with address %d on"
						 " MANU ID %d and bus patch %d that is not connected.",
						blockNumber,
						DataType2Text(block.fDataType).c_str(),
						block.fPtr,
						block.fSize,
						channelId,
						manuId,
						busPatchId
					);
					result = false;
					continue;
				}
				
				// Need to also load the correct segmentation to check the channel.
				const AliMpVSegmentation* seg =
					AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(
						busPatch->GetDEId(), manuId
					);
				if (seg == NULL)
				{
					HLTError("Could not load segmentation for detector element %d"
						 " and MANU ID %d.",
						busPatch->GetDEId(), manuId
					);
					result = false;
					continue;
				}
				
#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
				AliMpPad pad = seg->PadByLocation(manuId, channelId, warn);
#else // old AliMpPad functionality < r 31742
				AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId, channelId), warn);
#endif //HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
				if (not pad.IsValid())
				{
					HLTError("Problem found with data block %d, fDataType = '%s',"
						 " fPtr = %p and fSize = %u bytes."
						 " Assuming this is a tracker DDL raw data block."
						 " Problem: Found a channel with address %d on"
						 " MANU ID %d and bus patch %d that is not valid.",
						blockNumber,
						DataType2Text(block.fDataType).c_str(),
						block.fPtr,
						block.fSize,
						channelId,
						manuId,
						busPatchId
					);
					result = false;
					continue;
				}
			}
		}
		else
		{
			HLTWarning("Cannot check if the bus patch IDs, MANU ID and"
				" channel addresses for DDL raw data are valid without"
				" being able to load the mapping from CDB."
			);
			result = false;
		}
	}
	else
	{
		HLTError("Problem found with data block %d, fDataType = '%s',"
			 " fPtr = %p and fSize = %u bytes."
			 " Assuming this is a DDL raw data block."
			 " Problem: The specification does not contain a valid pattern,"
			 " received 0x%8.8X for the specification.",
			blockNumber,
			DataType2Text(block.fDataType).c_str(),
			block.fPtr,
			block.fSize,
			block.fSpecification
		);
		result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckTriggerRecordsBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a trigger records block.

	bool result = true;
	const char* name = "trigger records";
	
	if (not fIgnoreSpec)
	{
		if (not IsFromTriggerOnly(block, blockNumber, name))
			result = false;
	}
	
	AliHLTMUONTriggerRecordsBlockReader inblock(block.fPtr, block.fSize);
	if (not CheckBlockIntegrity(block, blockNumber, inblock, name))
		return false;
	
	bool ddl[22];
	AliHLTMUONUtils::UnpackSpecBits(block.fSpecification, ddl);
	
	// Min and max allowed chamber numbers for hits:
	Int_t minCh = AliMUONConstants::NCh() - AliMUONConstants::NTriggerCh();
	Int_t maxCh = AliMUONConstants::NCh() - 1;
	
	for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
	{
		// Check that each hit in each trigger record has a reasonable coordinate.
		AliHLTMUONParticleSign sign;
		bool hitset[4];
		AliHLTMUONUtils::UnpackTriggerRecordFlags(inblock[i].fFlags, sign, hitset);
	
		for (Int_t j = 0; j < 4; j++)  // loop over 4 trigger chamber hits.
		{
			if (not hitset[i]) continue; // ignore hits that are not initialised.
			bool hitOk = IsHitCoordinateOk(
					block, blockNumber, name, i, inblock[i].fHit[j],
					minCh, maxCh, j+10, ddl
				);
			if (not hitOk) result = false;
		}
		
		// We can also check the momentum vector.
		bool momOk = IsMomentumVectorOk(
				block, blockNumber, name, i,
				inblock[i].fPx, inblock[i].fPy, inblock[i].fPz
			);
		if (not momOk) result = false;
	}
	
	// Need to check that no entries have duplicated data but with a different
	// ID number.
	for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
	{
		AliHLTMUONTriggerRecordStruct ti = inblock[i];
		ti.fId = -1;
		for (AliHLTUInt32_t j = i+1; j < inblock.Nentries(); j++)
		{
			AliHLTMUONTriggerRecordStruct tj = inblock[j];
			tj.fId = ti.fId;
			
			if (ti == tj)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: The trigger records %d and %d contain the"
					" same data. The data might have been duplicated.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					i, j
				);
				result = false;
			}
		}
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckTrigRecsDebugBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a trigger records debug block.
	
	bool result = true;
	const char* name = "trigger records debug information";
	
	if (not fIgnoreSpec)
	{
		if (not IsFromTriggerOnly(block, blockNumber, name))
			result = false;
	}
	
	AliHLTMUONTrigRecsDebugBlockReader inblock(block.fPtr, block.fSize);
	if (not CheckBlockIntegrity(block, blockNumber, inblock, name))
		return false;
	
	bool ddl[22];
	AliHLTMUONUtils::UnpackSpecBits(block.fSpecification, ddl);
	
	// Check that each detector element ID is valid and the corresponding DDL
	// bit is set in the data block specification.
	if (FetchMappingStores() == 0)  // are stores loaded?
	{
		Bool_t warn = kFALSE;
		AliMpDEStore* store = AliMpDEStore::Instance(warn);
		for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
		for (AliHLTUInt32_t j = 0; j < 4; j++)
		{
			const AliHLTMUONTrigRecInfoStruct& trig = inblock[i];
			AliMpDetElement* de = store->GetDetElement(trig.fDetElemId[j], warn);
			if (de == NULL)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: The detector element number %d on chamber"
					 " %d for trigger record debug structure %d is not valid.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					trig.fDetElemId[j],
					j+11,
					i
				);
				result = false;
				continue;
			}
			
			// Check that the chamber number from the detector element number
			// has the expected value.
			Int_t chamber = AliMpDEManager::GetChamberId(trig.fDetElemId[j], warn);
			if (chamber != Int_t(j+10))
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: The detector element number %d for trigger"
					 " record debug structure %d, corresponds to chamber"
					 " %d, but we expected a hit for chamber %d.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					trig.fDetElemId[j],
					i,
					chamber+1,
					j+11
				);
				result = false;
			}
			
			if (fIgnoreSpec) continue;
			if (0 <= de->GetDdlId() and de->GetDdlId() < 22 and not ddl[de->GetDdlId()])
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: The detector element number %d for trigger"
					 " record %d corresponds to DDL number %d, but the"
					 " data block specification 0x%8.8X does not have the"
					 " corresponding bit set.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					trig.fDetElemId[j],
					i,
					de->GetDdlId(),
					block.fSpecification
				);
				result = false;
			}
		}
	}
	else
	{
		HLTWarning("Cannot check trigger record debug information without"
			" being able to load the mapping from CDB."
		);
		result = false;
	}
	
	// Need to check that no entries have duplicated data but with a different
	// ID number.
	for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
	{
		AliHLTMUONTrigRecInfoStruct ti = inblock[i];
		ti.fTrigRecId = -1;
		for (AliHLTUInt32_t j = i+1; j < inblock.Nentries(); j++)
		{
			AliHLTMUONTrigRecInfoStruct tj = inblock[j];
			tj.fTrigRecId = ti.fTrigRecId;
			
			if (ti == tj)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: The trigger record debug information"
					" structures %d and %d contain the same data."
					" The data might have been duplicated.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					i, j
				);
				result = false;
			}
		}
		
		// Can also check that the value of the fZmiddle and fBl.
		bool paramsOk = AreMomentumCalcParamsOk(
				block, blockNumber, name, i, ti.fZmiddle, ti.fBl
			);
		if (not paramsOk) result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckRecHitsBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a reconstructed hits block.

	bool result = true;
	const char* name = "reconstructed hits";
	
	if (not fIgnoreSpec)
	{
		if (not IsFromTrackerOnly(block, blockNumber, name))
			result = false;
	}
	
	AliHLTMUONRecHitsBlockReader inblock(block.fPtr, block.fSize);
	if (not CheckBlockIntegrity(block, blockNumber, inblock, name))
		return false;
	
	bool ddl[22];
	AliHLTMUONUtils::UnpackSpecBits(block.fSpecification, ddl);
	
	// Check that each hit has a reasonable coordinate.
	Int_t minCh = 0;
	Int_t maxCh = AliMUONConstants::NTrackingCh() - 1;
	for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
	{
		bool hitOk = IsHitCoordinateOk(
				block, blockNumber, name, i, inblock[i],
				minCh, maxCh, -1, ddl
			);
		if (not hitOk) result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckClustersBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a clusters block.

	bool result = true;
	const char* name = "clusters";
	
	if (not fIgnoreSpec)
	{
		if (not IsFromTrackerOnly(block, blockNumber, name))
			result = false;
	}
	
	AliHLTMUONClustersBlockReader inblock(block.fPtr, block.fSize);
	if (not CheckBlockIntegrity(block, blockNumber, inblock, name))
		return false;
	
	bool ddl[22];
	AliHLTMUONUtils::UnpackSpecBits(block.fSpecification, ddl);
	
	if (FetchMappingStores() == 0)  // are stores loaded?
	{
		Bool_t warn = kFALSE;
		AliMpDEStore* store = AliMpDEStore::Instance(warn);
		for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
		{
			const AliHLTMUONClusterStruct& cluster = inblock[i];
			
			// Check that the detector element ID is valid.
			AliMpDetElement* de = store->GetDetElement(cluster.fDetElemId, warn);
			if (de == NULL)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: The detector element number %d for cluster"
					 " %d is not valid.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					cluster.fDetElemId,
					i
				);
				result = false;
				continue;
			}
			
			// Check that the chamber number found from the hit coordinate and
			// that from the detector element number are the same.
			Int_t chamberFromHit = AliMUONConstants::ChamberNumber(cluster.fHit.fZ, warn);
			Int_t chamberFromDE = AliMpDEManager::GetChamberId(cluster.fDetElemId, warn);
			if (chamberFromHit != chamberFromDE)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: The detector element number %d for"
					 " cluster %d, corresponds to chamber %d, but"
					 " found a different chamber number %d for the"
					 " corresponding hit coordinate {x = %f, y = %f,"
					 " z = %f}.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					cluster.fDetElemId,
					i,
					chamberFromDE+1,
					chamberFromHit+1,
					cluster.fHit.fX,
					cluster.fHit.fY,
					cluster.fHit.fZ
				);
				result = false;
			}
			
			// Make sure the corresponding DDL bit is set in the data
			// block specification.
			if (fIgnoreSpec) continue;
			if (0 <= de->GetDdlId() and de->GetDdlId() < 22 and not ddl[de->GetDdlId()])
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: The detector element number %d for cluster"
					 " %d corresponds to DDL number %d, but the data"
					 " block specification 0x%8.8X does not have the"
					 " corresponding bit set.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					cluster.fDetElemId,
					i,
					de->GetDdlId(),
					block.fSpecification
				);
				result = false;
			}
			
			// Check that the total cluster charge is a reasonable value.
			if (cluster.fCharge < 0 and 1e4 < cluster.fCharge)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: The total charge %f for the cluster"
					 " %d is not in a reasonable range [0..1e4].",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					cluster.fCharge,
					i
				);
				result = false;
				continue;
			}
		}
	}
	else
	{
		HLTWarning("Cannot check cluster information without being able"
			" to load the mapping from CDB."
		);
		result = false;
	}
	
	// Min and max chamber numbers allowed for the cluster hits.
	Int_t minCh = 0;
	Int_t maxCh = AliMUONConstants::NTrackingCh() - 1;
	
	for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
	{
		// Need to check that no cluster data has duplicated data but with
		// a different ID number.
		AliHLTMUONClusterStruct ci = inblock[i];
		ci.fId = -1;
		for (AliHLTUInt32_t j = i+1; j < inblock.Nentries(); j++)
		{
			AliHLTMUONClusterStruct cj = inblock[j];
			cj.fId = ci.fId;
			
			if (ci == cj)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: The cluster structures %d and %d contain"
					" the same data. The data might have been duplicated.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					i, j
				);
				result = false;
			}
		}
		
		// Check that the hit structure in the cluster corresponds
		// to a tracker chamber.
		bool hitOk = IsHitCoordinateOk(
				block, blockNumber, name, i, ci.fHit,
				minCh, maxCh, -1, ddl
			);
		if (not hitOk) result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckChannelsBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a channels block.

	bool result = true;
	const char* name = "channels";
	
	if (not fIgnoreSpec)
	{
		if (not IsFromTrackerOnly(block, blockNumber, name))
			result = false;
	}
	
	AliHLTMUONChannelsBlockReader inblock(block.fPtr, block.fSize);
	if (not CheckBlockIntegrity(block, blockNumber, inblock, name))
		return false;
	
	bool ddl[22];
	AliHLTMUONUtils::UnpackSpecBits(block.fSpecification, ddl);
	
	if (FetchMappingStores() == 0)  // are stores loaded?
	{
		Bool_t warn = kFALSE;
		AliMpDDLStore* store = AliMpDDLStore::Instance(warn);
		
		for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
		{
			const AliHLTMUONChannelStruct& channel = inblock[i];
			
			// Check if the bus patch is valid.
			AliMpBusPatch* busPatch = store->GetBusPatch(channel.fBusPatch, warn);
			if (busPatch == NULL)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: Found a bus patch identifier %d that"
					" is not valid.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					channel.fBusPatch
				);
				result = false;
				continue;
			}
			
			// We can check that the bus patch is for a DDL
			// which is also indicated by the specification bits.
			if (not fIgnoreSpec and (
			     not (0 <= busPatch->GetDdlId() and busPatch->GetDdlId() < 20)
			     or  (0 <= busPatch->GetDdlId() and busPatch->GetDdlId() < 20
			          and not ddl[busPatch->GetDdlId()])
			   ))
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: Found a bus patch identifier %d for"
					 " DDL %d, but the data block specification 0x%8.8X"
					 " does not have the corresponding bit set.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					channel.fBusPatch,
					busPatch->GetDdlId(),
					block.fSpecification
				);
				result = false;
				continue;
			}
			
			// Check if the MANU ID is valid.
			if (not busPatch->HasManu(channel.fManu))
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: Found a MANU identifier %d on bus patch %d"
					 " that is not valid.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					channel.fManu,
					channel.fBusPatch
				);
				result = false;
				continue;
			}
			
			// Now try to fetch the detector element to check the MANU channel.
			AliMpDetElement* de = store->GetDetElement(busPatch->GetDEId(), warn);
			if (de == NULL)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: Found a bus patch identifier %d that"
					 " does not correspond to a detector element.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					channel.fBusPatch
				);
				result = false;
				continue;
			}
			
			if (not de->IsConnectedChannel(channel.fManu, channel.fChannelAddress))
			{
				// Just a warning because this is marked not
				// to be an error in the AliMUONDigitMaker.
				HLTWarning("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: Found a channel with address %d on"
					" MANU ID %d and bus patch %d that is not connected.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					channel.fChannelAddress,
					channel.fManu,
					channel.fBusPatch
				);
				result = false;
				continue;
			}
			
			// Need to also load the correct segmentation to check the channel.
			const AliMpVSegmentation* seg =
				AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(
					busPatch->GetDEId(), channel.fManu
				);
			if (seg == NULL)
			{
				HLTError("Could not load segmentation for detector element %d"
					 " and MANU ID %d.",
					busPatch->GetDEId(), channel.fManu
				);
				result = false;
				continue;
			}
			
			AliMpPad pad = seg->PadByLocation(
#ifndef HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					channel.fManu, channel.fChannelAddress,
#else // old AliMpPad functionality < r 31742
					AliMpIntPair(channel.fManu, channel.fChannelAddress),
#endif //HAVE_NOT_MUON_ALIMPPAD_GETPOSITION
					warn
				);
			if (not pad.IsValid())
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: Found a channel with address %d on"
					" MANU ID %d and bus patch %d that is not valid.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					channel.fChannelAddress,
					channel.fManu,
					channel.fBusPatch
				);
				result = false;
				continue;
			}
		}
	}
	else
	{
		HLTWarning("Cannot check channel information without being able"
			" to load the mapping from CDB."
		);
		result = false;
	}
	
	// Need to check that no channel data has duplicated data but with
	// a different cluster ID number.
	for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
	{
		AliHLTMUONChannelStruct ci = inblock[i];
		ci.fClusterId = -1;
		for (AliHLTUInt32_t j = i+1; j < inblock.Nentries(); j++)
		{
			AliHLTMUONChannelStruct cj = inblock[j];
			cj.fClusterId = ci.fClusterId;
			
			if (ci == cj)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: The channel structures %d and %d contain"
					" the same data. The data might have been duplicated.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					i, j
				);
				result = false;
			}
		}
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckMansoTracksBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a Manso tracks block.

	bool result = true;
	const char* name = "Manso tracks";
	
	if (not fIgnoreSpec)
	{
		if (not IsSpecificationValid(block, blockNumber, name))
			result = false;
	}
	
	AliHLTMUONMansoTracksBlockReader inblock(block.fPtr, block.fSize);
	if (not CheckBlockIntegrity(block, blockNumber, inblock, name))
		return false;
	
	bool ddl[22];
	AliHLTMUONUtils::UnpackSpecBits(block.fSpecification, ddl);
	
	for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
	{
		// Need to check that no entries have duplicated data but with
		// a different track ID number.
		AliHLTMUONMansoTrackStruct ti = inblock[i];
		ti.fId = -1;
		for (AliHLTUInt32_t j = i+1; j < inblock.Nentries(); j++)
		{
			AliHLTMUONMansoTrackStruct tj = inblock[j];
			tj.fId = ti.fId;
			
			if (ti == tj)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: The Manso tracks %d and %d contain the"
					" same data. The data might have been duplicated.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					i, j
				);
				result = false;
			}
		}
		
		bool trackOk = IsMansoTrackOk(block, blockNumber, name, i, ti, ddl);
		if (not trackOk) result = false;
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckMansoCandidatesBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a Manso candidates block.

	bool result = true;
	const char* name = "Manso track candidates";
	
	if (not fIgnoreSpec)
	{
		if (not IsSpecificationValid(block, blockNumber, name))
			result = false;
	}
	
	AliHLTMUONMansoCandidatesBlockReader inblock(block.fPtr, block.fSize);
	if (not CheckBlockIntegrity(block, blockNumber, inblock, name))
		return false;
	
	bool ddl[22];
	AliHLTMUONUtils::UnpackSpecBits(block.fSpecification, ddl);
	
	for (AliHLTUInt32_t i = 0; i < inblock.Nentries(); i++)
	{
		// Need to check that no entries have duplicated data but with a
		// different track ID number.
		AliHLTMUONMansoCandidateStruct ti = inblock[i];
		ti.fTrack.fId = -1;
		for (AliHLTUInt32_t j = i+1; j < inblock.Nentries(); j++)
		{
			AliHLTMUONMansoCandidateStruct tj = inblock[j];
			tj.fTrack.fId = ti.fTrack.fId;
			
			if (ti == tj)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: The Manso track candidates %d and %d"
					" contain the same data."
					" The data might have been duplicated.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					i, j
				);
				result = false;
			}
		}
		
		bool trackOk = IsMansoTrackOk(block, blockNumber, name, i, ti.fTrack, ddl);
		if (not trackOk) result = false;
		
		// Check that each ROI has a centre point somewhere on the correct
		// corresponding chamber and that the Radius is not bigger thant
		// the diameter of the chamber which would be pointless.
		for (AliHLTInt32_t j = 0; j < 4; j++)
		{
			if (ti.fRoI[j].fRadius == -1) continue; // Ignore invalid ROIs
			
			Int_t chamber = AliMUONConstants::ChamberNumber(
					ti.fRoI[j].fZ, false  // false = do not warn.
				);
			if (chamber != j+6)
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: The region of interest on chamber %d for"
					" Manso track candidate %d has a z-coordinate of %f"
					" cm that does not correspond to that chamber.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					j+6+1,
					i,
					ti.fRoI[j].fZ
				);
				result = false;
			}
			
			double x = ti.fRoI[j].fX;
			double y = ti.fRoI[j].fY;
			double r = sqrt(x*x + y*y);
			if (r > AliMUONConstants::Dmax((j+6)/2))
			{
				// Just a warning since this is not a data integrity problem
				// but rather just a data sanity problem.
				HLTWarning("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: The region of interest coordinate {x = %f,"
					" y = %f} cm on chamber %d for Manso track candidate %d"
					" does not correspond to that chamber.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					ti.fRoI[j].fX,
					ti.fRoI[j].fY,
					j+6+1,
					i
				);
				result = false;
			}
			
			if (ti.fRoI[j].fRadius > AliMUONConstants::Dmax((j+6)/2))
			{
				// Just a warning since this is not a data integrity problem
				// but rather just a data sanity problem.
				HLTWarning("Problem found with data block %d, fDataType = '%s',"
					" fPtr = %p and fSize = %u bytes."
					" Assuming this is a %s data block."
					" Problem: The region of interest radius of %f cm"
					" on chamber %d for Manso track candidate %d"
					" is bigger than the chamber diameter %f cm.",
					blockNumber,
					DataType2Text(block.fDataType).c_str(),
					block.fPtr,
					block.fSize,
					name,
					ti.fRoI[j].fRadius,
					j+6+1,
					i,
					AliMUONConstants::Dmax((j+6)/2)
				);
				result = false;
			}
		}
	}
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckSinglesDecisionBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a single tracks trigger decision block.

	bool result = true;
	const char* name = "singles decision";
	
	if (not fIgnoreSpec)
	{
		if (not IsSpecificationValid(block, blockNumber, name))
			result = false;
	}
	
	AliHLTMUONSinglesDecisionBlockReader inblock(block.fPtr, block.fSize);
	if (not CheckBlockIntegrity(block, blockNumber, inblock, name))
		return false;
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::CheckPairsDecisionBlock(
		const AliHLTComponentBlockData& block,
		AliHLTUInt32_t blockNumber
	) const
{
	/// Checks the validity of a track pairs trigger decision block.

	bool result = true;
	const char* name = "pairs decision";
	
	if (not fIgnoreSpec)
	{
		if (not IsSpecificationValid(block, blockNumber, name))
			result = false;
	}
	
	AliHLTMUONPairsDecisionBlockReader inblock(block.fPtr, block.fSize);
	if (not CheckBlockIntegrity(block, blockNumber, inblock, name))
		return false;
	
	return result;
}


bool AliHLTMUONDataCheckerComponent::AreMomentaCompatible(
		AliHLTFloat32_t px1,
		AliHLTFloat32_t py1,
		AliHLTFloat32_t pz1,
		AliHLTFloat32_t px2,
		AliHLTFloat32_t py2,
		AliHLTFloat32_t pz2
	) const
{
	/// Checks to see if the two momenta vectors are compatible or not.
	/// The vectors should not have an angle more than 10 degrees between
	/// them and their magnitudes should not be more different than 50%.
	
	double p1 = sqrt(double(px1)*double(px1) + double(py1)*double(py1) + double(pz1)*double(pz1));
	double p2 = sqrt(double(px2)*double(px2) + double(py2)*double(py2) + double(pz2)*double(pz2));
	if (p1 == 0 and p2 == 0) return true;
	if (fabs(p1 - p2) / ((p1 + p2)*0.5) > 0.5) return false;
	double denom = p1 * p2;
	if (denom == 0) return false;
	double ratio = (double(px1)*double(px2) + double(py1)*double(py2) + double(pz1)*double(pz2)) / denom;
	if (ratio < -1) return true;
	if (ratio > 1) return true;
	double angle = acos(ratio);
	if (angle > 3.14159265358979323846 * 10. / 180.) return false;
	return true;
}


bool AliHLTMUONDataCheckerComponent::IsScalarTooLarge(
		const AliHLTComponentBlockData* block,
		AliHLTUInt32_t blockNumber,
		const char* blockTypeName,
		const char* scalarName,
		AliHLTUInt32_t scalarValue,
		AliHLTUInt32_t totalTrackCount
	) const
{
	/// Checks if the scalar value is larger than the number of Manso
	/// tracks in the event.

	if (scalarValue > totalTrackCount)
	{
		HLTError("Problem found with %s trigger decision"
			" data block %d, fDataType = '%s', fPtr = %p and"
			" fSize = %u bytes."
			" Problem: The %s scalar with value %d is larger"
			" than the total number of Manso tracks found for the"
			" event (%d tracks).",
			blockTypeName,
			blockNumber,
			DataType2Text(block->fDataType).c_str(),
			block->fPtr,
			block->fSize,
			scalarName,
			scalarValue,
			totalTrackCount
		);
		return true;
	}
	else
	{
		return false;
	}
}


bool AliHLTMUONDataCheckerComponent::IsScalarALargerThanB(
		const AliHLTComponentBlockData* block,
		AliHLTUInt32_t blockNumber,
		const char* blockTypeName,
		const char* scalarAName,
		AliHLTUInt32_t scalarAValue,
		const char* scalarBName,
		AliHLTUInt32_t scalarBValue
	) const
{
	/// Checks if the scalar value is larger than the number of Manso
	/// tracks in the event.

	if (scalarAValue > scalarBValue)
	{
		HLTError("Problem found with %s trigger decision"
			" data block %d, fDataType = '%s', fPtr = %p and"
			" fSize = %u bytes."
			" Problem: The %s scalar with value %d is larger"
			" than scalar %s with value %d, but is should not be.",
			blockTypeName,
			blockNumber,
			DataType2Text(block->fDataType).c_str(),
			block->fPtr,
			block->fSize,
			scalarAName,
			scalarAValue,
			scalarBName,
			scalarBValue
		);
		return true;
	}
	else
	{
		return false;
	}
}


void AliHLTMUONDataCheckerComponent::MarkBlock(
		const AliHLTComponentBlockData* blocks,
		bool* blockOk,
		AliHLTUInt32_t blockCount,
		const AliHLTComponentBlockData* blockToMark
	) const
{
	/// Tries to find the 'blockToMark' in the list of blocks and sets the
	/// corresponding 'blockOk' flag to false.
	
	for (AliHLTUInt32_t i = 0; i < blockCount; i++)
	{
		if (&blocks[i] == blockToMark)
		{
			blockOk[i] = false;
			return;
		}
	}
}


void AliHLTMUONDataCheckerComponent::MakeGlobalChecks(
		const AliHLTComponentBlockData* blocks,
		bool* blockOk,
		AliHLTUInt32_t blockCount,
		const AliHLTComponentBlockData** trigRecBlocks,
		AliHLTUInt32_t trigRecBlocksCount,
		const AliHLTComponentBlockData** trigRecDebugBlocks,
		AliHLTUInt32_t trigRecDebugBlocksCount,
		const AliHLTComponentBlockData** hitBlocks,
		AliHLTUInt32_t hitBlocksCount,
		const AliHLTComponentBlockData** clusterBlocks,
		AliHLTUInt32_t clusterBlocksCount,
		const AliHLTComponentBlockData** channelBlocks,
		AliHLTUInt32_t channelBlocksCount,
		const AliHLTComponentBlockData** mansoTrackBlocks,
		AliHLTUInt32_t mansoTrackBlocksCount,
		const AliHLTComponentBlockData** mansoCandidateBlocks,
		AliHLTUInt32_t mansoCandidateBlocksCount,
		const AliHLTComponentBlockData** singleDecisionBlocks,
		AliHLTUInt32_t singleDecisionBlocksCount,
		const AliHLTComponentBlockData** pairDecisionBlocks,
		AliHLTUInt32_t pairDecisionBlocksCount
	) const
{
	/// The following set of global checks are performed:
	/// 1) Checks if all the ID numbers in all the blocks are unique.
	/// 2) Check if all the structures are unique up to their ID numbers,
	/// that it, make sure there are no structures with the same data but
	/// for a different ID number.
	/// 3) Check if the reference ID numbers are correct, i.e. are the
	/// trigger record ID numbers in the track structures found in any of
	/// the trigger record data blocks.
	/// 4) Do the number of channels claimed in the cluster correspond to
	/// the number of channel structures.
	/// 5) Check that the momentum vectors between the Manso tracks and
	/// the corresponding trigger record are compatible.
	/// 6) Check that the trigger decision scalars are reasonable.
	/// 7) Check that the detector element IDs are the same between rec
	/// hits and clusters / trigger record debug blocks.
	
	// Check if all the trigger record identifiers and data are unique.
	for (AliHLTUInt32_t bi = 0; bi < trigRecBlocksCount; bi++)
	{
		AliHLTMUONTriggerRecordsBlockReader inblocki(trigRecBlocks[bi]->fPtr, trigRecBlocks[bi]->fSize);
		if (not inblocki.BufferSizeOk()) continue;
		for (AliHLTUInt32_t bj = bi+1; bj < trigRecBlocksCount; bj++)
		{
			AliHLTMUONTriggerRecordsBlockReader inblockj(trigRecBlocks[bj]->fPtr, trigRecBlocks[bj]->fSize);
			if (not inblockj.BufferSizeOk()) continue;
			
			for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
			for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
			{
				if (inblocki[i].fId == inblockj[j].fId)
				{
					HLTError("Problem found with trigger record data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes,"
						" and trigger record data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: Trigger record %d in block %d and entry"
						" %d in block %d have the same identfier, but they"
						" should be unique.",
						bi,
						DataType2Text(trigRecBlocks[bi]->fDataType).c_str(),
						trigRecBlocks[bi]->fPtr,
						trigRecBlocks[bi]->fSize,
						bj,
						DataType2Text(trigRecBlocks[bj]->fDataType).c_str(),
						trigRecBlocks[bj]->fPtr,
						trigRecBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, trigRecBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, trigRecBlocks[bj]);
				}
				
				AliHLTMUONTriggerRecordStruct a = inblocki[i];
				AliHLTMUONTriggerRecordStruct b = inblockj[j];
				a.fId = b.fId = -1;
				if (a == b)
				{
					HLTError("Problem found with trigger record data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes,"
						" and trigger record data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: Trigger record %d in block %d and entry"
						" %d in block %d have the same data."
						" The data may have been duplicated.",
						bi,
						DataType2Text(trigRecBlocks[bi]->fDataType).c_str(),
						trigRecBlocks[bi]->fPtr,
						trigRecBlocks[bi]->fSize,
						bj,
						DataType2Text(trigRecBlocks[bj]->fDataType).c_str(),
						trigRecBlocks[bj]->fPtr,
						trigRecBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, trigRecBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, trigRecBlocks[bj]);
				}
			}
		}
	}
	
	for (AliHLTUInt32_t bi = 0; bi < trigRecDebugBlocksCount; bi++)
	{
		AliHLTMUONTrigRecsDebugBlockReader inblocki(trigRecDebugBlocks[bi]->fPtr, trigRecDebugBlocks[bi]->fSize);
		if (not inblocki.BufferSizeOk()) continue;
		
		for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
		{
			// Check if all the trigger record IDs in the debug information structures exist.
			bool found = false;
			
			for (AliHLTUInt32_t bj = 0; bj < trigRecBlocksCount and not found; bj++)
			{
				AliHLTMUONTriggerRecordsBlockReader inblockj(trigRecBlocks[bj]->fPtr, trigRecBlocks[bj]->fSize);
				if (not inblockj.BufferSizeOk()) continue;
				
				for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
				{
					if (inblocki[i].fTrigRecId == inblockj[j].fId)
					{
						found = true;
						
						// Since we found the corresponding trigger record,
						// check if the detector element IDs are the same.
						bool deOk = CheckDetElemIds(
								*trigRecDebugBlocks[bi], bi, i, inblocki[i],
								*trigRecBlocks[bj], bj, j, inblockj[j]
							);
						if (not deOk)
						{
							MarkBlock(blocks, blockOk, blockCount, trigRecDebugBlocks[bi]);
							MarkBlock(blocks, blockOk, blockCount, trigRecBlocks[bj]);
						}
						
						break;
					}
				}
			}
			
			if (not found)
			{
				HLTError("Problem found with trigger record debug information"
					" data block %d, fDataType = '%s', fPtr = %p and"
					" fSize = %u bytes."
					" Problem with entry %d in block: The trigger record"
					" identifier %d does not exist in any trigger record"
					" data block.",
					bi,
					DataType2Text(trigRecDebugBlocks[bi]->fDataType).c_str(),
					trigRecDebugBlocks[bi]->fPtr,
					trigRecDebugBlocks[bi]->fSize,
					i, inblocki[i].fTrigRecId
				);
				MarkBlock(blocks, blockOk, blockCount, trigRecDebugBlocks[bi]);
			}
		}
		
		// Check if all the trigger record debug information structures are unique.
		for (AliHLTUInt32_t bj = bi+1; bj < trigRecDebugBlocksCount; bj++)
		{
			AliHLTMUONTrigRecsDebugBlockReader inblockj(trigRecDebugBlocks[bj]->fPtr, trigRecDebugBlocks[bj]->fSize);
			if (not inblockj.BufferSizeOk()) continue;
			
			for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
			for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
			{
				AliHLTMUONTrigRecInfoStruct a = inblocki[i];
				AliHLTMUONTrigRecInfoStruct b = inblockj[j];
				a.fTrigRecId = b.fTrigRecId = -1;
				if (a == b)
				{
					HLTError("Problem found with trigger record debug information"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and trigger record debug"
						" information data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The trigger record debug inforamtion"
						" structure %d in block %d and entry"
						" %d in block %d have the same data."
						" The data may have been duplicated.",
						bi,
						DataType2Text(trigRecDebugBlocks[bi]->fDataType).c_str(),
						trigRecDebugBlocks[bi]->fPtr,
						trigRecDebugBlocks[bi]->fSize,
						bj,
						DataType2Text(trigRecDebugBlocks[bj]->fDataType).c_str(),
						trigRecDebugBlocks[bj]->fPtr,
						trigRecDebugBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, trigRecDebugBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, trigRecDebugBlocks[bj]);
				}
			}
		}
	}
	
	// Check that all the reconstructed hits are unique.
	for (AliHLTUInt32_t bi = 0; bi < hitBlocksCount; bi++)
	{
		AliHLTMUONRecHitsBlockReader inblocki(hitBlocks[bi]->fPtr, hitBlocks[bi]->fSize);
		if (not inblocki.BufferSizeOk()) continue;
		for (AliHLTUInt32_t bj = bi+1; bj < hitBlocksCount; bj++)
		{
			AliHLTMUONRecHitsBlockReader inblockj(hitBlocks[bj]->fPtr, hitBlocks[bj]->fSize);
			if (not inblockj.BufferSizeOk()) continue;
			
			for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
			for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
			{
				if (inblocki[i] == inblockj[j])
				{
					HLTError("Problem found with reconstructed hit data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes,"
						" and reconstructed hit data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: Reconstructed hit %d in block %d and entry"
						" %d in block %d are the same, but all hits"
						" should be unique.",
						bi,
						DataType2Text(hitBlocks[bi]->fDataType).c_str(),
						hitBlocks[bi]->fPtr,
						hitBlocks[bi]->fSize,
						bj,
						DataType2Text(hitBlocks[bj]->fDataType).c_str(),
						hitBlocks[bj]->fPtr,
						hitBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, hitBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, hitBlocks[bj]);
				}
			}
		}
	}
	
	for (AliHLTUInt32_t bi = 0; bi < clusterBlocksCount; bi++)
	{
		AliHLTMUONClustersBlockReader inblocki(clusterBlocks[bi]->fPtr, clusterBlocks[bi]->fSize);
		if (not inblocki.BufferSizeOk()) continue;
		
		for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
		{
			// Check if all the reconstructed hit coordinates in the cluster structures exist.
			bool found = false;
			
			for (AliHLTUInt32_t bj = 0; bj < hitBlocksCount and not found; bj++)
			{
				AliHLTMUONRecHitsBlockReader inblockj(hitBlocks[bj]->fPtr, hitBlocks[bj]->fSize);
				if (not inblockj.BufferSizeOk()) continue;
				
				for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
				{
					if (inblocki[i].fHit == inblockj[j])
					{
						found = true;
						
						// Since we found the corresponding cluster,
						// check if the detector element IDs are the same.
						bool deOk = CheckDetElemIds(
								*clusterBlocks[bi], bi, i, inblocki[i],
								*hitBlocks[bj], bj, j, inblockj[j]
							);
						if (not deOk)
						{
							MarkBlock(blocks, blockOk, blockCount, clusterBlocks[bi]);
							MarkBlock(blocks, blockOk, blockCount, hitBlocks[bj]);
						}
						
						break;
					}
				}
			}
			
			// If the hit was not found then it should be nil.
			if (not found and (inblocki[i].fHit != AliHLTMUONConstants::NilRecHitStruct()))
			{
				HLTError("Problem found with cluster data block %d,"
					" fDataType = '%s', fPtr = %p and fSize = %u bytes."
					" Problem with entry %d in block: The cluster hit"
					" coordinate {x = %f, y = %f, z = %f} does not exist"
					" in any reconstructed hit data block.",
					bi,
					DataType2Text(clusterBlocks[bi]->fDataType).c_str(),
					clusterBlocks[bi]->fPtr,
					clusterBlocks[bi]->fSize,
					i,
					inblocki[i].fHit.fX,
					inblocki[i].fHit.fY,
					inblocki[i].fHit.fZ
				);
				MarkBlock(blocks, blockOk, blockCount, clusterBlocks[bi]);
			}
			
			// Check that the fNchannels value is correct.
			AliHLTUInt32_t count = 0;
			for (AliHLTUInt32_t bj = 0; bj < channelBlocksCount and not found; bj++)
			{
				AliHLTMUONChannelsBlockReader inblockj(channelBlocks[bj]->fPtr, channelBlocks[bj]->fSize);
				if (not inblockj.BufferSizeOk()) continue;
				
				for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
				{
					if (inblocki[i].fId == inblockj[j].fClusterId)
					{
						count++;
					}
				}
			}
			
			if (inblocki[i].fNchannels != count)
			{
				HLTWarning("Problem found with cluster data block %d,"
					" fDataType = '%s', fPtr = %p and fSize = %u bytes."
					" Problem with entry %d in block: The number of"
					" channels in the cluster is reported as %d, but"
					" only %d channel structures were found.",
					bi,
					DataType2Text(clusterBlocks[bi]->fDataType).c_str(),
					clusterBlocks[bi]->fPtr,
					clusterBlocks[bi]->fSize,
					i,
					inblocki[i].fNchannels,
					count
				);
				MarkBlock(blocks, blockOk, blockCount, clusterBlocks[bi]);
			}
		}
		
		// Check if all the cluster structures are unique up to the identifier
		// and have unique identifiers.
		for (AliHLTUInt32_t bj = bi+1; bj < clusterBlocksCount; bj++)
		{
			AliHLTMUONClustersBlockReader inblockj(clusterBlocks[bj]->fPtr, clusterBlocks[bj]->fSize);
			if (not inblockj.BufferSizeOk()) continue;
			
			for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
			for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
			{
				if (inblocki[i].fId == inblockj[j].fId)
				{
					HLTError("Problem found with cluster"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and cluster data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The cluster %d in block %d and entry"
						" %d in block %d have the same identifier, but they"
						" should be unique.",
						bi,
						DataType2Text(clusterBlocks[bi]->fDataType).c_str(),
						clusterBlocks[bi]->fPtr,
						clusterBlocks[bi]->fSize,
						bj,
						DataType2Text(clusterBlocks[bj]->fDataType).c_str(),
						clusterBlocks[bj]->fPtr,
						clusterBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, clusterBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, clusterBlocks[bj]);
				}
				
				AliHLTMUONClusterStruct a = inblocki[i];
				AliHLTMUONClusterStruct b = inblockj[j];
				a.fId = b.fId = -1;
				if (a == b)
				{
					HLTError("Problem found with cluster"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and cluster data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The cluster %d in block %d and entry"
						" %d in block %d have the same data."
						" The data may have been duplicated.",
						bi,
						DataType2Text(clusterBlocks[bi]->fDataType).c_str(),
						clusterBlocks[bi]->fPtr,
						clusterBlocks[bi]->fSize,
						bj,
						DataType2Text(clusterBlocks[bj]->fDataType).c_str(),
						clusterBlocks[bj]->fPtr,
						clusterBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, clusterBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, clusterBlocks[bj]);
				}
			}
		}
	}
	
	for (AliHLTUInt32_t bi = 0; bi < channelBlocksCount; bi++)
	{
		AliHLTMUONChannelsBlockReader inblocki(channelBlocks[bi]->fPtr, channelBlocks[bi]->fSize);
		if (not inblocki.BufferSizeOk()) continue;
		
		// Check if all the cluster IDs in the channel structures exist.
		for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
		{
			bool found = false;
			
			for (AliHLTUInt32_t bj = 0; bj < clusterBlocksCount and not found; bj++)
			{
				AliHLTMUONClustersBlockReader inblockj(clusterBlocks[bj]->fPtr, clusterBlocks[bj]->fSize);
				if (not inblockj.BufferSizeOk()) continue;
				
				for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
				{
					if (inblocki[i].fClusterId == inblockj[j].fId)
					{
						found = true;
						break;
					}
				}
			}
			
			if (not found)
			{
				HLTError("Problem found with channel"
					" data block %d, fDataType = '%s', fPtr = %p and"
					" fSize = %u bytes."
					" Problem with entry %d in block: The cluster"
					" identifier %d does not exist in any cluster"
					" data block.",
					bi,
					DataType2Text(channelBlocks[bi]->fDataType).c_str(),
					channelBlocks[bi]->fPtr,
					channelBlocks[bi]->fSize,
					i, inblocki[i].fClusterId
				);
				MarkBlock(blocks, blockOk, blockCount, channelBlocks[bi]);
			}
		}
		
		// Check if all the channel structures are unique up to the cluster ID.
		for (AliHLTUInt32_t bj = bi+1; bj < channelBlocksCount; bj++)
		{
			AliHLTMUONChannelsBlockReader inblockj(channelBlocks[bj]->fPtr, channelBlocks[bj]->fSize);
			if (not inblockj.BufferSizeOk()) continue;
			
			for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
			for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
			{
				AliHLTMUONChannelStruct a = inblocki[i];
				AliHLTMUONChannelStruct b = inblockj[j];
				a.fClusterId = b.fClusterId = -1;
				if (a == b)
				{
					HLTError("Problem found with channel"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and channel data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The channel %d in block %d and entry"
						" %d in block %d have the same data."
						" The data may have been duplicated.",
						bi,
						DataType2Text(channelBlocks[bi]->fDataType).c_str(),
						channelBlocks[bi]->fPtr,
						channelBlocks[bi]->fSize,
						bj,
						DataType2Text(channelBlocks[bj]->fDataType).c_str(),
						channelBlocks[bj]->fPtr,
						channelBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, channelBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, channelBlocks[bj]);
				}
			}
		}
	}
	
	// Will need the total number of tracks later for comparison to trigger scalars.
	AliHLTUInt32_t totalTrackCount = 0;
	
	for (AliHLTUInt32_t bi = 0; bi < mansoTrackBlocksCount; bi++)
	{
		AliHLTMUONMansoTracksBlockReader inblocki(mansoTrackBlocks[bi]->fPtr, mansoTrackBlocks[bi]->fSize);
		if (not inblocki.BufferSizeOk()) continue;
		
		totalTrackCount += inblocki.Nentries();
		
		// Check if all the trigger record IDs in the Manso track structures exist.
		for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
		{
			bool found = false;
			
			for (AliHLTUInt32_t bj = 0; bj < trigRecBlocksCount and not found; bj++)
			{
				AliHLTMUONTriggerRecordsBlockReader inblockj(trigRecBlocks[bj]->fPtr, trigRecBlocks[bj]->fSize);
				if (not inblockj.BufferSizeOk()) continue;
				
				for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
				{
					if (inblocki[i].fTrigRec == inblockj[j].fId)
					{
						// At this point we can check if the momentum
						// is compatible with the trigger record.
						if (not AreMomentaCompatible(
								inblocki[i].fPx, inblocki[i].fPy, inblocki[i].fPz,
								inblockj[j].fPx, inblockj[j].fPy, inblockj[j].fPz
							)
						   )
						{
							HLTWarning("Problem found with Manso track"
								" data block %d, fDataType = '%s', fPtr = %p and"
								" fSize = %u bytes."
								" Problem with Manso track %d in block: The momentum"
								" vector of the track p = {%f, %f, %f} GeV/c is not"
								" compatible with the momentum vector of the trigger"
								" record with p = {%f, %f, %f} GeV/c.",
								bi,
								DataType2Text(mansoTrackBlocks[bi]->fDataType).c_str(),
								mansoTrackBlocks[bi]->fPtr,
								mansoTrackBlocks[bi]->fSize,
								i, inblocki[i].fPx, inblocki[i].fPy, inblocki[i].fPz,
								inblockj[j].fPx, inblockj[j].fPy, inblockj[j].fPz
							);
							MarkBlock(blocks, blockOk, blockCount, mansoTrackBlocks[bi]);
						}
						
						found = true;
						break;
					}
				}
			}
			
			if (not found)
			{
				HLTError("Problem found with Manso track"
					" data block %d, fDataType = '%s', fPtr = %p and"
					" fSize = %u bytes."
					" Problem with Manso track %d in block: The trigger"
					" record identifier %d does not exist in any trigger"
					" record data block.",
					bi,
					DataType2Text(mansoTrackBlocks[bi]->fDataType).c_str(),
					mansoTrackBlocks[bi]->fPtr,
					mansoTrackBlocks[bi]->fSize,
					i, inblocki[i].fTrigRec
				);
				MarkBlock(blocks, blockOk, blockCount, mansoTrackBlocks[bi]);
			}
		}
		
		// Check if all the hits in the Manso track structures exist.
		for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
		{
			AliHLTMUONParticleSign sign;
			bool hitset[4];
			AliHLTMUONUtils::UnpackMansoTrackFlags(inblocki[i].fFlags, sign, hitset);
			
			for (AliHLTUInt32_t n = 0; n < 4; n++)
			{
				if (not hitset[n]) continue;
				bool found = false;
				
				for (AliHLTUInt32_t bj = 0; bj < hitBlocksCount and not found; bj++)
				{
					AliHLTMUONRecHitsBlockReader inblockj(hitBlocks[bj]->fPtr, hitBlocks[bj]->fSize);
					if (not inblockj.BufferSizeOk()) continue;
					
					for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
					{
						if (inblocki[i].fHit[n] == inblockj[j])
						{
							found = true;
							break;
						}
					}
				}
				
				if (not found)
				{
					HLTError("Problem found with Manso track"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes."
						" Problem with Manso track %d in block: The hit"
						" for chamber %d does not exist in any"
						" reconstructed hits data block.",
						bi,
						DataType2Text(mansoTrackBlocks[bi]->fDataType).c_str(),
						mansoTrackBlocks[bi]->fPtr,
						mansoTrackBlocks[bi]->fSize,
						i, n+6+1
					);
					MarkBlock(blocks, blockOk, blockCount, mansoTrackBlocks[bi]);
				}
			}
		}
		
		// Check if all the Manso track structures are unique up to the ID and
		// have unique identifiers.
		for (AliHLTUInt32_t bj = bi+1; bj < mansoTrackBlocksCount; bj++)
		{
			AliHLTMUONMansoTracksBlockReader inblockj(mansoTrackBlocks[bj]->fPtr, mansoTrackBlocks[bj]->fSize);
			if (not inblockj.BufferSizeOk()) continue;
			
			for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
			for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
			{
				if (inblocki[i].fId == inblockj[j].fId)
				{
					HLTError("Problem found with Manso track"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and Manso track data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The Manso track %d in block %d and entry"
						" %d in block %d have the same identifier, but they"
						" should be unique.",
						bi,
						DataType2Text(mansoTrackBlocks[bi]->fDataType).c_str(),
						mansoTrackBlocks[bi]->fPtr,
						mansoTrackBlocks[bi]->fSize,
						bj,
						DataType2Text(mansoTrackBlocks[bj]->fDataType).c_str(),
						mansoTrackBlocks[bj]->fPtr,
						mansoTrackBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, mansoTrackBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, mansoTrackBlocks[bj]);
				}
				
				AliHLTMUONMansoTrackStruct a = inblocki[i];
				AliHLTMUONMansoTrackStruct b = inblockj[j];
				a.fId = b.fId = -1;
				if (a == b)
				{
					HLTError("Problem found with Manso track"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and Manso track data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The Manso track %d in block %d and entry"
						" %d in block %d have the same data."
						" The data may have been duplicated.",
						bi,
						DataType2Text(mansoTrackBlocks[bi]->fDataType).c_str(),
						mansoTrackBlocks[bi]->fPtr,
						mansoTrackBlocks[bi]->fSize,
						bj,
						DataType2Text(mansoTrackBlocks[bj]->fDataType).c_str(),
						mansoTrackBlocks[bj]->fPtr,
						mansoTrackBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, mansoTrackBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, mansoTrackBlocks[bj]);
				}
			}
		}
	}
	
	for (AliHLTUInt32_t bi = 0; bi < mansoCandidateBlocksCount; bi++)
	{
		AliHLTMUONMansoCandidatesBlockReader inblocki(mansoCandidateBlocks[bi]->fPtr, mansoCandidateBlocks[bi]->fSize);
		if (not inblocki.BufferSizeOk()) continue;
		
		// Check if all the trigger record IDs in the Manso track candidate structures exist.
		for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
		{
			bool found = false;
			
			for (AliHLTUInt32_t bj = 0; bj < trigRecBlocksCount and not found; bj++)
			{
				AliHLTMUONTriggerRecordsBlockReader inblockj(trigRecBlocks[bj]->fPtr, trigRecBlocks[bj]->fSize);
				if (not inblockj.BufferSizeOk()) continue;
				
				for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
				{
					if (inblocki[i].fTrack.fTrigRec == inblockj[j].fId)
					{
						// At this point we can check if the momentum
						// is compatible with the trigger record.
						if (not AreMomentaCompatible(
								inblocki[i].fTrack.fPx, inblocki[i].fTrack.fPy, inblocki[i].fTrack.fPz,
								inblockj[j].fPx, inblockj[j].fPy, inblockj[j].fPz
							)
						   )
						{
							HLTWarning("Problem found with Manso track candidate"
								" data block %d, fDataType = '%s', fPtr = %p and"
								" fSize = %u bytes."
								" Problem with track candidate %d in block: The momentum"
								" vector of the candidate p = {%f, %f, %f} GeV/c is not"
								" compatible with the momentum vector of the trigger"
								" record with p = {%f, %f, %f} GeV/c.",
								bi,
								DataType2Text(mansoTrackBlocks[bi]->fDataType).c_str(),
								mansoTrackBlocks[bi]->fPtr,
								mansoTrackBlocks[bi]->fSize,
								i,
								inblocki[i].fTrack.fPx, inblocki[i].fTrack.fPy, inblocki[i].fTrack.fPz,
								inblockj[j].fPx, inblockj[j].fPy, inblockj[j].fPz
							);
							MarkBlock(blocks, blockOk, blockCount, mansoTrackBlocks[bi]);
						}
						
						found = true;
						break;
					}
				}
			}
			
			if (not found)
			{
				HLTError("Problem found with Manso track candidate"
					" data block %d, fDataType = '%s', fPtr = %p and"
					" fSize = %u bytes."
					" Problem with track candidate %d in block: The trigger"
					" record identifier %d does not exist in any trigger"
					" record data block.",
					bi,
					DataType2Text(mansoCandidateBlocks[bi]->fDataType).c_str(),
					mansoCandidateBlocks[bi]->fPtr,
					mansoCandidateBlocks[bi]->fSize,
					i, inblocki[i].fTrack.fTrigRec
				);
				MarkBlock(blocks, blockOk, blockCount, mansoCandidateBlocks[bi]);
			}
		}
		
		// Check if all the hits in the Manso track candidate structures exist.
		for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
		{
			AliHLTMUONParticleSign sign;
			bool hitset[4];
			AliHLTMUONUtils::UnpackMansoTrackFlags(inblocki[i].fTrack.fFlags, sign, hitset);
			
			for (AliHLTUInt32_t n = 0; n < 4; n++)
			{
				if (not hitset[n]) continue;
				bool found = false;
				
				for (AliHLTUInt32_t bj = 0; bj < hitBlocksCount and not found; bj++)
				{
					AliHLTMUONRecHitsBlockReader inblockj(hitBlocks[bj]->fPtr, hitBlocks[bj]->fSize);
					if (not inblockj.BufferSizeOk()) continue;
					
					for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
					{
						if (inblocki[i].fTrack.fHit[n] == inblockj[j])
						{
							found = true;
							break;
						}
					}
				}
				
				if (not found)
				{
					HLTError("Problem found with Manso track candidate"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes."
						" Problem with track candidate %d in block: The hit"
						" for chamber %d does not exist in any"
						" reconstructed hits data block.",
						bi,
						DataType2Text(mansoTrackBlocks[bi]->fDataType).c_str(),
						mansoTrackBlocks[bi]->fPtr,
						mansoTrackBlocks[bi]->fSize,
						i, n+6+1
					);
					MarkBlock(blocks, blockOk, blockCount, mansoTrackBlocks[bi]);
				}
			}
		}
		
		// Check if all the Manso track candidate structures are unique up to the
		// track ID and have unique identifiers.
		for (AliHLTUInt32_t bj = bi+1; bj < mansoCandidateBlocksCount; bj++)
		{
			AliHLTMUONMansoCandidatesBlockReader inblockj(mansoCandidateBlocks[bj]->fPtr, mansoCandidateBlocks[bj]->fSize);
			if (not inblockj.BufferSizeOk()) continue;
			
			for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
			for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
			{
				if (inblocki[i].fTrack.fId == inblockj[j].fTrack.fId)
				{
					HLTError("Problem found with Manso track candidate"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and Manso track candidate data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The track candidate %d in block %d and entry"
						" %d in block %d have the same identifier, but they"
						" should be unique.",
						bi,
						DataType2Text(mansoCandidateBlocks[bi]->fDataType).c_str(),
						mansoCandidateBlocks[bi]->fPtr,
						mansoCandidateBlocks[bi]->fSize,
						bj,
						DataType2Text(mansoCandidateBlocks[bj]->fDataType).c_str(),
						mansoCandidateBlocks[bj]->fPtr,
						mansoCandidateBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, mansoCandidateBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, mansoCandidateBlocks[bj]);
				}
				
				AliHLTMUONMansoCandidateStruct a = inblocki[i];
				AliHLTMUONMansoCandidateStruct b = inblockj[j];
				a.fTrack.fId = b.fTrack.fId = -1;
				if (a == b)
				{
					HLTError("Problem found with Manso track candidate"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and Manso track candidate data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The track candidate %d in block %d and entry"
						" %d in block %d have the same data."
						" The data may have been duplicated.",
						bi,
						DataType2Text(mansoCandidateBlocks[bi]->fDataType).c_str(),
						mansoCandidateBlocks[bi]->fPtr,
						mansoCandidateBlocks[bi]->fSize,
						bj,
						DataType2Text(mansoCandidateBlocks[bj]->fDataType).c_str(),
						mansoCandidateBlocks[bj]->fPtr,
						mansoCandidateBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, mansoCandidateBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, mansoCandidateBlocks[bj]);
				}
			}
		}
	}
	
	for (AliHLTUInt32_t bi = 0; bi < singleDecisionBlocksCount; bi++)
	{
		AliHLTMUONSinglesDecisionBlockReader inblocki(singleDecisionBlocks[bi]->fPtr, singleDecisionBlocks[bi]->fSize);
		if (not inblocki.BufferSizeOk()) continue;
		
		// Check that the scalars are within reasonable limits.
		const AliHLTMUONSinglesDecisionBlockStruct& hdr = inblocki.BlockHeader();
		const AliHLTComponentBlockData* block = singleDecisionBlocks[bi];
		if (IsScalarTooLarge(block, bi, "single track", "fNlowPt", hdr.fNlowPt, totalTrackCount) or
		    IsScalarTooLarge(block, bi, "single track", "fNhighPt", hdr.fNhighPt, totalTrackCount) or
		    IsScalarALargerThanB(block, bi, "single track", "fNhighPt", hdr.fNhighPt, "fNlowPt", hdr.fNlowPt)
		   )
		{
			MarkBlock(blocks, blockOk, blockCount, block);
		}
		
		// Check if all the Manso track IDs in the trigger decision structures exist.
		for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
		{
			bool found = false;
			
			for (AliHLTUInt32_t bj = 0; bj < mansoTrackBlocksCount and not found; bj++)
			{
				AliHLTMUONMansoTracksBlockReader inblockj(mansoTrackBlocks[bj]->fPtr, mansoTrackBlocks[bj]->fSize);
				if (not inblockj.BufferSizeOk()) continue;
				
				for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
				{
					if (inblocki[i].fTrackId == inblockj[j].fId)
					{
						found = true;
						break;
					}
				}
			}
			
			if (not found)
			{
				HLTError("Problem found with single track trigger decision"
					" data block %d, fDataType = '%s', fPtr = %p and"
					" fSize = %u bytes."
					" Problem with decision %d in block: The track"
					" identifier %d does not exist in any Manso tracks"
					" data block.",
					bi,
					DataType2Text(singleDecisionBlocks[bi]->fDataType).c_str(),
					singleDecisionBlocks[bi]->fPtr,
					singleDecisionBlocks[bi]->fSize,
					i, inblocki[i].fTrackId
				);
				MarkBlock(blocks, blockOk, blockCount, singleDecisionBlocks[bi]);
			}
		}
		
		// Check if all the trigger decision structures are unique up to the ID and
		// have unique Manso track identifiers.
		for (AliHLTUInt32_t bj = bi+1; bj < singleDecisionBlocksCount; bj++)
		{
			AliHLTMUONSinglesDecisionBlockReader inblockj(singleDecisionBlocks[bj]->fPtr, singleDecisionBlocks[bj]->fSize);
			if (not inblockj.BufferSizeOk()) continue;
			
			for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
			for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
			{
				if (inblocki[i].fTrackId == inblockj[j].fTrackId)
				{
					HLTError("Problem found with single track trigger decision"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and single track trigger decision"
						" data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The trigger decision %d in block %d and entry"
						" %d in block %d have the same Manso track identifier,"
						" but they should be unique.",
						bi,
						DataType2Text(singleDecisionBlocks[bi]->fDataType).c_str(),
						singleDecisionBlocks[bi]->fPtr,
						singleDecisionBlocks[bi]->fSize,
						bj,
						DataType2Text(singleDecisionBlocks[bj]->fDataType).c_str(),
						singleDecisionBlocks[bj]->fPtr,
						singleDecisionBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, singleDecisionBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, singleDecisionBlocks[bj]);
				}
				
				AliHLTMUONTrackDecisionStruct a = inblocki[i];
				AliHLTMUONTrackDecisionStruct b = inblockj[j];
				a.fTrackId = b.fTrackId = -1;
				if (a == b)
				{
					HLTError("Problem found with single track trigger decision"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and single track trigger decision"
						" data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The trigger decision %d in block %d and entry"
						" %d in block %d have the same data."
						" The data may have been duplicated.",
						bi,
						DataType2Text(singleDecisionBlocks[bi]->fDataType).c_str(),
						singleDecisionBlocks[bi]->fPtr,
						singleDecisionBlocks[bi]->fSize,
						bj,
						DataType2Text(singleDecisionBlocks[bj]->fDataType).c_str(),
						singleDecisionBlocks[bj]->fPtr,
						singleDecisionBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, singleDecisionBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, singleDecisionBlocks[bj]);
				}
			}
		}
	}
	
	for (AliHLTUInt32_t bi = 0; bi < pairDecisionBlocksCount; bi++)
	{
		AliHLTMUONPairsDecisionBlockReader inblocki(pairDecisionBlocks[bi]->fPtr, pairDecisionBlocks[bi]->fSize);
		if (not inblocki.BufferSizeOk()) continue;
		
		AliHLTUInt32_t maxPairs = totalTrackCount * (totalTrackCount-1) / 2;
		const AliHLTMUONPairsDecisionBlockStruct& hdr = inblocki.BlockHeader();
		const AliHLTComponentBlockData* block = pairDecisionBlocks[bi];
		if (IsScalarTooLarge(block, bi, "track pair", "fNunlikeAnyPt", hdr.fNunlikeAnyPt, maxPairs) or
		    IsScalarTooLarge(block, bi, "track pair", "fNunlikeLowPt", hdr.fNunlikeLowPt, maxPairs) or
		    IsScalarTooLarge(block, bi, "track pair", "fNunlikeHighPt", hdr.fNunlikeHighPt, maxPairs) or
		    IsScalarTooLarge(block, bi, "track pair", "fNlikeAnyPt", hdr.fNlikeAnyPt, maxPairs) or
		    IsScalarTooLarge(block, bi, "track pair", "fNlikeLowPt", hdr.fNlikeLowPt, maxPairs) or
		    IsScalarTooLarge(block, bi, "track pair", "fNlikeHighPt", hdr.fNlikeHighPt, maxPairs) or
		    IsScalarTooLarge(block, bi, "track pair", "fNmassAny", hdr.fNmassAny, maxPairs) or
		    IsScalarTooLarge(block, bi, "track pair", "fNmassLow", hdr.fNmassLow, maxPairs) or
		    IsScalarTooLarge(block, bi, "track pair", "fNmassHigh", hdr.fNmassHigh, maxPairs) or
		    IsScalarALargerThanB(block, bi, "track pair", "fNunlikeHighPt", hdr.fNunlikeHighPt, "fNunlikeLowPt", hdr.fNunlikeLowPt) or
		    IsScalarALargerThanB(block, bi, "track pair", "fNunlikeLowPt", hdr.fNunlikeLowPt, "fNunlikeAnyPt", hdr.fNunlikeAnyPt) or
		    IsScalarALargerThanB(block, bi, "track pair", "fNlikeHighPt", hdr.fNlikeHighPt, "fNlikeLowPt", hdr.fNlikeLowPt) or
		    IsScalarALargerThanB(block, bi, "track pair", "fNlikeLowPt", hdr.fNlikeLowPt, "fNlikeAnyPt", hdr.fNlikeAnyPt) or
		    IsScalarALargerThanB(block, bi, "track pair", "fNmassHigh", hdr.fNmassHigh, "fNmassLow", hdr.fNmassLow) or
		    IsScalarALargerThanB(block, bi, "track pair", "fNmassLow", hdr.fNmassLow, "fNmassAny", hdr.fNmassAny)
		   )
		{
			MarkBlock(blocks, blockOk, blockCount, block);
		}
		
		// Check if all the Manso track IDs in the trigger decision structures exist.
		for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
		{
			bool found = false;
			
			for (AliHLTUInt32_t bj = 0; bj < mansoTrackBlocksCount and not found; bj++)
			{
				AliHLTMUONMansoTracksBlockReader inblockj(mansoTrackBlocks[bj]->fPtr, mansoTrackBlocks[bj]->fSize);
				if (not inblockj.BufferSizeOk()) continue;
				
				for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
				{
					if (inblocki[i].fTrackAId == inblockj[j].fId)
					{
						found = true;
						break;
					}
				}
			}
			
			if (not found)
			{
				HLTError("Problem found with track pair trigger decision"
					" data block %d, fDataType = '%s', fPtr = %p and"
					" fSize = %u bytes."
					" Problem with decision %d in block: The track"
					" identifier %d does not exist in any Manso tracks"
					" data block.",
					bi,
					DataType2Text(pairDecisionBlocks[bi]->fDataType).c_str(),
					pairDecisionBlocks[bi]->fPtr,
					pairDecisionBlocks[bi]->fSize,
					i, inblocki[i].fTrackAId
				);
				MarkBlock(blocks, blockOk, blockCount, pairDecisionBlocks[bi]);
			}
			
			found = false;
			
			for (AliHLTUInt32_t bj = 0; bj < mansoTrackBlocksCount and not found; bj++)
			{
				AliHLTMUONMansoTracksBlockReader inblockj(mansoTrackBlocks[bj]->fPtr, mansoTrackBlocks[bj]->fSize);
				if (not inblockj.BufferSizeOk()) continue;
				
				for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
				{
					if (inblocki[i].fTrackBId == inblockj[j].fId)
					{
						found = true;
						break;
					}
				}
			}
			
			if (not found)
			{
				HLTError("Problem found with track pair trigger decision"
					" data block %d, fDataType = '%s', fPtr = %p and"
					" fSize = %u bytes."
					" Problem with decision %d in block: The track"
					" identifier %d does not exist in any Manso tracks"
					" data block.",
					bi,
					DataType2Text(pairDecisionBlocks[bi]->fDataType).c_str(),
					pairDecisionBlocks[bi]->fPtr,
					pairDecisionBlocks[bi]->fSize,
					i, inblocki[i].fTrackBId
				);
				MarkBlock(blocks, blockOk, blockCount, pairDecisionBlocks[bi]);
			}
		}
		
		// Check if all the trigger decision structures are unique up to the ID and
		// have unique Manso track identifier pairs.
		for (AliHLTUInt32_t bj = bi+1; bj < pairDecisionBlocksCount; bj++)
		{
			AliHLTMUONPairsDecisionBlockReader inblockj(pairDecisionBlocks[bj]->fPtr, pairDecisionBlocks[bj]->fSize);
			if (not inblockj.BufferSizeOk()) continue;
			
			for (AliHLTUInt32_t i = 0; i < inblocki.Nentries(); i++)
			for (AliHLTUInt32_t j = 0; j < inblockj.Nentries(); j++)
			{
				if (inblocki[i].fTrackAId == inblockj[j].fTrackAId and
				    inblocki[i].fTrackBId == inblockj[j].fTrackBId
				   )
				{
					HLTError("Problem found with track pair trigger decision"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and track pair trigger decision"
						" data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The trigger decision %d in block %d and entry"
						" %d in block %d have the same Manso track identifier pair,"
						" but the pair should be unique.",
						bi,
						DataType2Text(pairDecisionBlocks[bi]->fDataType).c_str(),
						pairDecisionBlocks[bi]->fPtr,
						pairDecisionBlocks[bi]->fSize,
						bj,
						DataType2Text(pairDecisionBlocks[bj]->fDataType).c_str(),
						pairDecisionBlocks[bj]->fPtr,
						pairDecisionBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, pairDecisionBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, pairDecisionBlocks[bj]);
				}
				
				AliHLTMUONPairDecisionStruct a = inblocki[i];
				AliHLTMUONPairDecisionStruct b = inblockj[j];
				a.fTrackAId = a.fTrackBId = b.fTrackAId = b.fTrackBId = -1;
				if (a == b)
				{
					HLTError("Problem found with track pair trigger decision"
						" data block %d, fDataType = '%s', fPtr = %p and"
						" fSize = %u bytes, and track pair trigger decision"
						" data block %d,"
						" fDataType = '%s', fPtr = %p and fSize = %u bytes."
						" Problem: The trigger decision %d in block %d and entry"
						" %d in block %d have the same data."
						" The data may have been duplicated.",
						bi,
						DataType2Text(pairDecisionBlocks[bi]->fDataType).c_str(),
						pairDecisionBlocks[bi]->fPtr,
						pairDecisionBlocks[bi]->fSize,
						bj,
						DataType2Text(pairDecisionBlocks[bj]->fDataType).c_str(),
						pairDecisionBlocks[bj]->fPtr,
						pairDecisionBlocks[bj]->fSize,
						bi, i,
						bj, j
					);
					MarkBlock(blocks, blockOk, blockCount, pairDecisionBlocks[bi]);
					MarkBlock(blocks, blockOk, blockCount, pairDecisionBlocks[bj]);
				}
			}
		}
	}
}

