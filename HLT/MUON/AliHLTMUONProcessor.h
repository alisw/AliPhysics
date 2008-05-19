#ifndef ALIHLTMUONPROCESSOR_H
#define ALIHLTMUONPROCESSOR_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///
/// @file   AliHLTMUONProcessor.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   2007-12-12
/// @brief  Declaration of a common processor component abstract interface for dHLT components.
///

#include "AliHLTProcessor.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONUtils.h"

/**
 * @class AliHLTMUONProcessor
 * This component class is an abstract base class for dHLT components.
 * Some common methods useful to all dHLT specific components are implemented
 * by this class.
 */
class AliHLTMUONProcessor : public AliHLTProcessor
{
public:
	/// Default constructor.
	AliHLTMUONProcessor() : AliHLTProcessor() {}
	/// Default destructor.
	virtual ~AliHLTMUONProcessor() {}

protected:

	/**
	 * Method to check the block structure and log appropriate error messages.
	 * If a problem is found with the data block then an appropriate HLT error
	 * message is logged and the method returns false.
	 * \param block  The lightweight block reader whose data block should be checked.
	 * \param name  A string containing a descriptive name of the data block
	 *          type. This name is used in the logged error messages.
	 * \returns  true if the structure of the block looks OK and false otherwise.
	 * \note  The BlockType should be a class deriving from AliHLTMUONDataBlockReader.
	 */
	template <class BlockType>
	bool BlockStructureOk(const BlockType& block, const char* name) const;
	
	/// Checks the structure of a trigger records data block.
	bool BlockStructureOk(const AliHLTMUONTriggerRecordsBlockReader& block) const
	{
		return BlockStructureOk(block, "trigger records");
	}
	
	/// Checks the structure of a trigger records debug information data block.
	bool BlockStructureOk(const AliHLTMUONTrigRecsDebugBlockReader& block) const
	{
		return BlockStructureOk(block, "trigger records debug information");
	}

	/// Checks the structure of a reconstructed hits data block.
	bool BlockStructureOk(const AliHLTMUONRecHitsBlockReader& block) const
	{
		return BlockStructureOk(block, "reconstructed hits");
	}
	
	/// Checks the structure of a clusters data block.
	bool BlockStructureOk(const AliHLTMUONClustersBlockReader& block) const
	{
		return BlockStructureOk(block, "clusters");
	}
	
	/// Checks the structure of a ADC channels data block.
	bool BlockStructureOk(const AliHLTMUONChannelsBlockReader& block) const
	{
		return BlockStructureOk(block, "channels");
	}

	/// Checks the structure of a Manso tracks data block.
	bool BlockStructureOk(const AliHLTMUONMansoTracksBlockReader& block) const
	{
		return BlockStructureOk(block, "Manso tracks");
	}
	
	/// Checks the structure of a Manso track candidates data block.
	bool BlockStructureOk(const AliHLTMUONMansoCandidatesBlockReader& block) const
	{
		return BlockStructureOk(block, "Manso track candidates");
	}

	/// Checks the structure of a single track trigger decision data block.
	bool BlockStructureOk(const AliHLTMUONSinglesDecisionBlockReader& block) const
	{
		return BlockStructureOk(block, "singles decision");
	}

	/// Checks the structure of a track pairs trigger decision data block.
	bool BlockStructureOk(const AliHLTMUONPairsDecisionBlockReader& block) const
	{
		return BlockStructureOk(block, "pairs decision");
	}
	
private:

	// Do not allow copying of this class.
	AliHLTMUONProcessor(const AliHLTMUONProcessor& /*obj*/);
	AliHLTMUONProcessor& operator = (const AliHLTMUONProcessor& /*obj*/);
	
	ClassDef(AliHLTMUONProcessor, 0)  // Abstract base class for dHLT specific components.
};

//______________________________________________________________________________

template <class BlockType>
bool AliHLTMUONProcessor::BlockStructureOk(const BlockType& block, const char* name) const
{
	/// Performs basic checks to see if the input data block structure is OK,
	/// that it is not corrupt, too short etc...
	
	if (not block.BufferSizeOk())
	{
		size_t headerSize = sizeof(typename BlockType::HeaderType);
		if (block.BufferSize() < headerSize)
		{
			HLTError("Received a %s data block with a size of %d bytes,"
				" which is smaller than the minimum valid header size of %d bytes."
				" The block must be corrupt.",
				name, block.BufferSize(), headerSize
			);
			return false;
		}
		
		size_t expectedWidth = sizeof(typename BlockType::ElementType);
		if (block.CommonBlockHeader().fRecordWidth != expectedWidth)
		{
			HLTError("Received a %s data block with a record"
				" width of %d bytes, but the expected value is %d bytes."
				" The block might be corrupt.",
				name,
				block.CommonBlockHeader().fRecordWidth,
				expectedWidth
			);
			return false;
		}
		
		HLTError("Received a %s data block with a size of %d bytes,"
			" but the block header claims the block should be %d bytes."
			" The block might be corrupt.",
			name, block.BufferSize(), block.BytesUsed()
		);
		return false;
	}
	
	AliHLTMUONUtils::WhyNotValid reason;
	if (not AliHLTMUONUtils::HeaderOk(block.BlockHeader(), &reason))
	{
		HLTError("Received a %s data block which might be corrupt. %s",
			name, AliHLTMUONUtils::FailureReasonToMessage(reason)
		);
		return false;
	}
	
	return true;
}

#endif // ALIHLTMUONPROCESSOR_H
