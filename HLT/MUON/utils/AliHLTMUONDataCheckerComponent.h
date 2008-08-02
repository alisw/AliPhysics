#ifndef ALIHLTMUONDATACHECKERCOMPONENT_H
#define ALIHLTMUONDATACHECKERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id: AliHLTMUONDataCheckerComponent.h 26179 2008-05-29 22:27:27Z aszostak $ */

///
/// @file   AliHLTMUONDataCheckerComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   27 May 2008
/// @brief  Declaration of a component to check and validate the integrity of dHLT data.
///

#include "AliHLTMUONProcessor.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif

class AliMpDDLStore;
extern "C" struct AliHLTMUONRecHitStruct;
extern "C" struct AliHLTMUONClusterStruct;
extern "C" struct AliHLTMUONTriggerRecordStruct;
extern "C" struct AliHLTMUONTrigRecInfoStruct;
extern "C" struct AliHLTMUONMansoTrackStruct;

/**
 * @class AliHLTMUONDataCheckerComponent
 * This component is used to check and validate the integrity of dHLT internal
 * raw data blocks. If there are any problems then an appropriate error message
 * is logged.
 * The component should be rather used only for debugging and testing since the
 * validation procedure can be slow.<br>
 *
 * Component ID: \b MUONDataChecker <br>
 * Library: \b libAliHLTMUON.so  <br>
 *
 * Optional arguments:<br>
 * \li -ignoretype <br>
 *       Indicates if the data type of the raw data blocks as given by the HLT
 *       framework in the DoEvent method should be ignored.
 *       (default behaviour is not to ignore the type)<br>
 * \li -ignorespec <br>
 *       Indicates if the data specification of the raw data blocks as given by
 *       the HLT framework in the DoEvent method should be ignored.
 *       (default behaviour is not to ignore the specification)<br>
 * \li -dontforward <br>
 *       If set given then the data blocks are not forwarded as output from this
 *       component.
 *       (default behaviour is not to forward all data blocks)<br>
 * \li -filter <br>
 *       If specified then all the data blocks for which a problem was found are
 *       forwarded to the output and everything else is suppressed.
 *       (default behaviour is not to filter for bad data blocks)<br>
 * \li -no_global_check <br>
 *       If specified then no global checks between the data blocks are performed,
 *       but rather only the per data block checks.
 *       (default behaviour is to perform all per block and global checks)<br>
 * \li -warn_on_unexpected_block <br>
 *       Indicates if a warning message should be generated when this component
 *       receives an unknown data block type.
 *       (default behaviour is not to just log a debugging message)<br>
 * \li -return_error <br>
 *       Indicates if error codes should be returned from the DoEvent method which
 *       would tell the framework that processing of the event failed. Otherwise
 *       errors are just logged but the data is considered to be processed successfully.
 *       (default behaviour is not to return errors)<br>
 *
 * @ingroup alihlt_dimuon_component
 */
class AliHLTMUONDataCheckerComponent : public AliHLTMUONProcessor
{
public:
	AliHLTMUONDataCheckerComponent();
	virtual ~AliHLTMUONDataCheckerComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the component registration process.

	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();
	
protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component.

	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit();
	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr, 
			AliHLTUInt32_t& size,
			AliHLTComponentBlockDataList& outputBlocks
		);
	
	using AliHLTProcessor::DoEvent;
	
private:

	// Do not allow copying of this class.
	AliHLTMUONDataCheckerComponent(const AliHLTMUONDataCheckerComponent& /*obj*/);
	AliHLTMUONDataCheckerComponent& operator = (const AliHLTMUONDataCheckerComponent& /*obj*/);
	
	bool IsSpecificationValid(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber,
			const char* name
		) const;
	
	bool IsFromTrackerOnly(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber,
			const char* name
		) const;
	
	bool IsFromTriggerOnly(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber,
			const char* name
		) const;
	
	bool IsMomentumVectorOk(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber,
			const char* name,
			AliHLTUInt32_t entryNumber,
			AliHLTFloat32_t px,
			AliHLTFloat32_t py,
			AliHLTFloat32_t pz
		) const;
	
	bool AreMomentumCalcParamsOk(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber,
			const char* name,
			AliHLTUInt32_t entryNumber,
			AliHLTFloat32_t zmiddle,
			AliHLTFloat32_t bl
		) const;
	
	bool IsHitCoordinateOk(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber,
			const char* name,
			AliHLTUInt32_t entryNumber,
			const AliHLTMUONRecHitStruct& hit,
			AliHLTInt32_t minChamber,
			AliHLTInt32_t maxChamber,
			AliHLTInt32_t expectedChamber,
			bool ddl[22]
		) const;
	
	bool IsMansoTrackOk(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber,
			const char* name,
			AliHLTUInt32_t entryNumber,
			const AliHLTMUONMansoTrackStruct& track,
			bool ddl[22]
		) const;
	
	bool CheckDetElemIds(
			const AliHLTComponentBlockData& infoBlock,
			AliHLTUInt32_t infoBlockNumber,
			AliHLTUInt32_t infoEntryNumber,
			const AliHLTMUONTrigRecInfoStruct& info,
			const AliHLTComponentBlockData& trBlock,
			AliHLTUInt32_t trBlockNumber,
			AliHLTUInt32_t trEntryNumber,
			const AliHLTMUONTriggerRecordStruct& tr
		) const;
	
	bool CheckDetElemIds(
			const AliHLTComponentBlockData& clusterBlock,
			AliHLTUInt32_t clusterBlockNumber,
			AliHLTUInt32_t clusterEntryNumber,
			const AliHLTMUONClusterStruct& cluster,
			const AliHLTComponentBlockData& hitBlock,
			AliHLTUInt32_t hitBlockNumber,
			AliHLTUInt32_t hitEntryNumber,
			const AliHLTMUONRecHitStruct& hit
		) const;
	
	bool CheckRawDataBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	bool CheckTriggerRecordsBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	bool CheckTrigRecsDebugBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	bool CheckRecHitsBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	bool CheckClustersBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	bool CheckChannelsBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	bool CheckMansoTracksBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	bool CheckMansoCandidatesBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	bool CheckSinglesDecisionBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	bool CheckPairsDecisionBlock(
			const AliHLTComponentBlockData& block,
			AliHLTUInt32_t blockNumber
		) const;
	
	template <class DataBlock>
	bool CheckBlockHeaderOnly(
			const AliHLTComponentBlockData& descriptor,
			AliHLTUInt32_t blockNumber,
			DataBlock& block,
			const char* name
		) const;
	
	template <class DataBlock>
	bool CheckBlockIntegrity(
			const AliHLTComponentBlockData& descriptor,
			AliHLTUInt32_t blockNumber,
			DataBlock& block,
			const char* name
		) const;
	
	bool AreMomentaCompatible(
			AliHLTFloat32_t px1,
			AliHLTFloat32_t py1,
			AliHLTFloat32_t pz1,
			AliHLTFloat32_t px2,
			AliHLTFloat32_t py2,
			AliHLTFloat32_t pz2
		) const;
	
	bool IsScalarTooLarge(
			const AliHLTComponentBlockData* block,
			AliHLTUInt32_t blockNumber,
			const char* blockTypeName,
			const char* scalarName,
			AliHLTUInt32_t scalarValue,
			AliHLTUInt32_t totalTrackCount
		) const;
	
	bool IsScalarALargerThanB(
			const AliHLTComponentBlockData* block,
			AliHLTUInt32_t blockNumber,
			const char* blockTypeName,
			const char* scalarAName,
			AliHLTUInt32_t scalarAValue,
			const char* scalarBName,
			AliHLTUInt32_t scalarBValue
		) const;
	
	void MarkBlock(
			const AliHLTComponentBlockData* blocks,
			bool* blockOk,
			AliHLTUInt32_t blockCount,
			const AliHLTComponentBlockData* blockToMark
		) const;
	
	void MakeGlobalChecks(
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
		) const;
	
	bool fIgnoreType; ///< Flag indicating if we should ignore the data block type as given in DoEvent by the framework.
	bool fIgnoreSpec; ///< Flag indicating if we should ignore the data block specification as given in DoEvent by the framework.
	bool fDontForward; ///< Flag indicating if we should not forward the input data blocks as output.
	bool fFilterBadBlocks; ///< Flag indicating if we should pass through only bad blocks to output.
	bool fNoGlobalChecks;  ///< Flag indicating if we should perform global data consistancy checks between all the data blocks.
	bool fWarnForUnexpecedBlock; ///< Flag indicating if we should log a warning if we got a block of an unexpected type.
	bool fReturnError; ///< Flag indicating if we should return error codes from DoEvent.

	ClassDef(AliHLTMUONDataCheckerComponent, 0)  // dHLT raw internal data block checking component.
};

//______________________________________________________________________________

template <class DataBlock>
bool AliHLTMUONDataCheckerComponent::CheckBlockHeaderOnly(
		const AliHLTComponentBlockData& descriptor,
		AliHLTUInt32_t blockNumber,
		DataBlock& block,
		const char* name
	) const
{
	/// Method for checking the integrity of dHLT raw internal data blocks,
	/// which only need their headers checked.
	/// \returns  true if the block structure is OK and false otherwise.

	// 'false' set so that we do not check the common block header since it
	// will be done in IntegrityOk.
	if (not BlockStructureOk(block, false)) return false;
	
	AliHLTUInt32_t count = 256;
	AliHLTMUONUtils::WhyNotValid reason[256];
	if (not AliHLTMUONUtils::IntegrityOk(block.BlockHeader(), &reason[0], count))
	{
		for (AliHLTUInt32_t i = 0; i < count; i++)
		{
			HLTError("Problem found with data block %d, fDataType = '%s',"
				 " fPtr = %p and fSize = %u bytes."
				 " Assuming this is a %s data block. Problem: %s",
				blockNumber,
				DataType2Text(descriptor.fDataType).c_str(),
				descriptor.fPtr,
				descriptor.fSize,
				name,
				AliHLTMUONUtils::FailureReasonToMessage(reason[i])
			);
		}
	}
	
	return true;
}


template <class DataBlock>
bool AliHLTMUONDataCheckerComponent::CheckBlockIntegrity(
		const AliHLTComponentBlockData& descriptor,
		AliHLTUInt32_t blockNumber,
		DataBlock& block,
		const char* name
	) const
{
	/// Method for checking the integrity of dHLT raw internal data blocks.
	/// \returns  true if the block structure is OK and false otherwise.
	
	// 'false' set so that we do not check the common block header since it
	// will be done in IntegrityOk.
	if (not BlockStructureOk(block, false)) return false;
	
	AliHLTUInt32_t count = 256;
	AliHLTMUONUtils::WhyNotValid reason[256];
	AliHLTUInt32_t recordNum[256];
	if (not AliHLTMUONUtils::IntegrityOk(block.BlockHeader(), &reason[0], &recordNum[0], count))
	{
		for (AliHLTUInt32_t i = 0; i < count; i++)
		{
			if (AliHLTMUONUtils::RecordNumberWasSet(reason[i]))
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem with entry %d in data block: %s",
					blockNumber,
					DataType2Text(descriptor.fDataType).c_str(),
					descriptor.fPtr,
					descriptor.fSize,
					name,
					recordNum[i],
					AliHLTMUONUtils::FailureReasonToMessage(reason[i])
				);
			}
			else
			{
				HLTError("Problem found with data block %d, fDataType = '%s',"
					 " fPtr = %p and fSize = %u bytes."
					 " Assuming this is a %s data block."
					 " Problem: %s",
					blockNumber,
					DataType2Text(descriptor.fDataType).c_str(),
					descriptor.fPtr,
					descriptor.fSize,
					name,
					AliHLTMUONUtils::FailureReasonToMessage(reason[i])
				);
			}
		}
	}
	
	return true;
}

#endif // ALIHLTMUONDATACHECKERCOMPONENT_H

