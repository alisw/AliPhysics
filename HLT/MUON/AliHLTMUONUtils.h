#ifndef ALIHLTMUONUTILS_H
#define ALIHLTMUONUTILS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONUtils.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   17 May 2007
/// @brief  Class containing various dimuon HLT utility routines and macros.
///

#include "TObject.h"
#include "AliHLTMUONDataTypes.h"
#include <ostream>

// Forward declare structures.
extern "C" {
struct AliHLTMUONTriggerRecordStruct;
struct AliHLTMUONTriggerRecordsBlockStruct;
struct AliHLTMUONTrigRecInfoStruct;
struct AliHLTMUONTrigRecsDebugBlockStruct;
struct AliHLTMUONRecHitStruct;
struct AliHLTMUONRecHitsBlockStruct;
struct AliHLTMUONClusterStruct;
struct AliHLTMUONClustersBlockStruct;
struct AliHLTMUONChannelStruct;
struct AliHLTMUONChannelsBlockStruct;
struct AliHLTMUONMansoTrackStruct;
struct AliHLTMUONMansoTracksBlockStruct;
struct AliHLTMUONMansoCandidateStruct;
struct AliHLTMUONMansoCandidatesBlockStruct;
struct AliHLTMUONTrackDecisionStruct;
struct AliHLTMUONSinglesDecisionBlockStruct;
struct AliHLTMUONPairDecisionStruct;
struct AliHLTMUONPairsDecisionBlockStruct;
} // extern "C"

/**
 * AliHLTMUONUtils contains arbitrary utility methods to be used in various
 * parts of the dimuon HLT system.
 * These include methods to perform basic sanity checks on the integrity of
 * data blocks.
 */
class AliHLTMUONUtils
{
public:
	/**
	 * This packs the given parameters into the bits of a word appropriate
	 * for AliHLTMUONTriggerRecordStruct::fFlags.
	 * @param sign    The particle sign.
	 * @param hitset  Flags to indicate if the corresponding fHits[i] elements
	 *                was set/filled.
	 * @return  Returns the 32 bit packed word.
	 */
	static AliHLTUInt32_t PackTriggerRecordFlags(
			AliHLTMUONParticleSign sign, const bool hitset[4]
		);

	/**
	 * This unpacks the AliHLTMUONTriggerRecordStruct::fFlags bits into
	 * its component fields.
	 * @param flags  The flags from an AliHLTMUONTriggerRecordStruct structure.
	 * @param sign    Sets this to the particle sign.
	 * @param hitset  Sets the array elements to indicate if the corresponding
	 *                fHits[i] element was set/filled.
	 */
	static void UnpackTriggerRecordFlags(
			AliHLTUInt32_t flags, // [in]
			AliHLTMUONParticleSign& sign, // [out]
			bool hitset[4] // [out]
		);
	
	/**
	 * This packs the given parameters into the bits of a word appropriate
	 * for AliHLTMUONRecHitStruct::fFlags.
	 * @param chamber    The chamber number in the range [0..13].
	 * @param detElemId  Detector element ID number.
	 * @return  Returns the 32 bit packed word.
	 */
	static AliHLTUInt32_t PackRecHitFlags(
			AliHLTUInt8_t chamber, AliHLTUInt16_t detElemId
		);

	/**
	 * This unpacks the AliHLTMUONRecHitStruct::fFlags bits into
	 * its component fields.
	 * [in]  @param flags  The flags from an AliHLTMUONRecHitStruct structure.
	 * [out] @param chamber    Sets the chamber number in the range [0..13].
	 * [out] @param detElemId  Sets the detector element ID number.
	 */
	static void UnpackRecHitFlags(
			AliHLTUInt32_t flags, // [in]
			AliHLTUInt8_t& chamber, // [out]
			AliHLTUInt16_t& detElemId // [out]
		);

	/**
	 * This packs the given parameters into the bits of a word appropriate
	 * for AliHLTMUONMansoTrackStruct::fFlags.
	 * @param sign    The particle sign.
	 * @param hitset  Flags to indicate if the corresponding fHits[i] elements
	 *                was set/filled.
	 * @return  Returns the 32 bit packed word.
	 */
	static AliHLTUInt32_t PackMansoTrackFlags(
			AliHLTMUONParticleSign sign, const bool hitset[4]
		)
	{
		return PackTriggerRecordFlags(sign, hitset);
	}

	/**
	 * This unpacks the AliHLTMUONMansoTrackStruct::fFlags bits into
	 * its component fields.
	 * @param flags  The flags from an AliHLTMUONMansoTrackStruct structure.
	 * @param sign    Sets this to the particle sign.
	 * @param hitset  Sets the array elements to indicate if the corresponding
	 *                fHits[i] element was set/filled.
	 */
	static void UnpackMansoTrackFlags(
			AliHLTUInt32_t flags, // [in]
			AliHLTMUONParticleSign& sign, // [out]
			bool hitset[4] // [out]
		)
	{
		UnpackTriggerRecordFlags(flags, sign, hitset);
	}
	
	/**
	 * This packs the given parameters into the bits of a word appropriate
	 * for AliHLTMUONTrackDecisionStruct::fTriggerBits.
	 * @param highPt  Has the track passed the high pt cut.
	 * @param lowPt   Has the track passed the low pt cut.
	 * @return  Returns the 32 bit packed word.
	 */
	static AliHLTUInt32_t PackTrackDecisionBits(bool highPt, bool lowPt);
	
	/**
	 * This unpacks the AliHLTMUONTrackDecisionStruct::fTriggerBits bits into
	 * its component fields.
	 * @param bits  The trigger bits from an AliHLTMUONTrackDecisionStruct
	 *              structure.
	 * @param highPt Sets this to the value of the high pt cut bit.
	 * @param lowPt  Sets this to the value of the low pt cut bit.
	 */
	static void UnpackTrackDecisionBits(
			AliHLTUInt32_t bits, // [in]
			bool& highPt, // [out]
			bool& lowPt // [out]
		);
	
	/**
	 * This packs the given parameters into the bits of a word appropriate
	 * for AliHLTMUONPairDecisionStruct::fTriggerBits.
	 *
	 * @param highMass Has the track pair passed the high invariant mass cut.
	 * @param lowMass  Has the track pair passed the low invariant mass cut.
	 * @param unlike   Does the track pair have unlike signs.
	 * @param highPtCount The number of tracks that passed the high pt cut
	 *                    in the pair.
	 * @param lowPtCount  The number of tracks that passed the low pt cut
	 *                    in the pair.
	 * @return  Returns the 32 bit packed word.
	 *
	 * Note: Must have highPtCount + lowPtCount <= 2 and unlike == true if
	 * highMass or lowMass is true.
	 */
	static AliHLTUInt32_t PackPairDecisionBits(
			bool highMass, bool lowMass, bool unlike,
			AliHLTUInt8_t highPtCount, AliHLTUInt8_t lowPtCount
		);
	
	/**
	 * This unpacks the AliHLTMUONPairDecisionStruct::fTriggerBits bits into
	 * its component fields.
	 * @param bits  The trigger bits from an AliHLTMUONPairDecisionStruct
	 *              structure.
	 * @param highMass Sets this to the value of the high invariant mass cut bit.
	 * @param lowMass  Sets this to the value of the low invariant mass cut bit.
	 * @param unlike   Sets this if the pair is unlike sign.
	 * @param highPtCount Sets this to the high pt count bits.
	 * @param lowPtCount  Sets this to the low pt count bits.
	 */
	static void UnpackPairDecisionBits(
			AliHLTUInt32_t bits, // [in]
			bool& highMass, // [out]
			bool& lowMass, // [out]
			bool& unlike, // [out]
			AliHLTUInt8_t& highPtCount, // [out]
			AliHLTUInt8_t& lowPtCount // [out]
		);
	
	/**
	 * This packs the given parameters into the 32bit Pub/Sub specification
	 * word in the data block descriptor.
	 *
	 * @param ddl  The list of DDLs forming part of the readout. ddl[0]
	 *             indicates DDL number 2560, ddl[1] is for DDL 2561 and so
	 *             on up to ddl[19]. ddl[20] and ddl[21] will be for the
	 *             trigger DDLs 2816 and 2817 respectively.
	 * @return  Returns the 32 bit packed specification word.
	 */
	static AliHLTUInt32_t PackSpecBits(
			const bool ddl[22]
		);
	
	/**
	 * This unpacks the AliHLTMUONPairDecisionStruct::fTriggerBits bits into
	 * its component fields.
	 * @param bits  The Pub/Sub specification word from a data block descriptor.
	 * @param ddl  The output list of DDLs forming part of the readout. ddl[0]
	 *             indicates DDL number 2560, ddl[1] is for DDL 2561 and so
	 *             on up to ddl[19]. ddl[20] and ddl[21] will be for the
	 *             trigger DDLs 2816 and 2817 respectively.
	 */
	static void UnpackSpecBits(
			AliHLTUInt32_t bits, // [in]
			bool ddl[22] // [out]
		);

	/**
	 * This method converts the DDL number for the muon spectrometer in the
	 * range [0..21] to the equipment ID number.
	 * @param ddlNo  The DDL number in the range [0..21].
	 * @return  Returns the equipment ID number or -1 if ddlNo was invalid.
	 */
	static AliHLTInt32_t DDLNumberToEquipId(AliHLTInt32_t ddlNo);
	
	/**
	 * This method converts the equipment ID number for a muon spectrometer
	 * DDL to the DDL number in the range [0..21].
	 * @param id  The equipment ID of the DDL.
	 * @return  Returns the DDL number in the range [0..21] or -1 if the
	 *          equipment ID was invalid.
	 */
	static AliHLTInt32_t EquipIdToDDLNumber(AliHLTInt32_t id);
	
	/**
	 * This method converts a 32 bit data block specification for a MUON-HLT
	 * data block into its corresponding DDL equipment ID number.
	 * It is assumed that the specification is for a data block comming from
	 * a single DDL source. If more than one DDL contributed to the data block
	 * then -1 is returned.
	 * @param spec  The 32 bit specification for a data block.
	 * @return  Returns the equipment ID corresponding to the specification
	 *          or -1 if the specification was invalid.
	 */
	static AliHLTInt32_t SpecToEquipId(AliHLTUInt32_t spec);
	
	/**
	 * This method converts a equipment ID number for a DDL into its corresponding
	 * 32 bit data block specification for the MUON-HLT.
	 * @param id  The equipment ID number of the DDL.
	 * @return  Returns the 32 bit data block specification or 0x0 if the
	 *          equipment ID was invalid.
	 */
	static AliHLTUInt32_t EquipIdToSpec(AliHLTInt32_t id);
	
	/**
	 * This method converts a 32 bit data block specification for a MUON-HLT
	 * data block into its corresponding DDL number in the range [0..21].
	 * It is assumed that the specification is for a data block comming from
	 * a single DDL source. If more than one DDL contributed to the data block
	 * then -1 is returned.
	 * @param spec  The 32 bit specification for a data block.
	 * @return  Returns the corresponding DDL number for the specification
	 *          or -1 if the specification was invalid.
	 */
	static AliHLTInt32_t SpecToDDLNumber(AliHLTUInt32_t spec);
	
	/**
	 * This method converts a DDL number in the range [0..21] into its
	 * corresponding 32 bit data block specification for the MUON-HLT.
	 * @param ddlNo  The equipment ID number of the DDL.
	 * @return  Returns the 32 bit data block specification or 0x0 if the
	 *          DDL number was invalid (out of range).
	 */
	static AliHLTUInt32_t DDLNumberToSpec(AliHLTInt32_t ddlNo);

	/**
	 * Returns true if the given specification was for a single trigger DDL.
	 */
	static bool IsTriggerDDL(AliHLTUInt32_t spec)
	{
		AliHLTInt32_t ddl = SpecToDDLNumber(spec);
		return (20 <= ddl and ddl <= 21);
	}

	/**
	 * Returns true if the given specification was for a single tracker DDL.
	 */
	static bool IsTrackerDDL(AliHLTUInt32_t spec)
	{
		AliHLTInt32_t ddl = SpecToDDLNumber(spec);
		return (0 <= ddl and ddl <= 19);
	}

	/**
	 * Returns true if the given specification is in principle valid.
	 * It checks if the bits that should be zero are indeed zero.
	 */
	static bool IsSpecValid(AliHLTUInt32_t spec)
	{
		AliHLTUInt32_t mask = ~((1 << 22) - 1);  // First 22 bits indicate DDL number.
		return (spec & mask) == 0x0;
	}

	/**
	 * Returns true if the data specification indicates the data block contains
	 * information generated from a trigger DDL or data fragments thereof.
	 */
	static bool ContainsDataFromTrigger(AliHLTUInt32_t spec)
	{
		AliHLTUInt32_t mask = ((1 << 22) - 1) & ~((1 << 20) - 1);
		return (spec & mask) != 0x0;
	}

	/**
	 * Returns true if the data specification indicates the data block contains
	 * information generated from a tracker DDL or data fragments thereof.
	 */
	static bool ContainsDataFromTracker(AliHLTUInt32_t spec)
	{
		AliHLTUInt32_t mask = ((1 << 20) - 1);
		return (spec & mask) != 0x0;
	}
	
	/**
	 * Parses the string containing the type name of a dHLT data block and
	 * returns the corresponding AliHLTMUONDataBlockType value.
	 * \param  type  The string containing the type name.
	 * \returns  The data block type or kUnknownDataBlock if the type name
	 *      is invalid.
	 */
	static AliHLTMUONDataBlockType ParseCommandLineTypeString(const char* type);
	
	/**
	 * Converts a type ID to a type string to be used for the dHLT FilePublisher
	 * component configuration parameters for example.
	 */
	static const char* DataBlockTypeToString(AliHLTMUONDataBlockType type);
	
	/**
	 * These codes indicate the reason why a data block failed its
	 * validity check.
	 */
	enum WhyNotValid
	{
		kNoReason,   ///< There was no reason for failure.
		kHeaderContainsWrongType,  ///< The common header contains an incorrect type ID.
		kHeaderContainsWrongRecordWidth,  ///< The common header contains an incorrect data record width.
		kInvalidIdValue,  ///< The structure identifier does not have a valid value.
		kInvalidTriggerIdValue,  ///< The trigger structure identifier does not have a valid value.
		kInvalidTrackIdValue,  ///< The track structure identifier does not have a valid value.
		kReservedBitsNotZero,  ///< Reserved bits have not been set to zero.
		kParticleSignBitsNotValid,  ///< The particle sign bits are not a valid value.
		kHitNotMarkedAsNil,  ///< A hit was marked as not found, but the corresponding hit structure was not set to nil.
		kInvalidDetElementNumber,  ///< An invalid detector element ID was found.
		kInvalidChamberNumber,  ///< An invalid chamber number was found.
		kHitIsNil,  ///< The hit cannot be set to a nil value.
		kInvalidChannelCount,  ///< The number of channels indicated is zero or outside the valid range.
		kInvalidTotalCharge, ///< The total charge does not have a valid value.
		kInvalidBusPatchId,  ///< The bus patch ID is outside the valid range.
		kInvalidManuId,  ///< The MANU ID is outside the valid range.
		kInvalidChannelAddress,  ///< The MANU channel address is outside the valid range.
		kInvalidSignal,  ///< The ADC signal value is outside the valid range.
		kDataWordDifferent, ///< The raw data word is different from the unpacked values.
		kChiSquareInvalid,  ///< The chi squared value must be a positive value or -1 indicating a fitting error.
		kMomentumVectorNotZero, ///< The chi sqaured value is set to -1, but momentum vector not zero.
		kRoiRadiusInvalid, ///< The region of interest radius is invalid.
		kHitNotWithinRoi, ///< A tracks hit is not within the corresponding region of interest.
		kPtValueNotValid,  ///< The pT value is not positive nor -1 indicating an invalid value.
		kPairTrackIdsAreIdentical,  ///< The track IDs of the track pair are identical.
		kMassValueNotValid,  ///< The invariant mass value is not positive nor -1 indicating an invalid value.
		kLowPtCountInvalid,  ///< The low pT trigger count is greater than 2, which is invalid.
		kHighPtCountInvalid,  ///< The high pT trigger count is greater than 2, which is invalid.
		kFoundDuplicateIDs,  ///< Found duplicate identifiers, but they should all be unique.
		kFoundDuplicateHits,  ///< Found duplicate hits.
		kFoundDuplicateTriggers  ///< Found duplicate trigger decisions.
	};
	
	/**
	 * This method converts the WhyNotValid enumeration to a string representation.
	 */
	static const char* FailureReasonToString(WhyNotValid reason);
	
	/**
	 * This method returns a string containing a user readable message explaining
	 * the reason for failure described by the WhyNotValid enumeration.
	 */
	static const char* FailureReasonToMessage(WhyNotValid reason);

	/**
	 * Method used to check if the header information corresponds to the
	 * supposed type of the trigger records data block.
	 * This method will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason code.
	 * [in]  \param block  The data block to check.
	 * [out] \param reason  If this is not NULL, then the variable pointed to
	 *      by this pointer will be filled with the reason code describing why
	 *      the header is not valid, if and only if a problem is found with
	 *      the data.
	 * \returns  true if there is no problem with the header and false otherwise.
	 */
	static bool HeaderOk(const AliHLTMUONTriggerRecordsBlockStruct& block, WhyNotValid* reason = NULL)
	{
		AliHLTUInt32_t count = 1;
		return HeaderOk(block, reason, count);
	}
	
	/**
	 * Method used to check if the header information corresponds to the
	 * supposed type of the trigger debug information data block.
	 * This method will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason code.
	 * [in]  \param block  The data block to check.
	 * [out] \param reason  If this is not NULL, then the variable pointed to
	 *      by this pointer will be filled with the reason code describing why
	 *      the header is not valid, if and only if a problem is found with
	 *      the data.
	 * \returns  true if there is no problem with the header and false otherwise.
	 */
	static bool HeaderOk(const AliHLTMUONTrigRecsDebugBlockStruct& block, WhyNotValid* reason = NULL)
	{
		AliHLTUInt32_t count = 1;
		return HeaderOk(block, reason, count);
	}
	
	/**
	 * Method used to check if the header information corresponds to the
	 * supposed type of the reconstructed hits data block.
	 * This method will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason code.
	 * [in]  \param block  The data block to check.
	 * [out] \param reason  If this is not NULL, then the variable pointed to
	 *      by this pointer will be filled with the reason code describing why
	 *      the header is not valid, if and only if a problem is found with
	 *      the data.
	 * \returns  true if there is no problem with the header and false otherwise.
	 */
	static bool HeaderOk(const AliHLTMUONRecHitsBlockStruct& block, WhyNotValid* reason = NULL)
	{
		AliHLTUInt32_t count = 1;
		return HeaderOk(block, reason, count);
	}
	
	/**
	 * Method used to check if the header information corresponds to the
	 * supposed type of the clusters data block.
	 * This method will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason code.
	 * [in]  \param block  The data block to check.
	 * [out] \param reason  If this is not NULL, then the variable pointed to
	 *      by this pointer will be filled with the reason code describing why
	 *      the header is not valid, if and only if a problem is found with
	 *      the data.
	 * \returns  true if there is no problem with the header and false otherwise.
	 */
	static bool HeaderOk(const AliHLTMUONClustersBlockStruct& block, WhyNotValid* reason = NULL)
	{
		AliHLTUInt32_t count = 1;
		return HeaderOk(block, reason, count);
	}
	
	/**
	 * Method used to check if the header information corresponds to the
	 * supposed type of the channels data block.
	 * This method will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason code.
	 * [in]  \param block  The data block to check.
	 * [out] \param reason  If this is not NULL, then the variable pointed to
	 *      by this pointer will be filled with the reason code describing why
	 *      the header is not valid, if and only if a problem is found with
	 *      the data.
	 * \returns  true if there is no problem with the header and false otherwise.
	 */
	static bool HeaderOk(const AliHLTMUONChannelsBlockStruct& block, WhyNotValid* reason = NULL)
	{
		AliHLTUInt32_t count = 1;
		return HeaderOk(block, reason, count);
	}
	
	/**
	 * Method used to check if the header information corresponds to the
	 * supposed type of the Manso tracks data block.
	 * This method will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason code.
	 * [in]  \param block  The data block to check.
	 * [out] \param reason  If this is not NULL, then the variable pointed to
	 *      by this pointer will be filled with the reason code describing why
	 *      the header is not valid, if and only if a problem is found with
	 *      the data.
	 * \returns  true if there is no problem with the header and false otherwise.
	 */
	static bool HeaderOk(const AliHLTMUONMansoTracksBlockStruct& block, WhyNotValid* reason = NULL)
	{
		AliHLTUInt32_t count = 1;
		return HeaderOk(block, reason, count);
	}
	
	/**
	 * Method used to check if the header information corresponds to the
	 * supposed type of the Manso candidates data block.
	 * This method will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason code.
	 * [in]  \param block  The data block to check.
	 * [out] \param reason  If this is not NULL, then the variable pointed to
	 *      by this pointer will be filled with the reason code describing why
	 *      the header is not valid, if and only if a problem is found with
	 *      the data.
	 * \returns  true if there is no problem with the header and false otherwise.
	 */
	static bool HeaderOk(const AliHLTMUONMansoCandidatesBlockStruct& block, WhyNotValid* reason = NULL)
	{
		AliHLTUInt32_t count = 1;
		return HeaderOk(block, reason, count);
	}
	
	/**
	 * Method used to check if the header information corresponds to the
	 * supposed type of the single tracks dHLT trigger decision data block.
	 * This method will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason code.
	 * [in]  \param block  The data block to check.
	 * [out] \param reason  If this is not NULL, then the variable pointed to
	 *      by this pointer will be filled with the reason code describing why
	 *      the header is not valid, if and only if a problem is found with
	 *      the data.
	 * \returns  true if there is no problem with the header and false otherwise.
	 */
	static bool HeaderOk(const AliHLTMUONSinglesDecisionBlockStruct& block, WhyNotValid* reason = NULL)
	{
		AliHLTUInt32_t count = 1;
		return HeaderOk(block, reason, count);
	}
	
	/**
	 * Method used to check if the header information corresponds to the
	 * supposed type of the track pairs dHLT trigger decision data block.
	 * This method will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason code.
	 * [in]  \param block  The data block to check.
	 * [out] \param reason  If this is not NULL, then the variable pointed to
	 *      by this pointer will be filled with the reason code describing why
	 *      the header is not valid, if and only if a problem is found with
	 *      the data.
	 * \returns  true if there is no problem with the header and false otherwise.
	 */
	static bool HeaderOk(const AliHLTMUONPairsDecisionBlockStruct& block, WhyNotValid* reason = NULL)
	{
		AliHLTUInt32_t count = 1;
		return HeaderOk(block, reason, count);
	}

	/**
	 * Methods used to check if the header information corresponds to the
	 * supposed type of the data block.
	 * If the 'reason' parameter should point to an array which will store
	 * the reason codes indicating the problems with the data block.
	 * The 'reasonCount' parameter should initialy contain the number of
	 * elements that can be stored in reason. When the method exits it will
	 * store the number of elements in the 'reason' array actually filled.
	 */
	static bool HeaderOk(
			const AliHLTMUONTriggerRecordsBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool HeaderOk(
			const AliHLTMUONTrigRecsDebugBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool HeaderOk(
			const AliHLTMUONRecHitsBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool HeaderOk(
			const AliHLTMUONClustersBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool HeaderOk(
			const AliHLTMUONChannelsBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool HeaderOk(
			const AliHLTMUONMansoTracksBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool HeaderOk(
			const AliHLTMUONMansoCandidatesBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool HeaderOk(
			const AliHLTMUONSinglesDecisionBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool HeaderOk(
			const AliHLTMUONPairsDecisionBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * trigger record structure is OK and returns true in that case.
	 * [in] \param tr  The trigger record structure to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the structure is not valid, if and
	 *      only if a problem is found with the data.
	 * \returns  true if there is no problem with the structure and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONTriggerRecordStruct& tr,
			WhyNotValid* reason = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(tr, reason, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * dHLT raw internal data block is OK and returns true in that case.
	 * [in] \param block  The trigger record data block to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the data block is not valid, if and
	 *      only if a problem is found with the data.
	 * [out] \param recordNum  If this is not NULL, then it will be filled with
	 *      the number of the trigger record that had a problem. This value will
	 *      only contain a valid value if the method RecordNumberWasSet(*reason)
	 *      returns true. Thus, 'reason' must be set.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONTriggerRecordsBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(block, reason, recordNum, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * trigger record debug information structure is OK and returns true in that case.
	 * [in] \param trigInfo  The trigger record debug information structure to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the structure is not valid, if and
	 *      only if a problem is found with the data.
	 * \returns  true if there is no problem with the structure and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONTrigRecInfoStruct& trigInfo,
			WhyNotValid* reason = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(trigInfo, reason, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * dHLT raw internal data block is OK and returns true in that case.
	 * [in] \param block  The trigger record debugging information data block to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the data block is not valid, if and
	 *      only if a problem is found with the data.
	 * [out] \param recordNum  If this is not NULL, then it will be filled with
	 *      the number of the trigger record debug information structure that had
	 *      a problem. This value will only contain a valid value if the method
	 *      RecordNumberWasSet(*reason) returns true. Thus, 'reason' must be set.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONTrigRecsDebugBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(block, reason, recordNum, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * reconstructed hit structure is OK and returns true in that case.
	 * [in] \param hit  The reconstructed hit structure to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the structure is not valid, if and
	 *      only if a problem is found with the data.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONRecHitStruct& hit,
			WhyNotValid* reason = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(hit, reason, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * dHLT raw internal data block is OK and returns true in that case.
	 * [in] \param block  The reconstructed hits data block to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the data block is not valid, if and
	 *      only if a problem is found with the data.
	 * [out] \param recordNum  If this is not NULL, then it will be filled with
	 *      the number of the cluster structure that had a problem. This value
	 *      will only contain a valid value if the method
	 *      RecordNumberWasSet(*reason) returns true. Thus, 'reason' must be set.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONRecHitsBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(block, reason, recordNum, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * cluster data structure is OK and returns true in that case.
	 * [in] \param cluster  The cluster structure to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the structure is not valid, if and
	 *      only if a problem is found with the data.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONClusterStruct& cluster,
			WhyNotValid* reason = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(cluster, reason, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * dHLT raw internal data block is OK and returns true in that case.
	 * [in] \param block  The clusters data block to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the data block is not valid, if and
	 *      only if a problem is found with the data.
	 * [out] \param recordNum  If this is not NULL, then it will be filled with
	 *      the number of the cluster structure that had a problem. This value
	 *      will only contain a valid value if the method
	 *      RecordNumberWasSet(*reason) returns true. Thus, 'reason' must be set.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONClustersBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(block, reason, recordNum, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * channel data structure is OK and returns true in that case.
	 * [in] \param cluster  The channel structure to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the structure is not valid, if and
	 *      only if a problem is found with the data.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONChannelStruct& channel,
			WhyNotValid* reason = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(channel, reason, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * dHLT raw internal data block is OK and returns true in that case.
	 * [in] \param block  The ADC channels data block to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the data block is not valid, if and
	 *      only if a problem is found with the data.
	 * [out] \param recordNum  If this is not NULL, then it will be filled with
	 *      the number of the channel structure that had a problem. This value
	 *      will only contain a valid value if the method
	 *      RecordNumberWasSet(*reason) returns true. Thus, 'reason' must be set.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONChannelsBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(block, reason, recordNum, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * Manso track structure is OK and returns true in that case.
	 * [in] \param track  The Manso track structure to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the structure is not valid, if and
	 *      only if a problem is found with the data.
	 * \returns  true if there is no problem with the structure and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONMansoTrackStruct& track,
			WhyNotValid* reason = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(track, reason, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * dHLT raw internal data block is OK and returns true in that case.
	 * [in] \param block  The Manso track data block to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the data block is not valid, if and
	 *      only if a problem is found with the data.
	 * [out] \param recordNum  If this is not NULL, then it will be filled with
	 *      the number of the Manso track structure that had a problem.
	 *      This value will only contain a valid value if the method
	 *      RecordNumberWasSet(*reason) returns true. Thus, 'reason' must be set.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONMansoTracksBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(block, reason, recordNum, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * Manso track candidate structure is OK and returns true in that case.
	 * [in] \param candidate  The Manso track candidate structure to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the structure is not valid, if and
	 *      only if a problem is found with the data.
	 * \returns  true if there is no problem with the structure and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONMansoCandidateStruct& candidate,
			WhyNotValid* reason = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(candidate, reason, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * dHLT raw internal data block is OK and returns true in that case.
	 * [in] \param block  The Manso track candidate data block to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the data block is not valid, if and
	 *      only if a problem is found with the data.
	 * [out] \param recordNum  If this is not NULL, then it will be filled with
	 *      the number of the Manso track candidate structure that had a problem.
	 *      This value will only contain a valid value if the method
	 *      RecordNumberWasSet(*reason) returns true. Thus, 'reason' must be set.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONMansoCandidatesBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(block, reason, recordNum, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * single track trigger decision structure is OK and returns true in that case.
	 * [in] \param decision  The trigger decision structure to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the structure is not valid, if and
	 *      only if a problem is found with the data.
	 * \returns  true if there is no problem with the structure and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONTrackDecisionStruct& decision,
			WhyNotValid* reason = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(decision, reason, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * dHLT raw internal data block is OK and returns true in that case.
	 * [in] \param block  The single track trigger decision data block to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the data block is not valid, if and
	 *      only if a problem is found with the data.
	 * [out] \param recordNum  If this is not NULL, then it will be filled with
	 *      the number of the single track trigger decision structure that had
	 *      a problem. This value will only contain a valid value if the method
	 *      RecordNumberWasSet(*reason) returns true. Thus, 'reason' must be set.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONSinglesDecisionBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(block, reason, recordNum, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * track pair trigger decision structure is OK and returns true in that case.
	 * [in] \param decision  The trigger decision structure to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the structure is not valid, if and
	 *      only if a problem is found with the data.
	 * \returns  true if there is no problem with the structure and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONPairDecisionStruct& decision,
			WhyNotValid* reason = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(decision, reason, count);
	}
	
	/**
	 * This method is used to check more extensively if the integrity of the
	 * dHLT raw internal data block is OK and returns true in that case.
	 * [in] \param block  The track pair trigger decision data block to check.
	 * [out] \param reason  If this is not NULL, then it will be filled with
	 *      the reason code describing why the data block is not valid, if and
	 *      only if a problem is found with the data.
	 * [out] \param recordNum  If this is not NULL, then it will be filled with
	 *      the number of the track pairs trigger decision structure that had
	 *      a problem. This value will only contain a valid value if the method
	 *      RecordNumberWasSet(*reason) returns true. Thus, 'reason' must be set.
	 * \returns  true if there is no problem with the data and false otherwise.
	 */
	static bool IntegrityOk(
			const AliHLTMUONPairsDecisionBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		)
	{
		AliHLTUInt32_t count = 1;
		return IntegrityOk(block, reason, recordNum, count);
	}

	/**
	 * Methods used to check more extensively if the integrity of various
	 * types of data blocks are Ok and returns true in that case.
	 * These can be slow and should generally only be used for debugging.
	 * The methods are able to return multiple reasons for the problems related
	 * to the data block under test.
	 */
	static bool IntegrityOk(
			const AliHLTMUONTriggerRecordStruct& tr,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONTriggerRecordsBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t* recordNum,
			AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONTrigRecInfoStruct& trigInfo,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONTrigRecsDebugBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t* recordNum,
			AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONRecHitStruct& hit,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONRecHitsBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t* recordNum,
			AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONClusterStruct& cluster,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONClustersBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t* recordNum,
			AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONChannelStruct& channel,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONChannelsBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t* recordNum,
			AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONMansoTrackStruct& track,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONMansoTracksBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t* recordNum,
			AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONMansoCandidateStruct& candidate,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONMansoCandidatesBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t* recordNum,
			AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONTrackDecisionStruct& decision,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONSinglesDecisionBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t* recordNum,
			AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONPairDecisionStruct& decision,
			WhyNotValid* reason, AliHLTUInt32_t& reasonCount
		);
	
	static bool IntegrityOk(
			const AliHLTMUONPairsDecisionBlockStruct& block,
			WhyNotValid* reason, AliHLTUInt32_t* recordNum,
			AliHLTUInt32_t& reasonCount
		);
	
	/**
	 * Returns true if the \em recordNum in the corresponding IntegrityOk method
	 * would have been set, if it returned false and a reason was set.
	 * This helper method makes it easy to test if the \em recordNum parameter
	 * is filled with a valid value or not.
	 */
	static bool RecordNumberWasSet(WhyNotValid reason);

private:

	// Should never have to create, copy or destroy this object.
	AliHLTMUONUtils() {}
	AliHLTMUONUtils(const AliHLTMUONUtils& obj);
	virtual ~AliHLTMUONUtils() {}
	AliHLTMUONUtils& operator = (const AliHLTMUONUtils& obj);
	
	ClassDef(AliHLTMUONUtils, 0);  // Interface for helpful dHLT utility methods.
};

//_____________________________________________________________________________

inline std::ostream& operator << (std::ostream& stream, AliHLTMUONUtils::WhyNotValid reason)
{
	/// Stream operator for the WhyNotValid enumeration for usage with
	/// std::ostream classes. Allows usages such as:
	/// AliHLTMUONUtils::WhyNotValid r; std::cout << r;
	
	stream << AliHLTMUONUtils::FailureReasonToString(reason);
	return stream;
}

//_____________________________________________________________________________

// Since c++ is missing a finally "keyword" we define one. Its usage is identical
// to a try..finally statement in Java etc.. however, since it is officialy a macro
// one must use the ( ) brackets instead of { }
// If the compiler supports __finally use it otherwise make our own.
#if defined(__BORLANDC__)
#	define finally(str) __finally{str}
#else
#	define finally(code) \
		catch(...) \
		{ \
			code \
			throw; \
		}; \
		code
#endif // __BORLANDC__

// Here we define the DebugTrace(message) macro for easy embedding of debug
// information into the dimuon HLT code. Its usage is meant to be for generating
// traces of the program which are only useful during full scale debugging.
// Log messages should use the standard HLT logging mechanisms.
// The output is only generated in programs compiled with the DEBUG directive
// defined. Here is a usage example:
//
//  // statements...
//  DebugTrace("some debug information.");
//  // statements...
//
// One can also use C++ ostream operators and manipulators like so:
//
//  // statements...
//  int x, y;
//  DebugTrace("x = " << x << " and y = 0x" << std::hex << y );
//  // statements...
//
#ifdef DEBUG
#	include <iostream>
#	define DebugTrace(message) {std::cout << message << std::endl;}
#else // DEBUG
#	define DebugTrace(message)
#endif // DEBUG


#endif // ALIHLTMUONUTILS_H
