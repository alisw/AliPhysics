#ifndef ALIHLTMUONUTILS_H
#define ALIHLTMUONUTILS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONUtils.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   
/// @brief  Class containing various dimuon HLT utility routines and macros.
///

#include "AliHLTMUONDataTypes.h"
#include <ostream>

// Forward declare structures.
extern "C" {
struct AliHLTMUONTriggerRecordStruct;
struct AliHLTMUONTriggerRecordsBlockStruct;
struct AliHLTMUONTrigRecsDebugBlockStruct;
struct AliHLTMUONTriggerChannelsBlockStruct;
struct AliHLTMUONRecHitsBlockStruct;
struct AliHLTMUONClustersBlockStruct;
struct AliHLTMUONChannelsBlockStruct;
struct AliHLTMUONMansoTrackStruct;
struct AliHLTMUONMansoTracksBlockStruct;
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
	* Parses the string containing the type name of a dHLT data block and
	* returns the corresponding AliHLTMUONDataBlockType value.
	* \param  type  The string containing the type name.
	* \returns  The data block type or kUnknownDataBlock if the type name
	*      is invalid.
	*/
	static AliHLTMUONDataBlockType ParseCommandLineTypeString(const char* type);

	/**
	 * These codes indicate the reason why a data block failed its
	 * validity check.
	 */
	enum WhyNotValid
	{
		kNoReason,   ///< There was no reason for failure.
		kHeaderContainsWrongType,  ///< The common header contains an incorrect type ID.
		kHeaderContainsWrongRecordWidth,  ///< The common header contains an incorrect data record width.
		kReservedBitsNotZero,  ///< Reserved bits have not been set to zero.
		kParticleSignBitsNotValid,  ///< The particle sign bits are not a valid value.
		kHitNotMarkedAsNil,  ///< A hit was marked as not found, but the corresponding hit structure was not set to nil.
		kFoundDuplicateIDs,  ///< Found duplicate identifiers, but they should all be unique.
		kPtValueNotValid,  ///< The pT value is not positive nor -1 indicating an invalid value.
		kFoundDuplicateTriggers,  ///< Found duplicate trigger decisions.
		kPairTrackIdsAreIdentical,  ///< The track IDs of the track pair are identical.
		kMassValueNotValid,  ///< The invariant mass value is not positive nor -1 indicating an invalid value.
		kLowPtCountInvalid,  ///< The low pT trigger count is greater than 2, which is invalid.
		kHighPtCountInvalid  ///< The high pT trigger count is greater than 2, which is invalid.
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
	 * Methods used to check if the header information corresponds to the
	 * supposed type of the data block.
	 * If the 'reason' parameter is not NULL then these methods will fill the
	 * memory pointed to by reason with a code describing why the header is
	 * not valid, if and only if a problem is found with the data.
	 * These methods will return either kHeaderContainsWrongType or
	 * kHeaderContainsWrongRecordWidth as the reason.
	 */
	static bool HeaderOk(const AliHLTMUONTriggerRecordsBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONTrigRecsDebugBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONTriggerChannelsBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONRecHitsBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONClustersBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONChannelsBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONMansoTracksBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONMansoCandidatesBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONSinglesDecisionBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONPairsDecisionBlockStruct& block, WhyNotValid* reason = NULL);

	/**
	 * Methods used to check more extensively if the integrity of various
	 * types of data blocks are Ok and returns true in that case.
	 * These can be slow and should generally only be used for debugging.
	 */
	static bool IntegrityOk(const AliHLTMUONTriggerRecordStruct& tr, WhyNotValid* reason = NULL);
	
	static bool IntegrityOk(
			const AliHLTMUONTriggerRecordsBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		);
		
	static bool IntegrityOk(const AliHLTMUONTrigRecsDebugBlockStruct& block, WhyNotValid* reason = NULL);
	static bool IntegrityOk(const AliHLTMUONTriggerChannelsBlockStruct& block, WhyNotValid* reason = NULL);
	static bool IntegrityOk(const AliHLTMUONRecHitsBlockStruct& block, WhyNotValid* reason = NULL);
	static bool IntegrityOk(const AliHLTMUONClustersBlockStruct& block, WhyNotValid* reason = NULL);
	static bool IntegrityOk(const AliHLTMUONChannelsBlockStruct& block, WhyNotValid* reason = NULL);
	static bool IntegrityOk(const AliHLTMUONMansoTrackStruct& track, WhyNotValid* reason = NULL);
	
	static bool IntegrityOk(
			const AliHLTMUONMansoTracksBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		);
	
	static bool IntegrityOk(
			const AliHLTMUONMansoCandidatesBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		);
	
	static bool IntegrityOk(const AliHLTMUONTrackDecisionStruct& decision, WhyNotValid* reason = NULL);
	
	static bool IntegrityOk(
			const AliHLTMUONSinglesDecisionBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		);
	
	static bool IntegrityOk(const AliHLTMUONPairDecisionStruct& decision, WhyNotValid* reason = NULL);
	
	static bool IntegrityOk(
			const AliHLTMUONPairsDecisionBlockStruct& block,
			WhyNotValid* reason = NULL, AliHLTUInt32_t* recordNum = NULL
		);

private:
	// Should never have to create or destroy this object.
	AliHLTMUONUtils();
	~AliHLTMUONUtils();
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
