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

///
/// @file   AliHLTMUONUtils.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   17 May 2007
/// @brief  Implementation of AliHLTMUONUtils utility routines.
///

#include "AliHLTMUONUtils.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONTrigRecsDebugBlockStruct.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONClustersBlockStruct.h"
#include "AliHLTMUONChannelsBlockStruct.h"
#include "AliHLTMUONMansoTracksBlockStruct.h"
#include "AliHLTMUONMansoCandidatesBlockStruct.h"
#include "AliHLTMUONSinglesDecisionBlockStruct.h"
#include "AliHLTMUONPairsDecisionBlockStruct.h"
#include "AliMUONTrackerDDLDecoderEventHandler.h"
#include <cstring>
#include <cmath>
#include <cassert>

ClassImp(AliHLTMUONUtils);


AliHLTUInt32_t AliHLTMUONUtils::PackTriggerRecordFlags(
		AliHLTMUONParticleSign sign, const bool hitset[4]
	)
{
	///
	/// This packs the given parameters into the bits of a word appropriate
	/// for AliHLTMUONTriggerRecordStruct::fFlags.
	/// @param sign    The particle sign.
	/// @param hitset  Flags to indicate if the corresponding fHits[i] elements
	///                was set/filled.
	/// @return  Returns the 32 bit packed word.
	///
	
	AliHLTUInt32_t flags;
	switch (sign)
	{
	case kSignMinus: flags = 0x80000000; break;
	case kSignPlus:  flags = 0x40000000; break;
	default:         flags = 0x00000000; break;
	}

	return flags | (hitset[0] ? 0x1 : 0) | (hitset[1] ? 0x2 : 0)
		| (hitset[2] ? 0x4 : 0) | (hitset[3] ? 0x8 : 0);
}


void AliHLTMUONUtils::UnpackTriggerRecordFlags(
		AliHLTUInt32_t flags, AliHLTMUONParticleSign& sign, bool hitset[4]
	)
{
	///
	/// This unpacks the AliHLTMUONTriggerRecordStruct::fFlags bits into
	/// its component fields.
	/// @param flags  The flags from an AliHLTMUONTriggerRecordStruct structure.
	/// @param sign    Sets this to the particle sign.
	/// @param hitset  Sets the array elements to indicate if the corresponding
	///                fHits[i] element was set/filled.
	///
	
	AliHLTUInt32_t signbits = flags & 0xC0000000;
	switch (signbits)
	{
	case 0x80000000: sign = kSignMinus;   break;
	case 0x40000000: sign = kSignPlus;    break;
	default:         sign = kSignUnknown; break;
	}
	hitset[0] = (flags & 0x1) == 0x1;
	hitset[1] = (flags & 0x2) == 0x2;
	hitset[2] = (flags & 0x4) == 0x4;
	hitset[3] = (flags & 0x8) == 0x8;
}


AliHLTUInt32_t AliHLTMUONUtils::PackRecHitFlags(
		AliHLTUInt8_t chamber, AliHLTUInt16_t detElemId
	)
{
	/// This packs the given parameters into the bits of a word appropriate
	/// for AliHLTMUONRecHitStruct::fFlags.
	/// @param chamber    The chamber number in the range [0..13].
	/// @param detElemId  Detector element ID number.
	/// @return  Returns the 32 bit packed word.
	
	return ((chamber & 0xF) << 12) | (detElemId & 0xFFF);
}


void AliHLTMUONUtils::UnpackRecHitFlags(
		AliHLTUInt32_t flags, // [in]
		AliHLTUInt8_t& chamber, // [out]
		AliHLTUInt16_t& detElemId // [out]
	)
{
	/// This unpacks the AliHLTMUONRecHitStruct::fFlags bits into
	/// its component fields.
	/// [in]  @param flags  The flags from an AliHLTMUONRecHitStruct structure.
	/// [out] @param chamber    Sets the chamber number in the range [0..13].
	/// [out] @param detElemId  Sets the detector element ID number.
	
	chamber = (flags >> 12) & 0xF;
	detElemId = flags & 0xFFF;
}


AliHLTUInt32_t AliHLTMUONUtils::PackTrackDecisionBits(bool highPt, bool lowPt)
{
	///
	/// This packs the given parameters into the bits of a word appropriate
	/// for AliHLTMUONTrackDecisionStruct::fTriggerBits.
	/// @param highPt  Has the track passed the high pt cut.
	/// @param lowPt   Has the track passed the low pt cut.
	/// @return  Returns the 32 bit packed word.
	///
	
	return (highPt ? 0x2 : 0) | (lowPt ? 0x1 : 0);
}


void AliHLTMUONUtils::UnpackTrackDecisionBits(
		AliHLTUInt32_t bits, bool& highPt, bool& lowPt
	)
{
	///
	/// This unpacks the AliHLTMUONTrackDecisionStruct::fTriggerBits bits into
	/// its component fields.
	/// @param bits  The trigger bits from an AliHLTMUONTrackDecisionStruct
	///              structure.
	/// @param highPt Sets this to the value of the high pt cut bit.
	/// @param lowPt  Sets this to the value of the low pt cut bit.
	///
	
	lowPt  = (bits & 0x1) == 0x1;
	highPt = (bits & 0x2) == 0x2;
}


AliHLTUInt32_t AliHLTMUONUtils::PackPairDecisionBits(
		bool highMass, bool lowMass, bool unlike,
		AliHLTUInt8_t highPtCount, AliHLTUInt8_t lowPtCount
	)
{
	///
	/// This packs the given parameters into the bits of a word appropriate
	/// for AliHLTMUONPairDecisionStruct::fTriggerBits.
	///
	/// @param highMass Has the track pair passed the high invariant mass cut.
	/// @param lowMass  Has the track pair passed the low invariant mass cut.
	/// @param unlike   Does the track pair have unlike signs.
	/// @param highPtCount The number of tracks that passed the high pt cut
	///                    in the pair.
	/// @param lowPtCount  The number of tracks that passed the low pt cut
	///                    in the pair.
	/// @return  Returns the 32 bit packed word.
	///
	/// Note: Must have highPtCount <= 2, lowPtCount <= 2 and
	/// unlike == true if highMass or lowMass is true.
	///
	
	assert( lowPtCount <= 2 );
	assert( highPtCount <= 2 );
	// highMass and lowMass must be false if unlike is false:
	assert( not unlike ? (highMass == false and lowMass == false) : true );
	
	return (highMass ? 0x40 : 0) | (lowMass ? 0x20 : 0) | (unlike ? 0x10 : 0)
		| ((highPtCount & 0x3) << 2) | (lowPtCount & 0x3);
}


void AliHLTMUONUtils::UnpackPairDecisionBits(
		AliHLTUInt32_t bits, bool& highMass, bool& lowMass, bool& unlike,
		AliHLTUInt8_t& highPtCount, AliHLTUInt8_t& lowPtCount
	)
{
	///
	/// This unpacks the AliHLTMUONPairDecisionStruct::fTriggerBits bits into
	/// its component fields.
	/// @param bits  The trigger bits from an AliHLTMUONPairDecisionStruct
	///              structure.
	/// @param highMass Sets this to the value of the high invariant mass cut bit.
	/// @param lowMass  Sets this to the value of the low invariant mass cut bit.
	/// @param unlike   Sets this if the pair is unlike sign.
	/// @param highPtCount Sets this to the high pt count bits.
	/// @param lowPtCount  Sets this to the low pt count bits.
	///
	
	highMass = (bits & 0x40) == 0x40;
	lowMass  = (bits & 0x20) == 0x20;
	unlike   = (bits & 0x10) == 0x10;
	highPtCount = (bits & 0xC) >> 2;
	lowPtCount = bits & 0x3;
}


AliHLTUInt32_t AliHLTMUONUtils::PackSpecBits(
		const bool ddl[22]
	)
{
	///
	/// This packs the given parameters into the 32bit Pub/Sub specification
	/// word in the data block descriptor.
	///
	/// @param ddl  The list of DDLs forming part of the readout. ddl[0]
	///             indicates DDL number 2560, ddl[1] is for DDL 2561 and so
	///             on up to ddl[19]. ddl[20] and ddl[21] will be for the
	///             trigger DDLs 2816 and 2817 respectively.
	/// @return  Returns the 32 bit packed specification word.
	///
	
	// Pack the bits into the following format:
	//   bit:   [        31 - 22        ][     21     ][     20     ][  19 - 0 ]
	//   field: [ reserved, set to zero ][ TRGDDL2817 ][ TRGDDL2816 ][ TRKDDLS ]
	// Meaning of field acronyms:
	//   TRGDDL2816 - Trigger DDL number 2816.
	//   TRGDDL2817 - Trigger DDL number 2817.
	//   TRKDDLS - Tracking DDL flags where bit 0 will be for DDL number 2560,
	//             bit 1 for DDL no. 2561 etc. up to bit 19 which is for DDL 2579.
	AliHLTUInt32_t bits = 0;
	for (int i = 0; i < 22; i++)
		bits |= (ddl[i] ? 0x1 : 0x0) << i;
	return bits;
}


void AliHLTMUONUtils::UnpackSpecBits(
		AliHLTUInt32_t bits, bool ddl[22]
	)
{
	///
	/// This unpacks the AliHLTMUONPairDecisionStruct::fTriggerBits bits into
	/// its component fields.
	/// @param bits  The Pub/Sub specification word from a data block descriptor.
	/// @param ddl  The output list of DDLs forming part of the readout. ddl[0]
	///             indicates DDL number 2560, ddl[1] is for DDL 2561 and so
	///             on up to ddl[19]. ddl[20] and ddl[21] will be for the
	///             trigger DDLs 2816 and 2817 respectively.
	///
	
	// Perform the inverse operation of PackSpecBits.
	for (int i = 0; i < 22; i++)
		ddl[i] = ((bits >> i) & 0x1) == 1;
}


AliHLTInt32_t AliHLTMUONUtils::DDLNumberToEquipId(AliHLTInt32_t ddlNo)
{
	///
	/// This method converts the DDL number for the muon spectrometer in the
	/// range [0..21] to the equipment ID number.
	/// @param ddlNo  The DDL number in the range [0..21].
	/// @return  Returns the equipment ID number or -1 if ddlNo was invalid.
	///
	
	if (0 <= ddlNo and ddlNo <= 19)
	{
		return 2560 + ddlNo;
	}
	else if (20 <= ddlNo and ddlNo <= 21)
	{
		return 2816 + (ddlNo - 20);
	}
	else
	{
		return -1;
	}
}


AliHLTInt32_t AliHLTMUONUtils::EquipIdToDDLNumber(AliHLTInt32_t id)
{
	///
	/// This method converts the equipment ID number for a muon spectrometer
	/// DDL to the DDL number in the range [0..21].
	/// @param id  The equipment ID of the DDL.
	/// @return  Returns the DDL number in the range [0..21] or -1 if the
	///          equipment ID was invalid.
	///
	
	if (2560 <= id and id <= 2560+19)
	{
		return id - 2560;
	}
	else if (2816 <= id and id <= 2817)
	{
		return id - 2816 + 20;
	}
	else
	{
		return -1;
	}
}


AliHLTInt32_t AliHLTMUONUtils::SpecToEquipId(AliHLTUInt32_t spec)
{
	///
	/// This method converts a 32 bit data block specification for a MUON-HLT
	/// data block into its corresponding DDL equipment ID number.
	/// It is assumed that the specification is for a data block comming from
	/// a single DDL source. If more than one DDL contributed to the data block
	/// then -1 is returned.
	/// @param spec  The 32 bit specification for a data block.
	/// @return  Returns the equipment ID corresponding to the specification
	///          or -1 if the specification was invalid.
	///
	
	for (AliHLTInt32_t ddlNo = 0; ddlNo < 20; ddlNo++)
	{
		if (spec == AliHLTUInt32_t(0x1 << ddlNo))
			return ddlNo + 2560;
	}
	for (AliHLTInt32_t ddlNo = 20; ddlNo < 22; ddlNo++)
	{
		if (spec == AliHLTUInt32_t(0x1 << ddlNo))
			return ddlNo - 20 + 2816;
	}
	return -1;
}


AliHLTUInt32_t AliHLTMUONUtils::EquipIdToSpec(AliHLTInt32_t id)
{
	///
	/// This method converts a equipment ID number for a DDL into its corresponding
	/// 32 bit data block specification for the MUON-HLT.
	/// @param id  The equipment ID number of the DDL.
	/// @return  Returns the 32 bit data block specification or 0x0 if the
	///          equipment ID was invalid.
	///
	
	if (2560 <= id and id <= 2560+19)
	{
		return 0x1 << (id - 2560);
	}
	else if (2816 <= id and id <= 2817)
	{
		return 0x1 << (id - 2816 + 20);
	}
	else
	{
		return 0x0;
	}
}


AliHLTInt32_t AliHLTMUONUtils::SpecToDDLNumber(AliHLTUInt32_t spec)
{
	///
	/// This method converts a 32 bit data block specification for a MUON-HLT
	/// data block into its corresponding DDL number in the range [0..21].
	/// It is assumed that the specification is for a data block comming from
	/// a single DDL source. If more than one DDL contributed to the data block
	/// then -1 is returned.
	/// @param spec  The 32 bit specification for a data block.
	/// @return  Returns the corresponding DDL number for the specification
	///          or -1 if the specification was invalid.
	///
	
	for (AliHLTInt32_t ddlNo = 0; ddlNo < 22; ddlNo++)
	{
		if (spec == AliHLTUInt32_t(0x1 << ddlNo))
			return ddlNo;
	}
	return -1;
}


AliHLTUInt32_t AliHLTMUONUtils::DDLNumberToSpec(AliHLTInt32_t ddlNo)
{
	///
	/// This method converts a DDL number in the range [0..21] into its
	/// corresponding 32 bit data block specification for the MUON-HLT.
	/// @param ddlNo  The equipment ID number of the DDL.
	/// @return  Returns the 32 bit data block specification or 0x0 if the
	///          DDL number was invalid (out of range).
	///
	
	if (0 <= ddlNo and ddlNo <= 21)
	{
		return 0x1 << ddlNo;
	}
	else
	{
		return 0x0;
	}
}


AliHLTMUONDataBlockType AliHLTMUONUtils::ParseCommandLineTypeString(const char* type)
{
	/// Parses the string containing the type name of a dHLT data block and
	/// returns the corresponding AliHLTMUONDataBlockType value.
	/// \param  type  The string containing the type name.
	/// \returns  The data block type or kUnknownDataBlock if the type name
	///      is invalid.

	if (strcmp(type, "trigrecs") == 0)
	{
		return kTriggerRecordsDataBlock;
	}
	else if (strcmp(type, "trigrecsdebug") == 0)
	{
		return kTrigRecsDebugDataBlock;
	}
	else if (strcmp(type, "rechits") == 0)
	{
		return kRecHitsDataBlock;
	}
	else if (strcmp(type,"channels") == 0)
	{
		return kChannelsDataBlock;
	}
	else if (strcmp(type,"clusters") == 0)
	{
		return kClustersDataBlock;
	}
	else if (strcmp(type, "mansotracks") == 0)
	{
		return kMansoTracksDataBlock;
	}
	else if (strcmp(type, "mansocandidates") == 0)
	{
		return kMansoCandidatesDataBlock;
	}
	else if (strcmp(type, "singlesdecision") == 0)
	{
		return kSinglesDecisionDataBlock;
	}
	else if (strcmp(type, "pairsdecision") == 0)
	{
		return kPairsDecisionDataBlock;
	}
	
	return kUnknownDataBlock;
}


const char* AliHLTMUONUtils::DataBlockTypeToString(AliHLTMUONDataBlockType type)
{
	/// Converts a type ID to a type string compatible with
	/// HLT data types.
	
	static char str[kAliHLTComponentDataTypefIDsize+1];
	AliHLTComponentDataType t;
	switch (type)
	{
	case kTriggerRecordsDataBlock:
		t = AliHLTMUONConstants::TriggerRecordsBlockDataType();
		break;
	case kTrigRecsDebugDataBlock:
		t = AliHLTMUONConstants::TrigRecsDebugBlockDataType();
		break;
	case kRecHitsDataBlock:
		t = AliHLTMUONConstants::RecHitsBlockDataType();
		break;
	case kClustersDataBlock:
		t = AliHLTMUONConstants::ClusterBlockDataType();
		break;
	case kChannelsDataBlock:
		t = AliHLTMUONConstants::ChannelBlockDataType();
		break;
	case kMansoTracksDataBlock:
		t = AliHLTMUONConstants::MansoTracksBlockDataType();
		break;
	case kMansoCandidatesDataBlock:
		t = AliHLTMUONConstants::MansoCandidatesBlockDataType();
		break;
	case kSinglesDecisionDataBlock:
		t = AliHLTMUONConstants::SinglesDecisionBlockDataType();
		break;
	case kPairsDecisionDataBlock:
		t = AliHLTMUONConstants::PairsDecisionBlockDataType();
		break;
	default:
		return "UNKNOWN";
	}
	memcpy(&str, &t.fID, kAliHLTComponentDataTypefIDsize);
	// Must insert the NULL character to make this an ANSI C string.
	str[kAliHLTComponentDataTypefIDsize] = '\0';
	return &str[0];
}


const char* AliHLTMUONUtils::FailureReasonToString(WhyNotValid reason)
{
	/// This method converts the WhyNotValid enumeration to a string representation.
	
	switch (reason)
	{
	case kNoReason: return "kNoReason";
	case kHeaderContainsWrongType: return "kHeaderContainsWrongType";
	case kHeaderContainsWrongRecordWidth: return "kHeaderContainsWrongRecordWidth";
	case kInvalidIdValue: return "kInvalidIdValue";
	case kInvalidTriggerIdValue: return "kInvalidTriggerIdValue";
	case kInvalidTrackIdValue: return "kInvalidTrackIdValue";
	case kReservedBitsNotZero: return "kReservedBitsNotZero";
	case kParticleSignBitsNotValid: return "kParticleSignBitsNotValid";
	case kHitNotMarkedAsNil: return "kHitNotMarkedAsNil";
	case kInvalidDetElementNumber: return "kInvalidDetElementNumber";
	case kInvalidChamberNumber: return "kInvalidChamberNumber";
	case kHitIsNil: return "kHitIsNil";
	case kInvalidChannelCount: return "kInvalidChannelCount";
	case kInvalidTotalCharge: return "kInvalidTotalCharge";
	case kInvalidBusPatchId: return "kInvalidBusPatchId";
	case kInvalidManuId: return "kInvalidManuId";
	case kInvalidChannelAddress: return "kInvalidChannelAddress";
	case kInvalidSignal: return "kInvalidSignal";
	case kDataWordDifferent: return "kDataWordDifferent";
	case kChiSquareInvalid: return "kChiSquareInvalid";
	case kMomentumVectorNotZero: return "kMomentumVectorNotZero";
	case kRoiRadiusInvalid: return "kRoiRadiusInvalid";
	case kHitNotWithinRoi: return "kHitNotWithinRoi";
	case kPtValueNotValid: return "kPtValueNotValid";
	case kPairTrackIdsAreIdentical: return "kPairTrackIdsAreIdentical";
	case kMassValueNotValid: return "kMassValueNotValid";
	case kLowPtCountInvalid: return "kLowPtCountInvalid";
	case kHighPtCountInvalid: return "kHighPtCountInvalid";
	case kFoundDuplicateIDs: return "kFoundDuplicateIDs";
	case kFoundDuplicateHits: return "kFoundDuplicateHits";
	case kFoundDuplicateTriggers: return "kFoundDuplicateTriggers";
	default: return "INVALID";
	}
}


const char* AliHLTMUONUtils::FailureReasonToMessage(WhyNotValid reason)
{
	/// This method returns a string containing a user readable message explaining
	/// the reason for failure described by the WhyNotValid enumeration.
	
	switch (reason)
	{
	case kNoReason:
		return "There was no problem with the data block.";
	case kHeaderContainsWrongType:
		return "The common data header contains an incorrect type"
			" identifier.";
	case kHeaderContainsWrongRecordWidth:
		return "The common data header contains an incorrect data"
			" record width.";
	case kInvalidIdValue:
		return "The structure identifier does not have a valid value.";
	case kInvalidTriggerIdValue:
		return "The trigger structure identifier does not have a valid"
			" value.";
	case kInvalidTrackIdValue:
		return "The track structure identifier does not have a valid"
			" value.";
	case kReservedBitsNotZero:
		return "Reserved bits have not been set to zero.";
	case kParticleSignBitsNotValid:
		return "The particle sign bits are not a valid value.";
	case kHitNotMarkedAsNil:
		return "A hit was marked as not found, but the corresponding hit"
			" structure was not set to nil.";
	case kInvalidDetElementNumber:
		return "An invalid detector element ID was found.";
	case kInvalidChamberNumber:
		return "An invalid chamber number was found.";
	case kHitIsNil:
		return "The hit cannot be set to a nil value.";
	case kInvalidChannelCount:
		return "The number of channels indicated is zero or outside"
			" the valid range.";
	case kInvalidTotalCharge:
		return "The total charge does not have a valid value.";
	case kInvalidBusPatchId:
		return "The bus patch identifier is outside the valid range.";
	case kInvalidManuId:
		return "The MANU identifier is outside the valid range.";
	case kInvalidChannelAddress:
		return "The MANU channel address is outside the valid range.";
	case kInvalidSignal:
		return "The ADC signal value is outside the valid range.";
	case kDataWordDifferent:
		return "The raw data word is different from the unpacked values.";
	case kChiSquareInvalid:
		return "The chi squared value must be a positive value or -1"
			" indicating a fitting error.";
	case kMomentumVectorNotZero:
		return "The chi sqaured value is set to -1 indicating momentum"
			" was not fitted, but the momentum vector was not zero.";
	case kRoiRadiusInvalid:
		return "The region of interest radius is invalid.";
	case kHitNotWithinRoi:
		return "A tracks hit is not within the corresponding region"
			" of interest.";
	case kPtValueNotValid:
		return "The pT value is not positive, nor -1 indicating an"
			" invalid value.";
	case kPairTrackIdsAreIdentical:
		return "The track identifiers of the track pair are identical.";
	case kMassValueNotValid:
		return "The invariant mass value is not positive, nor -1"
			" indicating an invalid value.";
	case kLowPtCountInvalid:
		return "The low pT trigger count is greater than 2,"
			" which is invalid.";
	case kHighPtCountInvalid:
		return "The high pT trigger count is greater than 2,"
			" which is invalid.";
	case kFoundDuplicateIDs:
		return "Found duplicate data record identifiers, but they"
			" should all be unique.";
	case kFoundDuplicateHits:
		return "Found duplicate hit structures, but they should all"
			" be unique.";
	case kFoundDuplicateTriggers:
		return "Found duplicate trigger decisions.";
	default:
		return "UNKNOWN REASON CODE";
	}
}


bool AliHLTMUONUtils::RecordNumberWasSet(WhyNotValid reason)
{
	/// Returns true if the \em recordNum in the corresponding IntegrityOk method
	/// would have been set, if it returned false and a reason was set.
	/// This helper method makes it easy to test if the \em recordNum parameter
	/// is filled with a valid value or not.
	/// \param reason  The reason code as returned by the IntegrityOk method.
	/// \returns  true if the \em recordNum parameter was set for the given
	///      reason code.
	
	switch (reason)
	{
	case kInvalidIdValue:
	case kInvalidTriggerIdValue:
	case kInvalidTrackIdValue:
	case kReservedBitsNotZero:
	case kParticleSignBitsNotValid:
	case kHitNotMarkedAsNil:
	case kInvalidDetElementNumber:
	case kInvalidChamberNumber:
	case kHitIsNil:
	case kInvalidChannelCount:
	case kInvalidBusPatchId:
	case kInvalidManuId:
	case kInvalidChannelAddress:
	case kInvalidSignal:
	case kDataWordDifferent:
	case kChiSquareInvalid:
	case kPtValueNotValid:
	case kPairTrackIdsAreIdentical:
	case kMassValueNotValid:
	case kLowPtCountInvalid:
	case kHighPtCountInvalid:
		return true;
	default: return false;
	}
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONTriggerRecordsBlockStruct& block,
		WhyNotValid* reason, AliHLTUInt32_t& reasonCount
	)
{
	/// Method used to check if the header information corresponds to the
	/// supposed type of the raw dHLT data block.
	/// [in]  \param block  The data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the header is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the header and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The block must have the correct type.
	if (block.fHeader.fType != kTriggerRecordsDataBlock)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongType;
			reasonCount++;
		}
		result = false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTriggerRecordStruct))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongRecordWidth;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONTrigRecsDebugBlockStruct& block,
		WhyNotValid* reason, AliHLTUInt32_t& reasonCount
	)
{
	/// Method used to check if the header information corresponds to the
	/// supposed type of the raw dHLT data block.
	/// [in]  \param block  The data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the header is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the header and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The block must have the correct type.
	if (block.fHeader.fType != kTrigRecsDebugDataBlock)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongType;
			reasonCount++;
		}
		result = false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTrigRecInfoStruct))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongRecordWidth;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONRecHitsBlockStruct& block,
		WhyNotValid* reason, AliHLTUInt32_t& reasonCount
	)
{
	/// Method used to check if the header information corresponds to the
	/// supposed type of the raw dHLT data block.
	/// [in]  \param block  The data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the header is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the header and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The block must have the correct type.
	if (block.fHeader.fType != kRecHitsDataBlock)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongType;
			reasonCount++;
		}
		result = false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONRecHitStruct))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongRecordWidth;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONClustersBlockStruct& block,
		WhyNotValid* reason, AliHLTUInt32_t& reasonCount
	)
{
	/// Method used to check if the header information corresponds to the
	/// supposed type of the raw dHLT data block.
	/// [in]  \param block  The data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the header is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the header and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The block must have the correct type.
	if (block.fHeader.fType != kClustersDataBlock)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongType;
			reasonCount++;
		}
		result = false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONClusterStruct))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongRecordWidth;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONChannelsBlockStruct& block,
		WhyNotValid* reason, AliHLTUInt32_t& reasonCount
	)
{
	/// Method used to check if the header information corresponds to the
	/// supposed type of the raw dHLT data block.
	/// [in]  \param block  The data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the header is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the header and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The block must have the correct type.
	if (block.fHeader.fType != kChannelsDataBlock)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongType;
			reasonCount++;
		}
		result = false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONChannelStruct))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongRecordWidth;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONMansoTracksBlockStruct& block,
		WhyNotValid* reason, AliHLTUInt32_t& reasonCount
	)
{
	/// Method used to check if the header information corresponds to the
	/// supposed type of the raw dHLT data block.
	/// [in]  \param block  The data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the header is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the header and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The block must have the correct type.
	if (block.fHeader.fType != kMansoTracksDataBlock)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongType;
			reasonCount++;
		}
		result = false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONMansoTrackStruct))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongRecordWidth;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONMansoCandidatesBlockStruct& block,
		WhyNotValid* reason, AliHLTUInt32_t& reasonCount
	)
{
	/// Method used to check if the header information corresponds to the
	/// supposed type of the raw dHLT data block.
	/// [in]  \param block  The data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the header is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the header and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The block must have the correct type.
	if (block.fHeader.fType != kMansoCandidatesDataBlock)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongType;
			reasonCount++;
		}
		result = false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONMansoCandidateStruct))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongRecordWidth;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONSinglesDecisionBlockStruct& block,
		WhyNotValid* reason, AliHLTUInt32_t& reasonCount
	)
{
	/// Method used to check if the header information corresponds to the
	/// supposed type of the raw dHLT data block.
	/// [in]  \param block  The data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the header is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the header and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The block must have the correct type.
	if (block.fHeader.fType != kSinglesDecisionDataBlock)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongType;
			reasonCount++;
		}
		result = false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTrackDecisionStruct))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongRecordWidth;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONPairsDecisionBlockStruct& block,
		WhyNotValid* reason, AliHLTUInt32_t& reasonCount
	)
{
	/// Method used to check if the header information corresponds to the
	/// supposed type of the raw dHLT data block.
	/// [in]  \param block  The data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the header is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the header and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The block must have the correct type.
	if (block.fHeader.fType != kPairsDecisionDataBlock)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongType;
			reasonCount++;
		}
		result = false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONPairDecisionStruct))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHeaderContainsWrongRecordWidth;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONTriggerRecordStruct& tr,
		WhyNotValid* reason,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// trigger record structure is OK and returns true in that case.
	/// [in] \param tr  The trigger record structure to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the structure is
	///      not valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the structure and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// Check that the ID has a valid value.
	if (not (tr.fId >= 0 or tr.fId == -1))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidIdValue;
			reasonCount++;
		}
		result = false;
	}
	
	// Make sure that the reserved bits in the fFlags field are set
	// to zero.
	if ((tr.fFlags & 0x3FFFFFF0) != 0)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kReservedBitsNotZero;
			reasonCount++;
		}
		result = false;
	}

	// Make sure the sign is not invalid.
	if ((tr.fFlags & 0xC0000000) == 0xC0000000)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kParticleSignBitsNotValid;
			reasonCount++;
		}
		result = false;
	}

	// Check that fHit[i] is nil if the corresponding bit in the
	// flags word is zero.
	const AliHLTMUONRecHitStruct& nilhit
		= AliHLTMUONConstants::NilRecHitStruct();
	if ( ((tr.fFlags & 0x1) == 0 and tr.fHit[0] != nilhit) or
	     ((tr.fFlags & 0x2) == 0 and tr.fHit[1] != nilhit) or
	     ((tr.fFlags & 0x4) == 0 and tr.fHit[2] != nilhit) or
	     ((tr.fFlags & 0x8) == 0 and tr.fHit[3] != nilhit)
	   )
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHitNotMarkedAsNil;
			reasonCount++;
		}
		result = false;
	}
	
	// Check the individual hits
	for (int i = 0; i < 4; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(tr.fHit[i], reason + reasonCount, filledCount))
		{
			reasonCount += filledCount;
			result = false;
		}
	}

	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONTriggerRecordsBlockStruct& block,
		WhyNotValid* reason,
		AliHLTUInt32_t* recordNum,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// dHLT raw internal data block is OK and returns true in that case.
	/// [in] \param block  The trigger record data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the data block is
	///      not valid.
	/// [out] \param recordNum  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the number of the trigger record that had a problem.
	///      The value 'recordNum[i]' will only contain a valid value if
	///      the corresponding 'reason[i]' contains one of:
	///        - kInvalidIdValue
	///        - kReservedBitsNotZero
	///        - kParticleSignBitsNotValid
	///        - kHitNotMarkedAsNil
	///        - kInvalidDetElementNumber
	///        - kInvalidChamberNumber
	/// \note You can use RecordNumberWasSet(reason[i]) to check if 'recordNum[i]'
	///      was set and is valid or not.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason' and 'recordNum'. It will be filled
	///      with the number of items actually filled into the arrays upon exit
	///      from this method.
	/// \returns  true if there is no problem with the data and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = HeaderOk(block, reason, reasonCount);
	
	const AliHLTMUONTriggerRecordStruct* triggerRecord =
		reinterpret_cast<const AliHLTMUONTriggerRecordStruct*>(&block + 1);

	// Check if any ID is duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = triggerRecord[i].fId;
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (id == triggerRecord[j].fId)
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateIDs;
					reasonCount++;
				}
				result = false;
			}
		}
	}

	// Check integrity of individual trigger records.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(triggerRecord[i], reason+reasonCount, filledCount))
		{
			// reasons filled in IntegrityOk, now we just need to adjust
			// reasonCount and fill the recordNum values.
			if (recordNum != NULL)
			{
				for (AliHLTUInt32_t n = 0; n < filledCount; n++)
					recordNum[reasonCount + n] = i;
			}
			reasonCount += filledCount;
			result = false;
		}
	}

	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONTrigRecInfoStruct& trigInfo,
		WhyNotValid* reason,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// trigger record debug information structure is OK and returns true in that case.
	/// [in] \param trigInfo  The trigger record debug information structure to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the structure is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the structure and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// Check that the trigger ID has a valid value.
	if (not (trigInfo.fTrigRecId >= 0 or trigInfo.fTrigRecId == -1))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidIdValue;
			reasonCount++;
		}
		result = false;
	}

	// Check that the fDetElemId[i] numbers are valid.
	if ( not ((trigInfo.fDetElemId[0] >= 100 and trigInfo.fDetElemId[0] < 1500)
	          or trigInfo.fDetElemId[0] == -1)
	     or not ((trigInfo.fDetElemId[1] >= 100 and trigInfo.fDetElemId[1] < 1500)
	          or trigInfo.fDetElemId[1] == -1)
	     or not ((trigInfo.fDetElemId[2] >= 100 and trigInfo.fDetElemId[2] < 1500)
	          or trigInfo.fDetElemId[2] == -1)
	     or not ((trigInfo.fDetElemId[3] >= 100 and trigInfo.fDetElemId[3] < 1500)
	          or trigInfo.fDetElemId[3] == -1)
	   )
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidDetElementNumber;
			reasonCount++;
		}
		result = false;
	}

	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONTrigRecsDebugBlockStruct& block,
		WhyNotValid* reason,
		AliHLTUInt32_t* recordNum,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// dHLT raw internal data block is OK and returns true in that case.
	/// [in] \param block  The trigger record debugging information data block
	///      to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the data block is
	///      not valid.
	/// [out] \param recordNum  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the number of the trigger record debug information
	///      structure that had a problem.
	///      The value 'recordNum[i]' will only contain a valid value if
	///      the corresponding 'reason[i]' contains one of:
	///        - kInvalidIdValue
	///        - kInvalidDetElementNumber
	/// \note You can use RecordNumberWasSet(reason[i]) to check if 'recordNum[i]'
	///      was set and is valid or not.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason' and 'recordNum'. It will be filled
	///      with the number of items actually filled into the arrays upon exit
	///      from this method.
	/// \returns  true if there is no problem with the data and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = HeaderOk(block, reason, reasonCount);
	
	const AliHLTMUONTrigRecInfoStruct* triggerInfo =
		reinterpret_cast<const AliHLTMUONTrigRecInfoStruct*>(&block + 1);

	// Check if any trigger debug info structure has duplicated trigger IDs.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = triggerInfo[i].fTrigRecId;
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (id == triggerInfo[j].fTrigRecId)
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateIDs;
					reasonCount++;
				}
				result = false;
			}
		}
	}

	// Check integrity of individual trigger records.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(triggerInfo[i], reason+reasonCount, filledCount))
		{
			// reasons filled in IntegrityOk, now we just need to adjust
			// reasonCount and fill the recordNum values.
			if (recordNum != NULL)
			{
				for (AliHLTUInt32_t n = 0; n < filledCount; n++)
					recordNum[reasonCount + n] = i;
			}
			reasonCount += filledCount;
			result = false;
		}
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONRecHitStruct& hit,
		WhyNotValid* reason,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// reconstructed hit structure is OK and returns true in that case.
	/// [in] \param hit  The reconstructed hit structure to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the structure is
	///      not valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the structure and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// If this is a NIL hit then skip all other checks.
	if (hit == AliHLTMUONConstants::NilRecHitStruct())
	{
		return true;
	}
	
	// Make sure that the reserved bits in the fFlags field are set
	// to zero.
	if ((hit.fFlags & 0x3FFF0000) != 0)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kReservedBitsNotZero;
			reasonCount++;
		}
		result = false;
	}

	AliHLTUInt32_t detElemId = hit.fFlags & 0x00000FFF;
	AliHLTUInt32_t chamber = (hit.fFlags & 0x0000F000) >> 12;
	
	// Make sure the detector element ID number is valid.
	if (not (detElemId >= 100 and detElemId < 1500))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidDetElementNumber;
			reasonCount++;
		}
		result = false;
	}
	
	// Make sure the chamber number is valid.
	if (((detElemId / 100) - 1) != chamber or chamber >= 14)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidChamberNumber;
			reasonCount++;
		}
		result = false;
	}

	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONRecHitsBlockStruct& block,
		WhyNotValid* reason,
		AliHLTUInt32_t* recordNum,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// dHLT raw internal hits data block is OK and returns true in that case.
	/// [in] \param block  The reconstructed hits data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the data block is
	///      not valid.
	/// [out] \param recordNum  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the number of the reconstructed hits that had a problem.
	///      The value 'recordNum[i]' will only contain a valid value if
	///      the corresponding 'reason[i]' contains one of:
	///        - kReservedBitsNotZero
	///        - kInvalidDetElementNumber
	///        - kInvalidChamberNumber
	/// \note You can use RecordNumberWasSet(reason[i]) to check if 'recordNum[i]'
	///      was set and is valid or not.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason' and 'recordNum'. It will be filled
	///      with the number of items actually filled into the arrays upon exit
	///      from this method.
	/// \returns  true if there is no problem with the data and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = HeaderOk(block, reason, reasonCount);
	
	const AliHLTMUONRecHitStruct* hit =
		reinterpret_cast<const AliHLTMUONRecHitStruct*>(&block + 1);

	// Check if any hit structure has been duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		const AliHLTMUONRecHitStruct& h = hit[i];
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (h == hit[j])
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateHits;
					reasonCount++;
				}
				result = false;
			}
		}
	}

	// Check integrity of the individual hit structures.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(hit[i], reason+reasonCount, filledCount))
		{
			// reasons filled in IntegrityOk, now we just need to adjust
			// reasonCount and fill the recordNum values.
			if (recordNum != NULL)
			{
				for (AliHLTUInt32_t n = 0; n < filledCount; n++)
					recordNum[reasonCount + n] = i;
			}
			reasonCount += filledCount;
			result = false;
		}
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONClusterStruct& cluster,
		WhyNotValid* reason,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// cluster structure is OK and returns true in that case.
	/// [in] \param cluster  The cluster structure to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the structure is
	///      not valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the structure and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// Check that the cluster ID has a valid value.
	if (not (cluster.fId >= 0 or cluster.fId == -1))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidIdValue;
			reasonCount++;
		}
		result = false;
	}
	
	// Check that the cluster does not have a nil value for its hit.
	if (cluster.fHit == AliHLTMUONConstants::NilRecHitStruct())
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHitIsNil;
			reasonCount++;
		}
		result = false;
	}
	
	// Make sure the detector element is a valid value.
	if (not ((cluster.fDetElemId >= 100 and cluster.fDetElemId < 1500)
	    or cluster.fDetElemId == -1)
	   )
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidDetElementNumber;
			reasonCount++;
		}
		result = false;
	}
	
	// The number of channels should be in a reasonable range.
	// between 1 and the maximum number of channels per DDL.
	// 1<<17 taken from the 11 bits MANU ID + 6 bits channel address.
	if (cluster.fNchannels < 1 or (1<<17) < cluster.fNchannels)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidChannelCount;
			reasonCount++;
		}
		result = false;
	}
	
	// The charge must be a positive value of -1.
	if (not (cluster.fCharge >= 0 or cluster.fCharge == -1))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidTotalCharge;
			reasonCount++;
		}
		result = false;
	}

	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONClustersBlockStruct& block,
		WhyNotValid* reason,
		AliHLTUInt32_t* recordNum,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// dHLT internal clusters data block is OK and returns true in that case.
	/// [in] \param block  The clusters data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the data block is
	///      not valid.
	/// [out] \param recordNum  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the number of the cluster structure that had a problem.
	///      The value 'recordNum[i]' will only contain a valid value if
	///      the corresponding 'reason[i]' contains one of:
	///        - kInvalidIdValue
	///        - kHitIsNil
	///        - kInvalidDetElementNumber
	///        - kInvalidChannelCount
	/// \note You can use RecordNumberWasSet(reason[i]) to check if 'recordNum[i]'
	///      was set and is valid or not.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason' and 'recordNum'. It will be filled
	///      with the number of items actually filled into the arrays upon exit
	///      from this method.
	/// \returns  true if there is no problem with the data and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = HeaderOk(block, reason, reasonCount);

	const AliHLTMUONClusterStruct* cluster =
		reinterpret_cast<const AliHLTMUONClusterStruct*>(&block + 1);
	
	// Check if any ID is duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = cluster[i].fId;
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (id == cluster[j].fId)
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateIDs;
					reasonCount++;
				}
				result = false;
			}
		}
	}

	// Check if any hit structure has been duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		const AliHLTMUONRecHitStruct& h = cluster[i].fHit;
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (h == cluster[j].fHit)
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateHits;
					reasonCount++;
				}
				result = false;
			}
		}
	}

	// Check integrity of individual cluster structures.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(cluster[i], reason+reasonCount, filledCount))
		{
			// reasons filled in IntegrityOk, now we just need to adjust
			// reasonCount and fill the recordNum values.
			if (recordNum != NULL)
			{
				for (AliHLTUInt32_t n = 0; n < filledCount; n++)
					recordNum[reasonCount + n] = i;
			}
			reasonCount += filledCount;
			result = false;
		}
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONChannelStruct& channel,
		WhyNotValid* reason,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// channel structure is OK and returns true in that case.
	/// [in] \param cluster  The channel structure to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the structure is
	///      not valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the structure and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// Check that the channel ID has a valid value.
	if (not (channel.fClusterId >= 0 or channel.fClusterId == -1))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidIdValue;
			reasonCount++;
		}
		result = false;
	}
	
	// Check that the bus patch ID has a valid value, which fits into 12 bits.
	if ((channel.fBusPatch & (~0xFFF)) != 0)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidBusPatchId;
			reasonCount++;
		}
		result = false;
	}
	
	// Check that the MANU ID has a valid value, which fits into 11 bits.
	if ((channel.fManu & (~0x7FF)) != 0)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidManuId;
			reasonCount++;
		}
		result = false;
	}
	
	// Check that the channel address has a valid value, which fits into 6 bits.
	if ((channel.fChannelAddress & (~0x3F)) != 0)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidChannelAddress;
			reasonCount++;
		}
		result = false;
	}
	
	// Check that the ADC signal has a valid value, which fits into 12 bits.
	if ((channel.fSignal & (~0xFFF)) != 0)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidSignal;
			reasonCount++;
		}
		result = false;
	}
	
	// Check that the raw data word corresponds to the unpacked values for
	// the ADC signal, MANU ID and channel address.
	UShort_t manuId; UChar_t channelId; UShort_t adc;
	AliMUONTrackerDDLDecoderEventHandler::UnpackADC(
			channel.fRawDataWord, manuId, channelId, adc
		);
	if (manuId != channel.fManu or channelId != channel.fChannelAddress
	    or adc != channel.fSignal
	   )
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kDataWordDifferent;
			reasonCount++;
		}
		result = false;
	}

	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONChannelsBlockStruct& block,
		WhyNotValid* reason,
		AliHLTUInt32_t* recordNum,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// dHLT internal channels data block is OK and returns true in that case.
	/// [in] \param block  The channels data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the data block is
	///      not valid.
	/// [out] \param recordNum  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the number of the channel structure that had a problem.
	///      The value 'recordNum[i]' will only contain a valid value if
	///      the corresponding 'reason[i]' contains one of:
	///        - kInvalidIdValue
	///        - kInvalidBusPatchId
	///        - kInvalidManuId
	///        - kInvalidChannelAddress
	///        - kInvalidSignal
	///        - kDataWordDifferent
	/// \note You can use RecordNumberWasSet(reason[i]) to check if 'recordNum[i]'
	///      was set and is valid or not.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason' and 'recordNum'. It will be filled
	///      with the number of items actually filled into the arrays upon exit
	///      from this method.
	/// \returns  true if there is no problem with the data and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = HeaderOk(block, reason, reasonCount);
	
	const AliHLTMUONChannelStruct* channel =
		reinterpret_cast<const AliHLTMUONChannelStruct*>(&block + 1);
	
	// Check if any cluster ID is duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = channel[i].fClusterId;
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (id == channel[j].fClusterId)
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateIDs;
					reasonCount++;
				}
				result = false;
			}
		}
	}

	// Check integrity of individual channel structures.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(channel[i], reason+reasonCount, filledCount))
		{
			// reasons filled in IntegrityOk, now we just need to adjust
			// reasonCount and fill the recordNum values.
			if (recordNum != NULL)
			{
				for (AliHLTUInt32_t n = 0; n < filledCount; n++)
					recordNum[reasonCount + n] = i;
			}
			reasonCount += filledCount;
			result = false;
		}
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONMansoTrackStruct& track,
		WhyNotValid* reason,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// Manso track structure is OK and returns true in that case.
	/// [in] \param track  The track structure to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the structure is
	///      not valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the structure and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// Check that the Manso track ID has a valid value.
	if (not (track.fId >= 0 or track.fId == -1))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidIdValue;
			reasonCount++;
		}
		result = false;
	}
	
	// Check that the corresponding trigger record ID has a valid value.
	if (not (track.fTrigRec >= 0 or track.fTrigRec == -1))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidTriggerIdValue;
			reasonCount++;
		}
		result = false;
	}
	
	// Make sure that the reserved bits in the fFlags field are set
	// to zero.
	if ((track.fFlags & 0x3FFFFFF0) != 0)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kReservedBitsNotZero;
			reasonCount++;
		}
		result = false;
	}

	// Make sure the sign is not invalid.
	if ((track.fFlags & 0xC0000000) == 0xC0000000)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kParticleSignBitsNotValid;
			reasonCount++;
		}
		result = false;
	}

	// Check that fHit[i] is nil if the corresponding bit in the
	// flags word is zero.
	const AliHLTMUONRecHitStruct& nilhit
		= AliHLTMUONConstants::NilRecHitStruct();
	if ( ((track.fFlags & 0x1) == 0 and track.fHit[0] != nilhit) or
	     ((track.fFlags & 0x2) == 0 and track.fHit[1] != nilhit) or
	     ((track.fFlags & 0x4) == 0 and track.fHit[2] != nilhit) or
	     ((track.fFlags & 0x8) == 0 and track.fHit[3] != nilhit)
	   )
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHitNotMarkedAsNil;
			reasonCount++;
		}
		result = false;
	}

	// Check that the chi squared value is valid
	if (not (track.fChi2 >= 0 or track.fChi2 == -1))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kChiSquareInvalid;
			reasonCount++;
		}
		result = false;
	}

	// Check that if chi squared is -1 then the momentum vector is zero.
	if (track.fChi2 == -1 and
	    not (track.fPx == 0 and track.fPy == 0 and track.fPz == 0)
	   )
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kMomentumVectorNotZero;
			reasonCount++;
		}
		result = false;
	}
	
	// Check the individual hits
	for (int i = 0; i < 4; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(track.fHit[i], reason + reasonCount, filledCount))
		{
			reasonCount += filledCount;
			result = false;
		}
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONMansoTracksBlockStruct& block,
		WhyNotValid* reason,
		AliHLTUInt32_t* recordNum,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// dHLT internal Manso track data block is OK and returns true in that case.
	/// [in] \param block  The Manso track data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the data block is
	///      not valid.
	/// [out] \param recordNum  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the number of the Manso track that had a problem.
	///      The value 'recordNum[i]' will only contain a valid value if
	///      the corresponding 'reason[i]' contains one of:
	///        - kInvalidIdValue
	///        - kInvalidTriggerIdValue
	///        - kReservedBitsNotZero
	///        - kParticleSignBitsNotValid
	///        - kHitNotMarkedAsNil
	///        - kChiSquareInvalid
	///        - kMomentumVectorNotZero
	///        - kInvalidDetElementNumber
	///        - kInvalidChamberNumber
	/// \note You can use RecordNumberWasSet(reason[i]) to check if 'recordNum[i]'
	///      was set and is valid or not.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason' and 'recordNum'. It will be filled
	///      with the number of items actually filled into the arrays upon exit
	///      from this method.
	/// \returns  true if there is no problem with the data and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = HeaderOk(block, reason, reasonCount);
	
	const AliHLTMUONMansoTrackStruct* track =
		reinterpret_cast<const AliHLTMUONMansoTrackStruct*>(&block + 1);
	
	// Check if any track ID is duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = track[i].fId;
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (id == track[j].fId)
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateIDs;
					reasonCount++;
				}
				result = false;
			}
		}
	}
	
	// Check that all the tracks have integrity.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(track[i], reason+reasonCount, filledCount))
		{
			// reasons filled in IntegrityOk, now we just need to adjust
			// reasonCount and fill the recordNum values.
			if (recordNum != NULL)
			{
				for (AliHLTUInt32_t n = 0; n < filledCount; n++)
					recordNum[reasonCount + n] = i;
			}
			reasonCount += filledCount;
			result = false;
		}
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONMansoCandidateStruct& candidate,
		WhyNotValid* reason,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// Manso track candidate structure is OK and returns true in that case.
	/// [in] \param track  The track candidate structure to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the structure is
	///      not valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the structure and false otherwise.
	
	// First check the integrity of the candidate track structure.
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = IntegrityOk(candidate.fTrack, reason, reasonCount);
	
	// Now check that the ROIs are reasonable.
	// The radius must be positive or -1 indicating computation error and
	// the corresponding hit in the track must be within the ROI.
	for (AliHLTUInt32_t i = 0; i < 4; i++)
	{
		if (not (candidate.fRoI[i].fRadius >= 0 or candidate.fRoI[i].fRadius == -1))
		{
			if (reason != NULL and reasonCount < maxCount)
			{
				reason[reasonCount] = kRoiRadiusInvalid;
				reasonCount++;
			}
			result = false;
		}
		
		// Check if the corresponding hit was even found in the track.
		if ( (candidate.fTrack.fFlags & (0x1 << i)) == 0 ) continue;
		
		double dx = candidate.fRoI[i].fX - candidate.fTrack.fHit[i].fX;
		double dy = candidate.fRoI[i].fY - candidate.fTrack.fHit[i].fY;
		double dz = candidate.fRoI[i].fZ - candidate.fTrack.fHit[i].fZ;
		double r = sqrt(dx*dx + dy*dy);
		// Check if the projected distance between ROI centre and hit is
		// bigger than the ROI radius. Also the difference between z
		// coordinates should not exceed 20 cm.
		if (r > candidate.fRoI[i].fRadius or fabs(dz) > 20.)
		{
			if (reason != NULL and reasonCount < maxCount)
			{
				reason[reasonCount] = kHitNotWithinRoi;
				reasonCount++;
			}
			result = false;
		}
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONMansoCandidatesBlockStruct& block,
		WhyNotValid* reason,
		AliHLTUInt32_t* recordNum,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// dHLT internal Manso candidates data block is OK and returns true in
	/// that case.
	/// [in] \param block  The Manso track candidate data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the data block is
	///      not valid.
	/// [out] \param recordNum  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the number of the track candidate that had a problem.
	///      The value 'recordNum[i]' will only contain a valid value if
	///      the corresponding 'reason[i]' contains one of:
	///        - kInvalidIdValue
	///        - kInvalidTriggerIdValue
	///        - kReservedBitsNotZero
	///        - kParticleSignBitsNotValid
	///        - kHitNotMarkedAsNil
	///        - kChiSquareInvalid
	///        - kRoiRadiusInvalid
	///        - kHitNotWithinRoi
	/// \note You can use RecordNumberWasSet(reason[i]) to check if 'recordNum[i]'
	///      was set and is valid or not.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason' and 'recordNum'. It will be filled
	///      with the number of items actually filled into the arrays upon exit
	///      from this method.
	/// \returns  true if there is no problem with the data and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = HeaderOk(block, reason, reasonCount);

	const AliHLTMUONMansoCandidateStruct* candidate =
		reinterpret_cast<const AliHLTMUONMansoCandidateStruct*>(&block + 1);
	
	// Check if any candidate track ID is duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = candidate[i].fTrack.fId;
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (id == candidate[j].fTrack.fId)
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateIDs;
					reasonCount++;
				}
				result = false;
			}
		}
	}
	
	// Check that all the track candidates have integrity.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(candidate[i], reason+reasonCount, filledCount))
		{
			// reasons filled in IntegrityOk, now we just need to adjust
			// reasonCount and fill the recordNum values.
			if (recordNum != NULL)
			{
				for (AliHLTUInt32_t n = 0; n < filledCount; n++)
					recordNum[reasonCount + n] = i;
			}
			reasonCount += filledCount;
			result = false;
		}
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONTrackDecisionStruct& decision,
		WhyNotValid* reason,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// single track trigger decision structure is OK and returns true in that case.
	/// [in] \param decision  The trigger decision structure to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the structure is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the structure and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	// The track ID value must be positive or -1.
	if (not (decision.fTrackId >= 0 or decision.fTrackId == -1))
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidTrackIdValue;
			reasonCount++;
		}
		result = false;
	}
	
	// Make sure that the reserved bits in the fTriggerBits field are set
	// to zero.
	if ((decision.fTriggerBits & 0xFFFFFFFC) != 0)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kReservedBitsNotZero;
			reasonCount++;
		}
		result = false;
	}
	
	// The pT should be -1 or a positive number.
	if (decision.fPt != -1. and decision.fPt < 0.)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kPtValueNotValid;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONSinglesDecisionBlockStruct& block,
		WhyNotValid* reason,
		AliHLTUInt32_t* recordNum,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// dHLT internal single track trigger decision data block is OK and returns
	/// true in that case.
	/// [in] \param block  The single track trigger decision data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the data block is
	///      not valid.
	/// [out] \param recordNum  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the number of the trigger decision that had a problem.
	///      The value 'recordNum[i]' will only contain a valid value if
	///      the corresponding 'reason[i]' contains one of:
	///        - kInvalidTrackIdValue
	///        - kReservedBitsNotZero
	///        - kPtValueNotValid
	/// \note You can use RecordNumberWasSet(reason[i]) to check if 'recordNum[i]'
	///      was set and is valid or not.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason' and 'recordNum'. It will be filled
	///      with the number of items actually filled into the arrays upon exit
	///      from this method.
	/// \returns  true if there is no problem with the data and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = HeaderOk(block, reason, reasonCount);
	
	const AliHLTMUONTrackDecisionStruct* decision =
		reinterpret_cast<const AliHLTMUONTrackDecisionStruct*>(&block + 1);

	// Check that there are no duplicate trigger entries.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = decision[i].fTrackId;
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (id == decision[j].fTrackId)
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateTriggers;
					reasonCount++;
				}
				result = false;
			}
		}
	}
	
	// Check that the trigger bits for each track have integrity.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(decision[i], reason+reasonCount, filledCount))
		{
			// Reasons filled in IntegrityOk, now we just need to adjust
			// reasonCount and fill the recordNum values.
			if (recordNum != NULL)
			{
				for (AliHLTUInt32_t n = 0; n < filledCount; n++)
					recordNum[reasonCount + n] = i;
			}
			reasonCount += filledCount;
			result = false;
		}
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONPairDecisionStruct& decision,
		WhyNotValid* reason,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// track pair trigger decision structure is OK and returns true in that case.
	/// [in] \param decision  The trigger decision structure to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the structure is not
	///      valid.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason'. It will be filled with the number
	///      of items actually filled into the reason array upon exit from this
	///      method.
	/// \returns  true if there is no problem with the structure and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	reasonCount = 0;
	bool result = true;
	
	//kInvalidTrackIdValue
	
	// The track IDs must have a positive value or -1.
	if (not (decision.fTrackAId >= 0 or decision.fTrackAId == -1) or
	    not (decision.fTrackBId >= 0 or decision.fTrackBId == -1)
	   )
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kInvalidTrackIdValue;
			reasonCount++;
		}
		result = false;
	}
	
	// Make sure that the reserved bits in the fTriggerBits field are set
	// to zero.
	if ((decision.fTriggerBits & 0xFFFFFF80) != 0)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kReservedBitsNotZero;
			reasonCount++;
		}
		result = false;
	}
	
	// Check that the track IDs are not the same.
	if (decision.fTrackAId == decision.fTrackBId)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kPairTrackIdsAreIdentical;
			reasonCount++;
		}
		result = false;
	}
	
	// The invariant mass should be -1 or a positive number.
	if (decision.fInvMass != -1. and decision.fInvMass < 0.)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kMassValueNotValid;
			reasonCount++;
		}
		result = false;
	}
	
	// Neither the high pt (hipt) or low pt (lopt) count bits can be > 2.
	AliHLTUInt8_t lowPtCount = (decision.fTriggerBits & 0x00000003);
	AliHLTUInt8_t highPtCount = (decision.fTriggerBits & 0x0000000C) >> 2;
	if (lowPtCount > 2)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kLowPtCountInvalid;
			reasonCount++;
		}
		result = false;
	}
	if (highPtCount > 2)
	{
		if (reason != NULL and reasonCount < maxCount)
		{
			reason[reasonCount] = kHighPtCountInvalid;
			reasonCount++;
		}
		result = false;
	}
	
	return result;
}


bool AliHLTMUONUtils::IntegrityOk(
		const AliHLTMUONPairsDecisionBlockStruct& block,
		WhyNotValid* reason,
		AliHLTUInt32_t* recordNum,
		AliHLTUInt32_t& reasonCount
	)
{
	/// This method is used to check more extensively if the integrity of the
	/// dHLT internal track pair trigger decision data block is OK and returns
	/// true in that case.
	/// [in] \param block  The track pair trigger decision data block to check.
	/// [out] \param reason  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the reason codes describing why the data block is
	///      not valid.
	/// [out] \param recordNum  If this is not NULL, then it is assumed to point
	///      to an array of at least 'reasonCount' number of elements. It will
	///      be filled with the number of the trigger decision that had a problem.
	///      The value 'recordNum[i]' will only contain a valid value if
	///      the corresponding 'reason[i]' contains one of:
	///        - kInvalidTrackIdValue
	///        - kReservedBitsNotZero
	///        - kPairTrackIdsAreIdentical
	///        - kMassValueNotValid
	///        - kLowPtCountInvalid
	///        - kHighPtCountInvalid
	/// \note You can use RecordNumberWasSet(reason[i]) to check if 'recordNum[i]'
	///      was set and is valid or not.
	/// [in/out] \param reasonCount  This should initially specify the size of
	///      the array pointed to by 'reason' and 'recordNum'. It will be filled
	///      with the number of items actually filled into the arrays upon exit
	///      from this method.
	/// \returns  true if there is no problem with the data and false otherwise.
	
	AliHLTUInt32_t maxCount = reasonCount;
	bool result = HeaderOk(block, reason, reasonCount);

	const AliHLTMUONPairDecisionStruct* decision =
		reinterpret_cast<const AliHLTMUONPairDecisionStruct*>(&block + 1);
	
	// Check that there are no duplicate trigger entries.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t ta = decision[i].fTrackAId;
		AliHLTInt32_t tb = decision[i].fTrackBId;
		for (AliHLTUInt32_t j = i+1; j < block.fHeader.fNrecords; j++)
		{
			if (ta == decision[j].fTrackAId and tb == decision[j].fTrackBId)
			{
				if (reason != NULL and reasonCount < maxCount)
				{
					reason[reasonCount] = kFoundDuplicateTriggers;
					reasonCount++;
				}
				result = false;
			}
		}
	}
	
	// Check that the trigger bits for each track pair have integrity.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t filledCount = maxCount - reasonCount;
		if (not IntegrityOk(decision[i], reason+reasonCount, filledCount))
		{
			// Reasons filled in IntegrityOk, now we just need to adjust
			// reasonCount and fill the recordNum values.
			if (recordNum != NULL)
			{
				for (AliHLTUInt32_t n = 0; n < filledCount; n++)
					recordNum[reasonCount + n] = i;
			}
			reasonCount += filledCount;
			result = false;
		}
	}
	
	return result;
}

