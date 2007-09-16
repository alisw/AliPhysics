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
 * @file   AliHLTMUONUtils.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of AliHLTMUONUtils utility routines.
 */

#include "AliHLTMUONUtils.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONTrigRecsDebugBlockStruct.h"
#include "AliHLTMUONTriggerChannelsBlockStruct.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONClustersBlockStruct.h"
#include "AliHLTMUONChannelsBlockStruct.h"
#include "AliHLTMUONMansoTracksBlockStruct.h"
#include "AliHLTMUONMansoCandidatesBlockStruct.h"
#include "AliHLTMUONSinglesDecisionBlockStruct.h"
#include "AliHLTMUONPairsDecisionBlockStruct.h"
#include <cassert>


AliHLTUInt32_t AliHLTMUONUtils::PackTriggerRecordFlags(
		AliHLTMUONParticleSign sign, const bool hitset[4]
	)
{
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


AliHLTUInt32_t AliHLTMUONUtils::PackTrackDecisionBits(bool highPt, bool lowPt)
{
	return (highPt ? 0x2 : 0) | (lowPt ? 0x1 : 0);
}


void AliHLTMUONUtils::UnpackTrackDecisionBits(
		AliHLTUInt32_t bits, bool& highPt, bool& lowPt
	)
{
	lowPt  = (bits & 0x1) == 1;
	highPt = (bits & 0x2) == 1;
}


AliHLTUInt32_t AliHLTMUONUtils::PackPairDecisionBits(
		bool highMass, bool lowMass, bool unlike,
		AliHLTUInt8_t highPtCount, AliHLTUInt8_t lowPtCount
	)
{
	assert( highPtCount + lowPtCount <= 2 );
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
	highMass = (bits & 0x40) == 1;
	lowMass  = (bits & 0x20) == 1;
	unlike   = (bits & 0x10) == 1;
	highPtCount = (bits & 0xC) >> 2;
	lowPtCount = bits & 0x3;
}


AliHLTUInt32_t AliHLTMUONUtils::PackSpecBits(
		const bool ddl[22]
	)
{
	// Pack the bits into the following format:
	//   bit:   [        31 - 22        ][     21     ][     20     ][  19 - 0 ]
	//   field: [ reserved, set to zero ][ TRGDDL2816 ][ TRGDDL2817 ][ TRKDDLS ]
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
	// Perform the inverse operation of PackSpecBits.
	for (int i = 0; i < 22; i++)
		ddl[i] = ((bits >> i) & 0x1) == 1;
}


bool AliHLTMUONUtils::HeaderOk(
		const AliHLTMUONTriggerRecordsBlockStruct& block,
		WhyNotValid* reason
	)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kTriggerRecordsDataBlock)
	{
		if (reason != NULL) *reason = kHeaderContainsWrongType;
		return false;
	}
	
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTriggerRecordStruct))
	{
		if (reason != NULL) *reason = kHeaderContainsWrongRecordWidth;
		return false;
	}
	
	return true;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONTrigRecsDebugBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kTrigRecsDebugDataBlock) return false;
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTrigRecInfoStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONTriggerChannelsBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kTriggerChannelsDataBlock) return false;
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTriggerChannelStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONRecHitsBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kRecHitsDataBlock) return false;
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONRecHitStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONClustersBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kClustersDataBlock) return false;
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONClusterStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONChannelsBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kChannelsDataBlock) return false;
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONChannelStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONMansoTracksBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kMansoTracksDataBlock) return false;
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONMansoTrackStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONMansoCandidatesBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kMansoCandidatesDataBlock) return false;
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONMansoCandidateStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONSinglesDecisionBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kSinglesDecisionDataBlock) return false;
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTrackDecisionStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONPairsDecisionBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kPairsDecisionDataBlock) return false;
	// The block's record width must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONPairDecisionStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONTriggerRecordStruct& tr)
{
	// Make sure that the reserved bits in the fFlags field are set
	// to zero.
	if ((tr.fFlags & 0x3FFFFFF0) != 0) return false;

	// Make sure the sign is not invalid.
	if ((tr.fFlags & 0xC0000000) == 3) return false;

	// Check that fHit[i] is nil if the corresponding bit in the
	// flags word is zero.
	const AliHLTMUONRecHitStruct& nilhit
		= AliHLTMUONConstants::NilRecHitStruct();
	if ((tr.fFlags & 0x1) == 0 and tr.fHit[0] != nilhit) return false;
	if ((tr.fFlags & 0x2) == 0 and tr.fHit[1] != nilhit) return false;
	if ((tr.fFlags & 0x4) == 0 and tr.fHit[2] != nilhit) return false;
	if ((tr.fFlags & 0x8) == 0 and tr.fHit[3] != nilhit) return false;

	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONTriggerRecordsBlockStruct& block)
{
	if (not HeaderOk(block)) return false;

	// Check if any ID is duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = block.fTriggerRecord[i].fId;
		for (AliHLTUInt32_t j = i+1; i < block.fHeader.fNrecords; j++)
		{
			if (id == block.fTriggerRecord[j].fId)
				return false;
		}
	}

	// Check integrity of individual trigger records.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		if (not IntegrityOk(block.fTriggerRecord[i])) return false;
	}

	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONTrigRecsDebugBlockStruct& block)
{
	if (not HeaderOk(block)) return false;
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONTriggerChannelsBlockStruct& block)
{
	if (not HeaderOk(block)) return false;
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONRecHitsBlockStruct& block)
{
	if (not HeaderOk(block)) return false;
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONClustersBlockStruct& block)
{
	if (not HeaderOk(block)) return false;

	// Check if any ID is duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = block.fCluster[i].fId;
		for (AliHLTUInt32_t j = i+1; i < block.fHeader.fNrecords; j++)
		{
			if (id == block.fCluster[j].fId)
				return false;
		}
	}
	
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONChannelsBlockStruct& block)
{
	if (not HeaderOk(block)) return false;
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONMansoTrackStruct& track)
{
	// Make sure that the reserved bits in the fFlags field are set
	// to zero.
	if ((track.fFlags & 0x3FFFFFF0) != 0) return false;

	// Make sure the sign is not invalid.
	if ((track.fFlags & 0xC0000000) == 0xC0000000) return false;

	// Check that fHit[i] is nil if the corresponding bit in the
	// flags word is zero.
	const AliHLTMUONRecHitStruct& nilhit
		= AliHLTMUONConstants::NilRecHitStruct();
	if ((track.fFlags & 0x1) == 0 and track.fHit[0] != nilhit) return false;
	if ((track.fFlags & 0x2) == 0 and track.fHit[1] != nilhit) return false;
	if ((track.fFlags & 0x4) == 0 and track.fHit[2] != nilhit) return false;
	if ((track.fFlags & 0x8) == 0 and track.fHit[3] != nilhit) return false;
	
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONMansoTracksBlockStruct& block)
{
	if (not HeaderOk(block)) return false;

	// Check if any ID is duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTInt32_t id = block.fTrack[i].fId;
		for (AliHLTUInt32_t j = i+1; i < block.fHeader.fNrecords; j++)
		{
			if (id == block.fTrack[j].fId)
				return false;
		}
	}

	// Check that the tracks have integrity.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		if (not IntegrityOk(block.fTrack[i])) return false;
	}

	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONMansoCandidatesBlockStruct& block)
{
	if (not HeaderOk(block)) return false;

	// Check that the tracks have integrity.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		if (not IntegrityOk(block.fCandidate[i].fTrack)) return false;
	}
	
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONTrackDecisionStruct& decision)
{
	// Make sure that the reserved bits in the fTriggerBits field are set
	// to zero.
	if ((decision.fTriggerBits & 0xFFFFFFFC) != 0) return false;
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONSinglesDecisionBlockStruct& block)
{
	if (not HeaderOk(block)) return false;

	// Check that the trigger bits for each track have integrity.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		if (not IntegrityOk(block.fDecision[i])) return false;
	}
	
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONPairDecisionStruct& decision)
{
	// Make sure that the reserved bits in the fTriggerBits field are set
	// to zero.
	if ((decision.fTriggerBits & 0xFFFFFF80) != 0) return false;
	
	// The high mass or low mass bits can only be set if unlike bit is set.
	if ((decision.fTriggerBits & 0x00000010) == 0
	    and (decision.fTriggerBits & 0x00000060) != 0
	   )
		return false;
	
	// Neither the high pt (hipt) or low pt (lopt) count bits can be > 2.
	// And the sum must not be > 2.
	AliHLTUInt8_t lowPtCount = (decision.fTriggerBits & 0x00000003);
	AliHLTUInt8_t highPtCount = (decision.fTriggerBits & 0x0000000C) >> 2;
	if (lowPtCount + highPtCount > 2) return false;
	
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONPairsDecisionBlockStruct& block)
{
	if (not HeaderOk(block)) return false;

	// Check that the trigger bits for each track pair have integrity.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		if (not IntegrityOk(block.fDecision[i])) return false;
	}
	
	return true;
}
