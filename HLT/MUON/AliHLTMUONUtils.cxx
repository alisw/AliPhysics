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


AliHLTUInt32_t AliHLTMUONUtils::PackTriggerRecordFlags(
		AliHLTMUONParticleSign sign, bool hitset[4]
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
		AliHLTUInt32_t flags, // [in]
		AliHLTMUONParticleSign& sign, // [out]
		bool hitset[4] // [out]
	)
{
	AliHLTUInt32_t signbits = flags & 0xC0000000;
	switch (signbits)
	{
	case 0x80000000: sign = kSignMinus;   break;
	case 0x40000000: sign = kSignPlus;    break;
	default:         sign = kSignUnknown; break;
	}
	hitset[0] = (flags & 0x1) == 1;
	hitset[1] = (flags & 0x2) == 1;
	hitset[2] = (flags & 0x4) == 1;
	hitset[3] = (flags & 0x8) == 1;
}


bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONTriggerRecordsBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kTriggerRecordsDataBlock) return false;
	// The block's record widths must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTriggerRecordStruct))
		return false;
	return true;
}

bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONTrigRecsDebugBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kTrigRecsDebugDataBlock) return false;
	// The block's record widths must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTrigRecInfoStruct))
		return false;
	return true;
}

bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONTriggerChannelsBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kTriggerChannelsDataBlock) return false;
	// The block's record widths must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONTriggerChannelStruct))
		return false;
	return true;
}

bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONRecHitsBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kRecHitsDataBlock) return false;
	// The block's record widths must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONRecHitStruct))
		return false;
	return true;
}

bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONClustersBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kClustersDataBlock) return false;
	// The block's record widths must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONClusterStruct))
		return false;
	return true;
}

bool AliHLTMUONUtils::HeaderOk(const AliHLTMUONChannelsBlockStruct& block)
{
	// The block must have the correct type.
	if (block.fHeader.fType != kChannelsDataBlock) return false;
	// The block's record widths must be the correct size.
	if (block.fHeader.fRecordWidth != sizeof(AliHLTMUONChannelStruct))
		return false;
	return true;
}


bool AliHLTMUONUtils::IntegrityOk(const AliHLTMUONTriggerRecordsBlockStruct& block)
{
	if (not HeaderOk(block)) return false;

	// Check if any ID is duplicated.
	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		AliHLTUInt32_t id = block.fTriggerRecord[i].fId;
		for (AliHLTUInt32_t j = i+1; i < block.fHeader.fNrecords; j++)
		{
			if (id == block.fTriggerRecord[j].fId)
				return false;
		}
	}

	for (AliHLTUInt32_t i = 0; i < block.fHeader.fNrecords; i++)
	{
		const AliHLTMUONTriggerRecordStruct& tr = block.fTriggerRecord[i];

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
		AliHLTUInt32_t id = block.fCluster[i].fId;
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
