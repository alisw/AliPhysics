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
 * @file   AliHLTMUONTriggerChannelsBlockStruct.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of the stream and comparison operators.
 */

#include "AliHLTMUONTriggerChannelsBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <cassert>

std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTriggerChannelStruct& channel
	)
{
	stream	<< "{fTrigRecId = " << channel.fTrigRecId
		<< ", fChamber = " << channel.fChamber
		<< ", fSignal = " << channel.fSignal
		<< ", fRawDataWord = " << std::showbase << std::hex
		<< channel.fRawDataWord << std::dec << "}";
	return stream;
}


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONTriggerChannelsBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );

	stream 	<< "{fHeader = " << block.fHeader << ", fChannel[] = [";
	if (block.fHeader.fNrecords > 0) stream << block.fChannel[0];
	for (AliHLTUInt32_t i = 1; i < block.fHeader.fNrecords; i++)
		stream << ", " << block.fChannel[i];
	stream << "]}";
	return stream;
}


bool operator == (
		const AliHLTMUONTriggerChannelsBlockStruct& a,
		const AliHLTMUONTriggerChannelsBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );

	// First check if the blocks have the same header. If they do then check
	// if every channel is the same. In either case if we find a difference
	// return false.
	if (a.fHeader != b.fHeader) return false;
	for (AliHLTUInt32_t i = 0; i < a.fHeader.fNrecords; i++)
		if (a.fChannel[i] != b.fChannel[i]) return false;
	return true;
}
