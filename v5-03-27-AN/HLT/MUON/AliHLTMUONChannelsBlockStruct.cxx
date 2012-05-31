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

// $Id$

///
/// @file   AliHLTMUONChannelsBlockStruct.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   18 May 2007
/// @brief  Implementation of the stream and comparison operators.
///
/// The channels blocks are internal dimuon HLT block structures containing
/// debugging information about channels that belong to reconstructed
/// hit clusters.
///

#include "AliHLTMUONChannelsBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <cassert>

std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONChannelStruct& channel
	)
{
	std::ios::fmtflags oldflags = stream.flags();
	stream	<< "{fClusterId = " << channel.fClusterId
		<< ", fBusPatch = " << channel.fBusPatch
		<< ", fManu = " << channel.fManu
		<< ", fChannelAddress = " << channel.fChannelAddress
		<< ", fSignal = " << channel.fSignal
		<< ", fRawDataWord = " << std::showbase << std::hex
		<< channel.fRawDataWord << std::dec << "}";
	stream.flags(oldflags);
	return stream;
}


std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONChannelsBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );

	const AliHLTMUONChannelStruct* channel =
		reinterpret_cast<const AliHLTMUONChannelStruct*>(&block + 1);
	stream 	<< "{fHeader = " << block.fHeader << ", fChannel[] = [";
	if (block.fHeader.fNrecords > 0) stream << channel[0];
	for (AliHLTUInt32_t i = 1; i < block.fHeader.fNrecords; i++)
		stream << ", " << channel[i];
	stream << "]}";
	return stream;
}


bool operator == (
		const AliHLTMUONChannelsBlockStruct& a,
		const AliHLTMUONChannelsBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );
	
	const AliHLTMUONChannelStruct* channelA =
		reinterpret_cast<const AliHLTMUONChannelStruct*>(&a + 1);
	const AliHLTMUONChannelStruct* channelB =
		reinterpret_cast<const AliHLTMUONChannelStruct*>(&b + 1);

	// First check if the blocks have the same header. If they do then check
	// if every channel is the same. In either case if we find a difference
	// return false.
	if (a.fHeader != b.fHeader) return false;
	for (AliHLTUInt32_t i = 0; i < a.fHeader.fNrecords; i++)
		if (channelA[i] != channelB[i]) return false;
	return true;
}
