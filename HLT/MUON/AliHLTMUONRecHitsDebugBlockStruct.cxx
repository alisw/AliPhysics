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
 * @file   AliHLTMUONRecHitsDebugBlockStruct.cxx
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Implementation of the stream and comparison operators.
 */

#include "AliHLTMUONRecHitsDebugBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include <assert.h>


std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONClusterInfoStruct& cluster
	)
{
	stream	<< "{fClusterId = " << cluster.fClusterId
		<< ", fHit = " << cluster.fHit
		<< ", fDetElemId = " << cluster.fDetElemId
		<< ", fNchannels = " << cluster.fNchannels << "}";
	return stream;
}


std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONChannelInfoStruct& channel
	)
{
	stream	<< "{fClusterId = " << channel.fClusterId
		<< ", fManu = " << channel.fManu
		<< ", fChannel = " << channel.fChannel
		<< ", fSignal = " << channel.fSignal
		<< ", fRawDataWord = " << std::showbase << std::hex
		<< channel.fRawDataWord << std::dec << "}";
	return stream;
}


std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONRecHitsDebugBlockStruct& block
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(block) );

	const AliHLTMUONChannelInfoStruct* channel =
		reinterpret_cast<const AliHLTMUONChannelInfoStruct*>(
			& block.fCluster[block.fNclusters]
		);

	// Now write the block header information.
	stream	<< "{fNclusters = " << block.fNclusters
		<< ", fClustersStartOffset = " << block.fClustersStartOffset
		<< ", fClustersEndOffset = " << block.fClustersEndOffset
		<< ", fNchannels = " << block.fNchannels
		<< ", fChannelsStartOffset = " << block.fChannelsStartOffset
		<< ", fChannelsEndOffset = " << block.fChannelsEndOffset
		<< ", fCluster[] = ";

	// And finally write the arrays.
	if (block.fNclusters > 0) stream << block.fCluster[0];
	for (AliHLTUInt32_t i = 1; i < block.fNclusters; i++)
		stream << ", " << block.fCluster[i];

	stream << "], fChannels[] = ";
	if (block.fNchannels > 0) stream << channel[0];
	for (AliHLTUInt32_t i = 1; i < block.fNchannels; i++)
		stream << ", " << channel[i];
	stream << "]}";

	return stream;
}


bool operator == (
		const AliHLTMUONRecHitsDebugBlockStruct& a,
		const AliHLTMUONRecHitsDebugBlockStruct& b
	)
{
	assert( AliHLTMUONUtils::IntegrityOk(a) );
	assert( AliHLTMUONUtils::IntegrityOk(b) );

	const AliHLTMUONChannelInfoStruct* channelA =
		reinterpret_cast<const AliHLTMUONChannelInfoStruct*>(
			& a.fCluster[a.fNclusters]
		);
	const AliHLTMUONChannelInfoStruct* channelB =
		reinterpret_cast<const AliHLTMUONChannelInfoStruct*>(
			& b.fCluster[b.fNclusters]
		);

	// First check if the blocks have the same header fields. If they do then
	// check the cluster and channel arrays are the same. In either case if we
	// find a difference return false.
	if (
		a.fNclusters != b.fNclusters
		or a.fClustersStartOffset != b.fClustersStartOffset
		or a.fClustersEndOffset != b.fClustersEndOffset
		or a.fNchannels != b.fNchannels
		or a.fChannelsStartOffset != b.fChannelsStartOffset
		or a.fChannelsEndOffset != b.fChannelsEndOffset
	   )
		return false;

	for (AliHLTUInt32_t i = 0; i < a.fNclusters; i++)
		if (a.fCluster[i] != b.fCluster[i]) return false;
	for (AliHLTUInt32_t i = 0; i < a.fNchannels; i++)
		if (channelA[i] != channelB[i]) return false;

	return true;
}
