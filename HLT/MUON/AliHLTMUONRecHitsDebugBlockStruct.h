#ifndef ALIHLTMUONRECHITSDEBUGBLOCKSTRUCT_H
#define ALIHLTMUONRECHITSDEBUGBLOCKSTRUCT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 * @file   AliHLTMUONRecHitsDebugBlockStruct.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Definition of internal dimuon HLT block structure containing
 *         debugging information about reconstructed hits.
 * 
 * The structures are defined with C linkage since C generally gives us more
 * binary compatibility between compilers.
 */

#include "AliHLTMUONDataTypes.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include <ostream>

extern "C"
{

/**
 * Debugging information about a cluster and its reconstructed hit.
 */
struct AliHLTMUONClusterInfoStruct
{
	AliHLTUInt32_t fClusterId; // Unique ID for the cluster. It must be at
	                           // least unique for any given data block.

	AliHLTMUONRecHitStruct fHit; // Corresponding reconstructed hit.

	AliHLTInt32_t fDetElemId;  // Detector ID number from AliRoot geometry
	                           // on which the cluster was found.

	AliHLTUInt32_t fNchannels; // Number of channels/pads in the cluster.
};

/**
 * Gives the fired channel/pad information which was considered during
 * hit reconstruction and correlated into a common cluster.
 */
struct AliHLTMUONChannelInfoStruct
{
	AliHLTUInt32_t fClusterId;   // ID corresponding to the cluster this
	                             // channels is part of.

	AliHLTUInt16_t fManu;        // The MANU address on electronics.
	AliHLTUInt16_t fChannel;     // The channel address on electronics.
	AliHLTUInt16_t fSignal;      // ADC value of signal.
	AliHLTUInt32_t fRawDataWord; // The raw data word as found in the DDL stream.
};

/**
 * AliHLTMUONRecHitsDebugBlockStruct defines the format of the internal data
 * block which contains debugging information about reconstructed hits.
 */
struct AliHLTMUONRecHitsDebugBlockStruct
{
	AliHLTUInt32_t fNclusters; // Number of cluster records in this block.
	
	// Offset in bytes to the start and end of the cluster data from the
	// beginning of this AliHLTMUONRecHitsDebugBlockStruct structure.
	AliHLTUInt32_t fClustersStartOffset;
	AliHLTUInt32_t fClustersEndOffset;

	AliHLTUInt32_t fNchannels; // Number of channel records in this block.
	
	// Offset in bytes to the start and end of the channel data from the
	// beginning of this AliHLTMUONRecHitsDebugBlockStruct structure.
	AliHLTUInt32_t fChannelsStartOffset;
	AliHLTUInt32_t fChannelsEndOffset;

	// Array of debugging information for reconstructed hits.
	AliHLTMUONClusterInfoStruct fCluster[/*fNclusters*/];

	// The array of channel information would follow here.
	//AliHLTMUONChannelInfoStruct fChannel[fNchannels];
};

} // extern "C"


/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONClusterInfoStruct in the following format:
 *   {fClusterId = xx, fHit = yy, fDetElemId = zz, fNchannels = ww}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONClusterInfoStruct& cluster
	);

/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONChannelInfoStruct in the following format:
 *   {fClusterId = xx, fManu = yy, fChannel = zz, fSignal = ww, fRawDataWord = 0xXXXXXXXX}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONChannelInfoStruct& channel
	);


/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONRecHitsDebugBlockStruct in the following format:
 *   {fNclusters = xx, fClustersStartOffset = yy, fClustersEndOffset = zz,
 *    fNchannels = uu, fChannelsStartOffset = vv, fChannelsEndOffset = ww,
 *    fCluster[] = [{..}, {..}, ...], fChannel[] = [{..}, {..}, ...]}
 */
std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONRecHitsDebugBlockStruct& block
	);


inline bool operator == (
		const AliHLTMUONClusterInfoStruct& a,
		const AliHLTMUONClusterInfoStruct& b
	)
{
	return	a.fClusterId == b.fClusterId and a.fHit == b.fHit and
		a.fDetElemId == b.fDetElemId and a.fNchannels == b.fNchannels;
}

inline bool operator != (
		const AliHLTMUONClusterInfoStruct& a,
		const AliHLTMUONClusterInfoStruct& b
	)
{
	return not operator == (a, b);
}


inline bool operator == (
		const AliHLTMUONChannelInfoStruct& a,
		const AliHLTMUONChannelInfoStruct& b
	)
{
	return	a.fClusterId == b.fClusterId and a.fManu == b.fManu and
		a.fChannel == b.fChannel and a.fSignal == b.fSignal and
		a.fRawDataWord == b.fRawDataWord;
}

inline bool operator != (
		const AliHLTMUONChannelInfoStruct& a,
		const AliHLTMUONChannelInfoStruct& b
	)
{
	return not operator == (a, b);
}


bool operator == (
		const AliHLTMUONRecHitsDebugBlockStruct& a,
		const AliHLTMUONRecHitsDebugBlockStruct& b
	);

inline bool operator != (
		const AliHLTMUONRecHitsDebugBlockStruct& a,
		const AliHLTMUONRecHitsDebugBlockStruct& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONRECHITSDEBUGBLOCKSTRUCT_H
