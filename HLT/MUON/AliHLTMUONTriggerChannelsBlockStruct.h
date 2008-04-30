#ifndef ALIHLTMUONTRIGGERCHANNELSBLOCKSTRUCT_H
#define ALIHLTMUONTRIGGERCHANNELSBLOCKSTRUCT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 * @file   AliHLTMUONTriggerChannelsBlockStruct.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Definition of internal dimuon HLT block structure containing
 *         debugging information about channels that belong to trigger records.
 * 
 * The structures are defined with C linkage since C generally gives us more
 * binary compatibility between compilers.
 */

#include "AliHLTMUONDataTypes.h"
#include <ostream>

extern "C"
{

/**
 * Gives the fired channel/strip information which was considered during
 * hit reconstruction on trigger chambers and correlated into a common cluster.
 */
struct AliHLTMUONTriggerChannelStruct
{
	AliHLTInt32_t fTrigRecId; // ID corresponding to the trigger records this
	                          // channel is part of. -1 == invalid.
	
	AliHLTUInt16_t fChamber; // The chamber this hit corresponds to.
	                         // In the range [11..14].
	
	AliHLTUInt16_t fSignal;      // ADC value of signal.
	AliHLTUInt32_t fRawDataWord; // The raw data word as found in the DDL stream.
};

/**
 * AliHLTMUONTriggerChannelsBlockStruct defines the format of the internal
 * channel data block corresponding to clusters on the trigger chambers.
 */
struct AliHLTMUONTriggerChannelsBlockStruct
{
	AliHLTMUONDataBlockHeader fHeader; // Common data block header

	// Array of trigger channels/strips.
	//AliHLTMUONTriggerChannelStruct fChannel[/*fHeader.fNrecords*/];
};

} // extern "C"


/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONTriggerChannelStruct in the following format:
 *   {fTrigRecId = xx, fChamber = yy, fSignal = zz, fRawDataWord = 0xXXXXXXXX}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTriggerChannelStruct& channel
	);

/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONTriggerChannelsBlockStruct in the following format:
 *   {fHeader = xx, fChannel[] = [{..}, {..}, ...]}
 */
std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONTriggerChannelsBlockStruct& block
	);


inline bool operator == (
		const AliHLTMUONTriggerChannelStruct& a,
		const AliHLTMUONTriggerChannelStruct& b
	)
{
	return a.fTrigRecId == b.fTrigRecId and a.fChamber == b.fChamber and
		a.fSignal == b.fSignal and a.fRawDataWord == b.fRawDataWord;
}

inline bool operator != (
		const AliHLTMUONTriggerChannelStruct& a,
		const AliHLTMUONTriggerChannelStruct& b
	)
{
	return not operator == (a, b);
}


bool operator == (
		const AliHLTMUONTriggerChannelsBlockStruct& a,
		const AliHLTMUONTriggerChannelsBlockStruct& b
	);

inline bool operator != (
		const AliHLTMUONTriggerChannelsBlockStruct& a,
		const AliHLTMUONTriggerChannelsBlockStruct& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONTRIGGERCHANNELSBLOCKSTRUCT_H
