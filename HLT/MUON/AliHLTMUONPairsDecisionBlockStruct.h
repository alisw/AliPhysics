#ifndef ALIHLTMUONPAIRSDECISIONBLOCKSTRUCT_H
#define ALIHLTMUONPAIRSDECISIONBLOCKSTRUCT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/**
 * @file   AliHLTMUONPairsDecisionBlockStruct.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   21 May 2007
 * @brief  Definition of internal dimuon HLT trigger decision data structure containing decision information for pairs of tracks.
 *
 * The structures are defined with C linkage since C generally gives us more
 * binary compatibility between compilers.
 */

#include "AliHLTMUONDataTypes.h"

extern "C"
{

/**
 * Structure contains information about the trigger decision for a pair of tracks.
 */
struct AliHLTMUONPairDecisionStruct
{
	// The ID numbers of the tracks we are referring to.
	// These could also be trigger record ID numbers.
	AliHLTInt32_t fTrackAId;
	AliHLTInt32_t fTrackBId;
	
	// The trigger bits have the following meaning:
	// bit:  [31 --- 7][ 6 ][ 5 ][   4  ][3--2][1--0]
	// field: reserved  hiM  loM  unlike  hipt  lopt
	// Reserved bits should be set to zero.
	// hiM == Pair passed the high invariant mass cut.
	// loM == Pair passed the low invariant mass cut.
	// unlike == is the pair an unlike sign pair.
	// hipt: 0 if neither track has pt > high pt cut, 1 if either has
	// and 2 if both have. A value of 3 is invalid.
	// lopt: 0 if neither track has pt > lo pt cut, 1 if either has
	// and 2 if both have. A value of 3 is invalid.
	AliHLTUInt32_t fTriggerBits;
	
	AliHLTFloat32_t fInvMass; // Invariant mass of the track pair assuming
	                          // they are both muons in (GeV/c^2)
};

/**
 * AliHLTMUONPairsDecisionBlockStruct defines the format of the internal
 * dimuon HLT decision data block for pairs of tracks.
 */
struct AliHLTMUONPairsDecisionBlockStruct
{
	AliHLTMUONDataBlockHeader fHeader; // Common data block header.
	
	// Number of unlike sign pt triggers for both tracks having any pt.
	AliHLTUInt32_t fNunlikeAnyPt;
	
	// Number of unlike sign low pt triggers where both tracks have
	// pt > low cut.
	AliHLTUInt32_t fNunlikeLowPt;
	
	// Number of unlike sign high pt triggers where both tracks have
	// pt > high cut.
	AliHLTUInt32_t fNunlikeHighPt;
	
	// Number of like sign pt triggers where both tracks have any pt.
	AliHLTUInt32_t fNlikeAnyPt;
	
	// Number of like sign low pt triggers where both tracks have
	// pt > low cut.
	AliHLTUInt32_t fNlikeLowPt;
	
	// Number of like sign high pt triggers where both tracks have
	// pt > high cut.
	AliHLTUInt32_t fNlikeHighPt;
	
	// Number of pairs that have invariant mass > low mass cut, any pt
	// and unlike sign.
	AliHLTUInt32_t fNmassAny;
	
	// Number of pairs that have invariant mass > low mass cut,
	// pt > low pt cut and unlike sign.
	AliHLTUInt32_t fNmassLow;
	
	// Number of pairs that have invariant mass > high mass cut,
	// pt > high pt cut and unlike sign.
	AliHLTUInt32_t fNmassHigh;

	// Array of decisions for track pairs.
	//AliHLTMUONPairDecisionStruct fDecision[/*fHeader.fNrecords*/];
};

} // extern "C"


/**
 * Stream operator for usage with std::ostream classes which prints the
 * trigger decision information for pairs of tracks in the following format:
 *  {fTrackAId = xx, fTrackBId = yy, fTriggerBits = 0xZZ, fInvMass = ww}
 */
inline std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONPairDecisionStruct& trig
	)
{
	std::ios::fmtflags oldflags = stream.flags();
	stream	<< "{fTrackAId = " << trig.fTrackAId
		<< ", fTrackBId = " << trig.fTrackBId << ", fTriggerBits = "
		<< std::showbase << std::hex << trig.fTriggerBits << std::dec
		<< ", fInvMass = " << trig.fInvMass << "}";
	stream.flags(oldflags);
	return stream;
}

/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONPairsDecisionBlockStruct in the following format:
 *   {fHeader = aa, fNunlikeLowPt = bb, fNunlikeHighPt = cc, fNlikeLowPt = dd,
 *    fNlikeHighPt = ee, fNmassAny = ff, fNmassLow = gg, fNmassHigh = hh
 *    fDecision[] = [{..}, {..}, ...]}
 */
std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONPairsDecisionBlockStruct& block
	);


inline bool operator == (
		const AliHLTMUONPairDecisionStruct& a,
		const AliHLTMUONPairDecisionStruct& b
	)
{
	return a.fTrackAId == b.fTrackAId and a.fTrackBId == b.fTrackBId
		and a.fTriggerBits == b.fTriggerBits
		and a.fTriggerBits == b.fTriggerBits;
}

inline bool operator != (
		const AliHLTMUONPairDecisionStruct& a,
		const AliHLTMUONPairDecisionStruct& b
	)
{
	return not operator == (a, b);
}


bool operator == (
		const AliHLTMUONPairsDecisionBlockStruct& a,
		const AliHLTMUONPairsDecisionBlockStruct& b
	);

inline bool operator != (
		const AliHLTMUONPairsDecisionBlockStruct& a,
		const AliHLTMUONPairsDecisionBlockStruct& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONPAIRSDECISIONBLOCKSTRUCT_H
