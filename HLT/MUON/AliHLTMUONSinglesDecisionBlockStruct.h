#ifndef ALIHLTMUONSINGLESDECISIONBLOCKSTRUCT_H
#define ALIHLTMUONSINGLESDECISIONBLOCKSTRUCT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/**
 * @file   AliHLTMUONSinglesDecisionBlockStruct.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   21 May 2007
 * @brief  Definition of internal dimuon HLT trigger decision data structure
 *         containing decision information for single tracks.
 *
 * The structures are defined with C linkage since C generally gives us more
 * binary compatibility between compilers.
 */

#include "AliHLTMUONDataTypes.h"

extern "C"
{

/**
 * Structure contains information about the trigger decision for a single track.
 */
struct AliHLTMUONTrackDecisionStruct
{
	AliHLTInt32_t fTrackId; // The ID number of the track we are referring to.
	                        // This could also be a trigger record ID.
	
	// The trigger bits have the following meaning:
	// bit:  [31 --- 2][  1 ][  0 ]
	// field: reserved  hipt  lopt
	// Reserved bits should be set to zero.
	// hipt == passed high pt cut. lopt == passed low pt cut.
	AliHLTUInt32_t fTriggerBits;
	
	AliHLTFloat32_t fPt; // The calculated transverse momentum of the track in (GeV/c).
};

/**
 * AliHLTMUONSinglesDecisionBlockStruct defines the format of the internal
 * dimuon HLT decision data block for individual tracks.
 */
struct AliHLTMUONSinglesDecisionBlockStruct
{
	AliHLTMUONDataBlockHeader fHeader; // Common data block header.
	AliHLTUInt32_t fNlowPt;  // Number of low pt triggers.
	AliHLTUInt32_t fNhighPt; // Number of high pt triggers.

	// Array of decisions for individual tracks.
	//AliHLTMUONTrackDecisionStruct fDecision[/*fHeader.fNrecords*/];
};

} // extern "C"


/**
 * Stream operator for usage with std::ostream classes which prints the
 * trigger decision information for individual tracks in the following format:
 *  {fTrackId = xx, fTriggerBits = 0xYY}
 */
inline std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTrackDecisionStruct& trig
	)
{
	std::ios::fmtflags oldflags = stream.flags();
	stream	<< "{fTrackId = " << trig.fTrackId << ", fTriggerBits = "
		<< std::showbase << std::hex << trig.fTriggerBits << std::dec
		<< ", fPt = " << trig.fPt
		<< "}";
	stream.flags(oldflags);
	return stream;
}

/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONSinglesDecisionBlockStruct in the following format:
 *   {fHeader = xx, fNlowPt = yy, fNhighPt = zz, fDecision[] = [{..}, {..}, ...]}
 */
std::ostream& operator << (
		std::ostream& stream,
		const AliHLTMUONSinglesDecisionBlockStruct& block
	);


inline bool operator == (
		const AliHLTMUONTrackDecisionStruct& a,
		const AliHLTMUONTrackDecisionStruct& b
	)
{
	return a.fTrackId == b.fTrackId and a.fTriggerBits == b.fTriggerBits
		and a.fPt == b.fPt;
}

inline bool operator != (
		const AliHLTMUONTrackDecisionStruct& a,
		const AliHLTMUONTrackDecisionStruct& b
	)
{
	return not operator == (a, b);
}


bool operator == (
		const AliHLTMUONSinglesDecisionBlockStruct& a,
		const AliHLTMUONSinglesDecisionBlockStruct& b
	);

inline bool operator != (
		const AliHLTMUONSinglesDecisionBlockStruct& a,
		const AliHLTMUONSinglesDecisionBlockStruct& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONSINGLESDECISIONBLOCKSTRUCT_H
