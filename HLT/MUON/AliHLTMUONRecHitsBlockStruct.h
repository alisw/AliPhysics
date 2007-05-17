#ifndef ALIHLTMUONRECHITSBLOCKSTRUCT_H
#define ALIHLTMUONRECHITSBLOCKSTRUCT_H
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
 * @file   AliHLTMUONRecHitsBlockStruct.h
 * @author Indranil Das <indra.das@saha.ac.in>,
 *         Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Definition of internal dimuon HLT reconstructed hit block structure.
 * 
 * The structures are defined with C linkage since C generally gives us more
 * binary compatibility between compilers.
 */

#include "AliHLTMUONDataTypes.h"

extern "C"
{

/**
 * A 3D reconstructed hit point structure.
 * These are used to store the hits reconstructed on the tracking or trigger
 * chambers.
 * Reconstructed hit values of (0, 0, 0) indicate an invalid or nil hit.
 */
struct AliHLTMUONRecHitStruct
{
	AliHLTFloat32_t fX; // X coordinate.
	AliHLTFloat32_t fY; // Y coordinate.
	AliHLTFloat32_t fZ; // Z coordinate.
};

/**
 * AliHLTMUONRecHitsBlockStruct defines the format of the internal
 * reconstructed hit data block.
 */
struct AliHLTMUONRecHitsBlockStruct
{
	AliHLTUInt32_t fNhits;                   // Number of hits in this block.
	AliHLTMUONRecHitStruct fHit[/*fNhits*/]; // Array of reconstructed hits.
};

} // extern "C"


/**
 * Stream operator for usage with std::ostream classes which prints the
 * reconstructed hit in the following format: {fX = xx, fY = yy, fZ = zz}
 */
inline std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONRecHitStruct& hit
	)
{
	stream	<< "{fX = " << hit.fX << ", fY = " << hit.fY << ", fZ = "
		<< hit.fZ << "}";
	return stream;
}

/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONRecHitsBlockStruct in the following format:
 *   {fNhits = xx, fHit[] = [{..}, {..}, ...]}
 */
inline std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONRecHitsBlockStruct& block
	)
{
	stream << "{fNhits = " << block.fNhits << ", fHit[] = [";
	if (block.fNhits > 0) stream << block.fHit[0];
	for (AliHLTUInt32_t i = 1; i < block.fNhits; i++)
		stream << ", " << block.fHit[i];
	stream << "]}";
	return stream;
}


inline bool operator == (
		const AliHLTMUONRecHitStruct& a, const AliHLTMUONRecHitStruct& b
	)
{
	return a.fX == b.fX and a.fY == b.fY and a.fZ == b.fZ;
}

inline bool operator != (
		const AliHLTMUONRecHitStruct& a, const AliHLTMUONRecHitStruct& b
	)
{
	return not operator == (a, b);
}


inline bool operator == (
		const AliHLTMUONRecHitsBlockStruct& a, const AliHLTMUONRecHitsBlockStruct& b
	)
{
	// First check if the blocks have the same number of hits. If they do then
	// check if every hit is the same. In either case if we find a difference
	// return false.
	if (a.fNhits != b.fNhits) return false;
	for (AliHLTUInt32_t i = 0; i < a.fNhits; i++)
		if (a.fHit[i] != b.fHit[i]) return false;
	return true;
}

inline bool operator != (
		const AliHLTMUONRecHitsBlockStruct& a, const AliHLTMUONRecHitsBlockStruct& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONRECHITSBLOCKSTRUCT_H
