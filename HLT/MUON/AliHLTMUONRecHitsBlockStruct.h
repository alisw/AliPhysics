#ifndef ALIHLTMUONRECHITSBLOCKSTRUCT_H
#define ALIHLTMUONRECHITSBLOCKSTRUCT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 * @file   AliHLTMUONRecHitsBlockStruct.h
 * @author Indranil Das <indra.das@saha.ac.in>,
 *         Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Definition of internal dimuon HLT reconstructed hit data block structure.
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
	AliHLTMUONDataBlockHeader fHeader; // Common data block header

	// Array of reconstructed hits.
	//AliHLTMUONRecHitStruct fHit[/*fHeader.fNrecords*/];
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
 *   {fHeader = xx, fHit[] = [{..}, {..}, ...]}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONRecHitsBlockStruct& block
	);


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


bool operator == (
		const AliHLTMUONRecHitsBlockStruct& a,
		const AliHLTMUONRecHitsBlockStruct& b
	);

inline bool operator != (
		const AliHLTMUONRecHitsBlockStruct& a,
		const AliHLTMUONRecHitsBlockStruct& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONRECHITSBLOCKSTRUCT_H
