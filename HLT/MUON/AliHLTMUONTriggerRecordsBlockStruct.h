#ifndef ALIHLTMUONTRIGGERRECORDSBLOCKSTRUCT_H
#define ALIHLTMUONTRIGGERRECORDSBLOCKSTRUCT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 * @file   AliHLTMUONTriggerRecordsBlockStruct.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Definition of internal dimuon HLT trigger record data block structure.
 * 
 * The structures are defined with C linkage since C generally gives us more
 * binary compatibility between compilers.
 */

#include "AliHLTMUONRecHitsBlockStruct.h"

extern "C"
{

/**
 * Trigger record structures containing information decoded from the dimuon
 * spectrometer's hardware trigger electronics.
 */
struct AliHLTMUONTriggerRecordStruct
{
	AliHLTInt32_t fId; // Each trigger record should have an ID number unique
	                   // for a given event. -1 == invalid.

	// The flags word constains the following bit fields (bit 31 is most
	// significant):
	// bits:  [31][30][29 --- 4][  3   ][  2   ][  1   ][  0   ]
	// field:  -   +   reserved  hst[3]  hst[2]  hst[1]  hst[0]
	// Where fields hst[i] indicates if fHit[i] has been filled/set.
	// Reserved bits should be set to zero.
	// Particle sign: if '-' is set then particle has minus sign.
	//                if '+' is set then particle has negative sign.
	// Either '+' or '-' should be set and if neither then the particle
	// sign is unknown.
	AliHLTUInt32_t fFlags;

	AliHLTFloat32_t fPx; // Particle momentum X component in GeV/c.
	AliHLTFloat32_t fPy; // Particle momentum Y component in GeV/c.
	AliHLTFloat32_t fPz; // Particle momentum Z component in GeV/c.

	// Particle hit coordinates on trigger chambers 11 to 14.
	AliHLTMUONRecHitStruct fHit[4];
};

/**
 * AliHLTMUONTriggerRecordsBlockStruct defines the format of the internal
 * trigger records data block.
 */
struct AliHLTMUONTriggerRecordsBlockStruct
{
	AliHLTMUONDataBlockHeader fHeader; // Common data block header

	// Array of trigger records.
#ifndef __SUNPRO_CC
	AliHLTMUONTriggerRecordStruct fTriggerRecord[/*fHeader.fNrecords*/];
#else
	AliHLTMUONTriggerRecordStruct fTriggerRecord[1];
#endif
};

} // extern "C"


/**
 * Stream operator for usage with std::ostream classes which prints the trigger
 * record in the following format:
 *  {fId = xx, fFlags = 0xYY, fPx = uu, fPy = vv, fPz = ww, fHit[0] = {...},
 *   fHit[1] = {...}, fHit[2] = {...}, fHit[3] = {...}}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTriggerRecordStruct& trigrec
	);

/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONTriggerRecordsBlockStruct in the following format:
 *   {fHeader = xx, fTriggerRecord[] = [{..}, {..}, ...]}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTriggerRecordsBlockStruct& block
	);


inline bool operator == (
		const AliHLTMUONTriggerRecordStruct& a,
		const AliHLTMUONTriggerRecordStruct& b
	)
{
	return a.fId == b.fId and a.fFlags == b.fFlags
		and a.fPx == b.fPx and a.fPy == b.fPy and a.fPz == b.fPz
		and a.fHit[0] == b.fHit[0] and a.fHit[1] == b.fHit[1]
		and a.fHit[2] == b.fHit[2] and a.fHit[3] == b.fHit[3];
}

inline bool operator != (
		const AliHLTMUONTriggerRecordStruct& a,
		const AliHLTMUONTriggerRecordStruct& b
	)
{
	return not operator == (a, b);
}


bool operator == (
		const AliHLTMUONTriggerRecordsBlockStruct& a,
		const AliHLTMUONTriggerRecordsBlockStruct& b
	);

inline bool operator != (
		const AliHLTMUONTriggerRecordsBlockStruct& a,
		const AliHLTMUONTriggerRecordsBlockStruct& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONTRIGGERRECORDSBLOCKSTRUCT_H
