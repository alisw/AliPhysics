#ifndef ALIHLTMUONTRIGRECSDEBUGBLOCKSTRUCT_H
#define ALIHLTMUONTRIGRECSDEBUGBLOCKSTRUCT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/**
 * @file   AliHLTMUONTrigRecsDebugBlockStruct.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   19 May 2007
 * @brief  Definition of internal dimuon HLT data structures for storing debugging
 *         information about trigger records.
 *
 * The structures are defined with C linkage since C generally gives us more
 * binary compatibility between compilers.
 */

#include "AliHLTMUONDataTypes.h"
#include "AliMUONTriggerDDLDecoderEventHandler.h"

extern "C"
{

/**
 * Debug trigger record structures contain debug information about the trigger
 * records that were decoded from the dimuon hardware trigger electronics.
 */
struct AliHLTMUONTrigRecInfoStruct
{
	AliHLTInt32_t fTrigRecId; // The corresponding trigger record ID this
	                          // debug information corresponds to. -1 == invalid.

	AliHLTInt32_t fDetElemId[4];  // Detector ID numbers from AliRoot geometry
	                           // on which the hit clusters were found.
	                           // fDetElemId[0] is for chamber 11, fDetElemId[1]
	                           // for chamber 12 etc.

	// The parameters used for momentum estimation:
	AliHLTFloat32_t fZmiddle; // Z position in the middle of the magnetic field in centimetres.
	AliHLTFloat32_t fBl; // The integrated magnetic field in (T.m) tesla metres.
	
	AliMUONLocalInfoStruct fL0Struct;  // Copy of the L0 local trigger structure bits.
	AliMUONLocalInfoStruct fL0StructPrev;  // Copy of the previous L0 local trigger structure bits used. Can be zero if none used.
	AliMUONLocalInfoStruct fL0StructNext;  // Copy of the next L0 local trigger structure bits used. Can be zero if none used.
};

/**
 * AliHLTMUONTrigRecsDebugBlockStruct defines the format of the internal
 * data block for storing debugging information for trigger records.
 */
struct AliHLTMUONTrigRecsDebugBlockStruct
{
	AliHLTMUONDataBlockHeader fHeader; // Common data block header

	// Array of trigger records.
	//AliHLTMUONTrigRecInfoStruct fTrigRecInfo[/*fHeader.fNrecords*/];
};

} // extern "C"


/**
 * Stream operator for usage with std::ostream classes which prints the trigger
 * record debug information in the following format:
 *  {fTrigRecId = xx, fDetElemId = [aa, bb, cc, dd], fZmiddle = yy, fBl = zz}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTrigRecInfoStruct& info
	);

/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONTrigRecsDebugBlockStruct in the following format:
 *   {fHeader = xx, fTrigRecInfo[] = [{..}, {..}, ...]}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTrigRecsDebugBlockStruct& block
	);


bool operator == (
		const AliHLTMUONTrigRecInfoStruct& a,
		const AliHLTMUONTrigRecInfoStruct& b
	);

inline bool operator != (
		const AliHLTMUONTrigRecInfoStruct& a,
		const AliHLTMUONTrigRecInfoStruct& b
	)
{
	return not operator == (a, b);
}


bool operator == (
		const AliHLTMUONTrigRecsDebugBlockStruct& a,
		const AliHLTMUONTrigRecsDebugBlockStruct& b
	);

inline bool operator != (
		const AliHLTMUONTrigRecsDebugBlockStruct& a,
		const AliHLTMUONTrigRecsDebugBlockStruct& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONTRIGRECSDEBUGBLOCKSTRUCT_H
