#ifndef ALIHLTMUONTRACKSBLOCKSTRUCT_H
#define ALIHLTMUONTRACKSBLOCKSTRUCT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

// $Id: AliHLTMUONTracksBlockStruct.h 36627 2009-11-10 19:21:49Z aszostak $

///
/// @file   AliHLTMUONTracksBlockStruct.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   10 Feb 2010
/// @brief  Definition the internal dimuon HLT data block for full tracks.
///
/// The tracks data block is an internal dimuon HLT data block structure generated
/// by the tracker component to store found tracks.
/// The structures are defined with C linkage since C generally gives us more
/// binary compatibility between compilers.
///

#include "AliHLTMUONRecHitsBlockStruct.h"

extern "C"
{

/**
 * Track structure containing information about a track found by the tracker.
 */
struct AliHLTMUONTrackStruct
{
	AliHLTInt32_t fId; /// Each track should have an ID number unique for a given event. -1 == invalid.

	AliHLTInt32_t fTrigRec; /// The associated trigger record ID that was matched to the track.

	// The flags word constains the following bit fields (bit 31 is most
	// significant):
	// bits:  [31][30][29 -- 16][  15   ][  14   ]...[  1   ][  0   ]
	// field:  -   +   reserved  hst[15]  hst[14]     hst[1]  hst[0]
	// Where fields hst[i] indicates if fHit[i] has been filled/set.
	// Reserved bits should be set to zero.
	// Particle sign: if '-' is set then particle has minus sign.
	//                if '+' is set then particle has negative sign.
	// Either '+' or '-' should be set and if neither then the particle
	// sign is unknown.
	AliHLTUInt32_t fFlags;  /// Bit fields for the track structure.

	AliHLTFloat32_t fPx; /// Particle's momentum X component in GeV/c.
	AliHLTFloat32_t fPy; /// Particle's momentum Y component in GeV/c.
	AliHLTFloat32_t fPz; /// Particle's momentum Z component in GeV/c.
	
	AliHLTFloat32_t fInverseBendingMomentum; /// One over the momentum of the fitted track [GeV/c].
	AliHLTFloat32_t fThetaX; /// The slope of the fitted track in the non-bending plane.
	AliHLTFloat32_t fThetaY; /// The slope of the fitted track in the bending plane.
	AliHLTFloat32_t fX; /// Non-bending plane coordinate for the distance of closest approach (DCA) [cm].
	AliHLTFloat32_t fY; /// Bending plane coordinate for the DCA [cm].
	AliHLTFloat32_t fZ; /// Z coordinate for the DCA [cm].
	
	// Chi squared of the track fit.
	// If set to -1 then no fit was done and the momentum vector and DCA coordinate is invalid.
	AliHLTFloat32_t fChi2; /// The chi squared of the fit of fHit points to the track model.

	AliHLTMUONRecHitStruct fHit[16];  /// Particle hit coordinates found by the hit reconstruction stage.
};

/**
 * AliHLTMUONTracksBlockStruct defines the format of the internal tracks data block.
 */
struct AliHLTMUONTracksBlockStruct
{
	AliHLTMUONDataBlockHeader fHeader; // Common data block header.

	// Array of tracks.
	//AliHLTMUONTrackStruct fTrack[/*fHeader.fNrecords*/];
};

} // extern "C"


/**
 * Stream operator for usage with std::ostream classes which prints the track
 * information in the following format:
 *  {fId = xx, fTrigRec = yy, fFlags = 0xZZ, fPx = uu, fPy = vv, fPz = ww,
 *   fInverseBendingMomentum = pp, fThetaX = nn, fThetaY = mm, fX = kk, fY = ll, fZ = hh,
 *   fChi2 = ss, fHit[0] = {...}, fHit[1] = {...}, ... fHit[14] = {...}, fHit[15] = {...}}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTrackStruct& trigrec
	);

/**
 * Stream operator for usage with std::ostream classes which prints the
 * AliHLTMUONTracksBlockStruct in the following format:
 *   {fHeader = xx, fTrack[] = [{..}, {..}, ...]}
 */
std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONTracksBlockStruct& block
	);


bool operator == (
		const AliHLTMUONTrackStruct& a,
		const AliHLTMUONTrackStruct& b
	);

inline bool operator != (
		const AliHLTMUONTrackStruct& a,
		const AliHLTMUONTrackStruct& b
	)
{
	return not operator == (a, b);
}


bool operator == (
		const AliHLTMUONTracksBlockStruct& a,
		const AliHLTMUONTracksBlockStruct& b
	);

inline bool operator != (
		const AliHLTMUONTracksBlockStruct& a,
		const AliHLTMUONTracksBlockStruct& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONTRACKSBLOCKSTRUCT_H
