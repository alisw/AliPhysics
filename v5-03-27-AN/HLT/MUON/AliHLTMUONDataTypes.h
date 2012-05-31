#ifndef ALIHLTMUONDATATYPES_H
#define ALIHLTMUONDATATYPES_H
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

// $Id$

/**
 * @file   AliHLTMUONDataTypes.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   18 May 2007
 * @brief  Declaration of basic data types used in HLT/MUON module.
 *
 * The types and structs are defined with C linkage since C generally gives us
 * more binary compatibility between compilers.
 */

#include "AliHLTDataTypes.h"
#include <ostream>

extern "C"
{

/**
 * The common internal dimuon HLT data block header.
 * These headers help to identify the data block when it is written to disk and
 * helps to check the integrity of the data blocks in the system.
 */
struct AliHLTMUONDataBlockHeader
{
	AliHLTUInt16_t fType; // The type of the data block. Must contain a value
	                      // defined by AliHLTMUONDataBlockType.
	AliHLTUInt16_t fRecordWidth; // The number of bytes each record uses. 
	AliHLTUInt32_t fNrecords; // Number of records in this data block.
};

/**
 * The field structure of a single row of the AliHLTMUONHitReconstructor component's
 * lookup table.
 */
struct AliHLTMUONHitRecoLutRow
{
	AliHLTInt32_t fDetElemId;  // The detector element ID of the pad.
	AliHLTInt32_t fIX, fIY;    // The X,Y number of the pad.
	AliHLTFloat32_t fRealX, fRealY, fRealZ;  // The real coordinate of the pad.
	AliHLTFloat32_t fHalfPadSize; // half padsize in X for bending and Y for nonbending
	AliHLTFloat32_t fPadSizeXY; // padsize in Y for bending plane and X for nonbending
	AliHLTInt32_t fPlane;  // The plane and PCB zone ID numbers.
	AliHLTFloat32_t fPed, fSigma, fA0, fA1; // Calibration values
	AliHLTInt32_t fThres, fSat; //Calibration values
};

/**
 * The field structure of a single row of the AliHLTMUONTriggerReconstructor component's
 * lookup table.
 */
struct AliHLTMUONTriggerRecoLutRow
{
	AliHLTUInt32_t fIdFlags;  /// The chamber and detector element identifier packed in AliHLTMUONRecHitStruct::fFlags format.
	AliHLTFloat32_t fX;  // Global X coordinate of channel.
	AliHLTFloat32_t fY;  // Global Y coordinate of channel.
	AliHLTFloat32_t fZ;  // Global Z coordinate of channel.
};

/**
 * The lookup table structure for the AliHLTMUONTriggerReconstructor component.
 * The LUT is used for translating from channel addresses to geometrical positions
 * and other relevant information of the strips in the trigger chambers.
 */
struct AliHLTMUONTriggerRecoLookupTable
{
	// The dimentions of the LUT are as follows:
	// [trigger crate ID][local board ID][chamber][cathode - X/Y][bit set in bit pattern]
	// - The trigger crate ID comes form the regional headers in the DDL payload.
	// - Local board ID numbers come from the local trigger structures in the
	//   DDL payload.
	// - The chamber is for chambers 11 to 14 coded as [0..3], 0 for chamber 11,
	//   1 for chamber 12 etc.
	// - The cathode 0 is for the X trigger strips (bending plane) and 1 is
	//   for the Y strips (non-bending plane).
	// - The "bit set in pattern" indicates the bit number that was set in
	//   the strip patterns found in the local structures of the DDL payload.
	AliHLTMUONTriggerRecoLutRow fRow[16][16][4][2][16];
};

} // extern "C"


/**
 * The sign/charge of a particle.
 */
enum AliHLTMUONParticleSign
{
	kSignMinus   = -1,
	kSignUnknown = 0,
	kSignPlus    = 1
};

/**
 * The chamber names of the dimuon spectrometer.
 */
enum AliHLTMUONChamberName
{
	kUnknownChamber = -1,
	kChamber1 = 0,
	kChamber2 = 1,
	kChamber3 = 2,
	kChamber4 = 3,
	kChamber5 = 4,
	kChamber6 = 5,
	kChamber7 = 6,
	kChamber8 = 7,
	kChamber9 = 8,
	kChamber10 = 9,
	kChamber11 = 10,
	kChamber12 = 11,
	kChamber13 = 12,
	kChamber14 = 13
};

/**
 * The internal data block type codes.
 */
enum AliHLTMUONDataBlockType
{
	kUnknownDataBlock = 0,
	kTriggerRecordsDataBlock = 1000,
	kTrigRecsDebugDataBlock = 1001,
	kRecHitsDataBlock = 2000,
	kClustersDataBlock = 2001,
	kChannelsDataBlock = 2002,
	kMansoTracksDataBlock = 3000,
	kMansoCandidatesDataBlock = 3001,
	kTracksDataBlock = 3100,
	kSinglesDecisionDataBlock = 4000,
	kPairsDecisionDataBlock = 4001
};


/**
 * Stream operator to be able to print AliHLTMUONParticleSign with human
 * readable names to some stream object like cout (standard output).
 */
inline std::ostream& operator << (std::ostream& stream, AliHLTMUONParticleSign sign)
{
	switch (sign)
	{
	case kSignMinus:   stream << "kSignMinus";   break;
	case kSignUnknown: stream << "kSignUnknown"; break;
	case kSignPlus:    stream << "kSignPlus";    break;
	default:           stream << "INVALID";
	}
	return stream;
}

/**
 * Stream operator to be able to print AliHLTMUONChamberName with human
 * readable names to some stream object like cout (standard output).
 */
inline std::ostream& operator << (std::ostream& stream, AliHLTMUONChamberName chamber)
{
	switch (chamber)
	{
	case kUnknownChamber:  stream << "kUnknownChamber";  break;
	case kChamber1:        stream << "kChamber1";        break;
	case kChamber2:        stream << "kChamber2";        break;
	case kChamber3:        stream << "kChamber3";        break;
	case kChamber4:        stream << "kChamber4";        break;
	case kChamber5:        stream << "kChamber5";        break;
	case kChamber6:        stream << "kChamber6";        break;
	case kChamber7:        stream << "kChamber7";        break;
	case kChamber8:        stream << "kChamber8";        break;
	case kChamber9:        stream << "kChamber9";        break;
	case kChamber10:       stream << "kChamber10";       break;
	case kChamber11:       stream << "kChamber11";       break;
	case kChamber12:       stream << "kChamber12";       break;
	case kChamber13:       stream << "kChamber13";       break;
	case kChamber14:       stream << "kChamber14";       break;
	default:               stream << "INVALID";
	}
	return stream;
}

/**
 * Stream operator to be able to print AliHLTMUONDataBlockType with human
 * readable names to some stream object like cout (standard output).
 */
inline std::ostream& operator << (std::ostream& stream, AliHLTMUONDataBlockType type)
{
	switch (type)
	{
	case kUnknownDataBlock:         stream << "kUnknownDataBlock";         break;
	case kTriggerRecordsDataBlock:  stream << "kTriggerRecordsDataBlock";  break;
	case kTrigRecsDebugDataBlock:   stream << "kTrigRecsDebugDataBlock";   break;
	case kRecHitsDataBlock:         stream << "kRecHitsDataBlock";         break;
	case kClustersDataBlock:        stream << "kClustersDataBlock";        break;
	case kChannelsDataBlock:        stream << "kChannelsDataBlock";        break;
	case kMansoTracksDataBlock:     stream << "kMansoTracksDataBlock";     break;
	case kMansoCandidatesDataBlock: stream << "kMansoCandidatesDataBlock"; break;
	case kTracksDataBlock:          stream << "kTracksDataBlock";          break;
	case kSinglesDecisionDataBlock: stream << "kSinglesDecisionDataBlock"; break;
	case kPairsDecisionDataBlock:   stream << "kPairsDecisionDataBlock";   break;
	default:                        stream << "INVALID";
	}
	return stream;
}

/**
 * Stream operator for usage with std::ostream classes which prints the common
 * data block header in the following format:
 *  {fType = xx, fRecordWidth = yy, fNrecords = zz}
 */
inline std::ostream& operator << (
		std::ostream& stream, const AliHLTMUONDataBlockHeader& header
	)
{
	stream	<< "{fType = " << AliHLTMUONDataBlockType(header.fType)
		<< ", fRecordWidth = " << header.fRecordWidth
		<< ", fNrecords = " << header.fNrecords << "}";
	return stream;
}


inline bool operator == (
		const AliHLTMUONDataBlockHeader& a, const AliHLTMUONDataBlockHeader& b
	)
{
	return a.fType == b.fType and a.fRecordWidth == b.fRecordWidth
		and a.fNrecords == b.fNrecords;
}

inline bool operator != (
		const AliHLTMUONDataBlockHeader& a, const AliHLTMUONDataBlockHeader& b
	)
{
	return not operator == (a, b);
}

#endif // ALIHLTMUONDATATYPES_H
