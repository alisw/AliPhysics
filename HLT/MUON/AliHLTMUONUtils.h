#ifndef ALIHLTMUONUTILS_H
#define ALIHLTMUONUTILS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**
 * @file   AliHLTMUONUtils.h
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   
 * @brief  Class containing various dimuon HLT utility routines and macros.
 */

#include "AliHLTMUONDataTypes.h"

// Forward declare structures.
extern "C" {
struct AliHLTMUONTriggerRecordStruct;
struct AliHLTMUONTriggerRecordsBlockStruct;
struct AliHLTMUONTrigRecsDebugBlockStruct;
struct AliHLTMUONTriggerChannelsBlockStruct;
struct AliHLTMUONRecHitsBlockStruct;
struct AliHLTMUONClustersBlockStruct;
struct AliHLTMUONChannelsBlockStruct;
struct AliHLTMUONMansoTrackStruct;
struct AliHLTMUONMansoTracksBlockStruct;
struct AliHLTMUONMansoCandidatesBlockStruct;
struct AliHLTMUONTrackDecisionStruct;
struct AliHLTMUONSinglesDecisionBlockStruct;
struct AliHLTMUONPairDecisionStruct;
struct AliHLTMUONPairsDecisionBlockStruct;
} // extern "C"

/**
 * AliHLTMUONUtils contains arbitrary utility methods to be used in various
 * parts of the dimuon HLT system.
 * These include methods to perform basic sanity checks on the integrity of
 * data blocks.
 */
class AliHLTMUONUtils
{
public:
	/**
	 * This packs the given parameters into the bits of a word appropriate
	 * for AliHLTMUONTriggerRecordStruct::fFlags.
	 * @param sign    The particle sign.
	 * @param hitset  Flags to indicate if the corresponding fHits[i] elements
	 *                was set/filled.
	 * @return  Returns the 32 bit packed word.
	 */
	static AliHLTUInt32_t PackTriggerRecordFlags(
			AliHLTMUONParticleSign sign, const bool hitset[4]
		);

	/**
	 * This unpacks the AliHLTMUONTriggerRecordStruct::fFlags bits into
	 * its component fields.
	 * @param flags  The flags from an AliHLTMUONTriggerRecordStruct structure.
	 * @param sign    Sets this to the particle sign.
	 * @param hitset  Sets the array elements to indicate if the corresponding
	 *                fHits[i] element was set/filled.
	 */
	static void UnpackTriggerRecordFlags(
			AliHLTUInt32_t flags, // [in]
			AliHLTMUONParticleSign& sign, // [out]
			bool hitset[4] // [out]
		);

	/**
	 * This packs the given parameters into the bits of a word appropriate
	 * for AliHLTMUONMansoTrackStruct::fFlags.
	 * @param sign    The particle sign.
	 * @param hitset  Flags to indicate if the corresponding fHits[i] elements
	 *                was set/filled.
	 * @return  Returns the 32 bit packed word.
	 */
	static AliHLTUInt32_t PackMansoTrackFlags(
			AliHLTMUONParticleSign sign, const bool hitset[4]
		)
	{
		return PackTriggerRecordFlags(sign, hitset);
	}

	/**
	 * This unpacks the AliHLTMUONMansoTrackStruct::fFlags bits into
	 * its component fields.
	 * @param flags  The flags from an AliHLTMUONMansoTrackStruct structure.
	 * @param sign    Sets this to the particle sign.
	 * @param hitset  Sets the array elements to indicate if the corresponding
	 *                fHits[i] element was set/filled.
	 */
	static void UnpackMansoTrackFlags(
			AliHLTUInt32_t flags, // [in]
			AliHLTMUONParticleSign& sign, // [out]
			bool hitset[4] // [out]
		)
	{
		UnpackTriggerRecordFlags(flags, sign, hitset);
	}
	
	/**
	 * This packs the given parameters into the bits of a word appropriate
	 * for AliHLTMUONTrackDecisionStruct::fTriggerBits.
	 * @param highPt  Has the track passed the high pt cut.
	 * @param lowPt   Has the track passed the low pt cut.
	 * @return  Returns the 32 bit packed word.
	 */
	static AliHLTUInt32_t PackTrackDecisionBits(bool highPt, bool lowPt);
	
	/**
	 * This unpacks the AliHLTMUONTrackDecisionStruct::fTriggerBits bits into
	 * its component fields.
	 * @param bits  The trigger bits from an AliHLTMUONTrackDecisionStruct
	 *              structure.
	 * @param highPt Sets this to the value of the high pt cut bit.
	 * @param lowPt  Sets this to the value of the low pt cut bit.
	 */
	static void UnpackTrackDecisionBits(
			AliHLTUInt32_t bits, // [in]
			bool& highPt, // [out]
			bool& lowPt // [out]
		);
	
	/**
	 * This packs the given parameters into the bits of a word appropriate
	 * for AliHLTMUONPairDecisionStruct::fTriggerBits.
	 *
	 * @param highMass Has the track pair passed the high invariant mass cut.
	 * @param lowMass  Has the track pair passed the low invariant mass cut.
	 * @param unlike   Does the track pair have unlike signs.
	 * @param highPtCount The number of tracks that passed the high pt cut
	 *                    in the pair.
	 * @param lowPtCount  The number of tracks that passed the low pt cut
	 *                    in the pair.
	 * @return  Returns the 32 bit packed word.
	 *
	 * Note: Must have highPtCount + lowPtCount <= 2 and unlike == true if
	 * highMass or lowMass is true.
	 */
	static AliHLTUInt32_t PackPairDecisionBits(
			bool highMass, bool lowMass, bool unlike,
			AliHLTUInt8_t highPtCount, AliHLTUInt8_t lowPtCount
		);
	
	/**
	 * This unpacks the AliHLTMUONPairDecisionStruct::fTriggerBits bits into
	 * its component fields.
	 * @param bits  The trigger bits from an AliHLTMUONPairDecisionStruct
	 *              structure.
	 * @param highMass Sets this to the value of the high invariant mass cut bit.
	 * @param lowMass  Sets this to the value of the low invariant mass cut bit.
	 * @param unlike   Sets this if the pair is unlike sign.
	 * @param highPtCount Sets this to the high pt count bits.
	 * @param lowPtCount  Sets this to the low pt count bits.
	 */
	static void UnpackPairDecisionBits(
			AliHLTUInt32_t bits, // [in]
			bool& highMass, // [out]
			bool& lowMass, // [out]
			bool& unlike, // [out]
			AliHLTUInt8_t& highPtCount, // [out]
			AliHLTUInt8_t& lowPtCount // [out]
		);
	
	/**
	 * This packs the given parameters into the 32bit Pub/Sub specification
	 * word in the data block descriptor.
	 *
	 * @param ddl  The list of DDLs forming part of the readout. ddl[0]
	 *             indicates DDL number 2560, ddl[1] is for DDL 2561 and so
	 *             on up to ddl[19]. ddl[20] and ddl[21] will be for the
	 *             trigger DDLs 2816 and 2817 respectively.
	 * @return  Returns the 32 bit packed specification word.
	 */
	static AliHLTUInt32_t PackSpecBits(
			const bool ddl[22]
		);
	
	/**
	 * This unpacks the AliHLTMUONPairDecisionStruct::fTriggerBits bits into
	 * its component fields.
	 * @param bits  The Pub/Sub specification word from a data block descriptor.
	 * @param ddl  The output list of DDLs forming part of the readout. ddl[0]
	 *             indicates DDL number 2560, ddl[1] is for DDL 2561 and so
	 *             on up to ddl[19]. ddl[20] and ddl[21] will be for the
	 *             trigger DDLs 2816 and 2817 respectively.
	 */
	static void UnpackSpecBits(
			AliHLTUInt32_t bits, // [in]
			bool ddl[22] // [out]
		);

	/**
	 * These codes indicate the reason why a data block failed its
	 * validity check.
	 */
	enum WhyNotValid
	{
		kNoReason,
		kHeaderContainsWrongType,
		kHeaderContainsWrongRecordWidth,
	};

	/**
	 * Methods used to check if the header information corresponds to the
	 * supposed type of the data block.
	 * If the 'reason' parameter is not NULL then these methods will fill the
	 * memory pointed to by reason with a code describing of why the header
	 * is not valid, if and only if a problem is found with the data.
	 */
	static bool HeaderOk(const AliHLTMUONTriggerRecordsBlockStruct& block, WhyNotValid* reason = NULL);
	static bool HeaderOk(const AliHLTMUONTrigRecsDebugBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONTriggerChannelsBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONRecHitsBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONClustersBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONChannelsBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONMansoTracksBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONMansoCandidatesBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONSinglesDecisionBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONPairsDecisionBlockStruct& block);

	/**
	 * Methods used to check more extensively if the integrity of various
	 * types of data blocks are Ok and returns true in that case.
	 * These can be slow and should generally only be used for debugging.
	 */
	static bool IntegrityOk(const AliHLTMUONTriggerRecordStruct& tr);
	static bool IntegrityOk(const AliHLTMUONTriggerRecordsBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONTrigRecsDebugBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONTriggerChannelsBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONRecHitsBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONClustersBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONChannelsBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONMansoTrackStruct& track);
	static bool IntegrityOk(const AliHLTMUONMansoTracksBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONMansoCandidatesBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONTrackDecisionStruct& decision);
	static bool IntegrityOk(const AliHLTMUONSinglesDecisionBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONPairDecisionStruct& decision);
	static bool IntegrityOk(const AliHLTMUONPairsDecisionBlockStruct& block);

private:
	// Should never have to create or destroy this object.
	AliHLTMUONUtils();
	~AliHLTMUONUtils();
};

//_____________________________________________________________________________

// Since c++ is missing a finally "keyword" we define one. Its usage is identical
// to a try..finally statement in Java etc.. however, since it is officialy a macro
// one must use the ( ) brackets instead of { }
// If the compiler supports __finally use it otherwise make our own.
#if defined(__BORLANDC__)
#	define finally(str) __finally{str}
#else
#	define finally(code) \
		catch(...) \
		{ \
			code \
			throw; \
		}; \
		code
#endif // __BORLANDC__

// If we do not already have them, then define logical operators that are easier
// to read. 'and' = &&, 'or' = ||, 'not' = !
#if ! defined(__GNUC__) && ! defined(__CINT__)
// TODO: Should use iso646.h
#	define and &&
#	define or ||
#	define not !
#endif // __GNUC__ | __CINT__


// Here we define the DebugTrace(message) macro for easy embedding of debug
// information into the dimuon HLT code. Its usage is meant to be for generating
// traces of the program which are only useful during full scale debugging.
// Log messages should use the standard HLT logging mechanisms.
// The output is only generated in programs compiled with the DEBUG directive
// defined. Here is a usage example:
//
//  // statements...
//  DebugTrace("some debug information.");
//  // statements...
//
// One can also use C++ ostream operators and manipulators like so:
//
//  // statements...
//  int x, y;
//  DebugTrace("x = " << x << " and y = 0x" << std::hex << y );
//  // statements...
//
#ifdef DEBUG
#	include <iostream>
#	define DebugTrace(message) {std::cout << message << std::endl;}
#else // DEBUG
#	define DebugTrace(message)
#endif // DEBUG


#endif // ALIHLTMUONUTILS_H
