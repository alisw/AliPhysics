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

#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONClustersBlockStruct.h"
#include "AliHLTMUONChannelsBlockStruct.h"

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
			AliHLTMUONParticleSign sign, bool hitset[4]
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
			AliHLTMUONParticleSign sign, bool hitset[4]
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
	 * Methods used to check if the header information corresponds to the
	 * supposed type of the data block.
	 */
	static bool HeaderOk(const AliHLTMUONTriggerRecordsBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONRecHitsBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONClustersBlockStruct& block);
	static bool HeaderOk(const AliHLTMUONChannelsBlockStruct& block);

	/**
	 * Methods used to check extensively if the integrity of various types
	 * of data blocks are Ok and returns true in that case.
	 * These can be slow and should generally only be used for debugging.
	 */
	static bool IntegrityOk(const AliHLTMUONTriggerRecordsBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONRecHitsBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONClustersBlockStruct& block);
	static bool IntegrityOk(const AliHLTMUONChannelsBlockStruct& block);

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
//  DebugMsg("x = " << x << " and y = 0x" << std::hex << y );
//  // statements...
//
#ifdef DEBUG
#	include <ostream>
#	define DebugTrace(message) {std::cout << message << std::endl;}
#else // DEBUG
#	define DebugTrace(message)
#endif // DEBUG


#endif // ALIHLTMUONUTILS_H
