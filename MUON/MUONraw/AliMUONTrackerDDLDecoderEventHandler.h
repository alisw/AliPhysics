#ifndef ALIMUONTRACKERDDLDECODEREVENTHANDLER_H
#define ALIMUONTRACKERDDLDECODEREVENTHANDLER_H
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
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

///
/// \file   AliMUONTrackerDDLDecoderEventHandler.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   28-11-2007
/// \brief  Implementation of a high performance DDL decoder event handler 
/// for the muon tracking stations.
///

#include <cassert>
#include <ostream>
#include <Rtypes.h>


// We use C binding for the structures because C is more uniform with its application
// binary interface (ABI) between compilers.
extern "C"
{

// The following structures are the headers found in the DDL payload from the
// muon tracking chambers. The specification is defined in ALICE-INT-2005-012
// (https://edms.cern.ch/file/591904/1/ALICE-INT-2005-012.pdf)

/// The block header structure of the Tracker DDL payload.
struct AliMUONBlockHeaderStruct
{
	UInt_t     fDataKey;        ///< Data key word for CRT header 
	UInt_t     fTotalLength;    ///< total length of block structure (w/o padding word)
	UInt_t     fLength;         ///< length of raw data
	UInt_t     fDSPId;          ///< DSP id
	UInt_t     fL0Trigger;      ///< L0 trigger word
	UInt_t     fMiniEventId;    ///< Bunch Crossing for mini-event id (see TDR chapter 8)
	UInt_t     fEventId1;       ///< Event Id in bunch crossing
	UInt_t     fEventId2;       ///< Event Id in orbit number
};

/// The DSP header structure of the Tracker DDL payload.
struct AliMUONDSPHeaderStruct
{
	UInt_t     fDataKey;          ///< Data key word for FRT header 
	UInt_t     fTotalLength;      ///< total length of block structure
	UInt_t     fLength;           ///< length of raw data
	UInt_t     fDSPId;            ///< DSP id
	UInt_t     fBlkL1ATrigger;    ///< L1 accept in Block Structure (CRT)
	UInt_t     fMiniEventId;      ///< Mini Event Id in bunch crossing 
	UInt_t     fL1ATrigger;       ///< Number of L1 accept in DSP Structure (FRT)
	UInt_t     fL1RTrigger;       ///< Number of L1 reject in DSP Structure (FRT)
	UInt_t     fPaddingWord;      ///< padding dummy word for 64 bits transfer
	UInt_t     fErrorWord;        ///< Error word
};

/// The bus patch header structure of the Tracker DDL payload.
struct AliMUONBusPatchHeaderStruct
{
	UInt_t     fDataKey;       ///< Data key word for bus patch header 
	UInt_t     fTotalLength;   ///< total length of bus patch structure
	UInt_t     fLength;        ///< length of raw data
	UInt_t     fBusPatchId;    ///< bus patch id
};

} // extern "C"


/// \ingroup raw
/// \class AliMUONTrackerDDLDecoderEventHandler
/// \brief Callback event handler class for the AliMUONTrackerDDLDecoder.
///
/// This class is the base class defining what methods the event handler for the
/// high performance decoder should have. This handler actually does nothing.
/// The user of this decoder will have to derive from this class a custom event
/// handler that actually does something within the callback methods OnNewBusPatch,
/// OnData, OnError etc...
///
/// \author Artur Szostak <artursz@iafrica.com>

class AliMUONTrackerDDLDecoderEventHandler
{
public:

	/// The only reason for a virtual destructor is to make -Weffc++ shutup.
	/// This should not really be here since we do not need or use virtual methods.
	virtual ~AliMUONTrackerDDLDecoderEventHandler() {}

	/// All the possible error codes for the parsing.
	enum ErrorCode
	{
		kNoError = 0,                  /// Decoding was successful.
		// Offset our error codes to stay clear of any common codes in AliMUONRawStreamTracker:
		kBufferTooBig = 10,            /// The DDL raw data is larger than indicated by the headers; extra bytes are probably just garbage.
		kTooManyBlocks = 11,           /// Too many block structures found.
		kTooManyDSPs = 12,             /// Too many DSP structures found in the block.
		kTooManyBusPatches = 13,       /// Too many bus patch structures found in the DSP structure.
		kNoBlockHeader = 14,           /// Missing a block header.
		kBadBlockKey = 15,             /// The block header key word does not contain the correct value.
		kBadBlockLength = 16,          /// The block length field points past the end of the raw data size.
		kBadBlockTotalLength = 17,     /// The total block length field points past the end of the raw data size.
		kBlockLengthMismatch = 18,     /// The block length and total length fields do not correspond. One or both of these values is incorrect.
		kNoDSPHeader = 19,             /// Missing a DSP header.
		kBadDSPKey = 20,               /// The DSP header key word does not contain the correct value.
		kBadDSPLength = 21,            /// The DSP structure length field points past the end of the block structure.
		kBadDSPTotalLength = 22,       /// The total DSP structure length field points past the end of the block structure.
		kDSPLengthMismatch = 23,       /// The DSP structure length and total length fields do not correspond. One or both of these values is incorrect.
		kNoBusPatchHeader = 24,        /// Missing a bus patch header.
		kBadBusPatchKey = 25,          /// The bus patch header key word does not contain the correct value.
		kBadBusPatchLength = 26,       /// The bus patch length field points past the end of the DSP structure.
		kBadBusPatchTotalLength = 27,  /// The total bus patch length field points past the end of the DSP structure.
		kBusPatchLengthMismatch = 28,  /// The bus patch length and total length fields do not correspond. One or both of these values is incorrect.
		kNoDDLTrailerWords = 29,       /// No end of DDL markers found in the trailer words.
		kTooFewDDLTrailerWords = 30,   /// Only one end of DDL marker trailer word found but expected two.
		kUnknownDspError = 31,         /// The DSP error code is non-zero but of an unrecognised format.
		kTokenLost = 32,               /// The DSP contains a token lost error code that can affect the deadtime.
		// match up error codes with AliMUONRawStreamTracker:
		kGlitchFound = 1,              /// Found a glitch. This means a 1 byte word has been randomly inserted into the raw data by mistake.
		kBadPaddingWord = 2,           /// The padding word does not contain the correct value.
		kParityError = 3               /// Found a parity error in the data word.
	};

	// The following methods should be overridden for specific processing to
	// take place in your event handler.

	/// The OnNewBuffer method will be called whenever a new buffer containing
	/// a DDL payload is about to be processed.
	/// The default behaviour of this method is to do nothing.
	/// - param const void*  The pointer to the start of the memory buffer storing
	///                the DDL payload.
	/// - param UInt_t The size in bytes of the memory buffer.
	void OnNewBuffer(const void* /*buffer*/, UInt_t /*bufferSize*/) {}
	
	/// The OnEndOfBuffer method will be called whenever the buffer containing
	/// a DDL payload has been processed. For each OnNewBuffer method call a
	/// symmetric call to OnEndOfBuffer is made at the end of processing (after
	/// the last call to OnData)
	/// The default behaviour of this method is to do nothing.
	/// - param const void*  The pointer to the start of the memory buffer storing
	///                the DDL payload.
	/// - param UInt_t The size in bytes of the memory buffer.
	void OnEndOfBuffer(const void* /*buffer*/, UInt_t /*bufferSize*/) {}
	
	/// OnNewBlock is called whenever a new block header is found in the payload.
	/// The default behaviour of this method is to do nothing.
	/// - param const AliMUONBlockHeaderStruct* This is a pointer to the block header
	///                as found in the DDL payload.
	/// - param const void* This is a pointer to the start of the block's contents.
	/// Note: both pointers point into the memory buffer being parsed, so the
	/// contents must not be modified. On the other hand this is very efficient
	/// because no memory copying is required.
	void OnNewBlock(const AliMUONBlockHeaderStruct* /*header*/, const void* /*data*/) {}
	
	/// OnEndOfBlock is called whenever a block has been processed. Symmetric
	/// calls are made to OnEndOfBlock after each call to OnNewBlock. This happens
	/// once all DSP structures contained inside the current block have been
	/// processed.
	/// The default behaviour of this method is to do nothing.
	/// - param const AliMUONBlockHeaderStruct* This is a pointer to the processed
	///                block header as found in the DDL payload.
	/// - param const void* This is a pointer to the start of the block's contents.
	/// Note: both pointers point into the memory buffer being parsed, so the
	/// contents must not be modified. On the other hand this is very efficient
	/// because no memory copying is required.
	void OnEndOfBlock(const AliMUONBlockHeaderStruct* /*header*/, const void* /*data*/) {}
	
	/// OnNewDSP is called whenever a new DSP header is found in the payload.
	/// Every DSP header received by a call to OnNewDSP is associated to the
	/// block header received in the most recent call to OnNewBlock.
	/// The default behaviour of this method is to do nothing.
	/// - param const AliMUONDSPHeaderStruct*  This is a pointer to the DSP header
	///                as found in the DDL payload.
	/// - param const void*  This is a pointer to the start of the DSP's contents.
	/// Note: both pointers point into the memory buffer being parsed, so the
	/// contents must not be modified. On the other hand this is very efficient
	/// because no memory copying is required.
	void OnNewDSP(const AliMUONDSPHeaderStruct* /*header*/, const void* /*data*/) {}
	
	/// OnEndOfDSP is called whenever a DSP header has already been processed.
	/// For every call to OnNewDSP a symmetric call to OnEndOfDSP is made once
	/// all the bus patch structured contained in the DSP are processed.
	/// The default behaviour of this method is to do nothing.
	/// - param const AliMUONDSPHeaderStruct*  This is a pointer to the already
	///                processed DSP header as found in the DDL payload.
	/// - param const void*  This is a pointer to the start of the DSP's contents.
	/// Note: both pointers point into the memory buffer being parsed, so the
	/// contents must not be modified. On the other hand this is very efficient
	/// because no memory copying is required.
	void OnEndOfDSP(const AliMUONDSPHeaderStruct* /*header*/, const void* /*data*/) {}
	
	/// OnNewBusPatch is called whenever a new bus patch header is found in
	/// the payload. Every bus patch received by a call to OnNewBusPatch is
	/// associated to the DSP header received in the most recent call to OnNewDSP.
	/// The default behaviour of this method is to do nothing.
	/// - param const AliMUONBusPatchHeaderStruct*  This is a pointer to the bus patch
	///                header as found in the DDL payload.
	/// - param const void*  This is a pointer to the start of the bus patch's contents,
	///              specifically the raw data words.
	/// Note: both pointers point into the memory buffer being parsed, so the
	/// contents must not be modified. On the other hand this is very efficient
	/// because no memory copying is required.
	void OnNewBusPatch(const AliMUONBusPatchHeaderStruct* /*header*/, const void* /*data*/) {}
	
	/// OnEndOfBusPatch is called whenever a bus patch has been processed.
	/// For every call to OnNewBusPatch a symmetric call to OnEndOfBusPatch is
	/// made once the bus patch is completely processed (no more OnData calls).
	/// The default behaviour of this method is to do nothing.
	/// - param const AliMUONBusPatchHeaderStruct*  This is a pointer to the already
	///                processed bus patch header, as found in the DDL payload.
	/// - param const void*  This is a pointer to the start of the bus patch's contents,
	///              specifically the raw data words.
	/// Note: both pointers point into the memory buffer being parsed so the
	/// contents must not be modified. On the other hand this is very efficient
	/// because no memory copying is required.
	void OnEndOfBusPatch(const AliMUONBusPatchHeaderStruct* /*header*/, const void* /*data*/) {}
	
	/// OnData is called for every raw data word found within a bus patch.
	/// Every data ward received by a call to OnData is associated to the bus patch
	/// header received in the most recent call to OnNewBusPatch.
	/// The default behaviour of this method is to do nothing.
	/// - param UInt_t  This is the raw data word as found within the bus patch payload.
	/// - param bool  Flag indicating if the raw data word had a parity error.
	///       This will always be set to false if fSendDataOnParityError in the
	///       AliMUONTrackerDDLDecoder class was set to false.
	void OnData(UInt_t /*data*/, bool /*parityError*/) {}
	
	/// Whenever a parsing error of the DDL payload is encountered because of
	/// corruption of the raw data (eg. bit flips) the OnError method is called
	/// immediately at the point this error is discovered.
	/// The default behaviour of this method is to do nothing.
	/// - param ErrorCode  This is an error code indicating the kind of problem
	///               encountered with the DDL payload.
	/// - param const void*  This is a pointer into the DDL payload memory buffer
	///         indicating the exact location where the parsing error happened
	///         or i.e. the location of the corruption.
	/// Note that a relative offset in bytes from the start of the memory buffer
	/// can be calculated by: storing the buffer pointer received in OnNewBuffer
	/// earlier in fBufferStart for example, and then the offset is given by:
	///   offset = (unsigned long)location - (unsigned long)fBufferStart;
	void OnError(ErrorCode /*error*/, const void* /*location*/) {}
	
	/// This is a utility method which will unpack the MANU ID, channel ID and
	/// ADC signal value from a raw data word. It should normally be used in
	/// OnData() to unpack these fields.
	/// [in]  \param data  This is the raw data word found in the DDL payload.
	/// [out] \param manuId    This is filled with the unpacked MANU ID.
	/// [out] \param channelId This is filled with the unpacked MANU channel ID.
	/// [out] \param adc       This is filled with the unpacked ADC signal.
	static void UnpackADC(
			UInt_t data,
			UShort_t& manuId, UChar_t& channelId, UShort_t& adc
		)
	{
		manuId = (UShort_t)(data >> 18) & 0x7FF;
		channelId = (Char_t)(data >> 12) & 0x3F;
		adc = (UShort_t)(data & 0xFFF);
	}
	
	/// This is a utility method which converts an error code to a string
	/// representation for printing purposes.
	/// \param code  The error code as received in OnError for example.
	/// \return  An ANSI string containing the name of the error code symbol.
	static const char* ErrorCodeToString(ErrorCode code);
	
	/// This is a utility method which converts an error code to user friendly
	/// descriptive message useful for printing to the screen.
	/// \param code  The error code as received in OnError for example.
	/// \return  An ANSI string containing a descriptive message of the error.
	static const char* ErrorCodeToMessage(ErrorCode code);
};

//_____________________________________________________________________________

inline const char* AliMUONTrackerDDLDecoderEventHandler::ErrorCodeToString(ErrorCode code)
{
	/// This is a utility method which converts an error code to a string
	/// representation for printing purposes.
	/// \param code  The error code as received in OnError for example.
	/// \return  An ANSI string containing the name of the error code symbol.
	
	switch (code)
	{
	case kNoError: return "kNoError";
	case kBufferTooBig: return "kBufferTooBig";
	case kTooManyBlocks: return "kTooManyBlocks";
	case kTooManyDSPs: return "kTooManyDSPs";
	case kTooManyBusPatches: return "kTooManyBusPatches";
	case kNoBlockHeader: return "kNoBlockHeader";
	case kBadBlockKey: return "kBadBlockKey";
	case kBadBlockLength: return "kBadBlockLength";
	case kBadBlockTotalLength: return "kBadBlockTotalLength";
	case kBlockLengthMismatch: return "kBlockLengthMismatch";
	case kNoDSPHeader: return "kNoDSPHeader";
	case kBadDSPKey: return "kBadDSPKey";
	case kBadDSPLength: return "kBadDSPLength";
	case kBadDSPTotalLength: return "kBadDSPTotalLength";
	case kDSPLengthMismatch: return "kDSPLengthMismatch";
	case kNoBusPatchHeader: return "kNoBusPatchHeader";
	case kBadBusPatchKey: return "kBadBusPatchKey";
	case kBadBusPatchLength: return "kBadBusPatchLength";
	case kBadBusPatchTotalLength: return "kBadBusPatchTotalLength";
	case kBusPatchLengthMismatch: return "kBusPatchLengthMismatch";
	case kNoDDLTrailerWords: return "kNoDDLTrailerWords";
	case kTooFewDDLTrailerWords: return "kTooFewDDLTrailerWords";
	case kUnknownDspError: return "kUnknownDspError";
	case kTokenLost: return "kTokenLost";
	case kGlitchFound: return "kGlitchFound";
	case kBadPaddingWord: return "kBadPaddingWord";
	case kParityError: return "kParityError";
	default: return "INVALID";
	}
}


inline const char* AliMUONTrackerDDLDecoderEventHandler::ErrorCodeToMessage(ErrorCode code)
{
	/// This is a utility method which converts an error code to user friendly
	/// descriptive message useful for printing to the screen.
	/// \param code  The error code as received in OnError for example.
	/// \return  An ANSI string containing a descriptive message of the error.
	
	switch (code)
	{
	case kNoError:
		return "Decoding was successful.";
	case kBufferTooBig:
		return "The DDL raw data is larger than indicated by the headers;"
		       " extra bytes are probably just garbage.";
	case kTooManyBlocks:
		return "Too many block structures found.";
	case kTooManyDSPs:
		return "Too many DSP structures found in the block.";
	case kTooManyBusPatches:
		return "Too many bus patch structures found in the DSP structure.";
	case kNoBlockHeader:
		return "Missing a block header.";
	case kBadBlockKey:
		return "The block header key word does not contain the correct value.";
	case kBadBlockLength:
		return "The block length field points past the end of the raw data size.";
	case kBadBlockTotalLength:
		return "The total block length field points past the end of the"
		       " raw data size.";
	case kBlockLengthMismatch:
		return "The block length and total length fields do not correspond."
		       " One or both of these values is incorrect.";
	case kNoDSPHeader:
		return "Missing a DSP header.";
	case kBadDSPKey:
		return "The DSP header key word does not contain the correct value.";
	case kBadDSPLength:
		return "The DSP structure length field points past the end of the"
		       " block structure.";
	case kBadDSPTotalLength:
		return "The total DSP structure length field points past the end of"
		       " the block structure.";
	case kDSPLengthMismatch:
		return "The DSP structure length and total length fields do not"
		       " correspond. One or both of these values is incorrect.";
	case kNoBusPatchHeader:
		return "Missing a bus patch header.";
	case kBadBusPatchKey:
		return "The bus patch header key word does not contain the correct value.";
	case kBadBusPatchLength:
		return "The bus patch length field points past the end of the"
		       " DSP structure.";
	case kBadBusPatchTotalLength:
		return "The total bus patch length field points past the end of"
		       " the DSP structure.";
	case kBusPatchLengthMismatch:
		return "The bus patch length and total length fields do not correspond."
		       " One or both of these values is incorrect.";
	case kNoDDLTrailerWords:
		return "No end of DDL data key found in the trailer words.";
	case kTooFewDDLTrailerWords:
		return "Only one end of DDL data key word found in the trailer but expected two.";
	case kUnknownDspError:
		return "The DSP error code is non-zero but of an unrecognised format.";
	case kTokenLost:
		return "The DSP contains a token lost error code that can affect the deadtime.";
	case kGlitchFound:
		return "Found a glitch. This means a 1 byte word has been randomly"
		       " inserted into the raw data by mistake.";
	case kBadPaddingWord:
		return "The padding word does not contain the correct value.";
	case kParityError:
		return "Found a parity error in the data word.";
	default:
		return "Unknown error code!";
	}
}


inline std::ostream& operator << (std::ostream& os, AliMUONTrackerDDLDecoderEventHandler::ErrorCode code)
{
	/// This is the stream operator for std::ostream classes to be able to
	/// easily write the error messages associated with the error codes generated
	/// by the decoder to 'cout' or 'cerr' for example.
	
	os << AliMUONTrackerDDLDecoderEventHandler::ErrorCodeToMessage(code);
	return os;
}

#endif // ALIMUONTRACKERDDLDECODEREVENTHANDLER_H

