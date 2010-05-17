#ifndef ALIMUONTRACKERDDLDECODER_H
#define ALIMUONTRACKERDDLDECODER_H
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
/// \file   AliMUONTrackerDDLDecoder.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   28-11-2007
/// \brief  Implementation of a high performance DDL decoder for the muon tracking stations.
///
/// This file implementes the AliMUONTrackerDDLDecoder class, which contains
/// the core logic for decoding the payload in DDL streams coming from the muon
/// spectrometer's tracking chambers in a very efficient manner.
///
/// This implementation is derived from work done by Christian Finck for the
/// AliMUONPayloadTracker.
///
/// Note to maintainers: Please remember that this file is used by the online
/// dHLT system. As an online system, the dHLT requires the fastest code possible
/// in the decoders to satisfy its timing constraints. The performance impact
/// must be checked before any proposed modification is made to this file.
///

#include <string.h>
#include "AliMUONTrackerDDLDecoderEventHandler.h"

/// \ingroup raw
/// \class AliMUONTrackerDDLDecoder
/// \brief A high performance decoder class for MUON tracking DDL data.
///
/// This class implements a high performance decoder for decoding DDL payload
/// data coming from the muon spectrometers tracking chambers.
/// It has been implemented using the event driven paradigm with templates,
/// which allows us to minimise the number of method calls made in the inner
/// loops of the algorithm and minimise the memory footprint. At least for
/// optimised production compilations.
/// The decoder class only contains the basic decoding and error checking logic.
/// It calls methods such as OnNewBlock, OnNewBusPatch, OnData etc in
/// the event handler during the decoding to return the decoded data.
/// The event handler class is nothing more than a callback interface to deliver
/// the next chunks of decoded data.
/// To actually do something with the data, one needs to implement a custom
/// event handler (callback) class by inheriting from AliMUONTrackerDDLDecoderEventHandler
/// and overriding the callback methods like so:
/// \code
///  class MyCustomHandler : public AliMUONTrackerDDLDecoderEventHandler
///  {
///  public:
///     void OnData(UInt_t data, bool parityError)
///     {
///       // I can do something with 'data' here and check if there was
///       // a parity error with 'parityError'.
///     }
///  };
/// \endcode
///
/// Once the custom handler is written then the decoder is instantiated as
/// shown below, to use your new custom handler. Also to start decoding one needs
/// to call the Decode() method of the decoder.
/// \code
///  AliMUONTrackerDDLDecoder<MyCustomHandler> myDecoder;
///  muDecoder.Decoder(buffer, bufferSize);
/// \endcode
///
/// Note that this class was written as a template on purpose. To maximise the
/// compilers chance to make optimisations and inline the code we must use a template.
/// Depending on exactly what you do inside your handler, the decoder could be
/// significantly slower if run time polymorphism was used, i.e. making the class
/// AliMUONTrackerDDLDecoderEventHandler abstract and using virtual methods.
///
/// There has been a change to the data format that the real detector generates.
/// Two trailer words are added to the end of the DDL payload which indicated
/// the end of data. The decoder is initialised by default to automatically
/// check for these and deal with it correctly, if they exist or not.
/// However, if you want to override this behaviour then set the flag
/// fAutoDetectTrailer to false with AutoDetectTrailer(false). Then if you have
/// data with the old data format you should set fCheckForTrailer to false with
/// CheckForTrailer(false), otherwise for real data it should be
/// fCheckForTrailer = true. Only when fAutoDetectTrailer is true will the
/// fCheckForTrailer flag be ignored and no warnings will be generated for an
/// incorrect data format.
///
/// \author Artur Szostak <artursz@iafrica.com>

template <class EventHandler>
class AliMUONTrackerDDLDecoder
{
public:

	/// Default contructor.
	AliMUONTrackerDDLDecoder() :
		fExitOnError(true), fTryRecover(false),
		fSendDataOnParityError(false), fHadError(false),
		fAutoDetectTrailer(true), fCheckForTrailer(true),
		fMaxBlocks(2), fMaxDSPs(5), fMaxBusPatches(5), fHandler()
	{}
	
	/// Constant method to return the event handler instance.
	const EventHandler& GetHandler() const { return fHandler; }
	
	/// Returns the event handler instance.
	EventHandler& GetHandler() { return fHandler; }
	
	/// Returns the "exit on error" flag.
	/// i.e. should the decoder stop on the very first error found.
	bool ExitOnError() const { return fExitOnError; }
	
	/// Sets the "exit on error" flag.
	/// i.e. should the decoder stop on the very first error found.
	void ExitOnError(bool value) { fExitOnError = value; }
	
	/// Returns the "try to recover from errors" flag.
	/// i.e. should the decoder try to recover from errors found in the
	/// payload headers or trailers.
	bool TryRecover() const { return fTryRecover; }
	
	/// Sets the "try to recover from errors" flag.
	/// i.e. should the decoder try to recover from errors found in the
	/// payload headers or trailers.
	void TryRecover(bool value) { fTryRecover = value; }
	
	/// Returns the flag indicating if the raw data words in the bus patches
	/// that failed their parity tests (i.e. parity error / bit flip in the
	/// raw data word) will be sent to the event handler anyway through OnData.
	bool SendDataOnParityError() const { return fSendDataOnParityError; }
	
	/// Sets the flag indicating if the raw data words in the bus patches
	/// that failed their parity tests (i.e. parity error / bit flip in the
	/// raw data word) will be sent to the event handler anyway through OnData.
	void SendDataOnParityError(bool value) { fSendDataOnParityError = value; }
	
	/// Returns the maximum block count expected in the DDL payload.
	UInt_t MaxBlocks() const { return fMaxBlocks; }
	
	/// Sets the maximum block count expected in the DDL payload.
	void MaxBlocks(UInt_t n) { fMaxBlocks = n; }
	
	/// Returns the maximum DSP header count expected in any given block
	/// structure within the DDL payload.
	UInt_t MaxDSPs() const { return fMaxDSPs; }
	
	/// Sets the maximum DSP header count expected in any given block structure
	/// within the DDL payload.
	void MaxDSPs(UInt_t n) { fMaxDSPs = n; }
	
	/// Returns the maximum number of bus patches expected in any given DSP
	/// structure within the DDL payload.
	UInt_t MaxBusPatches() const { return fMaxBusPatches; }
	
	/// Sets the maximum number of bus patches expected in any given DSP
	/// structure within the DDL payload.
	void MaxBusPatches(UInt_t n) { fMaxBusPatches = n; }
	
	/// Returns the value of the auto-detect trailer flag.
	bool AutoDetectTrailer() const { return fAutoDetectTrailer; }
	
	/// Sets the value of the auto-detect trailer flag.
	void AutoDetectTrailer(bool value) { fAutoDetectTrailer = value; }
	
	/// Returns the value of the flag to check for the end of DDL trailer.
	bool CheckForTrailer() const { return fCheckForTrailer; }
	
	/// Sets the value of the flag to check for the end of DDL trailer.
	void CheckForTrailer(bool value) { fCheckForTrailer = value; }
	
	/// This method decodes the DDL payload contained in the buffer.
	bool Decode(const void* buffer, UInt_t bufferSize);
	
	/// First try fix data corruption and then decode the DDL payload.
	bool Decode(void* buffer, UInt_t bufferSize);
	
	/// Returns the block marker key.
	static UInt_t BlockDataKeyWord() { return fgkBlockDataKey; }
	
	/// Returns the DSP marker key.
	static UInt_t DspDataKeyWord() { return fgkDSPDataKey; }
	
	/// Returns the bus patch marker key.
	static UInt_t BusPatchDataKeyWord() { return fgkBusPatchDataKey; }
	
	/// Returns the expected padding word value.
	static UInt_t PaddingWord() { return fgkPaddingWord; }
	
	/// Returns the expected end of DDL marker.
	static UInt_t EndOfDDLWord() { return fgkEndOfDDL; }
	
private:

	bool fExitOnError; ///< Indicates if we should exit on the very first error.
	bool fTryRecover; ///< Indicates if we should try recover from a corrupt structure header or DDL trailer.
	bool fSendDataOnParityError; ///< If set to true then we issue a OnData() event even if the data word had a parity error.
	bool fHadError; ///< Indicates if we had an error decoding the data.
	bool fAutoDetectTrailer; ///< Indicates if we should automatically check for the end of DDL trailer (Default = true).
	bool fCheckForTrailer; ///< Indicates if we should check for the end of DDL trailer (Default = true). This flag is ignored if fAutoDetectTrailer is true.
	UInt_t fMaxBlocks; ///< Maximum number of block structures allowed in a DDL stream.
	UInt_t fMaxDSPs; ///< Maximum number of DSP structures allowed in a DDL stream.
	UInt_t fMaxBusPatches; ///< Maximum number of bus patch structures allowed in a DDL stream.
	EventHandler fHandler; ///< The event handler which deals with parsing events.

	void DecodeBuffer(const UChar_t* start, const UChar_t* end);
	
	bool DecodeBlockData(
			const AliMUONBlockHeaderStruct* blockHeader,
			const UChar_t* start, const UChar_t* end
		);

	bool DecodeDSPData(const UChar_t* start, const UChar_t* end);
	
	bool DecodeBusPatchData(const UChar_t* start, const UChar_t* end);

	/// Possible results that can be returned by the TryRecoverStruct method.
	enum RecoverResult
	{
		kRecoverFailed,        ///< The recovery failed. Cannot continue parsing.
		kStructRecovered,      ///< Indicates that we recovered from a corrupt structure header and can continue processing the given structure.
		kContinueToNextStruct  ///< Must continue parsing the next structure and ignore the current one.
	};

	RecoverResult TryRecoverStruct(
			UInt_t expectedKey,
			UInt_t headerSize,
			UInt_t totalLength,
			UInt_t length,
			const UChar_t* structStart,
			const UChar_t* bufferEnd,
			const UChar_t*& dataEnd,
			const UChar_t*& structEnd,
			const UChar_t*& current
		);
	
	const UChar_t* FindKey(
			UInt_t key, const UChar_t* start, const UChar_t* end
		);
	
	bool ParityIsOk(UInt_t data);
	
	static const UInt_t fgkBlockDataKey;     ///< The key word expected to identify block structure headers.
	static const UInt_t fgkDSPDataKey;       ///< The key word expected to identify DSP structure headers.
	static const UInt_t fgkBusPatchDataKey;  ///< The key word expected to identify bus patch headers.
	static const UInt_t fgkPaddingWord;      ///< The expected format of the padding word in the DDL payload.
	static const UInt_t fgkEndOfDDL;         ///< The end of DDL trailer word.
};

//_____________________________________________________________________________

// The following are the structure header keys which are used to identify the kind
// of structure header we are dealing with: block, DSP or bus patch header.
template <class EventHandler>
const UInt_t AliMUONTrackerDDLDecoder<EventHandler>::fgkBlockDataKey = 0xFC0000FC;
template <class EventHandler>
const UInt_t AliMUONTrackerDDLDecoder<EventHandler>::fgkDSPDataKey = 0xF000000F;
template <class EventHandler>
const UInt_t AliMUONTrackerDDLDecoder<EventHandler>::fgkBusPatchDataKey = 0xB000000B;
template <class EventHandler>
const UInt_t AliMUONTrackerDDLDecoder<EventHandler>::fgkPaddingWord = 0xBEEFFACE;
template <class EventHandler>
const UInt_t AliMUONTrackerDDLDecoder<EventHandler>::fgkEndOfDDL = 0xD000000D;


template <class EventHandler>
bool AliMUONTrackerDDLDecoder<EventHandler>::Decode(const void* buffer, UInt_t bufferSize)
{
	/// This method should be called to actually decode the DDL payload
	/// contained in a memory buffer. The payload should be for a muon tracking
	/// chamber DDL stream.
	/// As the decoder progresses it will make method calls to the event handler
	/// instance (which can be accessed with the GetHandler() method) to indicate
	/// the start of the new block, DSP and bus patch headers. For every raw
	/// data word the OnData method of the event handler is called.
	///
	/// If an error occurs during the parse because the data is corrupt then
	/// the OnError method is called indicating what the problem was.
	/// Decoding will stop at this point unless the fExitOnError flag is set
	/// to false. Also raw data words which contain a parity error are only
	/// sent to the event handler with OnData if the fSendDataOnParityError
	/// flag is set to true. There is also an optional flag fTryRecover which
	/// can enable logic which will attempt to recover the header structures found
	/// in the DDL payload if they are found to be inconsistent (assumed corrupt).
	/// fTryRecover set to true will also enable recovery from a corrupt
	/// DDL trailer marking the end of DDL payload.
	///
	/// \param buffer  This is the pointer to the start of the memory buffer
	///     containing the DDL payload. Remember that this must be the start of
	///     the payload and not the DDL stream. That is, this pointer should be
	///     equal to: DDL start pointer + 8 * sizeof(UInt_t).
	/// \param bufferSize  This is the pointer to the first byte just past the
	///     end of the block structure.
	/// \return Returns false if there was any problem with decoding the data,
	///     and true otherwise. Note: the data may have been partially decoded
	///     even if false was returned, which would be indicated by at least one
	///     call to the event handlers OnData method.
	
	assert( buffer != NULL );
	
	fHadError = false;
	
	// We are basically implementing something like a recursive decent parser.
	// So start by marking the start of buffer position and end of buffer.
	const UChar_t* start = reinterpret_cast<const UChar_t*>(buffer);
	const UChar_t* end = start + bufferSize;
	
	fHandler.OnNewBuffer(buffer, bufferSize);
	DecodeBuffer(start, end);
	fHandler.OnEndOfBuffer(buffer, bufferSize);
	return not fHadError;
}


template <class EventHandler>
bool AliMUONTrackerDDLDecoder<EventHandler>::Decode(void* buffer, UInt_t bufferSize)
{
	/// This decoding method performs a special checks to see if the DDL
	/// corruption is fixable and then fixes the buffer by modifying it directly.
	/// The fixes only apply if the fTryRecover flag is enabled.
	/// \note buffer might be modified.
	
	assert( buffer != NULL );
	if (fTryRecover)
	{
		///////////////////////////////////////////////////////
		// Trick to recover B off data Runs : 119041, 119047, 119055, 119057
		// A. Baldisseri May 2010
		// Check if 16 32-bit words from a buspatch have been inserted
		// incorrectly immediately at the beginning of the DDL payload.
		UChar_t* bufferChar = reinterpret_cast<UChar_t*>(buffer);
		UInt_t* bufferInt = reinterpret_cast<UInt_t*>(buffer);
		if (bufferSize > 18*4 // is the buffer large enough to check the header.
		    and bufferInt[0] != fgkBlockDataKey and bufferInt[16] == fgkBlockDataKey // was the header incorrectly moved.
		    and bufferSize > bufferInt[17]*4 + sizeof(AliMUONBlockHeaderStruct) // is the buffer large enough for the second block header.
		    and bufferInt[bufferInt[17]] == fgkBlockDataKey  // Check that the second header is in the correct location.
		    and bufferInt[17] == bufferInt[18] + 8  // Make sure that both lengths are correct.
		    and (bufferInt[17] + bufferInt[bufferInt[17]+1]+2)*4 == bufferSize // Check that both blocks will have the correct size if we make the fix.
		   )
		{
			UChar_t tmpbuf[16*4];
			memcpy(tmpbuf, bufferChar, 16*4);
			size_t sizeToMove = bufferInt[17]*4-16*4;
			memmove(bufferChar, bufferChar+16*4, sizeToMove);
			memcpy(bufferChar + sizeToMove, tmpbuf, 16*4);
		}
	}
	return Decode(reinterpret_cast<const void*>(buffer), bufferSize);
}


template <class EventHandler>
void AliMUONTrackerDDLDecoder<EventHandler>::DecodeBuffer(
		const UChar_t* start, const UChar_t* end
	)
{
	/// This method decodes the buffer's payload data. It unpacks the block
	/// structures contained inside and then for each block it calls the
	/// OnNewBlock method for the event handler to signal the start of each new
	/// block structure. OnEndOfBlock is called once each block is processed.
	/// \param start  This is the pointer to the start of the buffer.
	/// \param end  This is the pointer to the first byte just past the
	///             end of the buffer.
	/// fHadError is set to true if there were any errors decoding the buffer
	/// and the OnError method of the callback event handler is called for
	/// each error.
	
	const UChar_t* current = start;
	const UInt_t* bufferStart = reinterpret_cast<const UInt_t*>(start);
	const UInt_t* bufferEnd = reinterpret_cast<const UInt_t*>(end);
	bool problemWithTrailer = false;
	
	// The DDL payload normally has a 2 word trailer which contains the end of
	// DDL markers 0xD000000D. But this is not the case for older simulated
	// data so if we are auto-detecting the trailer then we need to carefully
	// check if these words are there or not.
	const UChar_t* endOfBlocks = end;
	const UInt_t* trailerWords = reinterpret_cast<const UInt_t*>(end) - 2;
	if (fAutoDetectTrailer)
	{
		if (trailerWords >= bufferStart and trailerWords < bufferEnd
		    and *trailerWords == fgkEndOfDDL and *(trailerWords+1) == fgkEndOfDDL
		   )
		{
			// Found the trailer so reposition the end of blocks marker.
			endOfBlocks = reinterpret_cast<const UChar_t*>(trailerWords);
		}
		// else assume we are dealing with the older data format.
	}
	else if (fCheckForTrailer)
	{
		if (trailerWords >= bufferStart and trailerWords < bufferEnd
		    and *trailerWords == fgkEndOfDDL and *(trailerWords+1) == fgkEndOfDDL
		   )
		{
			// Found the trailer so reposition the end of blocks marker.
			endOfBlocks = reinterpret_cast<const UChar_t*>(trailerWords);
		}
		else
		{
			if (trailerWords+1 >= bufferStart and trailerWords+1 < bufferEnd and *(trailerWords+1) == fgkEndOfDDL)
				fHandler.OnError(EventHandler::kTooFewDDLTrailerWords, trailerWords+1);
			else if (trailerWords >= bufferStart and trailerWords < bufferEnd and *(trailerWords) == fgkEndOfDDL)
				fHandler.OnError(EventHandler::kTooFewDDLTrailerWords, trailerWords);
			else
				fHandler.OnError(EventHandler::kNoDDLTrailerWords, end);
	
			// Stop the decoding if so requested by the user, otherwise
			// remember about the error so that we return false from the
			// Decode() method and continue decoding.
			fHadError = true;
			if (fExitOnError) return;
			
			// Mark that there was a problem with the trailer so that
			// for subsequent errors we try to deal with this better.
			problemWithTrailer = true;
			
			// We can also try figure out how many trailer words there
			// actually are and move the end of blocks marker back.
			
			if (fTryRecover)
			{
				trailerWords = bufferEnd;
				// There should only be a max of 2 trailer words.
				if (trailerWords-2 >= bufferStart and trailerWords-2 < bufferEnd and *(trailerWords-2) == fgkEndOfDDL)
					trailerWords -= 2;
				else if (trailerWords-1 >= bufferStart and trailerWords-1 < bufferEnd and *(trailerWords-1) == fgkEndOfDDL)
					trailerWords -= 1;
				endOfBlocks = reinterpret_cast<const UChar_t*>(trailerWords);
			}
		}
	}

	UInt_t blockCount = 0; // Indicates the number of blocks decoded.
	while (current < endOfBlocks)
	{
		// Mark the start of the block structure.
		const UChar_t* blockStart = current;
		
		// Get the block header, move the current pointer just past the end
		// of the header and check that we have not overflowed the buffer.
		const AliMUONBlockHeaderStruct* blockHeader
			= reinterpret_cast<const AliMUONBlockHeaderStruct*>(blockStart);
		current += sizeof(AliMUONBlockHeaderStruct);
		if (current > endOfBlocks or current < start)
		{
			// If we overflowed the pointer and already had an error then
			// we are clearly lost so just stop decoding before we segfault.
			if (current < start and fHadError) return;
			
			// We first check if we actually hit the end of DDL markers
			// If we did then either we did not/could not recover from
			// a corrupt trailer or we did not detect a correct trailer
			// in auto-detect mode.
			trailerWords = reinterpret_cast<const UInt_t*>(blockHeader);
			// The "trailerWords+1 <= bufferEnd" checks that we are
			// not reading beyond the end of the buffer.
			if (trailerWords+1 <= bufferEnd and *trailerWords == fgkEndOfDDL)
			{
				// If we aready knew the trailer was corrupt then just
				// return because the error was already announced.
				if (problemWithTrailer) return;
				
				if (fAutoDetectTrailer)
				{
					// If we got here then there is at least one correct trailer
					// word, but since we did not detect a correct trailer then
					// there must be only one. Announce the error and exit.
					fHandler.OnError(EventHandler::kTooFewDDLTrailerWords, trailerWords);
					fHadError = true;
					return;
				}
			}
			
			// So we only got part of a block header at the very end
			// of the buffer. Nothing to do but report the error and exit.
			if (blockCount == fMaxBlocks)
				// Special case where we got all the blocks we
				// expected, so the remaining data must be rubbish.
				fHandler.OnError(EventHandler::kBufferTooBig, blockHeader);
			else
				fHandler.OnError(EventHandler::kNoBlockHeader, blockHeader);
			fHadError = true;
			return;
		}
		
		// The header fits the buffer so we can mark the data start and
		// read from the header to find the end of data and block pointers.
		const UChar_t* dataStart = current;
		current += blockHeader->fLength * sizeof(UInt_t);
		const UChar_t* dataEnd = current;
		const UChar_t* blockEnd = blockStart
			+ blockHeader->fTotalLength * sizeof(UInt_t);
		
		// Now we need to check for the following things:
		// 1) Is the end of block or end of data pointer outside the buffer
		//    boundaries.
		// 2) Are the values for these pointers the same.
		// 3) Is the expected data key in the header present.
		// If any of the above fail then we know there is a problem with
		// the block header. It must be corrupted somehow.
		if (blockHeader->fDataKey != fgkBlockDataKey
		    or dataEnd > endOfBlocks or dataEnd < start
		    or blockEnd > endOfBlocks or blockEnd < start
		    or dataEnd != blockEnd)
		{
			// So let us see what exactly is wrong and report this.
			if (blockCount == fMaxBlocks)
			{
				// Special case where we got all the blocks we
				// expected, so the remaining data must be rubbish.
				// Don't even bother trying to recover the data.
				fHandler.OnError(EventHandler::kBufferTooBig, blockHeader);
				fHadError = true;
				return;
			}
			if (blockHeader->fDataKey != fgkBlockDataKey)
				fHandler.OnError(EventHandler::kBadBlockKey, &blockHeader->fDataKey);
			if (blockEnd > endOfBlocks or blockEnd < start)
				fHandler.OnError(EventHandler::kBadBlockLength, &blockHeader->fLength);
			if (dataEnd > endOfBlocks or dataEnd < start)
				fHandler.OnError(EventHandler::kBadBlockTotalLength, &blockHeader->fTotalLength);
			if (dataEnd != blockEnd)
				fHandler.OnError(EventHandler::kBlockLengthMismatch, blockHeader);
			
			// Stop the decoding if so requested by the user, otherwise
			// remember about the error so that we return false from the
			// Decode() method and continue decoding.
			fHadError = true;
			if (fExitOnError) return;
			
			// Try to recover from the corrupt header.
			RecoverResult result = TryRecoverStruct(
					fgkBlockDataKey, sizeof(AliMUONBlockHeaderStruct),
					blockHeader->fTotalLength, blockHeader->fLength,
					blockStart, endOfBlocks, dataEnd, blockEnd, current
				);
			if (result == kContinueToNextStruct)
				continue; // Try the next block at 'current'.
			if (result == kRecoverFailed) return;
		}
		
		// At this point we certainly have a valid block header, so we
		// need to check if we have more blocks than we expected. If not
		// then we can indicate we have another block and decode its data.
		if (++blockCount > fMaxBlocks)
		{
			fHandler.OnError(EventHandler::kTooManyBlocks, current);
			
			// In this case we stop the decoding because clearly
			// something is seriously wrong with the data if we are
			// getting more blocks than expected.
			fHadError = true;
			return;
		}
		
		fHandler.OnNewBlock(blockHeader, dataStart);
		if (not DecodeBlockData(blockHeader, dataStart, dataEnd))
		{
			// At this point we had a problem decoding the block structure's
			// data. Thus we should stop further decoding if so requested by
			// the user. Note the fHadError flag is already marked inside
			// DecodeBlockData.
			if (fExitOnError)
			{
				fHandler.OnEndOfBlock(blockHeader, dataStart);
				return;
			}
		}
		fHandler.OnEndOfBlock(blockHeader, dataStart);
	}
}


template <class EventHandler>
bool AliMUONTrackerDDLDecoder<EventHandler>::DecodeBlockData(
		const AliMUONBlockHeaderStruct* blockHeader,
		const UChar_t* start, const UChar_t* end
	)
{
	/// This method decodes a block structure's data payload. It unpacks the
	/// DSP structures contained inside and then for each DSP it calls the
	/// OnNewDSP method for the event handler to signal the start of each new
	/// DSP structure.
	/// \param blockHeader
	/// \param start  This is the pointer to the start of the block
	///               structure's data.
	/// \param end  This is the pointer to the first byte just past the
	///             end of the block structure.
	/// \return If the block structure's data was decoded without errors
	///      or we could recover from the errors, then true is returned.
	///      False is returned otherwise.
	
	const UChar_t* current = start;
	
	UInt_t dspCount = 0; // Indicates the number of DSPs decoded.
	while (current < end)
	{
		// Mark the start of the DSP structure.
		const UChar_t* dspStart = current;
		
		// Get the DSP header, move the current pointer just past the end
		// of the header and check that we have not overflowed the buffer.
		const AliMUONDSPHeaderStruct* dspHeader
			= reinterpret_cast<const AliMUONDSPHeaderStruct*>(dspStart);
		current += sizeof(AliMUONDSPHeaderStruct);
		if (current > end or current < start)
		{
			// If we overflowed the pointer and already had an error then
			// we are clearly lost so just stop decoding before we segfault.
			if (current < start and fHadError) return false;
			
			// So we only got part of a DSP header at the very end of
			// the block structure buffer. Nothing to do but report the
			// error and exit. Set fHadError in case of further decoding.
			fHandler.OnError(EventHandler::kNoDSPHeader, dspHeader);
			fHadError = true;
			return false;
		}
		
		// The header fits the buffer so we can mark the data start and
		// read from the header to find the end of data and DSP pointers.
		const UChar_t* dataStart = current;
		current += dspHeader->fLength * sizeof(UInt_t);
		const UChar_t* dataEnd = current;
		const UChar_t* dspEnd = dspStart + dspHeader->fTotalLength * sizeof(UInt_t);
		
		// Now we need to check for the following things:
		// 1) Is the end of DSP or end of data pointer outside the buffer
		//    boundaries.
		// 2) Are the values for these pointers the same.
		// 3) Is the expected data key in the header present.
		// If any of the above fail then we know there is a problem with
		// the DSP header. It must be corrupted somehow.
		if (dspHeader->fDataKey != fgkDSPDataKey
		    or dataEnd > end or dataEnd < start
		    or dspEnd > end or dspEnd < start
		    or dataEnd != dspEnd)
		{
			// So let us see what exactly is wrong and report this.
			if (dspHeader->fDataKey != fgkDSPDataKey)
				fHandler.OnError(EventHandler::kBadDSPKey, &dspHeader->fDataKey);
			if (dspEnd > end or dspEnd < start)
				fHandler.OnError(EventHandler::kBadDSPLength, &dspHeader->fLength);
			if (dataEnd > end or dataEnd < start)
				fHandler.OnError(EventHandler::kBadDSPTotalLength, &dspHeader->fTotalLength);
			if (dataEnd != dspEnd)
				fHandler.OnError(EventHandler::kDSPLengthMismatch, dspHeader);
			
			// Indicate we had and error and stop the decoding if so
			// requested by the user.
			fHadError = true;
			if (fExitOnError) return false;
			
			// Try to recover from the corrupt header.
			RecoverResult result = TryRecoverStruct(
					fgkDSPDataKey, sizeof(AliMUONDSPHeaderStruct),
					dspHeader->fTotalLength, dspHeader->fLength,
					dspStart, end, dataEnd, dspEnd, current
				);
			if (result == kContinueToNextStruct)
				continue; // Try the next DSP at 'current'.
			if (result == kRecoverFailed) return false;
		}
		
		// At this point we certainly have a valid DSP header, so we
		// need to check if we have more DSPs than we expected. If not
		// then we can indicate we have another DSP and decode its data.
		if (++dspCount > fMaxDSPs)
		{
			fHandler.OnError(EventHandler::kTooManyDSPs, current);
			
			// In this case we stop further decoding of the block
			// structure data because clearly something is seriously
			// wrong if we are getting more DSPs than expected.
			// Indicate that we had an error so the Decode() method
			// returns false.
			fHadError = true;
			return false;
		}
		
		fHandler.OnNewDSP(dspHeader, dataStart);
		
		// Check the error word in the header.
		if (dspHeader->fErrorWord != 0x0)
		{
			if (dspHeader->fErrorWord == (0x000000B1 | blockHeader->fDSPId)
			    or dspHeader->fErrorWord == (0x00000091 | blockHeader->fDSPId)
			   )
			{
				// An event with a glitch in the readout has been detected.
				// It means that somewhere a 1 byte word has been randomly
				// inserted and all the readout sequence is shifted until
				// the next event.
				fHandler.OnError(EventHandler::kGlitchFound, &dspHeader->fErrorWord);
			}
			else if ((dspHeader->fErrorWord & 0x0000FFF0) == 0x220)
			{
				// Detected a TOKEN_LOST error which can affect the dead time in the DAQ.
				fHandler.OnError(EventHandler::kTokenLost, &dspHeader->fErrorWord);
			}
			else
			{
				// The DSP error code is non-zero but has an unknown code.
				fHandler.OnError(EventHandler::kUnknownDspError, &dspHeader->fErrorWord);
			}
			
			fHadError = true;
			if (fExitOnError)
			{
				fHandler.OnEndOfDSP(dspHeader, dataStart);
				return false;
			}
			
			// Try recover by finding the very next DSP and continue
			// decoding from there. Note: to achieve all we have to do
			// is continue to the next iteration, because the logic
			// will land up calling the FindKey method within the
			// TryRecoverStruct method above.
			if (fTryRecover) continue;
		}
		
		// Check if we are padding. If we are, then the bus patch data is
		// actually 4 bytes smaller and the last word is a padding word.
		if (dspHeader->fPaddingWord == 1)
		{
			dataEnd -= sizeof(UInt_t);
			
			// Check the pad word is correct.
			const UInt_t* padWord = reinterpret_cast<const UInt_t*>(dataEnd);
			if (*padWord != fgkPaddingWord)
			{
				fHandler.OnError(EventHandler::kBadPaddingWord, padWord);
				fHadError = true;
				if (fExitOnError)
				{
					fHandler.OnEndOfDSP(dspHeader, dataStart);
					return false;
				}
			}
		}
		
		if (not DecodeDSPData(dataStart, dataEnd))
		{
			// At this point we had a problem decoding the DSP structure's
			// data, thus we should stop further decoding if so requested by
			// the user. Note the fHadError flag is already marked inside
			// DecodeDSPData.
			if (fExitOnError)
			{
				fHandler.OnEndOfDSP(dspHeader, dataStart);
				return false;
			}
		}
		fHandler.OnEndOfDSP(dspHeader, dataStart);
	}
	
	return true;
}


template <class EventHandler>
bool AliMUONTrackerDDLDecoder<EventHandler>::DecodeDSPData(
		const UChar_t* start, const UChar_t* end
	)
{
	/// This method decodes a DSP structure's data payload. It finds all the
	/// bus patches found inside and for each it calls the OnNewBusPatch method
	/// for the event handler to signal the start of each new bus patch.
	/// \param start  This is the pointer to the start of the DSP structure's data.
	/// \param end  This is the pointer to the first byte just past the
	///             end of the DSP structure.
	/// \return If the DSP structure's data was decoded without errors
	///      or we could recover from the errors, then true is returned.
	///      False is returned otherwise.
	
	const UChar_t* current = start;
	
	UInt_t busPatchCount = 0; // Indicates the number of bus patches decoded.
	while (current < end)
	{
		// Mark the start of the bus patch structure.
		const UChar_t* busPatchStart = current;
		
		// Get the bus patch header, move the current pointer just past
		// the end of the header and check that we have not overflowed
		// the buffer.
		const AliMUONBusPatchHeaderStruct* busPatchHeader
			= reinterpret_cast<const AliMUONBusPatchHeaderStruct*>(busPatchStart);
		current += sizeof(AliMUONBusPatchHeaderStruct);
		if (current > end or current < start)
		{
			// If we overflowed the pointer and already had an error then
			// we are clearly lost so just stop decoding before we segfault.
			if (current < start and fHadError) return false;
			
			// So we only got part of a bus patch header at the very
			// end of the DSP structure buffer. Nothing to do but
			// report the error and exit. Set fHadError in case of
			// further decoding.
			fHandler.OnError(EventHandler::kNoBusPatchHeader, busPatchHeader);
			fHadError = true;
			return false;
		}
		
		// The header fits the buffer so we can mark the data start and
		// read from the header to find the end of data and bus patch
		// structure pointers.
		const UChar_t* dataStart = current;
		current += busPatchHeader->fLength * sizeof(UInt_t);
		const UChar_t* dataEnd = current;
		const UChar_t* busPatchEnd = busPatchStart
			+ busPatchHeader->fTotalLength * sizeof(UInt_t);
		
		// Now we need to check for the following things:
		// 1) Is the end of bus patch structure or end of data pointer
		//    outside the buffer boundaries.
		// 2) Are the values for these pointers the same.
		// 3) Is the expected data key in the header present.
		// If any of the above fail then we know there is a problem with
		// the bus patch header. It must be corrupted somehow.
		if (busPatchHeader->fDataKey != fgkBusPatchDataKey
		    or dataEnd > end or dataEnd < start
		    or busPatchEnd > end or busPatchEnd < start
		    or dataEnd != busPatchEnd)
		{
			// So let us see what exactly is wrong and report this.
			if (busPatchHeader->fDataKey != fgkBusPatchDataKey)
				fHandler.OnError(EventHandler::kBadBusPatchKey, &busPatchHeader->fDataKey);
			if (busPatchEnd > end or busPatchEnd < start)
				fHandler.OnError(EventHandler::kBadBusPatchLength, &busPatchHeader->fLength);
			if (dataEnd > end or dataEnd < start)
				fHandler.OnError(EventHandler::kBadBusPatchTotalLength, &busPatchHeader->fTotalLength);
			if (dataEnd != busPatchEnd)
				fHandler.OnError(EventHandler::kBusPatchLengthMismatch, busPatchHeader);
			
			// Indicate we had and error and stop the decoding if so
			// requested by the user.
			fHadError = true;
			if (fExitOnError) return false;
			
			// Try to recover from the corrupt header.
			RecoverResult result = TryRecoverStruct(
					fgkBusPatchDataKey, sizeof(AliMUONBusPatchHeaderStruct),
					busPatchHeader->fTotalLength, busPatchHeader->fLength,
					busPatchStart, end, dataEnd, busPatchEnd, current
				);
			if (result == kContinueToNextStruct)
				continue; // Try the next bus patch at 'current'.
			if (result == kRecoverFailed) return false;
		}
		
		// At this point we certainly have a valid bus patch header, so
		// we need to check if we have more bus patches than we expected.
		// If not then we can indicate we have another bus patch and
		// decode its data.
		if (++busPatchCount > fMaxBusPatches)
		{
			fHandler.OnError(EventHandler::kTooManyBusPatches, current);
			
			// In this case we stop further decoding of the DSP
			// structure's data because clearly something is seriously
			// wrong if we are getting more bus patches than expected.
			// Indicate that we had an error so the Decode() method
			// returns false.
			fHadError = true;
			return false;
		}
		
		fHandler.OnNewBusPatch(busPatchHeader, dataStart);
		if (not DecodeBusPatchData(dataStart, dataEnd))
		{
			// At this point we had a problem decoding the bus patch data,
			// thus we should stop further decoding if so requested by the
			// user. Note the fHadError flag is already marked inside
			// DecodeBusPatchData.
			if (fExitOnError)
			{
				fHandler.OnEndOfBusPatch(busPatchHeader, dataStart);
				return false;
			}
		}
		fHandler.OnEndOfBusPatch(busPatchHeader, dataStart);
	}
	
	return true;
}


template <class EventHandler>
bool AliMUONTrackerDDLDecoder<EventHandler>::DecodeBusPatchData(
		const UChar_t* start, const UChar_t* end
	)
{
	/// This method decodes a single bus patch's data payload.
	/// It will check the parity of the raw data words and send them
	/// to the event handler instance with calls to OnData.
	/// \param start  This is the pointer to the start of the bus patch
	///               structure's data.
	/// \param end  This is the pointer to the first byte just past the
	///             end of the bus patch structure.
	/// \return If the bus patch's data was decoded without errors
	///      or we could recover from the errors, then true is returned.
	///      False is returned otherwise.

	// Assert that 'end' is always larger than start by n*sizeof(UInt_t)
	// where n is a positive integer. This should be the case because we
	// always add multiples of sizeof(UInt_t) to the 'current' pointer in
	// all the DecodeXYZ methods.
	assert( UInt_t(end - start) % 4 == 0 );
	
	// Now step through all the data words and issue OnData events.
	// We also need to check parity and signal OnError if it is not valid
	// for any of the data words.
	const UInt_t* data = reinterpret_cast<const UInt_t*>(start);
	const UInt_t* dataEnd = reinterpret_cast<const UInt_t*>(end);
	for (; data < dataEnd; data++)
	{
		if (ParityIsOk(*data))
		{
			fHandler.OnData(*data, false);
		}
		else
		{
			// Indicate we had a parity error and exit immediately
			// if the user so requested.
			fHandler.OnError(EventHandler::kParityError, data);
			fHadError = true;
			if (fExitOnError) return false;
			
			if (fSendDataOnParityError)
				fHandler.OnData(*data, true);
		}
	}
	
	return true;
}


template <class EventHandler>
typename AliMUONTrackerDDLDecoder<EventHandler>::RecoverResult
AliMUONTrackerDDLDecoder<EventHandler>::TryRecoverStruct(
		UInt_t expectedKey,
		UInt_t headerSize,
		UInt_t totalLength,
		UInt_t length,
		const UChar_t* structStart,
		const UChar_t* bufferEnd,
		const UChar_t*& dataEnd,
		const UChar_t*& structEnd,
		const UChar_t*& current
	)
{
	/// This method attempts to recover from a corrupt structure header by
	/// figuring out which of the structure size indicators is correct.
	/// This is possible because each header has some redundant information.
	/// The recovery procedure is only attempted if fTryRecover was set to
	/// true. If the recovery procedure is successful then this method will
	/// also update the pointers indicating the start of data, end of structure
	/// and current parsing position with the correct values.
	///
	/// [in]  \param expectedKey This is the expected block key for the header
	///           currently being processed.
	/// [in]  \param headerSize  The expected header size as given by the sizeof
	///           operator for example.
	/// [in]  \param totalLength The total length as given by the fTotalLength
	///           field in the current header being handled.
	/// [in]  \param length  The data length as given by the fLength field
	///           in the current header being handled.
	/// [in]  \param structStart A pointer to the start of the structure header.
	/// [in]  \param bufferEnd A pointer to the first byte just past the end
	///           of the buffer. This could be the pointer to the first byte
	///           just past the end of the parent structure if we are dealing
	///           with a DSP structure or bus patch. The parent structure for
	///           the DSP is a block structure and for a bus patch it is a DSP.
	/// [out] \param dataEnd This is the pointer to the first byte just past
	///           the end of the structure being processed. It should be equal to
	///           structStart + sizeof(structure header) + fLength, where fLength
	///           is the field found in the structure's header itself. This value
	///           will be corrected and updated if we could recover from the
	///           corruption in the header.
	/// [out] \param structEnd A pointer to the first byte just past the end of
	///           the structure. This value should be set equal to
	///           structStart + fTotalLength * sizeof(UInt_t), where fTotalLength
	///           is the field found in the structure's header itself. This value
	///           will be corrected and updated if we could recover from the
	///           corruption in the header.
	/// [out] \param current This is the pointer to the current location in
	///           the DDL payload being parsed. It should in principle point
	///           to the start of the structures data. This value will be
	///           corrected and updated if we could recover from the corruption
	///           in the header.
	///
	/// \return Returns the result of the recovery attempt, which can be one
	///    of the following:
	///      kRecoverFailed - The recovery failed completely so the caller
	///           cannot continue parsing any more structures. If the failure
	///           is within a DSP then one could still continue parsing
	///           from the next block. Similarly for bus patches, parsing could
	///           continue from the next DSP structure.
	///      kStructRecovered - Indicates that we recovered from a corrupt
	///           structure header and can continue processing the data of the
	///           structure in question.
	///      kContinueToNextStruct - Either fTryRecover was set to false or we
	///           could not recover from the corrupt header but we did find the
	///           start of another header matching the expected key so parsing
	///           can continue from the updated current position.

	// Check if the user wants us to try and recover from a corrupt header.
	if (not fTryRecover) return kContinueToNextStruct;
	
	// If the user wants us to try recover, then try to recover what the
	// correct values for dataEnd, structEnd and current were supposed to be.
	// The recovery procedure is as follows: We have 4 conditions for a correct
	// header:
	//   1) The header key is what we expect.
	//   2) The totalLength equals length + headerSize.
	//   3) The word at dataEnd contains a valid key. (implies length is
	//      correct.)
	//   4) The word at structEnd contains a valid key. (implies totalLength
	//      is correct.)
	// If any 2 of these conditions hold then we know that only one of the
	// header fields is corrupt and we have enough information to reconstruct
	// the third field. Note that if conditions 3 and 4 are true then this
	// implies 2 is also true. (not necessarily the other way around though.)
	// The valid key mentioned above at dataEnd and structEnd should be:
	//   a) A bus patch key, DSP key or end of buffer if expectedKey indicates
	//      a buspatch.
	//   b) A DSP key, block structure key or end of buffer if expectedKey
	//      indicates a DSP.
	//   c) A block structure key or end of buffer if expectedKey indicates
	//      a DSP.
	const UInt_t* headerKey = reinterpret_cast<const UInt_t*>(structStart);
	bool headerKeyOk = (expectedKey == *headerKey);
	
	bool lengthsMatch = (totalLength*4 == length*4 + headerSize);
	
	bool lengthIsCorrect = false;
	bool totalLengthIsCorrect = false;
	const UInt_t* keyAtDataEnd = reinterpret_cast<const UInt_t*>(dataEnd);
	const UInt_t* keyAtStructEnd = reinterpret_cast<const UInt_t*>(structEnd);
	

        if ( expectedKey == fgkBlockDataKey )
        {
		if (dataEnd == bufferEnd)
		{
			// Are we at the end of the buffer?
			lengthIsCorrect = true;
		}
		else
		{
			// Must check that we can read another 4 bytes before
			// checking the key at dataEnd.
			if (dataEnd + sizeof(UInt_t) <= bufferEnd and dataEnd + sizeof(UInt_t) > structStart)
			{
				if (*keyAtDataEnd == fgkBlockDataKey)
					lengthIsCorrect = true;
			}
		}
		
		if (structEnd == bufferEnd)
		{
			// Are we at the end of the buffer?
			totalLengthIsCorrect = true;
		}
		else
		{
			// Must check that we can read another 4 bytes before
			// checking the key at structEnd.
			if (structEnd + sizeof(UInt_t) <= bufferEnd and structEnd + sizeof(UInt_t) > structStart)
			{
				if (*keyAtStructEnd == fgkBlockDataKey)
					totalLengthIsCorrect = true;
			}
		}
        }        
			
        else if ( expectedKey == fgkDSPDataKey )
        {
		if (dataEnd == bufferEnd)
		{
			// Are we at the end of the buffer?
			lengthIsCorrect = true;
		}
		else
		{
			// Must check that we can read another 4 bytes before
			// checking the key at dataEnd.
			if (dataEnd + sizeof(UInt_t) <= bufferEnd and dataEnd + sizeof(UInt_t) > structStart)
			{
				if (*keyAtDataEnd == fgkBlockDataKey
				    or *keyAtDataEnd == fgkDSPDataKey)
					lengthIsCorrect = true;
			}
		}
		
		if (structEnd == bufferEnd)
		{
			// Are we at the end of the buffer?
			totalLengthIsCorrect = true;
		}
		else
		{
			// Must check that we can read another 4 bytes before
			// checking the key at structEnd.
			if (structEnd + sizeof(UInt_t) <= bufferEnd and structEnd + sizeof(UInt_t) > structStart)
			{
				if (*keyAtStructEnd == fgkBlockDataKey
				    or *keyAtStructEnd == fgkDSPDataKey)
					totalLengthIsCorrect = true;
			}
		}
        }        
        else if ( expectedKey == fgkBusPatchDataKey )
        {
		if (dataEnd == bufferEnd)
		{
			// Are we at the end of the buffer?
			lengthIsCorrect = true;
		}
		else
		{
			// Must check that we can read another 4 bytes before
			// checking the key at dataEnd.
			if (dataEnd + sizeof(UInt_t) <= bufferEnd and dataEnd + sizeof(UInt_t) > structStart)
			{
				if (*keyAtDataEnd == fgkDSPDataKey
				    or *keyAtDataEnd == fgkBusPatchDataKey)
					lengthIsCorrect = true;
			}
		}
		
		if (structEnd == bufferEnd)
		{
			// Are we at the end of the buffer?
			totalLengthIsCorrect = true;
		}
		else
		{
			// Must check that we can read another 4 bytes before
			// checking the key at structEnd.
			if (structEnd + sizeof(UInt_t) <= bufferEnd and structEnd + sizeof(UInt_t) > structStart)
			{
				if (*keyAtStructEnd == fgkDSPDataKey
				    or *keyAtStructEnd == fgkBusPatchDataKey)
					totalLengthIsCorrect = true;
			}
		}
        }
	
	if (headerKeyOk and lengthIsCorrect)
	{
		// totalLength was wrong, dataEnd is correct.
		structEnd = dataEnd;
		current = dataEnd;
		return kStructRecovered;
	}
	if (headerKeyOk and totalLengthIsCorrect)
	{
		// Length was wrong, structEnd is correct.
		dataEnd = structEnd;
		current = structEnd;
		return kStructRecovered;
	}
	if (lengthsMatch and lengthIsCorrect and totalLengthIsCorrect)
	{
		// The header's key was wrong but the lengths and pointers are OK.
		return kStructRecovered;
	}
	
	// Could not recover the header from the available information, so find
	// the next key in the stream that is the same as the currently expected
	// one and continue decoding from there.
	const UChar_t* location = FindKey(
			expectedKey, structStart + sizeof(UInt_t), bufferEnd
		);
	if (location != NULL)
	{
		current = location;
		return kContinueToNextStruct;
	}

	return kRecoverFailed;
}


template <class EventHandler>
const UChar_t* AliMUONTrackerDDLDecoder<EventHandler>::FindKey(
		UInt_t key, const UChar_t* start, const UChar_t* end
	)
{
	/// Searches for the first occurrence of the key value in the buffer marked by
	/// 'start' and 'end'. 'start' should point to the start of the buffer and 'end'
	/// should point to 'start + bufferSize', i.e. just past the last byte of the
	/// buffer. If the key was found then the pointer to that location is returned
	/// otherwise NULL is returned.
	
	if (end + sizeof(UInt_t) < start) return NULL;  // check for pointer overflow.
	const UChar_t* current = start;
	while (current + sizeof(UInt_t) <= end)
	{
		UInt_t data = * reinterpret_cast<const UInt_t*>(current);
		if (data == key) return current;
		current++;
	}
	return NULL;
}


template <class EventHandler>
bool AliMUONTrackerDDLDecoder<EventHandler>::ParityIsOk(UInt_t data)
{
	/// Optimised parity check addapted from:
	/// http://graphics.stanford.edu/~seander/bithacks.html#ParityParallel
	
	// parity of the 32 bits must be zero if the last bit is equal
	// to the parity of the first 31 bits.
	// Reason: the parity bit xor the parity of the first 31 bits must give
	// zero, unless there was a bit error.
	data ^= data >> 16;
	data ^= data >> 8;
	data ^= data >> 4;
	data &= 0xf;
	data = ((0x6996 >> data) & 1);
	return data == 0;
}

#endif // ALIMUONTRACKERDDLDECODER_H
