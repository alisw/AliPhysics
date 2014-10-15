#ifndef ALIMUONTRIGGERDDLDECODER_H
#define ALIMUONTRIGGERDDLDECODER_H
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
/// \file   AliMUONTriggerDDLDecoder.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   28-11-2007
/// \brief  Implementation of a high performance DDL decoder for the muon trigger.
///
/// This file implementes the AliMUONTriggerDDLDecoder class, which contains
/// the core logic for decoding the payload in DDL streams comming from the muon
/// spectrometer's hardware trigger in a very efficient manner.
///
/// This implementation is derived from work done by Christian Finck for the
/// AliMUONPayloadTrigger class.
///
/// Note to maintainers: Please remember that this file is used by the online
/// dHLT system. As an online system, the dHLT requires the fastest code possible
/// in the decoders to satisfy its timing constraints. The performance impact
/// must be checked before any proposed modification is made to this file.
///

#include "AliMUONTriggerDDLDecoderEventHandler.h"

/// \ingroup raw
/// \class AliMUONTriggerDDLDecoder
/// \brief A high performance decoder class for MUON trigger DDL data.
///
/// This class implements a high performance decoder for decoding DDL payload
/// data coming from the muon spectrometers trigger stations.
/// It has been implemented using the event driven paradigm with templates,
/// which allows us to minimise the number of method calls made in the inner
/// loops of the algorithm and minimise the memory footprint.
/// At least for optimised production compilations.
/// The decoder class only contains the basic decoding and error checking logic.
/// It calls methods such as OnNewBuffer, OnDarcHeader, OnNewRegionalStruct,
/// OnLocalStruct etc in the event handler during the decoding to return the
/// decoded data and headers.
/// The event handler class is nothing more than a callback interface to deliver
/// the next chunks of decoded data.
/// To actually do something with the data, one needs to implement a custom
/// event handler (callback) class by inheriting from AliMUONTriggerDDLDecoderEventHandler
/// and overriding the callback methods like so:
/// \code
///  class MyCustomHandler : public AliMUONTriggerDDLDecoderEventHandler
///  {
///  public:
///     void OnLocalStruct(const AliMUONLocalInfoStruct* localStruct,
///                        const AliMUONLocalScalarsStruct* scalars)
///     {
///       // I can do something with the data in 'localStruct' here.
///       // and the 'scalars' if they are not NULL.
///     }
///  };
/// \endcode
///
/// Once the custom handler is written then the decoder is instantiated as
/// shown below, to use your new custom handler. Also to start decoding one needs
/// to call the Decode() method of the decoder.
/// \code
///  AliMUONTriggerDDLDecoder<MyCustomHandler> myDecoder;
///  muDecoder.Decoder(buffer, bufferSize);
/// \endcode
///
/// Note that this class was written as a template on purpose. To maximise the
/// compilers chance to make optimisations and inline the code we must use a template.
/// Depending on exactly what you do inside your handler, the decoder could be
/// significantly slower if run time polymorphism was used, i.e. making the class
/// AliMUONTriggerDDLDecoderEventHandler abstract and using virtual methods.
///
template <class EventHandler>
class AliMUONTriggerDDLDecoder
{
public:

	/// Default contructor.
	AliMUONTriggerDDLDecoder() :
		fExitOnError(true), fTryRecover(false),
		fAutoDetectScalars(false), fHadError(false),
		fNoRegionals(0), fMaxRegionals(8), fMaxLocals(16),
		fHandler()
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
	/// payload headers.
	bool TryRecover() const { return fTryRecover; }
	
	/// Sets the "try to recover from errors" flag.
	/// i.e. should the decoder try to recover from errors found in the
	/// payload headers.
	void TryRecover(bool value) { fTryRecover = value; }
	
	/// Returns the flag indicating if we should check for the scalars of
	/// software triggers automatically or not.
	bool AutoDetectScalars() const { return fAutoDetectScalars; }
	
	/// Sets the flag indicating if we should check for the scalars of
	/// software triggers automatically or not. If set to true then the
	/// scalarEvent parameter of the Decode method is ignored.
	void AutoDetectScalars(bool value) { fAutoDetectScalars = value; }
	
	/// Returns the maximum regional structure count expected in the DDL payload.
	UInt_t MaxRegionals() const { return fMaxRegionals; }
	
	/// Sets the maximum regional structure count expected in the DDL payload.
	void MaxRegionals(UInt_t n) { fMaxRegionals = n; }
	
	/// Returns the number of regional structures we actually attempted to decode.
	UInt_t RegionalsDecoded() const { return fNoRegionals; }
	
	/// Returns the maximum local structure count expected in any given regional
	/// card structure within the DDL payload.
	UInt_t MaxLocals() const { return fMaxLocals; }
	
	/// Sets the maximum local structure count expected in any given regional
	/// card structure within the DDL payload.
	void MaxLocals(UInt_t n) { fMaxLocals = n; }
	
	/// This method decodes the DDL payload contained in the buffer.
	bool Decode(const void* buffer, UInt_t bufferSize, bool scalarEvent = false);
	
	/// Returns the end of DARC marker key.
	static UInt_t EndOfDarcWord() { return fgkEndOfDarc; }
	
	/// Returns the end of global header marker key.
	static UInt_t EndOfGlobalWord() { return fgkEndOfGlobal; }
	
	/// Returns the end of regional structure marker key.
	static UInt_t EndOfRegionalWord() { return fgkEndOfReg; }
	
	/// Returns the regional error word.
	static UInt_t RegionalErrorWord() { return fgkErrorWord; }
	
	/// Returns the end of local structure marker key.
	static UInt_t EndOfLocalWord() { return fgkEndOfLocal; }
	
	/// Returns the local card disable word.
	static UInt_t LocalDisableWord() { return fgkDisableWord; }
	
	/// Returns value of default DARC type.
	static UInt_t DarcDefaultType() { return fgkDarcDefaultType; }
	
	/// Returns value of Vadorh DARC type.
	static UInt_t DarcVadorhType()  { return fgkDarcVadorhType; }
	
private:

	bool fExitOnError; ///< Indicates if we should exit on the very first error.
	bool fTryRecover; ///< Indicates if we should try recover from a corrupt structures.
	bool fAutoDetectScalars; ///< Flag to indicate if we should auto-detect if there are scalars in the data.
	bool fHadError; ///< Indicates if we had an error decoding the data.
	UInt_t fNoRegionals; ///< The number of regional card structures actually decoded.
	UInt_t fMaxRegionals; ///< Maximum number of regional card structures allowed in a DDL stream.
	UInt_t fMaxLocals; ///< Maximum number of local card structures per regional structure allowed in a DDL stream.
	EventHandler fHandler; ///< The event handler which deals with the generated parsing events.

	void DecodeBuffer(const UChar_t* start, const UChar_t* end, bool scalarEvent);
	void DecodeRegionalStructs(const UChar_t* start, const UChar_t* end, bool scalarEvent);
	const UChar_t* DecodeLocalStructs(const UChar_t* start, const UChar_t* end, bool scalarEvent);
	
	const UChar_t* FindKey(
			UInt_t key, const UChar_t* start, const UChar_t* end
		);

	static const UInt_t fgkEndOfDarc;   ///< Indicates the end of the DARC header.
	static const UInt_t fgkEndOfGlobal; ///< Indicates the end of the global header just after the DARC header.
	static const UInt_t fgkEndOfReg;    ///< Indicates the end of a regional card structure.
	static const UInt_t fgkErrorWord;   ///< The error word when a regional board is missing.
	static const UInt_t fgkEndOfLocal;  ///< Indicates the end of a local card structure.
	static const UInt_t fgkDisableWord; ///< Word used to fill "empty" slots.
	static const UInt_t fgkDarcDefaultType;   ///< default type for DARC def.
	static const UInt_t fgkDarcVadorhType;    ///< default type for DARC vadorh
};

//_____________________________________________________________________________

// The following are the structure key words which are used to demarcate (identify
// the end) of the structures within the raw trigger data.
template <class EventHandler>
const UInt_t AliMUONTriggerDDLDecoder<EventHandler>::fgkEndOfDarc    = 0xDEADFACE;
template <class EventHandler>
const UInt_t AliMUONTriggerDDLDecoder<EventHandler>::fgkEndOfGlobal  = 0xDEADBEEF;
template <class EventHandler>
const UInt_t AliMUONTriggerDDLDecoder<EventHandler>::fgkEndOfReg     = 0xBEEFFACE;
template <class EventHandler>
const UInt_t AliMUONTriggerDDLDecoder<EventHandler>::fgkErrorWord    = 0xCAFEDEAD;
template <class EventHandler>
const UInt_t AliMUONTriggerDDLDecoder<EventHandler>::fgkEndOfLocal   = 0xCAFEFADE;
template <class EventHandler>
const UInt_t AliMUONTriggerDDLDecoder<EventHandler>::fgkDisableWord  = 0x10CADEAD;
template <class EventHandler>
const UInt_t AliMUONTriggerDDLDecoder<EventHandler>::fgkDarcDefaultType = 0x6;
template <class EventHandler>
const UInt_t AliMUONTriggerDDLDecoder<EventHandler>::fgkDarcVadorhType  = 0x4;


template <class EventHandler>
bool AliMUONTriggerDDLDecoder<EventHandler>::Decode(
		const void* buffer, UInt_t bufferSize, bool scalarEvent
	)
{
	/// This method should be called to actually decode the DDL payload
	/// contained in a memory buffer. The payload should be in the muon trigger
	/// system DDL stream format otherwise it will be recognised as corrupt.
	/// As the decoder progresses it will make method calls to the event handler
	/// instance (which can be accessed with the GetHandler() method) to indicate
	/// the start of the DARC header, global header, new regional structures and
	/// local card structures.
	///
	/// If an error occurs during the parse, because the data is corrupt, then
	/// the OnError method is called indicating what the problem was.
	/// Decoding will stop at this point unless the fExitOnError flag is set
	/// to false. There is an optional flag fTryRecover which can enable logic
	/// which will attempt to recover from corruption of the data structures
	/// in the DDL payload, if they are found to be inconsistent (assumed corrupt).
	/// The recovery procedure simply involves trying to find in the DDL stream
	/// the location of the expected end of header/structure marker/key or the
	/// next expected key, and then continue decoding from there. This kind of
	/// recovery logic is (and should be) turned off by default and only
	/// enabled when trying to recover corrupted data.
	///
	/// \param buffer  This is the pointer to the start of the memory buffer
	///     containing the DDL payload. Remember that this must be the start of
	///     the payload and not the DDL stream. That is, this pointer should be
	///     equal to: DDL start pointer + 8 * sizeof(UInt_t).
	/// \param bufferSize  This is the pointer to the first byte just past the
	///     end of the block structure.
	/// \param scalarEvent  Set to true if this DDL contains a scalar event
	///     and false if it is a normal physics event. If the fAutoDetectScalars
	///     flag is true then we ignore this parameter.
	/// \return Returns false if there was any problem with decoding the data,
	///     and true otherwise. Note: the data may have been partially decoded
	///     even if false was returned, which would be indicated by at least one
	///     call to the event handler's OnLocalStruct or OnNewRegionalStruct methods.
	
	assert( buffer != NULL );
	
	fHadError = false;
	
	// We are basically implementing something like a recursive decent parser.
	// So start by marking the start of buffer position and end of buffer.
	const UChar_t* start = reinterpret_cast<const UChar_t*>(buffer);
	const UChar_t* end = start + bufferSize;
	
	fHandler.OnNewBuffer(buffer, bufferSize);
	DecodeBuffer(start, end, scalarEvent);
	fHandler.OnEndOfBuffer(buffer, bufferSize);
	return not fHadError;
}


template <class EventHandler>
void AliMUONTriggerDDLDecoder<EventHandler>::DecodeBuffer(
		const UChar_t* start, const UChar_t* end, bool scalarEvent
	)
{
	/// This method does the work to decode the buffer's payload data. It
	/// unpacks the DARC header, the global header, then calls the
	/// DecodeRegionalStructs() method to decode the regional structures.
	/// For the DARC and global headers the methods OnDarcHeader and
	/// OnGlobalHeader of the event handler object are called respectively.
	/// \param start  This is the pointer to the start of the buffer to decode.
	/// \param end  This is the pointer to the first byte just past the
	///             end of the buffer.
	/// \param scalarEvent  Set to true if this DDL contains a scalar event
	///     and false if it is a normal physics event.
	///
	/// fHadError is set to true if there were any errors decoding the buffer
	/// and the OnError method of the callback event handler is called for
	/// each error.
	
	const UChar_t* current = start;
	
	// Mark the DARC header, but check that we do not overrun the buffer.
	const UInt_t* darcHeader = reinterpret_cast<const UInt_t*>(current);
	current += sizeof(UInt_t);
	if (current > end or current < start)
	{
		// Indicate we had an error and stop the decoding because we
		// hit the end of the buffer already.
		fHandler.OnError(EventHandler::kNoDarcHeader, darcHeader);
		fHadError = true;
		return;
	}
	
	// Check if we need to figure out if this is a scalar event.
	// If we do, then we do this by checking the event type in the DARC header
	// and double checking this by checking if the fgkEndOfDarc key is in the
	// expected location for a scalar event DDL payload.
	if (fAutoDetectScalars)
	{
		const UInt_t* expectedEndOfDarc = reinterpret_cast<const UInt_t*>(
				current + sizeof(AliMUONDarcScalarsStruct)
			);
		
		// First make sure not to read past the end of buffer. Then check
		// the value of the event type in the DARC header.
		// Physics events are indicated by the two trigger bits of the
		// DARC header set to 01b. Everything else is a software trigger
		// of some sort.
		if (reinterpret_cast<const UChar_t*>(expectedEndOfDarc+1) <= end and
		    reinterpret_cast<const UChar_t*>(expectedEndOfDarc+1) > start and
		    EventHandler::GetDarcEventType(*darcHeader) != 0x1
		    and *expectedEndOfDarc == fgkEndOfDarc
		   )
		{
			scalarEvent = true;
		}
		else
		{
			scalarEvent = false;
		}
	}
	
	// Detect how many regional blocks we expect. If we have no idea then
	// just use what the maximum setting is.
	UInt_t darkType = EventHandler::GetDarcType(*darcHeader);
	if (darkType == fgkDarcVadorhType)
	{
		fNoRegionals = 1;
	}
	else if (darkType == fgkDarcDefaultType)
	{
		fNoRegionals = 8;
	}
	else
	{
		fNoRegionals = fMaxRegionals;
	}
	
	// Check if the DARC header indicates we expect more regionals than we
	// are allowed to decode according to our max regional structures count.
	// If we do then this is an error and we should indicate it and exit
	// if so requested by the user. Also we can fix the number of regionals
	// to expect if we are trying to recover from errors.
	if (fNoRegionals > fMaxRegionals)
	{
		fHandler.OnError(EventHandler::kTooManyRegionals, darcHeader);
		fHadError = true;
		if (fExitOnError) return;
		if (fTryRecover)
		{
			fNoRegionals = fMaxRegionals;
		}
	}
	
	// Check that the DARC header indicates correctly if we are a scalar event or not.
	bool darcShowsScalar = (EventHandler::GetDarcEventType(*darcHeader) != 0x1);
	if (darcShowsScalar != scalarEvent)
	{
		// Indicate we had an error and stop the decoding if so requested
		// by the user.
		fHandler.OnError(EventHandler::kWrongEventType, darcHeader);
		fHadError = true;
		if (fExitOnError) return;
	}
	
	// Decode the DARC scalars if this is a scalar event.
	const AliMUONDarcScalarsStruct* darcScalars = NULL;
	if (scalarEvent)
	{
		darcScalars = reinterpret_cast<const AliMUONDarcScalarsStruct*>(current);
		current += sizeof(AliMUONDarcScalarsStruct);
		if (current > end or current < start)
		{
			// If we overflowed the pointer and already had an error then
			// we are clearly lost so just stop decoding before we segfault.
			if (current < start and fHadError) return;
			
			// Indicate we had an error and stop the decoding because we
			// hit the end of the buffer already.
			fHandler.OnError(EventHandler::kNoDarcScalars, darcScalars);
			fHadError = true;
			return;
		}
	}
	
	// Now check that the end of DARC header marker is OK.
	const UInt_t* endOfDarc = reinterpret_cast<const UInt_t*>(current);
	current += sizeof(UInt_t);
	if (current > end or current < start)
	{
		// If we overflowed the pointer and already had an error then
		// we are clearly lost so just stop decoding before we segfault.
		if (current < start and fHadError) return;
		
		// Indicate we had an error and stop the decoding because we
		// hit the end of the buffer already.
		fHandler.OnError(EventHandler::kNoEndOfDarc, endOfDarc);
		fHadError = true;
		return;
	}
	if (*endOfDarc != fgkEndOfDarc)
	{
		// Indicate we had an error and stop the decoding if so requested
		// by the user.
		fHandler.OnError(EventHandler::kBadEndOfDarc, endOfDarc);
		fHadError = true;
		if (fExitOnError) return;
		
		// If the user requested for us to try and recover from structure
		// errors then we need to try locate the key in the data stream
		// and continue decoding from there.
		if (fTryRecover)
		{
			const UChar_t* keypos = FindKey(fgkEndOfDarc,
					reinterpret_cast<const UChar_t*>(darcHeader),
					end
				);
			if (keypos != NULL)
			{
				// remember to continue decoding just past the key.
				current = keypos + sizeof(UInt_t);
			}
		}
	}
	
	fHandler.OnDarcHeader(*darcHeader, darcScalars, current);
	
	// Next, we mark the Global header and check that we do not overrun the buffer.
	const AliMUONGlobalHeaderStruct* globalHeader =
		reinterpret_cast<const AliMUONGlobalHeaderStruct*>(current);
	current += sizeof(AliMUONGlobalHeaderStruct);
	if (current > end or current < start)
	{
		// If we overflowed the pointer and already had an error then
		// we are clearly lost so just stop decoding before we segfault.
		if (current < start and fHadError) return;
		
		// Indicate we had an error and stop the decoding because we
		// hit the end of the buffer already.
		fHandler.OnError(EventHandler::kNoGlobalHeader, globalHeader);
		fHadError = true;
		return;
	}
	
	// Decode the global scalars if this is a scalar event.
	const AliMUONGlobalScalarsStruct* globalScalars = NULL;
	if (scalarEvent)
	{
		globalScalars = reinterpret_cast<const AliMUONGlobalScalarsStruct*>(current);
		current += sizeof(AliMUONGlobalScalarsStruct);
		if (current > end or current < start)
		{
			// If we overflowed the pointer and already had an error then
			// we are clearly lost so just stop decoding before we segfault.
			if (current < start and fHadError) return;
			
			// Indicate we had an error and stop the decoding because we
			// hit the end of the buffer already.
			fHandler.OnError(EventHandler::kNoGlobalScalars, globalScalars);
			fHadError = true;
			return;
		}
	}
	
	// Now check that the end of global header marker is OK.
	const UInt_t* endOfGlobal = reinterpret_cast<const UInt_t*>(current);
	current += sizeof(UInt_t);
	if (current > end or current < start)
	{
		// If we overflowed the pointer and already had an error then
		// we are clearly lost so just stop decoding before we segfault.
		if (current < start and fHadError) return;
		
		// Indicate we had an error and stop the decoding because we
		// hit the end of the buffer already.
		fHandler.OnError(EventHandler::kNoEndOfGlobal, endOfGlobal);
		fHadError = true;
		return;
	}
	if (*endOfGlobal != fgkEndOfGlobal)
	{
		// Indicate we had an error and stop the decoding if so requested
		// by the user.
		fHandler.OnError(EventHandler::kBadEndOfGlobal, endOfGlobal);
		fHadError = true;
		if (fExitOnError) return;
		
		// If the user requested for us to try and recover from structure
		// errors then we need to try locate the key in the data stream
		// and continue decoding from there.
		if (fTryRecover)
		{
			const UChar_t* keypos = FindKey(fgkEndOfGlobal,
					reinterpret_cast<const UChar_t*>(globalHeader),
					end
				);
			if (keypos != NULL)
			{
				// remember to continue decoding just past the key.
				current = keypos + sizeof(UInt_t);
			}
		}
	}
	
	fHandler.OnGlobalHeader(globalHeader, globalScalars, current);
	
	DecodeRegionalStructs(current, end, scalarEvent);
}


template <class EventHandler>
void AliMUONTriggerDDLDecoder<EventHandler>::DecodeRegionalStructs(
		const UChar_t* start, const UChar_t* end, bool scalarEvent
	)
{
	/// This method decodes the regional structures in the DDL data.
	/// For each regional header found, the OnNewRegionalStruct method of
	/// the event handler is called, to signal the start of each
	/// new regional structure. The DecodeLocalStructs() method is then
	/// called to decode the local trigger structures. Finally, a symmetrical
	/// call to OnEndOfRegionalStruct of the event handler is called to
	/// signal that the regional structure has been decoded.
	/// \param start  This is the pointer to the start of the buffer.
	/// \param end  This is the pointer to the first byte just past the
	///             end of the buffer.
	/// \param scalarEvent  Set to true if this DDL contains a scalar event
	///     and false if it is a normal physics event.
	///
	/// fHadError is set to true if there were any errors decoding the buffer
	/// and the OnError method of the callback event handler is called for
	/// each error.
	
	const UChar_t* current = start;
	
	for (UInt_t iReg = 0; iReg < fNoRegionals; iReg++)
	{
		const AliMUONRegionalHeaderStruct* regionalHeader =
			reinterpret_cast<const AliMUONRegionalHeaderStruct*>(current);
		current += sizeof(AliMUONRegionalHeaderStruct);
		
		if (current > end or current < start)
		{
			// If we overflowed the pointer and already had an error then
			// we are clearly lost so just stop decoding before we segfault.
			if (current < start and fHadError) return;
			
			// So we only got part of a regional header, nothing to do but
			// report the error and exit.
			fHandler.OnError(EventHandler::kNoRegionalHeader, regionalHeader);
			fHadError = true;
			return;
		}
		
		// Skip empty regional board (not connected or with error reading).
		if (regionalHeader->fDarcWord == fgkErrorWord)
		{
			//current += sizeof(AliMUONRegionalHeaderStruct);  // already done above
			if (scalarEvent)
				current += sizeof(AliMUONRegionalScalarsStruct);
			current += sizeof(UInt_t); // skip the end of regional structure key.
			
			// now also skip the local structure part:
			current += fMaxLocals * sizeof(AliMUONLocalInfoStruct);
			if (scalarEvent)
				current += fMaxLocals * sizeof(AliMUONLocalScalarsStruct);
			current += fMaxLocals * sizeof(UInt_t); // skip all the end of local structure keys.
			
			continue;
		}
		
		// Decode the regional scalar words if this is a scalar event.
		const AliMUONRegionalScalarsStruct* regionalScalars = NULL;
		if (scalarEvent)
		{
			regionalScalars = reinterpret_cast<const AliMUONRegionalScalarsStruct*>(current);
			current += sizeof(AliMUONRegionalScalarsStruct);
			if (current > end or current < start)
			{
				// If we overflowed the pointer and already had an error then
				// we are clearly lost so just stop decoding before we segfault.
				if (current < start and fHadError) return;
				
				// Indicate we had an error and stop the decoding because we
				// hit the end of the buffer already.
				fHandler.OnError(EventHandler::kNoRegionalScalars, regionalScalars);
				fHadError = true;
				return;
			}
		}
		
		// Now check that the end of regional header marker is OK.
		const UInt_t* endOfRegional = reinterpret_cast<const UInt_t*>(current);
		current += sizeof(UInt_t);
		if (current > end or current < start)
		{
			// If we overflowed the pointer and already had an error then
			// we are clearly lost so just stop decoding before we segfault.
			if (current < start and fHadError) return;
			
			// Indicate we had an error and stop the decoding because we
			// hit the end of the buffer already.
			fHandler.OnError(EventHandler::kNoEndOfRegional, endOfRegional);
			fHadError = true;
			return;
		}
		if (*endOfRegional != fgkEndOfReg)
		{
			// Indicate we had an error and stop the decoding if so requested
			// by the user.
			fHandler.OnError(EventHandler::kBadEndOfRegional, endOfRegional);
			fHadError = true;
			if (fExitOnError) return;
			
			// If the user requested for us to try and recover from structure
			// errors then we need to try locate the key in the data stream
			// and continue decoding from there.
			if (fTryRecover)
			{
				const UChar_t* keypos = FindKey(fgkEndOfReg,
						reinterpret_cast<const UChar_t*>(regionalHeader),
						end
					);
				
				// If the fgkEndOfReg key was found exactly one regional
				// structure later then we should not adjust the current
				// decoding position because it is more likely that the
				// end of regional structure key was just corrupt rather
				// than offset.
				// If we could not find another good key then just continue.
				size_t sizeOfRegional = sizeof(AliMUONRegionalHeaderStruct) + sizeof(UInt_t)
					+ fMaxLocals * (sizeof(AliMUONLocalInfoStruct) + sizeof(UInt_t));
				if (scalarEvent)
				{
					sizeOfRegional += sizeof(AliMUONRegionalScalarsStruct)
						+ fMaxLocals * sizeof(AliMUONLocalScalarsStruct);
				}
				
				if (keypos != NULL and keypos != current + sizeOfRegional)
				{
					current = keypos + sizeof(UInt_t);
				}
			}
		}
		
		// Tell the handler that we have a new regional block and decode it.
		// When done, tell the handler again.
		// We call both versions of OnNewRegionalStruct so that user event
		// handlers can implement the one they prefer to use.
		const UChar_t* startOfLocals = current;
		fHandler.OnNewRegionalStruct(regionalHeader, regionalScalars, startOfLocals);
		fHandler.OnNewRegionalStructV2(iReg, regionalHeader, regionalScalars, startOfLocals);
		current = DecodeLocalStructs(current, end, scalarEvent);
		fHandler.OnEndOfRegionalStruct(regionalHeader, regionalScalars, startOfLocals);
		fHandler.OnEndOfRegionalStructV2(iReg, regionalHeader, regionalScalars, startOfLocals);
	}
	
	// Now just check that there is no extra rubbish at the end of the DDL.
	if (current != end)
	{
		fHandler.OnError(EventHandler::kBufferTooBig, current);
		fHadError = true;
	}
}


template <class EventHandler>
const UChar_t* AliMUONTriggerDDLDecoder<EventHandler>::DecodeLocalStructs(
		const UChar_t* start, const UChar_t* end, bool scalarEvent
	)
{
	/// This method decodes the local structures in the DDL data for a
	/// single regional structure. For each local trigger structure found,
	/// the OnLocalStruct method of the event handler is called.
	/// \param start  This is the pointer to the start of the regional structure
	///               payload (The pointer just past the regional header key).
	/// \param end  This is the pointer to the first byte just past the
	///             end of the buffer.
	/// \param scalarEvent  Set to true if this DDL contains a scalar event
	///     and false if it is a normal physics event.
	/// \returns  The position in the buffer where this method stopped decoding.
	///
	/// fHadError is set to true if there were any errors decoding the buffer
	/// and the OnError method of the callback event handler is called for
	/// each error.
	
	const UChar_t* current = start;
	
	for (UInt_t iLocal = 0; iLocal < fMaxLocals; iLocal++)
	{
		const AliMUONLocalInfoStruct* localStruct =
			reinterpret_cast<const AliMUONLocalInfoStruct*>(current);
		current += sizeof(AliMUONLocalInfoStruct);
		
		if (current > end or current < start)
		{
			// If we overflowed the pointer and already had an error then
			// we are clearly lost so just stop decoding before we segfault.
			if (current < start and fHadError) return end;
			
			// So we only got part of a local structure, nothing to do but
			// report the error and exit.
			fHandler.OnError(EventHandler::kNoLocalStruct, localStruct);
			fHadError = true;
			return end;
		}
		
		// Skip empty local board if card not notified.
		if (localStruct->fX2X1 == fgkDisableWord and
		    localStruct->fX4X3 == fgkDisableWord and
		    localStruct->fY2Y1 == fgkDisableWord and
		    localStruct->fY4Y3 == fgkDisableWord and
		    localStruct->fTriggerBits == fgkDisableWord
		   )
		{
			//current += sizeof(AliMUONLocalInfoStruct); // already done above
			if (scalarEvent)
				current += sizeof(AliMUONLocalScalarsStruct);
			current += sizeof(UInt_t); // skip the end of local structure key.
			continue;
		}
		
		// Decode the regional scalar words if this is a scalar event.
		const AliMUONLocalScalarsStruct* localScalars = NULL;
		if (scalarEvent)
		{
			localScalars = reinterpret_cast<const AliMUONLocalScalarsStruct*>(current);
			current += sizeof(AliMUONLocalScalarsStruct);
			if (current > end or current < start)
			{
				// If we overflowed the pointer and already had an error then
				// we are clearly lost so just stop decoding before we segfault.
				if (current < start and fHadError) return end;
				
				// Indicate we had an error and stop the decoding because we
				// hit the end of the buffer already.
				fHandler.OnError(EventHandler::kNoLocalScalars, localScalars);
				fHadError = true;
				return end;
			}
		}
		
		// Now check that the end of regional header marker is OK.
		const UInt_t* endOfLocal = reinterpret_cast<const UInt_t*>(current);
		current += sizeof(UInt_t);
		if (current > end or current < start)
		{
			// If we overflowed the pointer and already had an error then
			// we are clearly lost so just stop decoding before we segfault.
			if (current < start and fHadError) return end;
			
			// Indicate we had an error and stop the decoding because we
			// hit the end of the buffer already. We can however signal that
			// we got a local structure, because the buffer contains enough
			// data to potencially contain the real data of the structure.
			fHandler.OnError(EventHandler::kNoEndOfLocal, endOfLocal);
			if (not fExitOnError)
			{
				fHandler.OnLocalStruct(localStruct, localScalars);
			}
			fHadError = true;
			return end;
		}
		if (*endOfLocal != fgkEndOfLocal)
		{
			// Indicate we had an error and stop the decoding if so requested
			// by the user.
			fHandler.OnError(EventHandler::kBadEndOfLocal, endOfLocal);
			fHadError = true;
			if (fExitOnError) return current;
			
			// If the user requested for us to try and recover from structure
			// errors then we need to try locate the key in the data stream
			// and continue decoding from there.
			if (fTryRecover)
			{
				const UChar_t* searchPos = reinterpret_cast<const UChar_t*>(localStruct);
				const UChar_t* firstLocalKey = FindKey(fgkEndOfLocal, searchPos, end);
				const UChar_t* firstRegionalKey = FindKey(fgkEndOfReg, searchPos, end);
				
				// If a regional key was found first, then give up on
				// anymore local structures from this regional block and
				// continue decoding from the next regional block.
				// Also if the fgkEndOfLocal key was found exactly one
				// local structure later then we should not adjust the
				// current decoding position because it is more likely that
				// the end of local structure key was just corrupt rather
				// than offset.
				if (firstLocalKey != NULL and firstRegionalKey != NULL)
				{
					if (firstLocalKey < firstRegionalKey)
					{
						size_t sizeOflocalStruct = sizeof(AliMUONLocalInfoStruct) + sizeof(UInt_t);
						if (scalarEvent)
							sizeOflocalStruct += sizeof(AliMUONLocalScalarsStruct);
						
						if (firstLocalKey != current + sizeOflocalStruct)
							current = firstLocalKey + sizeof(UInt_t);
					}
					else
					{
						// Adjust back to the start of the regional header.
						current = firstRegionalKey - sizeof(AliMUONRegionalHeaderStruct);
						if (scalarEvent)
							current -= sizeof(AliMUONRegionalScalarsStruct);
						break;
					}
				}
				
				// If we could not find another good key then just continue.
			}
		}
		
		// Call both handlers for local structures so that the user decoder event
		// handler class can implement the one it prefers to use.
		fHandler.OnLocalStruct(localStruct, localScalars);
		fHandler.OnLocalStructV2(iLocal, localStruct, localScalars);
	}
	
	return current;
}


template <class EventHandler>
const UChar_t* AliMUONTriggerDDLDecoder<EventHandler>::FindKey(
		UInt_t key, const UChar_t* start, const UChar_t* end
	)
{
	/// Searches for the first occurrence of the key value in the buffer
	/// marked by 'start' and 'end'. 'start' should point to the start of
	/// the buffer and 'end' should point to 'start + bufferSize'.
	/// \param key  The 32 bit word to look for.
	/// \param start  The start location to begin the search at.
	/// \param end  The pointer to the first byte just past end of the
	///             buffer to check.
	/// \returns  The location of the first occurance of key in the buffer,
	///           otherwise NULL is returned if nothing was found.
 
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

#endif // ALIMUONTRIGGERDDLDECODER_H
