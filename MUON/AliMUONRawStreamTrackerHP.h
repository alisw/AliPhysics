#ifndef ALIMUONRAWSTREAMTRACKERHP_H
#define ALIMUONRAWSTREAMTRACKERHP_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

///
/// \file   AliMUONRawStreamTrackerHP.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   29-11-2007
/// \brief  Declaration of the high performance decoder for muon trigger chamber raw streams.
///

#include "AliMUONVRawStreamTracker.h"
#include "AliMUONTrackerDDLDecoder.h"

class AliMUONRawStreamTrackerHP : public AliMUONVRawStreamTracker
{
public:

	/// Default constructor.
	AliMUONRawStreamTrackerHP();
	
	/// Constructor for setting the raw reader.
	AliMUONRawStreamTrackerHP(AliRawReader* rawReader);
	
	/// Default destructor.
	virtual ~AliMUONRawStreamTrackerHP();
	
	// The following public methods are all inherited from AliMUONVRawStreamTracker:
	
	/// Initialize iterator
	virtual void First();
	
	/// DDL iterator 
	virtual Bool_t NextDDL();
	
	/// Whether the iteration is finished or not
	virtual Bool_t IsDone() const;
	
	/// Nothing is actually done in the AddErrorMessage method because we log
	/// the error messages as we find them in AliDecoderEventHandler::OnError().
	virtual void AddErrorMessage() {};
	
	/// Advance one step in the iteration. Returns false if finished.
	virtual Bool_t Next(Int_t& busPatchId,
				UShort_t& manuId, UChar_t& manuChannel,
				UShort_t& adc);
	
	/// Returns the next batch of decoded channel data.
	virtual UInt_t Next(const AliChannelInfo*& channels);
	
	/// Return maximum number of blocks per DDL allowed.
	virtual Int_t GetMaxBlock() const { return (Int_t) fDecoder.MaxBlocks(); }
	/// Return maximum number of Dsp per block allowed.
	virtual Int_t GetMaxDsp() const { return (Int_t) fDecoder.MaxDSPs(); }
	/// Return maximum number of Buspatch per Dsp allowed.
	virtual Int_t GetMaxBus() const { return (Int_t) fDecoder.MaxBusPatches(); }
	
	/// Set maximum number of blocks per DDL allowed.
	virtual void SetMaxBlock(Int_t blk) { fDecoder.MaxBlocks( (UInt_t) blk ); }
	/// Set maximum number of Dsp per block allowed.
	virtual void SetMaxDsp(Int_t dsp) { fDecoder.MaxDSPs( (UInt_t) dsp ); }
	/// Set maximum number of Buspatch per Dsp allowed.
	virtual void SetMaxBus(Int_t bus) { fDecoder.MaxBusPatches( (UInt_t) bus ); }
	
	/// Return number of the current DDL.
	virtual Int_t GetDDL() const { return fDDL - 1; }
	
	/// check error/Warning presence
	virtual Bool_t IsErrorMessage() const { return fHadError; }

private:

	// Do not allow copying of this class.
        /// Not implemented
	AliMUONRawStreamTrackerHP(const AliMUONRawStreamTrackerHP& stream);
        /// Not implemented
	AliMUONRawStreamTrackerHP& operator = (const AliMUONRawStreamTrackerHP& stream);
	
	/// This is the custom event handler (callback interface) class which
	/// unpacks raw data words and fills an internal buffer with decoded digits
	/// as they are decoded by the high performance decoder.
	/// Any errors are logged to the parent AliMUONVRawStreamTracker, so one
	/// must set this pointer appropriately before decoding and DDL payload.
	class AliDecoderEventHandler : public AliMUONTrackerDDLDecoderEventHandler
	{
	public:
	
		/// Default constructor.
		AliDecoderEventHandler();
		/// Default destructor.
		virtual ~AliDecoderEventHandler();
		
		/// Sets the raw stream object which should be the parent of this class.
		void SetRawStream(AliMUONVRawStreamTracker* rawStream) { fRawStream = rawStream; }
		
		/// Return the number of channels in the buffer returned by Channels().
		UInt_t ChannelCount() const { return fChannelCount; }
		
		/// Return the buffer of decoded channel data.
		const AliChannelInfo* Channels() const { return fChannelBuffer; }
		
		// The following methods are inherited from AliMUONTrackerDDLDecoderEventHandler:
		
		/// New buffer handler.
		void OnNewBuffer(const void* buffer, UInt_t bufferSize);
		
		/// New bus patch handler.
		void OnNewBusPatch(const AliMUONBusPatchHeaderStruct* header, const void* data);
		
		/// Raw data word handler.
		void OnData(UInt_t data);
		
		/// Error handler.
		void OnError(ErrorCode error, const void* location);
	
	private:
	
		// Do not allow copying of this class.
                /// Not implemented
		AliDecoderEventHandler(const AliDecoderEventHandler& /*obj*/);
                /// Not implemented
		AliDecoderEventHandler& operator = (const AliDecoderEventHandler& /*obj*/);

		Int_t fBusPatchId;     //!< The bus patch ID of the current bus patch being decoded.
		UInt_t fChannelCount;  //!< Number of elements in fChannelBuffer.
		UInt_t fMaxChannels;   //!< Maximum number of elements that can be stored in fChannelBuffer.
		AliChannelInfo* fChannelBuffer;  //!< Buffer of decoded channel structures.
		AliMUONVRawStreamTracker* fRawStream; //!< Pointer to the parent raw stream object.
		const void* fBufferStart;   //!< Pointer to the start of the current DDL payload buffer.
	};
	
	AliMUONTrackerDDLDecoder<AliDecoderEventHandler> fDecoder;  //!< The decoder for the DDL payload.
	Int_t fDDL;         //!< The current DDL number being handled.
	UInt_t fCurrentChannel;  //!< The current channel to return by Next().
	Int_t fBufferSize;  //!< This is the buffer size in bytes of fBuffer.
	UChar_t* fBuffer;   //!< This is the buffer in which we store the DDL payload read from AliRawReader.
	Bool_t fHadError;   //!< Flag indicating if there was a decoding error or not.

	ClassDef(AliMUONRawStreamTrackerHP, 0) // High performance decoder for reading MUON raw digits from tracking chamber DDL data.
};

#endif  // ALIMUONRAWSTREAMTRACKERHP_H
