#ifndef ALIMUONVRAWSTREAMTRACKER_H
#define ALIMUONVRAWSTREAMTRACKER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

///
/// \file   AliMUONVRawStreamTracker.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   28-11-2007
/// \brief  Declaration of the abstract base class for muon trigger chamber raw stream decoders.
///

#include "AliMUONRawStream.h"

class AliMUONVRawStreamTracker : public AliMUONRawStream
{
public:

	/// Default constructor.
	AliMUONVRawStreamTracker();
	
	/// Constructor setting the raw reader.
	AliMUONVRawStreamTracker(AliRawReader* rawReader);
	
	/// Default destructor.
	virtual ~AliMUONVRawStreamTracker();
	
	/// Advance one step in the iteration. Returns false if finished.
	virtual Bool_t Next(Int_t& busPatchId,
				UShort_t& manuId, UChar_t& manuChannel,
				UShort_t& adc) = 0;
	
	/// Class used in the following Next() method to return blocks of decoded
	/// channel data. This is better because we generate a lot fewer method calls.
	class AliChannelInfo
	{
	public:
		/// Default constructor.
		AliChannelInfo(Int_t busPatchId = 0, UShort_t manuId = 0, UChar_t channelId = 0, UShort_t adc = 0) :
			fBusPatchId(busPatchId), fManuId(manuId), fChannelId(channelId), fADC(adc)
		{}
		
		/// Returns the bus patch ID.
		Int_t BusPatchId() const { return fBusPatchId; }
		/// Returns the MANU ID.
		UShort_t ManuId() const { return fManuId; }
		/// Returns the channel ID.
		UShort_t ChannelId() const { return fChannelId; }
		/// ADC signal.
		UShort_t ADC() const { return fADC; }
	
	private:
		Int_t fBusPatchId;  //!< The bus patch ID for this channel.
		UShort_t fManuId;   //!< MANU ID.
		UChar_t fChannelId; //!< MANU channel ID.
		UShort_t fADC;      //!< ADC signal.
	};
	
	/// Returns the next batch of decoded channel data.
	/// [out] \param channels  This is filled with the pointer to the array
	///                        containing the channel information.
	/// \return The number of elements in the 'channels' array is returned.
	///     Zero is returned if there are no more digits to be fetched.
	virtual UInt_t Next(const AliChannelInfo*& channels) = 0;
	
	/// Return maximum number of DDLs
	static Int_t GetMaxDDL() { return fgkMaxDDL; };
	
	/// Return maximum number of blocks per DDL allowed.
	virtual Int_t GetMaxBlock() const = 0;
	/// Return maximum number of Dsp per block allowed.
	virtual Int_t GetMaxDsp()   const = 0;
	/// Return maximum number of Buspatch per Dsp allowed.
	virtual Int_t GetMaxBus()   const = 0;
	
	/// Set maximum number of blocks per DDL allowed.
	virtual void SetMaxBlock(Int_t blk) = 0;
	/// Set maximum number of Dsp per block allowed.
	virtual void SetMaxDsp(Int_t dsp) = 0;
	/// Set maximum number of Buspatch per Dsp allowed.
	virtual void SetMaxBus(Int_t bus) = 0;
	
	/// Return number of the current DDL.
	virtual Int_t GetDDL() const = 0;
	
	/// check error/Warning presence
	virtual Bool_t IsErrorMessage() const = 0;
	
	/// error numbers
	enum rawStreamTrackerError
	{
		kGlitchErr      = 1, ///< glitch error
		kPaddingWordErr = 2, ///< padding word error
		kParityErr      = 3  ///< parity error
	};

private:

	// Do not allow copying of this class.
        /// Not implemented
	AliMUONVRawStreamTracker(const AliMUONVRawStreamTracker& stream);
        /// Not implemented
	AliMUONVRawStreamTracker& operator = (const AliMUONVRawStreamTracker& stream);

	static const Int_t fgkMaxDDL;   //!< maximum number of DDLs

	ClassDef(AliMUONVRawStreamTracker, 0) // Base class for reading MUON raw digits from tracking chambers.
};

#endif  // ALIMUONVRAWSTREAMTRACKER_H
