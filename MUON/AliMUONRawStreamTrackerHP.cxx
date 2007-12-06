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

/* $Id$*/

///
/// \file   AliMUONRawStreamTrackerHP.cxx
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   29-11-2007
/// \brief  Implementation of the the high performance decoder AliMUONRawStreamTrackerHP.
///

//-----------------------------------------------------------------------------
/// \ingroup raw
/// \class AliMUONRawStreamTrackerHP
/// \brief A high performance stream decoder for muon tracking DDL streams.
///
/// This is the raw stream class which interfaces between the high performance
/// core decoder and the AliRawReader class.
/// To gain the most out of the decoder, the Next() method which returns batches
/// of decoded digit / channel information should be used. That is:
/// \code
///   UInt_t Next(const AliChannelInfo*& channels);
/// \endcode
///
/// \author Artur Szostak <artursz@iafrica.com>
//-----------------------------------------------------------------------------

#include "AliMUONRawStreamTrackerHP.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include <cassert>

/// \cond CLASSIMP
ClassImp(AliMUONRawStreamTrackerHP)
/// \endcond


AliMUONRawStreamTrackerHP::AliMUONRawStreamTrackerHP() :
	AliMUONVRawStreamTracker(),
	fDecoder(),
	fDDL(0),
	fCurrentChannel(0),
	fBufferSize(8192),
	fBuffer(new UChar_t[8192]),
	fHadError(kFALSE)
{
	///
	/// Default constructor.
	///
}


AliMUONRawStreamTrackerHP::AliMUONRawStreamTrackerHP(AliRawReader* rawReader) :
	AliMUONVRawStreamTracker(rawReader),
	fDecoder(),
	fDDL(0),
	fCurrentChannel(0),
	fBufferSize(8192),
	fBuffer(new UChar_t[8192]),
	fHadError(kFALSE)
{
	///
	/// Constructor with AliRawReader as argument.
	///
	
	fDecoder.GetHandler().SetRawStream(this);
}


AliMUONRawStreamTrackerHP::~AliMUONRawStreamTrackerHP()
{
	///
	/// Default destructor.
	///
	
	if (fBuffer != NULL)
	{
		delete [] fBuffer;
	}
}


void AliMUONRawStreamTrackerHP::First()
{
	/// Initialise or reset the iterator.
	/// The first DDL will be found and decoded.
	
	assert( GetReader() != NULL );
	
	fDDL = 0;
	NextDDL();
}


Bool_t AliMUONRawStreamTrackerHP::NextDDL()
{
	/// Reading the next tracker DDL and decode the payload with the 
	/// high performance decoder.
	/// \return kTRUE if the next DDL was successfully read and kFALSE otherwise.

	assert( GetReader() != NULL );
	
	if (IsDone()) return kFALSE;
	
	do
	{
		GetReader()->Reset();
		GetReader()->Select("MUONTRK", fDDL, fDDL);  // Select the DDL file to be read.
		if (GetReader()->ReadHeader()) break;
		AliDebug(3, Form("Skipping DDL %d which does not seem to be there", fDDL+1));
		fDDL++;
	}
	while (fDDL < GetMaxDDL());
	
	AliDebug(3, Form("DDL Number %d\n", fDDL));
	
	Int_t dataSize = GetReader()->GetDataSize(); // in bytes
	// Check if we have enough buffer space already in fBuffer. If we do then
	// just continue reading otherwise we need to resize the buffer.
	if (fBufferSize < dataSize)
	{
		if (fBuffer != NULL)
		{
			delete [] fBuffer;
			fBuffer = NULL;
			fBufferSize = 0;
		}
		try
		{
			fBuffer = new UChar_t[dataSize];
			fBufferSize = dataSize;
		}
		catch (const std::bad_alloc&)
		{
			AliError("Could not allocate more buffer space. Cannot decode DDL.");
			return kFALSE;
		}
	}
	
	if (not GetReader()->ReadNext(fBuffer, dataSize))
	{
		return kFALSE;
	}
	
#ifndef R__BYTESWAP
	Swap(fBuffer, dataSize); // Swap needed for mac power pc.
#endif
	
	bool result = false;
	try
	{
		// Since we might allocate memory inside OnNewBuffer in the event
		// handler we need to trap any memory allocation exception to be robust.
		result = fDecoder.Decode(fBuffer, dataSize);
	}
	catch (const std::bad_alloc&)
	{
		AliError("Could not allocate more buffer space. Cannot decode DDL.");
		return kFALSE;
	}
	
	fCurrentChannel = 0;
	fDDL++; // Remember to increment index to next DDL.
	return result;
}


Bool_t AliMUONRawStreamTrackerHP::IsDone() const
{
	/// Indicates whether the iteration is finished or not.
	/// \return kTRUE if we already read all the digits and kFALSE if not.
	
	return fDDL == GetMaxDDL();
}


Bool_t AliMUONRawStreamTrackerHP::Next(
		Int_t& busPatchId, UShort_t& manuId, UChar_t& manuChannel,
		UShort_t& adc
	)
{
	/// Advance one step in the iteration. Returns false if finished.
	/// [out] \param busPatchId  This is filled with the bus patch ID of the digit.
	/// [out] \param manuId      This is filled with the MANU ID of the digit.
	/// [out] \param manuChannel This is filled with the MANU channel ID of the digit.
	/// [out] \param adc         This is filled with the ADC signal value of the digit.
	/// \return kTRUE if we read another digit and kFALSE if we have read all the
	///    digits already, i.e. at the end of the iteration.
	
	if (fCurrentChannel < fDecoder.GetHandler().ChannelCount())
	{
		const AliChannelInfo& channel = fDecoder.GetHandler().Channels()[fCurrentChannel];
		busPatchId = channel.BusPatchId();
		manuId = channel.ManuId();
		manuChannel = channel.ChannelId();
		adc = channel.ADC();
		fCurrentChannel++;
		return kTRUE;
	}
	else
	{
		if (NextDDL())
		{
			if (fCurrentChannel < fDecoder.GetHandler().ChannelCount())
			{
				const AliChannelInfo& channel = fDecoder.GetHandler().Channels()[fCurrentChannel];
				busPatchId = channel.BusPatchId();
				manuId = channel.ManuId();
				manuChannel = channel.ChannelId();
				adc = channel.ADC();
				fCurrentChannel++;
				return kTRUE;
			}
		}
	}
	return kFALSE;
}


UInt_t AliMUONRawStreamTrackerHP::Next(const AliChannelInfo*& channels)
{
	/// Returns the next batch of decoded channel data.
	/// [out] \param channels This is filled with a pointer to an array of
	///          channels / digits that were decoded. This method does not
	///          modify 'channels' if zero is returned.
	/// \return The number of elements in the array pointed to by 'channels'
	///    is returned. If zero is returned then there are no more channels to read.
	
	// Check if we are already at the end of the channels array. If so then we
	// need to fetch the next non empty DDL.
	if (fCurrentChannel >= fDecoder.GetHandler().ChannelCount())
	{
		do
		{
			if (not NextDDL()) return 0;
		}
		// Make sure to keep going even for empty DDL payloads:
		while (fDecoder.GetHandler().ChannelCount() == 0);
	}
	channels = fDecoder.GetHandler().Channels() + fCurrentChannel;
	return fDecoder.GetHandler().ChannelCount() - fCurrentChannel;
}

///////////////////////////////////////////////////////////////////////////////

AliMUONRawStreamTrackerHP::AliDecoderEventHandler::AliDecoderEventHandler() :
	fBusPatchId(0),
	fChannelCount(0),
	fMaxChannels(8192),
	fChannelBuffer(new AliChannelInfo[8192]),
	fRawStream(NULL),
	fBufferStart(NULL)
{
	/// Default constructor initialises the internal buffer to store decoded
	/// channel / digit data to 8192 elements.
	/// This array will grow dynamically if needed.
}


AliMUONRawStreamTrackerHP::AliDecoderEventHandler::~AliDecoderEventHandler()
{
	/// Default destructor cleans up the allocated memory.
	
	if (fChannelBuffer != NULL)
	{
		delete [] fChannelBuffer;
	}
}

void AliMUONRawStreamTrackerHP::AliDecoderEventHandler::OnNewBuffer(
		const void* buffer, UInt_t bufferSize
	)
{
	/// This is called by the high performance decoder when a new DDL payload
	/// is about to be decoded.
	/// \param buffer  The pointer to the buffer storing the DDL payload.
	/// \param bufferSize  The size of the buffer in bytes.

	// remember the start of the buffer to be used in OnError.
	fBufferStart = buffer;
	
	// Reset the number of channels found.
	fChannelCount = 0;
	
	// Check if we will have enough space in the fChannelBuffer array.
	// If we do not then we need to resize the array.
	// bufferSize / sizeof(UInt_t) will be a safe over estimate of the
	// number of channels that we will find.
	UInt_t maxPossibleChannels = bufferSize / sizeof(UInt_t);
	if (maxPossibleChannels > fMaxChannels)
	{
		if (fChannelBuffer != NULL)
		{
			delete [] fChannelBuffer;
			fChannelBuffer = NULL;
			fMaxChannels = 0;
		}
		fChannelBuffer = new AliChannelInfo[maxPossibleChannels];
		fMaxChannels = maxPossibleChannels;
	}
}


void AliMUONRawStreamTrackerHP::AliDecoderEventHandler::OnNewBusPatch(
		const AliMUONBusPatchHeaderStruct* header, const void* /*data*/
	)
{
	/// This is called by the high performance decoder when a new bus patch
	/// is found within the DDL payload. All we do is remember the bus patch ID.
	/// \param header  The bus patch header structure.
	/// \param data  The bus patch data (not used in this method).
	
	fBusPatchId = header->fBusPatchId;
}


void AliMUONRawStreamTrackerHP::AliDecoderEventHandler::OnData(UInt_t data)
{
	/// This is called by the high performance decoder when a new bus patch
	/// is found within the DDL payload. All we do is remember the bus patch ID.
	/// \param data  The bus patch data (not used in this method).
	
	assert( fChannelCount < fMaxChannels );
	
	UShort_t manuId; UChar_t channelId; UShort_t adc;
	UnpackADC(data, manuId, channelId, adc);
	fChannelBuffer[fChannelCount] = AliChannelInfo(fBusPatchId, manuId, channelId, adc);
	fChannelCount++;
}


void AliMUONRawStreamTrackerHP::AliDecoderEventHandler::OnError(
		ErrorCode error, const void* location
	)
{
	/// This is called by the high performance decoder when a error occurs
	/// when trying to decode the DDL payload. This indicates corruption in
	/// the data. This method converts the error code to a descriptive message
	/// and log this with the raw reader.
	/// \param error  The error code indicating the problem.
	/// \param location  A pointer to the location within the DDL payload buffer
	///              being decoded where the problem with the data was found.

	assert( fRawStream != NULL );
	assert( fRawStream->GetReader() != NULL );
	
	const char* msg = ErrorCodeToMessage(error);
	unsigned long pos = (unsigned long)location - (unsigned long)fBufferStart;
	
	switch (error)
	{
	case kBadPaddingWord:
	case kParityError:
		fRawStream->GetReader()->AddMinorErrorLog(error, Form("%s [At byte: %d]", msg, pos));
		break;
	default:
		fRawStream->GetReader()->AddMajorErrorLog(error, Form("%s [At byte: %d]", msg, pos));
		break;
	}
}

