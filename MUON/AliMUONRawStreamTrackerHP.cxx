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
///   const AliBusPatch* Next();
/// \endcode
///
/// This decoder tries to implement as similar an interface as possible to
/// AliMUONRawStreamTracker where possible. However certain constructs which
/// would slow us down too much are avoided.
///
/// \author Artur Szostak <artursz@iafrica.com>
//-----------------------------------------------------------------------------

#include "AliMUONRawStreamTrackerHP.h"
#include "AliMUONTrackerDDLDecoder.h"
#include "AliMUONDspHeader.h"
#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONDDLTracker.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include <cassert>
#include <iostream>
using std::cout;
using std::endl;
using std::hex;
using std::dec;

/// \cond CLASSIMP
ClassImp(AliMUONRawStreamTrackerHP)
/// \endcond

const Int_t AliMUONRawStreamTrackerHP::fgkMaxDDL = 20;

AliMUONRawStreamTrackerHP::AliMUONRawStreamTrackerHP() :
	TObject(),
fEnableErrorLogger(kFALSE),
	fDecoder(),
	fDDL(0),
	fBufferSize(8192),
	fBuffer(new UChar_t[8192]),
	fkCurrentBusPatch(NULL),
	fkCurrentData(NULL),
	fkEndOfData(NULL),
	fHadError(kFALSE),
	fDone(kFALSE),
	fDDLObject(NULL),
fTotalNumberOfGlitchErrors(0),
fTotalNumberOfParityErrors(0),
fTotalNumberOfPaddingErrors(0)
{
	///
	/// Default constructor.
	///
	
	// Must set this flag to get all information about parity errors though
	// the OnData method. OnError gets them either way.
	fDecoder.ExitOnError(false);
	fDecoder.SendDataOnParityError(true);

	fDecoder.GetHandler().SetMaxStructs(
			fDecoder.MaxBlocks(),
			fDecoder.MaxDSPs(),
			fDecoder.MaxBusPatches()
		);
}


AliMUONRawStreamTrackerHP::AliMUONRawStreamTrackerHP(AliRawReader* rawReader) :
TObject(),
fEnableErrorLogger(kFALSE),
	fDecoder(),
	fDDL(0),
	fBufferSize(8192),
	fBuffer(new UChar_t[8192]),
	fkCurrentBusPatch(NULL),
	fkCurrentData(NULL),
	fkEndOfData(NULL),
	fHadError(kFALSE),
	fDone(kFALSE),
fDDLObject(NULL),
fTotalNumberOfGlitchErrors(0),
fTotalNumberOfParityErrors(0),
fTotalNumberOfPaddingErrors(0)
{
	///
	/// Constructor with AliRawReader as argument.
	///
	
	// Must set this flag to get all information about parity errors though
	// the OnData method. OnError gets them either way.
	fDecoder.ExitOnError(false);
	fDecoder.SendDataOnParityError(true);

	fDecoder.GetHandler().SetMaxStructs(
			fDecoder.MaxBlocks(),
			fDecoder.MaxDSPs(),
			fDecoder.MaxBusPatches()
		);
  
  fDecoder.GetHandler().SetReader(rawReader);
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
	if (fDDLObject != NULL)
	{
		delete fDDLObject;
	}
}


void AliMUONRawStreamTrackerHP::First()
{
	/// Initialise or reset the iterator.
	/// The first DDL will be found and decoded.
	
	assert( GetReader() != NULL );
	
	fDDL = 0;
	fDone = kFALSE;
	NextDDL();
  fTotalNumberOfGlitchErrors = 0;
  fTotalNumberOfPaddingErrors = 0;
  fTotalNumberOfParityErrors = 0;
}


Bool_t AliMUONRawStreamTrackerHP::NextDDL()
{
	/// Reading the next tracker DDL and decode the payload with the 
	/// high performance decoder.
	/// \return kTRUE if the next DDL was successfully read and kFALSE otherwise.

	assert( GetReader() != NULL );
	
	// The temporary object if generated in GetDDLTracker, is now stale,
	// so delete it.
	if (fDDLObject != NULL)
	{
		delete fDDLObject;
		fDDLObject = NULL;
	}
	
	// Better to reset these pointers.
	fkCurrentBusPatch = NULL;
	fkCurrentData = NULL;
	fkEndOfData = NULL;
	
	while (fDDL < GetMaxDDL())
	{
		GetReader()->Reset();
		GetReader()->Select("MUONTRK", fDDL, fDDL);  // Select the DDL file to be read.
		if (GetReader()->ReadHeader()) break;
		AliDebug(3, Form("Skipping DDL %d which does not seem to be there", fDDL+1));
		fDDL++;
	}

	// If we reach the end of the DDL list for this event then reset the
	// DDL counter, mark the iteration as done and exit.
	if (fDDL >= GetMaxDDL())
	{
		fDDL = 0;
		fDone = kTRUE;
		return kFALSE;
	}
	else
	{
		fDone = kFALSE;
	}

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
	Swap(reinterpret_cast<UInt_t*>(fBuffer), dataSize / sizeof(UInt_t)); // Swap needed for mac power pc.
#endif
	
	bool result = false;
	try
	{
		// Since we might allocate memory inside OnNewBuffer in the event
		// handler we need to trap any memory allocation exception to be robust.
		result = fDecoder.Decode(fBuffer, dataSize);
		fHadError = (result == true ? kFALSE : kTRUE);
    fTotalNumberOfGlitchErrors += fDecoder.GetHandler().GlitchErrorCount();
    fTotalNumberOfParityErrors += fDecoder.GetHandler().ParityErrorCount();
		fTotalNumberOfPaddingErrors += fDecoder.GetHandler().PaddingErrorCount();
	}
	catch (const std::bad_alloc&)
	{
		AliError("Could not allocate more buffer space. Cannot decode DDL.");
		return kFALSE;
	}

	// Update the current bus patch pointers.
	fkCurrentBusPatch = fDecoder.GetHandler().FirstBusPatch();
	if (fkCurrentBusPatch != fDecoder.GetHandler().EndOfBusPatch())
	{
		fkCurrentData = fkCurrentBusPatch->GetData();
		fkEndOfData = fkCurrentData + fkCurrentBusPatch->GetDataCount();
	}
	else
	{
		// If the DDL did not have any bus patches then mark both fCurrentData
		// and fEndOfData as NULL so that in Next() we are forced to find the
		// first non empty DDL.
		fkCurrentData = fkEndOfData = NULL;
	}

	fDDL++; // Remember to increment index to next DDL.
	return kTRUE;
}


Bool_t AliMUONRawStreamTrackerHP::IsDone() const
{
	/// Indicates whether the iteration is finished or not.
	/// \return kTRUE if we already read all the digits and kFALSE if not.
	
	return fDone;
}


Bool_t AliMUONRawStreamTrackerHP::Next(Int_t& busPatchId, 
                                       UShort_t& manuId, 
                                       UChar_t& manuChannel,
                                       UShort_t& adc, 
                                       Bool_t skipParityErrors)
{
	/// Advance one step in the iteration. Returns false if finished.
	/// [out] \param busPatchId  This is filled with the bus patch ID of the digit.
	/// [out] \param manuId      This is filled with the MANU ID of the digit.
	/// [out] \param manuChannel This is filled with the MANU channel ID of the digit.
	/// [out] \param adc         This is filled with the ADC signal value of the digit.
  /// [in] \param skipParityErrors If this is kTRUE, we'll skip the buspatches that
  ///                              have some parity errors
	/// \return kTRUE if we read another digit and kFALSE if we have read all the
	///    digits already, i.e. at the end of the iteration.
	
	if (fkCurrentData == NULL) return kFALSE;
	
	fDecoder.SendDataOnParityError(skipParityErrors);

retry:
	// Check if we still have data to be returned for the current bus patch.
	if (fkCurrentData != fkEndOfData)
	{
		busPatchId = fkCurrentBusPatch->GetBusPatchId();
		AliMUONTrackerDDLDecoderEventHandler::UnpackADC(*fkCurrentData, manuId, manuChannel, adc);
		fkCurrentData++;
		return kTRUE;
	}
	else
	{
		// We hit the end of the current bus patch so check if we have any more
		// bus patches to process for the current DDL. If we do, then increment
		// the current bus patch, make sure it is not the last one and then try
		// reading the first element again.
		if (fkCurrentBusPatch != fDecoder.GetHandler().EndOfBusPatch())
		{
			fkCurrentBusPatch++;
			if (fkCurrentBusPatch != fDecoder.GetHandler().EndOfBusPatch())
			{
				fkCurrentData = fkCurrentBusPatch->GetData();
				fkEndOfData = fkCurrentData + fkCurrentBusPatch->GetDataCount();
        if ( skipParityErrors ) 
        {
          Bool_t ok(kTRUE);
          for ( Int_t i = 0; i < fkCurrentBusPatch->GetLength() && ok; ++ i )
          {
            ok = fkCurrentBusPatch->IsParityOk(i);
          }
          if (!ok) fkCurrentData = fkEndOfData;
        }
				goto retry;
			}
		}

		// This was the last bus patch in the DDL so read in the next one and
		// try reading the first data element again.
		// Note: fCurrentBusPatch is set inside NextDDL().
		if (NextDDL()) goto retry;
	}
	return kFALSE;
}


AliMUONDDLTracker* AliMUONRawStreamTrackerHP::GetDDLTracker() const
{
	/// Construct and return a pointer to the DDL payload object.
	/// \return Pointer to internally constructed AliMUONDDLTracker object.
	///         The object is owned by this class and should not be deleted
	///         by the caller.
	///
	/// \note This method should not be used just to gain access to the DDL
	/// payload, unless there is a good reason to have the AliMUONDDLTracker
	/// object. For example, if you want to modify the data and then save it
	/// to another DDL stream. Otherwise it can be an order of magnitude
	/// faster to access the DDL headers and data with the GetBlockHeader,
	/// GetDspHeader and GetBusPatch methods for example.
	/// Refer to the MUONRawStreamTracker.C macro to see how to use the fast
	/// decoder interface optimally.
	
	if (fDDLObject != NULL) return fDDLObject;
	
	fDDLObject = new AliMUONDDLTracker;
	for (Int_t iBlock = 0; iBlock < (Int_t)GetBlockCount(); iBlock++)
	{
		AliMUONBlockHeader blockHeader;
		AliMUONDspHeader dspHeader;
		AliMUONBusStruct busPatch;
		
		const AliBlockHeader* bh = GetBlockHeader(iBlock);
		// Copy block header and add it to the DDL object.
		memcpy(blockHeader.GetHeader(), bh->GetHeader(), sizeof(AliMUONBlockHeaderStruct));
		fDDLObject->AddBlkHeader(blockHeader);
		
		for (Int_t iDsp = 0; iDsp < (Int_t)bh->GetDspCount(); iDsp++)
		{
			const AliDspHeader* dh = bh->GetDspHeader(iDsp);
			// Copy DSP header and add it to the DDL object.
			memcpy(dspHeader.GetHeader(), dh->GetHeader(), sizeof(AliMUONDSPHeaderStruct));
			fDDLObject->AddDspHeader(dspHeader, iBlock);
			
			const AliBusPatch* bp = dh->GetFirstBusPatch();
			while (bp != NULL)
			{
				// Copy bus patch header, data and add everything into DDL object.
				memcpy(busPatch.GetHeader(), bp->GetHeader(), sizeof(AliMUONBusPatchHeaderStruct));
				busPatch.SetAlloc(bp->GetLength());
				memcpy(busPatch.GetData(), bp->GetData(), bp->GetDataCount()*sizeof(UInt_t));
				busPatch.SetBlockId(iBlock);
				busPatch.SetDspId(iDsp);
				fDDLObject->AddBusPatch(busPatch, iBlock, iDsp);
				bp = bp->Next();
			}
		}
	}
	
	return fDDLObject;
}


void AliMUONRawStreamTrackerHP::SetMaxBlock(Int_t blk)
{
	/// Set maximum number of blocks per DDL allowed.
	fDecoder.MaxBlocks( (UInt_t) blk );
	
	fDecoder.GetHandler().SetMaxStructs(
			fDecoder.MaxBlocks(),
			fDecoder.MaxDSPs(),
			fDecoder.MaxBusPatches()
		);
}


void AliMUONRawStreamTrackerHP::SetMaxDsp(Int_t dsp)
{
	/// Set maximum number of Dsp per block allowed.
	fDecoder.MaxDSPs( (UInt_t) dsp );
	
	fDecoder.GetHandler().SetMaxStructs(
			fDecoder.MaxBlocks(),
			fDecoder.MaxDSPs(),
			fDecoder.MaxBusPatches()
		);
}


void AliMUONRawStreamTrackerHP::SetMaxBus(Int_t bus)
{
	/// Set maximum number of Buspatch per Dsp allowed.
	fDecoder.MaxBusPatches( (UInt_t) bus );
	
	fDecoder.GetHandler().SetMaxStructs(
			fDecoder.MaxBlocks(),
			fDecoder.MaxDSPs(),
			fDecoder.MaxBusPatches()
		);
}

AliRawReader* AliMUONRawStreamTrackerHP::GetReader()
{
  return fDecoder.GetHandler().GetReader();
}

void AliMUONRawStreamTrackerHP::SetReader(AliRawReader* reader)
{
  fDecoder.GetHandler().SetReader(reader);
}

///////////////////////////////////////////////////////////////////////////////

void AliMUONRawStreamTrackerHP::AliBlockHeader::Print() const
{
	/// Print header to screen.
	
	cout << "CRT info"        << endl;
	if (fHeader == NULL)
	{
		cout << "Header is NULL" << endl;
		return;
	}
	cout << "DataKey: 0x"     << hex << fHeader->fDataKey << dec << endl;
	cout << "TotalLength: "   << fHeader->fTotalLength << endl;
	cout << "Length: "        << fHeader->fLength << endl;
	cout << "DspId: "         << fHeader->fDSPId << endl;
	cout << "L0Trigger: "     << fHeader->fL0Trigger << endl;
	cout << "MiniEventId: "   << fHeader->fMiniEventId<< endl; 
	cout << "EventId1: "      << fHeader->fEventId1 << endl;
	cout << "EventId2: "      << fHeader->fEventId2 << endl;
}


void AliMUONRawStreamTrackerHP::AliDspHeader::Print() const
{
	/// Print header to screen.
	
	cout << "FRT info"        << endl;
	if (fHeader == NULL)
	{
		cout << "Header is NULL" << endl;
		return;
	}
	cout << "DataKey: 0x"     << hex << fHeader->fDataKey << dec << endl;
	cout << "TotalLength: "   << fHeader->fTotalLength << endl;
	cout << "Length : "       << fHeader->fLength << endl;
	cout << "DspId: "         << fHeader->fDSPId << endl;
	cout << "BlkL1ATrigger: " << fHeader->fBlkL1ATrigger << endl;
	cout << "MiniEventId: "   << fHeader->fMiniEventId << endl;
	cout << "L1ATrigger: "    << fHeader->fL1ATrigger << endl;
	cout << "L1RTrigger: "    << fHeader->fL1RTrigger << endl;
	cout << "PaddingWord: "   << fHeader->fPaddingWord << endl;
	cout << "ErrorWord: "     << fHeader->fErrorWord << endl;
}


void AliMUONRawStreamTrackerHP::AliBusPatch::Print(const Option_t* opt) const
{
	/// Print header to screen.
	cout << "Bus patch info" << endl;
	if (fHeader == NULL)
	{
		cout << "Header is NULL" << endl;
		return;
	}
	cout << "DataKey: 0x"    << hex << fHeader->fDataKey << dec << endl;
	cout << "fTotalLength: " << fHeader->fTotalLength << endl;
	cout << "fLength: "      << fHeader->fLength << endl;
	cout << "fBusPatchId: "  << fHeader->fBusPatchId << endl;

	if (TString(opt).Contains("all"))
	{
		for (UInt_t i = 0; i < fHeader->fLength; ++i)
			cout << "Data["<< i << "] = " << fData[i] << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

AliMUONRawStreamTrackerHP::AliDecoderEventHandler::AliDecoderEventHandler() :
fRawReader(0x0),
	fBufferStart(NULL),
	fBlockCount(0),
	fBlocks(NULL),
	fDSPs(NULL),
	fBusPatches(NULL),
	fEndOfBusPatches(NULL),
	fMaxChannels(8192),
	fParityOk(new Bool_t[8192]),
	fCurrentBlock(NULL),
	fCurrentDSP(NULL),
	fCurrentBusPatch(NULL),
	fCurrentParityOkFlag(NULL),
	fParityErrors(0),
	fGlitchErrors(0),
	fPaddingErrors(0),
	fWarnings(kTRUE),
  fMaxBlocks(),
  fMaxDsps(),
  fMaxBusPatches()
{
	/// Default constructor initialises the internal parity flags buffer to
	/// store 8192 elements. This array will grow dynamically if needed.
}


AliMUONRawStreamTrackerHP::AliDecoderEventHandler::~AliDecoderEventHandler()
{
	/// Default destructor cleans up the allocated memory.
	
	if (fParityOk != NULL) delete [] fParityOk;
	if (fBlocks != NULL) delete [] fBlocks;
	if (fDSPs != NULL) delete [] fDSPs;
	if (fBusPatches != NULL) delete [] fBusPatches;
}


void AliMUONRawStreamTrackerHP::AliDecoderEventHandler::SetMaxStructs(
		UInt_t maxBlocks, UInt_t maxDsps, UInt_t maxBusPatches
	)
{
	/// Sets the maximum number of structures allowed.
	
	// Start by clearing the current arrays.
	if (fBlocks != NULL)
	{
		delete [] fBlocks;
		fBlocks = NULL;
	}
	if (fDSPs != NULL)
	{
		delete [] fDSPs;
		fDSPs = NULL;
	}
	if (fBusPatches != NULL)
	{
		delete [] fBusPatches;
		fBusPatches = NULL;
	}
	fCurrentBlock = NULL;
	fCurrentDSP = NULL;
	fCurrentBusPatch = NULL;
	
	// Allocate new memory.
	fBlocks = new AliBlockHeader[maxBlocks];
	fDSPs = new AliDspHeader[maxBlocks*maxDsps];
	fBusPatches = new AliBusPatch[maxBlocks*maxDsps*maxBusPatches];
	fEndOfBusPatches = fBusPatches;

  fMaxBlocks = maxBlocks;
  fMaxDsps = maxDsps;
  fMaxBusPatches = maxBusPatches;
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

	// Reset error counters.
	fParityErrors = 0;
	fGlitchErrors = 0;
	fPaddingErrors = 0;

	// Check if we will have enough space in the fParityOk array.
	// If we do not then we need to resize the array.
	// bufferSize / sizeof(UInt_t) will be a safe over estimate of the
	// number of channels that we will find.
	UInt_t maxChannelsPossible = bufferSize / sizeof(UInt_t);
	if (maxChannelsPossible > fMaxChannels)
	{
		if (fParityOk != NULL)
		{
			delete [] fParityOk;
			fParityOk = NULL;
			fMaxChannels = 0;
		}
		fParityOk = new Bool_t[maxChannelsPossible];
		fMaxChannels = maxChannelsPossible;
	}
	
	// Reset the current pointers which will be used to track where we need to
	// fill fBlocks, fDSPs, fBusPatches and the parity flag. We have to subtract
	// one space because we will increment the pointer the first time in the
	// OnNewXZY methods.
	fCurrentBlock = fBlocks-1;
	fCurrentDSP = fDSPs-1;
	fCurrentBusPatch = fBusPatches-1;
	fCurrentParityOkFlag = fParityOk-1;
	fBlockCount = 0;
}


void AliMUONRawStreamTrackerHP::AliDecoderEventHandler::OnError(
		ErrorCode error, const void* location
	)
{
	/// This is called by the high performance decoder when a error occurs
	/// when trying to decode the DDL payload. This indicates corruption in
	/// the data. This method converts the error code to a descriptive message
	/// and logs this with the raw reader.
	/// \param error  The error code indicating the problem.
	/// \param location  A pointer to the location within the DDL payload buffer
	///              being decoded where the problem with the data was found.

	assert( fRawReader != NULL );
	
	Char_t* message = NULL;
	UInt_t word = 0;

	switch (error)
	{
	case kGlitchFound:
		fGlitchErrors++;
		message = Form(
			"Glitch error detected in DSP %d, skipping event ",
			fCurrentBlock->GetDspId()
		);
		fRawReader->AddMajorErrorLog(error, message);
		break;

	case kBadPaddingWord:
		fPaddingErrors++;
		// We subtract 1 from the current numbers of blocks, DSPs
		// and bus patches to get the indices.
		message = Form(
			"Padding word error for iBlock %d, iDsp %d, iBus %d\n", 
			fBlockCount-1,
			fCurrentBlock->GetDspCount()-1,
			fCurrentDSP->GetBusPatchCount()-1
		);
		fRawReader->AddMinorErrorLog(error, message);
		break;

	case kParityError:
		fParityErrors++;
		// location points to the incorrect data word and
		// fCurrentBusPatch->GetData() returns a pointer to the start of
		// bus patches data, so the difference divided by 4 gives the 32
		// bit word number.
		word = ((unsigned long)location - (unsigned long)fCurrentBusPatch->GetData())
				/ sizeof(UInt_t);
		message = Form(
			"Parity error in word %d for manuId %d and channel %d in buspatch %d\n", 
			word,
			fCurrentBusPatch->GetManuId(word),
			fCurrentBusPatch->GetChannelId(word),
			fCurrentBusPatch->GetBusPatchId()
		);
		fRawReader->AddMinorErrorLog(error, message);
		break;

	default:
		message = Form(
			"%s (At byte %d in DDL.)",
			ErrorCodeToMessage(error),
			(unsigned long)location - (unsigned long)fBufferStart + sizeof(AliRawDataHeader)
		);
		fRawReader->AddMajorErrorLog(error, message);
		break;
	}

	if (fWarnings)
	{
		AliWarningGeneral(
				"AliMUONRawStreamTrackerHP::AliDecoderEventHandler",
				message
			);
	}
}

