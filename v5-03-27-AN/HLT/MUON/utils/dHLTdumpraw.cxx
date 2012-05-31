/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author:                                                                *
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
 * @file   dHLTdumpraw.cxx
 * @author Artur Szostak <artursz@iafrica.com>,
 *         Seforo Mohlalisi <seforomohlalisi@yahoo.co.uk>
 * @date   1 July 2007
 * @brief  Command line utility to dump dHLT's internal raw data blocks.
 */

// We define NDEBUG for the AliHLTMUONDataBlockReader.h header file since this
// program by definition handles corrupt data. So we do not need the assertions
// in the AliHLTMUONDataBlockReader class to be checked.
#ifndef NDEBUG
#define NDEBUG
#endif
#include "AliHLTMUONDataBlockReader.h"
#if defined(DEBUG) && defined(NDEBUG)
#undef NDEBUG
#endif

#include "AliHLTMUONUtils.h"
#include "AliHLTMUONConstants.h"
#include "Rtypes.h"
#include "AliRawDataHeader.h"
#include "AliMUONTrackerDDLDecoder.h"
#include "AliMUONTriggerDDLDecoder.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include "AliLog.h"
#include "TClassTable.h"
#include "TString.h"
#include "TRegexp.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <new>
#include <fstream>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::showbase;
using std::noshowbase;
using std::hex;
using std::dec;

#include <iomanip>
using std::setw;
using std::left;
using std::right;
using std::internal;


#define CMDLINE_ERROR 1
#define PARSE_ERROR 2
#define SYSTEM_ERROR 3
#define FATAL_ERROR 4
#define HLTSYSTEM_ERROR 5


// Adding enum types for extending AliHLTMUONDataBlockType with the
// raw DDL data types.
enum AliHLTMUONDDLRawDataType
{
	kTrackerDDLRawData = 10,
	kTriggerDDLRawData = 11
};


void PrintRubbishData(AliHLTUInt32_t offset, const char* padByte, AliHLTUInt32_t padCount)
{
	if (padCount == 0) return;
	
	cerr << "ERROR: Found the following unexpected rubbish data at the"
		" end of the data block:" << endl;
	cerr << "Byte #\tValue\tCharacter" << endl;
	for (AliHLTUInt32_t i = 0; i < padCount; i++)
	{
		short value = short(padByte[i]) & 0xFF;
		char character = padByte[i];
		cerr << offset + i + 1 << "\t"
			<< noshowbase << hex << "0x" << value << dec << "\t"
			<< character << endl;
	}
}


void PrintBitPattern(AliHLTUInt32_t value, int bitcount = 32)
{
	// Prints a bit pattern to cout.
	
	for (int i = bitcount-1; i >= 0; i--)
	{
		if ( ((value >> i) & 0x1) == 1 )
			cout << "1";
		else
			cout << "0";
	}
}


template <typename FieldType>
int CheckHeaderField(
		FieldType& field, const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	const char* fieldptr = reinterpret_cast<const char*>(&field);
	const char* endptr = buffer + bufferSize;
	AliHLTUInt32_t bufferRemaining = endptr > fieldptr ? endptr - fieldptr : 0;
	
	if (bufferRemaining < sizeof(field))
	{
		cout << "..." << endl; // We may be half way through printing a line so end it.
		cerr << "ERROR: The data block is too short. The header is corrupt." << endl;
		if (continueParse)
		{
			AliHLTUInt32_t offset = fieldptr - buffer;
			PrintRubbishData(offset, fieldptr, bufferRemaining);
		}
		return PARSE_ERROR;
	}
	return EXIT_SUCCESS;
}


template <typename FieldType>
int CheckField(
		FieldType& field, const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	const char* fieldptr = reinterpret_cast<const char*>(&field);
	const char* endptr = buffer + bufferSize;
	AliHLTUInt32_t bufferRemaining = endptr > fieldptr ? endptr - fieldptr : 0;
	
	if (bufferRemaining < sizeof(field))
	{
		cout << "..." << endl; // We may be half way through printing a line so end it.
		cerr << "ERROR: The data block is too short. The data is corrupt." << endl;
		if (continueParse)
		{
			AliHLTUInt32_t offset = fieldptr - buffer;
			PrintRubbishData(offset, fieldptr, bufferRemaining);
		}
		return PARSE_ERROR;
	}
	return EXIT_SUCCESS;
}


template <typename BlockType>
int CheckCommonHeader(
		BlockType& block, const char* /*buffer*/, unsigned long bufferSize,
		bool continueParse
	)
{
	int result = EXIT_SUCCESS;

	// Check the fRecordWidth field in the common header.
	if (block.CommonBlockHeader().fRecordWidth !=
		sizeof(typename BlockType::ElementType))
	{
		cerr << "ERROR: The record width found in the header is incorrect."
			" Found a record width of "
			<< block.CommonBlockHeader().fRecordWidth << " bytes, but expected"
			" a value of " << sizeof(typename BlockType::ElementType)
			<< " bytes." << endl;
		result = PARSE_ERROR;
		if (not continueParse) return result;
	}
	
	if (not block.BufferSizeOk())
	{
		cerr << "ERROR: The size of the file is incorrect. It is "
			<< bufferSize << " bytes big, but according"
			" to the data block header it should be " << block.BytesUsed()
			<< " bytes." << endl;
		result = PARSE_ERROR;
		if (not continueParse) return result;
	}
	
	return result;
}


template <typename BlockType>
AliHLTUInt32_t CalculateNEntries(BlockType& block, unsigned long bufferSize)
{
	// Calculate how many entries we can display. If the buffer size is correct
	// we just use the number of entries the block specifies. Otherwise we need
	// to calculate it from the buffer size.
	AliHLTUInt32_t nentries;
	if (block.BytesUsed() == bufferSize)
	{
		nentries = block.Nentries();
	}
	else
	{
		AliHLTInt32_t dataSize = bufferSize
			- sizeof(typename BlockType::HeaderType);
		nentries = dataSize / sizeof(typename BlockType::ElementType);
		if (dataSize % sizeof(typename BlockType::ElementType) > 0)
			nentries++;
	}
	return nentries;
}


namespace
{
	/**
	 * Common methods for DDL decoder event handlers.
	 */
	class AliDecoderHandler
	{
	public:
		AliDecoderHandler() :
			fBufferStart(NULL),
			fDumpStart(NULL),
			fDumpData(false)
		{
		}
		
		virtual ~AliDecoderHandler() {}
		
	protected:
		
		// Do not allow copying of this class.
		AliDecoderHandler(const AliDecoderHandler& obj);
		AliDecoderHandler& operator = (const AliDecoderHandler& obj);
		
		void HandleError(
				const char* errorMessage, int errorCode,
				const char* errorCodeString, const void* location
			);
	
		void TryDumpCorruptData(const void* dumpEnd);
		
		const void* fBufferStart;  ///< Start location of buffer.
		const void* fDumpStart;  ///< Start location of corrupt data to dump.
		bool fDumpData;  ///< Flag indicating if fDumpStart points to corrupt data and should be dumped.
	};
	
	
	void AliDecoderHandler::HandleError(
			const char* errorMessage, int errorCode,
			const char* errorCodeString, const void* location
		)
	{
		unsigned long offset = (unsigned long)location - (unsigned long)fBufferStart
			+ sizeof(AliRawDataHeader);
		
		cerr << "ERROR: " << errorMessage
			<< " [Error code = " << errorCode << " ("
			<< errorCodeString << "), at byte "
			<< offset << " (" << noshowbase << hex << "0x"
			<< offset << dec << ")]" << endl;
		
		if (fDumpStart == NULL) fDumpStart = location;
		fDumpData = true;
	}
	
	
	void AliDecoderHandler::TryDumpCorruptData(const void* dumpEnd)
	{
		if (dumpEnd < fDumpStart) return;
		if (not fDumpData) return;
		
		unsigned long startOffset = (unsigned long)fDumpStart - (unsigned long)fBufferStart
			+ sizeof(AliRawDataHeader);
		unsigned long endOffset = (unsigned long)dumpEnd - (unsigned long)fBufferStart
			+ sizeof(AliRawDataHeader);
		if (endOffset - startOffset > 264)
		{
			endOffset = startOffset + 264;
			dumpEnd = reinterpret_cast<const char*>(fBufferStart) + endOffset;
		}
		cerr << "Dumping corrupt data words from byte " << startOffset
			<< " (" << noshowbase << hex << "0x" << startOffset
			<< dec << "), to byte " << endOffset << " (" << noshowbase
			<< hex << "0x" << endOffset << dec << "):" << endl;
		const UInt_t* start = reinterpret_cast<const UInt_t*>(fDumpStart);
		const UInt_t* end = reinterpret_cast<const UInt_t*>(dumpEnd);
		cerr << "     Start byte     | Data words" << endl;
		for (const UInt_t* current = start; current < end; current++)
		{
			unsigned long currentByte = (unsigned long)current
				- (unsigned long)fBufferStart + sizeof(AliRawDataHeader);
			cerr << right << setw(9) << dec << currentByte << setw(0)
				<< " 0x" << left << setw(7) << noshowbase << hex
				<< currentByte << setw(0) << right << " | ";
			char fillChar = cerr.fill();
			cerr.fill('0');
			for (int i = 0; i < 4 and current < end; i++, current++)
			{
				cerr << noshowbase << hex << "0x" << setw(8)
					<< (*current) << setw(0) << dec << " ";
			}
			cerr.fill(fillChar);
			cerr << endl;
		}
		fDumpStart = NULL;
		fDumpData = false;
	}

	/**
	 * Event handler for the tracker DDL decoder.
	 * It simply prints the structure to standard output.
	 */
	class AliTrackerDecoderHandler :
		public AliMUONTrackerDDLDecoderEventHandler, public AliDecoderHandler
	{
	public:
		AliTrackerDecoderHandler() :
			AliMUONTrackerDDLDecoderEventHandler(),
			AliDecoderHandler(),
			fBlockNum(0),
			fDSPNum(0),
			fBusPatchNum(0)
		{}
		
		virtual ~AliTrackerDecoderHandler() {}
	
		void OnNewBuffer(const void* buffer, UInt_t /*bufferSize*/)
		{
			fBufferStart = buffer;
			fBlockNum = fDSPNum = fBusPatchNum = 0;
		}
		
		void OnEndOfBuffer(const void* buffer, UInt_t bufferSize)
		{
			const char* bufferEnd =
				reinterpret_cast<const char*>(buffer) + bufferSize;
			TryDumpCorruptData(bufferEnd);
		}
		
		void OnNewBlock(const AliMUONBlockHeaderStruct* header, const void* /*data*/);
		
		void OnNewDSP(const AliMUONDSPHeaderStruct* header, const void* /*data*/);
		
		void OnNewBusPatch(const AliMUONBusPatchHeaderStruct* header, const void* /*data*/);
		
		void OnData(UInt_t data, bool parityError);
		
		void OnError(ErrorCode error, const void* location);
		
	private:
	
		UInt_t fBlockNum;     ///< Number of block being processed [1..maxBlock].
		UInt_t fDSPNum;       ///< Number of DSP being processed [1..maxDSP].
		UInt_t fBusPatchNum;  ///< Number of bus patch being processed [1..maxBusPatch].
	};
	
	
	void AliTrackerDecoderHandler::OnNewBlock(const AliMUONBlockHeaderStruct* header, const void* /*data*/)
	{
		TryDumpCorruptData(header);
		fBlockNum++;
		char fillChar = cout.fill();  // store current fill char to restore it.
		cout << "================================ Block header "
			<< setw(3) << fBlockNum << setw(0)
			<< " ================================" << endl;
		cout << "   Data key word for block : 0x" << noshowbase << setw(8) << setfill('0')
			<< hex << header->fDataKey << dec << setfill(fillChar) << setw(0) << endl;
		cout << "     Total Length of block : " << header->fTotalLength << endl;
		cout << "        Length of raw data : " << header->fLength << endl;
		cout << "                    DSP ID : " << header->fDSPId << endl;
		cout << "           L0 trigger word : " << header->fL0Trigger << " (0x"
			<< noshowbase << setw(8) << setfill('0') << hex << header->fL0Trigger
			<< dec << setw(0) << setfill(fillChar) << ")" << endl;
		cout << "             Mini event ID : " << header->fMiniEventId << " (0x"
			<< noshowbase << setw(8) << setfill('0') << hex << header->fMiniEventId
			<< dec << setw(0) << setfill(fillChar) << ")" << endl;
		cout << "Event ID In Bunch Crossing : " << header->fEventId1 << " (0x"
			<< noshowbase << setw(8) << setfill('0') << hex << header->fEventId1
			<< dec << setw(0) << setfill(fillChar) << ")" << endl;
		cout << "  Event ID In Orbit Number : " << header->fEventId2 << " (0x"
			<< noshowbase << setw(8) << setfill('0') << hex << header->fEventId2
			<< dec << setw(0) << setfill(fillChar) << ")" << endl;
	}
	
	
	void AliTrackerDecoderHandler::OnNewDSP(const AliMUONDSPHeaderStruct* header, const void* /*data*/)
	{
		TryDumpCorruptData(header);
		fDSPNum++;
		char fillChar = cout.fill();  // store current fill char to restore it.
		cout << "================================= DSP header "
			<< setw(3) << fDSPNum << setw(0)
			<< " =================================" << endl;
		cout << "                     Data key word for FRT : 0x" << noshowbase
			<< setw(8) << setfill('0') << hex << header->fDataKey << dec
			<< setfill(fillChar) << setw(0) << endl;
		cout << "                 Total length of structure : " << header->fTotalLength << endl;
		cout << "                        Length of raw data : " << header->fLength << endl;
		cout << "                                    DSP ID : " << header->fDSPId << endl;
		cout << "        L1 accept in block Structure (CRT) : " << header->fBlkL1ATrigger
			<< " (0x" << noshowbase << setw(8) << setfill('0') << hex
			<< header->fBlkL1ATrigger << dec << setw(0) << setfill(fillChar) << ")" << endl;
		cout << "           Mini event ID in bunch crossing : " << header->fMiniEventId
			<< " (0x" << noshowbase << setw(8) << setfill('0') << hex
			<< header->fMiniEventId << dec << setw(0) << setfill(fillChar) << ")" << endl;
		cout << "Number of L1 accept in DSP Structure (FRT) : " << header->fL1ATrigger << endl;
		cout << "Number of L1 reject in DSP Structure (FRT) : " << header->fL1RTrigger << endl;
		const char* paddingWordValue = ") (INVALID)";
		if (header->fPaddingWord == 0) paddingWordValue = ") (false)";
		if (header->fPaddingWord == 1) paddingWordValue = ") (true)";
		cout << "     Padding word set for 64 bits transfer : " << header->fPaddingWord
			<< " (0x" << noshowbase << setw(8) << setfill('0') << hex
			<< header->fPaddingWord << setw(0) << setfill(fillChar) << dec
			<< paddingWordValue << endl;
		cout << "                                Error word : " << header->fErrorWord
			<< " (" << noshowbase << setw(8) << setfill('0') << hex
			<< header->fErrorWord << setw(0) << setfill(fillChar)
			<< dec << ")" << endl;
	}
	
	
	void AliTrackerDecoderHandler::OnNewBusPatch(const AliMUONBusPatchHeaderStruct* header, const void* /*data*/)
	{
		TryDumpCorruptData(header);
		fBusPatchNum++;
		char fillChar = cout.fill();  // store current fill char to restore it.
		cout << "============================== Bus patch header "
			<< setw(3) << fBusPatchNum << setw(0)
			<< " ==============================" << endl;
		cout << "Data key word for bus patch : 0x" << noshowbase << setw(8)
			<< setfill('0') << hex << header->fDataKey << dec
			<< setfill(fillChar) << setw(0) << endl;
		cout << "  Total length of structure : " << header->fTotalLength << endl;
		cout << "         Length of raw data : " << header->fLength << endl;
		cout << "               Bus patch ID : " << header->fBusPatchId << endl;
		if (header->fLength > 0)
		{
			cout << "    Raw bits |      Manu ID | Manu channel |   ADC Signal" << endl;
			cout << "----------------------------------------------------------" << endl;
		}
	}
	
	
	void AliTrackerDecoderHandler::OnData(UInt_t data, bool parityError)
	{
		UShort_t manuId, adc;
		UChar_t manuChannel;
		UnpackADC(data, manuId, manuChannel, adc);
		char fillChar = cout.fill();  // store current fill char to restore it.
		cout << noshowbase << "  0x" << setfill('0') << hex << setw(8) << data
			<< setw(0) << dec << setfill(fillChar) << "   "
			<< setw(12) << manuId << setw(0) << "   "
			<< setw(12) << (UInt_t)manuChannel << setw(0) << "   "
			<< setw(12) << adc << setw(0);
		if (parityError)
		{
			cout << " <= WARNING: Raw data word with parity error!" << endl;
		}
		else
		{
			cout << endl;
		}
	}
	
	
	void AliTrackerDecoderHandler::OnError(ErrorCode error, const void* location)
	{
		TryDumpCorruptData(location);
		HandleError(
			ErrorCodeToMessage(error), error,
			ErrorCodeToString(error), location
		);
	}

	/**
	 * Event handler for the trigger DDL decoder.
	 * It simply prints the structure to standard output.
	 */
	class AliTriggerDecoderHandler :
		public AliMUONTriggerDDLDecoderEventHandler, public AliDecoderHandler
	{
	public:
		AliTriggerDecoderHandler() :
			AliMUONTriggerDDLDecoderEventHandler(),
			AliDecoderHandler(),
			fRegNum(0),
			fLocNum(0)
		{}
		
		virtual ~AliTriggerDecoderHandler() {}
	
		void OnNewBuffer(const void* buffer, UInt_t /*bufferSize*/)
		{
			fBufferStart = buffer;
			fRegNum = fLocNum = 0;
		}
		
		void OnEndOfBuffer(const void* buffer, UInt_t bufferSize)
		{
			const char* bufferEnd =
				reinterpret_cast<const char*>(buffer) + bufferSize;
			TryDumpCorruptData(bufferEnd);
		}
		
		const char* EventTypeToString(UInt_t type);
		
		const char* DarcTypeToString(UInt_t type);
		
		void OnDarcHeader(
				UInt_t header,
				const AliMUONDarcScalarsStruct* scalars,
				const void* data
			);
		
		void OnGlobalHeader(
				const AliMUONGlobalHeaderStruct* header,
				const AliMUONGlobalScalarsStruct* scalars,
				const void* /*data*/
			);
		
		void OnNewRegionalStruct(
				const AliMUONRegionalHeaderStruct* regionalStruct,
				const AliMUONRegionalScalarsStruct* scalars,
				const void* /*data*/
			);
		
		void OnLocalStruct(
				const AliMUONLocalInfoStruct* localStruct,
				const AliMUONLocalScalarsStruct* scalars
			);
		
		void OnError(ErrorCode error, const void* location);
	
	private:
	
		UInt_t fRegNum;  ///< Number of block being processed [1..maxReg].
		UInt_t fLocNum;  ///< Number of DSP being processed [1..maxLoc].
	};
	
	
	const char* AliTriggerDecoderHandler::EventTypeToString(UInt_t type)
	{
		switch (type)
		{
		case 0: return "Other software trigger";
		case 1: return "Physics";
		case 2: return "Start of run";
		case 3: return "End of run";
		default: return "UNKNOWN";
		}
	}
	
	
	const char* AliTriggerDecoderHandler::DarcTypeToString(UInt_t type)
	{
		typedef AliMUONTriggerDDLDecoder<AliMUONTriggerDDLDecoderEventHandler> AliDec;
		if (type == AliDec::DarcDefaultType()) return "Default";
		if (type == AliDec::DarcVadorhType()) return "Vadorh";
		return "UNKNOWN";
	}
	
	
	void AliTriggerDecoderHandler::OnDarcHeader(
			UInt_t header,
			const AliMUONDarcScalarsStruct* scalars,
			const void* data
		)
	{
		if (scalars != NULL)
			TryDumpCorruptData(scalars);
		else
			TryDumpCorruptData(data);
		
		cout << "================================ DARC header =====================================" << endl;
		char fillChar = cout.fill();  // store current fill char to restore it.
		cout << noshowbase << "         Header bit pattern : 0x" << setfill('0') << hex << setw(8)
			<< header << setw(0) << dec << setfill(fillChar) << endl;
		UInt_t eventType = GetDarcEventType(header);
		cout << "                 Event type : " << eventType << " ("
			<< EventTypeToString(eventType) << ")" << endl;
		cout << "           Application type : " << ((header >> 27) & 0x7) << endl;
		UInt_t darcType = GetDarcType(header);
		cout << "                  DARC type : " << darcType << " ("
			<< DarcTypeToString(darcType) << ")" << endl;
		cout << "              Serial number : " << (UInt_t)GetDarcSerialNb(header) << endl;
		cout << "                    Version : " << (UInt_t)GetDarcVersion(header) << endl;
		cout << "           VME trigger used : " << boolalpha << GetDarcVMETrig(header) << endl;
		cout << "Global card data occurrence : " << boolalpha << GetDarcGlobalFlag(header) << endl;
		cout << "      CTP or LTU interfaced : " << boolalpha << GetDarcCTPTrig(header) << endl;
		cout << "             DAQ interfaced : " << boolalpha << GetDarcDAQFlag(header) << endl;
		cout << "  Regional cards occurrence : 0x" << noshowbase << hex << setw(2) << setfill('0')
			<< (UInt_t)GetDarcRegPattern(header) << dec << setfill(fillChar) << setw(0)
			<< " [bits: ";
		PrintBitPattern(((UInt_t)GetDarcRegPattern(header) >> 4) & 0xF, 4); cout << " ";
		PrintBitPattern(((UInt_t)GetDarcRegPattern(header) >> 0) & 0xF, 4);
		cout << "]" << endl;
		
		cout << "================================ DARC scalars ====================================" << endl;
		if (scalars != NULL)
		{
			cout << "    Trigger | Received |     Used" << endl;
			cout << "----------------------------------" << endl;
			cout << "         L0   " << setw(8) << ((scalars->fL0R >> 16) & 0xFFFF) << setw(0)
				<< "   " << setw(8) << ((scalars->fL0R >> 0) & 0xFFFF) << setw(0) << endl;
			cout << " L1 physics   " << setw(8) << ((scalars->fL1P >> 16) & 0xFFFF) << setw(0)
				<< "   " << setw(8) << ((scalars->fL1P >> 0) & 0xFFFF) << setw(0) << endl;
			cout << "L1 software   " << setw(8) << ((scalars->fL1S >> 16) & 0xFFFF) << setw(0)
				<< "   " << setw(8) << ((scalars->fL1S >> 0) & 0xFFFF) << setw(0) << endl;
			cout << "  L2 accept   " << setw(8) << ((scalars->fL2A >> 16) & 0xFFFF) << setw(0)
				<< "   " << setw(8) << ((scalars->fL2A >> 0) & 0xFFFF) << setw(0) << endl;
			cout << "  L2 reject   " << setw(8) << ((scalars->fL2R >> 16) & 0xFFFF) << setw(0)
				<< "   " << setw(8) << ((scalars->fL2R >> 0) & 0xFFFF) << setw(0) << endl;
			cout << "           Clock : " << scalars->fClk << endl;
			cout << "Hold (dead time) : " << scalars->fHold << endl;
			cout << "      Spare word : " << scalars->fSpare << " (0x"
				<< setw(8) << setfill('0') << hex << scalars->fSpare <<
				setw(0) << setfill(fillChar) << dec << ")" << endl;
		}
		else
		{
			cout << "(None found)" << endl;
		}
	}
	
	
	void AliTriggerDecoderHandler::OnGlobalHeader(
			const AliMUONGlobalHeaderStruct* header,
			const AliMUONGlobalScalarsStruct* scalars,
			const void* /*data*/
		)
	{
		TryDumpCorruptData(header);
		
		cout << "=============================== Global header ====================================" << endl;
		cout << "Global input from regional controllers:" << endl;
		cout << "      |                 |         High pT           |          Low pT          " << endl;
		cout << "      |    Bit pattern  | Single mu |    mu pair    | Single mu |    mu pair   " << endl;
		cout << "Input |  Hex |   Binary | -ve | +ve | unlike | like | -ve | +ve | unlike | like" << endl;
		cout << "--------------------------------------------------------------------------------" << endl;
		char fillChar = cout.fill();  // store current fill char to restore it.
		for (int i = 0; i < 4; i++)
		{
			union
			{
				UInt_t fWord;
				UChar_t fByte[4];
			} input;
			input.fWord = header->fInput[i];
			for (int j = 0; j < 4; j++)
			{
				cout << setw(5) << (i*4+j) << "   0x" << noshowbase << setw(2)
					<< setfill('0') << hex << (UInt_t)input.fByte[j] << setw(0)
					<< setfill(fillChar) << dec << "   ";
				PrintBitPattern(input.fByte[j], 8);
				cout << ((((input.fByte[j] >> 7) & 0x1) == 1) ? "   yes" : "    no");
				cout << ((((input.fByte[j] >> 6) & 0x1) == 1) ? "   yes" : "    no");
				cout << ((((input.fByte[j] >> 5) & 0x1) == 1) ? "     yes " : "      no ");
				cout << ((((input.fByte[j] >> 4) & 0x1) == 1) ? "    yes" : "     no");
				cout << ((((input.fByte[j] >> 3) & 0x1) == 1) ? "   yes" : "    no");
				cout << ((((input.fByte[j] >> 2) & 0x1) == 1) ? "   yes" : "    no");
				cout << ((((input.fByte[j] >> 1) & 0x1) == 1) ? "     yes " : "      no ");
				cout << ((((input.fByte[j] >> 0) & 0x1) == 1) ? "    yes" : "     no");
				cout << endl;
			}
		}
		cout << "        Global ouput : 0x" << noshowbase << setw(2) << setfill('0') << hex
			<< (UInt_t)GetGlobalOutput(header) << setw(0) << setfill(fillChar) << dec << " [Bits: ";
		PrintBitPattern(((UInt_t)GetGlobalOutput(header) >> 4) & 0xF, 4); cout << " ";
		PrintBitPattern(((UInt_t)GetGlobalOutput(header) >> 0) & 0xF, 4);
		cout << "]" << endl;
		cout << "          [ unlike sign |  like sign  | single muon ]" << endl;
		cout << "          [ high |  low | high |  low | high |  low ]" << endl;
		cout << "          [-----------------------------------------]" << endl;
		cout << "          [ ";
		cout << ((((GetGlobalOutput(header) >> 5) & 0x1) == 1) ? " yes" : "  no");
		cout << ((((GetGlobalOutput(header) >> 4) & 0x1) == 1) ? "    yes" : "     no");
		cout << ((((GetGlobalOutput(header) >> 3) & 0x1) == 1) ? "    yes" : "     no");
		cout << ((((GetGlobalOutput(header) >> 2) & 0x1) == 1) ? "    yes" : "     no");
		cout << ((((GetGlobalOutput(header) >> 1) & 0x1) == 1) ? "    yes" : "     no");
		cout << ((((GetGlobalOutput(header) >> 0) & 0x1) == 1) ? "    yes" : "     no");
		cout << " ]" << endl;
		cout << "Global configuration : 0x" << noshowbase << setw(4) << setfill('0') << hex
			<< GetGlobalConfig(header) << setw(0) << setfill(fillChar) << dec << " [Bits: ";
		PrintBitPattern(((UInt_t)GetGlobalConfig(header) >> 12) & 0xF, 4); cout << " ";
		PrintBitPattern(((UInt_t)GetGlobalConfig(header) >> 8) & 0xF, 4); cout << " ";
		PrintBitPattern(((UInt_t)GetGlobalConfig(header) >> 4) & 0xF, 4); cout << " ";
		PrintBitPattern(((UInt_t)GetGlobalConfig(header) >> 0) & 0xF, 4);
		cout << "]" << endl;
		
		cout << "=============================== Global scalars ===================================" << endl;
		if (scalars != NULL)
		{
			cout << "           Number of L0 triggers : " << scalars->fL0 << endl;
			cout << "          Number of clock cycles : " << scalars->fClk << endl;
			cout << " Number of unlike mu pair low pT : " << scalars->fScaler[0] << endl;
			cout << "Number of unlike mu pair high pT : " << scalars->fScaler[1] << endl;
			cout << "   Number of like mu pair low pT : " << scalars->fScaler[2] << endl;
			cout << "  Number of like mu pair high pT : " << scalars->fScaler[3] << endl;
			cout << "      Number of single mu low pT : " << scalars->fScaler[4] << endl;
			cout << "     Number of single mu high pT : " << scalars->fScaler[5] << endl;
			cout << "                Hold (dead time) : " << scalars->fHold << endl;
			cout << "                      Spare word : " << scalars->fSpare << " (0x"
				<< setw(8) << setfill('0') << hex << scalars->fSpare <<
				setw(0) << setfill(fillChar) << dec << ")" << endl;
		}
		else
		{
			cout << "(None found)" << endl;
		}
	}
	
	
	void AliTriggerDecoderHandler::OnNewRegionalStruct(
			const AliMUONRegionalHeaderStruct* regionalStruct,
			const AliMUONRegionalScalarsStruct* scalars,
			const void* /*data*/
		)
	{
		TryDumpCorruptData(regionalStruct);
		
		fRegNum++;
		cout << "========================= Regional structure header "
			<< setw(3) << fRegNum << setw(0)
			<< " ==========================" << endl;
		char fillChar = cout.fill();  // store current fill char to restore it.
		cout << "Darc word bit pattern : 0x" << noshowbase << setw(8) << setfill('0')
			<< hex << regionalStruct->fDarcWord << setw(0)
			<< setfill(fillChar) << dec << endl;
		UShort_t errBits = GetRegionalErrorBits(regionalStruct);
		cout << "  [           Error type : "
			<< (((errBits >> 7) & 0x1) == 1 ? "1 (Physics) " : "0 (Software)")
			<< " ]" << endl;
		cout << "  [       Regional error : " << ((errBits >> 6) & 0x1) << "            ]" << endl;
		cout << "  [           Full error : " << ((errBits >> 6) & 0x1) << "            ]" << endl;
		cout << "  [          Empty error : " << ((errBits >> 6) & 0x1) << "            ]" << endl;
		cout << "  [ DARC L2 reject error : " << ((errBits >> 6) & 0x1) << "            ]" << endl;
		cout << "  [        DARC L2 error : " << ((errBits >> 6) & 0x1) << "            ]" << endl;
		cout << "  [        DARC L1 error : " << ((errBits >> 6) & 0x1) << "            ]" << endl;
		cout << "  [        DARC L0 error : " << ((errBits >> 6) & 0x1) << "            ]" << endl;
		
		cout << "  [  FPGA number in card : "
			<< (UInt_t)GetRegionalFPGANumber(regionalStruct) << " (";
		PrintBitPattern((UInt_t)GetRegionalFPGANumber(regionalStruct), 3);
		cout << "b)     ]" << endl;
		
		cout << "  [      Physics trigger : " << setw(13) << boolalpha << left
			<< GetRegionalDarcPhysFlag(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [     Regional present : " << setw(13) << boolalpha << left
			<< GetRegionalPresentFlag(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [         RAM not full : " << setw(13) << boolalpha << left
			<< GetRegionalRamNotFullFlag(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [        RAM not empty : " << setw(13) << boolalpha << left
			<< GetRegionalRamNotEmptyFlag(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [          L2 rejected : " << setw(13) << boolalpha << left
			<< GetRegionalL2RejStatus(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [          L2 accepted : " << setw(13) << boolalpha << left
			<< GetRegionalL2AccStatus(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [             L1 found : " << setw(13) << boolalpha << left
			<< GetRegionalL1Status(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [             L0 found : " << setw(13) << boolalpha << left
			<< GetRegionalL0Status(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [ No. of events in RAM : " << setw(13) << left
			<< (UInt_t)GetRegionalEventInRam(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [            Busy word : 0x"
			<< (UInt_t)GetRegionalBusy(regionalStruct) << " (";
		PrintBitPattern((UInt_t)GetRegionalBusy(regionalStruct), 4);
		cout << "b)  ]" << endl;
		
		cout << "Regional word bit pattern : 0x" << noshowbase << setw(8) << setfill('0')
			<< hex << regionalStruct->fWord << setw(0)
			<< setfill(fillChar) << dec << endl;
		cout << "  [ Physics trigger occurred : " << setw(27) << boolalpha << left
			<< (UInt_t)GetRegionalPhysFlag(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [         Number of resets : " << setw(27) << left
			<< (UInt_t)GetRegionalResetNb(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [ Controller Serial number : " << setw(27) << left
			<< (UInt_t)GetRegionalSerialNb(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [        Regional crate ID : " << setw(27) << left
			<< (UInt_t)GetRegionalId(regionalStruct) << setw(0)
			<< right << "]" << endl;
		cout << "  [    FPGA software version : " << setw(27) << left
			<< (UInt_t)GetRegionalVersion(regionalStruct) << setw(0)
			<< right << "]" << endl;
		UInt_t output = GetRegionalOutput(regionalStruct);
		cout << "  [          Regional output : 0x" << setw(2) << setfill('0') << hex
			<< output << dec << setw(0) << setfill(fillChar) << setw(0) << right << " (";
		PrintBitPattern(output, 8);
		cout  << "b)           ]" << endl;
		
		cout << "  [         High pT           |          Low pT           ]" << endl;
		cout << "  [ Single mu |    mu pair    | Single mu |    mu pair    ]" << endl;
		cout << "  [ -ve | +ve | unlike | like | -ve | +ve | unlike | like ]" << endl;
		cout << "  [-------------------------------------------------------]" << endl;
		cout << ((((output >> 7) & 0x1) == 1) ? "  [ yes" : "  [  no");
		cout << ((((output >> 6) & 0x1) == 1) ? "   yes" : "    no");
		cout << ((((output >> 5) & 0x1) == 1) ? "     yes " : "      no ");
		cout << ((((output >> 4) & 0x1) == 1) ? "    yes" : "     no");
		cout << ((((output >> 3) & 0x1) == 1) ? "   yes" : "    no");
		cout << ((((output >> 2) & 0x1) == 1) ? "   yes" : "    no");
		cout << ((((output >> 1) & 0x1) == 1) ? "     yes " : "      no ");
		cout << ((((output >> 0) & 0x1) == 1) ? "    yes ]" : "     no ]");
		cout << endl;
		
		cout << "Regional input (low pT) bit pattern  : 0x" << noshowbase << setw(8)
			<< setfill('0') << hex << regionalStruct->fInput[0]
			<< setw(0) << setfill(fillChar) << dec << endl;
		cout << "Regional input (high pT) bit pattern : 0x" << noshowbase << setw(8)
			<< setfill('0') << hex << regionalStruct->fInput[1]
			<< setw(0) << setfill(fillChar) << dec << endl;
		cout << "Regional input logical dump: " << endl;
		cout << "    Local     | Local | low pT mu | high pT mu" << endl;
		cout << "structure no. | board | -ve | +ve |  -ve | +ve" << endl;
		cout << "-----------------------------------------------" << endl;
		for (int i = 15; i >= 0 ; i--)
		{
			cout << setw(13) << (fLocNum + 16 - i) << setw(0) << "   ";
			cout << setw(5) << (15 - i) << setw(0) << "   ";
			cout << ((((regionalStruct->fInput[0] >> (i*2+1)) & 0x1) == 1) ? "yes" : " no");
			cout << ((((regionalStruct->fInput[0] >> (i*2+0)) & 0x1) == 1) ? "   yes" : "    no");
			cout << ((((regionalStruct->fInput[1] >> (i*2+1)) & 0x1) == 1) ? "    yes" : "     no");
			cout << ((((regionalStruct->fInput[1] >> (i*2+0)) & 0x1) == 1) ? "   yes" : "    no");
			cout << endl;
		}
		
		UInt_t mask = GetRegionalMask(regionalStruct);
		cout << "Local mask : 0x" << noshowbase << setw(4) << setfill('0') << hex
			<< mask << dec << setw(0) << setfill(fillChar) << " [Bits: ";
		PrintBitPattern(mask >> 12, 4); cout << " ";
		PrintBitPattern(mask >> 8, 4); cout << " ";
		PrintBitPattern(mask >> 4, 4); cout << " ";
		PrintBitPattern(mask >> 0, 4); cout << "]" << endl;
		
		cout << "L0 counter : " << GetRegionalL0(regionalStruct) << endl;
		
		cout << "========================= Regional structure scalars "
			<< setw(3) << fRegNum << setw(0)
			<< " =========================" << endl;
		if (scalars != NULL)
		{
			cout << "            Number of clock cycles : " << scalars->fClk << endl;
			cout << "   Number of high pT +ve single mu : " << scalars->fScaler[0] << endl;
			cout << "   Number of high pT -ve single mu : " << scalars->fScaler[1] << endl;
			cout << "Number of high pT unlike sign pair : " << scalars->fScaler[2] << endl;
			cout << "  Number of high pT like sign pair : " << scalars->fScaler[3] << endl;
			cout << "    Number of low pT +ve single mu : " << scalars->fScaler[4] << endl;
			cout << "    Number of low pT -ve single mu : " << scalars->fScaler[5] << endl;
			cout << " Number of low pT unlike sign pair : " << scalars->fScaler[6] << endl;
			cout << "   Number of low pT like sign pair : " << scalars->fScaler[7] << endl;
			cout << "                  Hold (dead time) : " << scalars->fHold << endl;
		}
		else
		{
			cout << "(None found)" << endl;
		}
	}
	
	
	void AliTriggerDecoderHandler::OnLocalStruct(
			const AliMUONLocalInfoStruct* localStruct,
			const AliMUONLocalScalarsStruct* scalars
		)
	{
		TryDumpCorruptData(localStruct);
		
		fLocNum++;
		cout << "=========================== Local structure header "
			<< setw(3) << fLocNum << setw(0)
			<< " ===========================" << endl;
		char fillChar = cout.fill();  // store current fill char to restore it.
		cout << "L0 strip patterns:" << endl;
		cout << "Chamber |        X         |         Y        " << endl;
		cout << "----------------------------------------------" << endl;
		cout << "   11     ";
		PrintBitPattern(GetLocalX1(localStruct), 16);
		cout << "   ";
		PrintBitPattern(GetLocalY1(localStruct), 16);
		cout << endl;
		cout << "   12     ";
		PrintBitPattern(GetLocalX2(localStruct), 16);
		cout << "   ";
		PrintBitPattern(GetLocalY2(localStruct), 16);
		cout << endl;
		cout << "   13     ";
		PrintBitPattern(GetLocalX3(localStruct), 16);
		cout << "   ";
		PrintBitPattern(GetLocalY3(localStruct), 16);
		cout << endl;
		cout << "   14     ";
		PrintBitPattern(GetLocalX4(localStruct), 16);
		cout << "   ";
		PrintBitPattern(GetLocalY4(localStruct), 16);
		cout << endl;
		
		cout << "L0 trigger bits: (word = ";
		cout << showbase << hex << localStruct->fTriggerBits
			<< noshowbase << dec << ")" << endl;
		cout << "            ID |  Dec | TrigY | YPos | Sign XDev | XDev |  XPos " << endl;
		cout << "----------------------------------------------------------------" << endl;
		cout << "Decimal : " << setw(4) << (UInt_t)GetLocalId(localStruct) << setw(0) << "   ";
		cout << setw(4) << (UInt_t)GetLocalDec(localStruct) << setw(0) << "   ";
		cout << setw(3) << (UInt_t)GetLocalTrigY(localStruct) << setw(0) << "     ";
		cout << setw(4) << (UInt_t)GetLocalYPos(localStruct) << setw(0) << "   ";
		cout << setw(5) << (UInt_t)GetLocalSXDev(localStruct) << setw(0) << "       ";
		cout << setw(4) << (UInt_t)GetLocalXDev(localStruct) << setw(0) << "   ";
		cout << setw(4) << (UInt_t)GetLocalXPos(localStruct) << setw(0) << endl;
		cout << " Binary : ";
		PrintBitPattern(AliHLTUInt32_t(GetLocalId(localStruct)), 4);
		cout << "   ";
		PrintBitPattern(AliHLTUInt32_t(GetLocalDec(localStruct)), 4);
		cout << "     ";
		PrintBitPattern(AliHLTUInt32_t(GetLocalTrigY(localStruct)), 1);
		cout << "     ";
		PrintBitPattern(AliHLTUInt32_t(GetLocalYPos(localStruct)), 4);
		cout << "       ";
		PrintBitPattern(AliHLTUInt32_t(GetLocalSXDev(localStruct)), 1);
		cout << "       ";
		PrintBitPattern(AliHLTUInt32_t(GetLocalXDev(localStruct)), 4);
		cout << "   ";
		PrintBitPattern(AliHLTUInt32_t(GetLocalXPos(localStruct)), 5);
		cout << endl;
		
		cout << "L0 decision (Dec): [high pT: ";
		PrintBitPattern((UInt_t)GetLocalHpt(localStruct), 2);
		cout << "b (";
		switch ((UInt_t)GetLocalHpt(localStruct))
		{
		case 0: cout << "No trigger"; break;
		case 1: cout << "-ve particle"; break;
		case 2: cout << "+ve particle"; break;
		case 3: cout << "No deviation trigger"; break;
		default: cout << "UNKNOWN"; break;
		}
		cout << "), low pT: ";
		PrintBitPattern((UInt_t)GetLocalLpt(localStruct), 2);
		cout << "b (";
		switch ((UInt_t)GetLocalLpt(localStruct))
		{
		case 0: cout << "No trigger"; break;
		case 1: cout << "-ve particle"; break;
		case 2: cout << "+ve particle"; break;
		case 3: cout << "No deviation trigger"; break;
		default: cout << "UNKNOWN"; break;
		}
		cout << ")]" << endl;
		
		cout << "=========================== Local structure scalars "
			<< setw(3) << fLocNum << setw(0)
			<< " ==========================" << endl;
		if (scalars != NULL)
		{
			cout << "              Number of L0 triggers : " << scalars->fL0 << endl;
			cout << "                   Hold (dead time) : " << scalars->fHold << endl;
			cout << "             Number of clock cycles : " << scalars->fClk << endl;
			cout << "       Number of low pT no triggers : " << scalars->fLPtNTrig << endl;
			cout << "      Number of high pT no triggers : " << scalars->fHPtNTrig << endl;
			cout << "    Number of low pT right triggers : " << scalars->fLPtRTrig << endl;
			cout << "   Number of high pT right triggers : " << scalars->fHPtRTrig << endl;
			cout << "     Number of low pT left triggers : " << scalars->fLPtLTrig << endl;
			cout << "    Number of high pT left triggers : " << scalars->fHPtLTrig << endl;
			cout << " Number of low pT straight triggers : " << scalars->fLPtSTrig << endl;
			cout << "Number of high pT straight triggers : " << scalars->fHPtSTrig << endl;
			
			UInt_t xoryflag = GetLocalComptXY(scalars);
			if (xoryflag == 1)
			{
				cout << "Y strip counts:" << endl;
			}
			else
			{
				cout << "X strip counts:" << endl;
			}
			cout << "      |               Chamber            " << endl;
			cout << "Strip |     11 |     12 |     13 |     14" << endl;
			cout << "------------------------------------------" << endl;
			for (int i = 0; i < 16; i++)
			{
				cout << setw(5) << i << setw(0) << "   "
					<< setw(6) << (UInt_t)GetLocalXY1(scalars, i) << setw(0) << "   "
					<< setw(6) << (UInt_t)GetLocalXY2(scalars, i) << setw(0) << "   "
					<< setw(6) << (UInt_t)GetLocalXY3(scalars, i) << setw(0) << "   "
					<< setw(6) << (UInt_t)GetLocalXY4(scalars, i) << setw(0) << endl;
			}
			
			cout << "    EOS word : 0x" << setw(8) << setfill('0')
				<< hex << scalars->fEOS << setw(0) << setfill(fillChar) << dec << endl;
			cout << "    [ Switches : 0x" << setw(3)
				<< setfill('0') << hex << (UInt_t)GetLocalSwitch(scalars)
				<< setw(0) << setfill(fillChar) << dec << " (";
			PrintBitPattern((UInt_t)GetLocalSwitch(scalars) >> 8, 2); cout << " ";
			PrintBitPattern((UInt_t)GetLocalSwitch(scalars) >> 4, 4); cout << " ";
			PrintBitPattern((UInt_t)GetLocalSwitch(scalars) >> 0, 4);
			cout << "b)    ]" << endl;
			cout << "    [   X or Y : " << xoryflag
				<< ((xoryflag == 1) ? " (scalars for Y strips) ]" : " (scalars for X strips) ]")
				<< endl;
			
			cout << "Reset signal : " << scalars->fReset << endl;
		}
		else
		{
			cout << "(None found)" << endl;
		}
	}
	
	
	void AliTriggerDecoderHandler::OnError(ErrorCode error, const void* location)
	{
		TryDumpCorruptData(location);
		HandleError(
			ErrorCodeToMessage(error), error,
			ErrorCodeToString(error), location
		);
	}

} // end of namespace


int DumpRawDataHeader(
		const char* buffer, unsigned long bufferSize, const AliRawDataHeader* header,
		bool continueParse
	)
{
	// Dumps the common DDL raw data block header.
	
	cout << "*************************** Common DDL data header *******************************" << endl;
	char fillChar = cout.fill();  // remember fill char to set back to original later.
	int result = CheckHeaderField(header->fSize, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Size of the raw data in bytes : ";
	if (header->fSize != 0xFFFFFFFF)
		cout << header->fSize;
	else
		cout << "0xFFFFFFFF (unknown)";
	cout << endl;
	
	result = CheckHeaderField(header->fWord2, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "               Format version : " << UInt_t(header->GetVersion()) << endl;
	UInt_t l1msg = header->GetL1TriggerMessage();
	cout << "           L1 trigger message : 0x"
		<< noshowbase << hex << setfill('0') << setw(2) << l1msg
		<< setw(0) << setfill(fillChar) << dec
		<< " [Spare: " << ((l1msg >> 7) & 0x1)
		<< ", ClT: " << ((l1msg >> 6) & 0x1)
		<< ", RoC: 0x" << noshowbase << hex << ((l1msg >> 2) & 0x4) << dec
		<< ", ESR: " << ((l1msg >> 1) & 0x1)
		<< ", L1SwC: " << ((l1msg >> 0) & 0x1) << "]" << endl;
	cout << "   Bunch crossing (Event ID1) : " << header->GetEventID1() << " (0x"
		<< noshowbase << hex << setfill('0') << setw(3) << header->GetEventID1()
		<< dec << setw(0) << setfill(fillChar) << ")" << endl;
	
	result = CheckHeaderField(header->fEventID2, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "     Orbit number (Event ID2) : " << header->fEventID2 << " (0x"
		<< noshowbase << hex << setfill('0') << setw(6) << header->fEventID2
		<< dec << setw(0) << setfill(fillChar) << ")" << endl;
	
	result = CheckHeaderField(header->fAttributesSubDetectors, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "             Block attributes : " << UInt_t(header->GetAttributes()) << " (0x"
		<< noshowbase << hex << setfill('0') << setw(2) << UInt_t(header->GetAttributes())
		<< dec << setw(0) << setfill(fillChar) << ")" << endl;
	cout << "  Participating sub-detectors : 0x" << noshowbase << hex
		<< setfill('0') << setw(6) << header->GetSubDetectors() << dec
		<< setw(0) << setfill(fillChar) << "    [Bits: ";
	for (int i = 5; i >= 0; i--)
	{
		PrintBitPattern(header->GetSubDetectors() >> (i*4), 4);
		if (i != 0)
			cout << " ";
		else
			cout << "]" << endl;
	}
	
	result = CheckHeaderField(header->fStatusMiniEventID, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	UInt_t statusBits = header->GetStatus();
	cout << "          Status & error bits : 0x" << noshowbase << hex
		<< setfill('0') << setw(4) << statusBits << setw(0) << setfill(fillChar)
		<< dec << endl;
	cout << "          [               Original data bit : " << ((statusBits >> 15) & 0x1) << " ]" << endl;
	cout << "          [        Multi-event buffer error : " << ((statusBits >> 14) & 0x1) << " ]" << endl;
	cout << "          [        Trigger L1 missing error : " << ((statusBits >> 13) & 0x1) << " ]" << endl;
	cout << "          [           Trigger error (other) : " << ((statusBits >> 12) & 0x1) << " ]" << endl;
	cout << "          [                 Pre-pulse error : " << ((statusBits >> 11) & 0x1) << " ]" << endl;
	cout << "          [       Trigger L2 time violation : " << ((statusBits >> 10) & 0x1) << " ]" << endl;
	cout << "          [       Trigger L1 time violation : " << ((statusBits >> 9) & 0x1) << " ]" << endl;
	cout << "          [                DDG payload flag : " << ((statusBits >> 8) & 0x1) << " ]" << endl;
	cout << "          [                HLT payload flag : " << ((statusBits >> 7) & 0x1) << " ]" << endl;
	cout << "          [               HLT decision flag : " << ((statusBits >> 6) & 0x1) << " ]" << endl;
	cout << "          [                       FEE error : " << ((statusBits >> 5) & 0x1) << " ]" << endl;
	cout << "          [ Trigger information unavailable : " << ((statusBits >> 4) & 0x1) << " ]" << endl;
	cout << "          [            Control parity error : " << ((statusBits >> 3) & 0x1) << " ]" << endl;
	cout << "          [               Data parity error : " << ((statusBits >> 2) & 0x1) << " ]" << endl;
	cout << "          [           Trigger missing error : " << ((statusBits >> 1) & 0x1) << " ]" << endl;
	cout << "          [           Trigger overlap error : " << ((statusBits >> 0) & 0x1) << " ]" << endl;
	cout << "                Mini event ID : " << header->GetMiniEventID() << " (0x"
		<< noshowbase << hex << setfill('0') << setw(3) << header->GetMiniEventID()
		<< dec << setw(0) << setfill(fillChar) << ")" << endl;
	
	result = CheckHeaderField(header->fTriggerClassLow, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	result = CheckHeaderField(header->fROILowTriggerClassHigh, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	ULong64_t triggerClasses = header->GetTriggerClasses();
	cout << "              Trigger classes : 0x" << noshowbase << hex << setw(13)
		<< setfill('0') << triggerClasses << setw(0) << setfill(fillChar)
		<< dec << endl;
	cout << "                [Bits: ";
	PrintBitPattern(triggerClasses >> (12*4), 2);
	cout << " ";
	for (int i = 11; i >= 0; i--)
	{
		PrintBitPattern(triggerClasses >> (i*4), 4);
		if (i != 0)
			cout << " ";
		else
			cout << "]" << endl;
	}
	
	result = CheckHeaderField(header->fROIHigh, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	ULong64_t roiBits = header->GetROI();
	cout << "           Region of interest : 0x" << noshowbase << hex << setw(9)
		<< setfill('0') << roiBits << setw(0) << setfill(fillChar)
		<< dec << endl;
	cout << "             [Bits: ";
	for (int i = 8; i >= 0; i--)
	{
		PrintBitPattern(roiBits >> (i*4), 4);
		if (i != 0)
			cout << " ";
		else
			cout << "]" << endl;
	}
	cout << "**********************************************************************************" <<endl;
	return EXIT_SUCCESS;
}


int DumpTrackerDDLRawStream(
		char* buffer, unsigned long bufferSize,
		bool continueParse, bool tryrecover
	)
{
	// Dumps a tracker DDL raw stream data.
	
	AliRawDataHeader* header = reinterpret_cast<AliRawDataHeader*>(buffer);
	int result = DumpRawDataHeader(buffer, bufferSize, header, continueParse);
	if (result != EXIT_SUCCESS) return result;

	// Setup the decoder for the DDL payload.
	AliMUONTrackerDDLDecoder<AliTrackerDecoderHandler> decoder;
	decoder.ExitOnError(not continueParse);
	decoder.SendDataOnParityError(true);
	decoder.TryRecover(tryrecover);
	decoder.AutoDetectTrailer(true);
	decoder.CheckForTrailer(true);
	char* payload = buffer + sizeof(AliRawDataHeader);
	UInt_t payloadSize = bufferSize - sizeof(AliRawDataHeader);
	if (decoder.Decode(payload, payloadSize))
	{
		return EXIT_SUCCESS;
	}
	else
	{
		return PARSE_ERROR;
	}
}


int DumpTriggerDDLRawStream(
		const char* buffer, unsigned long bufferSize,
		bool continueParse, bool tryrecover
	)
{
	// Dumps a trigger DDL raw stream data.
	
	const AliRawDataHeader* header =
		reinterpret_cast<const AliRawDataHeader*>(buffer);
	int result = DumpRawDataHeader(buffer, bufferSize, header, continueParse);
	if(result != EXIT_SUCCESS) return result;
	bool scalarEvent = ((header->GetL1TriggerMessage() & 0x1) == 0x1);
	
	AliMUONTriggerDDLDecoder<AliTriggerDecoderHandler> decoder;
	decoder.ExitOnError(not continueParse);
	decoder.TryRecover(tryrecover);
	decoder.AutoDetectScalars(false);
	const char* payload = buffer + sizeof(AliRawDataHeader);
	UInt_t payloadSize = bufferSize - sizeof(AliRawDataHeader);
	if (decoder.Decode(payload, payloadSize, scalarEvent))
	{
		return EXIT_SUCCESS;
	}
	else
	{
		return PARSE_ERROR;
	}
}


int DumpRecHitStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONRecHitStruct* hit,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	
	int result = CheckField(hit->fFlags, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	AliHLTUInt8_t chamber = 0xFF;
	AliHLTUInt16_t detElemId = 0xFFFF;
	AliHLTMUONUtils::UnpackRecHitFlags(hit->fFlags, chamber, detElemId);
	if (chamber == 0 and detElemId == 0)
	{
		cout << setw(10) << left << (int)(chamber) << setw(0);
		cout << setw(12) << left << (int)detElemId << setw(0);
	}
	else
	{
		cout << setw(10) << left << (int)(chamber+1) << setw(0);
		cout << setw(12) << left << (int)detElemId << setw(0);
	}
	
	result = CheckField(hit->fX, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << hit->fX << setw(0);

	result = CheckField(hit->fY, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(13) << left << hit->fY << setw(0);

	result = CheckField(hit->fZ, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << hit->fZ << setw(0) << endl;

	return result;
}


int DumpRecHitsBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	// Dumps a reconstructed hits data block.
	
	int result = EXIT_SUCCESS;
	AliHLTMUONRecHitsBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << "Chamber | DetElemID | X (cm)     | Y (cm)     | Z (cm)" << endl;
	cout << "-----------------------------------------------------------" << endl;
	const AliHLTMUONRecHitStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		int subResult = DumpRecHitStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpTriggerRecordStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONTriggerRecordStruct* record,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(record->fId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Trigger Record ID: " << record->fId <<endl;
	
	result = CheckField(record->fFlags, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	ios::fmtflags oldflags = cout.flags();
	cout << "Flags: " << showbase << hex << record->fFlags << dec;
	cout.flags(oldflags);
		
	// Print the individual trigger bits.
	AliHLTMUONParticleSign sign;
	bool hitset[4];
	AliHLTMUONUtils::UnpackTriggerRecordFlags(record->fFlags, sign, hitset);
	cout << " [Sign: " << sign << ", Hits set on chambers: ";
	bool first = true;
	for (AliHLTUInt32_t i = 0; i < 4; i++)
	{
		if (hitset[i])
		{
			cout << (first ? "" : ", ") << i+11;
			first = false;
		}
	}
	cout << (first ? "none]" : "]") << endl;

	result = CheckField(record->fPx, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Momentum: (px = " << record->fPx << ", ";

	result = CheckField(record->fPy, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "py = " << record->fPy << ", ";

	result = CheckField(record->fPz, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "pz = " << record->fPz << ") GeV/c"<<endl;
	
	cout << "Hits on chambers:" << endl;
	cout << "Chamber | DetElemID | X (cm)     | Y (cm)     | Z (cm)" << endl;
	cout << "-----------------------------------------------------------" << endl;
	const AliHLTMUONRecHitStruct* hit = &record->fHit[0];
	for(AliHLTUInt32_t ch = 0; ch < 4; ch++)
	{
		result = DumpRecHitStruct(buffer, bufferSize, hit++, continueParse);
		if (result != EXIT_SUCCESS) return result;
	}

	return result;

}


int DumpTriggerRecordsBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	// Dumps a trigger records data block.
	
	AliHLTMUONTriggerRecordsBlockReader block(buffer, bufferSize);
	
	int result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONTriggerRecordStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << "============================== Trigger Record number " << i+1
			<< " of " << nentries << " ==============================" << endl;
		int subResult = DumpTriggerRecordStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return EXIT_SUCCESS;
}


int DumpLocalStruct(
		const char* buffer, unsigned long bufferSize,
		const AliMUONLocalInfoStruct* localStruct,
		bool continueParse,
		const char* title = ""
	)
{
	// Prints the fields of a L0 local structure as found in the buffer.
	
	typedef AliMUONTriggerDDLDecoderEventHandler AliH;
	
	cout << "L0 strip patterns" << title << ":" << endl;
	cout << "Chamber |        X         |         Y        " << endl;
	cout << "----------------------------------------------" << endl;
	int result = CheckField(localStruct->fX2X1, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "   11     ";
	PrintBitPattern(AliH::GetLocalX1(localStruct), 16);
	result = CheckField(localStruct->fY2Y1, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "   ";
	PrintBitPattern(AliH::GetLocalY1(localStruct), 16);
	cout << endl;
	cout << "   12     ";
	PrintBitPattern(AliH::GetLocalX2(localStruct), 16);
	cout << "   ";
	PrintBitPattern(AliH::GetLocalY2(localStruct), 16);
	cout << endl;
	
	result = CheckField(localStruct->fX4X3, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "   13     ";
	PrintBitPattern(AliH::GetLocalX3(localStruct), 16);
	result = CheckField(localStruct->fY4Y3, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "   ";
	PrintBitPattern(AliH::GetLocalY3(localStruct), 16);
	cout << endl;
	cout << "   12     ";
	PrintBitPattern(AliH::GetLocalX4(localStruct), 16);
	cout << "   ";
	PrintBitPattern(AliH::GetLocalY4(localStruct), 16);
	cout << endl;
	
	cout << "L0 trigger bits" << title << ": (word = ";
	result = CheckField(localStruct->fTriggerBits, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << showbase << hex << localStruct->fTriggerBits
		<< noshowbase << dec << ")" << endl;
	cout << "  ID |  Dec | TrigY | YPos | Sign XDev | XDev |  XPos " << endl;
	cout << "------------------------------------------------------" << endl;
	PrintBitPattern(AliHLTUInt32_t(AliH::GetLocalId(localStruct)), 4);
	cout << "   ";
	PrintBitPattern(AliHLTUInt32_t(AliH::GetLocalDec(localStruct)), 4);
	cout << "     ";
	PrintBitPattern(AliHLTUInt32_t(AliH::GetLocalTrigY(localStruct)), 1);
	cout << "     ";
	PrintBitPattern(AliHLTUInt32_t(AliH::GetLocalYPos(localStruct)), 4);
	cout << "       ";
	PrintBitPattern(AliHLTUInt32_t(AliH::GetLocalSXDev(localStruct)), 1);
	cout << "       ";
	PrintBitPattern(AliHLTUInt32_t(AliH::GetLocalXDev(localStruct)), 4);
	cout << "   ";
	PrintBitPattern(AliHLTUInt32_t(AliH::GetLocalXPos(localStruct)), 5);
	cout << endl;
	
	return result;
}


int DumpTrigRecInfoStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONTrigRecInfoStruct* debuginfo,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	
	int result = CheckField(debuginfo->fTrigRecId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Trigger Record ID: " << debuginfo->fTrigRecId << endl;

	cout << "Detector element IDs:" << endl;
	cout << "  Chamber :           11 |           12 |           13 |           14" << endl;
	cout << "       ID : ";
	for (int i = 0; i < 4; i++)
	{
		result = CheckField(debuginfo->fDetElemId[i], buffer, bufferSize, continueParse);
		if (result != EXIT_SUCCESS) return result;
		ios::fmtflags oldflags = cout.flags();
		cout << setw(12) << right << debuginfo->fDetElemId[i] << setw(0);
		cout.flags(oldflags);
		if (i != 3) cout << " | ";
	}
	cout << endl;
	
	cout << "Momentum estimation parameters:" << endl;
	cout << "  Parameter : Z_middle coordinate (cm) | Integrated magnetic field (T.m)" << endl;
	cout << "      Value : ";
	result = CheckField(debuginfo->fZmiddle, buffer, bufferSize, continueParse);
	if(result != EXIT_SUCCESS) return result;
	cout << setw(24) << right << debuginfo->fZmiddle << setw(0) << " | ";
	
	result = CheckField(debuginfo->fBl, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(31) << right << debuginfo->fBl << setw(0) << endl;
	
	result = DumpLocalStruct(buffer, bufferSize, &debuginfo->fL0Struct, continueParse, " for central local structure");
	if (result != EXIT_SUCCESS) return result;
	result = DumpLocalStruct(buffer, bufferSize, &debuginfo->fL0StructPrev, continueParse, " for previous local structure");
	if (result != EXIT_SUCCESS) return result;
	result = DumpLocalStruct(buffer, bufferSize, &debuginfo->fL0StructNext, continueParse, " for next local structure");
	if (result != EXIT_SUCCESS) return result;
	
	return result;
}


int DumpTrigRecsDebugBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	// Dumps the debugging information for trigger records.
	
	AliHLTMUONTrigRecsDebugBlockReader block(buffer, bufferSize);
	
	int result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONTrigRecInfoStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << "======================= Trigger Record debug data " << i+1
			<< " of " << nentries << " =======================" << endl;
		int subResult = DumpTrigRecInfoStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return EXIT_SUCCESS;
}


int DumpClusterStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONClusterStruct* cluster,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(cluster->fId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "                             Cluster ID: " << cluster->fId << endl;

	result = CheckField(cluster->fDetElemId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "                    Detector Element ID: " << cluster->fDetElemId << endl;

	result = CheckField(cluster->fNchannelsB, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "    Number of channels in bending plane: " << cluster->fNchannelsB << endl;

	result = CheckField(cluster->fNchannelsNB, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Number of channels in non-bending plane: " << cluster->fNchannelsNB << endl;
	
	result = CheckField(cluster->fChargeB, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "                Charge on bending plane: " << cluster->fChargeB << endl;

	result = CheckField(cluster->fChargeNB, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "            Charge on non bending plane: " << cluster->fChargeNB << endl;

	cout << "Corresponding Hit: "<< endl;
	cout << "Chamber | DetElemID | X (cm)     | Y (cm)     | Z (cm)" << endl;
	cout << "-----------------------------------------------------------" << endl;
	result = DumpRecHitStruct(buffer, bufferSize, &cluster->fHit, continueParse);

	return result;
}


int DumpClustersBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	// Dumps a clusters block structure.
	
        int result = EXIT_SUCCESS;
	AliHLTMUONClustersBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONClusterStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << "======================= Cluster Number "
			<< i+1 << " =======================" << endl;
		int subResult = DumpClusterStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}	

	return result;
}


int DumpChannelStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONChannelStruct* channel,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(channel->fClusterId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	ios::fmtflags oldflags = cout.flags();
	cout << setw(16) << left << channel->fClusterId << setw(0);
	cout.flags(oldflags);

	result = CheckField(channel->fBusPatch, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(16) << left << channel->fBusPatch << setw(0);
	cout.flags(oldflags);

	result = CheckField(channel->fManu, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(16) << left << channel->fManu << setw(0);
	cout.flags(oldflags);

	result = CheckField(channel->fChannelAddress, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(16) << left << channel->fChannelAddress << setw(0);
	cout.flags(oldflags);

	result = CheckField(channel->fSignal, buffer, bufferSize, continueParse);
	if(result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(16) << left << channel->fSignal << setw(0);
	cout.flags(oldflags);

	result = CheckField(channel->fRawDataWord, buffer, bufferSize, continueParse);
	if(result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << showbase << hex << channel->fRawDataWord << dec << setw(0) << endl;
	cout.flags(oldflags);

	return result;
}


int DumpChannelsBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	// Dumps a channels block structure.
	
        int result = EXIT_SUCCESS;
	AliHLTMUONChannelsBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << "Cluster ID    | Bus Patch     | Manu Address  | Channel Addr  | Signal Value  | Raw Data Word " << endl;
	cout << "----------------------------------------------------------------------------------------------" << endl;
	const AliHLTMUONChannelStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{ 
		int subResult = DumpChannelStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	return EXIT_SUCCESS;
}


int DumpMansoTrackStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONMansoTrackStruct* track,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(track->fId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Track ID: " << track->fId << "\t";

	result = CheckField(track->fTrigRec, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Trigger Record ID: " << track->fTrigRec << endl;
	
	result = CheckField(track->fFlags, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	ios::fmtflags oldflags = cout.flags();
	cout << "Flags: " << showbase << hex << track->fFlags << dec;
	cout.flags(oldflags);
	
	// Print the individual trigger bits.
	AliHLTMUONParticleSign sign;
	bool hitset[4];
	AliHLTMUONUtils::UnpackMansoTrackFlags(track->fFlags, sign, hitset);
	cout << " [Sign: " << sign << ", Hits set on chambers: ";
	bool first = true;
	for (AliHLTUInt32_t i = 0; i < 4; i++)
	{
		if (hitset[i])
		{
			cout << (first ? "" : ", ") << i+7;
			first = false;
		}
	}
	cout << (first ? "none]" : "]") << endl;

	result = CheckField(track->fPx, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Momentum: (px = " << track->fPx << ", ";

	result = CheckField(track->fPy, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "py = " << track->fPy << ", ";

	result = CheckField(track->fPz, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "pz = " << track->fPz << ") GeV/c\t";

	result = CheckField(track->fChi2, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Chi squared fit: " << track->fChi2 << endl;
	
	cout << "Track hits:" << endl;
	cout << "Chamber | DetElemID | X (cm)     | Y (cm)     | Z (cm)" << endl;
	cout << "-----------------------------------------------------------" << endl;
	const AliHLTMUONRecHitStruct* hit = &track->fHit[0];
	for(AliHLTUInt32_t ch = 0; ch < 4; ch++)
	{
		result = DumpRecHitStruct(buffer, bufferSize, hit++, continueParse);
		if (result != EXIT_SUCCESS) return result;
	}

	return result;
}


int DumpMansoTracksBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	// Dumps the Manso tracks block structure.
	
	int result = EXIT_SUCCESS;
	AliHLTMUONMansoTracksBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONMansoTrackStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << "============================== Manso track number " << i+1
			<< " of " << nentries << " ==============================" << endl;
		int subResult = DumpMansoTrackStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpMansoRoIStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONMansoRoIStruct* roi,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(roi->fX, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	ios::fmtflags oldflags = cout.flags();
	cout << setw(13) << left << roi->fX << setw(0);
	cout.flags(oldflags);

	result = CheckField(roi->fY, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(13) << left << roi->fY << setw(0);
	cout.flags(oldflags);

	result = CheckField(roi->fZ, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(13) << left << roi->fZ << setw(0);
	cout.flags(oldflags);

	result = CheckField(roi->fRadius, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << roi->fRadius << setw(0) << endl;
	cout.flags(oldflags);

	return result;
}


int DumpMansoCandidateStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONMansoCandidateStruct* candidate,
		bool continueParse
	)
{
	// Dumps the manso candidate structure.
	
	int result = DumpMansoTrackStruct(buffer, bufferSize, &candidate->fTrack, continueParse);
	if (result != EXIT_SUCCESS) return result;
	
	cout << "Regions of interest:" << endl;
	cout << "Chamber | X (cm)     | Y (cm)     | Z (cm)     | Radius (cm)" << endl;
	cout << "-------------------------------------------------------------" << endl;
	const AliHLTMUONMansoRoIStruct* roi = &candidate->fRoI[0];
	for(AliHLTUInt32_t ch = 0; ch < 4; ch++)
	{
		cout << setw(10) << ch + 7;
		result = DumpMansoRoIStruct(buffer, bufferSize, roi++, continueParse);
		if (result != EXIT_SUCCESS) return result;
	}
	
	cout << "Momentum estimation parameters:" << endl;
	cout << "  Parameter : Z_middle coordinate (cm) | Integrated magnetic field (T.m)" << endl;
	cout << "      Value : ";
	result = CheckField(candidate->fZmiddle, buffer, bufferSize, continueParse);
	if(result != EXIT_SUCCESS) return result;
	ios::fmtflags oldflags = cout.flags();
	cout << setw(24) << right << candidate->fZmiddle << setw(0) << " | ";
	cout.flags(oldflags);
	
	result = CheckField(candidate->fBl, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(31) << right << candidate->fBl << setw(0) << endl;
	cout.flags(oldflags);
	
	return result;
}


int DumpMansoCandidatesBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	// Dumps the manso candidates block structure.
	
	int result = EXIT_SUCCESS;
	AliHLTMUONMansoCandidatesBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONMansoCandidateStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << "=========================== Manso track candidate number " << i+1
			<< " of " << nentries << " ===========================" << endl;
		int subResult = DumpMansoCandidateStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpTrackStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONTrackStruct* track,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(track->fId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Track ID: " << track->fId << "\t";

	result = CheckField(track->fTrigRec, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Trigger Record ID: " << track->fTrigRec << endl;
	
	result = CheckField(track->fFlags, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	ios::fmtflags oldflags = cout.flags();
	cout << "Flags: " << showbase << hex << track->fFlags << dec;
	cout.flags(oldflags);
	
	// Print the individual trigger bits.
	AliHLTMUONParticleSign sign;
	bool hitset[16];
	AliHLTMUONUtils::UnpackTrackFlags(track->fFlags, sign, hitset);
	cout << " [Sign: " << sign << ", Hits set: ";
	bool first = true;
	for (AliHLTUInt32_t i = 0; i < 16; i++)
	{
		if (hitset[i])
		{
			cout << (first ? "" : ", ") << i;
			first = false;
		}
	}
	cout << (first ? "none]" : "]") << endl;

	result = CheckField(track->fPx, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Momentum: (px = " << track->fPx << ", ";

	result = CheckField(track->fPy, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "py = " << track->fPy << ", ";

	result = CheckField(track->fPz, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "pz = " << track->fPz << ") GeV/c" << endl;
	
	result = CheckField(track->fInverseBendingMomentum, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Inverse bending momentum: " << track->fInverseBendingMomentum << " c/GeV" << endl;
	
	result = CheckField(track->fThetaX, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Slope: (non-bending plane = " << track->fThetaX << ", ";
	
	result = CheckField(track->fThetaY, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "bending plane = " << track->fThetaY << ")" << endl;

	result = CheckField(track->fX, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "DCA vertex: (x = " << track->fX << ", ";

	result = CheckField(track->fY, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "y = " << track->fY << ", ";

	result = CheckField(track->fZ, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "z = " << track->fZ << ") cm" << endl;

	result = CheckField(track->fChi2, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Chi squared fit: " << track->fChi2 << endl;
	
	cout << "Track hits:" << endl;
	cout << "Chamber | DetElemID | X (cm)     | Y (cm)     | Z (cm)" << endl;
	cout << "-----------------------------------------------------------" << endl;
	const AliHLTMUONRecHitStruct* hit = &track->fHit[0];
	for(AliHLTUInt32_t ch = 0; ch < 16; ch++)
	{
		result = DumpRecHitStruct(buffer, bufferSize, hit++, continueParse);
		if (result != EXIT_SUCCESS) return result;
	}

	return result;
}


int DumpTracksBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	// Dumps the full tracks block structure.
	
	int result = EXIT_SUCCESS;
	AliHLTMUONTracksBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	const AliHLTMUONTrackStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		cout << "================================ track number " << i+1
			<< " of " << nentries << " ================================" << endl;
		int subResult = DumpTrackStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpSinglesDecisionBlockHeader(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONSinglesDecisionBlockStruct* header,
		bool continueParse
	)
{
	// Step through the header fields trying to print them.
	// At each step check if we have not overflowed the buffer, if we have
	// not then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckHeaderField(header->fNlowPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << " Number of low pt triggers: " << header->fNlowPt << endl;
	
	result = CheckHeaderField(header->fNhighPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Number of high pt triggers: " << header->fNhighPt << endl;

	return result;
}


int DumpTrackDecisionStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONTrackDecisionStruct* decision,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(decision->fTrackId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	ios::fmtflags oldflags = cout.flags();
	cout << setw(13) << left << decision->fTrackId << setw(0);
	cout.flags(oldflags);
	
	result = CheckField(decision->fTriggerBits, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(12) << left << showbase << hex << decision->fTriggerBits
		<< setw(0) << dec;
		
	// Print the individual trigger bits.
	bool highPt, lowPt;
	AliHLTMUONUtils::UnpackTrackDecisionBits(decision->fTriggerBits, highPt, lowPt);
	cout << setw(7) << left << (highPt ? "yes" : "no");
	cout << setw(8) << left << (lowPt ? "yes" : "no");
	cout.flags(oldflags);
	
	result = CheckField(decision->fPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << setw(0) << decision->fPt << endl;

	return result;
}


int DumpSinglesDecisionBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	int result = EXIT_SUCCESS;
	AliHLTMUONSinglesDecisionBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	// Dump the rest of the block header.
	const AliHLTMUONSinglesDecisionBlockStruct* header = &block.BlockHeader();
	int subResult = DumpSinglesDecisionBlockHeader(buffer, bufferSize, header, continueParse);
	if (subResult != EXIT_SUCCESS) return subResult;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << "           |        Trigger Bits      |" << endl;
	cout << "Track ID   | Raw         HighPt LowPt | pT" << endl;
	cout << "----------------------------------------------------" << endl;
	const AliHLTMUONTrackDecisionStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		subResult = DumpTrackDecisionStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpPairsDecisionBlockHeader(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONPairsDecisionBlockStruct* header,
		bool continueParse
	)
{
	// Step through the header fields trying to print them.
	// At each step check if we have not overflowed the buffer, if we have
	// not then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckHeaderField(header->fNunlikeAnyPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "      Number of unlike all pt triggers: " << header->fNunlikeAnyPt << endl;
	
	result = CheckHeaderField(header->fNunlikeLowPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "      Number of unlike low pt triggers: " << header->fNunlikeLowPt << endl;
	
	result = CheckHeaderField(header->fNunlikeHighPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "     Number of unlike high pt triggers: " << header->fNunlikeHighPt << endl;
	
	result = CheckHeaderField(header->fNlikeAnyPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "        Number of like any pt triggers: " << header->fNlikeAnyPt << endl;
	
	result = CheckHeaderField(header->fNlikeLowPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "        Number of like low pt triggers: " << header->fNlikeLowPt << endl;
	
	result = CheckHeaderField(header->fNlikeHighPt, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "       Number of like high pt triggers: " << header->fNlikeHighPt << endl;
	
	result = CheckHeaderField(header->fNmassAny, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << " Number of all invariant mass triggers: " << header->fNmassAny << endl;
	
	result = CheckHeaderField(header->fNmassLow, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << " Number of low invariant mass triggers: " << header->fNmassLow << endl;
	
	result = CheckHeaderField(header->fNmassHigh, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Number of high invariant mass triggers: " << header->fNmassHigh << endl;
	
	return result;
}


int DumpPairDecisionStruct(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONPairDecisionStruct* decision,
		bool continueParse
	)
{
	// Step through the fields trying to print them.
	// At each step check if we have not overflowed the buffer. If we have
	// not, then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckField(decision->fTrackAId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	ios::fmtflags oldflags = cout.flags();
	cout << setw(13) << left << decision->fTrackAId << setw(0);
	cout.flags(oldflags);
	
	result = CheckField(decision->fTrackBId, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(13) << left << decision->fTrackBId << setw(0);
	cout.flags(oldflags);
	
	result = CheckField(decision->fTriggerBits, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	oldflags = cout.flags();
	cout << setw(12) << left << showbase << hex << decision->fTriggerBits
		<< setw(0) << dec;
	cout.flags(oldflags);
		
	// Print the individual trigger bits.
	bool highMass, lowMass, unlike;
	AliHLTUInt8_t highPtCount, lowPtCount;
	AliHLTMUONUtils::UnpackPairDecisionBits(
			decision->fTriggerBits,
			highMass, lowMass, unlike, highPtCount, lowPtCount
		);
	oldflags = cout.flags();
	cout << setw(7) << left << (highMass ? "yes" : "no");
	cout << setw(7) << left << (lowMass ? "yes" : "no");
	cout << setw(7) << left << (unlike ? "yes" : "no");
	cout << setw(6) << left << AliHLTUInt16_t(highPtCount);
	cout << setw(8) << left << AliHLTUInt16_t(lowPtCount);
	cout << setw(0);
	cout.flags(oldflags);
	
	result = CheckField(decision->fInvMass, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << decision->fInvMass << endl;
	
	return EXIT_SUCCESS;
}


int DumpPairsDecisionBlock(
		const char* buffer, unsigned long bufferSize,
		bool continueParse
	)
{
	int result = EXIT_SUCCESS;
	AliHLTMUONPairsDecisionBlockReader block(buffer, bufferSize);
	
	result = CheckCommonHeader(block, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS and not continueParse) return result;
	
	// Dump the rest of the block header.
	const AliHLTMUONPairsDecisionBlockStruct* header = &block.BlockHeader();
	int subResult = DumpPairsDecisionBlockHeader(buffer, bufferSize, header, continueParse);
	if (subResult != EXIT_SUCCESS) return subResult;
	
	AliHLTUInt32_t nentries = CalculateNEntries(block, bufferSize);
	
	// Print the data block record entries.
	cout << "           |            |                Trigger Bits                  |" << endl;
	cout << "Track A ID | Track B ID | Raw         HiMass LoMass Unlike HiPt# LoPt# | Invariant mass" << endl;
	cout << "----------------------------------------------------------------------------------------" << endl;
	const AliHLTMUONPairDecisionStruct* entry = block.GetArray();
	for(AliHLTUInt32_t i = 0; i < nentries; i++)
	{
		subResult = DumpPairDecisionStruct(buffer, bufferSize, entry++, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	}
	
	return result;
}


int DumpCommonHeader(
		const char* buffer, unsigned long bufferSize,
		const AliHLTMUONDataBlockHeader* header, bool continueParse
	)
{
	// Step through the header fields trying to print them.
	// At each step check if we have not overflowed the buffer, if we have
	// not then we can print the field, otherwise we print the left over
	// bytes assumed to be corrupted rubbish.
	int result = CheckHeaderField(header->fType, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	AliHLTMUONDataBlockType type = AliHLTMUONDataBlockType(header->fType);
	cout << "       Block type: " << type << endl;
	
	result = CheckHeaderField(header->fRecordWidth, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "     Record width: " << header->fRecordWidth << endl;
	
	result = CheckHeaderField(header->fNrecords, buffer, bufferSize, continueParse);
	if (result != EXIT_SUCCESS) return result;
	cout << "Number of entries: " << header->fNrecords << endl;
	
	return result;
}

/**
 * Method to look for a certain data word key in the buffer.
 * \param buffer  The start location in the buffer to search.
 * \param end  The end of the buffer.
 * \returns  the pointer position just past the word found or
 *     NULL if the word was not found in the buffer.
 */
const char* FindDataWord(const char* buffer, const char* end, UInt_t word)
{
	for (const char* current = buffer; (current+1) <= end; current += sizeof(UInt_t))
	{
		const UInt_t* currentWord = reinterpret_cast<const UInt_t*>(current);
		if (*currentWord == word) return current + sizeof(UInt_t);
	}
	return NULL;
}

/**
 * Method to check if the data buffer is really a raw DDL stream.
 * \returns  kUnknownDataBlock if this does not look like a raw DDL stream.
 *     kTrackerDDLRawData if this looks like a raw DDL stream from the tracker.
 *     kTriggerDDLRawData if this looks like a raw DDL stream from the trigger.
 */
int CheckIfDDLStream(const char* buffer, unsigned long bufferSize)
{
	if (bufferSize < sizeof(AliRawDataHeader)) return kUnknownDataBlock;
	
	const AliRawDataHeader* cdhHeader =
		reinterpret_cast<const AliRawDataHeader*>(buffer);
	
	// Keep scores of indicators / tests that show this is a raw DDL stream
	// either from the trigger or tracker. We will decide if the stream is
	// indeed a raw DDL stream if the largest of the two scores is above a
	// minimum threshold.
	int trackerScore = 0;
	int triggerScore = 0;
	
	if (cdhHeader->fSize == UInt_t(-1) or cdhHeader->fSize == bufferSize)
	{
		trackerScore++;
		triggerScore++;
	}
	
	if (cdhHeader->GetVersion() == 2)
	{
		trackerScore++;
		triggerScore++;
	}
	
	const char* payload = buffer + sizeof(AliRawDataHeader);
	const char* payloadEnd = buffer + bufferSize;
	
	typedef AliMUONTrackerDDLDecoder<AliMUONTrackerDDLDecoderEventHandler> AliTrkDecoder;
	typedef AliMUONTriggerDDLDecoder<AliMUONTriggerDDLDecoderEventHandler> AliTrgDecoder;
	
	// See if the DDL data has a payload with data word keys as expected by
	// AliMUONTrackerDDLDecoder.
	const char* current = payload;
	while ( (current = FindDataWord(current, payloadEnd, AliTrkDecoder::BlockDataKeyWord())) != NULL )
	{
		trackerScore++;
	}
	current = payload;
	while ( (current = FindDataWord(current, payloadEnd, AliTrkDecoder::DspDataKeyWord())) != NULL )
	{
		trackerScore++;
	}
	current = payload;
	while ( (current = FindDataWord(current, payloadEnd, AliTrkDecoder::BusPatchDataKeyWord())) != NULL )
	{
		trackerScore++;
	}
	current = payload;
	while ( (current = FindDataWord(current, payloadEnd, AliTrkDecoder::EndOfDDLWord())) != NULL )
	{
		trackerScore++;
	}
	
	// See if the DDL data has a payload with data word keys as expected by
	// AliMUONTriggerDDLDecoder.
	current = payload;
	while ( (current = FindDataWord(current, payloadEnd, AliTrgDecoder::EndOfDarcWord())) != NULL )
	{
		triggerScore++;
	}
	current = payload;
	while ( (current = FindDataWord(current, payloadEnd, AliTrgDecoder::EndOfGlobalWord())) != NULL )
	{
		triggerScore++;
	}
	current = payload;
	while ( (current = FindDataWord(current, payloadEnd, AliTrgDecoder::EndOfRegionalWord())) != NULL )
	{
		triggerScore++;
	}
	current = payload;
	while ( (current = FindDataWord(current, payloadEnd, AliTrgDecoder::EndOfLocalWord())) != NULL )
	{
		triggerScore++;
	}
	
	if (triggerScore > trackerScore)
	{
		if (triggerScore >= 6) return kTriggerDDLRawData;
	}
	else
	{
		if (trackerScore >= 6) return kTrackerDDLRawData;
	}
	return kUnknownDataBlock;
}

/**
 * Parses the buffer and prints the contents to screen.
 * \param [in] buffer  The pointer to the buffer to parse.
 * \param [in] bufferSize  The size of the buffer in bytes.
 * \param [in] continueParse  If specified then the we try to continue parsing the
 *           buffer as much as possible.
 * \param [in] tryrecover Indicates if the DDL decoders should have special
 *           recovery logic enabled.
 * \param [in,out] type  Initialy this should indicate the type of the data block
 *           or kUnknownDataBlock if not known. On exit it will be filled with
 *           the type of the data block as discovered by this routine if type
 *           initially contained kUnknownDataBlock.
 * \returns  The error code indicating the problem. EXIT_SUCCESS is returned
 *           on success.
 */
int ParseBuffer(
		char* buffer, unsigned long bufferSize,
		bool continueParse, bool tryrecover, int& type
	)
{
	assert( buffer != NULL );
	int result = EXIT_SUCCESS;
	int subResult = EXIT_FAILURE;
	
	// If the -type|-t option was not used in the command line then we need to
	// figure out what type of data block this is from the data itself.
	bool ddlStream = false;
	if (type == kUnknownDataBlock)
	{
		// First check if this is a raw DDL stream, if not then assume it is
		// some kind of internal dHLT raw data block.
		int streamType = CheckIfDDLStream(buffer, bufferSize);
		if (streamType == kTrackerDDLRawData or streamType == kTriggerDDLRawData)
		{
			type = streamType;
			ddlStream = true;
		}
	}
	else if (type == kTrackerDDLRawData or type == kTriggerDDLRawData)
	{
		ddlStream = true;
	}
	
	if (not ddlStream)
	{
		if (bufferSize < sizeof(AliHLTMUONDataBlockHeader))
		{
			cerr << "ERROR: The size of the file is too small to contain a"
				" valid data block." << endl;
			result = PARSE_ERROR;
			if (not continueParse) return result;
		}
		const AliHLTMUONDataBlockHeader* header =
			reinterpret_cast<const AliHLTMUONDataBlockHeader*>(buffer);
	
		subResult = DumpCommonHeader(buffer, bufferSize, header, continueParse);
		if (subResult != EXIT_SUCCESS) return subResult;
	
		
		// Check if the block type in the header corresponds to the type given
		// by the '-type' command line parameter. If they do not then print an
		// error or big fat warning message and force interpretation of the data
		// block with the type given by '-type'.
		AliHLTMUONDataBlockType headerType = AliHLTMUONDataBlockType(header->fType);
		
		if (type == kUnknownDataBlock)
		{
			// -type not used in the command line so just use what is given
			// by the data block header.
			type = headerType;
		}
		else if (type != headerType)
		{
			cerr << "WARNING: The data block header indicates a type"
				" different from what was specified on the command line."
				" The data could be corrupt."
				<< endl;
			cerr << "WARNING: The type value in the file is "
				<< showbase << hex << header->fType
				<< " (" << headerType << "), but on the command line it is "
				<< showbase << hex << int(type) << dec
				<< " (" << type << ")."
				<< endl;
			cerr << "WARNING: Will force the interpretation of the data block"
				" with a type of " << type << "." << endl;
		}
	}
	
	// Now we know what type the data block is supposed to be, so we can
	// dump it to screen with the appropriate dump routine.
	switch (type)
	{
	case kTrackerDDLRawData:
		subResult = DumpTrackerDDLRawStream(buffer, bufferSize, continueParse, tryrecover);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kTriggerDDLRawData:
		subResult = DumpTriggerDDLRawStream(buffer, bufferSize, continueParse, tryrecover);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kTriggerRecordsDataBlock:
		subResult = DumpTriggerRecordsBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kTrigRecsDebugDataBlock:
		subResult = DumpTrigRecsDebugBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kRecHitsDataBlock:
		subResult = DumpRecHitsBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kClustersDataBlock:
		subResult = DumpClustersBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kChannelsDataBlock:
		subResult = DumpChannelsBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kMansoTracksDataBlock:
		subResult = DumpMansoTracksBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kMansoCandidatesDataBlock:
		subResult = DumpMansoCandidatesBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kTracksDataBlock:
		subResult = DumpTracksBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kSinglesDecisionDataBlock:
		subResult = DumpSinglesDecisionBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	case kPairsDecisionDataBlock:
		subResult = DumpPairsDecisionBlock(buffer, bufferSize, continueParse);
		if (subResult != EXIT_SUCCESS) result = subResult;
		break;
	default :
		cout << "ERROR: Unknown data block type. Found a type number of "
			<< showbase << hex << int(type) << noshowbase << dec
			<< " (" << int(type) << ")." << endl;
		result = PARSE_ERROR;
	}
	
	return result;
}

/**
 * Convert the type code to a string.
 */
const char* TypeToString(int type)
{
	if (type == kTrackerDDLRawData or type == kTriggerDDLRawData)
	{
		static char str[kAliHLTComponentDataTypefIDsize+1];
		AliHLTComponentDataType t = AliHLTMUONConstants::DDLRawDataType();
		memcpy(&str, &t.fID, kAliHLTComponentDataTypefIDsize);
		// Must insert the NULL character to make this an ANSI C string.
		str[kAliHLTComponentDataTypefIDsize] = '\0';
		return &str[0];
	}
	else
	{
		return AliHLTMUONUtils::DataBlockTypeToString(AliHLTMUONDataBlockType(type));
	}
}

/**
 * Find the data specification from the filename and return it in string format.
 */
const char* TryDecodeDataSpec(const char* filename)
{
	TString name = filename;
	
	TRegexp re1("MUONTR[GK]_[0123456789]+\\.ddl$");
	Ssiz_t length = 0;
	Ssiz_t pos = re1.Index(name, &length);
	if (pos != kNPOS)
	{
		TString substr;
		for (Ssiz_t i = 0; i < length; i++)
			substr += name[pos+i];
		TRegexp re("[0123456789]+");
		pos = re.Index(substr, &length);
		TString num;
		for (Ssiz_t i = 0; i < length; i++)
			num += substr[pos+i];
		AliHLTUInt32_t spec = AliHLTMUONUtils::EquipIdToSpec(num.Atoi());
		static char strbuf[32];
		sprintf(&strbuf[0], "0x%8.8X", spec);
		return &strbuf[0];
	}
	
	TRegexp re2("_MUON\\:.+_0x[0123456789abscefABCDEF]+\\.dat$");
	pos = re2.Index(name, &length);
	if (pos != kNPOS)
	{
		TString substr;
		for (Ssiz_t i = 0; i < length; i++)
			substr += name[pos+i];
		TRegexp re("0x[0123456789abscefABCDEF]+");
		pos = re.Index(substr, &length);
		TString num;
		for (Ssiz_t i = 0; i < length; i++)
			num += substr[pos+i];
		static TString result = num;
		return result.Data();
	}
	
	return NULL;
}

namespace
{
	// CDB path and run number to use.
	const char* gCDBPath = "local://$ALICE_ROOT/OCDB";
	Int_t gRunNumber = 0;
}

/**
 * Performs basic data integrity checks of the data block using the
 * AliHLTMUONDataCheckerComponent.
 * \param [in] sys  The HLT system framework.
 * \param [in] filename  The name of the file containing the data block to check.
 * \param [in] type  Must indicate the type of the data block.
 * \param [in] dataspec The data specification of the data block. NULL if none.
 * \param [in] maxLogging  If set to true then full logging is turned on for AliHLTSystem.
 * \returns  The error code indicating the problem. EXIT_SUCCESS is returned
 *           on success.
 */
int CheckDataIntegrity(
		AliHLTSystem& sys, const char* filename, int type,
		const char* dataspec, bool maxLogging
	)
{
	if (maxLogging)
	{
		AliLog::SetGlobalLogLevel(AliLog::kMaxType);
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	else
	{
		AliLog::SetGlobalLogLevel(AliLog::kWarning);
		int level = kHLTLogWarning | kHLTLogError | kHLTLogFatal;
		sys.SetGlobalLoggingLevel(AliHLTComponentLogSeverity(level));
	}
	
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTMUON.so");
	
	// Setup the component parameter lists and then the components.
	TString dcparams = "-return_error -warn_on_unexpected_block -no_global_check -cdbpath ";
	dcparams += gCDBPath;
	dcparams += " -run ";
	dcparams += gRunNumber;
	TString fpparams = "-datatype '";
	fpparams += TypeToString(type);
	fpparams += "' 'MUON'";
	if (dataspec != NULL)
	{
		fpparams += " -dataspec ";
		fpparams += dataspec;
	}
	else
	{
		const char* spec = TryDecodeDataSpec(filename);
		if (spec != NULL)
		{
			fpparams += " -dataspec ";
			fpparams += spec;
		}
		else
		{
			dcparams += " -ignorespec";
		}
	}
	fpparams += " -datafile ";
	fpparams += filename;
	TString fpname = "filePublisher_";
	fpname += filename;
	TString dcname = "checker_";
	dcname += filename;
	
	if (maxLogging)
	{
		cout << "DEBUG: Using the following flags for FilePublisher: \""
			<< fpparams.Data() << "\""<< endl;
		cout << "DEBUG: Using the following flags for "
			<< AliHLTMUONConstants::DataCheckerComponentId()
			<< ": \"" << dcparams.Data() << "\""<< endl;
	}
	
	AliHLTConfiguration(fpname.Data(), "FilePublisher", NULL, fpparams.Data());
	AliHLTConfiguration checker(
			dcname.Data(), AliHLTMUONConstants::DataCheckerComponentId(),
			fpname.Data(), dcparams.Data()
		);
	
	// Build and run the HLT tasks.
	if (sys.BuildTaskList(dcname.Data()) != 0) return HLTSYSTEM_ERROR;
	if (maxLogging) sys.PrintTaskList();
	if (sys.Run() != 1) return HLTSYSTEM_ERROR;
	
	// Now clean up.
	if (sys.CleanTaskList() != 0) return HLTSYSTEM_ERROR;

	return EXIT_SUCCESS;
}


/**
 * The caller is responsible for freeing memory allocated for buffer with a call
 * to delete [] buffer.
 */
int ReadFile(const char* filename, char*& buffer, unsigned long& bufferSize)
{
	assert( filename != NULL );
	
	// Open the file and find its size.
	fstream file;
	file.open(filename, ios::in);
	if (not file)
	{
		cerr << "ERROR: Could not open the file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	file.seekg(0, ios::end);
	if (not file)
	{
		cerr << "ERROR: Could not seek in the file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	bufferSize = file.tellg();
	if (not file)
	{
		cerr << "ERROR: Could not get file size for the file: " <<
			filename << endl;
		return SYSTEM_ERROR;
	}
	file.seekg(0, ios::beg);
	if (not file)
	{
		cerr << "ERROR: Could not seek in the file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	
	// Allocate the memory for the file.
	try
	{
		buffer = new char[bufferSize];
	}
	catch (const std::bad_alloc&)
	{
		cerr << "ERROR: Out of memory. Tried to allocate " << bufferSize
			<< " bytes." << endl;
		return SYSTEM_ERROR;
	}
	
	file.read(buffer, bufferSize);
	if (not file)
	{
		delete [] buffer;
		buffer = NULL;
		bufferSize = 0;
		cerr << "ERROR: Could not read from file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	
	file.close();
	if (not file)
	{
		delete [] buffer;
		buffer = NULL;
		bufferSize = 0;
		cerr << "ERROR: Could not close the file: " << filename << endl;
		return SYSTEM_ERROR;
	}
	
	return EXIT_SUCCESS;
}

/**
 * Prints the command line usage of this program to standard error.
 */
void PrintUsage(bool asError = true)
{
	std::ostream& os = asError ? cerr : cout;
	os << "Usage: dHLTdumpraw [-help|-h] [-continue|-c] [-type|-t <typename>] [-check|-k]" << endl;
	os << "         [-debug|-d] [-dataspec|-s <number>] [-cdbpath|-p <url>] [-run|-r <number>]" << endl;
	os << "         <filename> [<filename> ...]" << endl;
	os << "Where <filename> is the name of a file containing a raw data block." << endl;
	os << "Options:" << endl;
	os << " -help | -h" << endl;
	os << "       Displays this message." << endl;
	os << " -continue | -c" << endl;
	os << "       If specified, the program will try to continue parsing the data block" << endl;
	os << "       as much as possible rather than stopping at the first error." << endl;
	os << " -recover | -e" << endl;
	os << "       If specified, special DDL decoder recovery logic is enabled to handle" << endl;
	os << "       corrupt DDL data." << endl;
	os << " -type | -t <typename>" << endl;
	os << "       Forces the contents of the subsequent files specified on the command" << endl;
	os << "       line to be interpreted as a specific type of data block." << endl;
	os << "       Where <typename> can be one of:" << endl;
	os << "         rawtracker - raw DDL stream from tracker chambers." << endl;
	os << "         rawtrigger - raw DDL stream from trigger chambers." << endl;
	os << "         trigrecs - trigger records data." << endl;
	os << "         trigrecsdebug - debugging information about trigger records." << endl;
	os << "         rechits - reconstructed hits data." << endl;
	os << "         channels - channel debugging information from hit reconstruction." << endl;
	os << "         clusters - cluster debugging information from hit reconstruction." << endl;
	os << "         mansotracks - partial tracks from Manso algorithm." << endl;
	os << "         mansocandidates - track candidates considered in the Manso algorithm." << endl;
	os << "         tracks - tracks from full tracker component." << endl;
	os << "         singlesdecision - trigger decisions for single tracks." << endl;
	os << "         pairsdecision - trigger decisions for track pairs." << endl;
	os << "         autodetect - the type of the data block will be automatically" << endl;
	os << "                      detected." << endl;
	os << " -check | -k" << endl;
	os << "       If specified then data integrity checks are performed on the raw data." << endl;
	os << "       Warnings and errors are printed as problems are found with the data, but" << endl;
	os << "       the data will still be converted into ROOT objects as best as possible." << endl;
	os << " -debug | -d" << endl;
	os << "       If specified, then the all debug messages are printed by the AliHLTSystem." << endl;
	os << "       This is only useful if experiencing problems with the -check|-k option." << endl;
	os << " -dataspec | -s <number>" << endl;
	os << "       When specified, then <number> is used as the data specification for the" << endl;
	os << "       data file that follows. This option is only useful with the -check|-k option." << endl;
	os << " -cdbpath | -p <url>" << endl;
	os << "       The path to the CDB to use when running with the -check | -k option." << endl;
	os << " -run | -r <number>" << endl;
	os << "       The run number to use when running with the -check | -k option." << endl;
}

/**
 * Parses the command line.
 * @param argc  Number of arguments as given in main().
 * @param argv  Array of arguments as given in main().
 * @param filenames  Pointer to buffer storing file name strings.
 * @param filetypes  Array that receives the type of the data block expected, i.e.
 *                   the value of the -type flag for the corresponding file.
 * @param dataspecs  Data specifications to use for the data files.
 * @param numOfFiles  Receives the number of file name strings that were found
 *                    and added to 'filenames'.
 * @param continueParse  Set to true if the user requested to continue to parse
 *                      after errors.
 * @param tryrecover Set to true if the user wants DDL decoder recovery logic enabled.
 * @param checkData  Set to true if data integrity checking was requested.
 * @param maxLogging  Set to true if maximal logging was requested.
 * @return  A status flag suitable for returning from main(), containing either
 *          EXIT_SUCCESS or CMDLINE_ERROR.
 */
int ParseCommandLine(
		int argc,
		char** argv,
		const char** filenames,
		int* filetypes,
		const char** dataspecs,
		int& numOfFiles,
		bool& continueParse,
		bool& tryrecover,
		bool& checkData,
		bool& maxLogging
	)
{
	numOfFiles = 0;
	continueParse = false;
	tryrecover = false;
	maxLogging = false;
	checkData = false;
	int currentType = kUnknownDataBlock;
	const char* currentDataSpec = NULL;
	bool pathSet = false;
	bool runSet = false;

	// Parse the command line.
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-help") == 0 or strcmp(argv[i], "-h") == 0)
		{
			PrintUsage(false);
			return EXIT_SUCCESS;
		}
		else if (strcmp(argv[i], "-continue") == 0 or strcmp(argv[i], "-c") == 0)
		{
			continueParse = true;
		}
		else if (strcmp(argv[i], "-recover") == 0 or strcmp(argv[i], "-e") == 0)
		{
			tryrecover = true;
		}
		else if (strcmp(argv[i], "-type") == 0 or strcmp(argv[i], "-t") == 0)
		{
			if (++i >= argc)
			{
				cerr << "ERROR: Missing a type specifier." << endl << endl;
				PrintUsage();
				return CMDLINE_ERROR;
			}
			// Now we need to parse the typename in the command line.
			if (strcmp(argv[i], "autodetect") == 0)
			{
				currentType = kUnknownDataBlock;
			}
			else if (strcmp(argv[i], "rawtracker") == 0)
			{
				currentType = kTrackerDDLRawData;
			}
			else if (strcmp(argv[i], "rawtrigger") == 0)
			{
				currentType = kTriggerDDLRawData;
			}
			else
			{
				currentType = AliHLTMUONUtils::ParseCommandLineTypeString(argv[i]);
				if (currentType == kUnknownDataBlock)
				{
					cerr << "ERROR: Invalid type name '" << argv[i]
						<< "' specified for argument " << argv[i-1]
						<< "." << endl << endl;
					PrintUsage();
					return CMDLINE_ERROR;
				}
			}
		}
		else if (strcmp(argv[i], "-debug") == 0 or strcmp(argv[i], "-d") == 0)
		{
			maxLogging = true;
		}
		else if (strcmp(argv[i], "-check") == 0 or strcmp(argv[i], "-k") == 0)
		{
			checkData = true;
		}
		else if (strcmp(argv[i], "-dataspec") == 0 or strcmp(argv[i], "-s") == 0)
		{
			if (++i >= argc)
			{
				cerr << "ERROR: Missing a specification number." << endl << endl;
				PrintUsage();
				return CMDLINE_ERROR;
			}
			
			char* errPos = NULL;
			unsigned long num = strtoul(argv[i], &errPos, 0);
			if (errPos == NULL or *errPos != '\0')
			{
				cerr << "ERROR: Cannot convert '%s' to a data specification number."
					<< argv[i] << endl;
				return CMDLINE_ERROR;
			}
			if (not AliHLTMUONUtils::IsSpecValid(num))
			{
				cerr << "ERROR: The data specification number is not a valid format."
					<< endl;
				return CMDLINE_ERROR;
			}
			currentDataSpec = argv[i];
		}
		else if (strcmp(argv[i], "-cdbpath") == 0 or strcmp(argv[i], "-p") == 0)
		{
			if (pathSet)
			{
				cerr << "WARNING: Already used -cdbpath|-p with '" << gCDBPath
					<< "' before. Will override it with the last value specified with -cdbpath|-p."
					<< endl;
			}
			if (++i >= argc)
			{
				cerr << "ERROR: Missing the URL for the CDB path." << endl << endl;
				PrintUsage();
				return CMDLINE_ERROR;
			}
			gCDBPath = argv[i];
			pathSet = true;
		}
		else if (strcmp(argv[i], "-run") == 0 or strcmp(argv[i], "-r") == 0)
		{
			if (runSet)
			{
				cerr << "WARNING: Already used -run|-r with " << gRunNumber
					<< " before. Will override it with the last value specified with -run|-r."
					<< endl;
			}
			if (++i >= argc)
			{
				cerr << "ERROR: Missing the run number." << endl << endl;
				PrintUsage();
				return CMDLINE_ERROR;
			}
			
			char* cpErr = NULL;
			Int_t run = Int_t( strtol(argv[i], &cpErr, 0) );
			if (cpErr == NULL or *cpErr != '\0' or run < 0)
			{
				cerr << "ERROR: Cannot convert '" << argv[i] << "' to a valid run number."
					" Expected a positive integer value." << endl;
				return CMDLINE_ERROR;
			}

			gRunNumber = run;
			runSet = true;
		}
		else
		{
			assert( numOfFiles < argc );
			filenames[numOfFiles] = argv[i];
			filetypes[numOfFiles] = currentType;
			dataspecs[numOfFiles] = currentDataSpec;
			currentDataSpec = NULL;  // Reset because '-dataspec' option is only valid for one file.
			numOfFiles++;
		}
	}
	
	// Now check that we have at least one filename and all the flags we need.
	if (numOfFiles == 0)
	{
		cerr << "ERROR: Missing a file name. You must specify at least one file to process."
			<< endl << endl;
		PrintUsage();
		return CMDLINE_ERROR;
	}
	
	return EXIT_SUCCESS;
}


int main(int argc, char** argv)
{
	// Test endianess of this machine during runtime and print warning if it is not
	// little endian.
	union
	{
		int dword;
		char byte[4];
	} endianTest;
	endianTest.dword = 0x1;
	if (endianTest.byte[0] != 0x1)
	{
		cerr << "!!!! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !!!!!" << endl;
		cerr << "!!!! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !!!!!" << endl;
		cerr << "!!                                                                     !!" << endl;
		cerr << "!! This is not a little endian machine, but dHLT raw data is normally  !!" << endl;
		cerr << "!! generated in little endian format. Unless you are looking at localy !!" << endl;
		cerr << "!! created simulated data, then this program will not show you correct !!" << endl;
		cerr << "!! output.                                                             !!" << endl;
		cerr << "!!                                                                     !!" << endl;
		cerr << "!!!! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !!!!!" << endl;
		cerr << "!!!! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !! WARNING !!!!!" << endl;
		cerr << endl;
	}


	int numOfFiles = 0;
	bool continueParse = false;
	bool tryrecover = false;
	bool checkData = false;
	bool maxLogging = false;
	int returnCode = EXIT_SUCCESS;
	char* buffer = NULL;
	const char** filename = NULL;
	int* filetype = NULL;
	const char** dataspec = NULL;

	try
	{
		// Reduce logging now to get rid of informationals from AliHLTSystem constructor.
		AliLog::SetGlobalLogLevel(AliLog::kWarning);
		AliHLTSystem sys;
	
		// There will be at least 'argc' number of filenames.
		typedef const char* AnsiString;
		filename = new AnsiString[argc];
		filetype = new int[argc];
		dataspec = new AnsiString[argc];
		
		returnCode = ParseCommandLine(
				argc, argv, filename, filetype, dataspec, numOfFiles,
				continueParse, tryrecover, checkData, maxLogging
			);

		if (returnCode == EXIT_SUCCESS)
		{
			for (int i = 0; i < numOfFiles; i++)
			{
				unsigned long bufferSize = 0;
				returnCode = ReadFile(filename[i], buffer, bufferSize);
				if (returnCode != EXIT_SUCCESS) break;
				if (numOfFiles > 1)
				{
					cout << "########## Start of dump for file: "
						<< filename[i] << " ##########" << endl;
				}
				int result = ParseBuffer(buffer, bufferSize, continueParse, tryrecover, filetype[i]);
				if (buffer != NULL) delete [] buffer;
				if (result != EXIT_SUCCESS)
				{
					returnCode = result;
					if (not continueParse) break;
				}
				if (checkData)
				{
					result = CheckDataIntegrity(
							sys, filename[i], filetype[i],
							dataspec[i], maxLogging
						);
					if (result != EXIT_SUCCESS)
					{
						returnCode = result;
						if (not continueParse) break;
					}
				}
				if (numOfFiles > 1)
				{
					cout << "##########   End of dump for file: " <<
						filename[i] << " ##########" << endl;
				}
			}
		}
		
		delete [] filename;
		delete [] filetype;
		delete [] dataspec;
	}
	catch (...)
	{
		cerr << "FATAL ERROR: An unknown exception occurred!" << endl << endl;
		returnCode = FATAL_ERROR;
		if (buffer != NULL) delete [] buffer;
		if (filename != NULL) delete [] filename;
		if (filetype != NULL) delete [] filetype;
		if (dataspec != NULL) delete [] dataspec;
	}
	
	return returnCode;
}
