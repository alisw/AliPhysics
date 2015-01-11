#ifndef ALIMUONTRIGGERDDLDECODEREVENTHANDLER_H
#define ALIMUONTRIGGERDDLDECODEREVENTHANDLER_H
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
/// \file   AliMUONTriggerDDLDecoderEventHandler.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   28-11-2007
/// \brief  Implementation of the high performance trigger DDL decoder event handler.
///
/// This file implementes the AliMUONTriggerDDLDecoderEventHandler class,
/// which is the callback interface for the AliMUONTriggerDDLDecoder decoder class.
///

#include <cassert>
#include <ostream>
#include <Rtypes.h>


// We use C binding for the structures because C is more uniform with its application
// binary interface (ABI) between compilers.
extern "C"
{

// The following structures are the headers found in the DDL payload coming from
// the muon hardware trigger. The specification is given in ALICE-INT-2005-012
// (https://edms.cern.ch/file/591904/1/ALICE-INT-2005-012.pdf)

/// The optional DARC board scalars.
struct AliMUONDarcScalarsStruct
{
	UInt_t     fL0R;       ///< DARC L0 received and used
	UInt_t     fL1P;       ///< DARC L1 physics
	UInt_t     fL1S;       ///< DARC L1 software
	UInt_t     fL2A;       ///< DARC L2 accept
	UInt_t     fL2R;       ///< DARC L2 reject
	UInt_t     fClk;       ///< DARC clock
	UInt_t     fHold;      ///< DARC hold (dead time)
	UInt_t     fSpare;     ///< DARC Empty slot (for the moment)
};

/// The global input and output words just after the DARC header.
struct AliMUONGlobalHeaderStruct
{
	UInt_t     fInput[4];    ///< Global input. 8-bit words comming from the each of the 16 regional controlers.
	UInt_t     fOutput;      ///< Global ouput
};

/// The optional global card scalars.
struct AliMUONGlobalScalarsStruct
{
	UInt_t     fL0;         ///< Global number of L0 triggers
	UInt_t     fClk;        ///< Global number of clock cycles
	UInt_t     fScaler[6];  ///< Global card ouput scalars.
	UInt_t     fHold;       ///< Global number of hold (dead time)
	UInt_t     fSpare;      ///< Global spare word
};

/// Regional header
struct AliMUONRegionalHeaderStruct
{
	UInt_t    fDarcWord;        ///< darc word
	UInt_t    fWord;            ///< first reg word
	UInt_t    fInput[2];        ///< regional input
	UInt_t    fL0CountAndMask;  ///< L0 counter (16 bits) and local mask ("poids faible" 16 bits)
};

/// Optional regional card scalars.
struct AliMUONRegionalScalarsStruct
{
	UInt_t     fClk;        ///< Regional number of clock cycles.
	UInt_t     fScaler[8];  ///< Regional ouput scalars.
	UInt_t     fHold;       ///< Regional hold (dead time)
};


/// Local card trigger information.
struct AliMUONLocalInfoStruct
{
	UInt_t fX2X1;  ///< 16 bits X2 position in 16 most significant bits and 16 bits of X1 in least significant bits.
	UInt_t fX4X3;  ///< 16 bits X4 position in 16 most significant bits and 16 bits of X3 in least significant bits.
	UInt_t fY2Y1;  ///< 16 bits Y2 position in 16 most significant bits and 16 bits of Y1 in least significant bits.
	UInt_t fY4Y3;  ///< 16 bits Y4 position in 16 most significant bits and 16 bits of Y3 in least significant bits.
	UInt_t fTriggerBits;  ///< Trigger bits and deviation.
};

/// Local card trigger scalars.
struct AliMUONLocalScalarsStruct
{
	UInt_t     fL0;        ///< local number of L0 triggers.
	UInt_t     fHold;      ///< local hold (dead time)
	UInt_t     fClk;       ///< local number of clock cycles
	
	UInt_t     fLPtNTrig;  ///< local low Pt no trigger
	UInt_t     fHPtNTrig;  ///< local high Pt no trigger
	
	UInt_t     fLPtRTrig;  ///< local low Pt right trigger
	UInt_t     fHPtRTrig;  ///< local high Pt right trigger
	
	UInt_t     fLPtLTrig;  ///< local low Pt left trigger
	UInt_t     fHPtLTrig;  ///< local high Pt left trigger
	
	UInt_t     fLPtSTrig;  ///< local low Pt straight trigger
	UInt_t     fHPtSTrig;  ///< local high Pt straight trigger
	
	UInt_t     fScaler[8*4];   ///< local data
	UInt_t     fEOS;           ///< contains switches conf. & flag for reading X (0) or Y (1) in fScaler
	UInt_t     fReset;         ///< reset signal
};

} // extern "C"


/// \ingroup raw
/// \class AliMUONTriggerDDLDecoderEventHandler
/// \brief Callback event handler class for the AliMUONTriggerDDLDecoder.
/// This class is the base class defining what methods the event handler for the
/// high performance decoder should have. This handler actually does nothing.
/// The user of this decoder will have to derive from this class a custom event
/// handler that actually does something within the callback methods OnNewRegionalHeader,
/// OnLocalStruct, OnError etc...
///
class AliMUONTriggerDDLDecoderEventHandler
{
public:

	/// The only reason for a virtual destructor is to make -Weffc++ shutup.
	/// This should not really be here since we do not actually want or need
	/// run-time polymorphism.
	virtual ~AliMUONTriggerDDLDecoderEventHandler() {}

	/// All the possible error codes from the parsing.
	enum ErrorCode
	{
		kNoError = 0,              /// Decoding was successful.
		kTooManyRegionals = 1,     /// Too many regional card structures are expected in the DDL payload.
		kNoDarcHeader = 2,         /// The DARC header is missing. The DDL buffer is too short to hold a DARC header.
		kNoDarcScalars = 3,        /// The DARC scalars are missing or corrupt. The DDL buffer is too short to contain them.
		kWrongEventType = 4,       /// Wrong event type obtained from the Darc header.
		kNoEndOfDarc = 5,          /// The DDL buffer is too short to contain an end of DARC header key word.
		kBadEndOfDarc = 6,         /// End of DARC header key word is incorrect or corrupt.
		kNoGlobalHeader = 7,       /// The global header is missing. The DDL buffer is too short to hold a global header.
		kNoGlobalScalars = 8,      /// The global scalars are missing or corrupt. The DDL buffer is too short to contain them.
		kNoEndOfGlobal = 9,        /// The DDL buffer is too short to contain an end of global header key word.
		kBadEndOfGlobal = 10,      /// End of global header key word is incorrect or corrupt.
		kNoRegionalHeader = 11,    /// The regional header is missing. The DDL buffer is too short to hold another regional header.
		kNoRegionalScalars = 12,   /// The regional scalars are missing or corrupt. The DDL buffer is too short to contain them.
		kNoEndOfRegional = 13,     /// The DDL buffer is too short to contain an end of regional header key word.
		kBadEndOfRegional = 14,    /// End of regional header key word is incorrect or corrupt.
		kNoLocalStruct = 15,       /// The local structure is missing. The DDL buffer is too short to hold another local structure.
		kNoLocalScalars = 16,      /// The local scalars are missing or corrupt. The DDL buffer is too short to contain them.
		kNoEndOfLocal = 17,        /// The DDL buffer is too short to contain an end of local structure key word.
		kBadEndOfLocal = 18,       /// End of local structure key word is incorrect or corrupt.
		kBufferTooBig = 19         /// The DDL raw data is larger than indicated by the headers; extra bytes are probably just garbage.
	};

	// The following methods should be overridden for specific processing to
	// take place in your own event handler.

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
	/// the last call to OnLocalStruct)
	/// The default behaviour of this method is to do nothing.
	/// - param const void*  The pointer to the start of the memory buffer storing
	///                the DDL payload.
	/// - param UInt_t The size in bytes of the memory buffer.
	void OnEndOfBuffer(const void* /*buffer*/, UInt_t /*bufferSize*/) {}
	
	/// The OnDarcHeader method will be called when the DARC header has been
	/// found in the DDL payload.
	/// The default behaviour of this method is to do nothing.
	/// - param UInt_t  The DARC header word as found in the payload.
	/// - param const AliMUONDarcScalarsStruct*  The DARC scalars found in the
	///       raw data. If there are no scalars in the data then this pointer
	///       is set to NULL.
	/// - param const void*  A pointer to the remainder of the raw data after
	///       the DARC header and scalars.
	void OnDarcHeader(
			UInt_t /*header*/,
			const AliMUONDarcScalarsStruct* /*scalars*/,
			const void* /*data*/
		)
	{
	}
	
	/// The OnGlobalHeader method will be called when the global header has
	/// been found in the DDL payload.
	/// The default behaviour of this method is to do nothing.
	/// - param const AliMUONGlobalHeaderStruct*  A pointer to the global header.
	/// - param const AliMUONDarcScalarsStruct*  The global scalars found in the
	///       raw data. If there are no scalars in the data then this pointer
	///       is set to NULL.
	/// - param const void*  A pointer to the start of the regional data blocks.
	void OnGlobalHeader(
			const AliMUONGlobalHeaderStruct* /*header*/,
			const AliMUONGlobalScalarsStruct* /*scalars*/,
			const void* /*data*/
		)
	{
	}
	
	/// The OnNewRegionalStruct method will be called for each regional header
	/// found in the DDL payload.
	/// The default behaviour of this method is to do nothing.
	/// - param const AliMUONRegionalHeaderStruct*  A pointer to the regional
	///       structure header.
	/// - param const AliMUONRegionalScalarsStruct*  The regional scalars found
	///       in the raw data. If there are no scalars in the data then this
	///       pointer is set to NULL.
	/// - param const void*  A pointer to the start of the local trigger
	///       structures data for this regional block.
	void OnNewRegionalStruct(
			const AliMUONRegionalHeaderStruct* /*regionalStruct*/,
			const AliMUONRegionalScalarsStruct* /*scalars*/,
			const void* /*data*/
		)
	{
	}
	
	/// The OnNewRegionalStructV2 method will be called for each regional header
	/// found in the DDL payload.
	/// The default behaviour of this method is to do nothing.
	/// This method is an alternate version to OnNewRegionalStruct and an
	/// inheriting class needs only implement one or the other version.
	/// - param UInt_t  The structure index number of the regional structure.
	/// - param const AliMUONRegionalHeaderStruct*  A pointer to the regional
	///       structure header.
	/// - param const AliMUONRegionalScalarsStruct*  The regional scalars found
	///       in the raw data. If there are no scalars in the data then this
	///       pointer is set to NULL.
	/// - param const void*  A pointer to the start of the local trigger
	///       structures data for this regional block.
	void OnNewRegionalStructV2(
			UInt_t /*num*/,
			const AliMUONRegionalHeaderStruct* /*regionalStruct*/,
			const AliMUONRegionalScalarsStruct* /*scalars*/,
			const void* /*data*/
		)
	{
	}
	
	/// The OnEndOfRegionalStruct method will be called whenever a regional
	/// structure has been processed. For each OnNewRegionalStruct method
	/// call a symmetric call to OnEndOfRegionalStruct is made after processing
	/// of the regional structure is done (after the last call to OnLocalStruct
	/// for that regional structure).
	/// - param const AliMUONRegionalHeaderStruct*  A pointer to the regional
	///       structure header.
	/// - param const AliMUONRegionalScalarsStruct*  The regional scalars found
	///       in the raw data. If there are no scalars in the data then this
	///       pointer is set to NULL.
	/// - param const void*  A pointer to the start of the local trigger
	///       structures data for this regional block.
	void OnEndOfRegionalStruct(
			const AliMUONRegionalHeaderStruct* /*regionalStruct*/,
			const AliMUONRegionalScalarsStruct* /*scalars*/,
			const void* /*data*/
		)
	{
	}
	
	/// The OnEndOfRegionalStructV2 method will be called whenever a regional
	/// structure has been processed. For each OnNewRegionalStruct method
	/// call a symmetric call to OnEndOfRegionalStruct is made after processing
	/// of the regional structure is done (after the last call to OnLocalStruct
	/// for that regional structure).
	/// This method is an alternate version to OnEndOfRegionalStruct and an
	/// inheriting class needs only implement one or the other version.
	/// - param UInt_t  The structure index number of the regional structure.
	/// - param const AliMUONRegionalHeaderStruct*  A pointer to the regional
	///       structure header.
	/// - param const AliMUONRegionalScalarsStruct*  The regional scalars found
	///       in the raw data. If there are no scalars in the data then this
	///       pointer is set to NULL.
	/// - param const void*  A pointer to the start of the local trigger
	///       structures data for this regional block.
	void OnEndOfRegionalStructV2(
			UInt_t /*num*/,
			const AliMUONRegionalHeaderStruct* /*regionalStruct*/,
			const AliMUONRegionalScalarsStruct* /*scalars*/,
			const void* /*data*/
		)
	{
	}
	
	/// The OnLocalStruct method will be called for each local trigger
	/// structure found in the DDL payload. The user must overload this
	/// method to process the local structures as needed.
	/// The default behaviour of this method is to do nothing.
	/// - param const AliMUONRegionalHeaderStruct*  A pointer to the local
	///       trigger structure found.
	/// - param const AliMUONRegionalScalarsStruct*  The local scalars found
	///       in the raw data. If there are no scalars in the data then this
	///       pointer is set to NULL.
	void OnLocalStruct(
			const AliMUONLocalInfoStruct* /*localStruct*/,
			const AliMUONLocalScalarsStruct* /*scalars*/
		)
	{
	}
	
	/// The OnLocalStructV2 method will be called for each local trigger
	/// structure found in the DDL payload. The user must overload this
	/// method to process the local structures as needed.
	/// The default behaviour of this method is to do nothing.
	/// This method is an alternate version to OnLocalStruct and an
	/// inheriting class needs only implement one or the other version.
	/// - param UInt_t  The structure index number of the local structure.
	/// - param const AliMUONRegionalHeaderStruct*  A pointer to the local
	///       trigger structure found.
	/// - param const AliMUONRegionalScalarsStruct*  The local scalars found
	///       in the raw data. If there are no scalars in the data then this
	///       pointer is set to NULL.
	void OnLocalStructV2(
			UInt_t /*num*/,
			const AliMUONLocalInfoStruct* /*localStruct*/,
			const AliMUONLocalScalarsStruct* /*scalars*/
		)
	{
	}
	
	// The following static methods are helper routines for decoding the
	// DARC header bits.
	
	/// Return event type
	/// \param header  Should be the header as given by the OnDarkHeader() method.
	static UChar_t GetDarcEventType(UInt_t header) { return (UChar_t)(header >> 30) & 0x3; };
	
	/// Return Darc type
	/// \param header  Should be the header as given by the OnDarkHeader() method.
	static UChar_t GetDarcType(UInt_t header) { return (UChar_t)(header >> 24) & 0x7; }
	
	/// Return serial number
	/// \param header  Should be the header as given by the OnDarkHeader() method.
	static UChar_t GetDarcSerialNb(UInt_t header) { return (UChar_t)(header >> 20) & 0xF; }
	
	/// Return version
	/// \param header  Should be the header as given by the OnDarkHeader() method.
	static UChar_t GetDarcVersion(UInt_t header) { return (UChar_t)(header >> 12) & 0xFF; }
	
	/// Return VME trig
	/// \param header  Should be the header as given by the OnDarkHeader() method.
	static bool GetDarcVMETrig(UInt_t header) { return (header & 0x800); }
	
	/// Return global flag
	/// \param header  Should be the header as given by the OnDarkHeader() method.
	static bool GetDarcGlobalFlag(UInt_t header) { return (header & 0x400); }
	
	/// Return CPT trigger
	/// \param header  Should be the header as given by the OnDarkHeader() method.
	static bool GetDarcCTPTrig(UInt_t header) { return (header & 0x200); }
	
	/// Return DAQ flag
	/// \param header  Should be the header as given by the OnDarkHeader() method.
	static bool GetDarcDAQFlag(UInt_t header) { return (header & 0x100); }
	
	/// Return reg pattern
	/// \param header  Should be the header as given by the OnDarkHeader() method.
	static UChar_t GetDarcRegPattern(UInt_t header) { return (UChar_t)(header & 0xFF); }
	
	// The following static methods are helper routines for decoding the
	// regional structure header bits.
	
	/// Return L0
	static UShort_t GetRegionalL0(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fL0CountAndMask >> 16) & 0xFFFF;
	}
	
	/// Return mask
	static UShort_t GetRegionalMask(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return header->fL0CountAndMask & 0xFFFF;
	}
	
	/// Return RegPhysFlag
	static bool GetRegionalPhysFlag(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fWord & 0x80000000) != 0;
	}
	
	/// Return ResetNb
	static UChar_t GetRegionalResetNb(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (UChar_t)(header->fWord >> 25) &  0x3F;
	}
	
	/// Return SerialNb
	static UChar_t GetRegionalSerialNb(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (UChar_t)(header->fWord >> 20) &  0x1F;
	}
	
	/// Return Id
	static UChar_t GetRegionalId(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (UChar_t)(header->fWord >> 16) &  0x0F;
	}
	
	/// Return Version
	static UChar_t GetRegionalVersion(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (UChar_t)(header->fWord >> 8)  &  0xFF;
	}
	
	/// Return Output
	static UChar_t GetRegionalOutput(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (UChar_t)(header->fWord       &  0xFF);
	}
	
	/// Return ErrorBits
	static UShort_t GetRegionalErrorBits(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (UShort_t)(header->fDarcWord >> 22) &  0x3FF;
	}
	
	/// Return FPGANumber
	static UChar_t GetRegionalFPGANumber(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (UChar_t)  (header->fDarcWord >> 19) &  0x7;
	}
	
	/// Return DarcPhysFlag
	static bool GetRegionalDarcPhysFlag(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fDarcWord  &  0x8000) != 0;
	}
	
	/// Return PresentFlag
	static bool GetRegionalPresentFlag(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fDarcWord  &  0x4000) != 0;
	}
	
	/// Return RamNotFullFlag
	static bool GetRegionalRamNotFullFlag(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fDarcWord  &  0x2000) != 0;
	}
	
	/// Return RamNotEmptyFlag
	static bool GetRegionalRamNotEmptyFlag(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fDarcWord  &  0x1000) != 0;
	}
	
	/// Return L2RejStatus
	static bool GetRegionalL2RejStatus(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fDarcWord  &  0x800) != 0;
	}
	
	/// Return L2AccStatus
	static bool GetRegionalL2AccStatus(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fDarcWord  &  0x400) != 0;
	}
	
	/// Return L1Status
	static bool GetRegionalL1Status(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fDarcWord  &  0x200) != 0;
	}
	
	/// Return L0Status
	static bool GetRegionalL0Status(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (header->fDarcWord  &  0x100) != 0;
	}
	
	/// Return EventInRam
	static UChar_t GetRegionalEventInRam(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (UChar_t)  (header->fDarcWord >> 4)  &  0x4;
	}
	
	/// Return Busy
	static UChar_t GetRegionalBusy(const AliMUONRegionalHeaderStruct* header)
	{
		assert( header != NULL );
		return (UChar_t)  (header->fDarcWord)       &  0x4;
	}
	
	// The following static methods are helper routines for decoding the
	// global header bits.
	
	/// Return global output
	/// \param header  Should be the header as given by the OnGlobalHeader() method.
	static UChar_t GetGlobalOutput(const AliMUONGlobalHeaderStruct* header)
	{
		assert(header != NULL);
		return header->fOutput & 0xFF;
	}
	
	/// Return global config
	/// \param header  Should be the header as given by the OnGlobalHeader() method.
	static UShort_t GetGlobalConfig(const AliMUONGlobalHeaderStruct* header)
	{
		assert(header != NULL);
		return (header->fOutput >> 16) & 0xFFFF;
	}
	
	// The following static methods are helper routines for decoding the
	// local trigger structure and scalar bits.
	
	/// Return X2
	static UShort_t GetLocalX2(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fX2X1 >> 16) & 0xFFFF;
	}
	
	/// Return X1
	static UShort_t GetLocalX1(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fX2X1) & 0xFFFF;
	}
	
	/// Return X4
	static UShort_t GetLocalX4(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fX4X3 >> 16) & 0xFFFF;
	}
	
	/// Return X3
	static UShort_t GetLocalX3(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fX4X3) & 0xFFFF;
	}

	/// Return Y2
	static UShort_t GetLocalY2(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fY2Y1 >> 16) & 0xFFFF;
	}
	
	/// Return Y1
	static UShort_t GetLocalY1(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fY2Y1) & 0xFFFF;
	}
	
	/// Return Y4
	static UShort_t GetLocalY4(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fY4Y3 >> 16) & 0xFFFF;
	}
	
	/// Return Y3
	static UShort_t GetLocalY3(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fY4Y3) & 0xFFFF;
	}
	
	/// Return Id
	static UChar_t GetLocalId(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return local->fTriggerBits >> 19 &  0xF;
	}
	
	/// Return Dec
	static UChar_t GetLocalDec(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return local->fTriggerBits >> 15 &  0xF;
	}
	
	/// Return TrigY
	static bool GetLocalTrigY(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fTriggerBits >> 14 & 0x1);
	}
	
	/// Return TriggerY
	static bool GetLocalTriggerY(const AliMUONLocalInfoStruct* local)
	{
		return not (GetLocalTrigY(local) and GetLocalYPos(local) == 15);
	}
	
	/// Return Upos
	static UChar_t GetLocalYPos(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return local->fTriggerBits >> 10 & 0xF;
	}
	
	/// Get Sign of X deviation 
	static bool GetLocalSXDev(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return (local->fTriggerBits >> 9 & 0x1);
	}
	
	/// Get X deviation 
	static UChar_t GetLocalXDev(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return local->fTriggerBits >> 5 & 0xF;
	}
	
	/// Return TriggerX
	static bool GetLocalTriggerX(const AliMUONLocalInfoStruct* local)
	{
		return not (GetLocalSXDev(local) and (GetLocalXDev(local) == 0)
		            and GetLocalXPos(local) == 0);
	}
	
	/// Return Xpos
	static UChar_t GetLocalXPos(const AliMUONLocalInfoStruct* local)
	{
		assert(local != NULL);
		return local->fTriggerBits & 0x1F;
	}

	/// Return LPT
	static UChar_t GetLocalLpt(const AliMUONLocalInfoStruct* local) {return (GetLocalDec(local) & 0x3);}
	
	/// Return HPT
	static UChar_t GetLocalHpt(const AliMUONLocalInfoStruct* local) {return (GetLocalDec(local) >> 2) & 0x3;}
	
	/// Return switch
	static UShort_t GetLocalSwitch(const AliMUONLocalScalarsStruct* scalars)
	{
		assert(scalars != NULL);
		return (scalars->fEOS >> 1) & 0x3FF;
	}
	
	/// Return ComptXY
	static UChar_t GetLocalComptXY(const AliMUONLocalScalarsStruct* scalars)
	{
		assert(scalars != NULL);
		return scalars->fEOS & 0x1;
	}
	
	/// Return XY1
	static UShort_t GetLocalXY1(const AliMUONLocalScalarsStruct* scalars, UInt_t n)
	{
		assert(scalars != NULL and n < 16);
		return  (n % 2 == 1) ? (scalars->fScaler[(n/2)] & 0xFFFF)
		                     : ((scalars->fScaler[(n/2)] >> 16) &  0xFFFF);
	}
	
	/// Return XY2
	static UShort_t GetLocalXY2(const AliMUONLocalScalarsStruct* scalars, UInt_t n)
	{
		assert(scalars != NULL and n < 16);
		return  (n % 2 == 1) ? (scalars->fScaler[8 + (n/2)] & 0xFFFF)
		                     : ((scalars->fScaler[8 + (n/2)] >> 16) &  0xFFFF);
	}
	
	/// Return XY3
	static UShort_t GetLocalXY3(const AliMUONLocalScalarsStruct* scalars, UInt_t n)
	{
		assert(scalars != NULL and n < 16);
		return  (n % 2 == 1) ? (scalars->fScaler[8*2 + (n/2)] & 0xFFFF)
		                     : ((scalars->fScaler[8*2 + (n/2)] >> 16) &  0xFFFF);
	}
	
	/// Return XY4
	static UShort_t GetLocalXY4(const AliMUONLocalScalarsStruct* scalars, UInt_t n)
	{
		assert(scalars != NULL and n < 16);
		return  (n % 2 == 1) ? (scalars->fScaler[8*3 + (n/2)] & 0xFFFF)
		                     : ((scalars->fScaler[8*3 + (n/2)] >> 16) &  0xFFFF);
	}
	
	/// Whenever a parsing error of the DDL payload is encountered because of
	/// corruption of the raw data the OnError method is called immediately at
	/// the point this error is discovered.
	/// The default behaviour of this method is to do nothing.
	/// -param error  This is an error code indicating the kind of problem
	///               encountered with the DDL payload.
	/// -param location  This is a pointer into the DDL payload memory buffer
	///         indicating the exact location where the parsing error happened
	///         or i.e. the location of the corruption.
	/// Note that a relative offset in bytes from the start of the memory buffer
	/// can be calculated by: storing the buffer pointer recevied in OnNewBuffer
	/// earlier in fBufferStart for example, and then the offset is given by:
	///   offset = (unsigned long)location - (unsigned long)fBufferStart;
	void OnError(ErrorCode /*error*/, const void* /*location*/) {}
	
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

inline const char* AliMUONTriggerDDLDecoderEventHandler::ErrorCodeToString(ErrorCode code)
{
	/// This is a utility method which converts an error code to a string
	/// representation for printing purposes.
	/// \param code  The error code as received in OnError for example.
	/// \return  An ANSI string containing the name of the error code symbol.
	
	switch (code)
	{
	case kNoError: return "kNoError";
	case kTooManyRegionals: return "kTooManyRegionals";
	case kNoDarcHeader: return "kNoDarcHeader";
	case kNoDarcScalars: return "kNoDarcScalars";
	case kWrongEventType: return "kWrongEventType";
	case kNoEndOfDarc: return "kNoEndOfDarc";
	case kBadEndOfDarc: return "kBadEndOfDarc";
	case kNoGlobalHeader: return "kNoGlobalHeader";
	case kNoGlobalScalars: return "kNoGlobalScalars";
	case kNoEndOfGlobal: return "kNoEndOfGlobal";
	case kBadEndOfGlobal: return "kBadEndOfGlobal";
	case kNoRegionalHeader: return "kNoRegionalHeader";
	case kNoRegionalScalars: return "kNoRegionalScalars";
	case kNoEndOfRegional: return "kNoEndOfRegional";
	case kBadEndOfRegional: return "kBadEndOfRegional";
	case kNoLocalStruct: return "kNoLocalStruct";
	case kNoLocalScalars: return "kNoLocalScalars";
	case kNoEndOfLocal: return "kNoEndOfLocal";
	case kBadEndOfLocal: return "kBadEndOfLocal";
	case kBufferTooBig: return "kBufferTooBig";
	default: return "INVALID";
	}
}


inline const char* AliMUONTriggerDDLDecoderEventHandler::ErrorCodeToMessage(ErrorCode code)
{
	/// This is a utility method which converts an error code to user friendly
	/// descriptive message useful for printing to the screen.
	/// \param code  The error code as received in OnError for example.
	/// \return  An ANSI string containing a descriptive message of the error.
	
	switch (code)
	{
	case kNoError:
		return "Decoding was successful.";
	case kTooManyRegionals:
		return "Too many regional card structures are expected in the DDL payload.";
	case kNoDarcHeader:
		return "The DARC header is missing. The DDL buffer is too short"
			" to hold a DARC header.";
	case kNoDarcScalars:
		return "The DARC scalars are missing or corrupt."
			" The DDL buffer is too short to contain them.";
	case kWrongEventType:
		return "Wrong event type obtained from the Darc header.";
	case kNoEndOfDarc:
		return "The DDL buffer is too short to contain an end of DARC"
			" header key word.";
	case kBadEndOfDarc:
		return "End of DARC header key word is incorrect or corrupt.";
	case kNoGlobalHeader:
		return "The global header is missing. The DDL buffer is too"
			" short to hold a global header.";
	case kNoGlobalScalars:
		return "The global scalars are missing or corrupt. The DDL"
			" buffer is too short to contain them.";
	case kNoEndOfGlobal:
		return "The DDL buffer is too short to contain an end of global"
			" header key word.";
	case kBadEndOfGlobal:
		return "End of global header key word is incorrect or corrupt.";
	case kNoRegionalHeader:
		return "The regional header is missing. The DDL buffer is too"
			" short to hold another regional header.";
	case kNoRegionalScalars:
		return "The regional scalars are missing or corrupt. The DDL"
			" buffer is too short to contain them.";
	case kNoEndOfRegional:
		return "The DDL buffer is too short to contain an end of regional"
			" header key word.";
	case kBadEndOfRegional:
		return "End of regional header key word is incorrect or corrupt.";
	case kNoLocalStruct:
		return "The local structure is missing. The DDL buffer is too"
			" short to hold another local structure.";
	case kNoLocalScalars:
		return "The local scalars are missing or corrupt. The DDL buffer"
			" is too short to contain them.";
	case kNoEndOfLocal:
		return "The DDL buffer is too short to contain an end of local"
			" structure key word.";
	case kBadEndOfLocal:
		return "End of local structure key word is incorrect or corrupt.";
	case kBufferTooBig:
		return "The DDL raw data is larger than indicated by the headers;"
			" extra bytes are probably just garbage.";
	default:
		return "Unknown error code!";
	}
}


inline std::ostream& operator << (std::ostream& os, AliMUONTriggerDDLDecoderEventHandler::ErrorCode code)
{
	/// This is the stream operator for std::ostream classes to be able to
	/// easily write the error messages associated with the error codes generated
	/// by the decoder to 'cout' or 'cerr' for example.
	
	os << AliMUONTriggerDDLDecoderEventHandler::ErrorCodeToMessage(code);
	return os;
}

#endif // ALIMUONTRIGGERDDLDECODEREVENTHANDLER_H
