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

/// \class AliMUONRawStreamTriggerHP
///
/// Implementation of a streamer interface to the high performance trigger decoder.
/// This is the raw stream class which interfaces between the high performance
/// core decoder for MUON trigger chambers and the AliRawReader class.
/// To gain the most out of the decoder, the Next() method should be used,
/// for example:
/// \code
///   AliMUONRawStreamTriggerHP* rawStream;  // assume initialised
///   const AliMUONRawStreamTriggerHP::AliLocalStruct* localStruct;
///   while ((localStruct = rawStream->Next()) != NULL)
///   {
///      // Do something with localStruct here.
///   }
/// \endcode
///
/// This decoder tries to implement as similar an interface as possible to
/// AliMUONRawStreamTrigger where possible. However certain constructs which
/// would slow us down too much are avoided.
///
/// \author Artur Szostak <artursz@iafrica.com>

#include "AliMUONRawStreamTriggerHP.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include <cassert>
#include <iostream>
using std::cout;
using std::endl;
using std::hex;
using std::dec;

/// \cond CLASSIMP
ClassImp(AliMUONRawStreamTriggerHP)
/// \endcond

const Int_t AliMUONRawStreamTriggerHP::fgkMaxDDL = 2;
bool AliMUONRawStreamTriggerHP::AliLocalStruct::fgOverrideId = true;

const AliMUONRegionalHeaderStruct
AliMUONRawStreamTriggerHP::AliDecoderEventHandler::fgkEmptyHeader = {
	AliMUONTriggerDDLDecoder<AliDecoderEventHandler>::RegionalErrorWord(),
	0,
	{0, 0},
	0
};


AliMUONRawStreamTriggerHP::AliMUONRawStreamTriggerHP() :
	AliMUONVRawStreamTrigger(),
	fDecoder(),
	fDDL(0),
	fBufferSize(8192),
	fBuffer(new UChar_t[8192]),
	fkCurrentLocalStruct(NULL),
	fHadError(kFALSE),
	fDone(kFALSE),
	fDDLObject(NULL)
{
	///
	/// Default constructor.
	///
	
	fDecoder.ExitOnError(false);

	fDecoder.GetHandler().SetMaxStructs(
			fDecoder.MaxRegionals(),
			fDecoder.MaxLocals()
		);

	fDecoder.GetHandler().SetRawStream(this);
}


AliMUONRawStreamTriggerHP::AliMUONRawStreamTriggerHP(AliRawReader* rawReader) :
	AliMUONVRawStreamTrigger(rawReader),
	fDecoder(),
	fDDL(0),
	fBufferSize(8192),
	fBuffer(new UChar_t[8192]),
	fkCurrentLocalStruct(NULL),
	fHadError(kFALSE),
	fDone(kFALSE),
	fDDLObject(NULL)
{
	///
	/// Constructor with AliRawReader as argument.
	///
	
	fDecoder.ExitOnError(false);

	fDecoder.GetHandler().SetMaxStructs(
			fDecoder.MaxRegionals(),
			fDecoder.MaxLocals()
		);
	
	fDecoder.GetHandler().SetRawStream(this);
}


AliMUONRawStreamTriggerHP::~AliMUONRawStreamTriggerHP()
{
	///
	/// Default destructor which cleans up the memory allocated.
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


void AliMUONRawStreamTriggerHP::First()
{
	/// Initialise or reset the iterator.
	/// The first DDL will be found and decoded.
	
	assert( GetReader() != NULL );
	
	fDDL = 0;
	fDone = kFALSE;
	NextDDL();
}


Bool_t AliMUONRawStreamTriggerHP::NextDDL()
{
	/// Read in the next trigger DDL and decode the payload with the
	/// high performance decoder.
	/// \return kTRUE if the next DDL was successfully read and kFALSE
	///    otherwise.

	assert( GetReader() != NULL );
	
	// The temporary object if generated in GetDDLTracker, is now stale,
	// so delete it.
	if (fDDLObject != NULL)
	{
		delete fDDLObject;
		fDDLObject = NULL;
	}
	
	fkCurrentLocalStruct = NULL;
	
	while (fDDL < GetMaxDDL())
	{
		GetReader()->Reset();
		GetReader()->Select("MUONTRG", fDDL, fDDL);  // Select the DDL file to be read.
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
	
	// Check if this is a scalar event.
	bool scalerEvent = (GetReader()->GetDataHeader()->GetL1TriggerMessage() & 0x1) == 0x1;
	
	bool result = false;
	try
	{
		// Since we might allocate memory inside OnNewBuffer in the event
		// handler we need to trap any memory allocation exception to be robust.
		result = fDecoder.Decode(fBuffer, dataSize, scalerEvent);
		fHadError = (result == true ? kFALSE : kTRUE);
	}
	catch (const std::bad_alloc&)
	{
		AliError("Could not allocate more buffer space. Cannot decode DDL.");
		return kFALSE;
	}

	// Update the current local structure pointer.
	fkCurrentLocalStruct = fDecoder.GetHandler().FirstLocalStruct();

	fDDL++; // Remember to increment index to next DDL.
	return kTRUE;
}


Bool_t AliMUONRawStreamTriggerHP::IsDone() const
{
	/// Indicates whether the iteration is finished or not.
	/// \return kTRUE if we already read all the digits and kFALSE if not.
	
	return fDone;
}


Bool_t AliMUONRawStreamTriggerHP::Next(
		UChar_t& id,   UChar_t& dec,     Bool_t& trigY,
		UChar_t& yPos, UChar_t& sXDev,   UChar_t& xDev,
		UChar_t& xPos, Bool_t& triggerY, Bool_t& triggerX,
		TArrayS& xPattern, TArrayS& yPattern
	)
{
	/// Advance one step in the iteration. Returns kFALSE if finished.
	/// If kTRUE is returned then the output parameters are filled with
	/// the values found in the next local trigger circuit structure.
	
	const AliLocalStruct* localStruct = Next();
	if (localStruct == NULL) return kFALSE;
	
	id    = localStruct->GetId();
	dec   = localStruct->GetDec();
	trigY = localStruct->GetTrigY();
	yPos  = localStruct->GetYPos();
	sXDev = localStruct->GetSXDev();
	xDev  = localStruct->GetXDev();
	xPos  = localStruct->GetXPos();
	
	triggerX = localStruct->GetTriggerX();
	triggerY = localStruct->GetTriggerY();
	
	localStruct->GetXPattern(xPattern);
	localStruct->GetYPattern(yPattern);

	return kTRUE;
}


AliMUONDDLTrigger* AliMUONRawStreamTriggerHP::GetDDLTrigger() const
{
	/// Construct and return a pointer to the DDL payload object.
	/// \return Pointer to internally constructed AliMUONDDLTrigger object.
	///         The object is owned by this class and should not be deleted
	///         by the caller.
	///
	/// \note This method should not be used just to gain access to the DDL
	/// payload, unless there is a good reason to have the AliMUONDDLTrigger
	/// object. For example, if you want to modify the data and then save it
	/// to another DDL stream. Otherwise it can be an order of magnitude
	/// faster to access the DDL headers and data with the GetHeaders,
	/// GetRegionalHeader and GetLocalStruct methods for example.
	/// Refer to the MUONRawStreamTrigger.C macro to see how to use the fast
	/// decoder interface optimally.
	
	if (fDDLObject != NULL) return fDDLObject;
	
	fDDLObject = new AliMUONDDLTrigger;
	
	// Copy over all DARC, global headers and scalars.
	AliMUONDarcHeader* darcHeader = fDDLObject->GetDarcHeader();
	const AliHeader* hdr = GetHeaders();
	UInt_t word = hdr->GetDarcHeader();
	memcpy(darcHeader->GetHeader(), &word, sizeof(word));
	if (hdr->GetDarcScalars() != NULL)
	{
		memcpy(darcHeader->GetDarcScalers(), hdr->GetDarcScalars(), sizeof(AliMUONDarcScalarsStruct));
	}
	memcpy(darcHeader->GetGlobalInput(), hdr->GetGlobalHeader(), sizeof(AliMUONGlobalHeaderStruct));
	if (hdr->GetGlobalScalars() != NULL)
	{
		memcpy(darcHeader->GetGlobalScalers(), hdr->GetGlobalScalars(), sizeof(AliMUONGlobalScalarsStruct));
	}
	
	for (Int_t iReg = 0; iReg < (Int_t)GetRegionalHeaderCount(); iReg++)
	{
		AliMUONRegHeader regHeader;
		AliMUONLocalStruct localStruct;
		
		const AliRegionalHeader* rh = GetRegionalHeader(iReg);
		// Copy local structure and scalars and add everything into DDL object.
		memcpy(regHeader.GetHeader(), rh->GetHeader(), sizeof(AliMUONRegionalHeaderStruct));
		if (rh->GetScalars() != NULL)
		{
			memcpy(regHeader.GetScalers(), rh->GetScalars(), sizeof(AliMUONRegionalScalarsStruct));
		}
		fDDLObject->AddRegHeader(regHeader);
		
		const AliLocalStruct* lstruct = rh->GetFirstLocalStruct();
		while (lstruct != NULL)
		{
			// Copy local structure and scalars and add everything into DDL object.
			memcpy(localStruct.GetData(), lstruct->GetData(), sizeof(AliMUONLocalInfoStruct));
			if (lstruct->GetScalars() != NULL)
			{
				memcpy(localStruct.GetScalers(), lstruct->GetScalars(), sizeof(AliMUONLocalScalarsStruct));
			}
			if (AliMUONRawStreamTriggerHP::AliLocalStruct::GetOverrideIdFlag() == true)
			{
				// Since the override ID flag is set, we need to replace the
				// ID in the structure with the calculated one returned by GetId().
				AliMUONLocalInfoStruct* strptr = reinterpret_cast<AliMUONLocalInfoStruct*>( localStruct.GetData() );
				UInt_t word = strptr->fTriggerBits;
				word &= (0xF << 19);
				strptr->fTriggerBits = word | (lstruct->GetId() << 19);
			}
			fDDLObject->AddLocStruct(localStruct, iReg);
			lstruct = lstruct->Next();
		}
	}
	
	return fDDLObject;
}


void AliMUONRawStreamTriggerHP::SetMaxReg(Int_t reg)
{
	/// Set the maximum allowed number of regional cards in the DDL.
	
	fDecoder.MaxRegionals( (UInt_t) reg );
	
	fDecoder.GetHandler().SetMaxStructs(
			fDecoder.MaxRegionals(),
			fDecoder.MaxLocals()
		);
}


void AliMUONRawStreamTriggerHP::SetMaxLoc(Int_t loc)
{
	/// Sets the maximum number of local cards in the DDL.
	
	fDecoder.MaxLocals( (UInt_t) loc );
	
	fDecoder.GetHandler().SetMaxStructs(
			fDecoder.MaxRegionals(),
			fDecoder.MaxLocals()
		);
}

///////////////////////////////////////////////////////////////////////////////

void AliMUONRawStreamTriggerHP::AliHeader::Print() const
{
	/// Print DARC header, global header and global scalars to screen.
	
	cout << "===== DARC info =====" << endl;
	cout << "Header bits : 0x" << hex << fDarcHeader << dec << endl;
	if (fDarcScalars != NULL)
	{
		cout << "L0R :   " << fDarcScalars->fL0R << " (0x"
			<< hex << fDarcScalars->fL0R << dec << ")" << endl;
		cout << "L1P :   " << fDarcScalars->fL1P << " (0x"
			<< hex << fDarcScalars->fL1P << dec << ")" << endl;
		cout << "L1S :   " << fDarcScalars->fL1S << " (0x"
			<< hex << fDarcScalars->fL1S << dec << ")" << endl;
		cout << "L2A :   " << fDarcScalars->fL2A << " (0x"
			<< hex << fDarcScalars->fL2A << dec << ")" << endl;
		cout << "L2R :   " << fDarcScalars->fL2R << " (0x"
			<< hex << fDarcScalars->fL2R << dec << ")" << endl;
		cout << "Clock : " << fDarcScalars->fClk << " (0x"
			<< hex << fDarcScalars->fClk << dec << ")" << endl;
		cout << "Hold :  " << fDarcScalars->fHold << " (0x"
			<< hex << fDarcScalars->fHold << dec << ")" << endl;
		cout << "Spare : " << fDarcScalars->fSpare << " (0x"
			<< hex << fDarcScalars->fSpare << dec << ")" << endl;
	}
	else
	{
		cout << "Scalars == NULL" << endl;
	}
	
	cout << "===== Global info =====" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << "Input[" << i << "] : " << fGlobalHeader->fInput[i] << " (0x"
			<< hex << fGlobalHeader->fInput[i] << dec << ")" << endl;
	}
	cout << "Output :    " << fGlobalHeader->fOutput << " (0x"
		<< hex << fGlobalHeader->fOutput << dec << ")" << endl;
	if (fGlobalScalars != NULL)
	{
		cout << "L0 :         " << fGlobalScalars->fL0 << " (0x"
			<< hex << fGlobalScalars->fL0 << dec << ")" << endl;
		cout << "Clock :     " << fGlobalScalars->fClk << " (0x"
			<< hex << fGlobalScalars->fClk << dec << ")" << endl;
		for (int j = 0; j < 4; j++)
		{
			cout << "Scaler[" << j << "] : " << fGlobalScalars->fScaler[j] << " (0x"
				<< hex << fGlobalScalars->fScaler[j] << dec << ")" << endl;
		}
		cout << "Hold :      " << fGlobalScalars->fHold << " (0x"
			<< hex << fGlobalScalars->fHold << dec << ")" << endl;
		cout << "Spare :     " << fGlobalScalars->fSpare << " (0x"
			<< hex << fGlobalScalars->fSpare << dec << ")" << endl;
	}
	else
	{
		cout << "Scalars == NULL" << endl;
	}
}

void AliMUONRawStreamTriggerHP::AliRegionalHeader::Print() const
{
	/// Print the regional header and scalars to screen.
	
	cout << "===== Regional card info =====" << endl;
	cout << "DarcWord : " << fHeader->fDarcWord << " (0x"
		<< hex << fHeader->fDarcWord << dec << ")" << endl;
	cout << "Word :     " << fHeader->fWord << " (0x"
		<< hex << fHeader->fWord << dec << ")" << endl;
	cout << "Input[0] : " << fHeader->fInput[0] << " (0x"
		<< hex << fHeader->fInput[0] << dec << ")" << endl;
	cout << "Input[1] : " << fHeader->fInput[1] << " (0x"
		<< hex << fHeader->fInput[1] << dec << ")" << endl;
	cout << "L0/Mask :  " << fHeader->fL0CountAndMask << " (0x"
		<< hex << fHeader->fL0CountAndMask << dec << ")" << endl;
	if (fScalars != NULL)
	{
		cout << "Clock :     " << fScalars->fClk << " (0x"
			<< hex << fScalars->fClk << dec << ")" << endl;
		for (int i = 0; i < 8; i++)
		{
			cout << "Scaler[" << i << "] : " << fScalars->fScaler[i] << " (0x"
				<< hex << fScalars->fScaler[i] << dec << ")" << endl;
		}
		cout << "Hold :      " << fScalars->fHold << " (0x"
			<< hex << fScalars->fHold << dec << ")" << endl;
	}
	else
	{
		cout << "Scalars == NULL" << endl;
	}
}

void AliMUONRawStreamTriggerHP::AliLocalStruct::Print() const
{
	/// Print local trigger structure and scalars to screen.
	
	cout << "===== Local card info =====" << endl;
	cout << "X2X1 :         " << fLocalStruct->fX2X1 << " (0x"
		<< hex << fLocalStruct->fX2X1 << dec << ")" << endl;
	cout << "X4X3 :         " << fLocalStruct->fX4X3 << " (0x"
		<< hex << fLocalStruct->fX4X3 << dec << ")" << endl;
	cout << "Y2Y1 :         " << fLocalStruct->fY2Y1 << " (0x"
		<< hex << fLocalStruct->fY2Y1 << dec << ")" << endl;
	cout << "Y4Y3 :         " << fLocalStruct->fY4Y3 << " (0x"
		<< hex << fLocalStruct->fY4Y3 << dec << ")" << endl;
	cout << "Trigger bits : " << fLocalStruct->fTriggerBits << " (0x"
		<< hex << fLocalStruct->fTriggerBits << dec << ")" << endl;
	if (fScalars != NULL)
	{
		cout << "L0 :           " << fScalars->fL0 << " (0x"
			<< hex << fScalars->fL0 << dec << ")" << endl;
		cout << "Hold :         " << fScalars->fHold << " (0x"
			<< hex << fScalars->fHold << dec << ")" << endl;
		cout << "Clock :        " << fScalars->fClk << " (0x"
			<< hex << fScalars->fClk << dec << ")" << endl;
		cout << "LPtNTrig :     " << fScalars->fLPtNTrig << " (0x"
			<< hex << fScalars->fLPtNTrig << dec << ")" << endl;
		cout << "HPtNTrig :     " << fScalars->fHPtNTrig << " (0x"
			<< hex << fScalars->fHPtNTrig << dec << ")" << endl;
		cout << "LPtRTrig :     " << fScalars->fLPtRTrig << " (0x"
			<< hex << fScalars->fLPtRTrig << dec << ")" << endl;
		cout << "HPtRTrig :     " << fScalars->fHPtRTrig << " (0x"
			<< hex << fScalars->fHPtRTrig << dec << ")" << endl;
		cout << "LPtLTrig :     " << fScalars->fLPtLTrig << " (0x"
			<< hex << fScalars->fLPtLTrig << dec << ")" << endl;
		cout << "HPtLTrig :     " << fScalars->fHPtLTrig << " (0x"
			<< hex << fScalars->fHPtLTrig << dec << ")" << endl;
		cout << "LPtSTrig :     " << fScalars->fLPtSTrig << " (0x"
			<< hex << fScalars->fLPtSTrig << dec << ")" << endl;
		cout << "HPtSTrig :     " << fScalars->fHPtSTrig << " (0x"
			<< hex << fScalars->fHPtSTrig << dec << ")" << endl;
		for (int i = 0; i < 8*4; i++)
		{
			cout << "Scaler[" << i << "] :  " << fScalars->fScaler[i] << " (0x"
				<< hex << fScalars->fScaler[i] << dec << ")" << endl;
		}
		cout << "EOS :          " << fScalars->fEOS << " (0x"
			<< hex << fScalars->fEOS << dec << ")" << endl;
		cout << "Reset :        " << fScalars->fReset << " (0x"
			<< hex << fScalars->fReset << dec << ")" << endl;
	}
	else
	{
		cout << "Scalars == NULL" << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

AliMUONRawStreamTriggerHP::AliDecoderEventHandler::AliDecoderEventHandler() :
	fRawStream(NULL),
	fBufferStart(NULL),
	fDarcHeader(0),
	fDarcScalars(NULL),
	fHeaders(),
	fRegionalsCount(0),
	fRegionals(NULL),
	fLocals(NULL),
	fEndOfLocals(NULL),
	fCurrentRegional(NULL),
	fCurrentLocal(NULL),
	fDarcEoWErrors(0),
	fGlobalEoWErrors(0),
	fRegEoWErrors(0),
	fLocalEoWErrors(0),
	fWarnings(kTRUE)
{
	/// Default constructor
}


AliMUONRawStreamTriggerHP::AliDecoderEventHandler::~AliDecoderEventHandler()
{
	/// Default destructor cleans up the allocated memory.
	
	if (fRegionals != NULL) delete [] fRegionals;
	if (fLocals != NULL) delete [] fLocals;
}


void AliMUONRawStreamTriggerHP::AliDecoderEventHandler::SetMaxStructs(
		UInt_t maxRegionals, UInt_t maxLocals
	)
{
	/// Sets the maximum number of structures allowed.
	
	// Start by clearing the current arrays.
	if (fRegionals != NULL)
	{
		delete [] fRegionals;
		fRegionals = NULL;
	}
	if (fLocals != NULL)
	{
		delete [] fLocals;
		fLocals = NULL;
		fEndOfLocals = NULL;
	}
	fCurrentRegional = NULL;
	fCurrentLocal = NULL;
	
	// Allocate new memory.
	fRegionals = new AliRegionalHeader[maxRegionals];
	fLocals = new AliLocalStruct[maxRegionals*maxLocals];
	fEndOfLocals = fLocals;
	
	fRegionalsCount = maxRegionals;
}


void AliMUONRawStreamTriggerHP::AliDecoderEventHandler::OnNewBuffer(
		const void* buffer, UInt_t /*bufferSize*/
	)
{
	/// This is called by the high performance decoder when a new DDL payload
	/// is about to be decoded.

	assert( fRawStream != NULL );
	
	// remember the start of the buffer to be used in OnError.
	fBufferStart = buffer;

	// Reset error counters.
	fDarcEoWErrors = 0;
	fGlobalEoWErrors = 0;
	fRegEoWErrors = 0;
	fLocalEoWErrors = 0;
	
	// Reset the current local structure pointer which will be used to track
	// where we need to fill fLocals. We have to subtract one space because we
	// will increment the pointer the first time in the OnLocalStruct method.
	fCurrentLocal = fLocals-1;
	
	fCurrentRegional = NULL;
	
	// Reset and link up all the regional structures together.
	for (UInt_t i = 0; i+1 < fRegionalsCount; i++)
	{
		fRegionals[i] = AliRegionalHeader(fLocals, &fgkEmptyHeader, NULL);
		fRegionals[i].SetNext(&fRegionals[i+1]);
	}
	// Reset the last structure.
	fRegionals[fRegionalsCount-1] = AliRegionalHeader(fLocals, &fgkEmptyHeader, NULL);
}


void AliMUONRawStreamTriggerHP::AliDecoderEventHandler::OnError(
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

	assert( fRawStream != NULL );
	assert( fRawStream->GetReader() != NULL );
	
	Char_t* message = NULL;
	//UInt_t word = 0;

	switch (error)
	{
	case kWrongEventType:
		message = Form("Wrong event type obtained from the Darc header, take the one of CDH");
		break;
		
	case kBadEndOfDarc:
		fDarcEoWErrors++;
		message = Form(
			"Wrong end of Darc word %x instead of %x\n",
			*reinterpret_cast<const UInt_t*>(location),
			AliMUONTriggerDDLDecoder<AliMUONTriggerDDLDecoderEventHandler>::EndOfDarcWord()
		);
		fRawStream->GetReader()->AddMajorErrorLog(kDarcEoWErr, message);
		break;
		
	case kBadEndOfGlobal:
		fGlobalEoWErrors++;
		message = Form(
			"Wrong end of Global word %x instead of %x\n",
			*reinterpret_cast<const UInt_t*>(location),
			AliMUONTriggerDDLDecoder<AliMUONTriggerDDLDecoderEventHandler>::EndOfGlobalWord()
		);
		fRawStream->GetReader()->AddMajorErrorLog(kGlobalEoWErr, message);
		break;
		
	case kBadEndOfRegional:
		fRegEoWErrors++;
		message = Form(
			"Wrong end of Regional word %x instead of %x\n",
			*reinterpret_cast<const UInt_t*>(location),
			AliMUONTriggerDDLDecoder<AliMUONTriggerDDLDecoderEventHandler>::EndOfRegionalWord()
		);
		fRawStream->GetReader()->AddMajorErrorLog(kRegEoWErr, message);
		break;
		
	case kBadEndOfLocal:
		fLocalEoWErrors++;
		message = Form(
			"Wrong end of Local word %x instead of %x\n",
			*reinterpret_cast<const UInt_t*>(location),
			AliMUONTriggerDDLDecoder<AliMUONTriggerDDLDecoderEventHandler>::EndOfLocalWord()
		);
		fRawStream->GetReader()->AddMajorErrorLog(kLocalEoWErr, message);
		break;
		
	default:
		message = Form(
			"%s (At byte %d in DDL.)",
			ErrorCodeToMessage(error),
			(unsigned long)location - (unsigned long)fBufferStart + sizeof(AliRawDataHeader)
		);
		fRawStream->GetReader()->AddMajorErrorLog(error, message);
		break;
	}

	if (fWarnings)
	{
		AliWarningGeneral(
				"AliMUONRawStreamTriggerHP::AliDecoderEventHandler",
				message
			);
	}
}

