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

/// \ingroup macros
/// \file MUONTimeRawStreamTrigger.C
/// \brief Macro for checking the timing (speed) performace of the two different trigger decoders.
///
/// \author Artur Szostak <artursz@iafrica.com>
///
/// This macro is used to check the timing (speed) performance of the existing
/// offline decoder for the trigger DDLs and also for the new high performance
/// decoder. It can be invoked as follows:
/// 
///  $ aliroot
/// .L $ALICE_ROOT/MUON/MUONTimeRawStreamTrigger.C+
///  MUONTimeRawStreamTrigger(filename, maxEvent);
///
/// where \em filename is the name of a file containing the raw data, or alternatively
/// the directory containing rawX (X being an integer) paths with the raw DDL
/// data. The \em maxEvent value is the maximum event to process (default set to
/// 1000). Thus the macro will time the algorithm for all events in the range
/// [0 .. maxEvent-1].
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliCodeTimer.h"

// MUON includes
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONRawStreamTriggerHP.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"

// RAW includes
#include "AliRawReader.h"
#include "AliRawReaderMemory.h"
#include "AliRawDataHeader.h"

#include "TStopwatch.h"
#include "TMath.h"
#include "Riostream.h"

#endif


// Linked list node for buffer structures.
struct AliBufferInfo
{
	AliBufferInfo* fNext;
	Int_t fEquipId;
	UInt_t fBufferSize;
	UChar_t* fBuffer;
};


// local structure information to store
struct AliLocalStructInfo
{
	UChar_t fId;
	UChar_t fDec;
	Bool_t fTrigY;
	UChar_t fYPos;
	UChar_t fSXDev;
	UChar_t fXDev;
	UChar_t fXPos;
	Bool_t fTriggerY;
	Bool_t fTriggerX;
};


UInt_t LoadFiles(AliBufferInfo*& list, TString fileName = "./", Int_t maxEvent = 1000)
{
	/// Reads in the DDL files into memory buffers as a linked list.

	AliRawReader* rawReader = AliRawReader::Create(fileName.Data()); 

	if (rawReader == NULL)
	{
		cerr << "ERROR: Could not create AliRawReader." << endl;
		delete rawReader;
		return 0;
	}

	UInt_t count = 0;
	Int_t iEvent = 0;
	list = NULL;
	while (rawReader->NextEvent())
	{
		if (iEvent++ >= maxEvent) break;

		rawReader->Select("MUONTRG", 0, 1);
		while (rawReader->ReadHeader())
		{
			AliBufferInfo* info = new AliBufferInfo;
			if (info == NULL)
			{
				cerr << "ERROR: Out of memory, sorry. You should limit the number of events read in." << endl;
				delete rawReader;
				return count;
			}
			info->fNext = list;
			info->fEquipId = rawReader->GetEquipmentId();
			info->fBufferSize = rawReader->GetDataSize() + sizeof(AliRawDataHeader);
			info->fBuffer = new UChar_t[info->fBufferSize];
			if (info->fBuffer == NULL)
			{
				cerr << "ERROR: Out of memory, sorry. You should limit the number of events read in." << endl;
				delete rawReader;
				return count;
			}
			list = info;
			
			// Copy the header.
			memcpy(info->fBuffer, rawReader->GetDataHeader(), sizeof(AliRawDataHeader));

			// Now copy the payload.
			if (! rawReader->ReadNext(
						info->fBuffer + sizeof(AliRawDataHeader),
						info->fBufferSize - sizeof(AliRawDataHeader)
					)
			   )
			{
				cerr << "ERROR: Failed to read from AliRawReader." << endl;
			}
			count++;
		}
	}

	delete rawReader;
	return count;
}


UInt_t CountMaxStructs(AliBufferInfo* list)
{
	/// Counts the maximum number of local structures possible

	UInt_t total = 0;
	AliBufferInfo* current = list;
	while (current != NULL)
	{
		total += current->fBufferSize / 5;
		current = current->fNext;
	}
	return total;
}


void ReleaseBuffers(AliBufferInfo* list)
{
	/// Deletes the memory allocated for the linked list of buffers.

	AliBufferInfo* current = list;
	while (current != NULL)
	{
		AliBufferInfo* tmp = current;
		current = current->fNext;
		delete [] tmp->fBuffer;
		delete tmp;
	}
}


Double_t TimeUsingOldDecoder(AliBufferInfo* list, AliLocalStructInfo* buffer, UInt_t maxBufferSize)
{
	/// Perform a timing using the old decoder.

	AliRawReaderMemory rawReader;
	AliMUONRawStreamTrigger rawStream(&rawReader);
	rawReader.NextEvent();

	TStopwatch timer;
	timer.Start(kTRUE);

	UInt_t i = 0;
	AliBufferInfo* current = list;
	while (current != NULL)
	{
		rawReader.SetMemory(current->fBuffer, current->fBufferSize);
		rawReader.SetEquipmentID(current->fEquipId);
		rawReader.Reset();

		rawStream.First();
		
		TArrayS fXPattern;
		TArrayS fYPattern;

		while ( rawStream.Next(buffer[i].fId, buffer[i].fDec, buffer[i].fTrigY,
		                       buffer[i].fYPos, buffer[i].fSXDev, buffer[i].fXDev,
		                       buffer[i].fXPos, buffer[i].fTriggerY, buffer[i].fTriggerX,
		                       fXPattern, fYPattern)
		      )
		{
			if (i < maxBufferSize-1) i++;
		}

		current = current->fNext;
	}

	return timer.RealTime();
}


Double_t TimeUsingNewDecoder(AliBufferInfo* list, AliLocalStructInfo* buffer, UInt_t maxBufferSize)
{
	/// Perform a timing using the new decoder.

	AliRawReaderMemory rawReader;
	AliMUONRawStreamTriggerHP rawStream(&rawReader);
	rawReader.NextEvent();

	TStopwatch timer;
	timer.Start(kTRUE);

	UInt_t i = 0;
	AliBufferInfo* current = list;
	while (current != NULL)
	{
		rawReader.SetMemory(current->fBuffer, current->fBufferSize);
		rawReader.SetEquipmentID(current->fEquipId);
		rawReader.Reset();

		rawStream.First();

		const AliMUONRawStreamTriggerHP::AliLocalStruct* local = NULL;
		while ( (local = rawStream.Next()) != NULL )
		{
			if (i < maxBufferSize)
			{
				buffer[i].fId = local->GetId();
				buffer[i].fDec = local->GetDec();
				buffer[i].fTrigY = local->GetTrigY();
				buffer[i].fYPos = local->GetYPos();
				buffer[i].fSXDev = local->GetSXDev();
				buffer[i].fXDev = local->GetXDev();
				buffer[i].fXPos = local->GetXPos();
				buffer[i].fTriggerY = local->GetTriggerY();
				buffer[i].fTriggerX = local->GetTriggerX();
				i++;
			}
		}

		current = current->fNext;
	}

	return timer.RealTime();
}


Double_t TimeUsingNewDecoderOldInterface(AliBufferInfo* list, AliLocalStructInfo* buffer, UInt_t maxBufferSize)
{
	/// Perform a timing using the new decoder but the old Next() method
	/// as the interface.

	AliRawReaderMemory rawReader;
	AliMUONRawStreamTriggerHP rawStream(&rawReader);
	rawReader.NextEvent();

	TStopwatch timer;
	timer.Start(kTRUE);

	UInt_t i = 0;
	AliBufferInfo* current = list;
	while (current != NULL)
	{
		rawReader.SetMemory(current->fBuffer, current->fBufferSize);
		rawReader.SetEquipmentID(current->fEquipId);
		rawReader.Reset();

		rawStream.First();
		
		TArrayS fXPattern;
		TArrayS fYPattern;

		while ( rawStream.Next(buffer[i].fId, buffer[i].fDec, buffer[i].fTrigY,
		                       buffer[i].fYPos, buffer[i].fSXDev, buffer[i].fXDev,
		                       buffer[i].fXPos, buffer[i].fTriggerY, buffer[i].fTriggerX,
		                       fXPattern, fYPattern)
		      )
		{
			if (i < maxBufferSize-1) i++;
		}

		current = current->fNext;
	}

	return timer.RealTime();
}


void Loop(const char* filename, Bool_t newDecoder)
{
  AliCodeTimerAutoGeneral(Form("Loop %s",(newDecoder ? "NEW":"OLD")),);
  
  AliRawReader* reader = AliRawReader::Create(filename);
  
  AliMUONVRawStreamTrigger* stream;
  
  if ( newDecoder ) 
  {
    stream = new AliMUONRawStreamTriggerHP(reader);
  }
  else
  {
    stream = new AliMUONRawStreamTrigger(reader);
  }

  UChar_t id;
  UChar_t dec;
  Bool_t trigY;
  UChar_t yPos;
  UChar_t sXDev;
  UChar_t xDev;
  UChar_t xPos;
  Bool_t triggerY;
  Bool_t triggerX;
  TArrayS xPattern;
  TArrayS yPattern;
  
  while ( reader->NextEvent() )
  {
    stream->First();
    
    while ( stream->Next(id, dec, trigY, yPos, sXDev, xDev, xPos,
		         triggerY, triggerX, xPattern, yPattern) )
    {
      id *= 2;
    }
  }
  
  delete stream;
  delete reader;
}


void MUONTimeRawStreamTriggerDumb(TString fileName)
{
  AliCodeTimer::Instance()->Reset();
  
  // first check we can open the stream
  AliRawReader* reader = AliRawReader::Create(fileName.Data());
  if (!reader)
  {
    cerr << "Cannot create reader from " << fileName.Data() << endl;
    return;
  }
  delete reader;
  
  // now start the timing per se
  
  Loop(fileName,kFALSE);
  
  Loop(fileName,kTRUE);
  
  AliCodeTimer::Instance()->Print();
}


void MUONTimeRawStreamTrigger(TString fileName = "./", Int_t maxEvent = 1000)
{
	/// Performs a timing of old and new decoders and reports this.

	AliBufferInfo* list = NULL;
	UInt_t bufferCount = LoadFiles(list, fileName, maxEvent);
	if (bufferCount == 0)
	{
		cerr << "ERROR: No DDL files found or read in." << endl;
		return;
	}

	UInt_t maxBufferSize = CountMaxStructs(list);
	AliLocalStructInfo* buffer = new AliLocalStructInfo[maxBufferSize];
	if (buffer == NULL)
	{
		ReleaseBuffers(list);
		cerr << "ERROR: Out of memory, sorry. You should limit the number of events read in." << endl;
		return;
	}
	Double_t oldTimes[100];
	for (Int_t i = 0; i < 100; i++)
	{
		cout << "Timing old decoder: timing iteration " << i+1 << " of 100" << endl;
		oldTimes[i] = TimeUsingOldDecoder(list, buffer, maxBufferSize);
	}
	Double_t newTimes[100];
	for (Int_t i = 0; i < 100; i++)
	{
		cout << "Timing new decoder: timing iteration " << i+1 << " of 100" << endl;
		newTimes[i] = TimeUsingNewDecoder(list, buffer, maxBufferSize);
	}
	Double_t newTimes2[100];
	for (Int_t i = 0; i < 100; i++)
	{
		cout << "Timing new decoder with old interface: timing iteration " << i+1 << " of 100" << endl;
		newTimes2[i] = TimeUsingNewDecoderOldInterface(list, buffer, maxBufferSize);
	}

	ReleaseBuffers(list);
	delete buffer;

	Double_t oldTime = TMath::Mean(100, oldTimes) / Double_t(bufferCount);
	Double_t oldTimeErr = TMath::RMS(100, oldTimes) / Double_t(bufferCount);
	Double_t newTime = TMath::Mean(100, newTimes) / Double_t(bufferCount);
	Double_t newTimeErr = TMath::RMS(100, newTimes) / Double_t(bufferCount);
	Double_t newTime2 = TMath::Mean(100, newTimes2) / Double_t(bufferCount);
	Double_t newTime2Err = TMath::RMS(100, newTimes2) / Double_t(bufferCount);

	cout << "Average processing time per DDL for:" << endl;
	cout << "                   Old decoder = " << oldTime*1e6 << " +/- " << oldTimeErr*1e6/TMath::Sqrt(100) << " micro seconds" << endl;
	cout << "                   New decoder = " << newTime*1e6 << " +/- " << newTimeErr*1e6/TMath::Sqrt(100) << " micro seconds" << endl;
	cout << "New decoder with old interface = " << newTime2*1e6 << " +/- " << newTime2Err*1e6/TMath::Sqrt(100) << " micro seconds" << endl;
}

