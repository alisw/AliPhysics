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
/// \file MUONTimeRawStreamTracker.C
/// \brief Macro for checking the timing (speed) performance of the tracker decoder.
///
/// \author Artur Szostak <artursz@iafrica.com>
///
/// This macro is used to check the timing (speed) performance of the 
/// decoder for the tracker DDLs. It can be invoked as follows:
/// 
///  $ aliroot
/// .L $ALICE_ROOT/MUON/MUONTimeRawStreamTracker.C+
///  MUONTimeRawStreamTracker(filename, maxEvent);
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
#include "AliMUONRawStreamTrackerHP.h"
#include "AliMUONDspHeader.h"
#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONDDLTracker.h"

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


// Digit information to store
struct AliDigitInfo
{
	UShort_t fManuId;
	UShort_t fAdc;
	UChar_t fChannelId;
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

		rawReader->Select("MUONTRK", 0, 19);
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


UInt_t CountMaxDigits(AliBufferInfo* list)
{
	/// Counts the maximum number of digits possible in all the buffers.

	UInt_t total = 0;
	AliBufferInfo* current = list;
	while (current != NULL)
	{
		total += current->fBufferSize;
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


Double_t TimeDecoderBusPatchIteration(AliBufferInfo* list, AliDigitInfo* buffer, UInt_t maxBufferSize)
{
	/// Perform a timing using the new decoder using the "next bus patch" iteration

	AliRawReaderMemory rawReader;
	AliMUONRawStreamTrackerHP rawStream(&rawReader);
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

		UShort_t manuId, adc;
		UChar_t manuChannel;

		rawStream.First();

		const AliMUONRawStreamTrackerHP::AliBusPatch* buspatch = NULL;
		while ((buspatch = rawStream.Next()) != NULL)
		{
			for (UInt_t j = 0; j < buspatch->GetDataCount(); j++)
			{
				buspatch->GetData(j, manuId, manuChannel, adc);
				if (i < maxBufferSize)
				{
					buffer[i].fManuId = manuId;
					buffer[i].fAdc = adc;
					buffer[i].fChannelId = manuChannel;
					i++;
				}
			}
		}

		current = current->fNext;
	}

	return timer.RealTime();
}


Double_t TimeDecoderChannelIteration(AliBufferInfo* list, AliDigitInfo* buffer, UInt_t maxBufferSize, Bool_t skipParityErrors)
{
	/// Perform a timing using the "next channel" iteration
  
	AliRawReaderMemory rawReader;
	AliMUONRawStreamTrackerHP rawStream(&rawReader);
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

		Int_t busPatch;
		UShort_t manuId, adc;
		UChar_t manuChannel;

		rawStream.First();

		while ( rawStream.Next(busPatch,manuId,manuChannel,adc,skipParityErrors) )
		{
			if (i < maxBufferSize)
			{
				buffer[i].fManuId = manuId;
				buffer[i].fAdc = adc;
				buffer[i].fChannelId = manuChannel;
				i++;
			}
		}

		current = current->fNext;
	}

	return timer.RealTime();
}

void MUONTimeRawStreamTrackerDumb(TString fileName)
{
  AliCodeTimer::Instance()->Reset();
  
  // first check we can open the stream
  AliRawReader* reader = AliRawReader::Create(fileName.Data());
  if (!reader)
  {
    cerr << "Cannot create reader from " << fileName.Data() << endl;
    return;
  }
  
  AliMUONRawStreamTrackerHP stream(reader);
  
  Int_t busPatch;
  UShort_t manuId, adc;
  UChar_t manuChannel;
  
  while ( reader->NextEvent() ) 
  {
    stream.First();
    
    while ( stream.Next(busPatch,manuId,manuChannel,adc) ) 
    {
      adc *= 2;
    }
  }
  
  AliCodeTimer::Instance()->Print();
}


void MUONTimeRawStreamTracker(TString fileName = "./", Int_t maxEvent = 1000)
{
	/// Performs a timing of decoder

	AliBufferInfo* list = NULL;
	UInt_t bufferCount = LoadFiles(list, fileName, maxEvent);
	if (bufferCount == 0)
	{
		cerr << "ERROR: No DDL files found or read in." << endl;
		return;
	}

	UInt_t maxBufferSize = CountMaxDigits(list);
	AliDigitInfo* buffer = new AliDigitInfo[maxBufferSize];
	if (buffer == NULL)
	{
		ReleaseBuffers(list);
		cerr << "ERROR: Out of memory, sorry. You should limit the number of events read in." << endl;
		return;
	}
	Double_t bpTimes[100];
	for (Int_t i = 0; i < 100; i++)
	{
		cout << "Timing decoder: bus patch iteration " << i+1 << " of 100" << endl;
		bpTimes[i] = TimeDecoderBusPatchIteration(list, buffer, maxBufferSize);
	}
	Double_t channelTimes[100];
	for (Int_t i = 0; i < 100; i++)
	{
		cout << "Timing decoder: channel iteration w/ parity check" << i+1 << " of 100" << endl;
		channelTimes[i] = TimeDecoderChannelIteration(list, buffer, maxBufferSize,kTRUE);
	}
	Double_t channelTimes2[100];
	for (Int_t i = 0; i < 100; i++)
	{
		cout << "Timing decoder: channel iteration w/o parity check" << i+1 << " of 100" << endl;
		channelTimes2[i] = TimeDecoderChannelIteration(list, buffer, maxBufferSize,kFALSE);
	}
  
	ReleaseBuffers(list);
	delete buffer;

	Double_t bpTime = TMath::Mean(100, bpTimes) / Double_t(bufferCount);
	Double_t bpTimeErr = TMath::RMS(100, bpTimes) / Double_t(bufferCount);
	Double_t channelTime = TMath::Mean(100, channelTimes) / Double_t(bufferCount);
	Double_t channelTimeErr = TMath::RMS(100, channelTimes) / Double_t(bufferCount);
	Double_t channelTime2 = TMath::Mean(100, channelTimes2) / Double_t(bufferCount);
	Double_t channelTime2Err = TMath::RMS(100, channelTimes2) / Double_t(bufferCount);

	cout << "Average processing time per DDL for:" << endl;
	cout << "   bus patch iteration                    = " << bpTime*1e6 << " +/- " << bpTimeErr*1e6/TMath::Sqrt(100) << " micro seconds" << endl;
	cout << "   channel iteration with parity check    = " << channelTime*1e6 << " +/- " << channelTimeErr*1e6/TMath::Sqrt(100) << " micro seconds" << endl;
	cout << "   channel iteration without parity check = " << channelTime2*1e6 << " +/- " << channelTime2Err*1e6/TMath::Sqrt(100) << " micro seconds" << endl;
}

