////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* AliHLTMUONADCStreamSource is used to store the raw DDL data generated
   from an AliRoot simulation.
   This class is used as a storage class for the input dataset for
   AliHLTMUONMicrodHLT.
 */

#include "AliRoot/ADCStreamSource.hpp"
#include "AliRoot/Base.hpp"
#include "TSystem.h"
#include <stdio.h>

// TODO: Change all the Error message statements to AliError at some stage.

ClassImp(AliHLTMUONADCStreamSource)
ClassImp(AliHLTMUONADCStreamSource::AliDataBlock)


AliHLTMUONADCStreamSource::AliHLTMUONADCStreamSource() :
	TObject(), fCurrentStream(-1),
	fList(AliHLTMUONADCStreamSource::AliDataBlock::Class())
{
	fCurrentStream = -1;
}


AliHLTMUONADCStreamSource::~AliHLTMUONADCStreamSource()
{
	// Everything is cleaned up implicitly.
}


void AliHLTMUONADCStreamSource::FillFromFile(const TString& filename, Int_t eventnumber)
{
// Fills the internal data structures from the specified file.

	DebugMsg(1, "Entering FillFromFile, file = " << filename.Data()
		<< ", event number = " << eventnumber
	);
	FILE* file = fopen(filename, "r");
	try
	{
		Long_t id = -1;
		Long64_t size = -1;
		Long_t flags = -1;
		Long_t modtime = -1;
		if ( gSystem->GetPathInfo(filename, &id, &size, &flags, &modtime) == 0 )
		{
			DebugMsg(2, "Size of file: " << filename.Data() << " is " << size << " bytes");
			AliHLTMUONADCStream stream;
			stream.Size(size/sizeof(Int_t));
			size_t bytesread = fread(stream.Data(), 1, size, file);
			if (bytesread == (size_t)size)
				AddStream(stream, eventnumber);
			else
				Error("FillFromFile", "Could not read from file: %s", filename.Data());
		}
		else
			Error("FillFromFile", "Could not stat the file: %s", filename.Data());
	}
	finally
	(
		fclose(file);
	);
	
	DebugMsg(1, "Leaving FillFromFile");
}


void AliHLTMUONADCStreamSource::FillFrom(const TString& directory, Int_t eventnumber)
{
// Fills the internal data structures from the specified directory.
// FillFromFile is called for every file in the directory that is
// prefixed with MUON_ and ends in .ddl

	DebugMsg(1, "Entering FillFrom, directory = " << directory.Data()
		<< ", event number = " << eventnumber
	);
	
	void* dir = gSystem->OpenDirectory(directory);
	try
	{
		const char* entry = gSystem->GetDirEntry(dir);
		while (entry != NULL)
		{
			TString filename = entry;
			if (filename.BeginsWith("MUON"))
				FillFromFile(filename, eventnumber);
			entry = gSystem->GetDirEntry(dir);
		};
	}
	finally
	(
		gSystem->FreeDirectory(dir);
	);
	
	
	DebugMsg(1, "Leaving FillFrom");
}


void AliHLTMUONADCStreamSource::FillFrom(
		const TString& dirprefix, UInt_t firstevent, UInt_t lastevent
	)
{
// Same as the methods above except the directory name is created as
// dirprefix + eventnumber, where eventnumber is looped from firstevent
// to lastevent.

	DebugMsg(1, "Entering FillFrom");
	
	for (UInt_t i = firstevent; i <= lastevent; i++)
	{
		TString dirname = dirprefix;
		dirname += i;
		FillFrom(dirname, i);
	};
	
	DebugMsg(1, "Leaving FillFrom");
}


void AliHLTMUONADCStreamSource::Clear(Option_t* /*option*/)
{
// Clears all the internal arrays.

	fCurrentStream = -1;
	fList.Clear("C");
}


Int_t AliHLTMUONADCStreamSource::NumberOfStreams() const
{
// Returns the number of ADC streams stored.

	return fList.GetEntriesFast();
}


Bool_t AliHLTMUONADCStreamSource::GetStream(Int_t index) const
{
// Fetches the index'th ADC stream stored.
// kTRUE is returned if the stream was found, kFALSE otherwise.

	if ( 0 <= index && index < NumberOfStreams() )
	{
		fCurrentStream = index;
		return kTRUE;
	}
	else
	{
		// The index is out of bounds so inform the user.
		if (NumberOfStreams() > 0)
			Error(	"GetStream",
				"The ADC stream index (%d) is out of bounds. Valid range is [0, %d]",
				index, NumberOfStreams() - 1
			);
		else
			Error(	"GetStream",
				"The ADC stream index (%d) is out of bounds. No streams found.",
				index
			);
		return kFALSE;
	}
}


Bool_t AliHLTMUONADCStreamSource::FirstStream() const
{
// Fetches the first ADC stream stored.
// kTRUE is returned if the stream was found, kFALSE otherwise.

	if (NumberOfStreams() > 0)
	{
		fCurrentStream = 0;
		return kTRUE;
	}
	else
		return kFALSE;
}


Bool_t AliHLTMUONADCStreamSource::NextStream() const
{
// Fetches the next ADC stream stored following the currently selected one.
// kTRUE is returned if the stream was found, kFALSE otherwise.

	if ( 0 <= fCurrentStream && fCurrentStream < NumberOfStreams() - 1 )
	{
		fCurrentStream++;
		return kTRUE;
	}
	else
		return kFALSE;
};


Int_t AliHLTMUONADCStreamSource::EventNumber() const
{
// Returns the corresponding AliRoot event number for the current stream.
// -1 is returned if no stream is selected.

	if (fCurrentStream >= 0)
	{
		Assert( fCurrentStream < NumberOfStreams() );
		return ((AliDataBlock*)fList[fCurrentStream])->EventNumber();
	}
	else
	{
		Error("EventNumber", "No ADC stream selected.");
		return -1;
	}
}


Bool_t AliHLTMUONADCStreamSource::FetchStream(AliHLTMUONADCStream& stream) const
{
// Returns the current ADC stream selected.
// kFALSE is returned if there is no stream selected.

	if (fCurrentStream >= 0)
	{
		Assert( fCurrentStream < NumberOfStreams() );
		stream = ((AliDataBlock*)fList[fCurrentStream])->Stream();
		return kTRUE;
	}
	else
	{
		Error("FetchStream", "No ADC stream selected.");
		return kFALSE;
	}
}


Bool_t AliHLTMUONADCStreamSource::FetchStream(Int_t index, AliHLTMUONADCStream& stream) const
{
// Returns the index'th ADC stream.
// kTRUE is returned if the stream was found, kFALSE otherwise.

	if ( GetStream(index) )
		return FetchStream(stream);
	else
		return kFALSE;
}


const AliHLTMUONADCStream* AliHLTMUONADCStreamSource::FetchStream() const
{
// Returns the current ADC stream selected.
// A NULL pointer is returned if no ADC stream is selected.

	if (fCurrentStream >= 0)
	{
		Assert( fCurrentStream < NumberOfStreams() );
		return &( ((AliDataBlock*)fList[fCurrentStream])->Stream() );
	}
	else
	{
		Error("FetchStream", "No ADC stream selected.");
		return NULL;
	}
}


void AliHLTMUONADCStreamSource::AddStream(
		AliHLTMUONADCStream& stream, UInt_t eventnumber
	)
{
// Adds a new AliHLTMUONADCStream object to the internal arrays.

	DebugMsg(1, "Entering AddStream");
	
	fCurrentStream = fList.GetEntriesFast();
	new ( fList[fCurrentStream] ) AliDataBlock();
	((AliDataBlock*)fList[fCurrentStream])->EventNumber() = eventnumber;
	((AliDataBlock*)fList[fCurrentStream])->Stream() = stream;
	
	DebugMsg(1, "Leaving AddStream");
}

