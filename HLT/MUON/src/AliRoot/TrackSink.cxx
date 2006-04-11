////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/TrackSink.hpp"
#include "AliRoot/Base.hpp"

ClassImp(AliHLTMUONTrackSink)
ClassImp(AliHLTMUONTrackSink::EventData)


AliHLTMUONTrackSink::AliHLTMUONTrackSink() :
	TObject(), fEventList(AliHLTMUONTrackSink::EventData::Class())
{
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
}


AliHLTMUONTrackSink::~AliHLTMUONTrackSink()
{
	DebugMsg(1, "AliHLTMUONTrackSink::~AliHLTMUONTrackSink()");
	fEventList.Clear("C");
}


void AliHLTMUONTrackSink::AddEvent(Int_t eventnumber)
{
	DebugMsg(1, "AliHLTMUONTrackSink::AddEvent(" << eventnumber << ")");
	
	if (eventnumber < 0)
	{
		Error("AddEvent", "The event number must be positive, got: %d", eventnumber);
		return;
	};
	
	// Make sure that the event number is not already added to the list of events.
	TIter next(&fEventList);
	EventData* current;
	while ((current = (EventData*)next()))
	{
		if (current->fEventNumber == eventnumber)
		{
			Error("AddEvent", "The event number %d is already stored.", eventnumber);
			return;
		};
	};
	
	fEventIndex = fEventList.GetEntriesFast();
	new ( fEventList[fEventIndex] ) EventData(eventnumber); 
	fCurrentEvent = (EventData*) fEventList[fEventIndex];
	
	// Remember to reset the other pointers because the new event is empty.
	ResetBlockPointers();
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void AliHLTMUONTrackSink::AddBlock()
{
	DebugMsg(1, "AliHLTMUONTrackSink::AddBlock()");
	
	if (fCurrentEvent == NULL)
	{
		Error("AddBlock", "No event selected.");
		return;
	};
	
	fBlockIndex = fCurrentEvent->fBlocks.GetEntriesFast();
	new ( fCurrentEvent->fBlocks[fBlockIndex] ) TClonesArray(AliHLTMUONTrack::Class());
	fCurrentBlock = (TClonesArray*) fCurrentEvent->fBlocks[fBlockIndex];
	
	// Remember to reset the track pointer because the new block is empty.
	ResetTrackPointers();

	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void AliHLTMUONTrackSink::AddTrack(const AliHLTMUONTrack& track)
{
	DebugMsg(1, "AliHLTMUONTrackSink::AddTrack()");

	if (fCurrentBlock == NULL)
	{
		Error("AddTrack", "No block selected.");
		return;
	};
	
	fTrackIndex = fCurrentBlock->GetEntriesFast();
	new ( (*fCurrentBlock)[fTrackIndex] ) AliHLTMUONTrack(track);
	fCurrentTrack = (AliHLTMUONTrack*) (*fCurrentBlock)[fTrackIndex];
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


AliHLTMUONTrack* AliHLTMUONTrackSink::AddTrack()
{
	DebugMsg(1, "AliHLTMUONTrackSink::AddTrack()");

	if (fCurrentBlock == NULL)
	{
		Error("AddTrack", "No block selected.");
		return NULL;
	};
	
	fTrackIndex = fCurrentBlock->GetEntriesFast();
	fCurrentTrack = (AliHLTMUONTrack*) fCurrentBlock->New(fTrackIndex);
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
	return fCurrentTrack;
};


void AliHLTMUONTrackSink::AddTrack(
		Int_t triggerid, Int_t sign, Float_t momentum,
		Float_t pt, const AliHLTMUONPoint hits[10], const AliHLTMUONRegion regions[10]
	)
{
	DebugMsg(1, "AliHLTMUONTrackSink::AddTrack(" << triggerid << ", " << sign << ", "
		<< momentum << ", " << pt << ", " << (void*)(&hits[0]) << ", "
		<< (void*)(&regions[0]) << ")"
	);

	if (fCurrentBlock == NULL)
	{
		Error("AddTrack", "No block selected.");
		return;
	};
	
	fTrackIndex = fCurrentBlock->GetEntriesFast();
	new ( (*fCurrentBlock)[fTrackIndex] ) AliHLTMUONTrack(triggerid, sign, momentum, pt, hits, regions);
	fCurrentTrack = (AliHLTMUONTrack*) (*fCurrentBlock)[fTrackIndex];
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void AliHLTMUONTrackSink::SetNames(const AliHLTMUONTriggerSource* triggersource)
{
	fFilename = triggersource->FileName();
	fFoldername = triggersource->FolderName();
}


void AliHLTMUONTrackSink::Clear(Option_t* /*option*/)
{
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fEventList.Clear("C");
}


Bool_t AliHLTMUONTrackSink::GetEvent(Int_t eventnumber) const
{
	DebugMsg(1, "AliHLTMUONTrackSink::GetEvent(" << eventnumber << ")" );
	
	// Try find the corresponding event in the list of events.
	for (Int_t i = 0; i < fEventList.GetEntriesFast(); i++)
	{
		EventData* current = (EventData*) fEventList[i];
		if (current->fEventNumber == eventnumber)
		{
			fEventIndex = i;
			fCurrentEvent = current;
			GetFirstBlock();
			DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
				<< " , fTrackIndex = " << fTrackIndex
			);
			return kTRUE;
		}
	}
	return kFALSE;
}


Bool_t AliHLTMUONTrackSink::GetFirstEvent() const
{
	DebugMsg(1, "AliHLTMUONTrackSink::GetFirstEvent()");
	if (fEventList.GetEntriesFast() > 0)
	{
		fEventIndex = 0;
		fCurrentEvent = (EventData*) fEventList[0];
		GetFirstBlock();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return kTRUE;
	}
	else
	{
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return kFALSE;
	}
}


Bool_t AliHLTMUONTrackSink::MoreEvents() const
{
	return 0 <= fEventIndex && fEventIndex < fEventList.GetEntriesFast();
}


Bool_t AliHLTMUONTrackSink::GetNextEvent() const
{
	DebugMsg(1, "AliHLTMUONTrackSink::GetNextEvent()");
	if (fEventIndex < fEventList.GetEntriesFast() - 1)
	{
		fCurrentEvent = (EventData*) fEventList[ ++fEventIndex ];
		GetFirstBlock();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return kTRUE;
	}
	else
	{
		ResetAllPointers();
		return kFALSE;
	};
};


Int_t AliHLTMUONTrackSink::CurrentEvent() const
{
	if (fCurrentEvent != NULL)
		return fCurrentEvent->fEventNumber;
	else
		return -1;
}


Int_t AliHLTMUONTrackSink::NumberOfBlocks() const
{
	DebugMsg(1, "AliHLTMUONTrackSink::NumberOfBlocks()");
	if (fCurrentEvent == NULL)
	{
		Error("NumberOfBlocks", "No event selected.");
		return -1;
	}
	else
		return fCurrentEvent->fBlocks.GetEntriesFast();
}


Bool_t AliHLTMUONTrackSink::GetBlock(Int_t index) const
{
	DebugMsg(1, "AliHLTMUONTrackSink::GetBlock(" << index << ")");
	
	// Note NumberOfBlocks() also checks if the event was selected.
	Int_t numberofblocks = NumberOfBlocks();
	if (numberofblocks < 0) return kFALSE;

	if ( 0 <= index && index < numberofblocks )
	{
		fBlockIndex = index;
		fCurrentBlock = (TClonesArray*) fCurrentEvent->fBlocks[index];
		GetFirstTrack();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return kTRUE;
	}
	else
	{
		// The index is out of bounds so inform the user.
		if (numberofblocks > 0)
			Error(	"GetBlock",
				"The block index (%d) is out of bounds. Valid range is [0, %d]",
				index, numberofblocks - 1
			);
		else
			Error(	"GetBlock",
				"The block index (%d) is out of bounds. No blocks found.",
				index
			);
		return kFALSE;
	}
}


Bool_t AliHLTMUONTrackSink::GetFirstBlock() const
{
	DebugMsg(1, "AliHLTMUONTrackSink::GetFirstBlock()");
	// Note: NumberOfBlocks() also checks if fCurrentEvent != NULL.
	if (NumberOfBlocks() > 0)
	{
		fBlockIndex = 0;
		fCurrentBlock = (TClonesArray*) fCurrentEvent->fBlocks[fBlockIndex];
		GetFirstTrack();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return kTRUE;
	}
	else
		return kFALSE;
}


Bool_t AliHLTMUONTrackSink::MoreBlocks() const
{
	return 0 <= fBlockIndex && fBlockIndex < NumberOfBlocks();
}


Bool_t AliHLTMUONTrackSink::GetNextBlock() const
{
	DebugMsg(1, "AliHLTMUONTrackSink::GetNextBlock()");

	// Note: NumberOfBlocks() checks if fCurrentEvent != NULL. If it is then it returns -1
	// and since fBlockIndex is always >= -1 the if statement must go to the else part.
	if (fBlockIndex < NumberOfBlocks() - 1)
	{
		fCurrentBlock = (TClonesArray*) fCurrentEvent->fBlocks[ ++fBlockIndex ];
		GetFirstTrack();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return kTRUE;
	}
	else
	{
		ResetBlockPointers();
		return kFALSE;
	}
}


Int_t AliHLTMUONTrackSink::NumberOfTracks() const
{
	DebugMsg(1, "AliHLTMUONTrackSink::NumberOfTracks()");
	if (fCurrentBlock == NULL)
	{
		Error("NumberOfTracks", "No block selected.");
		return -1;
	}
	else
		return fCurrentBlock->GetEntriesFast();
};


const AliHLTMUONTrack* AliHLTMUONTrackSink::GetTrack(Int_t index) const
{
	DebugMsg(1, "AliHLTMUONTrackSink::GetTrack(" << index << ")");

	// Note NumberOfTracks() also checks if the event and block was selected.
	Int_t numberoftracks = NumberOfTracks();
	if (numberoftracks < 0) return NULL;
	
	if ( 0 <= index && index < numberoftracks )
	{
		fTrackIndex = index;
		fCurrentTrack = (AliHLTMUONTrack*) fCurrentBlock->At(index);
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return fCurrentTrack;
	}
	else
	{
		// The index is out of bounds so inform the user.
		if (numberoftracks > 0)
			Error(	"GetTrack",
				"The track index (%d) is out of bounds. Valid range is [0, %d]",
				index, numberoftracks - 1
			);
		else
			Error(	"GetTrack",
				"The track index (%d) is out of bounds. No tracks found.",
				index
			);
		return NULL;
	}
}


const AliHLTMUONTrack* AliHLTMUONTrackSink::GetFirstTrack() const
{
	DebugMsg(1, "AliHLTMUONTrackSink::GetFirstTrack()");
	// Note: NumberOfTracks() also checks if fCurrentBlock != NULL.
	if (NumberOfTracks() > 0)
	{
		fTrackIndex = 0;
		fCurrentTrack = (AliHLTMUONTrack*) fCurrentBlock->At(0);
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return fCurrentTrack;
	}
	else
		return NULL;
}


Bool_t AliHLTMUONTrackSink::MoreTracks() const
{
	return 0 <= fTrackIndex && fTrackIndex < NumberOfTracks();
}


const AliHLTMUONTrack* AliHLTMUONTrackSink::GetNextTrack() const
{
	DebugMsg(1, "AliHLTMUONTrackSink::GetNextTrack()");
	
	// Note: NumberOfTracks() checks if fCurrentBlock != NULL. If it is then it returns -1
	// and since fTrackIndex is always >= -1 the if statement must go to the else part.
	if (fTrackIndex < NumberOfTracks() - 1)
	{
		fCurrentTrack = (AliHLTMUONTrack*) fCurrentBlock->At( ++fTrackIndex );
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return fCurrentTrack;
	}
	else
	{
		ResetTrackPointers();
		return NULL;
	}
}


void AliHLTMUONTrackSink::ResetAllPointers() const
{
	fEventIndex = -1;
	fCurrentEvent = NULL;
	fBlockIndex = -1;
	fCurrentBlock = NULL;
	fTrackIndex = -1;
	fCurrentTrack = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void AliHLTMUONTrackSink::ResetBlockPointers() const
{
	fBlockIndex = -1;
	fCurrentBlock = NULL;
	fTrackIndex = -1;
	fCurrentTrack = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void AliHLTMUONTrackSink::ResetTrackPointers() const
{
	fTrackIndex = -1;
	fCurrentTrack = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


AliHLTMUONTrackSink::EventData::EventData() : fBlocks(TClonesArray::Class())
{
	fEventNumber = -1;
}


AliHLTMUONTrackSink::EventData::EventData(Int_t eventnumber)
	: fBlocks(TClonesArray::Class())
{
	fEventNumber = eventnumber;
	
	// If the following is not set then we do not write the fBlocks properly.
	fBlocks.BypassStreamer(kFALSE);
}


AliHLTMUONTrackSink::EventData::~EventData()
{
	DebugMsg(1, "AliHLTMUONTrackSink::EventData::~EventData()");
	fBlocks.Clear("C");
}

