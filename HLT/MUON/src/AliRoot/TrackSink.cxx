////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/TrackSink.hpp"
#include "AliRoot/Base.hpp"

ClassImp(AliMUONHLT::TrackSink)
ClassImp(AliMUONHLT::TrackSink::EventData)

namespace AliMUONHLT
{


TrackSink::TrackSink() :
	TObject(), fEventList(TrackSink::EventData::Class())
{
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
}


TrackSink::~TrackSink()
{
	DebugMsg(1, "TrackSink::~TrackSink()");
	fEventList.Clear("C");
}


void TrackSink::AddEvent(Int_t eventnumber)
{
	DebugMsg(1, "TrackSink::AddEvent(" << eventnumber << ")");
	
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


void TrackSink::AddBlock()
{
	DebugMsg(1, "TrackSink::AddBlock()");
	
	if (fCurrentEvent == NULL)
	{
		Error("AddBlock", "No event selected.");
		return;
	};
	
	fBlockIndex = fCurrentEvent->fBlocks.GetEntriesFast();
	new ( fCurrentEvent->fBlocks[fBlockIndex] ) TClonesArray(Track::Class());
	fCurrentBlock = (TClonesArray*) fCurrentEvent->fBlocks[fBlockIndex];
	
	// Remember to reset the track pointer because the new block is empty.
	ResetTrackPointers();

	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void TrackSink::AddTrack(const Track& track)
{
	DebugMsg(1, "TrackSink::AddTrack()");

	if (fCurrentBlock == NULL)
	{
		Error("AddTrack", "No block selected.");
		return;
	};
	
	fTrackIndex = fCurrentBlock->GetEntriesFast();
	new ( (*fCurrentBlock)[fTrackIndex] ) Track(track);
	fCurrentTrack = (Track*) (*fCurrentBlock)[fTrackIndex];
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


Track* TrackSink::AddTrack()
{
	DebugMsg(1, "TrackSink::AddTrack()");

	if (fCurrentBlock == NULL)
	{
		Error("AddTrack", "No block selected.");
		return NULL;
	};
	
	fTrackIndex = fCurrentBlock->GetEntriesFast();
	fCurrentTrack = (Track*) fCurrentBlock->New(fTrackIndex);
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
	return fCurrentTrack;
};


void TrackSink::AddTrack(
		Int_t triggerid, Int_t sign, Float_t momentum,
		Float_t pt, const Point hits[10], const Region regions[10]
	)
{
	DebugMsg(1, "TrackSink::AddTrack(" << triggerid << ", " << sign << ", "
		<< momentum << ", " << pt << ", " << (void*)(&hits[0]) << ", "
		<< (void*)(&regions[0]) << ")"
	);

	if (fCurrentBlock == NULL)
	{
		Error("AddTrack", "No block selected.");
		return;
	};
	
	fTrackIndex = fCurrentBlock->GetEntriesFast();
	new ( (*fCurrentBlock)[fTrackIndex] ) Track(triggerid, sign, momentum, pt, hits, regions);
	fCurrentTrack = (Track*) (*fCurrentBlock)[fTrackIndex];
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void TrackSink::SetNames(const TriggerSource* triggersource)
{
	fFilename = triggersource->FileName();
	fFoldername = triggersource->FolderName();
}


void TrackSink::Clear(Option_t* /*option*/)
{
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fEventList.Clear("C");
}


Bool_t TrackSink::GetEvent(Int_t eventnumber) const
{
	DebugMsg(1, "TrackSink::GetEvent(" << eventnumber << ")" );
	
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


Bool_t TrackSink::GetFirstEvent() const
{
	DebugMsg(1, "TrackSink::GetFirstEvent()");
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


Bool_t TrackSink::MoreEvents() const
{
	return 0 <= fEventIndex and fEventIndex < fEventList.GetEntriesFast();
}


Bool_t TrackSink::GetNextEvent() const
{
	DebugMsg(1, "TrackSink::GetNextEvent()");
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


Int_t TrackSink::CurrentEvent() const
{
	if (fCurrentEvent != NULL)
		return fCurrentEvent->fEventNumber;
	else
		return -1;
}


Int_t TrackSink::NumberOfBlocks() const
{
	DebugMsg(1, "TrackSink::NumberOfBlocks()");
	if (fCurrentEvent == NULL)
	{
		Error("NumberOfBlocks", "No event selected.");
		return -1;
	}
	else
		return fCurrentEvent->fBlocks.GetEntriesFast();
}


Bool_t TrackSink::GetBlock(Int_t index) const
{
	DebugMsg(1, "TrackSink::GetBlock(" << index << ")");
	
	// Note NumberOfBlocks() also checks if the event was selected.
	Int_t numberofblocks = NumberOfBlocks();
	if (numberofblocks < 0) return kFALSE;

	if ( 0 <= index and index < numberofblocks )
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


Bool_t TrackSink::GetFirstBlock() const
{
	DebugMsg(1, "TrackSink::GetFirstBlock()");
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


Bool_t TrackSink::MoreBlocks() const
{
	return 0 <= fBlockIndex and fBlockIndex < NumberOfBlocks();
}


Bool_t TrackSink::GetNextBlock() const
{
	DebugMsg(1, "TrackSink::GetNextBlock()");

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


Int_t TrackSink::NumberOfTracks() const
{
	DebugMsg(1, "TrackSink::NumberOfTracks()");
	if (fCurrentBlock == NULL)
	{
		Error("NumberOfTracks", "No block selected.");
		return -1;
	}
	else
		return fCurrentBlock->GetEntriesFast();
};


const Track* TrackSink::GetTrack(Int_t index) const
{
	DebugMsg(1, "TrackSink::GetTrack(" << index << ")");

	// Note NumberOfTracks() also checks if the event and block was selected.
	Int_t numberoftracks = NumberOfTracks();
	if (numberoftracks < 0) return NULL;
	
	if ( 0 <= index and index < numberoftracks )
	{
		fTrackIndex = index;
		fCurrentTrack = (Track*) fCurrentBlock->At(index);
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


const Track* TrackSink::GetFirstTrack() const
{
	DebugMsg(1, "TrackSink::GetFirstTrack()");
	// Note: NumberOfTracks() also checks if fCurrentBlock != NULL.
	if (NumberOfTracks() > 0)
	{
		fTrackIndex = 0;
		fCurrentTrack = (Track*) fCurrentBlock->At(0);
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTrackIndex = " << fTrackIndex
		);
		return fCurrentTrack;
	}
	else
		return NULL;
}


Bool_t TrackSink::MoreTracks() const
{
	return 0 <= fTrackIndex and fTrackIndex < NumberOfTracks();
}


const Track* TrackSink::GetNextTrack() const
{
	DebugMsg(1, "TrackSink::GetNextTrack()");
	
	// Note: NumberOfTracks() checks if fCurrentBlock != NULL. If it is then it returns -1
	// and since fTrackIndex is always >= -1 the if statement must go to the else part.
	if (fTrackIndex < NumberOfTracks() - 1)
	{
		fCurrentTrack = (Track*) fCurrentBlock->At( ++fTrackIndex );
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


void TrackSink::ResetAllPointers() const
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


void TrackSink::ResetBlockPointers() const
{
	fBlockIndex = -1;
	fCurrentBlock = NULL;
	fTrackIndex = -1;
	fCurrentTrack = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void TrackSink::ResetTrackPointers() const
{
	fTrackIndex = -1;
	fCurrentTrack = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


TrackSink::EventData::EventData() : fBlocks(TClonesArray::Class())
{
	fEventNumber = -1;
}


TrackSink::EventData::EventData(Int_t eventnumber)
	: fBlocks(TClonesArray::Class())
{
	fEventNumber = eventnumber;
	
	// If the following is not set then we do not write the fBlocks properly.
	fBlocks.BypassStreamer(kFALSE);
}


TrackSink::EventData::~EventData()
{
	DebugMsg(1, "TrackSink::EventData::~EventData()");
	fBlocks.Clear("C");
}


} // AliMUONHLT
