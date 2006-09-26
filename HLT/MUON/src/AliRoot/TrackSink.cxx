////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* AliHLTMUONTrackSink is used as the output target object by AliHLTMUONMicrodHLT.
   It is just a fancy buffer to store tracks that are found by the dHLT
   algorithm.
 */
 
#include "AliRoot/TrackSink.hpp"
#include "AliRoot/Base.hpp"

ClassImp(AliHLTMUONTrackSink)
ClassImp(AliHLTMUONTrackSink::AliEventData)


AliHLTMUONTrackSink::AliHLTMUONTrackSink() :
	TObject(),
	fFilename(""), fFoldername(""), fEventIndex(-1),
	fCurrentEvent(NULL), fBlockIndex(-1), fCurrentBlock(NULL),
	fTrackIndex(-1), fCurrentTrack(NULL),
	fEventList(AliHLTMUONTrackSink::AliEventData::Class())
{
// Default constructor

	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
}


AliHLTMUONTrackSink::~AliHLTMUONTrackSink()
{
// destructor.

	DebugMsg(1, "AliHLTMUONTrackSink::~AliHLTMUONTrackSink()");
	fEventList.Clear("C");
}


void AliHLTMUONTrackSink::AddEvent(Int_t eventnumber)
{
// Adds a new AliEventData block to the fEventList and updates internal pointers.
// Cannot have duplicate event numbers so this method will display an error
// message if one attempts to add the same event number more than once.

	DebugMsg(1, "AliHLTMUONTrackSink::AddEvent(" << eventnumber << ")");
	
	if (eventnumber < 0)
	{
		Error("AddEvent", "The event number must be positive, got: %d", eventnumber);
		return;
	};
	
	// Make sure that the event number is not already added to the list of events.
	TIter next(&fEventList);
	AliEventData* current;
	while ((current = (AliEventData*)next()))
	{
		if (current->EventNumber() == eventnumber)
		{
			Error("AddEvent", "The event number %d is already stored.", eventnumber);
			return;
		};
	};
	
	fEventIndex = fEventList.GetEntriesFast();
	new ( fEventList[fEventIndex] ) AliEventData(eventnumber); 
	fCurrentEvent = (AliEventData*) fEventList[fEventIndex];
	
	// Remember to reset the other pointers because the new event is empty.
	ResetBlockPointers();
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void AliHLTMUONTrackSink::AddBlock()
{
// Adds a new block to the current event and updates fCurrentBlock and
// fCurrentTrack.

	DebugMsg(1, "AliHLTMUONTrackSink::AddBlock()");
	
	if (fCurrentEvent == NULL)
	{
		Error("AddBlock", "No event selected.");
		return;
	};
	
	fBlockIndex = fCurrentEvent->Blocks().GetEntriesFast();
	new ( fCurrentEvent->Blocks()[fBlockIndex] ) TClonesArray(AliHLTMUONTrack::Class());
	fCurrentBlock = (TClonesArray*) fCurrentEvent->Blocks()[fBlockIndex];
	
	// Remember to reset the track pointer because the new block is empty.
	ResetTrackPointers();

	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


void AliHLTMUONTrackSink::AddTrack(const AliHLTMUONTrack& track)
{
// Adds the given track to the current block.
// If no current block is selected then an error message is displayed.

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
// Adds a new track to the current event and block, and returns a pointer
// to this track object to be filled by the caller.
// The fCurrentTrack is updated appropriately.
// If no current block is selected then NULL is returned.

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
// Adds the specified track parameters as a new track.
// The fCurrentTrack is updated appropriately.

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
// Sets the internal file and folder names from the trigger source.

	fFilename = triggersource->FileName();
	fFoldername = triggersource->FolderName();
}


void AliHLTMUONTrackSink::Clear(Option_t* /*option*/)
{
// Clears all the internal arrays.

	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fEventList.Clear("C");
}


Bool_t AliHLTMUONTrackSink::GetEvent(Int_t eventnumber) const
{
// Fetches the event data corresponding to the given event number.
// Sets the current block and track to the first block and track for 
// the newly selected event.
// If there are no blocks or tracks then these pointers are set to NULL.
// kTRUE is returned if the event was found. kFALSE is returned if the
// event was not found and the internal pointers left untouched.

	DebugMsg(1, "AliHLTMUONTrackSink::GetEvent(" << eventnumber << ")" );
	
	// Try find the corresponding event in the list of events.
	for (Int_t i = 0; i < fEventList.GetEntriesFast(); i++)
	{
		AliEventData* current = (AliEventData*) fEventList[i];
		if (current->EventNumber() == eventnumber)
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
// Fetches the first event stored in this AliHLTMUONTrackSink.
// Sets the current block and track to the first block and track of the
// first event.
// If there are no blocks or tracks then these pointers are set to NULL.
// kFALSE is returned if there are no events stored, kTRUE is returned
// on success.

	DebugMsg(1, "AliHLTMUONTrackSink::GetFirstEvent()");
	if (fEventList.GetEntriesFast() > 0)
	{
		fEventIndex = 0;
		fCurrentEvent = (AliEventData*) fEventList[0];
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
// Returns kTRUE if there are more events to iterate over, kFALSE otherwise.

	return 0 <= fEventIndex && fEventIndex < fEventList.GetEntriesFast();
}


Bool_t AliHLTMUONTrackSink::GetNextEvent() const
{
// Fetches the next event stored following the currently selected one.
// Sets the current block and track pointers to the first block and track
// in the newly selected event. These pointers are set to NULL if there
// are no blocks or tracks for this event.

	DebugMsg(1, "AliHLTMUONTrackSink::GetNextEvent()");
	if (fEventIndex < fEventList.GetEntriesFast() - 1)
	{
		fCurrentEvent = (AliEventData*) fEventList[ ++fEventIndex ];
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
// Returns the corresponding AliRoot event number for the current event.
// -1 is returned if no event is selected.

	if (fCurrentEvent != NULL)
		return fCurrentEvent->EventNumber();
	else
		return -1;
}


Int_t AliHLTMUONTrackSink::NumberOfBlocks() const
{
// Returns the number of track blocks in the current event.
// -1 is returned if no event is selected.

	DebugMsg(1, "AliHLTMUONTrackSink::NumberOfBlocks()");
	if (fCurrentEvent == NULL)
	{
		Error("NumberOfBlocks", "No event selected.");
		return -1;
	}
	else
		return fCurrentEvent->Blocks().GetEntriesFast();
}


Bool_t AliHLTMUONTrackSink::GetBlock(Int_t index) const
{
// Fetches the index'th block in the current event.
// Sets the current track to the first track in the block.
// If there are no tracks then the track pointers are reset.
// kTRUE is returned if the block was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONTrackSink::GetBlock(" << index << ")");
	
	// Note NumberOfBlocks() also checks if the event was selected.
	Int_t numberofblocks = NumberOfBlocks();
	if (numberofblocks < 0) return kFALSE;

	if ( 0 <= index && index < numberofblocks )
	{
		fBlockIndex = index;
		fCurrentBlock = (TClonesArray*) fCurrentEvent->Blocks()[index];
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
// Fetches the first block in the current event.
// Sets the current track to the first track in the block.
// If there are no tracks then the fCurrentTracks pointer is set to NULL.

	DebugMsg(1, "AliHLTMUONTrackSink::GetFirstBlock()");
	// Note: NumberOfBlocks() also checks if fCurrentEvent != NULL.
	if (NumberOfBlocks() > 0)
	{
		fBlockIndex = 0;
		fCurrentBlock = (TClonesArray*) fCurrentEvent->Blocks()[fBlockIndex];
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
// Returns kTRUE if there are more blocks to iterate over.

	return 0 <= fBlockIndex && fBlockIndex < NumberOfBlocks();
}


Bool_t AliHLTMUONTrackSink::GetNextBlock() const
{
// Fetches the next block in the current event.
// kTRUE is returned if the block was found, kFALSE otherwise.
// The current track pointers are reset if no more blocks are found.

	DebugMsg(1, "AliHLTMUONTrackSink::GetNextBlock()");

	// Note: NumberOfBlocks() checks if fCurrentEvent != NULL. If it is then it returns -1
	// and since fBlockIndex is always >= -1 the if statement must go to the else part.
	if (fBlockIndex < NumberOfBlocks() - 1)
	{
		fCurrentBlock = (TClonesArray*) fCurrentEvent->Blocks()[ ++fBlockIndex ];
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
// Returns the number of tracks in the current block.
// -1 is returned if no block is selected.

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
// Fetches the index'th track in the current block.
// NULL is returned if the track was not found.

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
// Fetches the first track in the current block.
// NULL is returned if the track was not found.

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
// Returns kTRUE if there are more tracks to iterate over in
// the current block. kFALSE is returned otherwise.

	return 0 <= fTrackIndex && fTrackIndex < NumberOfTracks();
}


const AliHLTMUONTrack* AliHLTMUONTrackSink::GetNextTrack() const
{
// Fetches the next track in the current block.
// NULL is returned if the track was not found.

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
// Sets all the current pointers to NULL and indices to -1.

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
// Sets the block and track pointers to NULL and indices to -1.

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
// Sets just the current track pointer to NULL and index to -1.

	fTrackIndex = -1;
	fCurrentTrack = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTrackIndex = " << fTrackIndex
	);
}


AliHLTMUONTrackSink::AliEventData::AliEventData()
	: fEventNumber(-1), fBlocks(TClonesArray::Class())
{
	fEventNumber = -1;
}


AliHLTMUONTrackSink::AliEventData::AliEventData(Int_t eventnumber)
	: fEventNumber(eventnumber), fBlocks(TClonesArray::Class())
{
// Creates a new event data block with the specified event number.

	fEventNumber = eventnumber;
	
	// If the following is not set then we do not write the fBlocks properly.
	fBlocks.BypassStreamer(kFALSE);
}


AliHLTMUONTrackSink::AliEventData::~AliEventData()
{
	DebugMsg(1, "AliHLTMUONTrackSink::AliEventData::~AliEventData()");
	fBlocks.Clear("C");
}

