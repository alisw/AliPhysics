////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* AliHLTMUONTriggerSource is used to extract L0 trigger information for
   the muon spectrometer from a simulated event stored in .root files by AliRoot.
   It is used by the AliHLTMUONMicrodHLT class as a input data set for the
   dHLT algorithm.
 */
 
#include "AliRoot/TriggerSource.hpp"
#include "AliRoot/Base.hpp"
#include "Tracking/Calculations.hpp"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliModule.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONHit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONDataInterface.h"
#include "TDatabasePDG.h"
#ifndef __alpha
#include <math.h>
#else
#include <float.h>
#endif

ClassImp(AliHLTMUONTriggerSource)
ClassImp(AliHLTMUONTriggerSource::AliEventData)


AliHLTMUONTriggerSource::AliHLTMUONTriggerSource()
	: TObject(),
	  fAreaToUse(kFromWholePlane), fDataToUse(kFromLocalTriggers),
	  fMaxBlockSize(0xFFFFFFFF), fUseLookupTable(kTRUE),
	  fFilename(""), fFoldername(""),
	  fEventIndex(-1), fCurrentEvent(NULL),
	  fBlockIndex(-1), fCurrentBlock(NULL),
	  fTriggerIndex(-1), fCurrentTrigger(NULL),
	  fEventList(AliHLTMUONTriggerSource::AliEventData::Class()),
	  fHadToLoadgAlice(kFALSE)
{
// Default constructor.

	fAreaToUse = kFromWholePlane;
	fDataToUse = kFromLocalTriggers;
	fMaxBlockSize = 0xFFFFFFFF;
	fUseLookupTable = kTRUE;
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fHadToLoadgAlice = kFALSE;
}


AliHLTMUONTriggerSource::AliHLTMUONTriggerSource(AliMUONDataInterface* data)
	: TObject(),
	  fAreaToUse(kFromWholePlane), fDataToUse(kFromLocalTriggers),
	  fMaxBlockSize(0xFFFFFFFF), fUseLookupTable(kTRUE),
	  fFilename(""), fFoldername(""),
	  fEventIndex(-1), fCurrentEvent(NULL),
	  fBlockIndex(-1), fCurrentBlock(NULL),
	  fTriggerIndex(-1), fCurrentTrigger(NULL),
	  fEventList(AliHLTMUONTriggerSource::AliEventData::Class()),
	  fHadToLoadgAlice(kFALSE)
{
// Creates a new trigger source object by filling data from the data interface.

	fAreaToUse = kFromWholePlane;
	fDataToUse = kFromLocalTriggers;
	fMaxBlockSize = 0xFFFFFFFF;
	fUseLookupTable = kTRUE;
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fHadToLoadgAlice = kFALSE;
	FillFrom(data);
}


AliHLTMUONTriggerSource::~AliHLTMUONTriggerSource()
{
	fEventList.Delete();
}


void AliHLTMUONTriggerSource::FillFrom(AliMUONDataInterface* data)
{
// Fills the internal data structures from the specified data interface
// for all the events found in AliMUONDataInterface.
	   
	DebugMsg(1, "FillFrom(AliMUONDataInterface*)");
	
	if (FileAndFolderOk(data))
	{
		AliMUON* module = NULL;
		if ( ! FetchAliMUON(module) ) return;
		
		for (Int_t i = 0; i < data->NumberOfEvents(); i++)
		{
			AddEventFrom(data, module, i);
		}
		
		FinishedWithAliMUON();
	}
}


void AliHLTMUONTriggerSource::FillFrom(AliMUONDataInterface* data, Int_t event)
{
// Fills the internal data structures from the specified data interface
// for the given event.

	DebugMsg(1, "FillFrom(AliMUONDataInterface*, Int_t)");
	
	if (FileAndFolderOk(data))
	{
		AliMUON* module = NULL;
		if ( ! FetchAliMUON(module) ) return;
		AddEventFrom(data, module, event);
		FinishedWithAliMUON();
	}
}


void AliHLTMUONTriggerSource::FillFrom(
		AliMUONDataInterface* data,
		Int_t event, Int_t trigger, Bool_t newblock
	)
{
// Fills the internal data structures from the specified data interface
// for the given event and trigger number.
// If 'newblock' is set to true then the new trigger record is added to 
// a new block. Otherwise the point is added to the current block.
// Note: This method ignores the fAreaToUse and fMaxBlockSize flags.
// This is very usefull for custom trigger source filling.
// For the case of adding data from AliMUONHit objects the 'trigger'
// parameter becomes the track number in TreeH and not the index of the
// AliMUONLocalTrigger object.

	DebugMsg(1, "FillFrom(AliMUONDataInterface*, Int_t, Int_t, Bool_t)");
	
	if (FileAndFolderOk(data))
	{
		data->GetEvent(event);
		AliMUON* module = NULL;
		if ( ! FetchAliMUON(module) ) return;

		// Check if the current event corresponds to the event number we are
		// attempting to add to. If they do not or no event is selected then
		// try find the event or create a new one.
		if ( fCurrentEvent == NULL )
		{
			Bool_t found = GetEvent(event);
			if ( ! found) AddEvent(event);
		}
		else
		{
			if (fCurrentEvent->EventNumber() != event)
			{
				Bool_t found = GetEvent(event);
				if ( ! found) AddEvent(event);
			}
		}
		
		if ( fCurrentBlock != NULL )
		{
			Assert( fCurrentEvent != NULL );
			// If the newblock flag is set then force a new block.
			if (newblock) AddBlock();
		}
		else
			AddBlock();  // No block selected so we need to create a new block.

		AddTriggerFrom(data, module, trigger);
		FinishedWithAliMUON();
	}
}


void AliHLTMUONTriggerSource::Clear(Option_t* /*option*/)
{
// Clears all the internal arrays.

	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fEventList.Clear("C");
}


Bool_t AliHLTMUONTriggerSource::GetEvent(Int_t eventnumber) const
{
// Fetches the specified event number stored in this AliHLTMUONTriggerSource.
// Sets the current block and trigger to the first block and trigger record in
// the event. If there are no blocks or trigger records then these pointers are
// set to NULL. kTRUE is returned if the event was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONTriggerSource::GetEvent(" << eventnumber << ")" );
	
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
				<< " , fTriggerIndex = " << fTriggerIndex
			);
			return kTRUE;
		}
	}
	return kFALSE;
}


Bool_t AliHLTMUONTriggerSource::GetFirstEvent() const
{
// Fetches the first event stored in this AliHLTMUONTriggerSource.
// Sets the current block and trigger record to the first block and trigger
// in the event. If there are no blocks or trigger records then these pointers
// are set to NULL. kTRUE is returned if the event was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONTriggerSource::GetFirstEvent()");
	if (fEventList.GetEntriesFast() > 0)
	{
		fEventIndex = 0;
		fCurrentEvent = (AliEventData*) fEventList[0];
		GetFirstBlock();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kTRUE;
	}
	else
	{
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kFALSE;
	}
}


Bool_t AliHLTMUONTriggerSource::MoreEvents() const
{
// Returns kTRUE if there are more events to iterate over.

	return 0 <= fEventIndex && fEventIndex < fEventList.GetEntriesFast();
}


Bool_t AliHLTMUONTriggerSource::GetNextEvent() const
{
// Fetches the next event stored following the currently selected one.
// kTRUE is returned if the event was found, kFALSE otherwise.
// The internal pointers are reset if we reached the last event.
	   
	DebugMsg(1, "AliHLTMUONTriggerSource::GetNextEvent()");
	if (fEventIndex < fEventList.GetEntriesFast() - 1)
	{
		fCurrentEvent = (AliEventData*) fEventList[ ++fEventIndex ];
		GetFirstBlock();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kTRUE;
	}
	else
	{
		ResetAllPointers();
		return kFALSE;
	}
}


Int_t AliHLTMUONTriggerSource::CurrentEvent() const
{
// Returns the corresponding AliRoot event number for the current event.
// -1 is returned if no event is selected.
	   
	if (fCurrentEvent != NULL)
		return fCurrentEvent->EventNumber();
	else
		return -1;
}


Int_t AliHLTMUONTriggerSource::NumberOfBlocks() const
{
// Returns the number of trigger record blocks in the current event.
// -1 is returned if no event is selected.
	   
	DebugMsg(1, "AliHLTMUONTriggerSource::NumberOfBlocks()");
	if (fCurrentEvent == NULL)
	{
		Error("NumberOfBlocks", "No event selected.");
		return -1;
	}
	else
		return fCurrentEvent->Blocks().GetEntriesFast();
}


Bool_t AliHLTMUONTriggerSource::GetBlock(Int_t index) const
{
// Fetches the index'th block in the current event.
// Sets the current trigger record to the first trigger in the block.
// If there are no trigger records then this pointer is set to NULL.
// kTRUE is returned if the block was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONTriggerSource::GetBlock(" << index << ")");
	
	// Note NumberOfBlocks() also checks if the event was selected.
	Int_t numberofblocks = NumberOfBlocks();
	if (numberofblocks < 0) return kFALSE;

	if ( 0 <= index && index < numberofblocks )
	{
		fBlockIndex = index;
		fCurrentBlock = (TClonesArray*) fCurrentEvent->Blocks()[index];
		GetFirstTrigger();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
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


Bool_t AliHLTMUONTriggerSource::GetFirstBlock() const
{
// Fetches the first block in the current event.
// Sets the current trigger record to the first trigger in the block.
// If there are no trigger records then this pointer is set to NULL.
// kTRUE is returned if the block was found, kFALSE otherwise.
	   
	DebugMsg(1, "AliHLTMUONTriggerSource::GetFirstBlock()");
	// Note: NumberOfBlocks() also checks if fCurrentEvent != NULL.
	if (NumberOfBlocks() > 0)
	{
		fBlockIndex = 0;
		fCurrentBlock = (TClonesArray*) fCurrentEvent->Blocks()[fBlockIndex];
		GetFirstTrigger();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kTRUE;
	}
	else
		return kFALSE;
}


Bool_t AliHLTMUONTriggerSource::MoreBlocks() const
{
// Returns kTRUE if there are more blocks to be traversed.

	return 0 <= fBlockIndex && fBlockIndex < NumberOfBlocks();
}


Bool_t AliHLTMUONTriggerSource::GetNextBlock() const
{
// Fetches the next block in the current event.
// kTRUE is returned if the block was found, kFALSE otherwise.
	   
	DebugMsg(1, "AliHLTMUONTriggerSource::GetNextBlock()");

	// Note: NumberOfBlocks() checks if fCurrentEvent != NULL. If it is then it returns -1
	// and since fBlockIndex is always >= -1 the if statement must go to the else part.
	if (fBlockIndex < NumberOfBlocks() - 1)
	{
		fCurrentBlock = (TClonesArray*) fCurrentEvent->Blocks()[ ++fBlockIndex ];
		GetFirstTrigger();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return kTRUE;
	}
	else
	{
		ResetBlockPointers();
		return kFALSE;
	}
}


Int_t AliHLTMUONTriggerSource::NumberOfTriggers() const
{
// Returns the number of trigger records in the current block.
// -1 is returned if no block is selected.

	DebugMsg(1, "AliHLTMUONTriggerSource::NumberOfTriggers()");
	if (fCurrentBlock == NULL)
	{
		Error("NumberOfTriggers", "No block selected.");
		return -1;
	}
	else
		return fCurrentBlock->GetEntriesFast();
}


const AliHLTMUONTriggerRecord* AliHLTMUONTriggerSource::GetTrigger(Int_t triggernumber) const
{
// Fetches the trigger record with the specified trigger number from
// the current block.
// NULL is returned if the record was not found.

	DebugMsg(1, "AliHLTMUONTriggerSource::GetTrigger(" << triggernumber << ")");

	if (fCurrentBlock == NULL)
	{
		Error("GetTrigger", "No block selected.");
		return NULL;
	}
	
	// Try find the corresponding trigger record in the list of events.
	for (Int_t i = 0; i < fCurrentBlock->GetEntriesFast(); i++)
	{
		AliHLTMUONTriggerRecord* current = (AliHLTMUONTriggerRecord*) fCurrentBlock->At(i);
		if (current->TriggerNumber() == triggernumber)
		{
			fTriggerIndex = i;
			fCurrentTrigger = current;
			DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
				<< " , fTriggerIndex = " << fTriggerIndex
			);
			return current;
		}
	}
	return NULL;
}


const AliHLTMUONTriggerRecord* AliHLTMUONTriggerSource::GetFirstTrigger() const
{
// Fetches the first trigger record in the current block.
// NULL is returned if the record was not found.

	DebugMsg(1, "AliHLTMUONTriggerSource::GetFirstTrigger()");
	// Note: NumberOfTriggers() also checks if fCurrentBlock != NULL.
	if (NumberOfTriggers() > 0)
	{
		fTriggerIndex = 0;
		fCurrentTrigger = (AliHLTMUONTriggerRecord*) fCurrentBlock->At(0);
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return fCurrentTrigger;
	}
	else
		return NULL;
}


Bool_t AliHLTMUONTriggerSource::MoreTriggers() const
{
// Returns kTRUE if there are more triggers to iterate over.

	return 0 <= fTriggerIndex && fTriggerIndex < NumberOfTriggers();
}


const AliHLTMUONTriggerRecord* AliHLTMUONTriggerSource::GetNextTrigger() const
{
// Fetches the next trigger record in the current block.
// NULL is returned if the record was not found.

	DebugMsg(1, "AliHLTMUONTriggerSource::GetNextTrigger()");
	
	// Note: NumberOfTriggers() checks if fCurrentBlock != NULL. If it is then it returns -1
	// and since fTriggerIndex is always >= -1 the if statement must go to the else part.
	if (fTriggerIndex < NumberOfTriggers() - 1)
	{
		fCurrentTrigger = (AliHLTMUONTriggerRecord*) fCurrentBlock->At( ++fTriggerIndex );
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fTriggerIndex = " << fTriggerIndex
		);
		return fCurrentTrigger;
	}
	else
	{
		ResetTriggerPointers();
		return NULL;
	}
}


Int_t AliHLTMUONTriggerSource::CurrentTrigger() const
{
// Returns the trigger record number for the currently selected trigger record.
// This number corresponds to the index'th AliMUONLocalTrigger object for the
// current event.
// -1 is returned if no trigger record is selected.

	if (fCurrentTrigger != NULL)
	{
		return fCurrentTrigger->TriggerNumber();
	}
	else
	{
		Error("CurrentTrigger", "No trigger record selected.");
		return -1;
	}
}


void AliHLTMUONTriggerSource::AddEvent(Int_t eventnumber)
{
// Adds a new AliEventData block to the fEventList and updates the fCurrentEvent,
// fCurrentBlock and fCurrentTrigger pointers.

	DebugMsg(1, "AliHLTMUONTriggerSource::AddEvent(" << eventnumber << ")");
	Assert( eventnumber >= 0 );

	// Assume the eventnumber does not already exist in the event list.
	fEventIndex = fEventList.GetEntriesFast();
	new ( fEventList[fEventIndex] ) AliEventData(eventnumber);
	fCurrentEvent = (AliEventData*) fEventList[fEventIndex];
	
	// Remember to reset the other pointers because the new event is empty.
	ResetBlockPointers();
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


void AliHLTMUONTriggerSource::AddBlock()
{
// Adds a new block to the current event and updates fCurrentBlock and fCurrentTrigger.

	DebugMsg(1, "AliHLTMUONTriggerSource::AddBlock()");
	
	if (fCurrentEvent == NULL)
	{
		Error("AddBlock", "No event selected.");
		return;
	}
	
	fBlockIndex = fCurrentEvent->Blocks().GetEntriesFast();
	new ( fCurrentEvent->Blocks()[fBlockIndex] ) TClonesArray(AliHLTMUONTriggerRecord::Class());
	fCurrentBlock = (TClonesArray*) fCurrentEvent->Blocks()[fBlockIndex];
	
	// Remember to reset the trigger pointer because the new block is empty.
	ResetTriggerPointers();

	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


void AliHLTMUONTriggerSource::AddTrigger(const AliHLTMUONTriggerRecord& data)
{
// Adds a new trigger record to the current event and block.
// The fCurrentTrigger is updated appropriately.

	DebugMsg(1, "AliHLTMUONTriggerSource::AddTrigger(" << (void*)&data << ")");

	if (fCurrentBlock == NULL)
	{
		Error("AddTrigger", "No block selected.");
		return;
	}
	
	fTriggerIndex = fCurrentBlock->GetEntriesFast();
	new ( (*fCurrentBlock)[fTriggerIndex] ) AliHLTMUONTriggerRecord(data);
	fCurrentTrigger = (AliHLTMUONTriggerRecord*) (*fCurrentBlock)[fTriggerIndex];
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


Bool_t AliHLTMUONTriggerSource::FileAndFolderOk(AliMUONDataInterface* data)
{
// Checks if the file and folder names correspond to this AliHLTMUONTriggerSource's 
// file and folder names. kTRUE is returned if they do.
// If the file and folder names are empty then they are assigned the names
// as found in the data interface and kTRUE is returned.

	if (fFilename == "")
	{
		// Nothing filled yet so set the file and folder names.
		fFilename = data->CurrentFile();
		fFoldername = data->CurrentFolder();
		return kTRUE;
	}

	if ( fFilename != data->CurrentFile() )
	{
		Error(	"FileAndFolderOk",
			"The Trigger source already contains data from file '%s', cannot add data from file '%s'",
			fFilename.Data(), data->CurrentFile().Data()
		);
		return kFALSE;
	}
	
	if ( fFoldername != data->CurrentFolder() )
	{
		Error(	"FileAndFolderOk",
			"The Trigger source already contains data from folder '%s', cannot add data from folder '%s'",
			fFoldername.Data(), data->CurrentFolder().Data()
		);
		return kFALSE;
	}
	
	return kTRUE;
}


void AliHLTMUONTriggerSource::AddEventFrom(AliMUONDataInterface* data, AliMUON* module, Int_t event)
{
// Adds the whole event from the data interface to the internal data structures.
// It is assumed that FileAndFolderOk(data) returns true just before calling
// this method.

	if ( data->GetEvent(event) )
	{
		AddEvent(event);

		AddBlock();
		UInt_t currentblocksize = 0;
		AliHLTMUONTriggerRecord trigdata;

		switch (fDataToUse)
		{
		case kFromHits:
			for (Int_t track = 0; track < data->NumberOfTracks(); track++)
			{
				if ( ! FillTriggerFromHits(data, track, trigdata) )
					continue;  // Continue if unable to find hits.

				if (InFillRegion(trigdata))
				{
					AddTrigger(trigdata);

					// Create a new block if we reached the maximum block size.
					if ( ++currentblocksize == fMaxBlockSize )
					{
						AddBlock();
						currentblocksize = 0;
					}
				}
			}
			break;

		case kFromLocalTriggers:
			Assert( module != NULL );
			DebugMsg(4, "Taking kFromLocalTriggers branch...");
			for (Int_t i = 0; i < data->NumberOfLocalTriggers(); i++)
			{
				DebugMsg(4, "for loop: i = " << i);
				AliMUONLocalTrigger* lt = data->LocalTrigger(i);
				FillTriggerFromLocalTrigger(lt, module, trigdata);
				trigdata.TriggerNumber(i);

				if (InFillRegion(trigdata))
				{
					AddTrigger(trigdata);

					// Create a new block if we reached the maximum block size.
					if ( ++currentblocksize == fMaxBlockSize )
					{
						AddBlock();
						currentblocksize = 0;
					}
				}
			}
			break;

		default:
			Error("AddChamberFrom", "fDataToUse is not set to a valid value.");
		}
	}  // Loop on events.
}


void AliHLTMUONTriggerSource::AddTriggerFrom(AliMUONDataInterface* data, AliMUON* module, Int_t trigger)
{
// Adds the specified trigger record from the given data interface.
// The data interface should be initialised correctly, that is the event
// should already be selected before calling this method.

	DebugMsg(1, "Entering AddTriggerFrom");

	AliHLTMUONTriggerRecord trigdata;

	switch (fDataToUse)
	{
	case kFromHits:
		{
		// Note: in this case we treat the trigger parameter as a track number.
		if ( ! FillTriggerFromHits(data, trigger, trigdata) )
			return;  // Continue if unable to find hits.
		}
		break;

	case kFromLocalTriggers:
		{
		Assert( module != NULL );
		AliMUONLocalTrigger* lt = data->LocalTrigger(trigger);
		FillTriggerFromLocalTrigger(lt, module, trigdata);
		trigdata.TriggerNumber(trigger);
		}
		break;

	default:
		Error("AddTriggerFrom", "fDataToUse is not set to a valid value.");
		return;
	}
	
	AddTrigger(trigdata);
	
	DebugMsg(1, "Leaving AddTriggerFrom");
}


Bool_t AliHLTMUONTriggerSource::InFillRegion(const AliHLTMUONTriggerRecord& data) const
{
// Checks to see if the specified trigger record is in the chamber region
// we want to fill from.
// kTRUE is returned if (x, y) is in the region, and kFALSE otherwise.

	switch (fAreaToUse)
	{
	case kFromWholePlane:     return kTRUE;
	case kFromLeftHalfPlane:  return data.Station1Point().X() <= 0;
	case kFromRightHalfPlane: return data.Station1Point().X() > 0;

	default:
		Error("InFillRegion", "fAreaToUse is not set to a valid value.");
		return kFALSE;
	}
}


void AliHLTMUONTriggerSource::FillTriggerFromLocalTrigger(
		AliMUONLocalTrigger* trigger, AliMUON* module, AliHLTMUONTriggerRecord& record
	)
{
// Fills the trigger data from the AliMUONLocalTrigger object.
// if the fUseLookupTable is set to true then we use the L0 lookup table to
// fill the Pt value otherwise we use the PtCal method in AliMUONTriggerCircuit.
// Note the fTriggerNumber parameter is not filled in to 'record'.

	DebugMsg(2, "Creating TriggerRecord from AliMUONLocalTrigger object: " << (void*)trigger );
	AliMUONTriggerCircuit& circuit = module->TriggerCircuit(trigger->LoCircuit());

	// Get the sign of the particle the sign of the muon.
	if (trigger->LoLpt() == 1 || trigger->LoHpt() == 1 || trigger->LoApt() == 1)
	{
		record.ParticleSign(-1);
	}
	else
	if (trigger->LoLpt() == 2 || trigger->LoHpt() == 2 || trigger->LoApt() == 2)
	{
		record.ParticleSign(+1);
	}
	else
	{
		record.ParticleSign(0);
	}
	DebugMsg(2, "Particle sign = " << record.ParticleSign() );

	// Compute the transverse momentum.
	if (fUseLookupTable)
	{
		// TODO: implement use of the L0 lookup table.
		Error("FillTriggerFromLocalTrigger", "Use of L0 lookup table is not yet implemented!");
	}
	else
	{
		Float_t pt = circuit.PtCal( trigger->LoStripX(), trigger->LoDev(), trigger->LoStripY() );
		record.Pt(pt);
	}
	DebugMsg(2, "Pt = " << record.Pt() );

	// Build the impact points.
	record.Station1Point().X() = circuit.GetX11Pos(trigger->LoStripY());
	record.Station1Point().Y() = circuit.GetY11Pos(trigger->LoStripX());
	record.Station2Point().Y() = circuit.GetY21Pos(trigger->LoStripX() + trigger->LoDev() + 1);  // Why + 1?
	record.Station2Point().X() = AliMUONConstants::DefaultChamberZ(12) * record.Station1Point().X() / AliMUONConstants::DefaultChamberZ(10);
	DebugMsg(2, "fStation1x = " << record.Station1Point().X());
	DebugMsg(2, "fStation1y = " << record.Station1Point().Y());
	DebugMsg(2, "fStation2x = " << record.Station2Point().X());
	DebugMsg(2, "fStation2y = " << record.Station2Point().Y());
}


Bool_t AliHLTMUONTriggerSource::FillTriggerFromHits(
		AliMUONDataInterface* data, Int_t track, AliHLTMUONTriggerRecord& record
	)
{
// Fills the TriggerRecord structure from AliMUONHit objects.
// The hits on the last 4 chambers are used (i.e. chambers 11 to 14).
// kTRUE is returned if the structure was filled successfully.

	DebugMsg(2, "Creating TriggerRecord from hits on track: " << track );
	
	Float_t x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
#ifndef __alpha
#ifndef __sun
	x1 = y1 = z1 = x2 = y2 = z2 = x3 = y3 = z3 = x4 = y4 = z4 = nanf("");
#else
	x1 = y1 = z1 = x2 = y2 = z2 = x3 = y3 = z3 = x4 = y4 = z4 = 0;
#endif
#else
	x1 = y1 = z1 = x2 = y2 = z2 = x3 = y3 = z3 = x4 = y4 = z4 = FLT_QNAN;
#endif
	// Find the hit that corresponds to chambers. 11 to 14. We can ignore any
	// hits above the first 14. If there are that many it means the particle
	// is cycling in the detector.
	for (Int_t i = 0; i < data->NumberOfHits(track) && i < 14; i++)
	{
		AliMUONHit* h = data->Hit(track, i);

		// Note AliMUONHit::Chamber() returns a value in the range 1..14
		// It is also important to have positive Z coordinates and under the
		// new version of Aliroot we use GEANT coordinates which return 
		// negatives. So we use the fabs routine.
		switch ( h->Chamber() )
		{
		case 11: x1 = h->X(); y1 = h->Y(); z1 = (Float_t)fabs(h->Z()); break;
		case 12: x2 = h->X(); y2 = h->Y(); z2 = (Float_t)fabs(h->Z()); break;
		case 13: x3 = h->X(); y3 = h->Y(); z3 = (Float_t)fabs(h->Z()); break;
		case 14: x4 = h->X(); y4 = h->Y(); z4 = (Float_t)fabs(h->Z()); break;
		}
	}
	DebugMsg(4, "Found: x1 = " << x1 << ", y1 = " << y1 << ", z1 = " << z1);
	DebugMsg(4, "Found: x2 = " << x2 << ", y2 = " << y2 << ", z2 = " << z2);
	DebugMsg(4, "Found: x3 = " << x3 << ", y3 = " << y3 << ", z3 = " << z3);
	DebugMsg(4, "Found: x4 = " << x4 << ", y4 = " << y4 << ", z4 = " << z4);

	// Get a coordinate for station 1, perferably from chamber 11 otherwise
	// use hits from chamber 12.
	if ( ! TMath::IsNaN(x1))
	{
		record.Station1Point().X() = x1;
		record.Station1Point().Y() = y1;
		DebugMsg(3, "Using value from chamber 11: x1 = " << x1 << ", y1 = " << y1 << ", z1 = " << z1 );
	}
	else if ( ! TMath::IsNaN(x2))
	{
		record.Station1Point().X() = x2;
		record.Station1Point().Y() = y2;
		z1 = z2;
		DebugMsg(3, "Using value from chamber 12: x2 = " << x2 << ", y2 = " << y2 << ", z2 = " << z2 );
	}
	else
	{
		// Return false if we could not find any hits on chambers 11 or 12.
		Warning("FillTriggerFromHits", "Could not find any hits on chambers 11 and 12.");
		return kFALSE;
	}

	// Get a coordinate for station 2, perferably from chamber 13 otherwise
	// use hits from chamber 14.
	if ( ! TMath::IsNaN(x3))
	{
		record.Station2Point().X() = x3;
		record.Station2Point().Y() = y3;
		z2 = z3;
		DebugMsg(3, "Using value from chamber 13: x3 = " << x3 << ", y3 = " << y3 << ", z3 = " << z3 );
	}
	else if ( ! TMath::IsNaN(x4))
	{
		record.Station2Point().X() = x4;
		record.Station2Point().Y() = y4;
		z2 = z4;
		DebugMsg(3, "Using value from chamber 14: x4 = " << x4 << ", y4 = " << y4 << ", z4 = " << z4 );
	}
	else
	{
		// Return false if we could not find any hits on chambers 13 or 14.
		Warning("FillTriggerFromHits", "Could not find any hits on chambers 13 and 14.");
		return kFALSE;
	}
	
	record.TriggerNumber(track);
	
	// Get the sign of the particle.
	Int_t particlecode = (Int_t) data->Hit(track, 0)->Particle();
	DebugMsg(3, "particle code = " << particlecode);
	TDatabasePDG* pdb = TDatabasePDG::Instance();
	TParticlePDG* pdata = pdb->GetParticle(particlecode);
	if (pdata->Charge() < 0)
		record.ParticleSign(-1);
	else if (pdata->Charge() > 0)
		record.ParticleSign(+1);
	else
		record.ParticleSign(0);
	DebugMsg(3, "Particle sign = " << record.ParticleSign());
	
	DebugMsg(3, "Calculating Pt: x1 = " << record.Station1Point().X()
			<< ", y1 = " << record.Station1Point().Y()
			<< ", y2 = " << record.Station2Point().Y()
			<< ", z1 = " << z1
			<< ", z2 = " << z2
		);
	// Calculate and assign the transverse momentum.
	Float_t pt = AliHLTMUONCoreCalculatePt(
				record.Station1Point().X(),
				record.Station1Point().Y(), record.Station2Point().Y(),
				z1, z2
			);
	record.Pt(pt);
	
	DebugMsg(3, "Pt = " << record.Pt());
	
	return kTRUE;
}


Bool_t AliHLTMUONTriggerSource::FetchAliMUON(AliMUON*& module)
{
// Fetches the AliMUON module from the AliRun global object. AliRun will be loaded
// by the runloader if it has not yet been loaded. In such a case the AliRun object
// will also we unloaded when we are done with it.
// kTRUE is returned if no error occured and kFALSE otherwise.
// Note that if fDataToUse is set to kFromHits then gAlice is not loaded and 'module'
// will be left untouched. The method will still return kTRUE however since this is
// not an error. We do not need the AliMUON module when filling from hits.

	// Check if we even need the MUON module. Not having to load it will
	// save a lot of loading time for AliRoot.
	if (fDataToUse == kFromHits)
	{
		// Make sure we do not attempt to unload gAlice in FinishedWithAliMUON,
		// by setting the fHadToLoadgAlice to false.
		fHadToLoadgAlice = kFALSE;
		return kTRUE;
	}
	
	AliRunLoader* runloader = AliRunLoader::GetRunLoader();
	if ( runloader == NULL )
	{
		Error("FetchAliMUON", "AliRunLoader not initialised!");
		return kFALSE;
	}

	// Try fetch the AliRun object. If it is not found then try load it using
	// the runloader.
	AliRun* alirun = runloader->GetAliRun();
	if (alirun == NULL)
	{
		if (runloader->LoadgAlice() != 0)
		{
			// Error.
			DebugMsg(1, "Leaving FillFrom(AliMUONDataInterface*)");
			return kFALSE;
		}
		fHadToLoadgAlice = kTRUE;
		alirun = runloader->GetAliRun();
	}
	else
		fHadToLoadgAlice = kFALSE;

	// Get the MUON module pointer and return it.
	module = dynamic_cast<AliMUON*>( alirun->GetModule("MUON") );
	return kTRUE;
}


void AliHLTMUONTriggerSource::FinishedWithAliMUON()
{
// After one is finished with the AliMUON object returned by GetAliMUON, one
// should call this method.
// If the gAlice object was loaded by GetAliMUON then it will be unloaded at
// this point, otherwise nothing is done.

	// Only unload the gAlice object if we had to load it ourselves.
	if (fHadToLoadgAlice)
		AliRunLoader::GetRunLoader()->UnloadgAlice();
}


void AliHLTMUONTriggerSource::ResetAllPointers() const
{
// Sets all the current pointers to NULL and indices to -1.

	fEventIndex = -1;
	fCurrentEvent = NULL;
	fBlockIndex = -1;
	fCurrentBlock = NULL;
	fTriggerIndex = -1;
	fCurrentTrigger = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


void AliHLTMUONTriggerSource::ResetBlockPointers() const
{
// Sets the block and trigger pointers to NULL and indices to -1.

	fBlockIndex = -1;
	fCurrentBlock = NULL;
	fTriggerIndex = -1;
	fCurrentTrigger = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


void AliHLTMUONTriggerSource::ResetTriggerPointers() const
{
// Sets just the current trigger record pointer to NULL and index to -1.

	fTriggerIndex = -1;
	fCurrentTrigger = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fTriggerIndex = " << fTriggerIndex
	);
}


AliHLTMUONTriggerSource::AliEventData::AliEventData() :
	fEventNumber(-1), fBlocks(TClonesArray::Class())
{
	fEventNumber = -1;
}


AliHLTMUONTriggerSource::AliEventData::AliEventData(Int_t eventnumber)
	: fEventNumber(eventnumber), fBlocks(TClonesArray::Class())
{
// Create a new event data block with specified event number.

	fEventNumber = eventnumber;

	// If the following is not set then we do not write the fBlocks properly.
	fBlocks.BypassStreamer(kFALSE);
}


AliHLTMUONTriggerSource::AliEventData::~AliEventData()
{
	//fBlocks.Clear("C");  // Done in fBlocks destructor 
}

