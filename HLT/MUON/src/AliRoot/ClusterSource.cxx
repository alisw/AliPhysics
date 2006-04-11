////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/ClusterSource.hpp"
#include "AliRoot/Base.hpp"
#include "AliMUONConstants.h"
#include "AliMUONHit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONDataInterface.h"
#ifndef __alpha
#include <math.h>
#else
#include <float.h>
#endif

ClassImp(AliHLTMUONClusterSource)
ClassImp(AliHLTMUONClusterSource::BlockData)
ClassImp(AliHLTMUONClusterSource::EventData)


AliHLTMUONClusterSource::AliHLTMUONClusterSource()
	: TObject(), fEventList(AliHLTMUONClusterSource::EventData::Class())
{
	fAreaToUse = FromWholePlane;
	fDataToUse = FromRawClusters;
	fMaxBlockSize = 0xFFFFFFFF;
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
}


AliHLTMUONClusterSource::AliHLTMUONClusterSource(AliMUONDataInterface* data)
	: TObject(), fEventList(AliHLTMUONClusterSource::EventData::Class())
{
	fAreaToUse = FromWholePlane;
	fDataToUse = FromRawClusters;
	fMaxBlockSize = 0xFFFFFFFF;
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	FillFrom(data);
}


AliHLTMUONClusterSource::~AliHLTMUONClusterSource()
{
	fEventList.Clear("C");
}


void AliHLTMUONClusterSource::FillFrom(AliMUONDataInterface* data)
{
	DebugMsg(1, "FillFrom(AliMUONDataInterface*)");
	
	if (FileAndFolderOk(data))
	{
		for (Int_t i = 0; i < data->NumberOfEvents(); i++)
		{
			AddEventFrom(data, i);
		}
	}
}


void AliHLTMUONClusterSource::FillFrom(AliMUONDataInterface* data, Int_t event)
{
	DebugMsg(1, "FillFrom(AliMUONDataInterface*, Int_t)");
	
	if (FileAndFolderOk(data))
		AddEventFrom(data, event);
}


void AliHLTMUONClusterSource::FillFrom(AliMUONDataInterface* data, Int_t event, Int_t chamber)
{
	DebugMsg(1, "FillFrom(AliMUONDataInterface*, Int_t)");
	
	if (FileAndFolderOk(data))
	{
		if ( data->GetEvent(event) )
		{
			AddEvent(event);
			AddChamberFrom(data, chamber);
		}
	}
}


void AliHLTMUONClusterSource::FillFrom(
		AliMUONDataInterface* data,
		Int_t event, Int_t chamber, Int_t cluster,
		Bool_t newblock
	)
{
	DebugMsg(1, "FillFrom(AliMUONDataInterface*, Int_t, Int_t, Int_t, Bool_t)");
	
	if (FileAndFolderOk(data))
	{
		data->GetEvent(event);

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
			if (fCurrentEvent->fEventNumber != event)
			{
				Bool_t found = GetEvent(event);
				if ( ! found) AddEvent(event);
			}
		}
		
		if ( fCurrentBlock != 0)
		{
			Assert( fCurrentEvent != NULL );
			
			if ( fCurrentBlock->fChamber != chamber)
			{
				// Create a new block if the current blocks chamber number does
				// not correspond to the specified chamber.
				AddBlock(chamber);
			}
			else
			{
				// If the newblock flag is set then force a new block.
				if (newblock) AddBlock(chamber);
			}
		}
		else
			AddBlock(chamber);  // No block selected so we need to create a new block.

		AddClusterFrom(data, chamber, cluster);
	}
}


void AliHLTMUONClusterSource::Clear(Option_t* /*option*/)
{
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fEventList.Delete();
}


Bool_t AliHLTMUONClusterSource::GetEvent(Int_t eventnumber) const
{
	DebugMsg(1, "AliHLTMUONClusterSource::GetEvent(" << eventnumber << ")" );
	
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
				<< " , fClusterIndex = " << fClusterIndex
			);
			return kTRUE;
		}
	}
	return kFALSE;
}


Bool_t AliHLTMUONClusterSource::GetFirstEvent() const
{
	DebugMsg(1, "AliHLTMUONClusterSource::GetFirstEvent()");
	if (fEventList.GetEntriesFast() > 0)
	{
		fEventIndex = 0;
		fCurrentEvent = (EventData*) fEventList[0];
		GetFirstBlock();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fClusterIndex = " << fClusterIndex
		);
		return kTRUE;
	}
	else
	{
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fClusterIndex = " << fClusterIndex
		);
		return kFALSE;
	}
}


Bool_t AliHLTMUONClusterSource::MoreEvents() const
{
	return 0 <= fEventIndex && fEventIndex < fEventList.GetEntriesFast();
}


Bool_t AliHLTMUONClusterSource::GetNextEvent() const
{
	DebugMsg(1, "AliHLTMUONClusterSource::GetNextEvent()");
	if (fEventIndex < fEventList.GetEntriesFast() - 1)
	{
		fCurrentEvent = (EventData*) fEventList[ ++fEventIndex ];
		GetFirstBlock();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fClusterIndex = " << fClusterIndex
		);
		return kTRUE;
	}
	else
	{
		ResetAllPointers();
		return kFALSE;
	}
}


Int_t AliHLTMUONClusterSource::CurrentEvent() const
{
	if (fCurrentEvent != NULL)
		return fCurrentEvent->fEventNumber;
	else
		return -1;
}


Int_t AliHLTMUONClusterSource::NumberOfBlocks() const
{
	DebugMsg(1, "AliHLTMUONClusterSource::NumberOfBlocks()");
	if (fCurrentEvent == NULL)
	{
		Error("NumberOfBlocks", "No event selected.");
		return -1;
	}
	else
		return fCurrentEvent->fBlocks.GetEntriesFast();
}


Bool_t AliHLTMUONClusterSource::GetBlock(Int_t index) const
{
	DebugMsg(1, "AliHLTMUONClusterSource::GetBlock(" << index << ")");
	
	// Note NumberOfBlocks() also checks if the event was selected.
	Int_t numberofblocks = NumberOfBlocks();
	if (numberofblocks < 0) return kFALSE;

	if ( 0 <= index && index < numberofblocks )
	{
		fBlockIndex = index;
		fCurrentBlock = (BlockData*) fCurrentEvent->fBlocks[index];
		GetFirstCluster();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fClusterIndex = " << fClusterIndex
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


Bool_t AliHLTMUONClusterSource::GetFirstBlock() const
{
	DebugMsg(1, "AliHLTMUONClusterSource::GetFirstBlock()");
	// Note: NumberOfBlocks() also checks if fCurrentEvent != NULL.
	if (NumberOfBlocks() > 0)
	{
		fBlockIndex = 0;
		fCurrentBlock = (BlockData*) fCurrentEvent->fBlocks[fBlockIndex];
		GetFirstCluster();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fClusterIndex = " << fClusterIndex
		);
		return kTRUE;
	}
	else
		return kFALSE;
}


Bool_t AliHLTMUONClusterSource::MoreBlocks() const
{
	return 0 <= fBlockIndex && fBlockIndex < NumberOfBlocks();
}


Bool_t AliHLTMUONClusterSource::GetNextBlock() const
{
	DebugMsg(1, "AliHLTMUONClusterSource::GetNextBlock()");

	// Note: NumberOfBlocks() checks if fCurrentEvent != NULL. If it is then it returns -1
	// and since fBlockIndex is always >= -1 the if statement must go to the else part.
	if (fBlockIndex < NumberOfBlocks() - 1)
	{
		fCurrentBlock = (BlockData*) fCurrentEvent->fBlocks[ ++fBlockIndex ];
		GetFirstCluster();
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fClusterIndex = " << fClusterIndex
		);
		return kTRUE;
	}
	else
	{
		ResetBlockPointers();
		return kFALSE;
	}
}


Int_t AliHLTMUONClusterSource::Chamber() const
{
	if (fCurrentBlock == NULL)
	{
		Error("Chamber", "No block selected.");
		return -1;
	}
	else
		return fCurrentBlock->fChamber;
}


Int_t AliHLTMUONClusterSource::NumberOfClusters() const
{
	DebugMsg(1, "AliHLTMUONClusterSource::NumberOfClusters()");
	if (fCurrentBlock == NULL)
	{
		Error("NumberOfClusters", "No block selected.");
		return -1;
	}
	else
		return fCurrentBlock->fClusters.GetEntriesFast();
}


const AliHLTMUONPoint* AliHLTMUONClusterSource::GetCluster(Int_t index) const
{
	DebugMsg(1, "AliHLTMUONClusterSource::GetCluster(" << index << ")");

	// Note NumberOfClusters() also checks if the event and block was selected.
	Int_t numberofclusters = NumberOfClusters();
	if (numberofclusters < 0) return NULL;
	
	if ( 0 <= index && index < numberofclusters )
	{
		fClusterIndex = index;
		fCurrentCluster = (AliHLTMUONPoint*) fCurrentBlock->fClusters[index];
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fClusterIndex = " << fClusterIndex
		);
		return fCurrentCluster;
	}
	else
	{
		// The index is out of bounds so inform the user.
		if (numberofclusters > 0)
			Error(	"GetCluster",
				"The cluster index (%d) is out of bounds. Valid range is [0, %d]",
				index, numberofclusters - 1
			);
		else
			Error(	"GetCluster",
				"The cluster index (%d) is out of bounds. No clusters found.",
				index
			);
		return NULL;
	}
}


const AliHLTMUONPoint* AliHLTMUONClusterSource::GetFirstCluster() const
{
	DebugMsg(1, "AliHLTMUONClusterSource::GetFirstCluster()");
	// Note: NumberOfClusters() also checks if fCurrentBlock != NULL.
	if (NumberOfClusters() > 0)
	{
		fClusterIndex = 0;
		fCurrentCluster = (AliHLTMUONPoint*) fCurrentBlock->fClusters[0];
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fClusterIndex = " << fClusterIndex
		);
		return fCurrentCluster;
	}
	else
		return NULL;
}


Bool_t AliHLTMUONClusterSource::MoreClusters() const
{
	return 0 <= fClusterIndex && fClusterIndex < NumberOfClusters();
}


const AliHLTMUONPoint* AliHLTMUONClusterSource::GetNextCluster() const
{
	DebugMsg(1, "AliHLTMUONClusterSource::GetNextCluster()");
	
	// Note: NumberOfClusters() checks if fCurrentBlock != NULL. If it is then it returns -1
	// and since fClusterIndex is always >= -1 the if statement must go to the else part.
	if (fClusterIndex < NumberOfClusters() - 1)
	{
		fCurrentCluster = (AliHLTMUONPoint*) fCurrentBlock->fClusters[ ++fClusterIndex ];
		DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
			<< " , fClusterIndex = " << fClusterIndex
		);
		return fCurrentCluster;
	}
	else
	{
		ResetClusterPointers();
		return NULL;
	}
}


Bool_t AliHLTMUONClusterSource::FetchCluster(Float_t& x, Float_t& y) const
{
	if (fCurrentCluster != NULL)
	{
		x = fCurrentCluster->fX;
		y = fCurrentCluster->fY;
		return kTRUE;
	}
	else
	{
		Error("FetchCluster", "No cluster point selected.");
		return kFALSE;
	}
}


void AliHLTMUONClusterSource::AddEvent(Int_t eventnumber)
{
	DebugMsg(1, "AliHLTMUONClusterSource::AddEvent(" << eventnumber << ")");
	Assert( eventnumber >= 0 );

	// Assume the eventnumber does not already exist in the event list.
	fEventIndex = fEventList.GetEntriesFast();
	new ( fEventList[fEventIndex] ) EventData(eventnumber);
	fCurrentEvent = (EventData*) fEventList[fEventIndex];
	
	// Remember to reset the other pointers because the new event is empty.
	ResetBlockPointers();
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


void AliHLTMUONClusterSource::AddBlock(Int_t chamber)
{
	DebugMsg(1, "AliHLTMUONClusterSource::AddBlock()");
	
	if (fCurrentEvent == NULL)
	{
		Error("AddBlock", "No event selected.");
		return;
	}
	
	fBlockIndex = fCurrentEvent->fBlocks.GetEntriesFast();
	new ( fCurrentEvent->fBlocks[fBlockIndex] ) BlockData(chamber);
	fCurrentBlock = (BlockData*) fCurrentEvent->fBlocks[fBlockIndex];
	
	// Remember to reset the trigger pointer because the new block is empty.
	ResetClusterPointers();

	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


void AliHLTMUONClusterSource::AddPoint(Float_t x, Float_t y)
{
	DebugMsg(1, "AliHLTMUONClusterSource::AddPoint(" << x << ", " << y << ")");

	if (fCurrentBlock == NULL)
	{
		Error("AddPoint", "No block selected.");
		return;
	}
	
	fClusterIndex = fCurrentBlock->fClusters.GetEntriesFast();
	new ( fCurrentBlock->fClusters[fClusterIndex] ) AliHLTMUONPoint(x, y);
	fCurrentCluster = (AliHLTMUONPoint*) fCurrentBlock->fClusters[fClusterIndex];
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


Bool_t AliHLTMUONClusterSource::FileAndFolderOk(AliMUONDataInterface* data)
{
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
			"The cluster source already contains data from file '%s', cannot add data from file '%s'",
			fFilename.Data(), data->CurrentFile().Data()
		);
		return kFALSE;
	}
	
	if ( fFoldername != data->CurrentFolder() )
	{
		Error(	"FileAndFolderOk",
			"The cluster source already contains data from folder '%s', cannot add data from folder '%s'",
			fFoldername.Data(), data->CurrentFolder().Data()
		);
		return kFALSE;
	}
	
	return kTRUE;
}


void AliHLTMUONClusterSource::AddEventFrom(AliMUONDataInterface* data, Int_t event)
{
	if ( data->GetEvent(event) )
	{
		AddEvent(event);
		for (Int_t chamber = 0; chamber < AliMUONConstants::NTrackingCh(); chamber++)
		{
			AddChamberFrom(data, chamber);
		}
	}
}


void AliHLTMUONClusterSource::AddChamberFrom(AliMUONDataInterface* data, Int_t chamber)
{
	DebugMsg(1, "Entering AddChamberFrom");
	
	AddBlock(chamber);
	UInt_t currentblocksize = 0;
#ifndef __alpha
#ifndef __sun
	Float_t x = nanf(""), y = nanf("");
#else
	Float_t x = 0, y = 0;
#endif
#else
	Float_t x = FLT_QNAN, y = FLT_QNAN;
#endif
	
	switch (fDataToUse)
	{
	case FromHits:
		for (Int_t track = 0; track < data->NumberOfTracks(); track++)
		{
			// Find the hit that corresponds to the current chamber number.
			Int_t i;
			for (i = 0; i < data->NumberOfHits(track); i++)
			{
				AliMUONHit* h = data->Hit(track, i);
				// Note AliMUONHit::Chamber() returns a value in the range 1..14
				if (h->Chamber() == chamber + 1)
				{
					x = h->X();
					y = h->Y();
					break;
				}
			}

			// Continue to the next track if we could not find a hit
			// on the current track and chamber.
			if (i >= data->NumberOfHits(track))
			{
				Warning("AddChamberFrom", "Could not find hit on chamber: %d , track: %d , event: %d", 
					chamber, track, data->CurrentEvent()
				);
				continue;
			}

			if (InFillRegion(x, y))
			{
				AddPoint(x, y);

				// Create a new block if we reached the maximum block size.
				if ( ++currentblocksize == fMaxBlockSize )
				{
					AddBlock(chamber);
					currentblocksize = 0;
				}
			}
		}
		break;

	case FromRawClusters:
		for (Int_t i = 0; i < data->NumberOfRawClusters(chamber); i++)
		{
			AliMUONRawCluster* rc = data->RawCluster(chamber, i);
			x = rc->GetX(0);
			y = rc->GetY(0);

			if (InFillRegion(x, y))
			{
				AddPoint(x, y);

				// Create a new block if we reached the maximum block size.
				if ( ++currentblocksize == fMaxBlockSize )
				{
					AddBlock(chamber);
					currentblocksize = 0;
				}
			}
		}
		break;

	default:
		Error("AddChamberFrom", "fDataToUse is not set to a valid value.");
	}
	
	DebugMsg(1, "Leaving AddChamberFrom");
}


void AliHLTMUONClusterSource::AddClusterFrom(
		AliMUONDataInterface* data, Int_t chamber, Int_t cluster
	)
{
	DebugMsg(1, "Entering AddClusterFrom");
#ifndef __alpha
#ifndef __sun	
	Float_t x = nanf(""), y = nanf("");
#else
	Float_t x = 0, y = 0;
#endif
#else
	Float_t x = FLT_QNAN, y = FLT_QNAN;
#endif

	switch (fDataToUse)
	{
	case FromHits:
		{
		Int_t i;
		// Note: cluster is now treated as the track number.
		for (i = 0; i < data->NumberOfHits(cluster); i++)
		{
			AliMUONHit* h = data->Hit(cluster, i);
			// Note AliMUONHit::Chamber() returns a value in the range 1..14
			if (h->Chamber() == chamber + 1)
			{
				x = h->X();
				y = h->Y();
				break;
			}
		}

		if (i >= data->NumberOfHits(cluster))
		{
			// Could not find a hit on the specified chamber so just return.
			Warning("AddClusterFrom", "Could not find hit on chamber: %d , track: %d , event: %d", 
				chamber, cluster, data->CurrentEvent()
			);
			DebugMsg(1, "Leaving AddClusterFrom");
			return;
		}
		}
		break;

	case FromRawClusters:
		{
		AliMUONRawCluster* rc = data->RawCluster(chamber, cluster);
		x = rc->GetX(0);
		y = rc->GetY(0);
		}
		break;

	default:
		Error("AddClusterFrom", "fDataToUse is not set to a valid value.");
		return;
	}
	
	AddPoint(x, y);
	
	DebugMsg(1, "Leaving AddClusterFrom");
}


Bool_t AliHLTMUONClusterSource::InFillRegion(Float_t x, Float_t /*y*/)
{
	switch (fAreaToUse)
	{
	case FromWholePlane:     return kTRUE;
	case FromLeftHalfPlane:  return x <= 0;
	case FromRightHalfPlane: return x > 0;

	default:
		Error("InFillRegion", "fAreaToUse is not set to a valid value.");
		return kFALSE;
	}
}


void AliHLTMUONClusterSource::ResetAllPointers() const
{
	fEventIndex = -1;
	fCurrentEvent = NULL;
	fBlockIndex = -1;
	fCurrentBlock = NULL;
	fClusterIndex = -1;
	fCurrentCluster = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


void AliHLTMUONClusterSource::ResetBlockPointers() const
{
	fBlockIndex = -1;
	fCurrentBlock = NULL;
	fClusterIndex = -1;
	fCurrentCluster = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


void AliHLTMUONClusterSource::ResetClusterPointers() const
{
	fClusterIndex = -1;
	fCurrentCluster = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


AliHLTMUONClusterSource::BlockData::BlockData() : fClusters(AliHLTMUONPoint::Class())
{
	fChamber = -1;
}

AliHLTMUONClusterSource::BlockData::BlockData(Int_t chamber) : fClusters(AliHLTMUONPoint::Class())
{
	fChamber = chamber;
}

AliHLTMUONClusterSource::BlockData::~BlockData()
{
	fClusters.Clear("C");
}

AliHLTMUONClusterSource::EventData::EventData() : fBlocks(AliHLTMUONClusterSource::BlockData::Class())
{
	fEventNumber = -1;
}

AliHLTMUONClusterSource::EventData::EventData(Int_t eventnumber)
	: fBlocks(AliHLTMUONClusterSource::BlockData::Class())
{
	fEventNumber = eventnumber;

	// If the following is not set then we do not write the fBlocks properly.
	fBlocks.BypassStreamer(kFALSE);
}

AliHLTMUONClusterSource::EventData::~EventData()
{
	fBlocks.Clear("C");
}

