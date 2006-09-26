////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* AliHLTMUONClusterSource is used to extract cluster points from the data
   stored in a root file for a AliRoot simulated event.
   It is used by AliHLTMUONMicrodHLT as the input data set object.
 */

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
ClassImp(AliHLTMUONClusterSource::AliBlockData)
ClassImp(AliHLTMUONClusterSource::AliEventData)


AliHLTMUONClusterSource::AliHLTMUONClusterSource()
	: TObject(),
	  fAreaToUse(kFromWholePlane), fDataToUse(kFromRawClusters),
	  fMaxBlockSize(0xFFFFFFFF), fFilename(""), fFoldername(""),
	  fEventIndex(-1), fCurrentEvent(NULL), fBlockIndex(-1),
	  fCurrentBlock(NULL), fClusterIndex(-1), fCurrentCluster(NULL),
	  fEventList(AliHLTMUONClusterSource::AliEventData::Class())
{
// Default contructor.

	fAreaToUse = kFromWholePlane;
	fDataToUse = kFromRawClusters;
	fMaxBlockSize = 0xFFFFFFFF;
	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
}


AliHLTMUONClusterSource::AliHLTMUONClusterSource(AliMUONDataInterface* data)
	: TObject(),
	  fAreaToUse(kFromWholePlane), fDataToUse(kFromRawClusters),
	  fMaxBlockSize(0xFFFFFFFF), fFilename(""), fFoldername(""),
	  fEventIndex(-1), fCurrentEvent(NULL), fBlockIndex(-1),
	  fCurrentBlock(NULL), fClusterIndex(-1), fCurrentCluster(NULL),
	  fEventList(AliHLTMUONClusterSource::AliEventData::Class())
{
// Creates a new cluster source object and fills it with data from 'data'.

	fAreaToUse = kFromWholePlane;
	fDataToUse = kFromRawClusters;
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
// Fills the internal data structures from the specified data interface
// for all the events found in AliMUONDataInterface.

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
// Fills the internal data structures from the specified data interface
// for the given event.

	DebugMsg(1, "FillFrom(AliMUONDataInterface*, Int_t)");
	
	if (FileAndFolderOk(data))
		AddEventFrom(data, event);
}


void AliHLTMUONClusterSource::FillFrom(AliMUONDataInterface* data, Int_t event, Int_t chamber)
{
// Fills the internal data structures from the specified data interface
// for the given event and chamber.

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
// Fills the internal data structures from the specified data interface
// for the given event, chamber and cluster number.
// If 'newblock' is set to true then the new cluster point is added to
// a new block. Otherwise the point is added to the current block.
// Note: This method ignores the fAreaToUse and fMaxBlockSize flags.
// This is very usefull for custom cluster source filling.
// For the case of adding data from AliMUONHit objects the 'cluster' parameter
// becomes the track number in TreeH and not the index of the AliMUONRawCluster
// object.

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
			if (fCurrentEvent->EventNumber() != event)
			{
				Bool_t found = GetEvent(event);
				if ( ! found) AddEvent(event);
			}
		}
		
		if ( fCurrentBlock != 0)
		{
			Assert( fCurrentEvent != NULL );
			
			if ( fCurrentBlock->Chamber() != chamber)
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
// Clears all the internal arrays.

	fFilename = "";
	fFoldername = "";
	ResetAllPointers();
	fEventList.Delete();
}


Bool_t AliHLTMUONClusterSource::GetEvent(Int_t eventnumber) const
{
// Fetches the specified event number stored in this AliHLTMUONClusterSource.
// Sets the current block and cluster point to the first block and cluster
// point respectively. If there are no blocks or clusters then these pointers
// are set to NULL.
// kTRUE is returned if the event was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONClusterSource::GetEvent(" << eventnumber << ")" );
	
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
				<< " , fClusterIndex = " << fClusterIndex
			);
			return kTRUE;
		}
	}
	return kFALSE;
}


Bool_t AliHLTMUONClusterSource::GetFirstEvent() const
{
// Fetches the first event stored in this AliHLTMUONClusterSource.
// Sets the current block and cluster point to the first block and cluster
// point respectively. If there are no blocks or clusters then these pointers
// are set to NULL.
// kTRUE is returned if the event was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONClusterSource::GetFirstEvent()");
	if (fEventList.GetEntriesFast() > 0)
	{
		fEventIndex = 0;
		fCurrentEvent = (AliEventData*) fEventList[0];
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
// Returns kTRUE if there are more events to interate over.

	return 0 <= fEventIndex && fEventIndex < fEventList.GetEntriesFast();
}


Bool_t AliHLTMUONClusterSource::GetNextEvent() const
{
// Fetches the next event stored following the currently selected one.
// kTRUE is returned if the event was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONClusterSource::GetNextEvent()");
	if (fEventIndex < fEventList.GetEntriesFast() - 1)
	{
		fCurrentEvent = (AliEventData*) fEventList[ ++fEventIndex ];
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
// Returns the corresponding AliRoot event number for the current event.
// -1 is returned if no event is selected.

	if (fCurrentEvent != NULL)
		return fCurrentEvent->EventNumber();
	else
		return -1;
}


Int_t AliHLTMUONClusterSource::NumberOfBlocks() const
{
// Returns the number of cluster blocks in the current event.
// -1 is returned if no event is selected.

	DebugMsg(1, "AliHLTMUONClusterSource::NumberOfBlocks()");
	if (fCurrentEvent == NULL)
	{
		Error("NumberOfBlocks", "No event selected.");
		return -1;
	}
	else
		return fCurrentEvent->Blocks().GetEntriesFast();
}


Bool_t AliHLTMUONClusterSource::GetBlock(Int_t index) const
{
// Fetches the index'th block in the current event.
// Sets the current cluster point to the first cluster point in the block.
// If there are no cluster points then this pointer is set to NULL.
// kTRUE is returned if the block was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONClusterSource::GetBlock(" << index << ")");
	
	// Note NumberOfBlocks() also checks if the event was selected.
	Int_t numberofblocks = NumberOfBlocks();
	if (numberofblocks < 0) return kFALSE;

	if ( 0 <= index && index < numberofblocks )
	{
		fBlockIndex = index;
		fCurrentBlock = (AliBlockData*) fCurrentEvent->Blocks()[index];
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
// Fetches the first block in the current event.
// Sets the current cluster point to the first cluster point in the block.
// If there are no cluster points then this pointer is set to NULL.
// kTRUE is returned if the block was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONClusterSource::GetFirstBlock()");
	// Note: NumberOfBlocks() also checks if fCurrentEvent != NULL.
	if (NumberOfBlocks() > 0)
	{
		fBlockIndex = 0;
		fCurrentBlock = (AliBlockData*) fCurrentEvent->Blocks()[fBlockIndex];
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
// Returns kTRUE if there are more blocks to interate over.

	return 0 <= fBlockIndex && fBlockIndex < NumberOfBlocks();
}


Bool_t AliHLTMUONClusterSource::GetNextBlock() const
{
// Fetches the next block in the current event.
// kTRUE is returned if the block was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONClusterSource::GetNextBlock()");

	// Note: NumberOfBlocks() checks if fCurrentEvent != NULL. If it is then it returns -1
	// and since fBlockIndex is always >= -1 the if statement must go to the else part.
	if (fBlockIndex < NumberOfBlocks() - 1)
	{
		fCurrentBlock = (AliBlockData*) fCurrentEvent->Blocks()[ ++fBlockIndex ];
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
// Returns the chamber number of the current block.
// -1 is returned if no block is selected.

	if (fCurrentBlock == NULL)
	{
		Error("Chamber", "No block selected.");
		return -1;
	}
	else
		return fCurrentBlock->Chamber();
}


Int_t AliHLTMUONClusterSource::NumberOfClusters() const
{
// Returns the number of cluster points in the current block.
// -1 is returned if no block is selected.

	DebugMsg(1, "AliHLTMUONClusterSource::NumberOfClusters()");
	if (fCurrentBlock == NULL)
	{
		Error("NumberOfClusters", "No block selected.");
		return -1;
	}
	else
		return fCurrentBlock->Clusters().GetEntriesFast();
}


const AliHLTMUONPoint* AliHLTMUONClusterSource::GetCluster(Int_t index) const
{
// Fetches the index'th cluster point in the current block.
// kTRUE is returned if the point was found, kFALSE otherwise.

	DebugMsg(1, "AliHLTMUONClusterSource::GetCluster(" << index << ")");

	// Note NumberOfClusters() also checks if the event and block was selected.
	Int_t numberofclusters = NumberOfClusters();
	if (numberofclusters < 0) return NULL;
	
	if ( 0 <= index && index < numberofclusters )
	{
		fClusterIndex = index;
		fCurrentCluster = (AliHLTMUONPoint*) fCurrentBlock->Clusters()[index];
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
// Fetches the first cluster point in the current block.
// NULL is returned if the point was not found.

	DebugMsg(1, "AliHLTMUONClusterSource::GetFirstCluster()");
	// Note: NumberOfClusters() also checks if fCurrentBlock != NULL.
	if (NumberOfClusters() > 0)
	{
		fClusterIndex = 0;
		fCurrentCluster = (AliHLTMUONPoint*) fCurrentBlock->Clusters()[0];
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
// Returns kTRUE if there are more cluster points to interate over.

	return 0 <= fClusterIndex && fClusterIndex < NumberOfClusters();
}


const AliHLTMUONPoint* AliHLTMUONClusterSource::GetNextCluster() const
{
// Fetches the next cluster point in the current block.
// NULL is returned if the point was not found.

	DebugMsg(1, "AliHLTMUONClusterSource::GetNextCluster()");
	
	// Note: NumberOfClusters() checks if fCurrentBlock != NULL. If it is then it returns -1
	// and since fClusterIndex is always >= -1 the if statement must go to the else part.
	if (fClusterIndex < NumberOfClusters() - 1)
	{
		fCurrentCluster = (AliHLTMUONPoint*) fCurrentBlock->Clusters()[ ++fClusterIndex ];
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
// Returns the x and y coordinate of the current cluster point.
// kFALSE is returned if there is no cluster point selected.

	if (fCurrentCluster != NULL)
	{
		x = fCurrentCluster->X();
		y = fCurrentCluster->Y();
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
// Adds a new AliEventData block to the fEventList and updates the fCurrentEvent,
// fCurrentBlock and fCurrentCluster pointers.

	DebugMsg(1, "AliHLTMUONClusterSource::AddEvent(" << eventnumber << ")");
	Assert( eventnumber >= 0 );

	// Assume the eventnumber does not already exist in the event list.
	fEventIndex = fEventList.GetEntriesFast();
	new ( fEventList[fEventIndex] ) AliEventData(eventnumber);
	fCurrentEvent = (AliEventData*) fEventList[fEventIndex];
	
	// Remember to reset the other pointers because the new event is empty.
	ResetBlockPointers();
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


void AliHLTMUONClusterSource::AddBlock(Int_t chamber)
{
// Adds a new block to the current event and updates fCurrentBlock and fCurrentCluster.
// The chamber number is assigned to the blocks fChamber value.

	DebugMsg(1, "AliHLTMUONClusterSource::AddBlock()");
	
	if (fCurrentEvent == NULL)
	{
		Error("AddBlock", "No event selected.");
		return;
	}
	
	fBlockIndex = fCurrentEvent->Blocks().GetEntriesFast();
	new ( fCurrentEvent->Blocks()[fBlockIndex] ) AliBlockData(chamber);
	fCurrentBlock = (AliBlockData*) fCurrentEvent->Blocks()[fBlockIndex];
	
	// Remember to reset the trigger pointer because the new block is empty.
	ResetClusterPointers();

	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


void AliHLTMUONClusterSource::AddPoint(Float_t x, Float_t y)
{
// Adds a new cluster point to the current event and block.
// The fCurrentCluster is updated appropriately.

	DebugMsg(1, "AliHLTMUONClusterSource::AddPoint(" << x << ", " << y << ")");

	if (fCurrentBlock == NULL)
	{
		Error("AddPoint", "No block selected.");
		return;
	}
	
	fClusterIndex = fCurrentBlock->Clusters().GetEntriesFast();
	new ( fCurrentBlock->Clusters()[fClusterIndex] ) AliHLTMUONPoint(x, y);
	fCurrentCluster = (AliHLTMUONPoint*) fCurrentBlock->Clusters()[fClusterIndex];
	
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


Bool_t AliHLTMUONClusterSource::FileAndFolderOk(AliMUONDataInterface* data)
{
// Checks if the file and folder names correspond to this AliHLTMUONClusterSource's
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
// Adds the whole event from the data interface to the internal data structures.
// It is assumed that FileAndFolderOk(data) returns true just before calling
// this method.

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
// Adds all cluster points found on the given chamber from the specified data
// interface. The data interface should be initialised correctly, that is the
// event should already be selected before calling this method.

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
	case kFromHits:
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

	case kFromRawClusters:
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
// Adds the cluster point from the specified data interface.
// The data interface should be initialised correctly, that is the event
// should already be selected before calling this method.

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
	case kFromHits:
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

	case kFromRawClusters:
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


Bool_t AliHLTMUONClusterSource::InFillRegion(Float_t x, Float_t /*y*/) const
{
// Checks to see if the x and y coordinate of the cluster point are in the
// chamber region we want to fill from.
// kTRUE is returned if (x, y) is in the region, and kFALSE otherwise.

	switch (fAreaToUse)
	{
	case kFromWholePlane:     return kTRUE;
	case kFromLeftHalfPlane:  return x <= 0;
	case kFromRightHalfPlane: return x > 0;

	default:
		Error("InFillRegion", "fAreaToUse is not set to a valid value.");
		return kFALSE;
	}
}


void AliHLTMUONClusterSource::ResetAllPointers() const
{
// Sets all the current pointers to NULL and indices to -1.

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
// Sets the block and trigger pointers to NULL and indices to -1.

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
// Sets just the current cluster point pointer to NULL and index to -1.

	fClusterIndex = -1;
	fCurrentCluster = NULL;
	DebugMsg(2, "\tfEventIndex = " << fEventIndex << " , fBlockIndex = " << fBlockIndex
		<< " , fClusterIndex = " << fClusterIndex
	);
}


AliHLTMUONClusterSource::AliBlockData::AliBlockData()
	: fChamber(-1), fClusters(AliHLTMUONPoint::Class())
{
	fChamber = -1;
}

AliHLTMUONClusterSource::AliBlockData::AliBlockData(Int_t chamber)
	: fChamber(chamber), fClusters(AliHLTMUONPoint::Class())
{
	fChamber = chamber;
}

AliHLTMUONClusterSource::AliBlockData::~AliBlockData()
{
	fClusters.Clear("C");
}

AliHLTMUONClusterSource::AliEventData::AliEventData()
	: fEventNumber(-1), fBlocks(AliHLTMUONClusterSource::AliBlockData::Class())
{
	fEventNumber = -1;
}

AliHLTMUONClusterSource::AliEventData::AliEventData(Int_t eventnumber)
	: fEventNumber(eventnumber), fBlocks(AliHLTMUONClusterSource::AliBlockData::Class())
{
// Creates a new event data block with specified event number.

	fEventNumber = eventnumber;

	// If the following is not set then we do not write the fBlocks properly.
	fBlocks.BypassStreamer(kFALSE);
}

AliHLTMUONClusterSource::AliEventData::~AliEventData()
{
	fBlocks.Clear("C");
}

