////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Tracking/EventHandler.hpp"
#include "Tracking/IOHandler.hpp"
#include "Error.hpp"

namespace dHLT
{
namespace Tracking
{

void EventHandler::AddTriggers(const TriggerRecord* triggers, const UInt count)
{
	Assert( iohandler != NULL );
	newtrackcount += count;
	for (UInt i = 0; i < count; i++)
	{
		Tracker* tracker = iohandler->NewTracker();
		tracker->SetCallback(this);
		trackers.Add(tracker);
		tracker->FindTrack(triggers[i]);
	};
};


void EventHandler::EndOfTriggers()
{
	endoftriggers = true;
	CheckForTotalEndOfRequests();
	CheckForEndOfTracks();
};


void EventHandler::AddClusters(const ROI region, const ClusterPoint* clusters, const UInt count)
{
	typedef RegionRegistrationTable::TrackerRecordList TrackerRecordList;
	
	ChamberID chamber; UChar level; UInt left, bottom;
	RegionOfInterest::Decode(region, chamber, level, left, bottom);

	TrackerRecordList foundlist;
	registration_table.AddClusters(clusters, count, chamber, left, bottom, level, foundlist);

	for (TrackerRecordList::Iterator rec = foundlist.First(); rec != foundlist.End(); rec++)
	{
		rec->tracker->ReturnClusters(rec->tag, clusters, count);
	};
};


void EventHandler::EndOfClusters(const ROI region)
{
	typedef RegionRegistrationTable::TrackerRecordList TrackerRecordList;

	ChamberID chamber; UChar level; UInt left, bottom;
	RegionOfInterest::Decode(region, chamber, level, left, bottom);

	TrackerRecordList foundlist;
	registration_table.MarkEndOfClusters(chamber, left, bottom, level, foundlist);

	for (TrackerRecordList::Iterator rec = foundlist.First(); rec != foundlist.End(); rec++)
	{
		rec->tracker->EndOfClusters(rec->tag);
	};
};


void EventHandler::RequestClusters(
		Tracker* tracker,
		const Float left, const Float right, const Float bottom, const Float top,
		const ChamberID chamber, const void* ctag
	)
{
	typedef RegionRegistrationTable::ClusterBlockList ClusterBlockList;

	Assert( chamber < 10 );
	Assert( left <= right );
	Assert( bottom <= top );
	Assert( iohandler != NULL );

	register void* tag = const_cast<void*>(ctag);

	RegionOfInterest region(left, right, bottom, top, chamber);
	UChar level; UInt l, b;
	ROI roi = region.Encode(level, l, b);

	ClusterBlockList foundlist;
	bool endofclusters = registration_table.AddRequest(tracker, tag, chamber, l, b, level, foundlist);

	if (foundlist.Empty())
	{
		iohandler->RequestClusters(eventid, roi);
	}
	else
	{
		for (ClusterBlockList::Iterator block = foundlist.First(); block != foundlist.End(); block++)
		{
			tracker->ReturnClusters(tag, block->clusters, block->count);
		};
		if (endofclusters)
			tracker->EndOfClusters(tag);
	};
};


void EventHandler::EndOfClusterRequests(Tracker* tracker)
{
	endofrequests_count++;
	// Might want to clean up some of the internal data here?
	CheckForTotalEndOfRequests();
};


void EventHandler::FoundTrack(Tracker* tracker)
{
	Assert( iohandler != NULL );
	if (trackblock == NULL)
	{
		UInt maxtracks = iohandler->GetMaxTracksInBlock();
		// maxtrackcount = min(maxtracks, newtrackcount)
		maxtrackcount = maxtracks < newtrackcount ? maxtracks : newtrackcount;
		Assert( maxtrackcount > 0 );
		trackblock = iohandler->AllocateTrackBlock( sizeof(Track) * maxtrackcount );
		currenttrackcount = 0;
		newtrackcount -= maxtrackcount;
	};

	tracker->FillTrackData(trackblock[currenttrackcount++]);

	if (currenttrackcount == maxtrackcount)
	{
		iohandler->ReturnTracks(eventid, trackblock, currenttrackcount);
		trackblock = NULL;
		CheckForEndOfTracks();
	};
};


void EventHandler::NoTrackFound(Tracker* tracker)
{
	// TODO
};


void EventHandler::CheckForTotalEndOfRequests()
{
	if ( endoftriggers and endofrequests_count == trackers.Count() )
	{
		Assert( iohandler != NULL );
		iohandler->EndOfClusterRequests(eventid);
	};
};


void EventHandler::CheckForEndOfTracks()
{
	Assert( iohandler != NULL );
	if (endoftriggers and newtrackcount == 0 )
		iohandler->EndOfTracks(eventid);
};


} // Tracking
} // dHLT
