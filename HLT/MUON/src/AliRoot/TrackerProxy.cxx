////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/TrackerProxy.hpp"
#include "AliRoot/convert.hpp"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONOutOfMemory.h"


AliHLTMUONTrackerProxy::AliHLTMUONTrackerProxy(AliHLTMUONTrackerInterface* client)
	: AliHLTMUONCoreTracker(), AliHLTMUONTrackerCallback(), fTracker(client)
{
	fTracker = client;
}


void AliHLTMUONTrackerProxy::FindTrack(const AliHLTMUONCoreTriggerRecord& trigger)
{
// Finds a track from the tigger seed by invoking the tracker object.

	AliHLTMUONTriggerRecord rec = AliHLTMUONConvert(trigger, 0);
	DebugMsg(6, "AliHLTMUONTrackerProxy::FindTrack : rec = " << rec);
	fTracker->FindTrack(rec);
}


void AliHLTMUONTrackerProxy::ReturnClusters(void* tag, const AliHLTMUONCoreClusterPoint* clusters, UInt count)
{
// Passes the cluster points the tracker requested to the tracker.

	AliHLTMUONPoint* points = new AliHLTMUONPoint[count];
	try
	{
		DebugMsg(6, "AliHLTMUONTrackerProxy::ReturnClusters");
		for (UInt i = 0; i < count; i++)
		{
			points[i] = AliHLTMUONConvert(clusters[i]);
			DebugMsg(6, "\tpoints[" << i << "] = " << points[i] );
		};
		fTracker->ReturnClusters(tag, points, count);
	}
	finally
	(
		delete [] points;
	)
}


void AliHLTMUONTrackerProxy::EndOfClusters(void* tag)
{
// Tells the tracker there are no more clusters.

	DebugMsg(6, "AliHLTMUONTrackerProxy::EndOfClusters");
	fTracker->EndOfClusters(tag);
}


void AliHLTMUONTrackerProxy::FillTrackData(AliHLTMUONCoreTrack& track)
{
// Fills the tracker data for the found track.

	AliHLTMUONTrack data;
	fTracker->FillTrackData(data);
	DebugMsg(6, "AliHLTMUONTrackerProxy::FillTrackData : data = " << data);
	track = AliHLTMUONConvert(data);
}


void AliHLTMUONTrackerProxy::Reset()
{
	DebugMsg(6, "AliHLTMUONTrackerProxy::Reset");
	fTracker->Reset();
}


void AliHLTMUONTrackerProxy::RequestClusters(
		Float_t left, Float_t right, Float_t bottom, Float_t top,
		Int_t chamber, const void* tag
	)
{
	DebugMsg(6, "AliHLTMUONTrackerProxy::RequestClusters");
	AliHLTMUONCoreTracker::RequestClusters(left, right, bottom, top, (AliHLTMUONCoreChamberID)chamber, tag);
}


void AliHLTMUONTrackerProxy::EndOfClusterRequests()
{
	DebugMsg(6, "AliHLTMUONTrackerProxy::EndOfClusterRequests");
	AliHLTMUONCoreTracker::EndOfClusterRequests();
}


void AliHLTMUONTrackerProxy::FoundTrack()
{
	DebugMsg(6, "AliHLTMUONTrackerProxy::FoundTrack");
	AliHLTMUONCoreTracker::FoundTrack();
}


void AliHLTMUONTrackerProxy::NoTrackFound()
{
	DebugMsg(6, "AliHLTMUONTrackerProxy::NoTrackFound");
	AliHLTMUONCoreTracker::NoTrackFound();
}

