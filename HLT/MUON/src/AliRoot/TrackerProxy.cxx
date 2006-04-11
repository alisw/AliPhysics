////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/TrackerProxy.hpp"
#include "AliRoot/convert.hpp"
#include "Utils.hpp"
#include "new.hpp"


AliHLTMUONTrackerProxy::AliHLTMUONTrackerProxy(AliHLTMUONTrackerInterface* client)
	: AliHLTMUONCoreTracker(), AliHLTMUONTrackerCallback()
{
	tracker = client;
}


void AliHLTMUONTrackerProxy::FindTrack(const AliHLTMUONCoreTriggerRecord& trigger)
{
	AliHLTMUONTriggerRecord rec = AliHLTMUONConvert(trigger, 0);
	DebugMsg(6, "AliHLTMUONTrackerProxy::FindTrack : rec = " << rec);
	tracker->FindTrack(rec);
}


void AliHLTMUONTrackerProxy::ReturnClusters(void* tag, const AliHLTMUONCoreClusterPoint* clusters, UInt count)
{
	AliHLTMUONPoint* points = new AliHLTMUONPoint[count];
	try
	{
		DebugMsg(6, "AliHLTMUONTrackerProxy::ReturnClusters");
		for (UInt i = 0; i < count; i++)
		{
			points[i] = AliHLTMUONConvert(clusters[i]);
			DebugMsg(6, "\tpoints[" << i << "] = " << points[i] );
		};
		tracker->ReturnClusters(tag, points, count);
	}
	finally
	(
		delete [] points;
	)
}


void AliHLTMUONTrackerProxy::EndOfClusters(void* tag)
{
	DebugMsg(6, "AliHLTMUONTrackerProxy::EndOfClusters");
	tracker->EndOfClusters(tag);
}


void AliHLTMUONTrackerProxy::FillTrackData(AliHLTMUONCoreTrack& track)
{
	AliHLTMUONTrack data;
	tracker->FillTrackData(data);
	DebugMsg(6, "AliHLTMUONTrackerProxy::FillTrackData : data = " << data);
	track = AliHLTMUONConvert(data);
}


void AliHLTMUONTrackerProxy::Reset()
{
	DebugMsg(6, "AliHLTMUONTrackerProxy::Reset");
	tracker->Reset();
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

