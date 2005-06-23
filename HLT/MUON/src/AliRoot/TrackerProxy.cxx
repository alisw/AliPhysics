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

namespace dHLT
{
namespace AliRoot
{


TrackerProxy::TrackerProxy(AliMUONHLT::TrackerInterface* client)
	: Tracking::Tracker(), AliMUONHLT::TrackerCallback()
{
	tracker = client;
}


void TrackerProxy::FindTrack(const TriggerRecord& trigger)
{
	AliMUONHLT::TriggerRecord rec = Convert(trigger, 0);
	DebugMsg(6, "TrackerProxy::FindTrack : rec = " << rec);
	tracker->FindTrack(rec);
}


void TrackerProxy::ReturnClusters(void* tag, const ClusterPoint* clusters, const UInt count)
{
	AliMUONHLT::Point* points = new AliMUONHLT::Point[count];
	try
	{
		DebugMsg(6, "TrackerProxy::ReturnClusters");
		for (UInt i = 0; i < count; i++)
		{
			points[i] = Convert(clusters[i]);
			DebugMsg(6, "\tpoints[" << i << "] = " << points[i] );
		};
		tracker->ReturnClusters(tag, points, count);
	}
	finally
	(
		delete [] points;
	)
}


void TrackerProxy::EndOfClusters(void* tag)
{
	DebugMsg(6, "TrackerProxy::EndOfClusters");
	tracker->EndOfClusters(tag);
}


void TrackerProxy::FillTrackData(Track& track)
{
	AliMUONHLT::Track data;
	tracker->FillTrackData(data);
	DebugMsg(6, "TrackerProxy::FillTrackData : data = " << data);
	track = Convert(data);
}


void TrackerProxy::Reset()
{
	DebugMsg(6, "TrackerProxy::Reset");
	tracker->Reset();
}


void TrackerProxy::RequestClusters(
		const Float_t left, const Float_t right, const Float_t bottom, const Float_t top,
		const Int_t chamber, const void* tag
	)
{
	DebugMsg(6, "TrackerProxy::RequestClusters");
	Tracking::Tracker::RequestClusters(left, right, bottom, top, (ChamberID)chamber, tag);
}


void TrackerProxy::EndOfClusterRequests()
{
	DebugMsg(6, "TrackerProxy::EndOfClusterRequests");
	Tracking::Tracker::EndOfClusterRequests();
}


void TrackerProxy::FoundTrack()
{
	DebugMsg(6, "TrackerProxy::FoundTrack");
	Tracking::Tracker::FoundTrack();
}


void TrackerProxy::NoTrackFound()
{
	DebugMsg(6, "TrackerProxy::NoTrackFound");
	Tracking::Tracker::NoTrackFound();
}


} // AliRoot
} // dHLT
