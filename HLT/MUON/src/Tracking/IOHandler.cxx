////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Tracking/IOHandler.hpp"
#include "Error.hpp"

#ifdef DEBUG
#	include "Debug/print.hpp"
#endif

namespace dHLT
{
namespace Tracking
{


IOHandler::IOHandler()
{
	framework = NULL;
	trackerstack = NULL;
};

////////////////////////////////////////////////////////////////////////////////
//                              IOInterface 

void IOHandler::AddTriggers(const EventID event, const TriggerRecord* triggers, const UInt count)
{
	EventHandler* handler = eventtable.Find(event);
	if ( handler == NULL )
	{
		handler = eventtable.Add(event);
		handler->SetIOHandler(this);
	};
	handler->AddTriggers(triggers, count);
};


void IOHandler::EndOfTriggers(const EventID event)
{
	EventHandler* handler = eventtable.Find(event);
	if ( handler != NULL )
	{
		handler->EndOfTriggers();
	}
	else
	{
		// handle error.
	};
};


void IOHandler::AddClusters(const EventID event, const ROI region, const ClusterPoint* clusters, const UInt count)
{
	EventHandler* handler = eventtable.Find(event);
	if ( handler != NULL )
	{
		handler->AddClusters(region, clusters, count);
	}
	else
	{
		// handle error.
	};
};


void IOHandler::EndOfClusters(const EventID event, const ROI region)
{
	EventHandler* handler = eventtable.Find(event);
	if ( handler != NULL )
	{
		handler->EndOfClusters(region);
	}
	else
	{
		// handle error.
	};
};


////////////////////////////////////////////////////////////////////////////////
//                         EventHandler callbacks

void IOHandler::RequestClusters(const EventID event, const ROI region)
{
	Assert( framework != NULL );
	framework->RequestClusters(event, region);
};

void IOHandler::EndOfClusterRequests(const EventID event)
{
	Assert( framework != NULL );
	framework->EndOfClusterRequests(event);
};

void IOHandler::ReturnTracks(const EventID event, Track* newtracks, const UInt count)
{
	Assert( framework != NULL );
	framework->ReturnTracks(event, newtracks, count);
}

void IOHandler::EndOfTracks(const EventID event)
{
	Assert( framework != NULL );
	eventtable.Remove(event);
	framework->EndOfTracks(event);
};


Track* IOHandler::AllocateTrackBlock(const UInt size)
{
	Assert( framework != NULL );
	return framework->AllocateTrackBlock(size);
};

void IOHandler::ReleaseTriggers(const TriggerRecord* triggers)
{
	Assert( framework != NULL );
	framework->ReleaseTriggers(triggers);
};

void IOHandler::ReleaseClusters(const ClusterPoint* clusters)
{
	Assert( framework != NULL );
	framework->ReleaseClusters(clusters);
};


Tracker* IOHandler::NewTracker()
{
	Assert( trackerstack != NULL );
	return trackerstack->NewTracker();
};

void IOHandler::FreeTracker(Tracker* tracker)
{
	Assert( trackerstack != NULL );
	trackerstack->FreeTracker(tracker);
};


} // Tracking
} // dHLT
