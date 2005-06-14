////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_TRACKING_IO_HANDLER_HPP
#define dHLT_TRACKING_IO_HANDLER_HPP

#include <vector>
#include <algorithm>

#include "Tracking/Tracker.hpp"
#include "Tracking/IOInterface.hpp"
#include "Tracking/EventHandler.hpp"

#include "Buffers/RegionTree.hpp"
#include "Buffers/List.hpp"
#include "Buffers/LookupTable.hpp"


/*
ostream& operator << (ostream& os, const dHLT::EventID& id)
{
	os << "[" << id.bunch << ", " << id.timestamp << "]";
	return os;
};
*/


namespace dHLT
{
namespace Tracking
{


using dHLT::Buffers::List;
using dHLT::Buffers::LookupTable;


class TrackerStack
{
public:

	virtual Tracker* NewTracker() = 0;
	virtual void FreeTracker(Tracker* tracker) = 0;
};


////////////////////////////////////////////////////////////////////////////////////////////


class IOHandler : public IOInterface
{
public:

	class EventLookupTable : public LookupTable<EventID, EventHandler> {};

	///////// IOInterface ///////////

	virtual void AddTriggers(const EventID event, const TriggerRecord* triggers, const UInt count);
	virtual void EndOfTriggers(const EventID event);
	virtual void AddClusters(const EventID event, const ROI region, const ClusterPoint* clusters, const UInt count);
	virtual void EndOfClusters(const EventID event, const ROI region);


	//////// EventHandler callbacks /////////////

	void RequestClusters(const EventID event, const ROI region);
	void EndOfClusterRequests(const EventID event);
	
	Track* AllocateTrackBlock(const UInt size);
	void ReturnTracks(const EventID event, Track* newtracks, const UInt count);
	void EndOfTracks(const EventID event);

	void ReleaseTriggers(const TriggerRecord* triggers);
	void ReleaseClusters(const ClusterPoint* clusters);

	Tracker* NewTracker();
	void FreeTracker(Tracker* tracker);

	////////////////////////////////////

	IOHandler();
	virtual ~IOHandler() {};

	void SetCallback(IOCallback* callback) { framework = callback; };
	void SetTrackerStack(TrackerStack* stack) { trackerstack = stack; };

	UInt GetMaxTracksInBlock()
	{
		return 100;
	};


private:

	IOCallback* framework;
	TrackerStack* trackerstack;

	EventLookupTable eventtable;

};

} // Tracking
} // dHLT

#endif // dHLT_TRACKING_IO_HANDLER_HPP
