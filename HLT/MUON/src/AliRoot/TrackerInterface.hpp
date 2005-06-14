////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_TRACKER_INTERFACE_HPP
#define dHLT_ALIROOT_TRACKER_INTERFACE_HPP

#include "Rtypes.h"


namespace AliMUONHLT
{

class Point;
class TriggerRecord;
class Track;
class Tracker;
class TrackerCallback;


class TrackerInterface
{
public:
	TrackerInterface(Tracker* tracker)
	{
		fTracker = tracker;
	};
	
	const Tracker* GetTracker() const
	{
		return fTracker;
	};
	
	void FindTrack(const TriggerRecord& trigger);
	void ReturnClusters(void* tag, const Point* clusters, const UInt_t count);
	void EndOfClusters(void* tag);
	void FillTrackData(Track& track);
	void Reset();
	void SetCallback(TrackerCallback* callback);

private:

	Tracker* fTracker;   //! Pointer to interpreted tracker class.
};


}; // AliMUONHLT

#endif // dHLT_ALIROOT_TRACKER_INTERFACE_HPP
