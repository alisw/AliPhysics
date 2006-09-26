////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONTRACKERPROXY_H
#define ALIHLTMUONTRACKERPROXY_H

#include "Tracking/Tracker.hpp"
#include "AliRoot/TrackerCallback.hpp"
#include "AliRoot/TrackerInterface.hpp"


class AliHLTMUONTrackerProxy : public AliHLTMUONCoreTracker, public AliHLTMUONTrackerCallback
{
public:

	AliHLTMUONTrackerProxy(AliHLTMUONTrackerInterface* client);
	virtual ~AliHLTMUONTrackerProxy() {};

	// inherited methods from Tracking::Tracker:
	virtual void FindTrack(const AliHLTMUONCoreTriggerRecord& trigger);
	virtual void ReturnClusters(void* tag, const AliHLTMUONCoreClusterPoint* clusters, UInt count);
	virtual void EndOfClusters(void* tag);
	virtual void FillTrackData(AliHLTMUONCoreTrack& track);
	virtual void Reset();

	// inherited methods from AliHLTMUONTrackerCallback:
	virtual void RequestClusters(
			Float_t left, Float_t right, Float_t bottom, Float_t top,
			Int_t chamber, const void* tag = NULL
		);
	virtual void EndOfClusterRequests();
	virtual void FoundTrack();
	virtual void NoTrackFound();

private:

	AliHLTMUONTrackerProxy(const AliHLTMUONTrackerProxy& /*object*/)
		: AliHLTMUONCoreTracker(), AliHLTMUONTrackerCallback(), fTracker(NULL)
	{}
	
	AliHLTMUONTrackerProxy& operator = (const AliHLTMUONTrackerProxy& /*object*/)
	{
		return *this;
	}


	AliHLTMUONTrackerInterface* fTracker;  // The tracker we are proxying for.
};


#endif // ALIHLTMUONTRACKERPROXY_H
