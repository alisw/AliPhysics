////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONTRACKERINTERFACE_H
#define ALIHLTMUONTRACKERINTERFACE_H

#include "Rtypes.h"

class AliHLTMUONPoint;
class AliHLTMUONTriggerRecord;
class AliHLTMUONTrack;
class AliHLTMUONTrackerCallback;
class AliHLTMUONDummyTracker;


class AliHLTMUONTrackerInterface
{
public:
	AliHLTMUONTrackerInterface(AliHLTMUONDummyTracker* tracker) : fTracker(tracker)
	{
		fTracker = tracker;
	};
	
	const AliHLTMUONDummyTracker* GetTracker() const
	{
		return fTracker;
	};
	
	void FindTrack(const AliHLTMUONTriggerRecord& trigger);
	void ReturnClusters(void* tag, const AliHLTMUONPoint* clusters, UInt_t count);
	void EndOfClusters(void* tag);
	void FillTrackData(AliHLTMUONTrack& track);
	void Reset();
	void SetCallback(AliHLTMUONTrackerCallback* callback);

private:
	AliHLTMUONTrackerInterface(const AliHLTMUONTrackerInterface& /*tracker*/)
		: fTracker(NULL)
	{}

	AliHLTMUONTrackerInterface& operator = (const AliHLTMUONTrackerInterface& /*tracker*/)
	{
		return *this;
	}

	AliHLTMUONDummyTracker* fTracker;   //! Pointer to interpreted tracker class.
};


#endif // ALIHLTMUONTRACKERINTERFACE_H
