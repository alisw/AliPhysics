////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_TRACKER_HPP
#define dHLT_ALIROOT_TRACKER_HPP

#ifndef __CINT__
#include "AliRoot/Point.hpp"
#include "AliRoot/TriggerRecord.hpp"
#include "AliRoot/Track.hpp"
#include "AliRoot/TrackerCallback.hpp"
#endif // __CINT__


namespace AliMUONHLT
{


class Tracker
{
public:

	Tracker() : fInterface(this)
	{
		fCallback = NULL;
	};

	/* Methods required to be implemented by the tracker.
	   These correspond to the dHLT::Tracker specification, refer to that
	   class for more information.
	 */
	virtual void FindTrack(const TriggerRecord& trigger) = 0;
	virtual void ReturnClusters(void* tag, const Point* clusters, const UInt_t count) = 0;
	virtual void EndOfClusters(void* tag) = 0;
	virtual void FillTrackData(Track& track) = 0;
	virtual void Reset() = 0;

	/* Set the callback for the tracker, so that the tracker can communicate
	   with the framework.
	 */
	void SetCallback(TrackerCallback* callback)
	{
		fCallback = callback;
	};
	
	/* Returns the TrackerInterface object to this tracker.
	   This is required by the MicrodHLT object.
	 */
	TrackerInterface* Interface() const
	{
		return &fInterface;
	};


protected:

	void RequestClusters(
			const Float_t left, const Float_t right, const Float_t bottom, const Float_t top,
			const Int_t chamber, const void* tag
		)
	{
		if (left > right)
			Error("RequestClusters", "The parameter left (%f) is larger than right (%f).",
				left, right
			);
		else if (bottom > top)
			Error("RequestClusters", "The parameter bottom (%f) is larger than top (%f).",
				bottom, top
			);
		else if (chamber < 0 or 9 < chamber)
			Error("RequestClusters", "The chamber parameter is out of range. Got: %d, expected a value in [0..9]",
				chamber
			);
		else if (fCallback != NULL)
			fCallback->RequestClusters(left, right, bottom, top, chamber, tag);
		else
			Error("RequestClusters", "Callback not set.");
	};


	void EndOfClusterRequests()
	{
		if (fCallback != NULL)
			fCallback->EndOfClusterRequests();
		else
			Error("EndOfClusterRequests", "Callback not set.");
	};


	void FoundTrack()
	{
		if (fCallback != NULL)
			fCallback->FoundTrack();
		else
			Error("FoundTrack", "Callback not set.");
	};


	void NoTrackFound()
	{
		if (fCallback != NULL)
			fCallback->NoTrackFound();
		else
			Error("NoTrackFound", "Callback not set.");
	};


private:

	TrackerInterface fInterface;  // The interface via which compiled code communicates with this object.
	TrackerCallback* fCallback;   // Callback interface to framework.
};


// Implementation of the TrackerInterface:
// This must come here so that it gets interpreted together with the rest
// of the AliMUONHLT::Tracker.

void TrackerInterface::FindTrack(const TriggerRecord& trigger)
{
	fTracker->FindTrack(trigger);
};

void TrackerInterface::ReturnClusters(void* tag, const Point* clusters, const UInt_t count)
{
	fTracker->ReturnClusters(tag, clusters, count);
};

void TrackerInterface::EndOfClusters(void* tag)
{
	fTracker->EndOfClusters(tag);
};

void TrackerInterface::FillTrackData(Track& track)
{
	fTracker->FillTrackData(track);
};

void TrackerInterface::Reset()
{
	fTracker->Reset();
};

void TrackerInterface::SetCallback(TrackerCallback* callback)
{
	fTracker->SetCallback(callback);
};


// Implementation of the SetTracker method which is undefined in MicrodHLT.
void MicrodHLT::SetTracker(Tracker* tracker)
{
	SetTracker(tracker->Interface());
};


}; // AliMUONHLT

#endif // dHLT_ALIROOT_TRACKER_HPP
