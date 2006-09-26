////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONDUMMYTRACKER_H
#define ALIHLTMUONDUMMYTRACKER_H

#ifndef __CINT__
#include "AliRoot/Point.hpp"
#include "AliRoot/TriggerRecord.hpp"
#include "AliRoot/Track.hpp"
#include "AliRoot/TrackerCallback.hpp"
#endif // __CINT__


class AliHLTMUONDummyTracker
{
public:

	AliHLTMUONDummyTracker() : fInterface(this), fCallback(NULL) {};

	AliHLTMUONDummyTracker(const AliHLTMUONDummyTracker& t)
		: fInterface(this), fCallback(t.fCallback)
	{};

	AliHLTMUONDummyTracker& operator = (const AliHLTMUONDummyTracker& t)
	{
		fCallback = t.fCallback;
		return *this;
	}

	virtual ~AliHLTMUONDummyTracker() {};

	/* Methods required to be implemented by the tracker.
	   These correspond to the AliHLTMUONCoreTracker specification, refer to that
	   class for more information.
	 */
	virtual void FindTrack(const AliHLTMUONTriggerRecord& trigger) = 0;
	virtual void ReturnClusters(void* tag, const AliHLTMUONPoint* clusters, const UInt_t count) = 0;
	virtual void EndOfClusters(void* tag) = 0;
	virtual void FillTrackData(AliHLTMUONTrack& track) = 0;
	virtual void Reset() = 0;

	/* Set the callback for the tracker, so that the tracker can communicate
	   with the framework.
	 */
	void SetCallback(AliHLTMUONTrackerCallback* callback)
	{
		fCallback = callback;
	};
	
	/* Returns the TrackerInterface object to this tracker.
	   This is required by the MicrodHLT object.
	 */
	AliHLTMUONTrackerInterface* Interface()
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

	AliHLTMUONTrackerInterface fInterface;  // The interface via which compiled code communicates with this object.
	AliHLTMUONTrackerCallback* fCallback;   // Callback interface to framework.
};


// Implementation of the TrackerInterface:
// This must come here so that it gets interpreted together with the rest
// of the AliMUONHLTTracker.

void AliHLTMUONTrackerInterface::FindTrack(const AliHLTMUONTriggerRecord& trigger)
{
	fTracker->FindTrack(trigger);
};

void AliHLTMUONTrackerInterface::ReturnClusters(void* tag, const AliHLTMUONPoint* clusters, const UInt_t count)
{
	fTracker->ReturnClusters(tag, clusters, count);
};

void AliHLTMUONTrackerInterface::EndOfClusters(void* tag)
{
	fTracker->EndOfClusters(tag);
};

void AliHLTMUONTrackerInterface::FillTrackData(AliHLTMUONTrack& track)
{
	fTracker->FillTrackData(track);
};

void AliHLTMUONTrackerInterface::Reset()
{
	fTracker->Reset();
};

void AliHLTMUONTrackerInterface::SetCallback(AliHLTMUONTrackerCallback* callback)
{
	fTracker->SetCallback(callback);
};


// Implementation of the SetTracker method which is undefined in MicrodHLT.
void AliHLTMUONMicrodHLT::SetTracker(AliHLTMUONDummyTracker* tracker)
{
	SetTracker(tracker->Interface());
};


#endif // ALIHLTMUONDUMMYTRACKER_H
