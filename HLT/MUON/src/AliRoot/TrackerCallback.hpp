////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_TRACKER_CALLBACK_HPP
#define dHLT_ALIROOT_TRACKER_CALLBACK_HPP

#include "TObject.h"


namespace AliMUONHLT
{


class TrackerCallback : public TObject
{
public:

	virtual void RequestClusters(
			const Float_t left, const Float_t right, const Float_t bottom, const Float_t top,
			const Int_t chamber, const void* tag = NULL
		) = 0;
	virtual void EndOfClusterRequests() = 0;
	virtual void FoundTrack() = 0;
	virtual void NoTrackFound() = 0;
	
	ClassDef(TrackerCallback, 0)  // Abstract tracker callback class.
};


} // AliMUONHLT

#endif // dHLT_ALIROOT_TRACKER_CALLBACK_HPP
